#include "integrals/response/ri_ao.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <cblas.h>
#include "parallel/openmp.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace xmvb::integrals {

RiAoResponse apply_ri_ao_response(
    const RiActiveTwoElectronResponseCache& active,
    const Eigen::Ref<const Eigen::MatrixXd>& inactive_duals,
    const Eigen::Ref<const Eigen::MatrixXd>& horizontal_direction,
    const Eigen::Ref<const Eigen::MatrixXd>& weighted_active,
    RiAoResponseWorkspace& workspace,std::size_t workspace_bytes) {
  const int n=active.n_basis_functions, na=active.n_active_orbitals;
  const int ni=inactive_duals.cols(), rank=ni+na, n_aux=active.n_auxiliary_functions;
  if (n<=0 || na<=0 || n_aux<=0 || inactive_duals.rows()!=n ||
      horizontal_direction.rows()!=n || horizontal_direction.cols()!=ni ||
      weighted_active.rows()!=n || weighted_active.cols()!=na ||
      !active.accepted_active_coefficients || !active.metric_whitened_ao_pair_factors)
    throw std::invalid_argument("invalid factor-native RI AO response dimensions");
  // Previous slabs are not inputs to this action. Release them before a
  // different rank/tile plan allocates its first returned transform.
  workspace.direction.resize(0,0);
  workspace.accepted.resize(0,0);
  const auto& c=*active.accepted_active_coefficients;
  // The response source is Z=U V^T+V U^T, with no numerical rank discovery.
  Eigen::MatrixXd u(n,rank),v(n,rank);
  u.leftCols(ni)=2*horizontal_direction;
  u.rightCols(na)=0.5*weighted_active;
  v.leftCols(ni)=inactive_duals;
  v.rightCols(na)=c;
  const int ao_workers=std::min(effective_openmp_thread_count(),n);
  const std::size_t limit=workspace_bytes/sizeof(double);
  const std::size_t fixed=2ULL*n*rank+n_aux+static_cast<std::size_t>(n)*(n+1)/2;
  const std::size_t per_worker=3ULL*n*na+2ULL*na*na;
  // Account for both live slabs and a returned transform while the previous
  // slab still exists, plus AO row slices owned by the shared transform team.
  const std::size_t per_aux=static_cast<std::size_t>(n)*(3*rank+ao_workers);
  if (fixed+per_worker+per_aux>limit)
    throw std::length_error("minimum factor-native RI AO tile exceeds workspace");
  const int workers=std::min<std::size_t>(std::min(ao_workers,n_aux),
      (limit-fixed-per_aux)/per_worker);
  const int tile=std::min<std::size_t>(n_aux,(limit-fixed-workers*per_worker)/per_aux);
  workspace.workers.resize(workers);
  if (ni) for (auto& worker:workspace.workers) {
    worker.left.resize(n,na); worker.right.resize(n,na);
    worker.projection.setZero(n,na);
    worker.core.resize(na,na); worker.mixed.resize(na,na);
  }
  workspace.coulomb.resize(n_aux);
  RiAoResponse result{Eigen::MatrixXd::Zero(n,n),Eigen::MatrixXd::Zero(n,na)};
  const int panels=std::min(n,static_cast<int>(std::ceil(
      (std::sqrt(1.0+8.0*workers)-1)/2)));
  std::vector<std::pair<int,int>> blocks;
  for (int row=0;row<panels;++row)
    for (int col=0;col<=row;++col) blocks.emplace_back(row,col);

  for (int first=0;first<n_aux;first+=tile) {
    const int count=std::min(tile,n_aux-first),leading=count*rank;
    workspace.direction=transform_ri_one_index(*active.metric_whitened_ao_pair_factors,u,first,count);
    if (active.accepted_one_index_factors) {
      workspace.accepted.resize(leading,n);
      if (ni)
        workspace.accepted.topRows(ni*count)=
            transform_ri_one_index(*active.metric_whitened_ao_pair_factors,inactive_duals,first,count);
      // Borrow the already-computed active transform, never recompute L C_a.
      for (int a=0;a<na;++a)
        workspace.accepted.middleRows((ni+a)*count,count)=
            active.accepted_one_index_factors->middleRows(a*n_aux+first,count);
    } else {
      workspace.accepted=transform_ri_one_index(*active.metric_whitened_ao_pair_factors,v,first,count);
    }
#pragma omp parallel num_threads(workers)
    {
      int thread=0;
#ifdef _OPENMP
      thread=omp_get_thread_num();
#endif
      auto& scratch=workspace.workers[thread];
      const int begin=thread*count/workers,end=(thread+1)*count/workers;
      auto charge=workspace.coulomb.segment(first+begin,end-begin);
      charge.setZero();
      // Contiguous auxiliary segments, not an auxiliary-major strided scan
      // through the same slab for every scalar charge.
      for (int mu=0;mu<n;++mu) for (int a=0;a<rank;++a)
        if (u(mu,a)!=0)
          charge.noalias()+=2*u(mu,a)*workspace.accepted.col(mu).segment(a*count+begin,end-begin);
      if (ni) {
#pragma omp for schedule(static)
      for (int aux=0;aux<count;++aux) {
        for (int mu=0;mu<n;++mu) {
          const double* accepted=workspace.accepted.col(mu).data()+aux;
          const double* direction=workspace.direction.col(mu).data()+aux;
          for (int a=0;a<na;++a) {
            scratch.left(mu,a)=direction[(ni+a)*count];
            scratch.right(mu,a)=accepted[(ni+a)*count];
          }
        }
        // Only the narrow G(K) C correction is needed, not the AO matrix G(K).
        scratch.core.noalias()=c.transpose()*scratch.right;
        scratch.mixed.noalias()=weighted_active.transpose()*scratch.right;
        const double weight=4*(weighted_active.array()*scratch.right.array()).sum();
        scratch.projection.noalias()+=weight*scratch.right;
        scratch.projection.noalias()-=2*scratch.left*scratch.core;
        scratch.projection.noalias()-=scratch.right*scratch.mixed;
      }
      }
    }
    // Every lower AO rectangle has one writer; no AO-square worker reductions.
#pragma omp parallel for schedule(dynamic,1) num_threads(workers)
    for (std::size_t block=0;block<blocks.size();++block) {
      const auto [row,col]=blocks[block];
      const int rb=row*n/panels,re=(row+1)*n/panels;
      const int cb=col*n/panels,ce=(col+1)*n/panels;
      double* output=result.combined.data()+static_cast<std::size_t>(cb)*n+rb;
      const double* ur=workspace.direction.data()+static_cast<std::size_t>(rb)*leading;
      const double* vr=workspace.accepted.data()+static_cast<std::size_t>(rb)*leading;
      const double* uc=workspace.direction.data()+static_cast<std::size_t>(cb)*leading;
      const double* vc=workspace.accepted.data()+static_cast<std::size_t>(cb)*leading;
      if (row==col) {
        cblas_dsyr2k(CblasColMajor,CblasLower,CblasTrans,re-rb,leading,
            -1,ur,leading,vr,leading,1,output,n);
      } else {
        cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,re-rb,ce-cb,leading,
            -1,ur,leading,vc,leading,1,output,n);
        cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,re-rb,ce-cb,leading,
            -1,vr,leading,uc,leading,1,output,n);
      }
    }
  }
  const Eigen::VectorXd coulomb=2*active.metric_whitened_ao_pair_factors->transpose()*
      workspace.coulomb;
  for (int mu=0;mu<n;++mu) for (int nu=0;nu<=mu;++nu) {
    result.combined(mu,nu)+=coulomb[static_cast<std::size_t>(mu)*(mu+1)/2+nu];
    result.combined(nu,mu)=result.combined(mu,nu);
  }
  if (ni) {
    for (const auto& worker:workspace.workers) result.projection-=0.25*worker.projection;
    result.projection.noalias()+=0.5*result.combined*c;
  }
  return result;
}

} // namespace xmvb::integrals
