#include "tnvb/gradient/adapter.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "core/pair_index.hpp"

namespace xmvb::tnvb {

void populate_active_space_adjoints(const WorkerResponse& response,
                                   vb::ActiveSpaceGradientResult* result) {
  if (result == nullptr) throw std::invalid_argument("TNVB adjoint result is null");
  const Eigen::MatrixXd& h = result->active_space_one_electron_result.h1e_act;
  const Eigen::MatrixXd& p = response.density_one;
  const int n_ao = h.rows();
  const std::size_t expected = std::size_t(n_ao)*n_ao*n_ao*n_ao;
  if (n_ao <= 0 || h.cols() != n_ao || p.rows() != n_ao || p.cols() != n_ao ||
      response.contracted_three.rows() != n_ao || response.contracted_three.cols() != n_ao ||
      response.density_two_buffer.size() != expected || !p.allFinite() ||
      !response.contracted_three.allFinite() || response.energies.size() == 0) {
    throw std::invalid_argument("TNVB density/adjoint dimensions or values are invalid");
  }
  for (double entry : response.density_two_buffer) {
    if (!std::isfinite(entry)) throw std::invalid_argument("nonfinite TNVB two-body density");
  }
  const auto pi = [&](int i, int j, int k, int l) {
    return response.density_two_buffer[i+std::size_t(n_ao)*(j+std::size_t(n_ao)*(k+std::size_t(n_ao)*l))];
  };
  // Norm derivative -Eact P excludes inactive and nuclear energies.
  Eigen::MatrixXd overlap_adjoint = 0.5*response.contracted_three-response.energies[0]*p;
  for (int j=0; j<n_ao; ++j) for (int i=0; i<n_ao; ++i) {
    for (int u=0; u<n_ao; ++u) for (int v=0; v<n_ao; ++v) {
      overlap_adjoint(i,j) += h(u,v)*pi(v,i,u,j);
    }
  }
  overlap_adjoint = (0.5*(overlap_adjoint+overlap_adjoint.transpose())).eval();
  result->active_orbital_overlap_gradient.assign(overlap_adjoint.data(),
      overlap_adjoint.data()+overlap_adjoint.size());
  result->active_one_electron_gradient.assign(p.data(), p.data()+p.size());
  const std::size_t n_pair=std::size_t(n_ao)*(n_ao+1)/2;
  result->packed_active_two_electron_gradient.assign(n_pair*(n_pair+1)/2, 0.0);
  // Enumerating ordered quadruples counts all packing permutations once,
  // including the reduced multiplicities when indices coincide.
  for (int l=0; l<n_ao; ++l) for (int k=0; k<n_ao; ++k)
    for (int j=0; j<n_ao; ++j) for (int i=0; i<n_ao; ++i) {
      const auto index=core::PairIndex::two_electron_storage_index(i,j,k,l);
      result->packed_active_two_electron_gradient[index] += 0.5*pi(l,j,k,i);
    }
}

}  // namespace xmvb::tnvb
