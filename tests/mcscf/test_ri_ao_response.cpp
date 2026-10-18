#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <Eigen/Core>
#include "integrals/response/ri_ao.hpp"
#include "parallel/openmp.hpp"

namespace {
using Matrix=Eigen::MatrixXd;

void close(const Matrix& actual,const Matrix& expected) {
  const double error=(actual-expected).cwiseAbs().maxCoeff();
  if (error>2e-11*std::max(1.0,expected.cwiseAbs().maxCoeff()))
    throw std::runtime_error("factor-native RI response differs from explicit contraction");
}

void check(int ni,bool cached,std::size_t budget) {
  constexpr int n=13,na=3,aux=19;
  Matrix factors(aux,n*(n+1)/2),c(n,na),a(n,ni),x(n,ni),w(n,na);
  for (int q=0;q<factors.rows();++q) for (int p=0;p<factors.cols();++p)
    factors(q,p)=0.1*std::sin(0.7*(p+1)*(q+2));
  for (int mu=0;mu<n;++mu) {
    for (int j=0;j<na;++j) {
      c(mu,j)=0.2*std::cos(0.3*(mu+1)*(j+1));
      w(mu,j)=0.1*std::sin(0.4*(mu+2)*(j+1));
    }
    for (int j=0;j<ni;++j) {
      a(mu,j)=0.3*std::sin(0.2*(mu+1)*(j+1));
      x(mu,j)=0.1*std::cos(0.7*(mu+2)*(j+1));
    }
  }
  xmvb::integrals::RiActiveTwoElectronResponseCache cache;
  cache.n_basis_functions=n; cache.n_active_orbitals=na; cache.n_auxiliary_functions=aux;
  cache.metric_whitened_ao_pair_factors=&factors;
  cache.accepted_active_coefficients=&c;
  std::shared_ptr<Matrix> resident;
  if (cached) resident=std::make_shared<Matrix>(aux*na,n);
  const Matrix d=x*a.transpose()+a*x.transpose(),k=w*c.transpose()+c*w.transpose();
  Matrix gd=Matrix::Zero(n,n),gk=gd;
  for (int q=0;q<aux;++q) {
    Matrix l(n,n);
    for (int mu=0;mu<n;++mu) for (int nu=0;nu<=mu;++nu)
      l(mu,nu)=l(nu,mu)=factors(q,mu*(mu+1)/2+nu);
    gd.noalias()+=2*(l.array()*d.array()).sum()*l-l*d*l;
    gk.noalias()+=2*(l.array()*k.array()).sum()*l-l*k*l;
    if (cached) {
      const Matrix lc=l*c;
      for (int j=0;j<na;++j)
        resident->row(j*aux+q)=lc.col(j).transpose();
    }
  }
  cache.accepted_one_index_factors=std::move(resident);
  xmvb::integrals::RiAoResponseWorkspace workspace;
  auto actual=xmvb::integrals::apply_ri_ao_response(cache,a,x,w,workspace,budget);
  close(actual.combined,2*gd+0.5*gk);
  close(actual.projection,gd*c);
  close(actual.combined,actual.combined.transpose());
  // Reuse must clear all per-worker projections and leave no stale charge.
  const auto zero=xmvb::integrals::apply_ri_ao_response(cache,a,0*x,0*w,workspace,budget);
  close(zero.combined,Matrix::Zero(n,n));
  close(zero.projection,Matrix::Zero(n,na));
}
} // namespace

int main() {
  for (int threads:{1,2,4}) {
#ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(threads);
#endif
    for (int ni:{0,2,5}) for (bool cached:{false,true})
      for (std::size_t budget:{std::size_t{8192},xmvb::integrals::ri_response_workspace_bytes})
        check(ni,cached,budget);
  }
  std::cout<<"ri_ao_response_exact_projection=PASS\n";
}
