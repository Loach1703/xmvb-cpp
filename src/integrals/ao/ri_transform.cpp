#include "integrals/ao/ri_transform.hpp"

#include <algorithm>
#include <stdexcept>
#include "parallel/openmp.hpp"

namespace xmvb::integrals {

Eigen::MatrixXd transform_ri_one_index(
    const Eigen::Ref<const Eigen::MatrixXd>& packed_factors,
    const Eigen::Ref<const Eigen::MatrixXd>& coefficients,
    int first_aux,int count) {
  const int n=coefficients.rows(),orbitals=coefficients.cols();
  if (n<=0 || first_aux<0 || count<=0 || first_aux+count>packed_factors.rows() ||
      packed_factors.cols()!=static_cast<Eigen::Index>(n)*(n+1)/2)
    throw std::invalid_argument("invalid RI one-index transform dimensions");
  Eigen::MatrixXd transformed(static_cast<Eigen::Index>(count)*orbitals,n);
  const int workers=std::min(effective_openmp_thread_count(),n);
#pragma omp parallel num_threads(workers)
  {
    Eigen::MatrixXd ao_slice(count,n);
#pragma omp for schedule(static)
    for (int mu=0;mu<n;++mu) {
      for (int nu=0;nu<n;++nu) {
        const int hi=std::max(mu,nu),lo=std::min(mu,nu);
        const Eigen::Index pair=static_cast<Eigen::Index>(hi)*(hi+1)/2+lo;
        ao_slice.col(nu)=packed_factors.col(pair).segment(first_aux,count);
      }
      Eigen::Map<Eigen::MatrixXd> output(transformed.col(mu).data(),count,orbitals);
      output.noalias()=ao_slice*coefficients;
    }
  }
  return transformed;
}

} // namespace xmvb::integrals
