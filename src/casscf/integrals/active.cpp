#include "casscf/integrals/active.hpp"

#include <cmath>
#include <stdexcept>
#include "core/pair_index.hpp"

namespace xmvb::casscf {

CasIntegrals make_active_integrals(
    const Eigen::Ref<const Eigen::MatrixXd>& h, const MoCoulomb& eri,
    int core, int active, double nuclear) {
  const int n=static_cast<int>(h.rows());
  if (core<0 || active<1 || core+active>n || h.cols()!=n ||
      !h.allFinite() || !eri || !std::isfinite(nuclear))
    throw std::invalid_argument("invalid orthonormal MO integrals or CAS partition");
  CasIntegrals result;
  result.core_energy=nuclear;
  for (int i=0; i<core; ++i) {
    result.core_energy+=2*h(i,i);
    for (int j=0; j<core; ++j)
      result.core_energy+=2*eri(i,i,j,j)-eri(i,j,j,i);
  }
  result.active.one_electron=h.block(core,core,active,active);
  for (int a=0; a<active; ++a)
    for (int b=0; b<active; ++b)
      for (int i=0; i<core; ++i)
        result.active.one_electron(a,b)+=
            2*eri(core+a,core+b,i,i)-eri(core+a,i,i,core+b);
  const int active_pairs=active*(active+1)/2;
  result.active.pair_kernel.resize(active_pairs,active_pairs);
  for (int a=0; a<active; ++a)
    for (int b=0; b<=a; ++b)
      for (int c=0; c<active; ++c)
        for (int d=0; d<=c; ++d)
          result.active.pair_kernel(core::PairIndex::packed_pair_index(a,b),
              core::PairIndex::packed_pair_index(c,d))=eri(core+a,core+b,core+c,core+d);
  return result;
}

}  // namespace xmvb::casscf
