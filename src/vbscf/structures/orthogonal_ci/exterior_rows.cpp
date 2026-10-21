#include "vbscf/structures/orthogonal_ci/exterior_rows.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <memory>

#include "core/timing.hpp"
#include "parallel/openmp.hpp"

namespace xmvb::vb {

ExteriorRowPlan::ExteriorRowPlan(int n_orbitals, int n_electrons,
    std::span<const std::uint64_t> masks)
    : n_(n_orbitals), electrons_(n_electrons),
      binomial_((n_+1)*(n_+1),0), column_ranks_(masks.size()) {
  for (int n=0;n<=n_;++n) {
    binomial_[n*(n_+1)]=1;
    for (int k=1;k<=n;++k)
      binomial_[n*(n_+1)+k]=choose(n-1,k-1)+choose(n-1,k);
  }
  for (std::size_t column=0;column<masks.size();++column) {
    auto mask=masks[column];
    std::uint64_t rank=0;
    for (int k=1;mask;++k) {
      const int orbital=std::countr_zero(mask);
      rank+=choose(orbital,k);
      mask&=mask-1;
    }
    column_ranks_[column]=static_cast<int>(rank);
  }
}

std::uint64_t ExteriorRowPlan::choose(int n,int k) const noexcept {
  return n<0 || k<0 || k>n ? 0:binomial_[n*(n_+1)+k];
}

ExteriorRowPlan::Rows ExteriorRowPlan::apply(const Eigen::MatrixXd& upper,
    std::span<const std::uint64_t> masks,std::span<const int> rows) const {
  profiling::ScopedPhase phase(profiling::Phase::CiTransformRight);
  const std::size_t dimension=masks.size();
  // A worker owns contiguous rows. Writing a column-major result here would
  // make neighboring workers repeatedly invalidate the same cache lines.
  Rows result(rows.size(),dimension);
  // This is derived scratch residency, not a change of representation.
  constexpr std::size_t scratch_bytes=64ULL*1024*1024;
  const int workers=bounded_openmp_thread_count(cooperative_openmp_thread_count(),
      std::min(rows.size(),std::max<std::size_t>(1,
          scratch_bytes/(2*sizeof(double)*dimension))));
  auto scratch=std::make_unique_for_overwrite<double[]>(2*dimension*workers);
  cooperative_openmp_for(workers,workers,[&](std::ptrdiff_t part) {
    std::array<int,63> occupied{},indices{};
    for (std::size_t row=rows.size()*part/workers;row<rows.size()*(part+1)/workers;++row) {
      auto mask=masks[rows[row]];
      for (int k=0;mask;++k) {
        occupied[k]=std::countr_zero(mask);
        mask&=mask-1;
      }
      double* previous=scratch.get()+2*dimension*part;
      double* current=previous+dimension;
      previous[0]=1.0;
      for (int k=1;k<=electrons_;++k) {
        const int base=electrons_-k;
        const auto count=choose(n_-base,k);
        std::fill_n(current,count,0.0);
        std::copy_n(occupied.data()+base,k,indices.data());
        for (;;) {
          std::uint64_t rank=0,suffix=0,prefix=0;
          for (int p=0;p<k;++p) {
            rank+=choose(indices[p]-base,p+1);
            if (p) suffix+=choose(indices[p]-base-1,p);
          }
          double value=0.0;
          // Laplace expansion along the newly prepended occupied row.
          for (int p=0;p<k;++p) {
            if (indices[0]!=base || p==0) {
              const double coefficient=upper(occupied[base],indices[p]);
              if (coefficient!=0.0)
                value+=(p%2 ? -coefficient:coefficient)*previous[prefix+suffix];
            }
            prefix+=choose(indices[p]-base-1,p+1);
            if (p+1<k) suffix-=choose(indices[p+1]-base-1,p+1);
          }
          current[rank]=value;
          int p=0;
          while (p<k && indices[p]+1==(p+1<k ? indices[p+1]:n_)) ++p;
          if (p==k) break;
          ++indices[p];
          std::copy_n(occupied.data()+base,p,indices.data());
        }
        std::swap(previous,current);
      }
      for (std::size_t column=0;column<dimension;++column)
        result(row,column)=previous[column_ranks_[column]];
    }
  });
  return result;
}

std::size_t ExteriorRowPlan::dynamic_bytes() const noexcept {
  return binomial_.capacity()*sizeof(std::uint64_t)+column_ranks_.capacity()*sizeof(int);
}

} // namespace xmvb::vb
