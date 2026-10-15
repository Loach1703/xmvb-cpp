#include "mcscf/ci/factors_internal.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <memory>
#include <stdexcept>

#include "core/openmp.hpp"
#include "core/timing.hpp"

namespace xmvb::mcscf {
namespace {

using Clock = std::chrono::steady_clock;

/** @brief Stream X B and X^T A through shared source-pair moments. */
class FactorMoments {
public:
  FactorMoments(const SpinFactorSource& source, int width, int columns)
      : source_(source), width_(width), a_(source.n_alpha,width),
        b_(source.n_beta,width), moment_(width,columns) {
    if (width_ == source_.count) source_.load(0,width_,a_,b_);
  }

  void add(const Eigen::Ref<const Eigen::MatrixXd>& alpha,
      const Eigen::Ref<const Eigen::MatrixXd>& beta,
      Eigen::Ref<Eigen::MatrixXd> result_alpha,
      Eigen::Ref<Eigen::MatrixXd> result_beta, double scale = 1.0) {
    for (int first = 0; first < source_.count; first += width_) {
      const int count = std::min(width_,source_.count-first);
      auto a = a_.leftCols(count), b = b_.leftCols(count);
      if (width_ != source_.count) source_.load(first,count,a,b);
      auto moment = moment_.topLeftCorner(count,alpha.cols());
      detail::transpose_product(b,beta,moment);
      if (scale != 1.0) moment *= scale;
      detail::add_product(a,moment,result_alpha);
      detail::transpose_product(a,alpha,moment);
      if (scale != 1.0) moment *= scale;
      detail::add_product(b,moment,result_beta);
    }
  }

private:
  const SpinFactorSource& source_;
  int width_;
  Eigen::MatrixXd a_, b_, moment_;
};

Eigen::MatrixXd generator_adjoint(const DirectCiSigmaAction& sigma,
    const SpinFactorSource& left, const SpinFactorSource& right,
    bool hamiltonian, double energy, std::size_t workspace_bytes) {
  profiling::ScopedPhase phase(profiling::Phase::CiGeneratorAdjoint);
  detail::validate_source(sigma,left); detail::validate_source(sigma,right);
  if (!std::isfinite(energy)) throw std::invalid_argument("nonfinite residual energy");
  const int n = sigma.n_orbitals(), d = sigma.channel_kernel().rows();
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(n,n);
  if (!left.count || !right.count) return result;
  // Only an identity bilinear permits exchanging sources followed by a
  // transpose. H on the left cannot be moved through an oriented generator.
  const bool transpose = !hamiltonian && left.count < right.count;
  const auto& bra = transpose ? right : left;
  const auto& ket = transpose ? left : right;
  const std::size_t rows = static_cast<std::size_t>(left.n_alpha) + left.n_beta;
  const std::size_t capacity = workspace_bytes / sizeof(double);
  const std::size_t reduction = static_cast<std::size_t>(
      effective_openmp_thread_count()+1) * n * n;
  const auto payload = [&](std::size_t l, std::size_t r, std::size_t p) {
    const long double source = static_cast<long double>(rows) * (l + (hamiltonian ? 3 : 2)*r);
    const long double moments = static_cast<long double>(l)*r*std::max<std::size_t>(1,p);
    const long double loading = std::max(static_cast<long double>(bra.workspace_rows)*l,
                                         static_cast<long double>(ket.workspace_rows)*r);
    const long double sweep = hamiltonian ? static_cast<long double>(rows)*r*(3*p+6) : 0;
    return source + moments + loading + sweep + reduction;
  };
  int p = hamiltonian ? d : 0;
  if (payload(1,1,hamiltonian ? 1 : 0) > capacity)
    throw std::length_error("factor generator workspace cannot hold one source/channel panel");
  if (hamiltonian && payload(1,1,p) > capacity) {
    int low = 1, high = d-1;
    while (low < high) {
      const int middle = low + (high-low+1)/2;
      if (payload(1,1,middle) <= capacity) low = middle; else high = middle-1;
    }
    p = low;
  }
  const bool resident = payload(bra.count,1,p) <= capacity;
  int low = 1, high = hamiltonian && p < d ? 1 : ket.count;
  while (low < high) {
    const int middle = low + (high-low+1)/2;
    if (payload(resident ? bra.count : middle,middle,p) <= capacity) low = middle;
    else high = middle-1;
  }
  const int nr = low;
  low = nr; high = bra.count;
  // H-left sources are not rank-sorted. A smaller bra is perfectly valid.
  low = std::min(low,high);
  while (low < high) {
    const int middle = low + (high-low+1)/2;
    if (payload(middle,nr,p) <= capacity) low = middle; else high = middle-1;
  }
  const int nl = low;
  const bool timing = profiling::phase_totals().enabled;
  const auto start = timing ? Clock::now() : Clock::time_point{};
  FactorMoments moments(bra,nl,nr*std::max(1,p));
  Eigen::MatrixXd a(left.n_alpha,nr), b(left.n_beta,nr);
  Eigen::MatrixXd wa(left.n_alpha,nr), wb(left.n_beta,nr);
  Eigen::MatrixXd ca(left.n_alpha,hamiltonian ? nr : 0), cb(left.n_beta,hamiltonian ? nr : 0);
  const SpinFactorSource shape{left.n_alpha,left.n_beta,nr,{}};
  std::unique_ptr<FactorActionWorkspace> sweep;
  if (hamiltonian)
    sweep = std::make_unique<FactorActionWorkspace>(shape,d,rows*nr*(3*p+6)*sizeof(double));
  for (int first = 0; first < ket.count; first += nr) {
    const int width = std::min(nr,ket.count-first);
    auto aa = a.leftCols(width), bb = b.leftCols(width);
    ket.load(first,width,aa,bb);
    auto weighted_a = wa.leftCols(width), weighted_b = wb.leftCols(width);
    weighted_a.setZero(); weighted_b.setZero();
    moments.add(aa,bb,weighted_a,weighted_b);
    if (!hamiltonian) {
      result += sigma.spin_generator_adjoint(false,weighted_a,aa);
      result += sigma.spin_generator_adjoint(true,weighted_b,bb);
      continue;
    }
    auto out_a = ca.leftCols(width), out_b = cb.leftCols(width);
    out_a = -energy * weighted_a; out_b = -energy * weighted_b;
    const SpinFactorSource weighted{left.n_alpha,left.n_beta,width,
        [&](int offset, int count, Eigen::Ref<Eigen::MatrixXd> x, Eigen::Ref<Eigen::MatrixXd> y) {
          x = weighted_a.middleCols(offset,count); y = weighted_b.middleCols(offset,count);
        }};
    sweep->apply(sigma,weighted,[](FactorCrossBlock) {}, [&](FactorImageBlock block) {
      out_a.middleCols(block.source,block.width) += block.same_alpha;
      out_b.middleCols(block.source,block.width) += block.same_beta;
    });
    const SpinFactorSource factors{left.n_alpha,left.n_beta,width,
        [&](int offset, int count, Eigen::Ref<Eigen::MatrixXd> x, Eigen::Ref<Eigen::MatrixXd> y) {
          x = aa.middleCols(offset,count); y = bb.middleCols(offset,count);
        }};
    sweep->apply(sigma,factors,[&](FactorCrossBlock block) {
      // Same-spin gather is complete: reuse raw channel slabs for X V_sigma.
      block.alpha_channels.setZero(); block.beta_channels.setZero();
      moments.add(block.mixed_alpha,block.mixed_beta,block.alpha_channels,block.beta_channels,2.0);
      sigma.gather_spin(false,block.alpha_channels,block.channel,block.channels,
          out_a.middleCols(block.source,block.width));
      sigma.gather_spin(true,block.beta_channels,block.channel,block.channels,
          out_b.middleCols(block.source,block.width));
    },[&](FactorImageBlock block) {
      moments.add(block.same_alpha,block.same_beta,
          out_a.middleCols(block.source,block.width),out_b.middleCols(block.source,block.width));
    });
    result += sigma.spin_generator_adjoint(false,out_a,aa);
    result += sigma.spin_generator_adjoint(true,out_b,bb);
  }
  if (timing)
    std::fprintf(stderr,"XMVB_FACTOR_GENERATOR,hamiltonian=%d,bra=%d,ket=%d,left_tile=%d,right_tile=%d,channels=%d,workspace_bytes=%.0Lf,seconds=%.9f\n",
        hamiltonian,bra.count,ket.count,nl,nr,p,sizeof(double)*payload(nl,nr,p),
        std::chrono::duration<double>(Clock::now()-start).count());
  if (transpose) return result.transpose();
  return result;
}

} // namespace

Eigen::MatrixXd factor_generator_adjoint(const DirectCiSigmaAction& sigma,
    const SpinFactorSource& left, const SpinFactorSource& right, std::size_t workspace_bytes) {
  return generator_adjoint(sigma,left,right,false,0.0,workspace_bytes);
}

Eigen::MatrixXd factor_residual_generator_adjoint(const DirectCiSigmaAction& sigma,
    const SpinFactorSource& left, const SpinFactorSource& right, double energy,
    std::size_t workspace_bytes) {
  return generator_adjoint(sigma,left,right,true,energy,workspace_bytes);
}

} // namespace xmvb::mcscf
