#include "mcscf/ci/factors.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <stdexcept>
#include <utility>

#include "core/openmp.hpp"
#include "core/timing.hpp"

namespace xmvb::mcscf {
namespace {

using Clock = std::chrono::steady_clock;

/** @brief Mix channel columns without replicating or coupling source labels. */
void mix_channels(const Eigen::Ref<const Eigen::MatrixXd>& channels,
    const Eigen::MatrixXd& kernel, int p, int q, int width,
    Eigen::Ref<Eigen::MatrixXd> image) {
  const int rows = channels.rows() * width;
  const int np = image.cols() / width, nq = channels.cols() / width;
  Eigen::Map<const Eigen::MatrixXd> input(channels.data(), rows, nq);
  Eigen::Map<Eigen::MatrixXd> output(image.data(), rows, np);
  const int workers = coarse_openmp_thread_count(effective_openmp_thread_count(),
      static_cast<std::size_t>(rows) * np * nq);
#pragma omp parallel for schedule(static) if(workers > 1) num_threads(workers)
  for (int part = 0; part < workers; ++part) {
    const int first = rows * part / workers, last = rows * (part + 1) / workers;
    if (last > first)
      output.middleRows(first, last - first).noalias() +=
          input.middleRows(first, last - first) * kernel.block(p, q, np, nq).transpose();
  }
}

/** @brief Sum two spin-channel Grams, with disjoint shared output row ownership. */
void add_channel_gram(const Eigen::Ref<const Eigen::MatrixXd>& ca,
    const Eigen::Ref<const Eigen::MatrixXd>& cb,
    const Eigen::Ref<const Eigen::MatrixXd>& qa,
    const Eigen::Ref<const Eigen::MatrixXd>& qb, int width,
    Eigen::Ref<Eigen::MatrixXd> result) {
  const int np = result.rows(), nq = result.cols();
  const int na = ca.rows() * width, nb = cb.rows() * width;
  const Eigen::Map<const Eigen::MatrixXd> a(ca.data(), na, np), b(cb.data(), nb, np);
  const Eigen::Map<const Eigen::MatrixXd> x(qa.data(), na, nq), y(qb.data(), nb, nq);
  const int workers = coarse_openmp_thread_count(effective_openmp_thread_count(),
      static_cast<std::size_t>(na + nb) * np * nq);
#pragma omp parallel for schedule(static) if(workers > 1) num_threads(workers)
  for (int part = 0; part < workers; ++part) {
    const int first = np * part / workers, last = np * (part + 1) / workers;
    if (last == first) continue;
    result.middleRows(first, last - first).noalias() += a.middleCols(first, last - first).transpose() * x;
    result.middleRows(first, last - first).noalias() += b.middleCols(first, last - first).transpose() * y;
  }
}

void validate_source(const DirectCiSigmaAction& sigma, const SpinFactorSource& source) {
  if (source.n_alpha != sigma.n_alpha_strings() || source.n_beta != sigma.n_beta_strings() ||
      source.count < 0 || !source.load)
    throw std::invalid_argument("invalid direct-CI spin-factor source");
}

/** @brief Accumulate one product by shared output rows, with single-thread BLAS. */
void add_product(const Eigen::Ref<const Eigen::MatrixXd>& left,
    const Eigen::Ref<const Eigen::MatrixXd>& right, Eigen::Ref<Eigen::MatrixXd> result) {
  const int workers = coarse_openmp_thread_count(effective_openmp_thread_count(),
      static_cast<std::size_t>(result.rows()) * result.cols() * left.cols());
#pragma omp parallel for schedule(static) if(workers > 1) num_threads(workers)
  for (int part = 0; part < workers; ++part) {
    const int first = result.rows() * part / workers, last = result.rows() * (part + 1) / workers;
    if (last > first)
      result.middleRows(first, last - first).noalias() +=
          left.middleRows(first, last - first) * right;
  }
}

/** @brief Source/channel moments in one product, partitioned by output rows. */
void transpose_product(const Eigen::Ref<const Eigen::MatrixXd>& channels,
    const Eigen::Ref<const Eigen::MatrixXd>& ket, Eigen::Ref<Eigen::MatrixXd> result) {
  const bool by_rows = result.rows() >= result.cols();
  const int extent = by_rows ? result.rows() : result.cols();
  const int workers = coarse_openmp_thread_count(effective_openmp_thread_count(),
      static_cast<std::size_t>(channels.rows()) * channels.cols() * ket.cols());
#pragma omp parallel for schedule(static) if(workers > 1) num_threads(workers)
  for (int part = 0; part < workers; ++part) {
    const int first = extent * part / workers, last = extent * (part + 1) / workers;
    if (last == first) continue;
    if (by_rows)
      result.middleRows(first, last - first).noalias() = channels.middleCols(first, last - first).transpose() * ket;
    else
      result.middleCols(first, last - first).noalias() = channels.transpose() * ket.middleCols(first, last - first);
  }
}

} // namespace

SpinFactorSource sum_factors(std::vector<SpinFactorSource> sources) {
  if (sources.empty()) throw std::invalid_argument("empty factor sum");
  SpinFactorSource result{sources.front().n_alpha, sources.front().n_beta, 0, {}};
  for (const auto& source : sources) {
    if (source.n_alpha != result.n_alpha || source.n_beta != result.n_beta ||
        source.count < 0 || !source.load)
      throw std::invalid_argument("incompatible factor sum");
    result.count += source.count;
    result.workspace_rows = std::max(result.workspace_rows, source.workspace_rows);
  }
  result.load = [sources = std::move(sources)](int first, int width,
      Eigen::Ref<Eigen::MatrixXd> a, Eigen::Ref<Eigen::MatrixXd> b) {
    int offset = 0;
    for (const auto& source : sources) {
      if (first >= source.count) { first -= source.count; continue; }
      const int count = std::min(width - offset, source.count - first);
      source.load(first, count, a.middleCols(offset, count), b.middleCols(offset, count));
      offset += count; first = 0;
      if (offset == width) break;
    }
  };
  return result;
}

SpinFactorSource generator_factors(const DirectCiSigmaAction& sigma,
    const SpinFactorSource& source, const Eigen::Ref<const Eigen::MatrixXd>& generator) {
  validate_source(sigma, source);
  return {source.n_alpha, source.n_beta, 2 * source.count,
      [&sigma, source, generator](int first, int width,
          Eigen::Ref<Eigen::MatrixXd> a, Eigen::Ref<Eigen::MatrixXd> b) {
        int offset = 0;
        while (offset < width) {
          const bool beta = first >= source.count;
          const int local = beta ? first - source.count : first;
          const int count = std::min(width - offset, source.count - local);
          auto aa = a.middleCols(offset, count), bb = b.middleCols(offset, count);
          source.load(local, count, aa, bb);
          Eigen::MatrixXd image(beta ? source.n_beta : source.n_alpha, count);
          sigma.apply_spin_generator(beta, beta ? bb : aa, generator, image);
          if (beta) bb = image; else aa = image;
          first += count; offset += count;
        }
      }, std::max(source.workspace_rows,
          static_cast<std::size_t>(std::max(source.n_alpha, source.n_beta)))};
}

FactorActionWorkspace::FactorActionWorkspace(const SpinFactorSource& source, int channels,
    std::size_t workspace_bytes, bool overlap_only)
    : n_alpha_(source.n_alpha), n_beta_(source.n_beta), channels_(channels),
      load_rows_(source.workspace_rows), overlap_only_(overlap_only) {
  if (n_alpha_ <= 0 || n_beta_ <= 0 || source.count < 0 || channels <= 0)
    throw std::invalid_argument("invalid factor workspace dimensions");
  const std::size_t unit_bytes = (static_cast<std::size_t>(n_alpha_) + n_beta_) * sizeof(double);
  const std::size_t load_bytes = source.workspace_rows * sizeof(double);
  if (workspace_bytes < unit_bytes * (overlap_only ? 2U : 9U) + load_bytes)
    throw std::length_error("source workspace cannot hold one single-spin channel tile");
  const int columns_per_source = overlap_only ? 2 : 3 * channels + 6;
  source_tile_ = static_cast<int>(std::max<std::size_t>(1,
      std::min<std::size_t>(source.count, workspace_bytes / (unit_bytes * columns_per_source + load_bytes))));
  channel_tile_ = overlap_only ? 0 : static_cast<int>(std::min<std::size_t>(channels,
      ((workspace_bytes / source_tile_ - load_bytes) / unit_bytes - 6) / 3));
  alpha_.resize(n_alpha_, source_tile_); beta_.resize(n_beta_, source_tile_);
  ha_.resize(n_alpha_, overlap_only ? 0 : source_tile_);
  hb_.resize(n_beta_, overlap_only ? 0 : source_tile_);
  ap_.resize(n_alpha_, source_tile_ * channel_tile_);
  bp_.resize(n_beta_, source_tile_ * channel_tile_);
  va_.resize(n_alpha_, source_tile_ * channel_tile_);
  vb_.resize(n_beta_, source_tile_ * channel_tile_);
  if (!overlap_only && channel_tile_ < channels) {
    aq_.resize(n_alpha_, source_tile_ * channel_tile_);
    bq_.resize(n_beta_, source_tile_ * channel_tile_);
  }
}

FactorActionWork FactorActionWorkspace::apply(const DirectCiSigmaAction& sigma,
    const SpinFactorSource& source, const std::function<void(FactorCrossBlock)>& cross,
    const std::function<void(FactorImageBlock)>& finish) {
  validate_source(sigma, source);
  if (source.n_alpha != n_alpha_ || source.n_beta != n_beta_ ||
      source.workspace_rows > load_rows_ || sigma.channel_kernel().rows() != channels_ ||
      !finish || (!overlap_only_ && !cross))
    throw std::invalid_argument("factor action does not match its workspace or consumers");
  FactorActionWork work;
  work.source_tile = source_tile_; work.channel_tile = channel_tile_;
  const bool timing = profiling::phase_totals().enabled;
  const auto tick = [&] { return timing ? Clock::now() : Clock::time_point{}; };
  const auto seconds = [](auto a, auto b) { return std::chrono::duration<double>(b-a).count(); };
  const auto& kernel = sigma.channel_kernel();
  for (int first = 0; first < source.count; first += source_tile_) {
    const int width = std::min(source_tile_, source.count - first);
    ++work.panels; work.sources += width;
    work.workspace_bound = std::max(work.workspace_bound, sizeof(double) *
        ((static_cast<std::size_t>(n_alpha_) + n_beta_) * width *
        (overlap_only_ ? 2 : 3 * channel_tile_ + 6) + source.workspace_rows * width));
    auto a = alpha_.leftCols(width), b = beta_.leftCols(width);
    auto ha = ha_.leftCols(overlap_only_ ? 0 : width);
    auto hb = hb_.leftCols(overlap_only_ ? 0 : width);
    source.load(first, width, a, b);
    ha.setZero(); hb.setZero();
    if (!overlap_only_) {
      for (int p = 0; p < channels_; p += channel_tile_) {
        const int np = std::min(channel_tile_, channels_ - p);
        auto ap = ap_.leftCols(np * width), bp = bp_.leftCols(np * width);
        auto va = va_.leftCols(np * width), vb = vb_.leftCols(np * width);
        const auto scattered = tick();
        sigma.scatter_spin(false, a, p, np, ap);
        sigma.scatter_spin(true, b, p, np, bp);
        work.scatter_seconds += seconds(scattered, tick());
        va.setZero(); vb.setZero();
        for (int q = 0; q < channels_; q += channel_tile_) {
          ++work.channel_blocks;
          const int nq = std::min(channel_tile_, channels_ - q);
          if (q == p) {
            const auto start = tick();
            mix_channels(ap, kernel, p, q, width, va);
            mix_channels(bp, kernel, p, q, width, vb);
            work.mix_seconds += seconds(start, tick());
          } else {
            const auto start = tick();
            auto aq = aq_.leftCols(nq * width), bq = bq_.leftCols(nq * width);
            sigma.scatter_spin(false, a, q, nq, aq);
            sigma.scatter_spin(true, b, q, nq, bq);
            const auto mixed = tick();
            work.scatter_seconds += seconds(start, mixed);
            mix_channels(aq, kernel, p, q, width, va);
            mix_channels(bq, kernel, p, q, width, vb);
            work.mix_seconds += seconds(mixed, tick());
          }
        }
        const auto gathered = tick();
        sigma.gather_spin(false, va, p, np, ha);
        sigma.gather_spin(true, vb, p, np, hb);
        work.gather_seconds += seconds(gathered, tick());
        cross({first, width, p, np, ap, vb});
      }
    }
    finish({first, width, a, b, ha, hb});
  }
  return work;
}

DirectCiIntegralAdjoint factor_integral_adjoint(const DirectCiSigmaAction& sigma,
    const SpinFactorSource& left, const SpinFactorSource& right, std::size_t workspace_bytes) {
  profiling::ScopedPhase phase(profiling::Phase::CiIntegralAdjoint);
  validate_source(sigma, left); validate_source(sigma, right);
  const int d = sigma.channel_kernel().rows();
  Eigen::MatrixXd gradient = Eigen::MatrixXd::Zero(d, d);
  if (!left.count || !right.count) return sigma.pullback_channel_kernel(gradient);
  // H is self-adjoint. Put the smaller factor family on the ket side, where
  // it controls the leading channel-Gram work, without changing the functional.
  const auto& bra = left.count >= right.count ? left : right;
  const auto& ket = left.count >= right.count ? right : left;
  const std::size_t rows = static_cast<std::size_t>(left.n_alpha) + left.n_beta;
  const std::size_t capacity = workspace_bytes / sizeof(double);
  const auto payload = [&](std::size_t l, std::size_t r, std::size_t p) {
    const int channel_families = p < static_cast<std::size_t>(d) ? 3 : 2;
    return static_cast<long double>(rows) * (l + 2 * r + channel_families * r * p) +
        static_cast<long double>(l) * r * p +
        std::max(static_cast<long double>(bra.workspace_rows) * l,
                 static_cast<long double>(ket.workspace_rows) * r);
  };
  if (payload(1, 1, 1) > capacity)
    throw std::length_error("factor adjoint workspace cannot hold one channel/source panel");
  int p = d;
  if (payload(1, 1, p) > capacity) {
    int low = 1, high = p - 1;
    while (low < high) {
      const int middle = low + (high - low + 1) / 2;
      if (payload(1, 1, middle) <= capacity) low = middle; else high = middle - 1;
    }
    p = low;
  }
  // When possible, retain the complete bra within the SAME slab budget.
  // This leaves the leading Gram/moment arithmetic unchanged but removes bra
  // regeneration across ket/channel panels. In particular a directional source
  // must not redo its exterior maps and one-body generators on every ket tile.
  const bool resident_bra = payload(bra.count, 1, p) <= capacity;
  int low = 1, high = ket.count;
  while (low < high) {
    const int middle = low + (high - low + 1) / 2;
    if (payload(resident_bra ? bra.count : middle, middle, p) <= capacity) low = middle;
    else high = middle - 1;
  }
  const int nr = low;
  low = nr; high = bra.count;
  while (low < high) {
    const int middle = low + (high - low + 1) / 2;
    if (payload(middle, nr, p) <= capacity) low = middle; else high = middle - 1;
  }
  const int nl = low;
  const bool timing = profiling::phase_totals().enabled;
  const auto tick = [&] { return timing ? Clock::now() : Clock::time_point{}; };
  const auto seconds = [](auto a, auto b) { return std::chrono::duration<double>(b-a).count(); };
  double load_time = 0, scatter_time = 0, covector_time = 0, gram_time = 0;
  Eigen::MatrixXd la(left.n_alpha, nl), lb(left.n_beta, nl);
  Eigen::MatrixXd ra(left.n_alpha, nr), rb(left.n_beta, nr);
  Eigen::MatrixXd wa(left.n_alpha, nr), wb(left.n_beta, nr);
  Eigen::MatrixXd rp_a(left.n_alpha, nr * p), rp_b(left.n_beta, nr * p);
  Eigen::MatrixXd rq_a(left.n_alpha, p < d ? nr * p : 0);
  Eigen::MatrixXd rq_b(left.n_beta, p < d ? nr * p : 0);
  Eigen::MatrixXd ca(left.n_alpha, nr * p), cb(left.n_beta, nr * p);
  Eigen::MatrixXd moments(nl, nr * p);
  // Retain a whole bra family only when it fits the same explicit slab budget.
  // Otherwise the provider streams it through the very same contraction.
  const bool bra_resident = nl == bra.count;
  if (bra_resident) {
    const auto start = tick();
    bra.load(0, nl, la, lb);
    load_time += seconds(start, tick());
  }
  for (int r = 0; r < ket.count; r += nr) {
    const int wr = std::min(nr, ket.count - r);
    auto a = ra.leftCols(wr), b = rb.leftCols(wr);
    const auto load_start = tick();
    ket.load(r, wr, a, b);
    load_time += seconds(load_start, tick());
    auto weighted_a = wa.leftCols(wr), weighted_b = wb.leftCols(wr);
    weighted_a.setZero(); weighted_b.setZero();
    // Move self-adjoint generators through the spectator contraction:
    // C_alpha,P = E_alpha,P(Y_left F_right) + Y_left(E_beta,P F_right).
    // The weighted bases are formed once, with no bra-channel scatter.
    for (int l = 0; l < bra.count; l += nl) {
      const int wl = std::min(nl, bra.count - l);
      auto u = la.leftCols(wl), v = lb.leftCols(wl);
      if (!bra_resident) {
        const auto start = tick(); bra.load(l, wl, u, v);
        load_time += seconds(start, tick());
      }
      const auto start = tick();
      auto moment = moments.topLeftCorner(wl, wr);
      transpose_product(v, b, moment); add_product(u, moment, weighted_a);
      transpose_product(u, a, moment); add_product(v, moment, weighted_b);
      covector_time += seconds(start, tick());
    }
    for (int first = 0; first < d; first += p) {
      const int np = std::min(p, d - first);
      auto bar_a = ca.leftCols(wr * np), bar_b = cb.leftCols(wr * np);
      auto pa = rp_a.leftCols(wr * np), pb = rp_b.leftCols(wr * np);
      const auto scattered = tick();
      sigma.scatter_spin(false, weighted_a, first, np, bar_a);
      sigma.scatter_spin(true, weighted_b, first, np, bar_b);
      sigma.scatter_spin(false, a, first, np, pa);
      sigma.scatter_spin(true, b, first, np, pb);
      scatter_time += seconds(scattered, tick());
      // Accumulate the entire bra before the Q-channel Gram: splitting bra
      // sources must not multiply the leading d^2 contraction by its tile count.
      for (int l = 0; l < bra.count; l += nl) {
        const int wl = std::min(nl, bra.count - l);
        auto u = la.leftCols(wl), v = lb.leftCols(wl);
        if (!bra_resident) {
          const auto start = tick(); bra.load(l, wl, u, v);
          load_time += seconds(start, tick());
        }
        const auto covectors = tick();
        auto moment = moments.topLeftCorner(wl, wr * np);
        transpose_product(v, pb, moment); add_product(u, moment, bar_a);
        transpose_product(u, pa, moment); add_product(v, moment, bar_b);
        covector_time += seconds(covectors, tick());
      }
      for (int second = 0; second < d; second += p) {
        const int nq = std::min(p, d - second);
        if (second == first) {
          const auto start = tick();
          add_channel_gram(bar_a, bar_b, pa, pb, wr, gradient.block(first, second, np, nq));
          gram_time += seconds(start, tick());
        } else {
          auto qa = rq_a.leftCols(wr * nq), qb = rq_b.leftCols(wr * nq);
          const auto scatter_start = tick();
          sigma.scatter_spin(false, a, second, nq, qa);
          sigma.scatter_spin(true, b, second, nq, qb);
          const auto gram_start = tick();
          scatter_time += seconds(scatter_start, gram_start);
          add_channel_gram(bar_a, bar_b, qa, qb, wr, gradient.block(first, second, np, nq));
          gram_time += seconds(gram_start, tick());
        }
      }
    }
  }
  if (timing)
    std::fprintf(stderr, "XMVB_FACTOR_ADJOINT,bra=%d,ket=%d,left_tile=%d,right_tile=%d,channel_tile=%d,workspace_bytes=%.0Lf,load=%.9f,scatter=%.9f,covectors=%.9f,gram=%.9f\n",
        bra.count, ket.count, nl, nr, p, sizeof(double) * payload(nl, nr, p),
        load_time, scatter_time, covector_time, gram_time);
  return sigma.pullback_channel_kernel(gradient);
}

} // namespace xmvb::mcscf
