#include "vbscf/structures/assembly/action.hpp"

#include <algorithm>
#include <cmath>
#include <chrono>
#include <stdexcept>

#include "core/openmp.hpp"
#include "core/timing.hpp"
#include "mcscf/ci/product_support.hpp"
#include "mcscf/ci/sigma.hpp"
#include "vbscf/structures/assembly/sources.hpp"
#include "vbscf/structures/orthogonal_ci/exterior_transform.hpp"

namespace xmvb::vb {
namespace {

/** @brief Mix channels without coupling or replicating the source columns. */
void mix_channels(const Eigen::Ref<const Eigen::MatrixXd>& channels, const Eigen::MatrixXd& kernel,
    int p, int q, int width, Eigen::Ref<Eigen::MatrixXd> image) {
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

} // namespace

bool StructureAction::prefer_source_contraction(const ExteriorOrbitalTransform& alpha,
    const ExteriorOrbitalTransform& beta, const mcscf::DirectCiSigmaAction& sigma) const {
  const long double r = sources_->representatives.size();
  const long double na = n_unique_alpha_, nb = n_unique_beta_;
  const long double pairs = sigma.channel_kernel().rows();
  const int n = alpha.n_orbitals(), ea = alpha.n_electrons(), eb = beta.n_electrons();
  // Arithmetic counts for two channel GEMMs, single-spin links, exterior
  // shears and sampling. This compares two exact contraction orders, not
  // two physical spaces; neither a memory budget nor a numerical threshold enters.
  const long double factor_work = r * (
      2 * (na + nb) * pairs * pairs +
      4 * (na * ea * (n - ea + 1) + nb * eb * (n - eb + 1)) +
      (pairs + 3) * (2 * (alpha.shear_pair_count() + beta.shear_pair_count()) + na + nb) +
      (2 * pairs + 6) * spin_products_.size());
  const long double full_rows = spin_exchange_parity_ == 0 ? na * nb :
      na * (na + spin_exchange_parity_) / 2;
  const long double channel_work = 2 * pairs * pairs;
  if (factor_work >= full_rows * channel_work) return false;
  if (!product_support_) return true;

  // Count only generated support until its GEMM work alone exceeds the
  // factor estimate. Stop at that proof; no full-product membership scan.
  auto cursor = product_support_->cursor(1);
  int a, b;
  bool supported;
  long double rows = 0;
  while (cursor.next(a, b, supported)) {
    if (spin_exchange_parity_ && (a > b || (a == b && spin_exchange_parity_ < 0))) continue;
    if (++rows * channel_work > factor_work) return true;
  }
  return false;
}

StructureActionResult StructureAction::apply_source_tiles(
    const Eigen::Ref<const Eigen::MatrixXd>& vectors,
    const ExteriorOrbitalTransform& alpha, const ExteriorOrbitalTransform& beta,
    const mcscf::DirectCiSigmaAction& sigma, bool overlap_only,
    std::size_t workspace_bytes) const {
  profiling::ScopedPhase phase(profiling::Phase::CiSigma);
  using Clock = std::chrono::steady_clock;
  const bool timing = profiling::phase_totals().enabled;
  const auto tick = [&] { return timing ? Clock::now() : Clock::time_point{}; };
  const auto seconds = [](auto first, auto last) { return std::chrono::duration<double>(last-first).count(); };
  double scatter_time = 0, mix_time = 0, gather_time = 0, pullback_time = 0, sample_time = 0;
  const auto& sources = *sources_;
  const auto& kernel = sigma.channel_kernel();
  const auto alpha_projection = alpha.projection(sources.alpha_rows);
  const auto beta_projection = beta.projection(sources.columns);
  const int pairs = kernel.rows();
  const std::size_t spin_rows = static_cast<std::size_t>(n_unique_alpha_) + n_unique_beta_;
  // At most three channel pairs plus six source/image pairs coexist. Tile
  // source columns first, then channels if even one full channel family is too big.
  const std::size_t unit_bytes = spin_rows * sizeof(double);
  if (workspace_bytes / unit_bytes < (overlap_only ? 2U : 9U))
    throw std::length_error("source workspace cannot hold one single-spin channel tile");
  const int columns_per_source = overlap_only ? 2 : 3 * pairs + 6;
  const int source_tile = static_cast<int>(std::max<std::size_t>(1,
      std::min<std::size_t>(sources.representatives.size(), workspace_bytes / unit_bytes / columns_per_source)));
  const int channel_tile = overlap_only ? 0 : static_cast<int>(std::min<std::size_t>(pairs,
      (workspace_bytes / unit_bytes / source_tile - 6) / 3));
  StructureActionResult result{
      Eigen::MatrixXd::Zero(n_structures_, vectors.cols()),
      Eigen::MatrixXd::Zero(n_structures_, vectors.cols())};
  // One shared allocation per slab. Active full-row views also handle final
  // partial source/channel tiles without resizing or refaulting large buffers.
  Eigen::MatrixXd a_buffer(n_unique_alpha_, source_tile), b_buffer(n_unique_beta_, source_tile);
  const int image_columns = overlap_only ? 0 : source_tile;
  Eigen::MatrixXd ha_buffer(n_unique_alpha_, image_columns), hb_buffer(n_unique_beta_, image_columns);
  Eigen::MatrixXd ap_buffer(n_unique_alpha_, source_tile * channel_tile);
  Eigen::MatrixXd bp_buffer(n_unique_beta_, source_tile * channel_tile);
  Eigen::MatrixXd va_buffer(n_unique_alpha_, source_tile * channel_tile);
  Eigen::MatrixXd vb_buffer(n_unique_beta_, source_tile * channel_tile);
  Eigen::MatrixXd aq_buffer, bq_buffer;
  if (!overlap_only && channel_tile < pairs) {
    aq_buffer.resize(n_unique_alpha_, source_tile * channel_tile);
    bq_buffer.resize(n_unique_beta_, source_tile * channel_tile);
  }
  for (int vector = 0; vector < vectors.cols(); ++vector) {
    Eigen::VectorXd coefficients = Eigen::VectorXd::Zero(spin_products_.size());
    Eigen::VectorXd h_products = coefficients, m_products = coefficients;
    ++result.tiles.expansions;
#pragma omp parallel for schedule(static) if(spin_products_.size() > 32)
    for (int group = 0; group < static_cast<int>(spin_products_.size()); ++group)
      for (std::size_t t = spin_term_offsets_[group]; t < spin_term_offsets_[group + 1]; ++t)
        coefficients[group] += spin_terms_[t].coefficient * vectors(spin_terms_[t].structure, vector);

    for (int first = 0; first < static_cast<int>(sources.representatives.size()); first += source_tile) {
      const int width = std::min<int>(source_tile, sources.representatives.size() - first);
      ++result.tiles.source_tiles;
      result.tiles.source_columns += width;
      result.tiles.source_workspace_bound = std::max(result.tiles.source_workspace_bound,
          unit_bytes * width * (overlap_only ? 2 : 3 * channel_tile + 6));
      auto a = a_buffer.leftCols(width), b = b_buffer.leftCols(width);
      auto ha = ha_buffer.leftCols(overlap_only ? 0 : width);
      auto hb = hb_buffer.leftCols(overlap_only ? 0 : width);
      a.setZero(); b.setZero(); ha.setZero(); hb.setZero();
      for (int j = 0; j < width; ++j) {
        const int source = sources.representatives[first + j];
        b(source, j) = 1.0;
        for (int group : sources.beta_groups[source]) {
          const int row = spin_products_[group] / n_unique_beta_;
          if (spin_exchange_parity_ && (row > source || (row == source && spin_exchange_parity_ < 0))) continue;
          a(row, j) = coefficients[group] * (spin_exchange_parity_ && row != source ? 2 : 1);
        }
      }
      alpha.apply_left_block(a);
      beta.apply_left_block(b);
      if (!overlap_only) {
        for (int p = 0; p < pairs; p += channel_tile) {
          const int np = std::min(channel_tile, pairs - p);
          auto ap = ap_buffer.leftCols(np * width), bp = bp_buffer.leftCols(np * width);
          auto va = va_buffer.leftCols(np * width), vb = vb_buffer.leftCols(np * width);
          const auto scatter_start = tick();
          sigma.scatter_spin(false, a, p, np, ap);
          sigma.scatter_spin(true, b, p, np, bp);
          scatter_time += seconds(scatter_start, tick());
          va.setZero(); vb.setZero();
          for (int q = 0; q < pairs; q += channel_tile) {
            ++result.tiles.channel_blocks;
            const int nq = std::min(channel_tile, pairs - q);
            if (q == p) {
              const auto start = tick();
              mix_channels(ap, kernel, p, q, width, va);
              mix_channels(bp, kernel, p, q, width, vb);
              mix_time += seconds(start, tick());
            } else {
              const auto start = tick();
              auto aq = aq_buffer.leftCols(nq * width), bq = bq_buffer.leftCols(nq * width);
              sigma.scatter_spin(false, a, q, nq, aq);
              sigma.scatter_spin(true, b, q, nq, bq);
              const auto mixed = tick();
              scatter_time += seconds(start, mixed);
              mix_channels(aq, kernel, p, q, width, va);
              mix_channels(bq, kernel, p, q, width, vb);
              mix_time += seconds(mixed, tick());
            }
          }
          const auto gathered = tick();
          sigma.gather_spin(false, va, p, np, ha);
          sigma.gather_spin(true, vb, p, np, hb);
          const auto pulling = tick();
          gather_time += seconds(gathered, pulling);
          // Commuting spin sectors give the factor two, since K is symmetric.
          // Sample only original products after the two single-spin adjoints.
          alpha.apply_adjoint_left_block(ap, alpha_projection);
          beta.apply_adjoint_left_block(vb, beta_projection);
          const auto sampled = tick();
          pullback_time += seconds(pulling, sampled);
#pragma omp parallel for schedule(static) if(spin_products_.size() > 32)
          for (int group = 0; group < static_cast<int>(spin_products_.size()); ++group) {
            const int i = spin_products_[group] / n_unique_beta_;
            const int j = spin_products_[group] % n_unique_beta_;
            h_products[group] += 2 * ap.row(i).dot(vb.row(j));
          }
          sample_time += seconds(sampled, tick());
        }
        alpha.apply_adjoint_left_block(ha, alpha_projection);
        beta.apply_adjoint_left_block(hb, beta_projection);
      }
      alpha.apply_adjoint_left_block(a, alpha_projection);
      beta.apply_adjoint_left_block(b, beta_projection);
#pragma omp parallel for schedule(static) if(spin_products_.size() > 32)
      for (int group = 0; group < static_cast<int>(spin_products_.size()); ++group) {
        const int i = spin_products_[group] / n_unique_beta_;
        const int j = spin_products_[group] % n_unique_beta_;
        m_products[group] += a.row(i).dot(b.row(j));
        if (!overlap_only)
          h_products[group] += ha.row(i).dot(b.row(j)) + a.row(i).dot(hb.row(j));
      }
    }
    ++result.tiles.structure_projections;
#pragma omp parallel for schedule(static) if(n_structures_ > 32)
    for (int structure = 0; structure < n_structures_; ++structure) {
      double h = 0, m = 0;
      for (std::size_t t = structure_term_offsets_[structure]; t < structure_term_offsets_[structure + 1]; ++t) {
        const auto& term = spin_terms_[structure_term_indices_[t]];
        const auto group = std::lower_bound(spin_products_.begin(), spin_products_.end(),
            term.spin_product) - spin_products_.begin();
        h += term.coefficient * h_products[group];
        m += term.coefficient * m_products[group];
      }
      result.hamiltonian(structure, vector) = h;
      result.overlap(structure, vector) = m;
    }
  }
  if (profiling::phase_totals().enabled)
    std::fprintf(stderr, "XMVB_SOURCES,structures=%d,sources=%zu,products=%zu,source_tile=%d,channel_tile=%d\n",
        n_structures_, sources.representatives.size(), spin_products_.size(), source_tile, channel_tile);
  if (timing)
    std::fprintf(stderr, "XMVB_SOURCE_TIME,scatter=%.9f,mix=%.9f,gather=%.9f,pullback=%.9f,sample=%.9f\n",
        scatter_time, mix_time, gather_time, pullback_time, sample_time);
  return result;
}

} // namespace xmvb::vb
