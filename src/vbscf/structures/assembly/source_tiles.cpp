#include "vbscf/structures/assembly/action.hpp"

#include <algorithm>
#include <chrono>

#include "core/openmp.hpp"
#include "core/timing.hpp"
#include "mcscf/ci/factors.hpp"
#include "mcscf/ci/product_support.hpp"
#include "mcscf/ci/sigma.hpp"
#include "vbscf/structures/assembly/sources.hpp"
#include "vbscf/structures/orthogonal_ci/exterior_transform.hpp"

namespace xmvb::vb {
Eigen::VectorXd StructureAction::expand_support(
    const Eigen::Ref<const Eigen::VectorXd>& vectors) const {
  Eigen::VectorXd coefficients = Eigen::VectorXd::Zero(spin_products_.size());
#pragma omp parallel for schedule(static) if(spin_products_.size() > 32)
  for (int group = 0; group < static_cast<int>(spin_products_.size()); ++group)
    for (std::size_t t = spin_term_offsets_[group]; t < spin_term_offsets_[group + 1]; ++t)
      coefficients[group] += spin_terms_[t].coefficient * vectors[spin_terms_[t].structure];
  return coefficients;
}

bool StructureAction::prefer_source_contraction(const ExteriorOrbitalTransform& alpha,
    const ExteriorOrbitalTransform& beta, const mcscf::DirectCiSigmaAction& sigma) const {
  const long double r = sources_->nodes.size();
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
  return apply_source_columns(vectors.cols(), [&](int column) {
    return expand_support(vectors.col(column));
  }, alpha, beta, sigma, overlap_only, workspace_bytes);
}

StructureActionResult StructureAction::apply_source_columns(int width,
    const std::function<Eigen::VectorXd(int)>& support_column,
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
  const auto& alpha_projection = sources.alpha_projection;
  const auto& beta_projection = sources.beta_projection;
  const mcscf::SpinFactorSource shape{n_unique_alpha_, n_unique_beta_,
      static_cast<int>(sources.nodes.size()), {}};
  mcscf::FactorActionWorkspace workspace(shape, sigma.channel_kernel().rows(), workspace_bytes, overlap_only);
  int source_tile = 0, channel_tile = 0;
  StructureActionResult result{
      Eigen::MatrixXd::Zero(n_structures_, width),
      Eigen::MatrixXd::Zero(n_structures_, width)};
  for (int vector = 0; vector < width; ++vector) {
    const Eigen::VectorXd coefficients = support_column(vector);
    Eigen::VectorXd h_products = Eigen::VectorXd::Zero(spin_products_.size()), m_products = h_products;
    ++result.tiles.expansions;

    const auto source = sources.orthogonal_source(spin_products_, coefficients, alpha, beta);
    const auto cross = [&](mcscf::FactorCrossBlock block) {
      const auto pulling = tick();
      alpha.apply_adjoint_left_block(block.alpha_channels, alpha_projection);
      beta.apply_adjoint_left_block(block.mixed_beta, beta_projection);
      const auto sampled = tick();
      pullback_time += seconds(pulling, sampled);
#pragma omp parallel for schedule(static) if(spin_products_.size() > 32)
      for (int group = 0; group < static_cast<int>(spin_products_.size()); ++group) {
        const int i = spin_products_[group] / n_unique_beta_;
        const int j = spin_products_[group] % n_unique_beta_;
        h_products[group] += 2 * block.alpha_channels.row(i).dot(block.mixed_beta.row(j));
      }
      sample_time += seconds(sampled, tick());
    };
    const auto finish = [&](mcscf::FactorImageBlock block) {
      if (!overlap_only) {
        alpha.apply_adjoint_left_block(block.same_alpha, alpha_projection);
        beta.apply_adjoint_left_block(block.same_beta, beta_projection);
      }
      alpha.apply_adjoint_left_block(block.alpha, alpha_projection);
      beta.apply_adjoint_left_block(block.beta, beta_projection);
#pragma omp parallel for schedule(static) if(spin_products_.size() > 32)
      for (int group = 0; group < static_cast<int>(spin_products_.size()); ++group) {
        const int i = spin_products_[group] / n_unique_beta_;
        const int j = spin_products_[group] % n_unique_beta_;
        m_products[group] += block.alpha.row(i).dot(block.beta.row(j));
        if (!overlap_only)
          h_products[group] += block.same_alpha.row(i).dot(block.beta.row(j)) +
              block.alpha.row(i).dot(block.same_beta.row(j));
      }
    };
    const auto work = workspace.apply(sigma, source, cross, finish);
    source_tile = work.source_tile; channel_tile = work.channel_tile;
    result.tiles.source_tiles += work.panels;
    result.tiles.source_columns += work.sources;
    result.tiles.channel_blocks += work.channel_blocks;
    result.tiles.source_workspace_bound = std::max(result.tiles.source_workspace_bound, work.workspace_bound);
    scatter_time += work.scatter_seconds; mix_time += work.mix_seconds; gather_time += work.gather_seconds;
    ++result.tiles.structure_projections;
#pragma omp parallel for schedule(static) if(n_structures_ > 32)
    for (int structure = 0; structure < n_structures_; ++structure) {
      double h = 0, m = 0;
      for (std::size_t t = structure_term_offsets_[structure]; t < structure_term_offsets_[structure + 1]; ++t) {
        const auto& term = spin_terms_[structure_term_indices_[t]];
        h += term.coefficient * h_products[term.group];
        m += term.coefficient * m_products[term.group];
      }
      result.hamiltonian(structure, vector) = h;
      result.overlap(structure, vector) = m;
    }
  }
  if (profiling::phase_totals().enabled)
    std::fprintf(stderr, "XMVB_SOURCES,structures=%d,sources=%zu,products=%zu,source_tile=%d,channel_tile=%d\n",
        n_structures_, sources.nodes.size(), spin_products_.size(), source_tile, channel_tile);
  if (timing)
    std::fprintf(stderr, "XMVB_SOURCE_TIME,scatter=%.9f,mix=%.9f,gather=%.9f,pullback=%.9f,sample=%.9f\n",
        scatter_time, mix_time, gather_time, pullback_time, sample_time);
  return result;
}

} // namespace xmvb::vb
