#include "vbscf/structures/assembly/action.hpp"

#include <algorithm>
#include <list>
#include <unordered_map>

#include "core/openmp.hpp"
#include "mcscf/ci/product_support.hpp"
#include "mcscf/ci/sigma.hpp"
#include "vbscf/structures/orthogonal_ci/exterior_transform.hpp"

namespace xmvb::vb {
namespace {

/** @brief Bounded LRU of column panels, without unbounded product-grid storage.
 * Eviction consumes an output panel before its storage is reused. Input
 * panels can be regenerated exactly; the budget changes residency, not H.
 */
class ColumnPanels {
public:
  using Build = std::function<Eigen::MatrixXd(int)>;
  using Consume = std::function<void(int, Eigen::MatrixXd&)>;
  using BuildMask = std::function<std::vector<std::uint8_t>(int)>;

  ColumnPanels(std::size_t bytes, Build build, Consume consume = {}, BuildMask mask = {})
      : budget_bytes_(std::max<std::size_t>(1, bytes)),
        build_(std::move(build)), consume_(std::move(consume)), mask_(std::move(mask)) {}

  Eigen::MatrixXd& get(int index) {
    auto found = index_.find(index);
    if (found != index_.end()) {
      panels_.splice(panels_.begin(), panels_, found->second);
      return found->second->values;
    }
    Eigen::MatrixXd values = build_(index);
    auto allowed = mask_ ? mask_(index) : std::vector<std::uint8_t>{};
    const auto bytes = static_cast<std::size_t>(values.size()) * sizeof(double) + allowed.capacity();
    while (!panels_.empty() &&
           (bytes > budget_bytes_ || resident_bytes_ > budget_bytes_ - bytes)) {
      auto& last = panels_.back();
      resident_bytes_ -= static_cast<std::size_t>(last.values.size()) * sizeof(double) +
                         last.allowed.capacity();
      if (consume_) consume_(last.index, last.values);
      index_.erase(last.index);
      panels_.pop_back();
      ++evictions;
    }
    panels_.push_front({index, std::move(values), std::move(allowed)});
    resident_bytes_ += bytes;
    peak_bytes = std::max(peak_bytes, resident_bytes_);
    index_[index] = panels_.begin();
    return panels_.front().values;
  }

  /** @brief Support of a resident panel; valid until the next cache mutation. */
  const std::vector<std::uint8_t>& mask(int index) const {
    return index_.at(index)->allowed;
  }

  void finish() {
    if (consume_)
      for (auto& panel : panels_) consume_(panel.index, panel.values);
    panels_.clear();
    index_.clear();
    resident_bytes_ = 0;
  }

  std::size_t evictions = 0, peak_bytes = 0;

private:
  struct Panel { int index; Eigen::MatrixXd values; std::vector<std::uint8_t> allowed; };
  std::size_t budget_bytes_;
  std::size_t resident_bytes_ = 0;
  Build build_;
  Consume consume_;
  BuildMask mask_;
  std::list<Panel> panels_;
  std::unordered_map<int, std::list<Panel>::iterator> index_;
};

}  // namespace

int StructureAction::detect_spin_exchange() const {
  int parity = 0;
  for (std::size_t group = 0; group < spin_products_.size(); ++group) {
    const auto begin = spin_term_offsets_[group], end = spin_term_offsets_[group + 1];
    if (begin == end) continue;
    const int product = spin_products_[group];
    const int exchanged = (product % n_unique_beta_) * n_unique_alpha_ +
                          product / n_unique_beta_;
    const auto found = std::lower_bound(spin_products_.begin(), spin_products_.end(), exchanged);
    if (found == spin_products_.end() || *found != exchanged) return 0;
    const auto other = static_cast<std::size_t>(found - spin_products_.begin());
    const auto paired = spin_term_offsets_[other];
    if (end - begin != spin_term_offsets_[other + 1] - paired) return 0;
    for (std::size_t offset = 0; offset < end - begin; ++offset) {
      const auto& a = spin_terms_[begin + offset];
      const auto& b = spin_terms_[paired + offset];
      if (a.structure != b.structure) return 0;
      const int sign = a.coefficient == b.coefficient ? 1 :
                      a.coefficient == -b.coefficient ? -1 : 0;
      if (!sign || (parity && parity != sign)) return 0;
      parity = sign;
    }
  }
  return parity;
}

StructureActionResult StructureAction::apply_orthogonal_tiled(
    const Eigen::Ref<const Eigen::MatrixXd>& vectors,
    const ExteriorOrbitalTransform& alpha_transform,
    const ExteriorOrbitalTransform& beta_transform,
    const mcscf::DirectCiSigmaAction& sigma,
    bool overlap_only,
    std::size_t workspace_bytes) const {
  StructureActionResult result{
      Eigen::MatrixXd::Zero(n_structures_, vectors.cols()),
      Eigen::MatrixXd::Zero(n_structures_, vectors.cols())};
  // The budget controls residency and RHS batching, never admission to CI.
  const std::size_t cache_bytes = std::max<std::size_t>(1, workspace_bytes / 4);
  const std::size_t pairs = static_cast<std::size_t>(alpha_transform.n_orbitals()) *
      (alpha_transform.n_orbitals() + 1) / 2;
  const std::size_t sigma_vector_bytes = sizeof(double) *
      static_cast<std::size_t>(n_unique_alpha_) * (2 * pairs + 4);
  // Keep a small carrier's RHS panels resident across the sigma halo sweeps.
  // A larger carrier still runs with one RHS and eviction; this never acts
  // as an admission test for the physical representation.
  const std::size_t resident_rhs = cache_bytes / sizeof(double) /
      static_cast<std::size_t>(n_unique_alpha_) / n_unique_beta_;
  const int batch_width = static_cast<int>(std::min<std::size_t>(vectors.cols(),
      std::max<std::size_t>(1, std::min(resident_rhs, cache_bytes / sigma_vector_bytes))));
  const std::size_t row_bytes = sizeof(double) *
      (static_cast<std::size_t>(n_unique_alpha_) * batch_width + n_unique_beta_);
  // Each worker gets at least a cache-line span in the right-transform
  // panel. A fixed 32-row panel gave one scalar row per worker at 32 cores.
  const std::size_t parallel_extent = static_cast<std::size_t>(
      effective_openmp_thread_count()) * (64 / sizeof(double));
  const int extent = static_cast<int>(std::min<std::size_t>(
      n_unique_beta_, std::max<std::size_t>(1,
          std::min(parallel_extent, cache_bytes / std::max<std::size_t>(1, 4 * row_bytes)))));
  const int count = (n_unique_beta_ + extent - 1) / extent;

  // Immutable accepted-point rows are shared across H/S and RHS batches.
  // Borrowed tiles remain valid even if the bounded transform cache evicts.
  const auto beta_rows = [&](int panel) {
    const int first = panel * extent;
    return beta_transform.row_block(first, std::min(extent, n_unique_beta_ - first));
  };

  const auto* support=product_support_.get();

  for (int vector_begin = 0; vector_begin < vectors.cols(); vector_begin += batch_width) {
    ++result.tiles.expansions;
    const int width = std::min<int>(batch_width, vectors.cols() - vector_begin);
    // Expand once per RHS, not once per beta panel. For a covalent space
    // these groups are exactly (A, complement(A)): only O(U) coefficients.
    Eigen::MatrixXd coefficients = Eigen::MatrixXd::Zero(spin_products_.size(), width);
#pragma omp parallel for schedule(static) if(spin_products_.size() > 32)
    for (int group = 0; group < static_cast<int>(spin_products_.size()); ++group)
      for (int vector = 0; vector < width; ++vector)
        for (std::size_t term = spin_term_offsets_[group];
             term < spin_term_offsets_[group + 1]; ++term)
          coefficients(group, vector) += spin_terms_[term].coefficient *
              vectors(spin_terms_[term].structure, vector_begin + vector);
    Eigen::MatrixXd h_products = Eigen::MatrixXd::Zero(spin_products_.size(), width);
    Eigen::MatrixXd m_products = Eigen::MatrixXd::Zero(spin_products_.size(), width);
    const auto project = [&](int panel, const Eigen::MatrixXd& rows,
                             Eigen::MatrixXd& image, Eigen::MatrixXd& target) {
      // C^T (T_alpha^T image T_beta): only touched original spin products
      // are sampled. This is exact for any selected structure subspace.
      if (support)
        alpha_transform.apply_supported_left(&image, *support,
            panel * extent, rows.rows(), true);
      else alpha_transform.apply_adjoint_left(&image);
#pragma omp parallel for schedule(static) if(spin_products_.size() > 32)
      for (int group = 0; group < static_cast<int>(spin_products_.size()); ++group) {
        const int alpha = spin_products_[group] / n_unique_beta_;
        const int beta = spin_products_[group] % n_unique_beta_;
        for (int vector = 0; vector < width; ++vector)
          target(group, vector) += image.row(alpha)
              .segment(vector * rows.rows(), rows.rows()).dot(rows.col(beta));
      }
    };
    const auto contract = [&]() {
      ++result.tiles.structure_projections;
      // Apply P^T only after all panels have contributed to each spin pair.
#pragma omp parallel for schedule(static) if(n_structures_ > 32)
      for (int structure = 0; structure < n_structures_; ++structure) {
        for (int vector = 0; vector < width; ++vector) {
          double h = 0.0, m = 0.0;
          for (std::size_t offset = structure_term_offsets_[structure];
               offset < structure_term_offsets_[structure + 1]; ++offset) {
            const auto& term = spin_terms_[structure_term_indices_[offset]];
            const auto group = std::lower_bound(spin_products_.begin(),
                spin_products_.end(), term.spin_product) - spin_products_.begin();
            h += term.coefficient * h_products(group, vector);
            m += term.coefficient * m_products(group, vector);
          }
          result.hamiltonian(structure, vector_begin + vector) = h;
          result.overlap(structure, vector_begin + vector) = m;
        }
      }
    };

    std::vector<bool> overlap_consumed(count, false);
    ColumnPanels input(cache_bytes, [&](int panel) -> Eigen::MatrixXd {
      ++result.tiles.input_panels;
      const auto row_tile = beta_rows(panel);
      const Eigen::MatrixXd& rows = *row_tile;
      Eigen::MatrixXd values = Eigen::MatrixXd::Zero(n_unique_alpha_, width * rows.rows());
      // Each alpha owns its row. CSR maps structures straight to strings;
      // no determinant vector or full raw spin-product grid is formed.
#pragma omp parallel for schedule(static) if(n_unique_alpha_ > 32)
      for (int alpha = 0; alpha < n_unique_alpha_; ++alpha) {
        const int first = alpha * n_unique_beta_;
        auto group = std::lower_bound(spin_products_.begin(), spin_products_.end(), first);
        for (; group != spin_products_.end() && *group < first + n_unique_beta_; ++group) {
          const auto id = static_cast<std::size_t>(group - spin_products_.begin());
          for (int vector = 0; vector < width; ++vector) {
            values.row(alpha).segment(vector * rows.rows(), rows.rows()).noalias() +=
                coefficients(id, vector) * rows.col(*group - first).transpose();
          }
        }
      }
      if (support)
        alpha_transform.apply_supported_left(&values, *support,
            panel * extent, rows.rows());
      else alpha_transform.apply_left(&values);
      if (!overlap_consumed[panel]) {
        ++result.tiles.overlap_projections;
        Eigen::MatrixXd overlap = values;
        project(panel, rows, overlap, m_products);
        overlap_consumed[panel] = true;
      }
      return values;
    });
    if (overlap_only) {
      for (int panel = 0; panel < count; ++panel) input.get(panel);
      result.tiles.input_evictions += input.evictions;
      result.tiles.input_resident_bytes = std::max(result.tiles.input_resident_bytes, input.peak_bytes);
      contract();
      continue;
    }
    // Prepare resident panels while outside the persistent sigma team, so
    // exterior transforms and overlap projection retain their own parallelism.
    // Larger carriers simply generate panels on demand; their space is unchanged.
    if (resident_rhs >= static_cast<std::size_t>(width))
      for (int panel = 0; panel < count; ++panel) input.get(panel);
    ColumnPanels output(cache_bytes,
        [&](int panel) -> Eigen::MatrixXd {
          return Eigen::MatrixXd::Zero(n_unique_alpha_,
              width * std::min(extent, n_unique_beta_ - panel * extent));
        },
        [&](int panel, Eigen::MatrixXd& image) {
          ++result.tiles.hamiltonian_projections;
          const auto rows = beta_rows(panel);
          project(panel, *rows, image, h_products);
        },
        [&](int panel) {
          // F^T = F^T P_S: only exact structural support S contributes to
          // the selected VB image. I\S is still kept in the sigma channels.
          // This byte mask belongs to one bounded output panel, not a global
          // product-space table, and is shared by every RHS in the panel.
          std::vector<std::uint8_t> allowed;
          if (!support) return allowed;
          const int first = panel * extent;
          const int columns = std::min(extent, n_unique_beta_ - first);
          allowed.resize(static_cast<std::size_t>(n_unique_alpha_) * columns, 0);
          int team = 1;
#ifdef _OPENMP
          team = omp_get_num_threads();
#endif
#pragma omp taskloop num_tasks(team) shared(allowed)
          for (int column = 0; column < columns; ++column) {
            auto paths = support->column(first + column, 0);
            int a, b;
            bool label;
            auto* mask = allowed.data() + static_cast<std::size_t>(column) * n_unique_alpha_;
            while (paths.next(a, b, label)) mask[a] = 1;
          }
          return allowed;
        });
    sigma.apply_tiled(width,
        [&](const std::vector<int>& betas, const std::vector<int>& alpha_rows,
            Eigen::Ref<Eigen::MatrixXd> columns) {
          int team = 1;
#ifdef _OPENMP
          team = omp_get_num_threads();
#endif
          for (int first = 0; first < static_cast<int>(betas.size());) {
            const int panel_id = betas[first] / extent;
            int last = first + 1;
            while (last < static_cast<int>(betas.size()) && betas[last] / extent == panel_id) ++last;
            const Eigen::MatrixXd& panel = input.get(panel_id);
            const int panel_width = static_cast<int>(panel.cols()) / width;
            const double* source = panel.data();
            double* target = columns.data();
            const Eigen::Index source_stride = panel.outerStride(), target_stride = columns.outerStride();
            // Implicit taskgroup completion pins this panel until all reads
            // finish; only then may the next get() mutate or evict the LRU.
#pragma omp taskloop num_tasks(team) shared(betas, alpha_rows)
            for (int packed = 0; packed < (last - first) * width; ++packed) {
              const int j = first + packed / width, v = packed % width;
              const double* from = source + (v * panel_width + betas[j] % extent) * source_stride;
              double* to = target + (j * width + v) * target_stride;
              for (int row = 0; row < static_cast<int>(alpha_rows.size()); ++row)
                to[row] = from[alpha_rows[row]];
            }
            first = last;
          }
        },
        [&](const std::vector<int>& betas, const std::vector<int>& alpha_rows,
            const mcscf::DirectCiColumnAccumulator& accumulate) {
          int team = 1;
#ifdef _OPENMP
          team = omp_get_num_threads();
#endif
          for (int first = 0; first < static_cast<int>(betas.size());) {
            const int panel_id = betas[first] / extent;
            int last = first + 1;
            while (last < static_cast<int>(betas.size()) && betas[last] / extent == panel_id) ++last;
            Eigen::MatrixXd& panel = output.get(panel_id);
            const auto& allowed = output.mask(panel_id);
            const auto* mask = allowed.empty() ? nullptr : allowed.data();
            const int panel_width = static_cast<int>(panel.cols()) / width;
            double* target = panel.data();
            const Eigen::Index target_stride = panel.outerStride();
            // Split alpha ranges as well as columns: the alpha-adjoint work
            // is concentrated in target beta columns. Tasks add directly to
            // disjoint destination ranges, without an image slab or atomics.
            const int chunks = std::min<int>(team, alpha_rows.size());
            const int column_count = (last - first) * width;
#pragma omp taskloop num_tasks(team) shared(betas, alpha_rows, accumulate)
            for (int packed = 0; packed < column_count * chunks; ++packed) {
              // Every task sees all beta columns over its alpha stripe.
              // Column-major task order would concentrate the expensive
              // alpha adjoint (target columns only) on a few workers.
              const int column = packed % column_count;
              const int chunk = packed / column_count;
              const int j = first + column / width, v = column % width;
              double* to = target + (v * panel_width + betas[j] % extent) * target_stride;
              const int begin = alpha_rows.size() * chunk / chunks;
              const int end = alpha_rows.size() * (chunk + 1) / chunks;
              const auto* column_mask = mask ? mask +
                  static_cast<std::size_t>(betas[j] % extent) * n_unique_alpha_ : nullptr;
              accumulate(j, v, begin, end, to, column_mask);
            }
            first = last;
          }
        }, support, spin_exchange_parity_);
    output.finish();
    // Exchange quotienting may never request some input columns. The overlap
    // is still the complete Gram action, independent of the Hamiltonian sweep.
    for (int panel = 0; panel < count; ++panel)
      if (!overlap_consumed[panel]) input.get(panel);
    result.tiles.input_evictions += input.evictions;
    result.tiles.output_evictions += output.evictions;
    result.tiles.input_resident_bytes = std::max(result.tiles.input_resident_bytes, input.peak_bytes);
    result.tiles.output_resident_bytes = std::max(result.tiles.output_resident_bytes, output.peak_bytes);
    contract();
  }
  return result;
}

}  // namespace xmvb::vb
