#include "vbscf/structures/assembly/action.hpp"

#include <algorithm>

#include "core/openmp.hpp"
#include "mcscf/ci/product_support.hpp"
#include "mcscf/ci/sigma.hpp"
#include "vbscf/structures/assembly/column_tiles.hpp"
#include "vbscf/structures/orthogonal_ci/exterior_transform.hpp"

namespace xmvb::vb {
using detail::ColumnPanels;
using detail::ColumnPattern;
using detail::ColumnTile;

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
  if (source_contraction_)
    return apply_source_tiles(vectors, alpha_transform, beta_transform, sigma,
        overlap_only, 64ULL * 1024ULL * 1024ULL);
  return contract_product_tiles(vectors,alpha_transform,beta_transform,sigma,
      overlap_only,workspace_bytes).images;
}

StructureResponseSweep StructureAction::contract_product_tiles(
    const Eigen::Ref<const Eigen::MatrixXd>& vectors,
    const ExteriorOrbitalTransform& alpha_transform,const ExteriorOrbitalTransform& beta_transform,
    const mcscf::DirectCiSigmaAction& sigma,bool overlap_only,std::size_t workspace_bytes,
    const StructureAdjointState* accepted,bool collect_images) const {
  const int n=sigma.n_orbitals(),d=sigma.channel_kernel().rows();
  StructureResponseSweep sweep{{Eigen::MatrixXd::Zero(collect_images ? n_structures_:0,vectors.cols()),
      Eigen::MatrixXd::Zero(collect_images ? n_structures_:0,vectors.cols())},
      Eigen::MatrixXd::Zero(accepted ? n:0,accepted ? n:0),
      Eigen::MatrixXd::Zero(accepted ? d:0,accepted ? d:0),
      Eigen::MatrixXd::Zero(accepted ? n:0,accepted ? n:0)};
  auto& result=sweep.images;
  if (accepted && (vectors.array()==0.0).all()) return sweep;
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
      cooperative_openmp_thread_count()) * (64 / sizeof(double));
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
  std::vector<std::uint8_t> panel_support(count,0);
  std::vector<int> reachable_panels;
  const auto include_panel=[&](int panel) {
    if (!panel_support[panel]) { panel_support[panel]=1; reachable_panels.push_back(panel); }
  };
  if (support && !support->full()) {
    auto columns=support->cursor(0);
    int beta;
    while (columns.next_column(beta)) include_panel(beta/extent);
    std::sort(reachable_panels.begin(),reachable_panels.end());
  } else for (int panel=0;panel<count;++panel) include_panel(panel);
  // Weak references share input/output patterns without an unbounded support
  // table. Their strong lifetime and byte charge belong to the bounded LRUs.
  std::vector<std::weak_ptr<const ColumnPattern>> patterns(count);
  const auto pattern_at=[&](int panel)->std::shared_ptr<const ColumnPattern> {
    if (auto pattern=patterns[panel].lock()) return pattern;
    auto pattern=std::make_shared<ColumnPattern>();
    const int cols=std::min(extent,n_unique_beta_-panel*extent);
    pattern->offsets.resize(cols+1);
    pattern->offsets[0]=0;
    pattern->full=!support || support->full();
    if (pattern->full)
      for (int b=0;b<cols;++b) pattern->offsets[b+1]=(b+1)*n_unique_alpha_;
    else {
      std::vector<std::vector<int>> rows(cols);
      const int workers=bounded_openmp_thread_count(cooperative_openmp_thread_count(),cols);
      cooperative_openmp_for(cols,workers,[&](std::ptrdiff_t b) {
        auto cursor=support->column(panel*extent+b,0);
        int a,beta; bool label;
        while (cursor.next(a,beta,label)) rows[b].push_back(a);
        if (!std::is_sorted(rows[b].begin(),rows[b].end())) std::sort(rows[b].begin(),rows[b].end());
      });
      for (int b=0;b<cols;++b) pattern->offsets[b+1]=pattern->offsets[b]+rows[b].size();
      pattern->rows.resize(pattern->offsets[cols]);
      for (int b=0;b<cols;++b)
        std::copy(rows[b].begin(),rows[b].end(),pattern->rows.data()+pattern->offsets[b]);
    }
    patterns[panel]=pattern;
    return pattern;
  };
  std::vector<int> alpha_position(n_unique_alpha_,-1);

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
    Eigen::MatrixXd h_products = Eigen::MatrixXd::Zero(collect_images ? spin_products_.size():0, width);
    Eigen::MatrixXd m_products = Eigen::MatrixXd::Zero(collect_images ? spin_products_.size():0, width);
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
      if (!collect_images) return;
      ++result.tiles.structure_projections;
      // Apply P^T only after all panels have contributed to each spin pair.
#pragma omp parallel for schedule(static) if(n_structures_ > 32)
      for (int structure = 0; structure < n_structures_; ++structure) {
        for (int vector = 0; vector < width; ++vector) {
          double h = 0.0, m = 0.0;
          for (std::size_t offset = structure_term_offsets_[structure];
               offset < structure_term_offsets_[structure + 1]; ++offset) {
            const auto& term = spin_terms_[structure_term_indices_[offset]];
            h += term.coefficient * h_products(term.group, vector);
            m += term.coefficient * m_products(term.group, vector);
          }
          result.hamiltonian(structure, vector_begin + vector) = h;
          result.overlap(structure, vector_begin + vector) = m;
        }
      }
    };

    std::vector<bool> overlap_consumed(count, false);
    ColumnPanels input(cache_bytes, [&](int panel) -> ColumnTile {
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
        if (accepted)
          for (int v=0;v<width;++v) {
            const int state=vector_begin+v;
            const double weight=accepted->weights[state];
            if (weight==0) continue;
            // D(r,q)=D(q,r)^T. Both right carriers are accepted data, so
            // the generated q panel needs neither neighbors nor another H.
            sweep.generator.noalias()+=2*weight*sigma.generator_adjoint_columns(
                panel*extent,rows.rows(),values.middleCols(v*rows.rows(),rows.rows()),
                {{accepted->orthogonal_coefficients[state],-accepted->energies[state],false},
                 {accepted->residuals[state],1.0,true}},true);
          }
        if (collect_images) {
          Eigen::MatrixXd overlap = values;
          project(panel, rows, overlap, m_products);
        }
        overlap_consumed[panel] = true;
      }
      return detail::pack_columns(std::move(values),pattern_at(panel));
    });
    if (overlap_only) {
      for (const int panel:reachable_panels) input.get(panel);
      result.tiles.input_evictions += input.evictions;
      result.tiles.input_resident_bytes = std::max(result.tiles.input_resident_bytes, input.peak_bytes);
      contract();
      continue;
    }
    // Prepare resident panels while outside the persistent sigma team, so
    // exterior transforms and overlap projection retain their own parallelism.
    // Larger carriers simply generate panels on demand; their space is unchanged.
    if (resident_rhs >= static_cast<std::size_t>(width))
      for (const int panel:reachable_panels) input.get(panel);
    ColumnPanels output(cache_bytes,
        [&](int panel) -> ColumnTile {
          auto pattern=pattern_at(panel);
          return {pattern,Eigen::MatrixXd::Zero(pattern->offsets[pattern->columns()],width)};
        },
        [&](int panel,ColumnTile& tile) {
          ++result.tiles.hamiltonian_projections;
          Eigen::MatrixXd image=tile.unpack(n_unique_alpha_);
          const int cols=tile.pattern->columns();
          if (accepted)
            for (int v=0;v<width;++v) {
              const int state=vector_begin+v;
              if (accepted->weights[state]==0) continue;
              // E_pq Y stays in S for p<=q. The S-masked image therefore
              // supplies the exact upper Cholesky cotangent, including when
              // sigma is an exchange representative of a common parity.
              sweep.generator.noalias()+=2*accepted->weights[state]*sigma.generator_adjoint_columns(
                  panel*extent,cols,image.middleCols(v*cols,cols),
                  {{accepted->orthogonal_coefficients[state],1.0,false}},true);
            }
          if (collect_images) {
            const auto rows=beta_rows(panel);
            project(panel,*rows,image,h_products);
          }
        });
    const mcscf::DirectCiTileReader read =
        [&](const std::vector<int>& betas, const std::vector<int>& alpha_rows,
            Eigen::Ref<Eigen::MatrixXd> columns) {
          int team = 1;
#ifdef _OPENMP
          team = omp_get_num_threads();
#endif
          for (int i=0;i<static_cast<int>(alpha_rows.size());++i) alpha_position[alpha_rows[i]]=i;
          for (int first = 0; first < static_cast<int>(betas.size());) {
            const int panel_id = betas[first] / extent;
            int last = first + 1;
            while (last < static_cast<int>(betas.size()) && betas[last] / extent == panel_id) ++last;
            if (!panel_support[panel_id]) {
              columns.middleCols(first*width,(last-first)*width).setZero();
              first=last;
              continue;
            }
            const ColumnTile& panel = input.get(panel_id);
            const auto* offsets=panel.pattern->offsets.data();
            const auto* indices=panel.pattern->rows.data();
            const bool full=panel.pattern->full;
            const double* source = panel.values.data();
            double* target = columns.data();
            const Eigen::Index source_stride = panel.values.outerStride(), target_stride = columns.outerStride();
            // Implicit taskgroup completion pins this panel until all reads
            // finish; only then may the next get() mutate or evict the LRU.
#pragma omp taskloop num_tasks(team) shared(betas, alpha_rows)
            for (int packed = 0; packed < (last - first) * width; ++packed) {
              const int j = first + packed / width, v = packed % width;
              const int beta=betas[j]%extent;
              const double* from=source+v*source_stride;
              double* to = target + (j * width + v) * target_stride;
              if (full)
                for (int row=0;row<static_cast<int>(alpha_rows.size());++row)
                  to[row]=from[offsets[beta]+alpha_rows[row]];
              else {
                std::fill_n(to,alpha_rows.size(),0.0);
                for (int i=offsets[beta];i<offsets[beta+1];++i) {
                  const int row=alpha_position[indices[i]];
                  if (row>=0) to[row]=from[i];
                }
              }
            }
            first = last;
          }
          for (const int a:alpha_rows) alpha_position[a]=-1;
        };
    const mcscf::DirectCiTileWriter add =
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
            if (!panel_support[panel_id]) { first=last; continue; }
            ColumnTile& panel = output.get(panel_id);
            const auto* offsets=panel.pattern->offsets.data();
            const auto* indices=panel.pattern->rows.data();
            const bool full=panel.pattern->full;
            double* target = panel.values.data();
            const Eigen::Index target_stride = panel.values.outerStride();
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
              const int beta=betas[j]%extent,offset=offsets[beta],length=offsets[beta+1]-offset;
              if (length==0) continue;
              double* to=target+v*target_stride+offset;
              const int begin = alpha_rows.size() * chunk / chunks;
              const int end = alpha_rows.size() * (chunk + 1) / chunks;
              accumulate(j,v,begin,end,to,full ? nullptr:indices+offset,full ? 0:length);
            }
            first = last;
          }
        };
    if (accepted) {
      const mcscf::DirectCiTileReader right = [&](const std::vector<int>& betas,
          const std::vector<int>& rows,Eigen::Ref<Eigen::MatrixXd> values) {
        for (int j=0;j<static_cast<int>(betas.size());++j)
          for (int v=0;v<width;++v) {
            const int state=vector_begin+v;
            const double weight=accepted->weights[state];
            if (weight==0) values.col(j*width+v).setZero();
            else
              for (int a=0;a<static_cast<int>(rows.size());++a)
                values(a,j*width+v)=weight*accepted->orthogonal_coefficients[state](rows[a],betas[j]);
          }
      };
      const auto adjoint=sigma.apply_tiled_with_adjoint(width,read,add,right,support,spin_exchange_parity_);
      ++result.tiles.integral_adjoint_sweeps;
      sweep.one_electron.noalias()+=2*adjoint.one_electron;
      sweep.pair_kernel.noalias()+=2*adjoint.pair_kernel;
    } else sigma.apply_tiled(width,read,add,support,spin_exchange_parity_);
    output.finish();
    // Exchange quotienting may never request some input columns. The overlap
    // is still the complete Gram action, independent of the Hamiltonian sweep.
    for (const int panel:reachable_panels)
      if (!overlap_consumed[panel]) input.get(panel);
    result.tiles.input_evictions += input.evictions;
    result.tiles.output_evictions += output.evictions;
    result.tiles.input_resident_bytes = std::max(result.tiles.input_resident_bytes, input.peak_bytes);
    result.tiles.output_resident_bytes = std::max(result.tiles.output_resident_bytes, output.peak_bytes);
    contract();
  }
  return sweep;
}

}  // namespace xmvb::vb
