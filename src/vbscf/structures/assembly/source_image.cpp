#include "vbscf/structures/assembly/action.hpp"

#include <algorithm>
#include <stdexcept>
#include <utility>

#include "parallel/openmp.hpp"
#include "mcscf/ci/factors_internal.hpp"

namespace xmvb::vb {

StructureAction::SourceImageAdjoint StructureAction::project_source_adjoint(
    const mcscf::SpinFactorSource& source, const mcscf::SpinFactorSource& right,
    const ExteriorOrbitalTransform& alpha, const ExteriorOrbitalTransform& beta,
    const mcscf::DirectCiSigmaAction& sigma, const Eigen::MatrixXd* output_generator,
    bool project_output, double residual_energy) const {
  // Image and covector slabs coexist, so partition one budget, not two copies.
  constexpr std::size_t image_bytes = mcscf::kFactorWorkspaceBytes/2;
  mcscf::FactorActionWorkspace workspace(source,sigma.channel_kernel().rows(),
      image_bytes,false,project_output && output_generator ? 2 : 0);
  const int columns = workspace.source_tile()*std::max(1,workspace.channel_tile());
  const std::size_t rows = static_cast<std::size_t>(n_unique_alpha_)+n_unique_beta_;
  const int n = sigma.n_orbitals();
  const std::size_t reductions = static_cast<std::size_t>(effective_openmp_thread_count()+1)*n*n;
  const std::size_t capacity = (mcscf::kFactorWorkspaceBytes-image_bytes)/sizeof(double);
  const std::size_t per_column = 2*rows+columns+right.workspace_rows;
  if (capacity < reductions+per_column)
    throw std::length_error("source-image adjoint cannot hold one right-factor panel");
  const int tile = std::min<std::size_t>(right.count,(capacity-reductions)/per_column);
  SourceImageAdjoint result{{Eigen::MatrixXd::Zero(project_output ? n_structures_ : 0,1),
                            Eigen::MatrixXd::Zero(project_output ? n_structures_ : 0,1)},Eigen::MatrixXd::Zero(n,n)};
  if (!tile) return result;
  Eigen::MatrixXd a(n_unique_alpha_,tile), b(n_unique_beta_,tile);
  Eigen::MatrixXd za(n_unique_alpha_,tile), zb(n_unique_beta_,tile), moments(columns,tile);
  for (int first = 0; first < right.count; first += tile) {
    const int count = std::min(tile,right.count-first);
    auto ra = a.leftCols(count), rb = b.leftCols(count);
    auto wa = za.leftCols(count), wb = zb.leftCols(count);
    right.load(first,count,ra,rb); wa.setZero(); wb.setZero();
    const auto consume = [&](const Eigen::Ref<const Eigen::MatrixXd>& x,
        const Eigen::Ref<const Eigen::MatrixXd>& y, double scale) {
      auto moment = moments.topLeftCorner(x.cols(),count);
      const auto add = [&](const Eigen::Ref<const Eigen::MatrixXd>& left,
          const Eigen::Ref<const Eigen::MatrixXd>& spectator, double weight) {
        mcscf::detail::transpose_product(spectator,rb,moment); moment *= weight;
        mcscf::detail::add_product(left,moment,wa);
        mcscf::detail::transpose_product(left,ra,moment); moment *= weight;
        mcscf::detail::add_product(spectator,moment,wb);
      };
      // Only the left H image needs physical parity reconstruction. Generators
      // commute with common spin exchange; no selected-space closure is used.
      add(x,y,scale*(spin_exchange_parity_ ? 0.5 : 1));
      if (spin_exchange_parity_) add(y,x,scale*0.5*spin_exchange_parity_);
    };
    // A smaller RHS residency can replay source panels, but projection is done
    // only on the first pass. It is never repeated for the extra adjoint tiles.
    auto images = project_factor_source(source,workspace,alpha,beta,sigma,false,
        output_generator,consume,project_output && first == 0,residual_energy);
    if (first == 0) result.images = std::move(images);
    else {
      result.images.tiles.source_tiles += images.tiles.source_tiles;
      result.images.tiles.source_columns += images.tiles.source_columns;
      result.images.tiles.channel_blocks += images.tiles.channel_blocks;
    }
    result.generator += sigma.spin_generator_adjoint(false,wa,ra);
    result.generator += sigma.spin_generator_adjoint(true,wb,rb);
  }
  result.images.tiles.source_workspace_bound += sizeof(double)*
      (per_column*tile+reductions);
  return result;
}

} // namespace xmvb::vb
