#pragma once

#include "mcscf/ci/factors.hpp"

namespace xmvb::mcscf::detail {

enum class FactorContraction { Integral, Generator, Residual };

/** @brief Joint source/channel schedule; counts work, never changes the operator. */
struct FactorPanels {
  int left = 0, right = 0, channels = 0;
  std::size_t workspace_bytes = 0;
  long double work = 0, source_columns = 0;
};

FactorPanels plan_factor_panels(const DirectCiSigmaAction& sigma,
    const SpinFactorSource& left, const SpinFactorSource& right,
    std::size_t workspace_bytes, FactorContraction contraction);

void validate_source(const DirectCiSigmaAction& sigma, const SpinFactorSource& source);
void add_product(const Eigen::Ref<const Eigen::MatrixXd>& left,
    const Eigen::Ref<const Eigen::MatrixXd>& right, Eigen::Ref<Eigen::MatrixXd> result);
void transpose_product(const Eigen::Ref<const Eigen::MatrixXd>& left,
    const Eigen::Ref<const Eigen::MatrixXd>& right, Eigen::Ref<Eigen::MatrixXd> result);

} // namespace xmvb::mcscf::detail
