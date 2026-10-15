#pragma once

#include "mcscf/ci/factors.hpp"

namespace xmvb::mcscf::detail {

void validate_source(const DirectCiSigmaAction& sigma, const SpinFactorSource& source);
void add_product(const Eigen::Ref<const Eigen::MatrixXd>& left,
    const Eigen::Ref<const Eigen::MatrixXd>& right, Eigen::Ref<Eigen::MatrixXd> result);
void transpose_product(const Eigen::Ref<const Eigen::MatrixXd>& left,
    const Eigen::Ref<const Eigen::MatrixXd>& right, Eigen::Ref<Eigen::MatrixXd> result);

} // namespace xmvb::mcscf::detail
