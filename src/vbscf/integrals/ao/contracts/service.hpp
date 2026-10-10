#pragma once

#include <vector>
#include <Eigen/Core>
#include "vbscf/integrals/ao/contracts/input.hpp"
#include "vbscf/integrals/ao/libcint/input.hpp"

namespace xmvb::vb {

using PairMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

/** @brief Integral-only worker transport; never owns optimizer/structure state. */
class IntegralService {
public:
  virtual ~IntegralService() = default;
  virtual AoIntegralInput build(const LibcintInput&) = 0;
  virtual bool dispatching() const noexcept = 0;
  virtual std::vector<double> h1e(const AoIntegralInput&, const double*, bool transpose) = 0;
  virtual void fused(const AoIntegralInput&, const double*, const double*,
                     std::vector<double>*, std::vector<double>*) = 0;
  virtual void pair(const AoPairGraph&, int n_bf, const Eigen::Ref<const PairMatrix>&,
                    PairMatrix*) = 0;
  /** @brief Set by the MPI session only; null in ordinary single-process runs. */
  inline static IntegralService* current = nullptr;
};

}  // namespace xmvb::vb
