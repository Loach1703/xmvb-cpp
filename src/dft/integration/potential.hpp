#pragma once

#include <Eigen/Core>
#include <memory>

namespace xmvb::integrals {
struct LibcintInput;
class BasisEvaluator;
}  // namespace xmvb::integrals

namespace xmvb::dft {
struct Grid;
class LdaFunctional;

struct XcResult {
  double energy = 0.0;
  double electron_count = 0.0;
  Eigen::MatrixXd potential;
};

/** Integrates an unpolarized LDA using the total restricted AO density.
 * Input, grid, and functional must outlive the integrator and remain unmodified.
 */
class LdaIntegrator {
 public:
  LdaIntegrator(const integrals::LibcintInput& input,
                const Grid& grid,
                const LdaFunctional& functional);
  ~LdaIntegrator();
  LdaIntegrator(const LdaIntegrator&) = delete;
  LdaIntegrator& operator=(const LdaIntegrator&) = delete;

  XcResult evaluate(
      const Eigen::Ref<const Eigen::MatrixXd>& density) const;

 private:
  const Grid& grid_;
  const LdaFunctional& functional_;
  std::unique_ptr<integrals::BasisEvaluator> basis_;
};

}  // namespace xmvb::dft
