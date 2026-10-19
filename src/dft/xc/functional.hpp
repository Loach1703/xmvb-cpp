#pragma once

#include <Eigen/Core>
#include <memory>
#include <string>

namespace xmvb::dft {

/** Unpolarized LDA components named by Libxc, separated by commas.
 * `exc` is energy per electron; `vrho` is d(rho * exc)/d(rho).
 * Libxc's own per-functional low-density thresholds remain in effect.
 */
class LdaFunctional {
 public:
  explicit LdaFunctional(const std::string& name = "LDA_X,LDA_C_VWN");
  ~LdaFunctional();
  LdaFunctional(const LdaFunctional&) = delete;
  LdaFunctional& operator=(const LdaFunctional&) = delete;
  LdaFunctional(LdaFunctional&&) noexcept;
  LdaFunctional& operator=(LdaFunctional&&) noexcept;

  static std::string library_version();

  void evaluate(const Eigen::Ref<const Eigen::VectorXd>& density,
                Eigen::VectorXd* exc, Eigen::VectorXd* vrho) const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace xmvb::dft
