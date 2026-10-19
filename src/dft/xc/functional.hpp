#pragma once

#include <Eigen/Core>

#include <memory>
#include <string>

namespace xmvb::dft {

struct XcValues {
  Eigen::VectorXd exc;     // Energy per electron.
  Eigen::MatrixXd vrho;    // Derivative of rho * exc with respect to spin density.
  Eigen::MatrixXd vsigma;  // Derivative with respect to gradient dot products.
};

/** A sum of pure Libxc LDA and GGA exchange/correlation components.
 * rho is N x n_spin: total density for n_spin=1, alpha/beta for n_spin=2.
 * For GGA, sigma is N x 1 or N x 3 (alpha-alpha, alpha-beta, beta-beta).
 * For LDA, sigma may have zero or the corresponding GGA number of columns;
 * its values are ignored. Libxc's own low-density thresholds remain in effect.
 *
 * evaluate allocates, in addition to the XcValues results, bulk Libxc buffers
 * of N * (1 + 2*n_spin) doubles for LDA or
 * N * (1 + 2*n_spin + 2*n_sigma) doubles for GGA. It allocates no storage
 * per point or per component. n_sigma is 1 or 3 according to n_spin.
 */
class Functional {
 public:
  explicit Functional(const std::string& name = "LDA_X,LDA_C_VWN",
                      int n_spin = 1);
  ~Functional();
  Functional(const Functional&) = delete;
  Functional& operator=(const Functional&) = delete;
  Functional(Functional&&) noexcept;
  Functional& operator=(Functional&&) noexcept;

  static std::string library_version();
  int n_spin() const;
  bool is_gga() const;

  void evaluate(const Eigen::Ref<const Eigen::MatrixXd>& rho,
                const Eigen::Ref<const Eigen::MatrixXd>& sigma,
                XcValues* out) const;

  /** Apply the pointwise XC potential Hessian to a density/sigma direction.
   * Inputs have the same shapes and low-density handling as evaluate();
   * delta_rho and delta_sigma have the shapes of rho and sigma, respectively.
   * For LDA, sigma and delta_sigma are ignored and delta_vsigma is N x 0.
   * For GGA, delta_vrho and delta_vsigma contain the directional changes of
   * the derivatives of rho * exc with respect to rho and sigma. The direction
   * is applied at the (possibly roundoff-clamped) input; it is not clamped.
   * Inputs and outputs must not alias. Temporary Libxc storage is limited to
   * 128 points, so callers may supply point tiles without an N-scaled copy.
   */
  void response(const Eigen::Ref<const Eigen::MatrixXd>& rho,
                const Eigen::Ref<const Eigen::MatrixXd>& sigma,
                const Eigen::Ref<const Eigen::MatrixXd>& delta_rho,
                const Eigen::Ref<const Eigen::MatrixXd>& delta_sigma,
                Eigen::MatrixXd* delta_vrho,
                Eigen::MatrixXd* delta_vsigma) const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace xmvb::dft
