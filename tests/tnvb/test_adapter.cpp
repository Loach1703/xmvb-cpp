#include "tnvb/gradient/adapter.hpp"
#include <cmath>
#include <stdexcept>
#include <limits>
#include "acceptance.hpp"

int main() {
  // Preserve benzene's strict matrix gate while applying the separately
  // confirmed C14 policy; nonfinite diagnostics can never pass either profile.
  const tnvb_test::ComparisonErrors measured{3.23e-11,3.18e-9,1.524e-10,7.05e-10,9.11e-10,1.75e-9};
  if (tnvb_test::acceptance_profile(6).passed(measured) ||
      !tnvb_test::acceptance_profile(14).passed(measured) ||
      tnvb_test::AcceptanceProfile::strict_matrices_passed(measured)) {
    throw std::runtime_error("benzene/C14 acceptance scopes were mixed");
  }
  auto invalid=measured;
  invalid.s=std::numeric_limits<double>::quiet_NaN();
  if (tnvb_test::acceptance_profile(14).passed(invalid)) {
    throw std::runtime_error("nonfinite C14 matrix diagnostic bypassed validation");
  }
  // An orthonormal two-electron closed-shell determinant on orbital zero.
  // E=2h00+g0000. At fixed h/g its norm derivative is -2h00-2g0000.
  xmvb::tnvb::WorkerResponse response;
  response.energies=Eigen::VectorXd::Constant(1,-1.7);
  response.density_one=Eigen::MatrixXd::Zero(2,2);
  response.density_one(0,0)=2.0;
  response.density_two_buffer.assign(16,0.0);
  response.density_two_buffer[0]=2.0;
  response.contracted_three=Eigen::MatrixXd::Zero(2,2);
  xmvb::vb::ActiveSpaceGradientResult result;
  result.active_space_one_electron_result.h1e_act=Eigen::MatrixXd::Zero(2,2);
  result.active_space_one_electron_result.h1e_act(0,0)=-1.0;
  xmvb::tnvb::populate_active_space_adjoints(response,&result);
  if (std::abs(result.active_orbital_overlap_gradient[0]-1.4)>1e-14 ||
      result.active_one_electron_gradient[0]!=2.0 ||
      result.packed_active_two_electron_gradient[0]!=1.0) {
    throw std::runtime_error("two-electron TNVB adjoint convention failed");
  }
  return 0;
}
