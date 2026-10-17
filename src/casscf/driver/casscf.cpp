#include "casscf/driver/casscf.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include "casscf/derivatives/coupled.hpp"
#include "mcscf/ci/strings.hpp"
#include "mcscf/optimization/problem.hpp"
#include "mcscf/optimization/trace.hpp"

namespace xmvb::casscf {
namespace {
using Clock = std::chrono::steady_clock;

Eigen::MatrixXd initial_ci(int nactive, int electrons, int spin,
    std::size_t workspace_bytes) {
  if (nactive<1 || nactive>63 || electrons<1 || electrons>2*nactive ||
      spin<0 || spin>std::min(electrons,2*nactive-electrons) || (electrons-spin)%2)
    throw std::invalid_argument("invalid complete-CAS electron/spin space");
  const auto na=mcscf::string_count(nactive,(electrons+spin)/2);
  const auto nb=mcscf::string_count(nactive,(electrons-spin)/2);
  if (!na || !nb || na>static_cast<std::size_t>(std::numeric_limits<int>::max())/nb ||
      na>workspace_bytes/sizeof(double)/nb)
    throw std::length_error("initial CAS coefficient carrier exceeds workspace");
  // Nested alpha/beta occupations make the first determinant a pure highest-Ms
  // state. Davidson subsequently finds the lowest root in this total-spin sector.
  Eigen::MatrixXd c=Eigen::MatrixXd::Zero(na,nb);
  c(0,0)=1;
  return c;
}
}

CasscfResult solve_casscf(const integrals::AoIntegralInput& ao,
    const Eigen::MatrixXd& orbitals, int ncore, int nactive, int electrons,
    int spin, double nuclear, const CasscfOptions& options) {
  if (!(options.energy_tolerance>0) || !(options.gradient_tolerance>0) ||
      !std::isfinite(options.energy_tolerance) || !std::isfinite(options.gradient_tolerance) ||
      options.max_iterations<1 || options.ci.n_roots!=1)
    throw std::invalid_argument("CASSCF requires positive tolerances and one CI root");
  CasscfResult result;
  const auto prepare=[&](const Eigen::MatrixXd& u,const Eigen::MatrixXd& c) {
    auto point=std::make_unique<FixedCiPoint>(ao,u,ncore,nactive,nuclear,electrons,spin,c);
    const auto ci=solve_casci(point->sigma(),point->spin_space(),
        point->integrals().active.one_electron,spin,point->core_energy(),options.ci,c);
    result.ci_actions+=ci.active_ci.block_actions;
    const auto& values=ci.active_ci.eigenpairs.eigenvector_matrix;
    point->set_coefficients(Eigen::Map<const Eigen::MatrixXd>(values.data(),c.rows(),c.cols()));
    return point;
  };
  auto point=prepare(orbitals,initial_ci(nactive,electrons,spin,options.ci.workspace_bytes));
  result.initial_energy=point->energy();
  const RotationChart chart(ao.n_basis_functions,ncore,nactive);
  const auto gradient_norm=[&](const FixedCiPoint& value) {
    const double orbital=chart.pullback(value.orbitals(),value.orbital_gradient_raw(),
        Eigen::VectorXd::Zero(chart.size())).norm();
    return std::hypot(orbital,std::sqrt(2.0)*value.ci_residual().norm());
  };
  result.initial_gradient_norm=gradient_norm(*point);

  if (chart.size()==0) {
    result.gradient_norm=std::sqrt(2.0)*point->ci_residual().norm();
    result.converged=result.gradient_norm<options.gradient_tolerance;
    result.termination_reason=result.converged ? "complete_casci" : "ci_residual_tolerance";
    result.energy=point->energy();
    result.orbitals=point->orbitals();
    result.coefficients=point->coefficients();
    return result;
  }

  for (int iteration=1; iteration<=options.max_iterations; ++iteration) {
    const auto start=Clock::now();
    const int ci_before=result.ci_actions;
    const double energy=point->energy();
    CoupledModel model(*point,ao,ncore,nactive);
    const auto action=[&](const Eigen::VectorXd& direction) {
      ++result.hvp_actions;
      return model.apply(direction);
    };
    const auto norm=[](const Eigen::VectorXd& vector) { return vector.norm(); };
    const auto identity=[](const Eigen::VectorXd& vector) -> Eigen::VectorXd { return vector; };
    mcscf::QuadraticProblem problem(model.gradient(),action,identity,norm);
    mcscf::CiahOptions inner;
    inner.coordinate_norm=norm;
    inner.residual_preconditioner=[&](const Eigen::VectorXd& r,double shift) {
      return model.precondition(r,shift-1e-8);
    };
    inner.residual_direction=identity;
    mcscf::CiahWorkspace workspace(problem,-model.precondition(model.gradient(),0),inner);
    mcscf::CiahMicroCallbacks callbacks;
    callbacks.coordinate_norm=norm;
    callbacks.metric_action=identity;
    callbacks.maximum_component=[&](const Eigen::VectorXd& vector) {
      return model.maximum_component(vector);
    };
    callbacks.orbital_size=model.orbital_size();
    callbacks.orbital_coordinate_norm=norm;
    callbacks.evaluate_keyframe=[&](const Eigen::VectorXd& displacement) {
      FixedCiPoint trial(ao,model.orbitals(displacement),ncore,nactive,nuclear,
          electrons,spin,model.coefficients(displacement));
      return model.keyframe(trial,displacement);
    };
    const int hvp_before=result.hvp_actions;
    auto micro=mcscf::run_ciah_microsteps(workspace,model.gradient(),energy,
        options.gradient_tolerance,callbacks);

    // Physical macro acceptance: relax CI in the finite rotated orbitals.
    // The small AH space remains available for an origin descent/curvature
    // restart; this uses only its existing H images, never another HVP.
    Eigen::VectorXd displacement=micro.displacement;
    bool restarted=false;
    const auto restart_origin=[&]() {
      if (restarted) return false;
      restarted=true;
      const auto ray=workspace.origin_ray(model.gradient());
      if (!ray) return false;
      const double largest=model.maximum_component(ray->step);
      if (!(largest>0)) return false;
      const double scale=largest>mcscf::ciah_max_component*ray->reference_component
          ? mcscf::ciah_max_component/largest : 1/ray->reference_component;
      displacement=scale*ray->step;
      return true;
    };
    std::unique_ptr<FixedCiPoint> accepted;
    double accepted_norm=0;
    while (true) {
      auto trial=prepare(model.orbitals(displacement),model.coefficients(displacement));
      const double change=trial->energy()-energy;
      if (change<0 || std::abs(change)<options.energy_tolerance) {
        accepted_norm=gradient_norm(*trial);
        if (change<0 || accepted_norm<options.gradient_tolerance) {
          accepted=std::move(trial);
          break;
        }
      }
      if (micro.keyframes.available()) {
        displacement=micro.keyframes.take().displacement;
        continue;
      }
      if ((micro.energy_model_stopped || model.gradient().dot(displacement)>=0) &&
          restart_origin()) continue;
      displacement*=0.5;
      if (displacement.norm()<=std::numeric_limits<double>::epsilon()) {
        if (restart_origin()) continue;
        result.termination_reason="ciah_macro_energy_rejected";
        break;
      }
    }
    if (!accepted) break;
    const double change=accepted->energy()-energy;
    result.iterations.push_back({iteration,result.hvp_actions-hvp_before,
        micro.keyframe_evaluations,result.ci_actions-ci_before,micro.micro_iterations,
        accepted->energy(),change,accepted_norm,
        std::chrono::duration<double>(Clock::now()-start).count()});
    if (options.observe_iteration) options.observe_iteration(result.iterations.back());
    mcscf::diagnostics::ciah_trace("casscf_macro",
        "iteration=%d,energy=%.17g,de=%.17g,gradient=%.17g,hvp=%d,ci_actions=%d",
        iteration,accepted->energy(),change,accepted_norm,
        result.hvp_actions-hvp_before,result.ci_actions-ci_before);
    point=std::move(accepted);
    if (std::abs(change)<options.energy_tolerance && accepted_norm<options.gradient_tolerance) {
      result.converged=true;
      result.termination_reason="ciah_gradient_tolerance";
      break;
    }
  }
  if (result.termination_reason.empty()) result.termination_reason="maximum_iterations";
  result.energy=point->energy();
  result.gradient_norm=gradient_norm(*point);
  result.orbitals=point->orbitals();
  result.coefficients=point->coefficients();
  return result;
}

}  // namespace xmvb::casscf
