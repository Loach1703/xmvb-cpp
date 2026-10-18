#include "mcscf/optimization/ciah/driver.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "core/timing.hpp"
#include "mcscf/optimization/trace.hpp"

namespace xmvb::mcscf {
namespace {

// Initial policies follow PySCF newton_casscf.update_orb_ci. They govern
// approximate microsteps, not the caller's outer accuracy contract.
constexpr double gradient_trust = 3.0;
constexpr int ah_start_cycle = 3;
constexpr int keyframe_interval = 5;
constexpr int base_micro_cycles = 10;

}  // namespace

CiahMicroResult run_ciah_microsteps(CiahWorkspace& workspace,
    const Eigen::VectorXd& initial_gradient, double energy,
    double outer_gradient_tolerance, const CiahMicroCallbacks& callbacks) {
  if (!callbacks.coordinate_norm || !callbacks.metric_action ||
      !callbacks.maximum_component || !callbacks.evaluate_keyframe ||
      callbacks.orbital_size < 0 || callbacks.orbital_size > initial_gradient.size())
    throw std::invalid_argument("CIAH microstep callbacks are incomplete");
  CiahMicroResult result(energy, initial_gradient);
  auto& displacement = result.displacement;
  auto& hessian_displacement = result.hessian_displacement;
  auto& predicted_gradient = result.predicted_gradient;
  auto& last_step = result.last_step;
  auto& candidate = result.candidate;
  auto& energy_model = result.energy_model;
  auto& keyframes = result.keyframes;
  const auto& coordinate_norm = callbacks.coordinate_norm;
  const Eigen::Index no = callbacks.orbital_size;
  const Eigen::Index ns = initial_gradient.size() - no;
  const auto orbital_norm = [&](const Eigen::VectorXd& g) {
    return callbacks.orbital_coordinate_norm
        ? callbacks.orbital_coordinate_norm(g.head(no)) : g.head(no).norm();
  };
  // Retain the finite Ritz residual rather than assuming exact CI stationarity.
  double keyframe_norm = coordinate_norm(initial_gradient);
  // newton_casscf.kernel supplies 0.3 times the outer gradient tolerance.
  const double micro_gradient_tolerance = 0.3 * outer_gradient_tolerance;
  int since_keyframe = 0;
  bool audit_next_keyframe_step = false;
  bool audited_keyframe_step = false;
  bool audited_origin_step = false;
  bool audited_keyframe = false;
  CiahStepPath path;
  while (true) {
    bool replayed = path.pending();
    const double rhs_norm = coordinate_norm(predicted_gradient);
    double replay_residual=0.0;
    if (replayed) {
      replay_residual=coordinate_norm(path.residual_at(predicted_gradient));
      if (!(replay_residual<=rhs_norm/gradient_trust)) {
        diagnostics::ciah_trace("ciah_replay_reset",
            "basis=%d,remaining_fraction=%.17g,current_residual=%.17g,target=%.17g",
            candidate.basis_size,path.remaining_fraction(),replay_residual,rhs_norm/gradient_trust);
        path.clear();
        // The gradient changed after consuming a clipped prefix. Reproject
        // it in the existing span before considering another physical HVP.
        workspace.reset_rhs();
        replayed=false;
      }
    }
    if (!replayed) {
      candidate = workspace.next(predicted_gradient);
      ++result.micro_iterations;
    }
    // A cached-only candidate is usable when the current AH equation is
    // resolved. Otherwise enrich before consuming an approximate step.
    const bool forcing_resolved = replayed ||
        candidate.resolved(rhs_norm, 1.0 / gradient_trust);
    const bool ready = replayed || candidate.converged || candidate.exhaustion ||
        (candidate.basis_size >= ah_start_cycle &&
         (forcing_resolved || candidate.new_hvp_actions != 0));
    if (!ready) continue;

    if (!replayed) path.set(candidate);
    else candidate.converged = false; // Source Ritz root is not the current RHS.
    double largest = callbacks.maximum_component(path.direction());
    largest *= path.remaining_fraction();
    const double step_scale = largest > ciah_max_component ? ciah_max_component / largest : 1.0;
    auto [step, image] = path.take(step_scale);
    if (step_scale < 1.0) candidate.converged = false;
    if (callbacks.audit_step) {
      callbacks.audit_step("micro_origin", displacement, step, result.micro_steps + 1);
      callbacks.audit_step("micro_increment", displacement, step, result.micro_steps + 1);
    }
    if (!audited_origin_step && displacement.isZero() && callbacks.audit_origin) {
      callbacks.audit_origin(displacement, step, image, predicted_gradient);
      audited_origin_step = true;
    }
    if (audit_next_keyframe_step && callbacks.audit_resumed) {
      callbacks.audit_resumed(displacement, step, image, predicted_gradient,
          candidate.basis_size);
      audit_next_keyframe_step = false;
      audited_keyframe_step = true;
    }
    const Eigen::VectorXd previous_gradient = predicted_gradient;
    const Eigen::VectorXd previous_step = last_step;
    displacement += step;
    hessian_displacement += image;
    energy_model.add(step, image);
    predicted_gradient += image;
    last_step = step;
    ++result.micro_steps;
    ++since_keyframe;
    const double predicted_norm = coordinate_norm(predicted_gradient);
    diagnostics::ciah_trace("ciah_micro",
        "basis=%d,micro=%d,omega=%.17g,alpha=%.17g,residual=%.17g,gradient=%.17g,step_scale=%.17g,max_component=%.17g,orbital_residual_euclidean=%.17g,structure_residual_euclidean=%.17g,rhs_norm=%.17g,forcing_target=%.17g,forcing_resolved=%d,replayed=%d,new_hvp=%d,remaining_fraction=%.17g,residual_rhs_current=%d",
        candidate.basis_size, result.micro_steps, candidate.eigenvalue,
        candidate.reference_component, candidate.residual_norm, predicted_norm,
        step_scale, largest, candidate.residual.head(no).norm(), candidate.residual.tail(ns).norm(),
        rhs_norm,rhs_norm / gradient_trust,static_cast<int>(forcing_resolved),
        static_cast<int>(replayed),replayed ? 0 : candidate.new_hvp_actions,
        path.remaining_fraction(),static_cast<int>(!replayed));
    if (result.micro_steps > 3 && predicted_norm > gradient_trust * keyframe_norm) {
      displacement -= step;
      hessian_displacement -= image;
      energy_model.remove(step, image);
      predicted_gradient = previous_gradient;
      last_step = previous_step;
      ++result.rejected_trial_count;
      result.stop_reason = "gradient_growth";
      break;
    }
    const int micro_limit = std::max(base_micro_cycles,
        base_micro_cycles - static_cast<int>(2.0 * std::log(keyframe_norm + 1.0e-7)));
    if (result.micro_steps >= micro_limit ||
        predicted_norm < 0.3 * micro_gradient_tolerance || candidate.exhaustion) {
      result.stop_reason = predicted_norm < 0.3 * micro_gradient_tolerance ? "micro_gradient" :
          candidate.exhaustion ? "subspace_exhaustion" : "micro_budget";
      diagnostics::ciah_trace("ciah_micro_stop",
          "reason=%s,micro=%d,budget=%d,basis=%d,keyframe_norm=%.17g,predicted=%.17g,target=%.17g,residual=%.17g,ah_converged=%d",
          result.stop_reason,result.micro_steps,micro_limit,candidate.basis_size,keyframe_norm,
          predicted_norm,0.3*micro_gradient_tolerance,candidate.residual_norm,
          static_cast<int>(candidate.converged));
      break;
    }

    // The interval uses displacement since the last committed keyframe.
    const auto& keyframe_step = energy_model.displacement();
    const double length = std::sqrt(std::max(0.0,
        keyframe_step.dot(callbacks.metric_action(keyframe_step))));
    const double interval = std::max<double>(keyframe_interval,
        keyframe_interval - std::log(length + 1.0e-7));
    if (since_keyframe >= interval || predicted_norm < keyframe_norm / gradient_trust) {
      if (callbacks.audit_step)
        callbacks.audit_step("keyframe", displacement, displacement, result.micro_steps);
      const CiahKeyframe exact = callbacks.evaluate_keyframe(displacement);
      ++result.keyframe_evaluations;
      if (!audited_keyframe && callbacks.audit_keyframe) {
        callbacks.audit_keyframe(displacement, hessian_displacement, exact);
        audited_keyframe = true;
      }
      const double actual_norm = coordinate_norm(exact.gradient);
      const Eigen::VectorXd correction = exact.gradient - predicted_gradient;
      const double mismatch_norm = coordinate_norm(correction);
      const bool small_gradient = actual_norm < gradient_trust * micro_gradient_tolerance;
      const bool mismatch_small = mismatch_norm < gradient_trust * predicted_norm;
      if (diagnostics::ciah_event_trace_enabled()) {
        const double captured = workspace.projected_dual_norm(correction);
        const double secant = correction.dot(keyframe_step);
        diagnostics::ciah_trace("ciah_keyframe_secant",
            "mismatch_norm=%.17g,captured=%.17g,displacement_norm=%.17g,secant=%.17g,mismatch_small=%d",
            mismatch_norm, captured, length, secant, static_cast<int>(mismatch_small));
      }
      const bool keep_keyframe = small_gradient || mismatch_small;
      // Retain evaluated lower-energy checkpoints even if the model is stale.
      keyframes.save(exact.energy,displacement,hessian_displacement,last_step,candidate);
      diagnostics::ciah_trace("ciah_keyframe",
          "actual=%.17g,predicted=%.17g,mismatch=%.17g,mismatch_small=%d,energy=%.17g,physical_energy=%.17g,accepted=%d",
          actual_norm, predicted_norm, mismatch_norm, static_cast<int>(mismatch_small),
          exact.energy, exact.physical_energy, static_cast<int>(keep_keyframe));
      if (diagnostics::ciah_event_trace_enabled()) {
        const double actual_orbital = orbital_norm(exact.gradient);
        const double predicted_orbital = orbital_norm(predicted_gradient);
        const double mismatch_orbital = orbital_norm(correction);
        diagnostics::ciah_trace("ciah_keyframe_blocks",
            "actual_orbital=%.17g,actual_structure_squared=%.17g,predicted_orbital=%.17g,predicted_structure_squared=%.17g,mismatch_orbital=%.17g,mismatch_structure_squared=%.17g",
            actual_orbital, std::max(0.0,actual_norm*actual_norm-actual_orbital*actual_orbital),
            predicted_orbital, std::max(0.0,predicted_norm*predicted_norm-predicted_orbital*predicted_orbital),
            mismatch_orbital, std::max(0.0,mismatch_norm*mismatch_norm-mismatch_orbital*mismatch_orbital));
      }
      if (keep_keyframe) {
        if (energy_model.descent_mismatch(exact.energy)) {
          result.energy_model_stopped = true;
          result.stop_reason = "keyframe_energy_model";
          diagnostics::ciah_trace("ciah_keyframe_energy_model",
              "anchor=%.17g,predicted=%.17g,actual=%.17g,stop=1",
              energy_model.anchor_energy(),energy_model.energy(),exact.energy);
          break;
        }
        xmvb::profiling::ScopedPhase refresh_phase(xmvb::profiling::Phase::CiahRefresh);
        predicted_gradient = exact.gradient;
        energy_model.refresh(exact.energy, exact.gradient);
        path.clear();
        diagnostics::ciah_trace("ciah_model_refresh",
            "energy=%.17g,source_energy=%.17g,anchored_prediction=%.17g,origin_prediction=%.17g",
            exact.energy, energy, energy - exact.energy,
            -initial_gradient.dot(displacement) -
                0.5 * displacement.dot(hessian_displacement));
        keyframe_norm = actual_norm;
        since_keyframe = 0;
        audit_next_keyframe_step = !audited_keyframe_step &&
            static_cast<bool>(callbacks.audit_resumed);
      } else {
        displacement -= step;
        hessian_displacement -= image;
        energy_model.remove(step, image);
        predicted_gradient = previous_gradient;
        last_step = previous_step;
        ++result.rejected_trial_count;
        result.stop_reason = "keyframe_mismatch";
        break;
      }
    }
    if (!path.pending()) workspace.reset_rhs();
  }

  if (diagnostics::ciah_event_trace_enabled()) {
    const Eigen::VectorXd frozen_gradient = initial_gradient + hessian_displacement;
    diagnostics::ciah_trace("ciah_micro_end",
        "reason=%s,source_coordinate_norm=%.17g,predicted_coordinate_norm=%.17g,frozen_model_coordinate_norm=%.17g,keyframe_correction_norm=%.17g,last_ritz=%.17g",
        result.stop_reason, coordinate_norm(initial_gradient), coordinate_norm(predicted_gradient),
        coordinate_norm(frozen_gradient), coordinate_norm(predicted_gradient - frozen_gradient),
        candidate.residual_norm);
  }
  return result;
}

}  // namespace xmvb::mcscf
