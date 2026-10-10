#pragma once

#include <array>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <ctime>

namespace xmvb::vb::diagnostics {

enum class Phase {
  Chart, Projection, Curvature, BlockFactor, BlockEigen,
  TrialForward, TrialGradient, ArmijoAccepted, ArmijoRejected,
  NeoSolve, HvpApply, HvpSetup, AoResponse, ActiveDirection,
  SelectedDirection, ActiveBackward, OrbitalPullback, TwoElectronResponse,
  FixedPullback, StructureSolve, StructureAdjoint,
  NeoProjectedSolve, NeoPreconditioner, NeoBasis, NeoMetric, Count
};

/** @brief Opt-in wall/process-CPU totals. Nested scopes are not additive.
 * Process CPU includes all worker threads, including spin/wait time. Scopes
 * belong at calling-thread phase boundaries, not inside parallel worker loops.
 */
class PhaseTotals {
 public:
  PhaseTotals() {
    const char* value = std::getenv("XMVB_PHASE_TIMING");
    enabled = value != nullptr && value[0] == '1';
  }

  ~PhaseTotals() {
    if (!enabled) return;
    constexpr const char* names[] = {
        "chart", "projection", "curvature", "block_factor", "block_eigen",
        "trial_forward", "trial_gradient", "armijo_accepted", "armijo_rejected",
        "neo_solve", "hvp_apply", "hvp_setup", "ao_response", "active_direction",
        "selected_direction", "active_backward", "orbital_pullback", "two_electron_response",
        "fixed_pullback", "structure_solve", "structure_adjoint",
        "neo_projected_solve", "neo_preconditioner", "neo_basis", "neo_metric"};
    static_assert(sizeof(names) / sizeof(names[0]) == static_cast<std::size_t>(Phase::Count));
    std::fprintf(stderr, "XMVB_PHASE,stage,calls,wall_seconds\n");
    for (std::size_t i = 0; i < calls.size(); ++i)
      std::fprintf(stderr, "XMVB_PHASE,%s,%llu,%.9f\n", names[i],
                   calls[i].load(std::memory_order_relaxed),
                   nanoseconds[i].load(std::memory_order_relaxed) * 1e-9);
    std::fprintf(stderr, "XMVB_CPU,stage,cpu_seconds\n");
    for (std::size_t i = 0; i < calls.size(); ++i)
      std::fprintf(stderr, "XMVB_CPU,%s,%.9f\n", names[i],
                   cpu_ticks[i].load(std::memory_order_relaxed) /
                       static_cast<double>(CLOCKS_PER_SEC));
    std::fprintf(stderr, "XMVB_WORKSPACE,projector_peak_payload_bytes,%llu\n",
                 projector_bytes.load(std::memory_order_relaxed));
  }

  bool enabled = false;
  std::atomic<unsigned long long> projector_bytes{0};
  std::array<std::atomic<unsigned long long>, static_cast<std::size_t>(Phase::Count)> calls{};
  std::array<std::atomic<unsigned long long>, static_cast<std::size_t>(Phase::Count)> nanoseconds{};
  std::array<std::atomic<unsigned long long>, static_cast<std::size_t>(Phase::Count)> cpu_ticks{};
};

inline PhaseTotals& phase_totals() {
  static PhaseTotals totals;
  return totals;
}

/** @brief Records an event count without changing the optimization path. */
inline void count_phase(Phase phase) {
  auto& totals = phase_totals();
  if (totals.enabled)
    totals.calls[static_cast<std::size_t>(phase)].fetch_add(1, std::memory_order_relaxed);
}

/** @brief Largest projector-owned matrix payload, excluding referenced inputs. */
inline void record_projector_bytes(unsigned long long bytes) {
  auto& totals = phase_totals();
  if (!totals.enabled) return;
  auto peak = totals.projector_bytes.load(std::memory_order_relaxed);
  while (peak < bytes && !totals.projector_bytes.compare_exchange_weak(
             peak, bytes, std::memory_order_relaxed)) {}
}

/** @brief Times a scope only when XMVB_PHASE_TIMING=1. */
class ScopedPhase {
 public:
  explicit ScopedPhase(Phase phase) : phase_(phase), enabled_(phase_totals().enabled) {
    if (enabled_) {
      start_ = std::chrono::steady_clock::now();
      cpu_start_ = std::clock();
    }
  }
  ScopedPhase(const ScopedPhase&) = delete;
  ScopedPhase& operator=(const ScopedPhase&) = delete;
  ~ScopedPhase() { stop(); }

  void stop() {
    if (!enabled_) return;
    const auto cpu_end = std::clock();
    const auto duration = std::chrono::steady_clock::now() - start_;
    const auto i = static_cast<std::size_t>(phase_);
    auto& totals = phase_totals();
    totals.calls[i].fetch_add(1, std::memory_order_relaxed);
    totals.nanoseconds[i].fetch_add(
        std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count(),
        std::memory_order_relaxed);
    if (cpu_start_ != static_cast<std::clock_t>(-1) && cpu_end >= cpu_start_)
      totals.cpu_ticks[i].fetch_add(cpu_end - cpu_start_, std::memory_order_relaxed);
    enabled_ = false;
  }

 private:
  Phase phase_;
  bool enabled_;
  std::chrono::steady_clock::time_point start_;
  std::clock_t cpu_start_{};
};

}  // namespace xmvb::vb::diagnostics
