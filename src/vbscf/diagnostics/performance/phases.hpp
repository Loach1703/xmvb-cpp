#pragma once

#include <array>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>

namespace xmvb::vb::diagnostics {

enum class Phase {
  Chart, Projection, Curvature, BlockFactor, BlockEigen,
  TrialForward, TrialGradient, ArmijoAccepted, ArmijoRejected, Count
};

/** @brief Opt-in wall-clock totals; nested phases are inclusive, not additive. */
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
        "trial_forward", "trial_gradient", "armijo_accepted", "armijo_rejected"};
    std::fprintf(stderr, "XMVB_PHASE,stage,calls,wall_seconds\n");
    for (std::size_t i = 0; i < calls.size(); ++i)
      std::fprintf(stderr, "XMVB_PHASE,%s,%llu,%.9f\n", names[i],
                   calls[i].load(std::memory_order_relaxed),
                   nanoseconds[i].load(std::memory_order_relaxed) * 1e-9);
    std::fprintf(stderr, "XMVB_WORKSPACE,projector_peak_payload_bytes,%llu\n",
                 projector_bytes.load(std::memory_order_relaxed));
  }

  bool enabled = false;
  std::atomic<unsigned long long> projector_bytes{0};
  std::array<std::atomic<unsigned long long>, static_cast<std::size_t>(Phase::Count)> calls{};
  std::array<std::atomic<unsigned long long>, static_cast<std::size_t>(Phase::Count)> nanoseconds{};
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
    if (enabled_) start_ = std::chrono::steady_clock::now();
  }
  ScopedPhase(const ScopedPhase&) = delete;
  ScopedPhase& operator=(const ScopedPhase&) = delete;
  ~ScopedPhase() { stop(); }

  void stop() {
    if (!enabled_) return;
    const auto duration = std::chrono::steady_clock::now() - start_;
    const auto i = static_cast<std::size_t>(phase_);
    auto& totals = phase_totals();
    totals.calls[i].fetch_add(1, std::memory_order_relaxed);
    totals.nanoseconds[i].fetch_add(
        std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count(),
        std::memory_order_relaxed);
    enabled_ = false;
  }

 private:
  Phase phase_;
  bool enabled_;
  std::chrono::steady_clock::time_point start_;
};

}  // namespace xmvb::vb::diagnostics
