#pragma once

#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace xmvb::mcscf::diagnostics {

/** @brief Full NEO diagnostics, including explicit spectral/model audits. */
inline bool neo_trace_enabled() {
  static const bool enabled = [] {
    const char* value = std::getenv("XMVB_NEO_TRACE");
    return value && value[0] == '1' && value[1] == '\0';
  }();
  return enabled;
}

/** @brief Existing events only; scalar mode leaves the native numerical path unchanged. */
inline bool neo_event_trace_enabled() {
  static const bool enabled = [] {
    const char* value = std::getenv("XMVB_NEO_TRACE");
    return value && (std::strcmp(value, "1") == 0 ||
                     std::strcmp(value, "scalar") == 0);
  }();
  return enabled;
}

/** @brief Expensive nonlinear candidate audit; never a production policy. */
inline bool neo_candidate_audit_enabled() {
  static const bool enabled = [] {
    const char* value = std::getenv("XMVB_NEO_CANDIDATE_AUDIT");
    return neo_trace_enabled() && value && value[0] == '1' && value[1] == '\0';
  }();
  return enabled;
}

/** @brief Identifies an accepted-point solve and its radius retry. */
struct NeoTraceContext {
  int keyframe = 0;
  int attempt = 0;
  const char* source = "initial";
};

inline thread_local NeoTraceContext neo_trace_context;

/** @brief Writes existing scalar diagnostics without retaining vectors. */
template<class... Args>
inline void neo_trace(const char* event, const char* format, Args... args) {
  if (!neo_event_trace_enabled()) return;
  std::fprintf(stderr, "NEO_DIAG,%d,%d,%s,",
               neo_trace_context.keyframe, neo_trace_context.attempt, event);
  std::fprintf(stderr, format, args...);
  std::fputc('\n', stderr);
}

}  // namespace xmvb::mcscf::diagnostics
