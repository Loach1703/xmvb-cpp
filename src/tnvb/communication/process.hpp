#pragma once

#include <filesystem>
#include <string>
#include "tnvb/communication/response.hpp"
#include "vbscf/core/contracts/input.hpp"
#include "vbscf/integrals/active/preparation/space.hpp"

namespace xmvb::tnvb {

/** Unique calculation/request scratch owner; never uses the old XMVB worker. */
class WorkerProcess {
 public:
  WorkerProcess();
  ~WorkerProcess();
  WorkerProcess(const WorkerProcess&) = delete;
  WorkerProcess& operator=(const WorkerProcess&) = delete;
  WorkerResponse evaluate(const vb::VbScfInput&,
                          const vb::PreparedActiveSpaceContext&,
                          bool gradient) const;
 private:
  std::filesystem::path scratch_;
};

}  // namespace xmvb::tnvb
