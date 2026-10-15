#pragma once

#include <filesystem>
#include <fstream>

namespace xmvb::vb {
struct NeoIterationRecord;
}

namespace xmvb::output {

/** Write accepted NEO iteration diagnostics as tab-separated rows. */
class NeoTraceWriter {
 public:
  explicit NeoTraceWriter(const std::filesystem::path& output_path);

  void write(const vb::NeoIterationRecord& record);

 private:
  std::ofstream stream_;
  std::filesystem::path output_path_;
};

}  // namespace xmvb::output
