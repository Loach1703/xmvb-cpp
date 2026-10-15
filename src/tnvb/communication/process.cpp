#include "tnvb/communication/process.hpp"
#include "tnvb/config.hpp"

#include <array>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <spawn.h>
#include <stdexcept>
#include <system_error>
#include <sys/wait.h>
#include <unistd.h>

extern char** environ;

namespace xmvb::tnvb {
namespace {

std::filesystem::path unique_directory(const std::filesystem::path& parent) {
  std::string pattern = (parent / "xmvb-cpp-tnvb-XXXXXX").string();
  if (::mkdtemp(pattern.data()) == nullptr) {
    throw std::system_error(errno, std::generic_category(), "TNVB scratch creation");
  }
  return pattern;
}

template <class T>
void write_values(std::ostream& stream, const T* data, std::size_t count) {
  stream.write(reinterpret_cast<const char*>(data), count * sizeof(T));
  if (!stream) throw std::runtime_error("TNVB request write failed");
}

template <class T>
void read_values(std::istream& stream, T* data, std::size_t count) {
  stream.read(reinterpret_cast<char*>(data), count * sizeof(T));
  if (!stream) throw std::runtime_error("truncated TNVB worker response");
}

std::filesystem::path worker_path() {
  // Prefer the matched installed worker relative to the executable; builds
  // use their own staged copy. Never fall back to xmvb4.0 source files.
  std::array<char, 4096> executable{};
  const ssize_t length = ::readlink("/proc/self/exe", executable.data(), executable.size()-1);
  if (length > 0) {
    const auto candidate = std::filesystem::path(std::string(executable.data(), length))
        .parent_path().parent_path() / "libexec/xmvb-cpp/tnvb/evaluate.py";
    if (std::filesystem::exists(candidate)) return candidate;
  }
  if (std::filesystem::exists(XMVB_TNVB_WORKER)) return XMVB_TNVB_WORKER;
  throw std::runtime_error("matched independent TNVB worker is missing");
}

void execute_worker(const std::filesystem::path& request,
                    const std::filesystem::path& response,
                    const std::filesystem::path& cache) {
  std::map<std::string, std::string> environment;
  for (char** item = environ; *item != nullptr; ++item) {
    const std::string entry(*item);
    const auto equals = entry.find('=');
    if (equals != std::string::npos) {
      environment.emplace(entry.substr(0, equals), entry.substr(equals + 1));
    }
  }
  environment["PYTHONPATH"] = std::string(XMVB_TNVB_BLOCK2_MODULE) + ":" +
      XMVB_TNVB_BLOCK2_SOURCE;
  const auto previous_preload=environment.find("LD_PRELOAD");
  const std::string inherited_preload=previous_preload==environment.end()
      ? std::string{} : previous_preload->second;
  environment["LD_PRELOAD"] = std::string(XMVB_TNVB_LIBRARY_DIR) +
      "/libstdc++.so.6:" + XMVB_TNVB_LIBRARY_DIR + "/libgomp.so.1" +
      (inherited_preload.empty() ? std::string{} : ":" + inherited_preload);
  environment["XMVB_BLOCK2_WARM_CACHE"] = cache.string();
  environment["XMVB_BLOCK2_SCRATCH_PARENT"] = request.parent_path().string();
  environment["XMVB_BLOCK2_RDM_CHECK_TOL"] = "1e-7";
  environment["OPENBLAS_NUM_THREADS"] = "1";
  environment["OMP_DYNAMIC"] = "FALSE";
  environment["OMP_MAX_ACTIVE_LEVELS"] = "1";
  std::vector<std::string> entries;
  entries.reserve(environment.size());
  for (const auto& item : environment) entries.push_back(item.first + "=" + item.second);
  std::vector<char*> envp;
  envp.reserve(entries.size() + 1);
  for (auto& item : entries) envp.push_back(item.data());
  envp.push_back(nullptr);
  std::string python = XMVB_TNVB_PYTHON;
  std::string worker = worker_path().string();
  std::string in = request.string(), out = response.string();
  std::array<char*, 5> argv{python.data(), worker.data(), in.data(), out.data(), nullptr};
  // Flush parent reports before the child writes to the inherited stdout.
  // Otherwise redirected logs can splice worker summaries into an orbital row.
  std::cout.flush();
  std::cerr.flush();
  pid_t child = -1;
  const int error = ::posix_spawn(&child, python.c_str(), nullptr, nullptr,
                                 argv.data(), envp.data());
  if (error != 0) {
    throw std::system_error(error, std::generic_category(), "TNVB worker launch");
  }
  int status = 0;
  while (::waitpid(child, &status, 0) == -1) {
    if (errno != EINTR) {
      throw std::system_error(errno, std::generic_category(), "TNVB worker wait");
    }
  }
  if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
    throw std::runtime_error("TNVB worker failed (wait status " + std::to_string(status) + ")");
  }
}

}  // namespace

WorkerProcess::WorkerProcess() {
  const std::uint16_t endian = 1;
  if (*reinterpret_cast<const unsigned char*>(&endian) != 1) {
    throw std::runtime_error("TNVB binary protocol requires a little-endian host");
  }
  const char* parent = std::getenv("XMVB_BLOCK2_SCRATCH_PARENT");
  scratch_ = unique_directory(parent && *parent ? parent : std::filesystem::temp_directory_path());
}

WorkerProcess::~WorkerProcess() {
  std::error_code ignored;
  std::filesystem::remove_all(scratch_, ignored);
}

WorkerResponse WorkerProcess::evaluate(const vb::VbScfInput& input,
    const vb::PreparedActiveSpaceContext& prepared, bool gradient) const {
  const auto directory = unique_directory(scratch_);
  const auto request = directory / "request.bin", response = directory / "response.bin";
  const auto& topology = input.raw_structure_data;
  const int n_str = topology.n_structures;
  const int n_ao = input.orbital_preparation_input.n_active_orbitals;
  const std::size_t n_pair = std::size_t(n_ao) * (n_ao + 1) / 2;
  const std::size_t n_eri = n_pair * (n_pair + 1) / 2;
  if (n_str <= 0 || n_ao <= 0 ||
      topology.raw_structure_orbitals.size() != std::size_t(n_str)*topology.n_total_electrons) {
    throw std::invalid_argument("TNVB raw structure dimensions are inconsistent");
  }
  std::vector<double> packed = prepared.active_space_two_electron_result.packed_active_two_electron_integrals;
  if (packed.empty() && prepared.active_space_two_electron_result.representation ==
      vb::ActiveSpaceTwoElectronRepresentation::ResolutionOfIdentity) {
    const auto& factors = prepared.active_space_two_electron_result.ri_active_pair_factors;
    if (factors.cols() != static_cast<Eigen::Index>(n_pair)) {
      throw std::invalid_argument("TNVB RI active-pair dimensions are inconsistent");
    }
    packed.resize(n_eri);
    for (std::size_t q = 0; q < n_pair; ++q) {
      for (std::size_t p = 0; p <= q; ++p) {
        packed[q*(q+1)/2+p] = factors.col(p).dot(factors.col(q));
      }
    }
  }
  if (packed.size() != n_eri) throw std::invalid_argument("TNVB active ERI size mismatch");
  {
    std::ofstream stream(request, std::ios::binary);
    stream.write("CPTNREQ1", 8);
    const std::array<std::int32_t, 8> header{1, n_str, n_ao, topology.n_active_electrons,
        topology.n_total_electrons, input.tnvb.structure_bond, input.tnvb.total_bond,
        gradient ? 1 : 0};
    write_values(stream, header.data(), header.size());
    // RawStructureData is already structure-major, i.e. Fortran (NEL,NSTR).
    for (const int orbital : topology.raw_structure_orbitals) {
      const std::int32_t value = orbital;
      write_values(stream, &value, 1);
    }
    const auto& overlap = prepared.orbital_result.active_orbital_overlap_matrix;
    const auto& h = prepared.active_space_one_electron_result.h1e_act;
    if (overlap.size() != std::size_t(n_ao)*n_ao || h.size() != n_ao*n_ao) {
      throw std::invalid_argument("TNVB active matrix size mismatch");
    }
    write_values(stream, overlap.data(), overlap.size());
    write_values(stream, h.data(), h.size());
    write_values(stream, packed.data(), packed.size());
  }
  execute_worker(request, response, scratch_ / "warm");
  WorkerResponse result;
  std::ifstream stream(response, std::ios::binary);
  std::array<char, 8> magic{};
  read_values(stream, magic.data(), magic.size());
  std::array<std::int32_t, 8> header{};
  read_values(stream, header.data(), header.size());
  if (std::memcmp(magic.data(), "CPTNRES1", 8) || header[0] != 1 || header[1] != 0 ||
      header[2] != n_str || header[3] != n_ao || header[4] != (gradient ? 1 : 0) ||
      header[5] != n_str || header[6] || header[7]) {
    throw std::runtime_error("TNVB response header mismatch");
  }
  result.hamiltonian.resize(n_str, n_str); result.overlap.resize(n_str, n_str);
  result.energies.resize(n_str); result.coefficients.resize(n_str, n_str);
  read_values(stream, result.hamiltonian.data(), result.hamiltonian.size());
  read_values(stream, result.overlap.data(), result.overlap.size());
  read_values(stream, result.energies.data(), result.energies.size());
  read_values(stream, result.coefficients.data(), result.coefficients.size());
  std::array<double, 4> checks{};
  read_values(stream, checks.data(), checks.size());
  result.fit_energy_error=checks[0]; result.rdm_identity_error=checks[1];
  result.hs_wall_seconds=checks[2]; result.rdm_wall_seconds=checks[3];
  if (gradient) {
    result.density_one.resize(n_ao, n_ao); result.contracted_three.resize(n_ao, n_ao);
    result.density_two_buffer.resize(std::size_t(n_ao)*n_ao*n_ao*n_ao);
    read_values(stream, result.density_one.data(), result.density_one.size());
    read_values(stream, result.density_two_buffer.data(), result.density_two_buffer.size());
    read_values(stream, result.contracted_three.data(), result.contracted_three.size());
  }
  if (stream.peek() != std::char_traits<char>::eof()) {
    throw std::runtime_error("trailing TNVB response data");
  }
  std::error_code ignored;
  std::filesystem::remove_all(directory, ignored);
  return result;
}

}  // namespace xmvb::tnvb
