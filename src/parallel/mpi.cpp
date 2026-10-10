#include "parallel/mpi.hpp"

#include <mpi.h>
#include <sys/resource.h>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <omp.h>

#include "core/openmp.hpp"
#include "libcint/materialized_provider.hpp"
#include "vbscf/diagnostics/performance/phases.hpp"
#include "vbscf/integrals/ao/pairs/two_electron_index.hpp"
#include "vbscf/integrals/active/two_electron/transformation/ao_pair_operator.hpp"

namespace xmvb::parallel {
namespace {
enum Command { Stop, Build, Pair, Forward, Transpose, Fused };
constexpr std::size_t tile_bytes = 8 * 1024 * 1024;
constexpr std::size_t h1e_private_bytes = 32 * 1024 * 1024;

void main_thread() {
  int main = 0;
  MPI_Is_thread_main(&main);
  if (omp_in_parallel() || !main)
    throw std::runtime_error("MPI integral actions require the main thread outside OpenMP regions");
}

template <typename T>
void broadcast(std::vector<T>& data, MPI_Datatype type, int rank) {
  std::uint64_t count = data.size();
  MPI_Bcast(&count, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
  if (count > static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
    throw std::overflow_error("MPI molecule metadata exceeds int count");
  if (rank != 0) data.resize(count);
  MPI_Bcast(data.data(), static_cast<int>(count), type, 0, MPI_COMM_WORLD);
}

void metadata(vb::LibcintInput& input, int rank) {
  MPI_Bcast(&input.n_atoms, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&input.n_shells, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&input.n_gaussian_primitives, 1, MPI_INT, 0, MPI_COMM_WORLD);
  broadcast(input.atm, MPI_INT, rank);
  broadcast(input.bas, MPI_INT, rank);
  broadcast(input.basidx, MPI_INT, rank);
  broadcast(input.env, MPI_DOUBLE, rank);
}

void reduce(std::vector<double>& values, int rank) {
  MPI_Reduce(rank == 0 ? MPI_IN_PLACE : values.data(),
             rank == 0 ? values.data() : nullptr,
             static_cast<int>(values.size()), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}
}  // namespace

MpiSession::MpiSession(int& argc, char**& argv) {
  int provided = 0;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &size_);
  // Root owns the SCF phase report; worker totals must not look like root totals.
  if (rank_ != 0) vb::diagnostics::phase_totals().enabled = false;
  MPI_Comm local;
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank_, MPI_INFO_NULL, &local);
  int local_size = 0;
  MPI_Comm_size(local, &local_size);
  MPI_Comm_free(&local);
  int invalid = provided < MPI_THREAD_FUNNELED || local_size != 1;
  MPI_Allreduce(MPI_IN_PLACE, &invalid, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (invalid) {
    if (rank_ == 0) std::cerr << "MPI requires FUNNELED support and one rank per node\n";
    MPI_Abort(MPI_COMM_WORLD, 2);
  }
  if (size_ > 1 && rank_ == 0) vb::IntegralService::current = this;
}

MpiSession::~MpiSession() {
  if (size_ > 1 && rank_ == 0) {
    int header[3] = {Stop, 0, 0};
    MPI_Bcast(header, 3, MPI_INT, 0, MPI_COMM_WORLD);
  }
  vb::IntegralService::current = nullptr;
  rusage usage{};
  getrusage(RUSAGE_SELF, &usage);
  std::cerr << "XMVB_MPI_RANK," << rank_ << ",integrals," << owned_integrals_
            << ",peak_rss_kib," << usage.ru_maxrss << '\n';
  if (rank_ == 0 && size_ > 1)
    std::cerr << "XMVB_MPI_ACTIONS," << commands_ << ",seconds," << action_seconds_
              << ",peak_tile_bytes," << peak_tile_bytes_ << '\n';
  MPI_Finalize();
}

vb::AoIntegralInput MpiSession::generate(vb::LibcintInput input) {
  metadata(input, rank_);
  vb::LibcintMaterializedIntegralProviderOptions options;
  options.shard_index = rank_;
  options.shard_count = size_;
  auto buffers = vb::LibcintMaterializedIntegralProvider().build(input, options);
  owned_integrals_ = buffers.two_electron_values.size();
  vb::AoIntegralInput ao;
  ao.n_basis_functions = buffers.n_basis_functions;
  ao.ao_core_hamiltonian_matrix = std::move(buffers.ao_core_hamiltonian_matrix);
  ao.pair_graph = vb::build_ao_pair_graph(std::move(buffers.left_pair_indices),
      std::move(buffers.right_pair_indices), std::move(buffers.two_electron_values),
      ao.n_basis_functions);
  return ao;
}

vb::AoIntegralInput MpiSession::build(const vb::LibcintInput& input) {
  main_thread();
  send(Build, 0, 0);
  try {
    auto ao = generate(input);
    ao.pair_graph.service = this;
    return ao;
  } catch (...) { MPI_Abort(MPI_COMM_WORLD, 2); throw; }
}

void MpiSession::send(int command, int rows, int columns) {
  main_thread();
  int header[3] = {command, rows, columns};
  MPI_Bcast(header, 3, MPI_INT, 0, MPI_COMM_WORLD);
}

bool MpiSession::run_workers() {
  if (rank_ == 0) return false;
  try {
    while (true) {
      int header[3];
      MPI_Bcast(header, 3, MPI_INT, 0, MPI_COMM_WORLD);
      if (header[0] == Stop) return true;
      if (header[0] == Build) worker_ = generate({});
      else {
        receive(header[0], header[1], header[2]);
        if (header[0] == Pair)
          execute_pair(worker_.pair_graph, worker_.n_basis_functions, header[1], header[2]);
        else execute(header[0], worker_, header[1], header[2]);
      }
    }
  } catch (const std::exception& error) {
    std::cerr << "MPI worker " << rank_ << ": " << error.what() << '\n';
    MPI_Abort(MPI_COMM_WORLD, 2);
    throw;
  }
}

void MpiSession::receive(int command, int rows, int columns) {
  if (rows <= 0 || columns <= 0 ||
      static_cast<std::size_t>(rows) * columns > tile_bytes / sizeof(double))
    throw std::invalid_argument("MPI contraction tile exceeds memory bound");
  const int count = rows * columns;
  source_.resize(count);
  MPI_Bcast(source_.data(), count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (command == Fused) {
    adjoint_.resize(count);
    MPI_Bcast(adjoint_.data(), count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  peak_tile_bytes_ = std::max(peak_tile_bytes_, source_.size() * sizeof(double));
}

void MpiSession::execute(int command, const vb::AoIntegralInput& ao, int rows, int columns) {
  dispatching_ = true;
  const int threads = std::max(1, std::min(xmvb::effective_openmp_thread_count(),
      static_cast<int>(h1e_private_bytes / (2 * sizeof(double) * rows))));
  try {
    if (command == Fused) {
      vb::apply_ao_h1e_fused(source_.data(), adjoint_.data(), ao, threads,
                           &h1e_workspace_, &forward_, &transpose_);
      reduce(forward_, rank_);
      reduce(transpose_, rank_);
    } else if (command == Forward || command == Transpose) {
      forward_ = command == Forward ? vb::apply_ao_h1e(source_.data(), ao, threads) :
                                     vb::apply_ao_h1e_transpose(source_.data(), ao, threads);
      reduce(forward_, rank_);
    } else throw std::invalid_argument("unknown MPI integral command");
    dispatching_ = false;
  } catch (...) { MPI_Abort(MPI_COMM_WORLD, 2); throw; }
}

std::vector<double> MpiSession::h1e(const vb::AoIntegralInput& ao, const double* input, bool transpose) {
  main_thread();
  const std::size_t count = static_cast<std::size_t>(ao.n_basis_functions) * ao.n_basis_functions;
  if (!input || count == 0 || count > tile_bytes / sizeof(double))
    throw std::invalid_argument("MPI AO-H1E input exceeds supported tile size");
  const double start = MPI_Wtime();
  source_.assign(input, input + count);
  const int command = transpose ? Transpose : Forward;
  send(command, count, 1);
  receive(command, count, 1);
  execute(command, ao, count, 1);
  ++commands_;
  action_seconds_ += MPI_Wtime() - start;
  return forward_;
}

void MpiSession::fused(const vb::AoIntegralInput& ao, const double* input, const double* adjoint,
                       std::vector<double>* forward, std::vector<double>* transpose) {
  main_thread();
  const std::size_t count = static_cast<std::size_t>(ao.n_basis_functions) * ao.n_basis_functions;
  if (!input || !adjoint || !forward || !transpose || count == 0 || count > tile_bytes / sizeof(double))
    throw std::invalid_argument("invalid or oversized MPI fused AO-H1E input");
  const double start = MPI_Wtime();
  source_.assign(input, input + count);
  adjoint_.assign(adjoint, adjoint + count);
  send(Fused, count, 1);
  receive(Fused, count, 1);
  execute(Fused, ao, count, 1);
  *forward = forward_;
  *transpose = transpose_;
  ++commands_;
  action_seconds_ += MPI_Wtime() - start;
}

void MpiSession::execute_pair(const vb::AoPairGraph& graph, int n_bf, int rows, int columns) {
  dispatching_ = true;
  try {
    const auto source = Eigen::Map<const vb::PairMatrix>(source_.data(), rows, columns);
    const int saved_threads = omp_get_max_threads();
    if (graph.has_compact_integral_stream()) omp_set_num_threads(1);
    vb::detail::apply_exact_ao_pair_graph(graph, source, n_bf, columns, &pair_image_);
    omp_set_num_threads(saved_threads);
    MPI_Reduce(rank_ == 0 ? MPI_IN_PLACE : pair_image_.data(),
        rank_ == 0 ? pair_image_.data() : nullptr, rows * columns,
        MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    dispatching_ = false;
  } catch (...) { MPI_Abort(MPI_COMM_WORLD, 2); throw; }
}

void MpiSession::pair(const vb::AoPairGraph& graph, int n_bf, const Eigen::Ref<const vb::PairMatrix>& input,
                      vb::PairMatrix* result) {
  main_thread();
  const Eigen::Index rows = static_cast<Eigen::Index>(n_bf) * (n_bf + 1) / 2;
  if (!result || input.rows() != rows || rows <= 0 ||
      rows > static_cast<Eigen::Index>(tile_bytes / sizeof(double)))
    throw std::invalid_argument("invalid or oversized MPI pair input");
  result->resize(input.rows(), input.cols());
  const Eigen::Index width = tile_bytes / sizeof(double) / rows;
  for (Eigen::Index begin = 0; begin < input.cols(); begin += width) {
    const double start = MPI_Wtime();
    const Eigen::Index columns = std::min(width, input.cols() - begin);
    source_.resize(rows * columns);
    Eigen::Map<vb::PairMatrix>(source_.data(), rows, columns) = input.middleCols(begin, columns);
    send(Pair, rows, columns);
    receive(Pair, rows, columns);
    execute_pair(graph, n_bf, rows, columns);
    result->middleCols(begin, columns) = pair_image_;
    ++commands_;
    action_seconds_ += MPI_Wtime() - start;
  }
}

}  // namespace xmvb::parallel
