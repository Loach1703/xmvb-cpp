#include <mpi.h>
#include <sys/resource.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

#include "core/openmp.hpp"
#include "input/deck/model.hpp"
#include "input/deck/primary_basis.hpp"
#include "libcint/materialized_provider.hpp"
#include "vbscf/integrals/ao/pairs/two_electron_index.hpp"
#include "vbscf/integrals/ao/one_electron/direct_operator.hpp"
#include "vbscf/integrals/active/two_electron/transformation/ao_pair_operator.hpp"

namespace {

/** @brief Broadcast molecule metadata, never ERIs or optimizer state. */
template <typename T>
void broadcast(std::vector<T>& values, MPI_Datatype type, int rank) {
  std::uint64_t size = values.size();
  MPI_Bcast(&size, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
  if (size > static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
    throw std::overflow_error("molecule metadata exceeds MPI count limit");
  if (rank != 0) values.resize(size);
  MPI_Bcast(values.data(), static_cast<int>(size), type, 0, MPI_COMM_WORLD);
}

std::uint64_t mix(std::uint64_t value) {
  value ^= value >> 30;
  value *= UINT64_C(0xbf58476d1ce4e5b9);
  value ^= value >> 27;
  value *= UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31);
}

template <typename T>
std::uint64_t bytes(const std::vector<T>& vector) {
  return vector.capacity() * sizeof(T);
}

double probe(std::size_t index, int channel) {
  return std::sin((index + 1) * (channel + 1) * 0.0031);
}

}  // namespace

/**
 * @brief MPI generation/linearity check; deliberately not an SCF driver.
 *
 * Only rank zero parses the deck. Every rank generates its own shell-quartet
 * shard, without scatter/gather of integral values. Scalar fingerprints and
 * independent bilinear operator probes are compared with a one-rank run.
 */
int main(int argc, char** argv) {
  int provided = 0;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  int rank = 0, size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  try {
    if (provided < MPI_THREAD_FUNNELED)
      throw std::runtime_error("MPI does not provide FUNNELED thread support");
    if (argc != 2) throw std::invalid_argument("usage: check_integral_shards input.xmi");
    xmvb::vb::LibcintInput input;
    if (rank == 0) {
      const auto deck = xmvb::vb::parse_input_deck_model(argv[1]);
      input = xmvb::vb::build_input_deck_basis(deck, deck.metadata.basis_name).libcint_input;
    }
    MPI_Bcast(&input.n_atoms, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&input.n_shells, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&input.n_gaussian_primitives, 1, MPI_INT, 0, MPI_COMM_WORLD);
    broadcast(input.atm, MPI_INT, rank);
    broadcast(input.bas, MPI_INT, rank);
    broadcast(input.basidx, MPI_INT, rank);
    broadcast(input.env, MPI_DOUBLE, rank);
    xmvb::vb::LibcintMaterializedIntegralProviderOptions options;
    options.shard_index = rank;
    options.shard_count = size;
    options.include_core_hamiltonian = false;
    MPI_Barrier(MPI_COMM_WORLD);
    const double start = MPI_Wtime();
    auto buffers = xmvb::vb::LibcintMaterializedIntegralProvider().build(input, options);
    const double generation_seconds = MPI_Wtime() - start;
    if (buffers.ao_core_hamiltonian_matrix.size() != 0)
      throw std::runtime_error("worker unexpectedly allocated the core Hamiltonian");
    const std::uint64_t count = buffers.two_electron_values.size();
    std::uint64_t hash_sum = 0, hash_xor = 0;
    for (std::size_t index = 0; index < count; ++index) {
      std::uint64_t bits = 0;
      static_assert(sizeof(bits) == sizeof(double));
      std::memcpy(&bits, &buffers.two_electron_values[index], sizeof(bits));
      const auto key = (static_cast<std::uint64_t>(buffers.left_pair_indices[index]) << 32) |
                       static_cast<std::uint64_t>(buffers.right_pair_indices[index]);
      const auto hash = mix(key) ^ mix(bits);
      hash_sum += hash;
      hash_xor ^= hash;
    }
    const std::uint64_t payload_bytes = bytes(buffers.two_electron_values) +
        bytes(buffers.left_pair_indices) + bytes(buffers.right_pair_indices);
    rusage generation_usage{};
    getrusage(RUSAGE_SELF, &generation_usage);
    const int n_bf = buffers.n_basis_functions;
    xmvb::vb::AoIntegralInput ao;
    ao.n_basis_functions = n_bf;
    ao.pair_graph = xmvb::vb::build_ao_pair_graph(
        std::move(buffers.left_pair_indices), std::move(buffers.right_pair_indices),
        std::move(buffers.two_electron_values), n_bf);
    const auto& graph = ao.pair_graph;
    const std::uint64_t graph_bytes = bytes(graph.values) + bytes(graph.columns) +
        bytes(graph.row_offsets) + bytes(graph.integral_rows) + bytes(graph.integral_edges) +
        bytes(graph.integral_columns) + bytes(graph.integral_values) +
        bytes(graph.pair_first) + bytes(graph.pair_second);
    // These diagnostic probes have a strict 8 MiB input/image bound. They are
    // not a production plan that replicates arbitrary active-space matrices.
    const std::size_t n_pairs = static_cast<std::size_t>(n_bf) * (n_bf + 1) / 2;
    if (n_pairs * 3 * sizeof(double) > 8 * 1024 * 1024)
      throw std::runtime_error("diagnostic pair probe exceeds its workspace bound");
    xmvb::vb::ExactCtxPairMatrix source(n_pairs, 3), image;
    for (Eigen::Index row = 0; row < source.rows(); ++row)
      for (int channel = 0; channel < 3; ++channel)
        source(row, channel) = probe(row, channel);
    xmvb::vb::detail::apply_exact_ao_pair_kernel(ao, source, n_bf, 3, &image);
    double contractions[5]{};
    for (int channel = 0; channel < 3; ++channel)
      contractions[channel] = source.col(channel).dot(image.col(channel));
    std::vector<double> density(static_cast<std::size_t>(n_bf) * n_bf);
    for (std::size_t index = 0; index < density.size(); ++index)
      density[index] = probe(index, 4);
    const auto h1e = xmvb::vb::apply_ao_h1e(
        density.data(), ao, xmvb::effective_openmp_thread_count());
    const auto transpose = xmvb::vb::apply_ao_h1e_transpose(
        density.data(), ao, xmvb::effective_openmp_thread_count());
    for (std::size_t index = 0; index < density.size(); ++index) {
      contractions[3] += probe(index, 5) * h1e[index];
      contractions[4] += probe(index, 5) * transpose[index];
    }
    double totals[5]{};
    MPI_Reduce(contractions, totals, 5, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    std::uint64_t total_count = 0, total_sum = 0, total_xor = 0;
    MPI_Reduce(&count, &total_count, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&hash_sum, &total_sum, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&hash_xor, &total_xor, 1, MPI_UINT64_T, MPI_BXOR, 0, MPI_COMM_WORLD);
    double max_seconds = 0;
    MPI_Reduce(&generation_seconds, &max_seconds, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    rusage usage{};
    getrusage(RUSAGE_SELF, &usage);
    // Rank-labelled output gives unambiguous memory and ownership evidence;
    // a launcher may deliver different ranks' stdout out of order.
    for (int owner = 0; owner < size; ++owner) {
      MPI_Barrier(MPI_COMM_WORLD);
      if (rank == owner) {
        char name[MPI_MAX_PROCESSOR_NAME];
        int length = 0;
        MPI_Get_processor_name(name, &length);
        std::cout << std::setprecision(17) << "rank=" << rank << " node="
                  << std::string(name, length) << " integrals=" << count
                  << " payload_bytes=" << payload_bytes << " graph_bytes=" << graph_bytes
                  << " generation_seconds=" << generation_seconds
                  << " generation_rss_kib=" << generation_usage.ru_maxrss
                  << " peak_rss_kib=" << usage.ru_maxrss << std::endl;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
      std::cout << "ranks=" << size << " total_integrals=" << total_count
                << " hash_sum=" << total_sum << " hash_xor=" << total_xor
                << " max_generation_seconds=" << max_seconds << '\n';
      for (int channel = 0; channel < 5; ++channel)
        std::cout << "probe_" << channel << '=' << totals[channel] << '\n';
    }
    MPI_Finalize();
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "rank " << rank << ": " << error.what() << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 2);
    return 2;
  }
}
