#include <cmath>
#include <algorithm>
#include <atomic>
#include <iostream>
#include <map>
#include <stdexcept>

extern "C" {
#include "cint.h"
}

#include "libcint/ri_provider.hpp"
#include "libcint/materialized_provider.hpp"
#include "libcint/direct_shell.hpp"
#include "libcint/shell_normalization.hpp"

namespace {

void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

/** @brief One normalized s Gaussian; duplicate shells make the metric singular. */
xmvb::vb::LibcintInput make_input(int n_shells) {
  xmvb::vb::LibcintInput input;
  input.n_atoms = 1;
  input.n_shells = n_shells;
  input.n_gaussian_primitives = n_shells;
  input.atm.assign(ATM_SLOTS, 0);
  input.atm[CHARGE_OF] = 1;
  input.atm[PTR_COORD] = PTR_ENV_START;
  input.env.assign(PTR_ENV_START + 3, 0.0);
  input.bas.assign(n_shells * BAS_SLOTS, 0);
  input.basidx.resize(2 * n_shells);
  for (int shell = 0; shell < n_shells; ++shell) {
    const int offset = static_cast<int>(input.env.size());
    input.env.push_back(1.0);
    input.env.push_back(xmvb::vb::shell_normalization(0, 1.0));
    input.bas[shell * BAS_SLOTS + NPRIM_OF] = 1;
    input.bas[shell * BAS_SLOTS + NCTR_OF] = 1;
    input.bas[shell * BAS_SLOTS + PTR_EXP] = offset;
    input.bas[shell * BAS_SLOTS + PTR_COEFF] = offset + 1;
    input.basidx[2 * shell] = shell;
    input.basidx[2 * shell + 1] = 1;
  }
  return input;
}

/** @brief Independent legacy integral reference, including Cartesian d shells. */
void check_shell_workspace() {
  auto input = make_input(3);
  const int counts[] = {1, 3, 6};
  int offset = 0;
  for (int shell = 0; shell < 3; ++shell) {
    input.bas[shell * BAS_SLOTS + ANG_OF] = shell;
    input.env[input.bas[shell * BAS_SLOTS + PTR_COEFF]] =
        xmvb::vb::shell_normalization(shell, 1.0);
    input.basidx[2 * shell] = offset;
    input.basidx[2 * shell + 1] = counts[shell];
    offset += counts[shell];
  }
  xmvb::vb::LibcintDirectShellEvaluator evaluator(input);
  const auto& norm = evaluator.ao_normalization();
  std::vector<std::vector<double>> references(81);
  for (int task = 0; task < 81; ++task) {
    FINT shells[] = {task / 27, (task / 9) % 3, (task / 3) % 3, task % 3};
    auto& values = references[task];
    values.resize(counts[shells[0]] * counts[shells[1]] *
                  counts[shells[2]] * counts[shells[3]]);
    cint2e_cart(values.data(), shells, input.atm.data(), input.n_atoms,
                input.bas.data(), input.n_shells, input.env.data(), nullptr);
    int index = 0;
    for (int l = 0; l < counts[shells[3]]; ++l)
      for (int k = 0; k < counts[shells[2]]; ++k)
        for (int j = 0; j < counts[shells[1]]; ++j)
          for (int i = 0; i < counts[shells[0]]; ++i, ++index)
            values[index] *= norm[input.basidx[2 * shells[0]] + i] *
                norm[input.basidx[2 * shells[1]] + j] *
                norm[input.basidx[2 * shells[2]] + k] *
                norm[input.basidx[2 * shells[3]] + l];
  }
  std::atomic<bool> valid(true);
#pragma omp parallel
  {
    xmvb::vb::LibcintShellWorkspace workspace;
#pragma omp for schedule(dynamic, 1)
    for (int task = 0; task < 81; ++task) {
      const auto& block = evaluator.evaluate_two_electron_shell_quartet(
          task / 27, (task / 9) % 3, (task / 3) % 3, task % 3, workspace);
      for (std::size_t i = 0; i < block.values.size(); ++i)
        if (std::abs(block.values[i] - references[task][i]) > 1e-12)
          valid.store(false);
    }
  }
  require(valid.load(), "shared libcint optimizer differs from legacy reference");
  xmvb::vb::LibcintShellWorkspace workspace;
  evaluator.evaluate_two_electron_shell_quartet(2, 2, 2, 2, workspace);
  const double* values = workspace.quartet.values.data();
  const double* cache = workspace.cache.data();
  for (int repeat = 0; repeat < 20; ++repeat) {
    evaluator.evaluate_two_electron_shell_quartet(2, 2, 2, 2, workspace);
    require(values == workspace.quartet.values.data() && cache == workspace.cache.data(),
            "same-size shell evaluation reallocates workspace");
  }
}

/** @brief Each canonical ERI has exactly one owner, including empty shards. */
void check_integral_shards() {
  auto input = make_input(3);
  const int counts[] = {1, 3, 6};
  int start = 0;
  for (int shell = 0; shell < 3; ++shell) {
    input.bas[shell * BAS_SLOTS + ANG_OF] = shell;
    input.env[input.bas[shell * BAS_SLOTS + PTR_COEFF]] =
        xmvb::vb::shell_normalization(shell, 1.0);
    input.basidx[2 * shell] = start;
    input.basidx[2 * shell + 1] = counts[shell];
    start += counts[shell];
  }
  xmvb::vb::LibcintMaterializedIntegralProvider provider;
  xmvb::vb::LibcintMaterializedIntegralProviderOptions full_options;
  full_options.use_historical_shell_pair_prescreen = false;
  const auto full = provider.build(input, full_options);
  using Key = std::pair<int, int>;
  std::map<Key, double> reference;
  for (std::size_t i = 0; i < full.two_electron_values.size(); ++i)
    require(reference.emplace(Key{full.left_pair_indices[i], full.right_pair_indices[i]},
                              full.two_electron_values[i]).second,
            "duplicate canonical integral in reference");
  // Independent complete Cartesian shell traversal: no provider ownership or
  // shell-order filtering. This catches shell-vs-AO canonical-order mistakes.
  xmvb::vb::LibcintDirectShellEvaluator evaluator(input);
  std::map<Key, double> direct;
  xmvb::vb::LibcintShellWorkspace workspace;
  for (int si = 0; si < 3; ++si)
    for (int sj = 0; sj < 3; ++sj)
      for (int sk = 0; sk < 3; ++sk)
        for (int sl = 0; sl < 3; ++sl) {
          const auto& quartet = evaluator.evaluate_two_electron_shell_quartet(si, sj, sk, sl, workspace);
          for (int li = 0; li < counts[si]; ++li)
            for (int lj = 0; lj < counts[sj]; ++lj)
              for (int lk = 0; lk < counts[sk]; ++lk)
                for (int ll = 0; ll < counts[sl]; ++ll) {
                  const int i = input.basidx[2 * si] + li;
                  const int j = input.basidx[2 * sj] + lj;
                  const int k = input.basidx[2 * sk] + lk;
                  const int l = input.basidx[2 * sl] + ll;
                  if (j > i || k > i || l > k || (i == k && l > j)) continue;
                  const auto index = li + counts[si] * (lj + counts[sj] * (lk + counts[sk] * ll));
                  const double value = quartet.values[index];
                  if (std::abs(value) < full_options.integral_tolerance) continue;
                  require(direct.emplace(Key{i * (i + 1) / 2 + j, k * (k + 1) / 2 + l}, value).second,
                          "duplicate integral in independent Cartesian traversal");
                }
        }
  require(direct == reference, "provider differs from complete Cartesian reference");
  for (int n_shards : {1, 2, 3, 16}) {
    std::map<Key, double> combined;
    for (int owner = 0; owner < n_shards; ++owner) {
      xmvb::vb::LibcintMaterializedIntegralProviderOptions options;
      options.shard_count = n_shards;
      options.shard_index = owner;
      options.include_core_hamiltonian = false;
      options.use_historical_shell_pair_prescreen = false;
      const auto shard = provider.build(input, options);
      require(shard.ao_core_hamiltonian_matrix.size() == 0,
              "worker retained core Hamiltonian");
      for (std::size_t i = 0; i < shard.two_electron_values.size(); ++i)
        require(combined.emplace(Key{shard.left_pair_indices[i], shard.right_pair_indices[i]},
                                 shard.two_electron_values[i]).second,
                "canonical integral has multiple shard owners");
    }
    require(combined == reference, "shards do not exactly reproduce all reference ERIs");
  }
  for (const auto bounds : {Key{-1, 2}, Key{2, 2}, Key{0, 0}}) {
    xmvb::vb::LibcintMaterializedIntegralProviderOptions options;
    options.shard_index = bounds.first;
    options.shard_count = bounds.second;
    bool rejected = false;
    try { provider.build(input, options); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "invalid shard ownership was accepted");
  }
  xmvb::vb::LibcintMaterializedIntegralProviderOptions options;
  options.shard_count = 2;
  bool rejected = false;
  try {
    xmvb::vb::build_materialized_ao_integral_input(provider.build(input, options));
  } catch (const std::invalid_argument&) { rejected = true; }
  require(rejected, "partial ERIs were accepted by the complete-input loader");
  options.shard_index = 1;
  require(provider.build(input, options).ao_core_hamiltonian_matrix.size() == 0,
          "non-root shard allocated a core Hamiltonian by default");
}

}  // namespace

int main() {
  try {
    check_shell_workspace();
    check_integral_shards();
    constexpr double pi = 3.14159265358979323846;
    using xmvb::vb::radial_cartesian_scale;
    require(radial_cartesian_scale(0,0,0) == 1.0, "s normalization");
    require(radial_cartesian_scale(1,0,0) == 1.0, "p normalization");
    require(std::abs(std::pow(radial_cartesian_scale(2,0,0),2) - 4*pi/5) < 1e-14,
            "Cartesian dxx radial norm");
    require(std::abs(std::pow(radial_cartesian_scale(1,1,0),2) - 4*pi/15) < 1e-14,
            "Cartesian dxy radial norm");
    require(std::abs(std::pow(radial_cartesian_scale(3,0,0),2) - 4*pi/7) < 1e-14,
            "Cartesian fxxx radial norm");

    const auto primary = make_input(1);
    xmvb::vb::LibcintRiIntegralProvider provider;
    const auto positive = provider.build(primary, make_input(1), {});
    const auto singular = provider.build(primary, make_input(3), {});
    require(positive.n_auxiliary_functions == 1, "positive metric rank");
    require(singular.n_auxiliary_functions == 1, "duplicate auxiliary shells must be removed");
    require(singular.auxiliary_metric_matrix.rows() == 3, "original metric size");
    require(singular.metric_whitened_ao_pair_factors.rows() == 1, "rectangular factor rank");
    const double a = positive.metric_whitened_ao_pair_factors.squaredNorm();
    const double b = singular.metric_whitened_ao_pair_factors.squaredNorm();
    require(std::abs(a-b) < 1e-12, "Cholesky and spectral RI contractions differ");
    std::cout << "RI metric normalization and rank regression passed\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
