#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <sys/resource.h>

#include <Eigen/Core>

#include "core/pair_index.hpp"
#include "mcscf/ci/sigma.hpp"
#include "mcscf/ci/product_support.hpp"
#include "vbscf/determinants/pairs/same_spin_cache.hpp"
#include "vbscf/structures/assembly/action.hpp"
#include "vbscf/structures/orthogonal_ci/exterior_transform.hpp"
#include "vbscf/structures/orthogonal_ci/integrals.hpp"
#include "../../mcscf/product_paths.hpp"

namespace {

using Strings = std::vector<std::vector<int>>;
using Pairing = std::vector<std::pair<int, int>>;
using Terms = std::vector<std::vector<xmvb::vb::StructureExpansionTerm>>;
using Clock = std::chrono::steady_clock;

struct CovalentExpansion {
  Terms terms;
  Strings product_alpha;
  Strings product_beta;
};

void require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

void append_combinations(int n, int k, int first, std::vector<int>* occupied,
                         Strings* strings) {
  if (static_cast<int>(occupied->size()) == k) {
    strings->push_back(*occupied);
    return;
  }
  const int remaining = k - static_cast<int>(occupied->size());
  for (int orbital = first; orbital <= n - remaining; ++orbital) {
    occupied->push_back(orbital);
    append_combinations(n, k, orbital + 1, occupied, strings);
    occupied->pop_back();
  }
}

Strings complete_strings(int n, int k) {
  Strings result;
  std::vector<int> occupied;
  append_combinations(n, k, 0, &occupied, &result);
  return result;
}

// A first site's mate splits a noncrossing matching into two even intervals.
std::vector<Pairing> rumer_pairings(int first, int past_last) {
  if (first == past_last) return {Pairing{}};
  std::vector<Pairing> result;
  for (int mate = first + 1; mate < past_last; mate += 2) {
    const auto inside = rumer_pairings(first + 1, mate);
    const auto outside = rumer_pairings(mate + 1, past_last);
    for (const auto& left : inside) {
      for (const auto& right : outside) {
        Pairing pairing{{first, mate}};
        pairing.insert(pairing.end(), left.begin(), left.end());
        pairing.insert(pairing.end(), right.begin(), right.end());
        result.push_back(std::move(pairing));
      }
    }
  }
  return result;
}

int occupation_mask(const std::vector<int>& occupied) {
  int mask = 0;
  for (int orbital : occupied) mask |= 1 << orbital;
  return mask;
}

// Each bond is (a_alpha^+ b_beta^+ - a_beta^+ b_alpha^+)/sqrt(2).
// Count the creator swaps needed to reach the alpha-then-beta determinant
// convention. This matters even when every site is singly occupied.
CovalentExpansion covalent_rumer_terms(
    int n, const Strings& alpha, const Strings& beta,
    const std::vector<Pairing>& pairings) {
  const int n_beta = static_cast<int>(beta.size());
  std::map<int, std::vector<xmvb::vb::StructureExpansionTerm>> sparse_terms;
  std::map<int, int> alpha_index;
  std::map<int, int> beta_index;
  for (int i = 0; i < static_cast<int>(alpha.size()); ++i)
    alpha_index.emplace(occupation_mask(alpha[i]), i);
  for (int i = 0; i < n_beta; ++i)
    beta_index.emplace(occupation_mask(beta[i]), i);
  const double scale = std::pow(2.0, -0.5 * (n / 2));
  for (int structure = 0; structure < static_cast<int>(pairings.size());
       ++structure) {
    const auto& pairing = pairings[structure];
    require(static_cast<int>(pairing.size()) == n / 2,
            "Rumer pairing has the wrong number of bonds");
    for (int choice = 0; choice < (1 << (n / 2)); ++choice) {
      std::vector<int> creators;
      int alpha_mask = 0;
      int beta_mask = 0;
      int sign = 1;
      for (int bond = 0; bond < n / 2; ++bond) {
        const auto [a, b] = pairing[bond];
        const bool reverse = (choice & (1 << bond)) != 0;
        const int alpha_site = reverse ? b : a;
        const int beta_site = reverse ? a : b;
        alpha_mask |= 1 << alpha_site;
        beta_mask |= 1 << beta_site;
        creators.push_back(a == alpha_site ? a : n + a);
        creators.push_back(b == alpha_site ? b : n + b);
        if (reverse) sign = -sign;
      }
      for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
          if (creators[i] > creators[j]) sign = -sign;
      const int product = alpha_index.at(alpha_mask) * n_beta +
                          beta_index.at(beta_mask);
      sparse_terms[product].push_back({structure, sign * scale});
    }
  }
  CovalentExpansion expansion;
  expansion.terms.reserve(sparse_terms.size());
  expansion.product_alpha.reserve(sparse_terms.size());
  expansion.product_beta.reserve(sparse_terms.size());
  for (auto& [product, product_terms] : sparse_terms) {
    expansion.product_alpha.push_back(alpha[product / n_beta]);
    expansion.product_beta.push_back(beta[product % n_beta]);
    expansion.terms.push_back(std::move(product_terms));
  }
  return expansion;
}

struct Fixture {
  int n = 0;
  Strings alpha;
  Strings beta;
  std::vector<double> overlap;
  Eigen::MatrixXd one_electron;
  xmvb::vb::ActiveSpaceTwoElectronResult two_electron;
};

Fixture make_fixture(int n, int n_alpha, int n_beta) {
  Fixture fixture;
  fixture.n = n;
  fixture.alpha = complete_strings(n, n_alpha);
  fixture.beta = complete_strings(n, n_beta);
  Eigen::MatrixXd upper = Eigen::MatrixXd::Identity(n, n);
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < j; ++i)
      upper(i, j) = 0.004 * std::cos(0.31 * (i + 2) * (j + 1));
  const Eigen::MatrixXd metric = upper.transpose() * upper;
  fixture.overlap.assign(metric.data(), metric.data() + metric.size());
  fixture.one_electron.resize(n, n);
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i <= j; ++i) {
      const double value = i == j ? -0.8 + 0.06 * i
                                  : 0.015 * std::sin(0.19 * (i + 1) * (j + 3));
      fixture.one_electron(i, j) = value;
      fixture.one_electron(j, i) = value;
    }
  }
  const int n_pairs = xmvb::vb::packed_active_pair_count(n);
  Eigen::MatrixXd factor(4, n_pairs);
  for (int j = 0; j < n_pairs; ++j)
    for (int i = 0; i < factor.rows(); ++i)
      factor(i, j) = 0.018 * std::cos(0.13 * (i + 2) * (j + 1));
  const Eigen::MatrixXd kernel = factor.transpose() * factor;
  auto& packed = fixture.two_electron.packed_active_two_electron_integrals;
  packed.resize(xmvb::vb::packed_active_two_electron_integral_count(n));
  for (int j = 0; j < n_pairs; ++j)
    for (int i = 0; i <= j; ++i)
      packed[xmvb::core::PairIndex::packed_pair_of_pairs_index(i, j)] =
          kernel(i, j);
  return fixture;
}

xmvb::vb::SameSpinPairCacheContext topology(
    const Fixture& fixture, const Strings& product_alpha,
    const Strings& product_beta) {
  return xmvb::vb::build_orthogonal_spin_topology(
      product_alpha, product_beta, fixture.n);
}

std::pair<Strings, Strings> full_products(const Fixture& fixture) {
  Strings product_alpha;
  Strings product_beta;
  const std::size_t products = fixture.alpha.size() * fixture.beta.size();
  product_alpha.reserve(products);
  product_beta.reserve(products);
  for (const auto& a : fixture.alpha)
    for (const auto& b : fixture.beta) {
      product_alpha.push_back(a);
      product_beta.push_back(b);
    }
  return {std::move(product_alpha), std::move(product_beta)};
}

Eigen::MatrixXd test_vectors(int n_structures, int width) {
  Eigen::MatrixXd vectors(n_structures, width);
  for (int column = 0; column < width; ++column)
    for (int row = 0; row < n_structures; ++row)
      vectors(row, column) =
          std::cos(0.17 * (row + 1) * (column + 2)) +
          0.2 * std::sin(0.11 * (row + 3) * (column + 1));
  return vectors;
}

// Reference independently expands C into the complete orthonormal carrier,
// uses public direct-CI sigma, then applies exterior adjoints and C^T.
xmvb::vb::StructureActionResult reference_images(
    const Fixture& fixture, const Terms& terms,
    const Strings& product_alpha, const Strings& product_beta,
    const Eigen::MatrixXd& vectors) {
  const auto integrals = xmvb::vb::orthogonalize_active_integrals(
      fixture.overlap, fixture.one_electron, fixture.two_electron, fixture.n);
  const xmvb::mcscf::DirectCiSigmaAction sigma(
      fixture.alpha, fixture.beta, integrals);
  const xmvb::vb::ExteriorOrbitalTransform alpha_transform(
      fixture.alpha, integrals.orbital_transform);
  const xmvb::vb::ExteriorOrbitalTransform beta_transform(
      fixture.beta, integrals.orbital_transform);
  std::map<int, int> alpha_index;
  std::map<int, int> beta_index;
  for (int i = 0; i < static_cast<int>(fixture.alpha.size()); ++i)
    alpha_index.emplace(occupation_mask(fixture.alpha[i]), i);
  for (int i = 0; i < static_cast<int>(fixture.beta.size()); ++i)
    beta_index.emplace(occupation_mask(fixture.beta[i]), i);
  const int n_beta = static_cast<int>(fixture.beta.size());
  Eigen::MatrixXd carrier = Eigen::MatrixXd::Zero(
      fixture.alpha.size(), vectors.cols() * n_beta);
  for (int product = 0; product < static_cast<int>(terms.size()); ++product) {
    const int a = alpha_index.at(occupation_mask(product_alpha[product]));
    const int b = beta_index.at(occupation_mask(product_beta[product]));
    for (const auto& term : terms[product])
      for (int column = 0; column < vectors.cols(); ++column)
        carrier(a, column * n_beta + b) +=
            term.coefficient * vectors(term.structure_index, column);
  }
  alpha_transform.apply_left(&carrier);
  for (int column = 0; column < vectors.cols(); ++column) {
    auto block = carrier.middleCols(column * n_beta, n_beta);
    beta_transform.apply_right_block(block);
  }
  Eigen::MatrixXd hamiltonian = sigma.apply(carrier);
  for (int column = 0; column < vectors.cols(); ++column) {
    auto h = hamiltonian.middleCols(column * n_beta, n_beta);
    auto s = carrier.middleCols(column * n_beta, n_beta);
    beta_transform.apply_adjoint_right_block(h);
    beta_transform.apply_adjoint_right_block(s);
  }
  alpha_transform.apply_adjoint_left(&hamiltonian);
  alpha_transform.apply_adjoint_left(&carrier);
  xmvb::vb::StructureActionResult result{
      Eigen::MatrixXd::Zero(vectors.rows(), vectors.cols()),
      Eigen::MatrixXd::Zero(vectors.rows(), vectors.cols())};
  for (int product = 0; product < static_cast<int>(terms.size()); ++product) {
    const int a = alpha_index.at(occupation_mask(product_alpha[product]));
    const int b = beta_index.at(occupation_mask(product_beta[product]));
    for (const auto& term : terms[product])
      for (int column = 0; column < vectors.cols(); ++column) {
        result.hamiltonian(term.structure_index, column) +=
            term.coefficient * hamiltonian(a, column * n_beta + b);
        result.overlap(term.structure_index, column) +=
            term.coefficient * carrier(a, column * n_beta + b);
      }
  }
  return result;
}

void check_images(const xmvb::vb::StructureActionResult& actual,
                  const xmvb::vb::StructureActionResult& expected) {
  const double h_error = (actual.hamiltonian - expected.hamiltonian).norm();
  const double s_error = (actual.overlap - expected.overlap).norm();
  require(h_error <= 2e-10 * std::max(1.0, expected.hamiltonian.norm()),
          "structure Hamiltonian differs from direct-CI reference");
  require(s_error <= 2e-11 * std::max(1.0, expected.overlap.norm()),
          "structure overlap differs from carrier reference");
}

void check_rumer_case(int n) {
  const Fixture fixture = make_fixture(n, n / 2, n / 2);
  const auto pairings = rumer_pairings(0, n);
  const CovalentExpansion expansion = covalent_rumer_terms(
      n, fixture.alpha, fixture.beta, pairings);
  const int expected = n == 4 ? 2 : 5;
  require(static_cast<int>(pairings.size()) == expected,
          "Rumer pairing count is incorrect");
  int term_count = 0;
  for (const auto& product : expansion.terms) term_count += product.size();
  require(term_count == expected * (1 << (n / 2)),
          "Rumer expansion is not complete");
  require(expansion.terms.size() == fixture.alpha.size(),
          "singly occupied Rumer products were not compacted");
  const auto cache = topology(
      fixture, expansion.product_alpha, expansion.product_beta);
  const xmvb::vb::StructureAction action(
      expansion.terms, expected, cache, fixture.overlap, fixture.one_electron,
      fixture.two_electron, n, nullptr,
      xmvb::vb::StructureActionMode::ActionOnly);
  require(action.storage().orthogonal_direct_ci,
          "complete spin strings failed to select direct CI");
  const Eigen::MatrixXd vectors = test_vectors(expected, 3);
  const auto actual = action.apply(vectors);
  // Force one-column panels and eviction/recomputation in both caches.
  // The resource choice must not change the operator or selected subspace.
  check_images(action.apply_streamed(vectors, 512), actual);
  check_images(actual, reference_images(
      fixture, expansion.terms, expansion.product_alpha,
      expansion.product_beta, vectors));
  require((actual.overlap.array() * vectors.array()).colwise().sum().minCoeff()
              > 0.0,
          "Rumer structure metric is not positive on test vectors");
}

void check_selected_unequal_spin_case() {
  const Fixture fixture = make_fixture(5, 2, 1);
  Strings product_alpha;
  Strings product_beta;
  for (int a = 0; a < 7; ++a)
    for (int b = 0; b < 4; ++b) {
      product_alpha.push_back(fixture.alpha[a]);
      product_beta.push_back(fixture.beta[b]);
    }
  const int products = static_cast<int>(product_alpha.size());
  constexpr int n_structures = 7;
  Terms terms(products);
  for (int product = 0; product < products; ++product)
    for (int structure = 0; structure < n_structures; ++structure)
      if ((product + 3 * structure) % 11 < 3)
        terms[product].push_back({structure,
            std::sin(0.23 * (product + 1) * (structure + 2))});
  const auto cache = topology(fixture, product_alpha, product_beta);
  const xmvb::vb::StructureAction action(
      terms, n_structures, cache, fixture.overlap, fixture.one_electron,
      fixture.two_electron, fixture.n, nullptr,
      xmvb::vb::StructureActionMode::ActionOnly);
  require(action.storage().orthogonal_direct_ci,
          "selected spin strings did not select completed direct CI");
  const Eigen::MatrixXd vectors = test_vectors(n_structures, 4);
  check_images(action.apply(vectors), reference_images(
      fixture, terms, product_alpha, product_beta, vectors));
  check_images(action.apply_streamed(vectors, 512), reference_images(
      fixture, terms, product_alpha, product_beta, vectors));
}

void check_complete_unequal_spin_case() {
  const Fixture fixture = make_fixture(5, 2, 1);
  const auto [product_alpha, product_beta] = full_products(fixture);
  constexpr int n_structures = 7;
  Terms terms(product_alpha.size());
  for (int product = 0; product < static_cast<int>(terms.size()); ++product)
    for (int structure = 0; structure < n_structures; ++structure)
      if ((product + 3 * structure) % 11 < 3)
        terms[product].push_back({structure,
            std::sin(0.23 * (product + 1) * (structure + 2))});
  const auto cache = topology(fixture, product_alpha, product_beta);
  const xmvb::vb::StructureAction action(
      terms, n_structures, cache, fixture.overlap, fixture.one_electron,
      fixture.two_electron, fixture.n, nullptr,
      xmvb::vb::StructureActionMode::ActionOnly);
  require(action.storage().orthogonal_direct_ci,
          "unequal complete carrier did not select direct CI");
  const Eigen::MatrixXd vectors = test_vectors(n_structures, 3);
  check_images(action.apply(vectors), reference_images(
      fixture, terms, product_alpha, product_beta, vectors));
}

void check_sigma_columns(const Fixture& fixture, int width) {
  const auto integrals = xmvb::vb::orthogonalize_active_integrals(
      fixture.overlap, fixture.one_electron, fixture.two_electron, fixture.n);
  const xmvb::mcscf::DirectCiSigmaAction sigma(
      fixture.alpha, fixture.beta, integrals);
  const xmvb::vb::ExteriorOrbitalTransform transform(
      fixture.alpha, integrals.orbital_transform);
  const auto rows = transform.row_block(1, 3);
  require(rows.get() == transform.row_block(1, 3).get(),
          "accepted-point exterior rows were regenerated");
  Eigen::MatrixXd selected = Eigen::MatrixXd::Zero(3, fixture.alpha.size());
  selected.middleCols(1, 3).setIdentity();
  transform.apply_adjoint_right_block(selected);
  require((*rows - selected).norm() < 1e-14,
          "cached exterior rows differ from the shear transform");
  const int n_beta = static_cast<int>(fixture.beta.size());
  Eigen::MatrixXd carrier(fixture.alpha.size(), width * n_beta);
  for (int j = 0; j < carrier.cols(); ++j)
    for (int i = 0; i < carrier.rows(); ++i)
      carrier(i, j) = std::sin(0.14 * (i + 2) * (j + 1));
  Eigen::MatrixXd tiled = Eigen::MatrixXd::Zero(carrier.rows(), carrier.cols());
  const auto read = [&](int beta, const std::vector<int>& alpha_rows, Eigen::Ref<Eigen::MatrixXd> block) {
        for (int vector = 0; vector < width; ++vector)
          for (int row = 0; row < static_cast<int>(alpha_rows.size()); ++row)
            block(row, vector) = carrier(alpha_rows[row], vector * n_beta + beta);
      };
  const auto add = [&](int beta, const std::vector<int>& alpha_rows, const Eigen::Ref<const Eigen::MatrixXd>& block) {
        for (int vector = 0; vector < width; ++vector)
          for (int row = 0; row < static_cast<int>(alpha_rows.size()); ++row)
            tiled(alpha_rows[row], vector * n_beta + beta) += block(row, vector);
      };
  sigma.apply_tiled(width, read, add);
  const Eigen::MatrixXd ordinary = sigma.apply(carrier);
  require((tiled - ordinary).norm() <=
              2e-11 * std::max(1.0, ordinary.norm()),
          "multi-column tiled direct-CI sigma differs from ordinary sigma");
  if (fixture.n == 5) {
    bool read_failed = false, write_failed = false;
    try {
      sigma.apply_tiled(width,
          [](int, const std::vector<int>&, Eigen::Ref<Eigen::MatrixXd>) {
            throw std::runtime_error("reader failure");
          }, add);
    } catch (const std::runtime_error& error) {
      read_failed = std::string(error.what()) == "reader failure";
    }
    try {
      sigma.apply_tiled(width, read,
          [](int, const std::vector<int>&, const Eigen::Ref<const Eigen::MatrixXd>&) {
            throw std::runtime_error("writer failure");
          });
    } catch (const std::runtime_error& error) {
      write_failed = std::string(error.what()) == "writer failure";
    }
    require(read_failed && write_failed,
            "tile callback exceptions did not propagate out of the OpenMP team");
  }
}

void check_prefix_words() {
  for (const int n : {4, 5, 12, 13, 16, 32, 63}) {
    const int na = n / 2, nb = n - na;
    const auto full = (std::uint64_t{1} << n) - 1;
    const auto low_a = (std::uint64_t{1} << na) - 1;
    const auto low_b = (std::uint64_t{1} << nb) - 1;
    const std::vector<std::uint64_t> alpha{low_a, full ^ low_b, full & 0xAAAAAAAAAAAAAAAAULL};
    const std::vector<std::uint64_t> beta{low_b, full ^ low_a, full ^ alpha.back()};
    const xmvb::mcscf::ProductSupport support(n, alpha, beta);
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b)
        for (int depth = 0; depth <= 1; ++depth) {
          int count = 0;
          bool expected = true;
          for (int k = 0; k < n; ++k) {
            count += ((alpha[a] >> k) & 1) + ((beta[b] >> k) & 1);
            expected = expected && count + depth >= k + 1;
          }
          require(support.contains(a, b, depth) == expected,
                  "packed prefix comparison differs from scalar occupation counts");
        }
  }
}

int string_index(const Strings& strings, const std::vector<int>& occupied) {
  const auto found = std::find(strings.begin(), strings.end(), occupied);
  require(found != strings.end(), "test spin string is absent");
  return static_cast<int>(found - strings.begin());
}

bool dominates_seed(const std::vector<int>& target,
                    const std::vector<int>& seed) {
  if (target.size() != seed.size()) return false;
  for (std::size_t i = 0; i < target.size(); ++i)
    if (target[i] > seed[i]) return false;
  return true;
}

void check_product_support_oracle(
    const Fixture& fixture, const std::vector<int>& products) {
  std::vector<std::uint64_t> alpha, beta;
  for (const auto& occupied : fixture.alpha) alpha.push_back(occupation_mask(occupied));
  for (const auto& occupied : fixture.beta) beta.push_back(occupation_mask(occupied));
  const int na = static_cast<int>(alpha.size());
  const int nb = static_cast<int>(beta.size());
  const xmvb::mcscf::ProductSupport support(fixture.n, alpha, beta, products);
  const auto direct = [&](int a, int b) {
    for (const int product : products)
      if (dominates_seed(fixture.alpha[a], fixture.alpha[product / nb]) &&
          dominates_seed(fixture.beta[b], fixture.beta[product % nb]))
        return true;
    return false;
  };
  std::vector<xmvb::test::ProductPaths::Seed> seeds;
  for (const int product : products)
    seeds.emplace_back(alpha[product / nb], beta[product % nb]);
  const xmvb::test::ProductPaths paths(fixture.n, fixture.alpha[0].size(),
                                      fixture.beta[0].size(), seeds);
  // Full masks exist only in this small independent oracle, not in paths.
  std::vector<unsigned char> generated(na * nb, 0), generated_halo(na * nb, 0);
  const auto check_paths = [&](int depth, std::vector<unsigned char>& found) {
    paths.visit(depth, [&](std::uint64_t a, std::uint64_t b, bool in_support) {
      const auto ai = std::find(alpha.begin(), alpha.end(), a) - alpha.begin();
      const auto bi = std::find(beta.begin(), beta.end(), b) - beta.begin();
      require(ai < na && bi < nb, "generated string is outside fixed-spin sector");
      require(found[ai * nb + bi]++ == 0, "duplicate product path");
      require(in_support == direct(ai, bi), "incorrect generated support label");
    });
  };
  check_paths(0, generated);
  check_paths(1, generated_halo);
  for (int a = 0; a < na; ++a)
    for (int b = 0; b < nb; ++b) {
      const bool expected = direct(a, b);
      require(support.contains(a, b) == expected,
              "paired product downset differs from sorted-occupation oracle");
      require(static_cast<bool>(generated[a * nb + b]) == expected,
              "generated paths differ from sorted-occupation oracle");
      bool neighbor = expected;
      for (int other = 0; other < na; ++other)
        if (__builtin_popcountll(alpha[a] ^ alpha[other]) == 2)
          neighbor = neighbor || direct(other, b);
      for (int other = 0; other < nb; ++other)
        if (__builtin_popcountll(beta[b] ^ beta[other]) == 2)
          neighbor = neighbor || direct(a, other);
      require(support.contains(a, b, 1) == neighbor,
              "paired product generator closure differs from explicit neighbors");
      require(static_cast<bool>(generated_halo[a * nb + b]) == neighbor,
              "generated halo differs from explicit one-electron neighbors");
    }
}

void check_generic_product_support() {
  const auto check = [](const Fixture& fixture,
                        const std::vector<std::pair<std::vector<int>, std::vector<int>>>& seeds) {
    std::vector<int> products;
    for (const auto& [alpha, beta] : seeds)
      products.push_back(string_index(fixture.alpha, alpha) * fixture.beta.size() +
                         string_index(fixture.beta, beta));
    check_product_support_oracle(fixture, products);
  };
  check(make_fixture(4, 2, 2), {{{1, 2}, {1, 2}}}); // Ionic seed.
  check(make_fixture(5, 2, 1),
        {{{0, 3}, {2}}, {{1, 4}, {0}}});
  const Fixture cross = make_fixture(6, 3, 3);
  const int a = string_index(cross.alpha, {0, 4, 5});
  const int b = string_index(cross.beta, {1, 2, 3});
  const int c = string_index(cross.alpha, {1, 2, 3});
  const int d = string_index(cross.beta, {0, 4, 5});
  check_product_support_oracle(cross,
      {a * static_cast<int>(cross.beta.size()) + b,
       c * static_cast<int>(cross.beta.size()) + d});
  std::vector<std::uint64_t> alpha, beta;
  for (const auto& occupied : cross.alpha) alpha.push_back(occupation_mask(occupied));
  for (const auto& occupied : cross.beta) beta.push_back(occupation_mask(occupied));
  const xmvb::mcscf::ProductSupport support(cross.n, alpha, beta,
      {a * static_cast<int>(cross.beta.size()) + b,
       c * static_cast<int>(cross.beta.size()) + d});
  require(!support.contains(a, d) && !support.contains(a, d, 1),
          "alpha and beta matched different support seeds");
  check_product_support_oracle(make_fixture(4, 2, 2), {});
  check_product_support_oracle(make_fixture(4, 2, 2), {0, 0, 35});
  check_product_support_oracle(make_fixture(4, 0, 2), {3});
  // No SpinStrings tables: a huge universe with tiny support is generated
  // directly. This is an ionic product test, not a 30-orbital covalent VB run.
  using Paths = xmvb::test::ProductPaths;
  const auto low = (std::uint64_t{1} << 15) - 1;
  const Paths sparse(30, 15, 15, {{low, low}});
  const auto sparse0 = sparse.visit(0, [](auto, auto, bool) {});
  const auto sparse1 = sparse.visit(1, [](auto, auto, bool) {});
  require(sparse0.products == 1 && sparse1.products == 451 && sparse1.support == 1,
          "large-universe sparse path count is wrong");
  require(sparse1.prefixes <= 31 * sparse1.products,
          "prefix traversal is not output bounded");
  const Paths bit63(63, 1, 0, {{std::uint64_t{1} << 62, 0}});
  require(bit63.visit(1, [](auto, auto, bool) {}).products == 63,
          "63-orbital path generation failed");
}

void check_supported_sigma(int n, int n_alpha, int n_beta_electrons, int width,
                           const std::vector<int>& products = {}) {
  const auto fixture = make_fixture(n, n_alpha, n_beta_electrons);
  const auto integrals = xmvb::vb::orthogonalize_active_integrals(
      fixture.overlap, fixture.one_electron, fixture.two_electron, n);
  const xmvb::mcscf::DirectCiSigmaAction sigma(fixture.alpha, fixture.beta, integrals);
  std::vector<std::uint64_t> alpha, beta;
  for (const auto& occupied : fixture.alpha) alpha.push_back(occupation_mask(occupied));
  for (const auto& occupied : fixture.beta) beta.push_back(occupation_mask(occupied));
  const xmvb::mcscf::ProductSupport support = products.empty()
      ? xmvb::mcscf::ProductSupport(n, alpha, beta)
      : xmvb::mcscf::ProductSupport(n, alpha, beta, products);
  const int na = static_cast<int>(alpha.size()), nb = static_cast<int>(beta.size());
  Eigen::MatrixXd carrier = Eigen::MatrixXd::Zero(na, width * nb);
  for (int a = 0; a < na; ++a)
    for (int b = 0; b < nb; ++b) {
      // Independently enumerate all one-electron neighbors to certify that
      // prefix deficit <=1 is the closure, not just an assumed input mask.
      bool reachable = support.contains(a, b);
      for (int other = 0; other < na; ++other)
        if (__builtin_popcountll(alpha[a] ^ alpha[other]) == 2)
          reachable = reachable || support.contains(other, b);
      for (int other = 0; other < nb; ++other)
        if (__builtin_popcountll(beta[b] ^ beta[other]) == 2)
          reachable = reachable || support.contains(a, other);
      require(reachable == support.contains(a, b, 1), "one-generator closure is incorrect");
      if (support.contains(a, b))
        for (int v = 0; v < width; ++v)
          carrier(a, v * nb + b) = std::sin(0.17 * (a + 1) * (b + v + 2));
    }
  Eigen::MatrixXd actual = Eigen::MatrixXd::Zero(na, width * nb);
  sigma.apply_tiled(width,
      [&](int b, const std::vector<int>& rows, Eigen::Ref<Eigen::MatrixXd> values) {
        for (int v = 0; v < width; ++v)
          for (int i = 0; i < static_cast<int>(rows.size()); ++i)
            values(i, v) = carrier(rows[i], v * nb + b);
      },
      [&](int b, const std::vector<int>& rows, const Eigen::Ref<const Eigen::MatrixXd>& values) {
        for (int v = 0; v < width; ++v)
          for (int i = 0; i < static_cast<int>(rows.size()); ++i) {
            require(support.contains(rows[i], b), "sigma wrote outside the output support");
            actual(rows[i], v * nb + b) += values(i, v);
          }
      }, &support);
  auto expected = sigma.apply(carrier);
  for (int a = 0; a < na; ++a)
    for (int b = 0; b < nb; ++b)
      if (!support.contains(a, b))
        for (int v = 0; v < width; ++v) expected(a, v * nb + b) = 0;
  require((actual - expected).norm() <= 2e-11 * std::max(1.0, expected.norm()),
          "support/closure sigma differs from independently projected full sigma");
  Eigen::MatrixXd generator = Eigen::MatrixXd::Ones(n, n);
  const auto one_body = sigma.apply_one_body_generator(carrier, generator);
  double outside_input = 0.0;
  for (int a = 0; a < na; ++a)
    for (int b = 0; b < nb; ++b)
      if (!support.contains(a, b)) outside_input += one_body(a, b) * one_body(a, b);
  require(outside_input > 1e-10, "fixture does not exercise intermediate-only support");

  // Independently contract the ordered E_pq E_rs Hamiltonian using only
  // generated I paths. Dense generator images are small-test oracles, not
  // proposed production storage. This checks paths that leave S and return.
  std::vector<xmvb::test::ProductPaths::Seed> seeds;
  if (products.empty()) {
    const std::uint64_t full = (std::uint64_t{1} << n) - 1;
    for (const auto a : alpha) seeds.emplace_back(a, full ^ a);
  } else {
    for (const int product : products)
      seeds.emplace_back(alpha[product / nb], beta[product % nb]);
  }
  const xmvb::test::ProductPaths paths(n, n_alpha, n_beta_electrons, seeds);
  std::vector<Eigen::MatrixXd> moments;
  for (int p = 0; p < n; ++p)
    for (int q = 0; q < n; ++q) {
      Eigen::MatrixXd e = Eigen::MatrixXd::Zero(n, n);
      e(p, q) = 1.0;
      moments.push_back(sigma.apply_one_body_generator(carrier, e));
    }
  const auto eri = [&](int p, int q, int r, int s) {
    return integrals.pair_kernel(xmvb::core::PairIndex::packed_pair_index(p, q),
                                xmvb::core::PairIndex::packed_pair_index(r, s));
  };
  Eigen::MatrixXd h = integrals.one_electron;
  for (int p = 0; p < n; ++p)
    for (int s = 0; s < n; ++s)
      for (int q = 0; q < n; ++q) h(p, s) -= 0.5 * eri(p, q, q, s);
  Eigen::VectorXd contracted = Eigen::VectorXd::Zero(width);
  Eigen::VectorXd halo_contribution = Eigen::VectorXd::Zero(width);
  paths.visit(1, [&](std::uint64_t a, std::uint64_t b, bool in_support) {
    const int ai = std::find(alpha.begin(), alpha.end(), a) - alpha.begin();
    const int bi = std::find(beta.begin(), beta.end(), b) - beta.begin();
    for (int v = 0; v < width; ++v) {
      const int col = v * nb + bi;
      for (int p = 0; p < n; ++p)
        for (int q = 0; q < n; ++q) {
          if (in_support)
            contracted[v] += carrier(ai, col) * h(p, q) * moments[p * n + q](ai, col);
          for (int r = 0; r < n; ++r)
            for (int s = 0; s < n; ++s) {
              const double term = 0.5 * eri(p, q, r, s) *
                  moments[q * n + p](ai, col) * moments[r * n + s](ai, col);
              contracted[v] += term;
              if (!in_support) halo_contribution[v] += term;
            }
        }
    }
  });
  require(halo_contribution.cwiseAbs().maxCoeff() > 1e-10,
          "fixture does not require leave-and-return Hamiltonian paths");
  double maximum_error = 0.0;
  for (int v = 0; v < width; ++v) {
    const double reference = carrier.middleCols(v * nb, nb)
        .cwiseProduct(expected.middleCols(v * nb, nb)).sum();
    require(std::abs(contracted[v] - reference) <= 3e-11 * std::max(1.0, std::abs(reference)),
            "generated two-generator contraction differs from direct-CI Hamiltonian");
    maximum_error = std::max(maximum_error, std::abs(contracted[v] - reference));
  }
  std::cout << std::setprecision(17) << "generated_contraction n=" << n
            << " na=" << n_alpha << " nb=" << n_beta_electrons
            << " absolute_error=" << maximum_error
            << " outside_support_contribution=" << halo_contribution.cwiseAbs().maxCoeff()
            << '\n';
}

void check_covalent_subset(int selected_count, bool add_ionic) {
  const auto fixture = make_fixture(8, 4, 4);
  const auto pairings = rumer_pairings(0, 8);
  std::vector<Pairing> selected;
  const int selected_indices[] = {0, 5, 13};
  for (int i = 0; i < selected_count; ++i)
    selected.push_back(pairings[selected_indices[i]]);
  auto expansion = covalent_rumer_terms(
      8, fixture.alpha, fixture.beta, selected);
  if (add_ionic) {
    expansion.product_alpha.push_back({0, 1, 2, 3});
    expansion.product_beta.push_back({0, 4, 5, 6});
    expansion.terms.push_back({{selected_count, 0.3}});
  }
  const int n_structures = selected_count + static_cast<int>(add_ionic);
  const auto cache = topology(
      fixture, expansion.product_alpha, expansion.product_beta);
  const xmvb::vb::StructureAction action(expansion.terms, n_structures, cache,
      fixture.overlap, fixture.one_electron, fixture.two_electron, 8, nullptr,
      xmvb::vb::StructureActionMode::ActionOnly);
  require(action.storage().orthogonal_direct_ci, "selected covalent fixture missed direct CI");
  const auto vectors = test_vectors(n_structures, 2);
  const auto expected = reference_images(fixture, expansion.terms,
      expansion.product_alpha, expansion.product_beta, vectors);
  check_images(action.apply(vectors), expected);
  check_images(action.apply_streamed(vectors, 512), expected);
  require((action.apply_overlap(vectors) - expected.overlap).norm() <=
              2e-11 * std::max(1.0, expected.overlap.norm()),
          "selected covalent overlap action differs from reference");
}

void benchmark_selected12(int selected_count, bool compare_reference, bool paths_only = false) {
  require(selected_count == 1 || selected_count == 3 ||
          selected_count == 20 || selected_count == 132,
          "selected12 count must be 1, 3, 20, or 132");
  const Fixture fixture = make_fixture(12, 6, 6);
  const auto pairings = rumer_pairings(0, 12);
  const std::vector<Pairing> selected(pairings.begin(),
                                      pairings.begin() + selected_count);
  const auto expansion = covalent_rumer_terms(
      12, fixture.alpha, fixture.beta, selected);
  std::size_t term_count = 0;
  for (const auto& terms : expansion.terms) term_count += terms.size();
  std::map<int, int> alpha_index, beta_index;
  std::vector<std::uint64_t> alpha_masks, beta_masks;
  for (int i = 0; i < static_cast<int>(fixture.alpha.size()); ++i) {
    const int mask = occupation_mask(fixture.alpha[i]);
    alpha_index.emplace(mask, i);
    alpha_masks.push_back(mask);
  }
  for (int i = 0; i < static_cast<int>(fixture.beta.size()); ++i) {
    const int mask = occupation_mask(fixture.beta[i]);
    beta_index.emplace(mask, i);
    beta_masks.push_back(mask);
  }
  std::vector<int> products;
  products.reserve(expansion.terms.size());
  for (int i = 0; i < static_cast<int>(expansion.terms.size()); ++i)
    products.push_back(alpha_index.at(occupation_mask(expansion.product_alpha[i])) *
        static_cast<int>(fixture.beta.size()) +
        beta_index.at(occupation_mask(expansion.product_beta[i])));
  std::vector<xmvb::test::ProductPaths::Seed> seeds;
  for (const int product : products)
    seeds.emplace_back(alpha_masks[product / fixture.beta.size()],
                       beta_masks[product % fixture.beta.size()]);
  const auto paths_begin = Clock::now();
  const xmvb::test::ProductPaths paths(12, 6, 6, seeds);
  const auto paths_built = Clock::now();
  const auto generated0 = paths.visit(0, [](auto, auto, bool) {});
  const auto paths_zero = Clock::now();
  const auto generated1 = paths.visit(1, [](auto, auto, bool) {});
  const auto paths_done = Clock::now();
  require(generated0.products == generated1.support,
          "generated halo lost the input support");
  const int index = selected_count == 1 ? 0 : selected_count == 3 ? 1 : selected_count == 20 ? 2 : 3;
  const std::size_t expected0[] = {136374, 157108, 192380, 226512};
  const std::size_t expected1[] = {523307, 559367, 603701, 608751};
  require(generated0.products == expected0[index] && generated1.products == expected1[index],
          "generated 12-orbital supports differ from the earlier exhaustive oracle");

  const auto seconds = [](Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration<double>(b - a).count();
  };
  std::cout << std::setprecision(17)
            << "case=CAS(12,12)_selected\n"
            << "structures=" << selected_count << '\n'
            << "products=" << expansion.terms.size() << '\n'
            << "terms=" << term_count << '\n'
            << "support0_products=" << generated0.products << '\n'
            << "support1_products=" << generated1.products << '\n'
            << "paths_build_seconds=" << seconds(paths_begin, paths_built) << '\n'
            << "paths_support_seconds=" << seconds(paths_built, paths_zero) << '\n'
            << "paths_halo_seconds=" << seconds(paths_zero, paths_done) << '\n'
            << "paths_support_prefixes=" << generated0.prefixes << '\n'
            << "paths_halo_prefixes=" << generated1.prefixes << '\n'
            << "paths_halo_word_updates=" << generated1.words << '\n'
            << "paths_metadata_bytes=" << paths.metadata_bytes() << '\n'
            << "paths_workspace_bytes=" << paths.workspace_bytes(1) << '\n';
  if (paths_only) return;

  const auto cache = topology(fixture, expansion.product_alpha, expansion.product_beta);
  const xmvb::vb::StructureAction action(
      expansion.terms, selected_count, cache, fixture.overlap,
      fixture.one_electron, fixture.two_electron, 12, nullptr,
      xmvb::vb::StructureActionMode::ActionOnly);
  require(action.storage().orthogonal_direct_ci,
          "selected12 did not retain direct CI");

  const Eigen::MatrixXd vectors = test_vectors(selected_count, 1);
  const auto start = Clock::now();
  const auto actual = action.apply(vectors);
  const auto cold = Clock::now();
  const auto warm = action.apply(vectors);
  const auto warmed = Clock::now();
  check_images(warm, actual);
  xmvb::vb::StructureActionResult reference;
  if (compare_reference) {
    reference = reference_images(fixture, expansion.terms,
        expansion.product_alpha, expansion.product_beta, vectors);
    check_images(actual, reference);
  }
  const double denominator = vectors.col(0).dot(actual.overlap.col(0));
  require(denominator > 0.0, "selected12 Rayleigh denominator is nonpositive");
  const double energy =
      vectors.col(0).dot(actual.hamiltonian.col(0)) / denominator;
  rusage usage{};
  require(getrusage(RUSAGE_SELF, &usage) == 0, "getrusage failed");
  std::cout << "cold_apply_seconds=" << seconds(start, cold) << '\n'
            << "warm_apply_seconds=" << seconds(cold, warmed) << '\n'
            << "max_rss_kib=" << usage.ru_maxrss << '\n'
            << "rayleigh_energy=" << energy << '\n'
            << "reference_checked=" << static_cast<int>(compare_reference) << '\n';
  if (compare_reference)
    std::cout << "hamiltonian_reference_error="
              << (actual.hamiltonian - reference.hamiltonian).norm() << '\n'
              << "overlap_reference_error="
              << (actual.overlap - reference.overlap).norm() << '\n';
}

void benchmark12(bool compare_reference) {
  const auto start = Clock::now();
  const Fixture fixture = make_fixture(12, 6, 6);
  const auto pairings = rumer_pairings(0, 12);
  require(pairings.size() == 132, "CAS(12,12) Rumer count is not 132");
  require(fixture.alpha.size() == 924 && fixture.beta.size() == 924,
          "CAS(12,12) complete spin-string count is not 924");
  const CovalentExpansion expansion = covalent_rumer_terms(
      12, fixture.alpha, fixture.beta, pairings);
  int term_count = 0;
  for (const auto& product : expansion.terms) term_count += product.size();
  require(term_count == 132 * 64, "CAS(12,12) expansion is not 64 terms per structure");
  require(expansion.terms.size() == 924,
          "CAS(12,12) covalent support is not 924 spin products");
  const auto cache = topology(
      fixture, expansion.product_alpha, expansion.product_beta);
  const auto prepared = Clock::now();
  const xmvb::vb::StructureAction action(
      expansion.terms, 132, cache, fixture.overlap, fixture.one_electron,
      fixture.two_electron, 12, nullptr,
      xmvb::vb::StructureActionMode::ActionOnly);
  require(action.storage().orthogonal_direct_ci,
          "CAS(12,12) did not retain direct CI");
  const auto constructed = Clock::now();
  const Eigen::MatrixXd vectors = test_vectors(132, 1);
  const auto actual = action.apply(vectors);
  const auto applied = Clock::now();
  const auto warm = action.apply(vectors);
  const auto warmed = Clock::now();
  check_images(warm, actual);
  xmvb::vb::StructureActionResult reference;
  if (compare_reference) {
    reference = reference_images(
        fixture, expansion.terms, expansion.product_alpha,
        expansion.product_beta, vectors);
    check_images(actual, reference);
  }
  const auto referenced = Clock::now();
  const double denominator = vectors.col(0).dot(actual.overlap.col(0));
  require(denominator > 0.0, "CAS(12,12) Rayleigh denominator is nonpositive");
  const double energy =
      vectors.col(0).dot(actual.hamiltonian.col(0)) / denominator;
  require(std::isfinite(energy), "CAS(12,12) Rayleigh energy is nonfinite");
  rusage usage{};
  std::vector<std::uint64_t> alpha_masks, beta_masks;
  for (const auto& a : fixture.alpha) alpha_masks.push_back(occupation_mask(a));
  for (const auto& b : fixture.beta) beta_masks.push_back(occupation_mask(b));
  const xmvb::mcscf::ProductSupport support(12, alpha_masks, beta_masks);
  std::size_t supported = 0, closure = 0;
  for (int a = 0; a < 924; ++a)
    for (int b = 0; b < 924; ++b) {
      supported += support.contains(a, b);
      closure += support.contains(a, b, 1);
    }
  require(supported == 226512 && closure == 608751,
          "CAS(12,12) support/closure count differs from the reflection formula");
  require(getrusage(RUSAGE_SELF, &usage) == 0, "getrusage failed");
  const auto seconds = [](Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration<double>(b - a).count();
  };
  std::cout << std::setprecision(17)
            << "case=CAS(12,12)_fixed_orbital_covalent_kernel\n"
            << "alpha_strings=924 beta_strings=924 structures=132 products=924 terms=8448\n"
            << "supported_products=" << supported << '\n'
            << "generator_closure_products=" << closure << '\n'
            << "prepare_seconds=" << seconds(start, prepared) << '\n'
            << "construct_seconds=" << seconds(prepared, constructed) << '\n'
            << "apply_seconds=" << seconds(constructed, applied) << '\n'
            << "warm_apply_seconds=" << seconds(applied, warmed) << '\n'
            << "reference_seconds="
            << (compare_reference ? seconds(warmed, referenced) : 0.0) << '\n'
            << "reference_checked=" << (compare_reference ? 1 : 0) << '\n'
            << "max_rss_kib=" << usage.ru_maxrss << '\n'
            << "rayleigh_energy=" << energy << '\n';
  if (compare_reference)
    std::cout << "hamiltonian_reference_error="
              << (actual.hamiltonian - reference.hamiltonian).norm() << '\n'
              << "overlap_reference_error="
              << (actual.overlap - reference.overlap).norm() << '\n';
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc == 2 && std::string(argv[1]) == "--benchmark12") {
      benchmark12(true);
    } else if (argc == 2 &&
               std::string(argv[1]) == "--benchmark12-no-reference") {
      benchmark12(false);
    } else if (argc == 3 && std::string(argv[1]) == "--selected12") {
      benchmark_selected12(std::stoi(argv[2]), false);
    } else if (argc == 3 && std::string(argv[1]) == "--selected12-parity") {
      benchmark_selected12(std::stoi(argv[2]), true);
    } else if (argc == 3 && std::string(argv[1]) == "--paths12") {
      benchmark_selected12(std::stoi(argv[2]), false, true);
    } else if (argc == 1) {
      check_rumer_case(4);
      check_rumer_case(6);
      check_selected_unequal_spin_case();
      check_complete_unequal_spin_case();
      check_sigma_columns(make_fixture(5, 2, 1), 3);
      check_sigma_columns(make_fixture(10, 5, 4), 2);
      check_prefix_words();
      check_supported_sigma(4, 2, 2, 2);
      check_supported_sigma(5, 3, 2, 3);
      check_supported_sigma(6, 3, 3, 1);
      check_generic_product_support();
      check_supported_sigma(5, 2, 1, 2, {3 * 5 + 2, 6 * 5 + 0});
      check_covalent_subset(1, false);
      check_covalent_subset(3, false);
      check_covalent_subset(3, true);
      std::cout << "structure tiles: PASS\n";
    } else {
      throw std::invalid_argument(
          "usage: test_structure_tiles [--benchmark12|--benchmark12-no-reference|"
          "--selected12 N|--selected12-parity N|--paths12 N]");
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "structure tiles: FAIL: " << error.what() << '\n';
    return 1;
  }
}
