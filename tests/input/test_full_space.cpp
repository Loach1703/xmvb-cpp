#include <algorithm>
#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "input/loading/generated_structures.hpp"
#include "vbscf/structures/analysis/spin.hpp"
#include "vbscf/structures/expansion/expander.hpp"

namespace {

struct Case {
  const char* name;
  int orbitals;
  int electrons;
  int multiplicity;
  int inactive;
  int structures;
  int determinants;
};

struct Sector {
  std::vector<int> structures;
  std::vector<int> determinants;
};

void require(bool condition, const std::string& message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

int choose(int n, int k) {
  if (k < 0 || k > n) return 0;
  k = std::min(k, n - k);
  int result = 1;
  for (int i = 1; i <= k; ++i) {
    result = result * (n - k + i) / i;
  }
  return result;
}

std::vector<int> raw_occupation(
    const xmvb::vb::RawStructureData& raw,
    int structure,
    int n_active_orbitals,
    int n_inactive) {
  std::vector<int> occupation(n_active_orbitals, 0);
  const int* orbital = raw.structure_orbitals_data(structure);
  for (int electron = 2 * n_inactive;
       electron < raw.n_total_electrons; ++electron) {
    const int active = orbital[electron] - n_inactive - 1;
    require(active >= 0 && active < n_active_orbitals,
            "raw structure contains an out-of-range active orbital");
    require(++occupation[active] <= 2,
            "raw structure contains more than two electrons in one orbital");
  }
  return occupation;
}

std::vector<int> determinant_occupation(
    const xmvb::vb::FullDeterminantStructureData& expansion,
    int determinant,
    int n_active_orbitals) {
  std::vector<int> occupation(n_active_orbitals, 0);
  for (const int orbital : expansion.alpha_det[determinant]) {
    require(orbital >= 0 && orbital < n_active_orbitals,
            "alpha string contains an out-of-range orbital");
    ++occupation[orbital];
  }
  for (const int orbital : expansion.beta_det[determinant]) {
    require(orbital >= 0 && orbital < n_active_orbitals,
            "beta string contains an out-of-range orbital");
    ++occupation[orbital];
  }
  require(std::all_of(occupation.begin(), occupation.end(),
                      [](int count) { return count <= 2; }),
          "determinant occupation exceeds two electrons");
  return occupation;
}

// Exact rank over a prime field. All current expansion coefficients are integers;
// a full modular rank certifies full real rank because every nonzero minor modulo
// the prime is also a nonzero integer minor.
int modular_rank(std::vector<std::vector<std::int64_t>> matrix) {
  constexpr std::int64_t prime = 1000000007;
  const int rows = static_cast<int>(matrix.size());
  const int columns = rows == 0 ? 0 : static_cast<int>(matrix.front().size());
  auto power = [](std::int64_t base, std::int64_t exponent) {
    std::int64_t result = 1;
    while (exponent > 0) {
      if (exponent & 1) result = result * base % prime;
      base = base * base % prime;
      exponent >>= 1;
    }
    return result;
  };
  int rank = 0;
  for (int column = 0; column < columns && rank < rows; ++column) {
    int pivot = rank;
    while (pivot < rows && matrix[pivot][column] == 0) ++pivot;
    if (pivot == rows) continue;
    std::swap(matrix[pivot], matrix[rank]);
    const std::int64_t inverse = power(matrix[rank][column], prime - 2);
    for (int row = rank + 1; row < rows; ++row) {
      if (matrix[row][column] == 0) continue;
      const std::int64_t factor = matrix[row][column] * inverse % prime;
      for (int j = column; j < columns; ++j) {
        matrix[row][j] =
            (matrix[row][j] - factor * matrix[rank][j] % prime + prime) % prime;
      }
    }
    ++rank;
  }
  return rank;
}

void check_case(const Case& test) {
  xmvb::vb::InputDeckMetadata metadata;
  metadata.structure_class_keyword = "FULL";
  metadata.declared_active_orbitals = test.orbitals;
  metadata.declared_active_electrons = test.electrons;
  metadata.declared_spin_multiplicity = test.multiplicity;
  const auto raw = xmvb::vb::build_generated_raw_structures_from_structure_class(
      metadata, test.electrons + 2 * test.inactive);
  require(raw.n_structures == test.structures,
          std::string(test.name) + ": wrong FULL structure count");
  const auto expansion =
      xmvb::vb::FullDeterminantStructureExpander().expand(raw);
  require(static_cast<int>(expansion.alpha_det.size()) == test.determinants,
          std::string(test.name) + ": wrong determinant count");
  require(expansion.alpha_det.size() == expansion.beta_det.size() &&
              expansion.alpha_det.size() ==
                  expansion.determinant_to_structure_terms.size(),
          std::string(test.name) + ": inconsistent expansion dimensions");

  const auto spin =
      xmvb::vb::analyze_structure_spin(expansion, test.multiplicity);
  require(spin.spin_adapted(),
          std::string(test.name) + ": a structure is not highest weight");

  std::map<std::vector<int>, Sector> sectors;
  for (int structure = 0; structure < raw.n_structures; ++structure) {
    sectors[raw_occupation(raw, structure, test.orbitals, test.inactive)]
        .structures.push_back(structure);
  }
  for (int determinant = 0; determinant < test.determinants; ++determinant) {
    const auto occupation =
        determinant_occupation(expansion, determinant, test.orbitals);
    const auto sector = sectors.find(occupation);
    require(sector != sectors.end(),
            std::string(test.name) + ": determinant has no structure sector");
    require(!expansion.determinant_to_structure_terms[determinant].empty(),
            std::string(test.name) + ": determinant has no expansion term");
    sector->second.determinants.push_back(determinant);
    for (const auto& term :
         expansion.determinant_to_structure_terms[determinant]) {
      require(term.structure_index >= 0 &&
                  term.structure_index < raw.n_structures,
              std::string(test.name) + ": invalid structure index");
      require(raw_occupation(raw, term.structure_index, test.orbitals,
                             test.inactive) == occupation,
              std::string(test.name) + ": term crosses occupation sectors");
    }
  }

  int certified_dimension = 0;
  for (const auto& [occupation, sector] : sectors) {
    const int n_open = static_cast<int>(std::count(
        occupation.begin(), occupation.end(), 1));
    const int spin_excess = test.multiplicity - 1;
    require(n_open >= spin_excess && (n_open - spin_excess) % 2 == 0,
            std::string(test.name) + ": invalid sector spin parity");
    const int n_alpha_open = (n_open + spin_excess) / 2;
    const int determinant_dimension = choose(n_open, n_alpha_open);
    const int spin_dimension = determinant_dimension -
        choose(n_open, n_alpha_open + 1);
    require(static_cast<int>(sector.determinants.size()) ==
                determinant_dimension,
            std::string(test.name) + ": incomplete determinant sector");
    require(static_cast<int>(sector.structures.size()) == spin_dimension,
            std::string(test.name) + ": wrong spin-adapted sector count");

    std::map<int, int> local_structure;
    for (int column = 0; column < spin_dimension; ++column) {
      local_structure.emplace(sector.structures[column], column);
    }
    std::vector<std::vector<std::int64_t>> incidence(
        determinant_dimension,
        std::vector<std::int64_t>(spin_dimension, 0));
    for (int row = 0; row < determinant_dimension; ++row) {
      const int determinant = sector.determinants[row];
      for (const auto& term :
           expansion.determinant_to_structure_terms[determinant]) {
        const auto column = local_structure.find(term.structure_index);
        require(column != local_structure.end(),
                std::string(test.name) + ": term crosses sector columns");
        require(term.coefficient == 1.0 || term.coefficient == -1.0,
                std::string(test.name) + ": incidence is not integral unit sign");
        constexpr std::int64_t prime = 1000000007;
        incidence[row][column->second] =
            (incidence[row][column->second] +
             static_cast<std::int64_t>(term.coefficient) + prime) % prime;
      }
    }
    require(modular_rank(std::move(incidence)) == spin_dimension,
            std::string(test.name) + ": deficient spin-adapted sector rank");
    certified_dimension += spin_dimension;
  }
  require(certified_dimension == test.structures,
          std::string(test.name) + ": certified dimensions do not sum to CAS");
  std::cout << test.name << ": " << sectors.size() << " sectors, "
            << certified_dimension << " independent spin-adapted structures\n";
}

}  // namespace

int main() {
  check_case({"FeCl2 CAS(10,8) quintet", 8, 10, 5, 0, 420, 448});
  check_case({"pyrrole CAS(11,10) doublet", 10, 11, 2, 0, 27720, 52920});
  check_case({"singlet CAS(4,4)", 4, 4, 1, 0, 20, 36});
  check_case({"singlet CAS(2,3)", 3, 2, 1, 0, 6, 9});
  check_case({"all-alpha holes CAS(3,5)", 5, 3, 4, 2, 10, 10});
}
