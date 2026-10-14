#include "output/orbitals/writer.hpp"

#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <stdexcept>

namespace xmvb::output {

std::filesystem::path write_orbital_file(
    const std::filesystem::path& input_path,
    const vb::OrbitalPreparationInput& orbitals) {
  auto path = std::filesystem::absolute(input_path);
  path.replace_extension(".orb");
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot open orbital output: " + path.string());
  out.imbue(std::locale::classic());
  for (std::size_t i = 0; i < orbitals.n_orbitals; ++i) {
    out << std::setw(6) << vb::stored_sparse_orbital_coefficient_count(orbitals, i);
    if ((i + 1) % 20 == 0 || i + 1 == orbitals.n_orbitals) out << '\n';
  }
  // Scientific notation retains tiny coefficients and always has a decimal
  // point, which the GUS reader uses to identify the first coefficient line.
  out << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10 - 1);
  for (std::size_t i = 0; i < orbitals.n_orbitals; ++i) {
    const int count = vb::stored_sparse_orbital_coefficient_count(orbitals, i);
    out << "# ORBITAL " << i + 1 << " NAO = " << count << '\n';
    for (int j = 0; j < count; ++j) {
      const auto offset = i * orbitals.n_basis_functions + j;
      out << std::setw(25) << orbitals.orbital_value_table[offset]
          << ' ' << std::setw(6) << orbitals.orbital_basis_index_table[offset];
      if ((j + 1) % 4 == 0 || j + 1 == count) out << '\n';
      else out << ' ';
    }
  }
  out.close();
  if (!out) throw std::runtime_error("cannot write orbital output: " + path.string());
  return path;
}

}  // namespace xmvb::output
