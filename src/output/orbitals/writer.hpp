#pragma once

#include <filesystem>
#include "vbscf/orbitals/preparation/input.hpp"

namespace xmvb::output {

/** @brief Write converged orbitals beside the input as .orb, in GUS body format.
 * Preserves the stored AO indices and coefficients, including allowed zeros.
 * No GUS/END delimiters, normalization or orbital transformation is applied.
 */
std::filesystem::path write_orbital_file(
    const std::filesystem::path& input_path,
    const vb::OrbitalPreparationInput& orbitals);

}  // namespace xmvb::output
