#include "vbscf/optimization/preconditioners/hessian_diagonal.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>

#include "core/openmp.hpp"
#include "vbscf/core/contracts/input.hpp"
#include "vbscf/derivatives/hessian/context/accepted_point.hpp"
#include "vbscf/integrals/active/two_electron/construction/indexer.hpp"
#include "vbscf/integrals/active/two_electron/response/types.hpp"
#include "vbscf/integrals/active/two_electron/transformation/ao_pair_operator.hpp"
#include "vbscf/integrals/active/two_electron/transformation/packed_pair_map.hpp"
#include "vbscf/integrals/ao/ri/factorization.hpp"
#include "vbscf/integrals/ao/one_electron/direct_operator.hpp"
#include "vbscf/integrals/ao/one_electron/ri_operator.hpp"

namespace xmvb::vb {
namespace {

std::size_t tensor4_index(
    int p, int q, int r, int s, int dimension) noexcept {
  return (((static_cast<std::size_t>(p) * dimension + q) * dimension + r) *
          dimension + s);
}

std::size_t slice_index(int p, int a, int b, int n_active) noexcept {
  return ((static_cast<std::size_t>(p) * n_active + a) * n_active + b);
}

std::size_t packed_eri_count(int n_active) noexcept {
  const std::size_t n_pairs =
      static_cast<std::size_t>(n_active) * (n_active + 1) / 2;
  return n_pairs * (n_pairs + 1) / 2;
}

void validate(
    const AnalyticDiagonalIntermediates& x,
    const Eigen::Ref<const Eigen::MatrixXd>& dm1,
    const std::vector<double>& packed_gradient) {
  const int n_mo = static_cast<int>(x.h_core.rows());
  if (x.n_core < 0 || x.n_active <= 0 || n_mo < x.n_core + x.n_active ||
      x.h_core.cols() != n_mo || x.v_core.rows() != n_mo ||
      x.v_core.cols() != n_mo || x.j_pc.rows() != n_mo ||
      x.j_pc.cols() != x.n_core || x.k_pc.rows() != n_mo ||
      x.k_pc.cols() != x.n_core || dm1.rows() != x.n_active ||
      dm1.cols() != x.n_active ||
      x.j_active.size() != static_cast<std::size_t>(n_mo) *
                             x.n_active * x.n_active ||
      x.k_active.size() != x.j_active.size() ||
      x.active_eri.size() != static_cast<std::size_t>(x.n_active) *
                             x.n_active * x.n_active * x.n_active ||
      packed_gradient.size() != packed_eri_count(x.n_active) ||
      !x.h_core.allFinite() || !x.v_core.allFinite() ||
      !x.j_pc.allFinite() || !x.k_pc.allFinite() || !dm1.allFinite()) {
    throw std::invalid_argument("invalid VBSCF analytic Hessian-diagonal input");
  }
  if (!std::all_of(x.j_active.begin(), x.j_active.end(), [](double value) {
        return std::isfinite(value);
      }) ||
      !std::all_of(x.k_active.begin(), x.k_active.end(), [](double value) {
        return std::isfinite(value);
      }) ||
      !std::all_of(x.active_eri.begin(), x.active_eri.end(), [](double value) {
        return std::isfinite(value);
      }) ||
      !std::all_of(packed_gradient.begin(), packed_gradient.end(), [](double value) {
        return std::isfinite(value);
      })) {
    throw std::invalid_argument("non-finite VBSCF analytic Hessian-diagonal tensor");
  }
}

/** @brief Diagonal MO image of G[C W C^T] for symmetric source weights.
 *
 * Since G[D]=2J[D]-K[D], a single occupied-density contraction provides
 * every external-orbital exchange diagonal. No external-by-occupied pair
 * transform or full four-index MO tensor is required. The exact integral
 * sweeps use contracted density channels rather than external-by-occupied
 * mixed-pair channels. Only AO-sized worker outputs
 * are used, with their aggregate storage bounded before choosing a team.
 */
Eigen::VectorXd density_operator_diagonal(
    const VbScfInput& input, const Eigen::MatrixXd& mo,
    const Eigen::Ref<const Eigen::MatrixXd>& coefficients,
    const Eigen::Ref<const Eigen::MatrixXd>& weights) {
  const int n_bf = static_cast<int>(mo.rows());
  Eigen::MatrixXd image;
  if (input.standard_two_electron_mode ==
      StandardTwoElectronMode::ResolutionOfIdentity) {
    if (!input.ri_factorization) {
      throw std::invalid_argument("RI diagonal requires accepted AO factors");
    }
    AoEffectiveOneElectronRiLowRankFactors factors;
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> spectrum(weights);
    if (spectrum.info() != Eigen::Success || !spectrum.eigenvalues().allFinite()) {
      throw std::runtime_error("failed to factor active diagonal density source");
    }
    factors.scaled_factor_matrix.resize(n_bf,weights.rows());
    int positive = 0;
    for (Eigen::Index j = 0; j < weights.rows(); ++j)
      if (spectrum.eigenvalues()[j] >= 0.0) ++positive;
    factors.n_positive_components = positive;
    factors.n_negative_components = weights.rows()-positive;
    int positive_column = 0, negative_column = positive;
    for (Eigen::Index j = 0; j < weights.rows(); ++j) {
      const double eigenvalue = spectrum.eigenvalues()[j];
      const int column = eigenvalue >= 0.0 ? positive_column++ : negative_column++;
      factors.scaled_factor_matrix.col(column).noalias() =
          std::sqrt(std::abs(eigenvalue)) * coefficients * spectrum.eigenvectors().col(j);
    }
    const auto storage = apply_ao_effective_one_electron_ri_operator(
        factors, *input.ri_factorization, n_bf);
    image = Eigen::Map<const Eigen::MatrixXd>(storage.data(), n_bf, n_bf);
  } else {
    const Eigen::MatrixXd density = coefficients * weights * coefficients.transpose();
    constexpr std::size_t worker_budget = 64 * 1024 * 1024;
    const std::size_t image_bytes =
        sizeof(double) * static_cast<std::size_t>(n_bf) * n_bf;
    const int threads = std::min(
        effective_openmp_thread_count(),
        static_cast<int>(std::max<std::size_t>(1, worker_budget / image_bytes)));
    auto storage = apply_ao_h1e(density.data(), input.ao_integral_input, threads);
    Eigen::Map<Eigen::MatrixXd> triangular(storage.data(), n_bf, n_bf);
    for (int mu = 0; mu < n_bf; ++mu)
      for (int nu = 0; nu <= mu; ++nu) {
        triangular(mu, nu) += triangular(nu, mu);
        triangular(nu, mu) = triangular(mu, nu);
      }
    image = triangular;
  }
  const Eigen::MatrixXd products = image * mo;
  return (mo.array() * products.array()).colwise().sum().transpose();
}

struct DiagonalContractions {
  /** @brief Total effective-potential diagonal, including inactive reference. */
  Eigen::VectorXd potential;
  /** @brief Accepted orbital-gradient diagonal in the orthonormal frame. */
  Eigen::VectorXd gradient;
  /** @brief Contracted one-RDM Coulomb/exchange correction for occupied rows. */
  Eigen::MatrixXd one_rdm;
  /** @brief Contracted two-RDM correction for all frame rows. */
  Eigen::MatrixXd two_rdm;
};

/** @brief Assemble the same signed diagonal from already contracted sources. */
Eigen::MatrixXd assemble_rotation_diagonal(
    const AnalyticDiagonalIntermediates& x,
    const Eigen::Ref<const Eigen::MatrixXd>& dm1_active,
    const DiagonalContractions& c) {
  const int n_mo = x.h_core.rows();
  const int n_core = x.n_core;
  const int n_active = x.n_active;
  const int n_occ = n_core+n_active;
  Eigen::MatrixXd dm1 = Eigen::MatrixXd::Zero(n_mo,n_mo);
  dm1.diagonal().head(n_core).setConstant(2.0);
  dm1.block(n_core,n_core,n_active,n_active) = dm1_active;
  Eigen::MatrixXd diagonal = Eigen::MatrixXd::Zero(n_mo,n_mo);
  for (int p = 0; p < n_mo; ++p) {
    for (int q = 0; q < n_mo; ++q) {
      diagonal(p,q) = x.h_core(p,p)*dm1(q,q)-x.h_core(p,q)*dm1(p,q)+
          x.h_core(q,q)*dm1(p,p)-x.h_core(q,p)*dm1(q,p)-c.gradient[p]-c.gradient[q];
    }
    diagonal(p,p) += 2.0*c.gradient[p];
  }
  for (int p = 0; p < n_mo; ++p)
    for (int core = 0; core < n_core; ++core) {
      diagonal(p,core) += 2.0*c.potential[p];
      diagonal(core,p) += 2.0*c.potential[p];
    }
  for (int core = 0; core < n_core; ++core)
    diagonal(core,core) -= 4.0*c.potential[core];
  for (int p = 0; p < n_mo; ++p)
    for (int a = 0; a < n_active; ++a) {
      const double value = x.v_core(p,p)*dm1_active(a,a);
      diagonal(p,n_core+a) += value;
      diagonal(n_core+a,p) += value;
    }
  for (int a = 0; a < n_active; ++a)
    for (int b = 0; b < n_active; ++b) {
      const double value = -x.v_core(n_core+a,n_core+b)*dm1_active(a,b);
      diagonal(n_core+a,n_core+b) += value;
      diagonal(n_core+b,n_core+a) += value;
    }
  for (int p = n_core; p < n_mo; ++p)
    for (int core = 0; core < n_core; ++core) {
      const double value = 6.0*x.k_pc(p,core)-2.0*x.j_pc(p,core);
      diagonal(p,core) += value;
      diagonal(core,p) += value;
    }
  for (int p = 0; p < n_occ; ++p)
    for (int a = 0; a < n_active; ++a) {
      diagonal(p,n_core+a) -= c.one_rdm(p,a);
      diagonal(n_core+a,p) -= c.one_rdm(p,a);
    }
  for (int p = 0; p < n_mo; ++p)
    for (int a = 0; a < n_active; ++a) {
      diagonal(n_core+a,p) += c.two_rdm(p,a);
      diagonal(p,n_core+a) += c.two_rdm(p,a);
    }
  diagonal *= 2.0;
  diagonal = (0.5*(diagonal+diagonal.transpose())).eval();
  if (!diagonal.allFinite())
    throw std::runtime_error("non-finite VBSCF orbital Hessian diagonal");
  return diagonal;
}

/** @brief Evaluate only the density span consumed by the analytic diagonal.
 *
 * Untruncated Householder QR reduces the number of integral sweeps to
 * min(n_active_pairs,2*n_active+1). There is no numerical rank cutoff,
 * division by a small singular value, or approximate density compression.
 */
Eigen::MatrixXd contracted_density_images(
    const VbScfInput& input, const Eigen::MatrixXd& mo, int n_core,
    const Eigen::MatrixXd& dm1, const std::vector<double>& dm2) {
  const int n = dm1.rows();
  const int n_pairs = n*(n+1)/2;
  const int n_sources = 2*n+1;
  Eigen::MatrixXd sources = Eigen::MatrixXd::Zero(n_pairs,n_sources);
  for (int w = 0; w < n; ++w)
    for (int z = 0; z <= w; ++z) {
      const int pair = TwoElectronIndexer::packed_pair_index(w,z);
      sources(pair,0) = w == z ? dm1(w,w) : dm1(w,z)+dm1(z,w);
      for (int a = 0; a < n; ++a) {
        sources(pair,1+a) = w == z
            ? (w == a ? dm1(a,w) : 0.0)
            : (w == a ? dm1(a,z) : 0.0)+(z == a ? dm1(a,w) : 0.0);
        const double wz = dm2[tensor4_index(a,w,z,a,n)]+
                          dm2[tensor4_index(w,a,z,a,n)];
        const double zw = dm2[tensor4_index(a,z,w,a,n)]+
                          dm2[tensor4_index(z,a,w,a,n)];
        sources(pair,1+n+a) = w == z ? wz : wz+zw;
      }
    }
  const int width = std::min(n_pairs,n_sources);
  const Eigen::HouseholderQR<Eigen::MatrixXd> qr(sources);
  const Eigen::MatrixXd basis = qr.householderQ()*Eigen::MatrixXd::Identity(n_pairs,width);
  const Eigen::MatrixXd reconstruction =
      qr.matrixQR().topRows(width).template triangularView<Eigen::Upper>();
  Eigen::MatrixXd images(mo.cols(),width);
  const auto active = mo.middleCols(n_core,n);
  for (int j = 0; j < width; ++j) {
    Eigen::MatrixXd weights(n,n);
    for (int w = 0; w < n; ++w)
      for (int z = 0; z <= w; ++z) {
        const int pair = TwoElectronIndexer::packed_pair_index(w,z);
        weights(w,z) = weights(z,w) = basis(pair,j)*(w == z ? 1.0 : 0.5);
      }
    images.col(j) = density_operator_diagonal(input,mo,active,weights);
  }
  return images*reconstruction;
}

Eigen::VectorXd diagonal_gradient(
    const AnalyticDiagonalIntermediates& x,
    const Eigen::Ref<const Eigen::MatrixXd>& dm1,
    const std::vector<double>& dm2, const Eigen::VectorXd& potential) {
  const int n = x.n_active;
  const int core = x.n_core;
  Eigen::VectorXd gradient = Eigen::VectorXd::Zero(x.h_core.rows());
  gradient.head(core) = 2.0*(x.h_core.diagonal().head(core)+potential.head(core));
  for (int a = 0; a < n; ++a)
    for (int u = 0; u < n; ++u) {
      gradient[core+a] += (x.h_core(core+a,core+u)+x.v_core(core+a,core+u))*dm1(u,a);
      for (int w = 0; w < n; ++w)
        for (int z = 0; z < n; ++z)
          gradient[core+a] += x.active_eri[tensor4_index(a,u,w,z,n)]*
              dm2[tensor4_index(w,z,u,a,n)];
    }
  return gradient;
}

Eigen::MatrixXd build_pair_transform(
    const Eigen::Ref<const Eigen::MatrixXd>& orbital_transform) {
  const int n = static_cast<int>(orbital_transform.rows());
  if (n <= 0 || orbital_transform.cols() != n) {
    throw std::invalid_argument("orbital pair transform must be square");
  }
  const int n_pairs = n * (n + 1) / 2;
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(n_pairs, n_pairs);
  for (int new_p = 0; new_p < n; ++new_p) {
    for (int new_q = 0; new_q <= new_p; ++new_q) {
      const int new_pair = TwoElectronIndexer::packed_pair_index(new_p, new_q);
      for (int old_p = 0; old_p < n; ++old_p) {
        for (int old_q = 0; old_q <= old_p; ++old_q) {
          const int old_pair = TwoElectronIndexer::packed_pair_index(old_p, old_q);
          double value =
              orbital_transform(old_p, new_p) *
              orbital_transform(old_q, new_q);
          if (old_p != old_q) {
            value +=
                orbital_transform(old_q, new_p) *
                orbital_transform(old_p, new_q);
          }
          result(old_pair, new_pair) = value;
        }
      }
    }
  }
  return result;
}

ExactCtxPairMatrix apply_ao_pair_operator(
    const VbScfInput& input,
    const Eigen::Ref<const ExactCtxPairMatrix>& coefficients) {
  ExactCtxPairMatrix result;
  const int n_bf = input.orbital_preparation_input.n_basis_functions;
  if (input.standard_two_electron_mode ==
      StandardTwoElectronMode::ResolutionOfIdentity) {
    if (!input.ri_factorization) {
      throw std::invalid_argument("RI diagonal requires accepted AO factors");
    }
    const auto& factors = input.ri_factorization->metric_whitened_ao_pair_factors;
    result.noalias() = factors.transpose() * (factors * coefficients);
  } else {
    detail::apply_exact_ao_pair_kernel(
        input.ao_integral_input,
        coefficients,
        n_bf,
        static_cast<std::size_t>(coefficients.cols()),
        &result);
  }
  return result;
}

std::vector<double> transform_dm2(
    const std::vector<double>& source,
    const Eigen::Ref<const Eigen::MatrixXd>& transform) {
  const int n = static_cast<int>(transform.rows());
  std::vector<double> current = source;
  std::vector<double> next(current.size());
  for (int axis = 0; axis < 4; ++axis) {
    std::fill(next.begin(), next.end(), 0.0);
    for (int p = 0; p < n; ++p)
      for (int q = 0; q < n; ++q)
        for (int r = 0; r < n; ++r)
          for (int s = 0; s < n; ++s) {
            const int output[4] = {p, q, r, s};
            double value = 0.0;
            for (int old = 0; old < n; ++old) {
              int input_index[4] = {p, q, r, s};
              input_index[axis] = old;
              value += transform(output[axis], old) * current[tensor4_index(
                  input_index[0], input_index[1], input_index[2],
                  input_index[3], n)];
            }
            next[tensor4_index(p, q, r, s, n)] = value;
          }
    current.swap(next);
  }
  return current;
}

std::vector<double> pack_symmetric_dm2(
    const std::vector<double>& dm2,
    int n) {
  std::vector<double> packed(packed_eri_count(n), 0.0);
  for (int p = 0; p < n; ++p) {
    for (int q = 0; q <= p; ++q) {
      const int pq = TwoElectronIndexer::packed_pair_index(p, q);
      for (int r = 0; r < n; ++r) {
        for (int s = 0; s <= r; ++s) {
          const int rs = TwoElectronIndexer::packed_pair_index(r, s);
          if (rs > pq) continue;
          const int multiplicity =
              (p == q ? 1 : 2) * (r == s ? 1 : 2) * (pq == rs ? 1 : 2);
          packed[TwoElectronIndexer::packed_pair_of_pairs_index(pq, rs)] =
              0.5 * multiplicity * dm2[tensor4_index(p, q, r, s, n)];
        }
      }
    }
  }
  return packed;
}

}  // namespace

std::vector<double> unpack_symmetric_active_two_rdm(
    const std::vector<double>& packed_eri_gradient,
    int n_active) {
  if (n_active <= 0 || packed_eri_gradient.size() != packed_eri_count(n_active)) {
    throw std::invalid_argument("invalid packed active ERI gradient");
  }
  std::vector<double> dm2(
      static_cast<std::size_t>(n_active) * n_active * n_active * n_active,
      0.0);
  for (int p = 0; p < n_active; ++p) {
    for (int q = 0; q < n_active; ++q) {
      const int pq = TwoElectronIndexer::packed_pair_index(p, q);
      for (int r = 0; r < n_active; ++r) {
        for (int s = 0; s < n_active; ++s) {
          const int rs = TwoElectronIndexer::packed_pair_index(r, s);
          const int packed =
              TwoElectronIndexer::packed_pair_of_pairs_index(pq, rs);
          const int multiplicity =
              (p == q ? 1 : 2) * (r == s ? 1 : 2) * (pq == rs ? 1 : 2);
          dm2[tensor4_index(p, q, r, s, n_active)] =
              2.0 * packed_eri_gradient[packed] / multiplicity;
        }
      }
    }
  }
  return dm2;
}

Eigen::MatrixXd build_analytic_rotation_hessian_diagonal(
    const AnalyticDiagonalIntermediates& x,
    const Eigen::Ref<const Eigen::MatrixXd>& dm1_active,
    const std::vector<double>& packed_gradient) {
  validate(x, dm1_active, packed_gradient);
  const int n_mo = static_cast<int>(x.h_core.rows());
  const int n_core = x.n_core;
  const int n_active = x.n_active;
  const int n_occ = n_core + n_active;
  const std::vector<double> dm2 =
      unpack_symmetric_active_two_rdm(packed_gradient, n_active);

  Eigen::VectorXd v_total_diagonal = x.v_core.diagonal();
  for (int p = 0; p < n_mo; ++p) {
    double value = 0.0;
    for (int a = 0; a < n_active; ++a) {
      for (int b = 0; b < n_active; ++b) {
        value +=
            x.j_active[slice_index(p, a, b, n_active)] * dm1_active(a, b) -
            0.5 * x.k_active[slice_index(p, a, b, n_active)] * dm1_active(a, b);
      }
    }
    v_total_diagonal[p] += value;
  }
  Eigen::VectorXd gradient_diagonal =
      diagonal_gradient(x,dm1_active,dm2,v_total_diagonal);

  Eigen::MatrixXd jkcaa = Eigen::MatrixXd::Zero(n_occ, n_active);
  Eigen::MatrixXd v_dm2 = Eigen::MatrixXd::Zero(n_mo, n_active);
  for (int p = 0; p < n_mo; ++p) {
    for (int a = 0; a < n_active; ++a) {
      double diagonal_dm2 = 0.0;

      for (int w = 0; w < n_active; ++w) {
        for (int z = 0; z < n_active; ++z) {
          diagonal_dm2 +=
              x.j_active[slice_index(p, w, z, n_active)] *
              dm2[tensor4_index(w, z, a, a, n_active)];
          diagonal_dm2 +=
              x.k_active[slice_index(p, w, z, n_active)] *
              (dm2[tensor4_index(a, w, z, a, n_active)] +
               dm2[tensor4_index(w, a, z, a, n_active)]);
        }
      }
      v_dm2(p, a) = diagonal_dm2;
      if (p < n_occ) {
        double value = 0.0;
        for (int v = 0; v < n_active; ++v) {
          value +=
              (6.0 * x.k_active[slice_index(p, a, v, n_active)] -
               2.0 * x.j_active[slice_index(p, a, v, n_active)]) *
              dm1_active(a, v);
        }
        jkcaa(p, a) = value;
      }
    }
  }

  return assemble_rotation_diagonal(x, dm1_active,
      {std::move(v_total_diagonal), std::move(gradient_diagonal),
       std::move(jkcaa), std::move(v_dm2)});
}

AnalyticOrbitalDiagonal build_analytic_orbital_diagonal(
    const VbScfInput& input,
    const AcceptedPointContext& accepted) {
  const auto& orbital_input = input.orbital_preparation_input;
  const auto& prepared = accepted.prepared_active_space;
  const auto& orbitals = prepared.orbital_result;
  const int n_bf = orbital_input.n_basis_functions;
  const int n_core = prepared.n_inactive_doubly_occupied_orbitals;
  const int n_active = accepted.n_active_orbitals;
  const int n_occ = n_core + n_active;
  if (n_bf < n_occ ||
      orbitals.auxiliary_orbital_matrix.rows() != n_bf ||
      orbitals.auxiliary_orbital_matrix.cols() != n_bf ||
      accepted.active_one_electron_gradient.size() !=
          static_cast<std::size_t>(n_active) * n_active) {
    throw std::invalid_argument(
        "analytic orbital diagonal requires a complete accepted AO frame");
  }
  const Eigen::Map<const Eigen::MatrixXd> overlap(
      orbital_input.ao_overlap_matrix.data(), n_bf, n_bf);
  const Eigen::Map<const Eigen::MatrixXd> active_overlap(
      orbitals.active_orbital_overlap_matrix.data(), n_active, n_active);
  Eigen::LLT<Eigen::MatrixXd> active_llt(active_overlap);
  if (active_llt.info() != Eigen::Success) {
    throw std::runtime_error("active overlap factorization failed");
  }
  const Eigen::MatrixXd active_transform = active_llt.matrixU();
  const Eigen::MatrixXd inverse_active_transform =
      active_transform.template triangularView<Eigen::Upper>().solve(
          Eigen::MatrixXd::Identity(n_active, n_active));

  Eigen::MatrixXd mo = orbitals.auxiliary_orbital_matrix;
  mo.leftCols(n_core) =
      orbitals.physical_orbital_frame.inactive_orthonormal_orbital_matrix;
  mo.middleCols(n_core, n_active) =
      orbitals.auxiliary_orbital_matrix.middleCols(n_core, n_active) *
      inverse_active_transform;
  if ((mo.transpose() * overlap * mo -
       Eigen::MatrixXd::Identity(n_bf, n_bf)).cwiseAbs().maxCoeff() > 2.0e-8) {
    throw std::runtime_error("analytic orbital frame is not AO-metric orthonormal");
  }

  const Eigen::Map<const Eigen::MatrixXd> original_dm1(
      accepted.active_one_electron_gradient.data(), n_active, n_active);
  const Eigen::MatrixXd dm1 =
      active_transform * original_dm1 * active_transform.transpose();
  const std::vector<double> packed_dm2 = pack_symmetric_dm2(
      transform_dm2(
          unpack_symmetric_active_two_rdm(
              accepted.packed_active_two_electron_gradient, n_active),
          active_transform),
      n_active);

  AnalyticDiagonalIntermediates x;
  x.n_core = n_core;
  x.n_active = n_active;
  x.h_core = mo.transpose() *
      input.ao_integral_input.ao_core_hamiltonian_matrix * mo;
  x.v_core = mo.transpose() *
      prepared.ao_effective_one_electron_result.ao_coulomb_exchange_matrix * mo;

  const auto& active_two_electron = prepared.active_space_two_electron_result;
  const Eigen::MatrixXd active_pair_transform =
      build_pair_transform(inverse_active_transform);
  ExactCtxPairMatrix active_pair_products;
  if (active_two_electron.representation ==
      ActiveSpaceTwoElectronRepresentation::ResolutionOfIdentity) {
    if (active_two_electron.ri_active_pair_factors.size() == 0 ||
        !input.ri_factorization) {
      throw std::invalid_argument(
          "RI analytic Hessian diagonal requires accepted active-pair factors");
    }
    const Eigen::MatrixXd orthogonal_active_factors =
        active_two_electron.ri_active_pair_factors * active_pair_transform;
    active_pair_products.noalias() =
        input.ri_factorization->metric_whitened_ao_pair_factors.transpose() *
        orthogonal_active_factors;
  } else if (active_two_electron.dense_ao_pair_products.size() != 0) {
    active_pair_products.noalias() =
        active_two_electron.dense_ao_pair_products * active_pair_transform;
  } else {
    // The memory-bounded exact forward path deliberately does not retain
    // G*Q_active. Rebuild only that absent accepted tensor from the AO graph.
    ExactCtxPairMatrix active_pair_map;
    build_packed_orbital_pair_map(
        mo.middleCols(n_core, n_active), &active_pair_map);
    active_pair_products = apply_ao_pair_operator(input, active_pair_map);
  }
  const int n_active_pairs = n_active * (n_active + 1) / 2;
  x.j_active.resize(static_cast<std::size_t>(n_bf) * n_active * n_active);
  // Only the diagonal MO pair (pp|ab) enters the analytic diagonal.
  // Accumulate directly: each worker owns O(n_active^2), not an AO-pair map.
#pragma omp parallel for schedule(static)
  for (int p = 0; p < n_bf; ++p) {
    Eigen::RowVectorXd contraction = Eigen::RowVectorXd::Zero(n_active_pairs);
    Eigen::Index pair = 0;
    for (int mu = 0; mu < n_bf; ++mu)
      for (int nu = 0; nu <= mu; ++nu, ++pair) {
        const double coefficient =
            mo(mu, p) * mo(nu, p) * (mu == nu ? 1.0 : 2.0);
        contraction.noalias() += coefficient * active_pair_products.row(pair);
      }
    for (int a = 0; a < n_active; ++a)
      for (int b = 0; b < n_active; ++b)
        x.j_active[slice_index(p, a, b, n_active)] =
            contraction[TwoElectronIndexer::packed_pair_index(a, b)];
  }
  ExactCtxPairMatrix active_pair_map;
  build_packed_orbital_pair_map(
      mo.middleCols(n_core, n_active), &active_pair_map);
  const Eigen::MatrixXd active_pairs =
      active_pair_map.transpose() * active_pair_products;
  x.active_eri.resize(static_cast<std::size_t>(n_active) *
                      n_active * n_active * n_active);
  for (int a = 0; a < n_active; ++a)
    for (int b = 0; b < n_active; ++b)
      for (int c = 0; c < n_active; ++c)
        for (int d = 0; d < n_active; ++d)
          x.active_eri[tensor4_index(a, b, c, d, n_active)] = active_pairs(
              TwoElectronIndexer::packed_pair_index(a, b),
              TwoElectronIndexer::packed_pair_index(c, d));
  const auto dm2 = unpack_symmetric_active_two_rdm(packed_dm2,n_active);
  const Eigen::MatrixXd density_images =
      contracted_density_images(input,mo,n_core,dm1,dm2);
  DiagonalContractions contractions;
  contractions.potential = x.v_core.diagonal()+0.5*density_images.col(0);
  contractions.gradient = diagonal_gradient(x,dm1,dm2,contractions.potential);
  contractions.one_rdm.resize(n_occ,n_active);
  contractions.two_rdm.resize(n_bf,n_active);
  for (int p = 0; p < n_bf; ++p)
    for (int a = 0; a < n_active; ++a) {
      double coulomb = 0.0;
      for (int w = 0; w < n_active; ++w)
        for (int z = 0; z < n_active; ++z) {
          const double weight = dm2[tensor4_index(w,z,a,a,n_active)]+
              2.0*(dm2[tensor4_index(a,w,z,a,n_active)]+
                   dm2[tensor4_index(w,a,z,a,n_active)]);
          coulomb += x.j_active[slice_index(p,w,z,n_active)]*weight;
        }
      contractions.two_rdm(p,a) = coulomb-density_images(p,1+n_active+a);
      if (p < n_occ) {
        double one_coulomb = 0.0;
        for (int v = 0; v < n_active; ++v)
          one_coulomb += x.j_active[slice_index(p,a,v,n_active)]*dm1(a,v);
        contractions.one_rdm(p,a) = 10.0*one_coulomb-6.0*density_images(p,1+a);
      }
    }
  x.j_pc = Eigen::MatrixXd::Zero(n_bf, n_core);
  x.k_pc = Eigen::MatrixXd::Zero(n_bf, n_core);
  ExactCtxPairMatrix core_diagonal_map(
      static_cast<Eigen::Index>(n_bf) * (n_bf + 1) / 2, n_core);
  Eigen::Index ao_pair = 0;
  for (int mu = 0; mu < n_bf; ++mu)
    for (int nu = 0; nu <= mu; ++nu, ++ao_pair)
      for (int c = 0; c < n_core; ++c)
        core_diagonal_map(ao_pair, c) =
            mo(mu, c) * mo(nu, c) * (mu == nu ? 1.0 : 2.0);
  const ExactCtxPairMatrix core_diagonal_products =
      apply_ao_pair_operator(input, core_diagonal_map);
  for (int p = 0; p < n_bf; ++p) {
    Eigen::VectorXd diagonal_pair(core_diagonal_map.rows());
    ao_pair = 0;
    for (int mu = 0; mu < n_bf; ++mu)
      for (int nu = 0; nu <= mu; ++nu, ++ao_pair)
        diagonal_pair[ao_pair] =
            mo(mu, p) * mo(nu, p) * (mu == nu ? 1.0 : 2.0);
    x.j_pc.row(p) = diagonal_pair.transpose() * core_diagonal_products;
  }

  for (int c = 0; c < n_core; ++c)
    x.k_pc.col(c) =
        2.0 * x.j_pc.col(c) - density_operator_diagonal(
            input,mo,mo.col(c),Eigen::MatrixXd::Ones(1,1));
  AnalyticOrbitalDiagonal result;
  result.mo_coefficients = std::move(mo);
  result.rotation_diagonal =
      assemble_rotation_diagonal(x, dm1, contractions);
  result.inactive_inverse_transform =
      orbitals.physical_orbital_frame
          .localized_representative_selector
          .inactive_inverse_transpose_right_transform.transpose();
  result.active_inverse_transform = inverse_active_transform;
  return result;
}

}  // namespace xmvb::vb
