#include "vbscf/structures/assembly/action.hpp"

#include <algorithm>
#include <stdexcept>

#include "core/openmp.hpp"
#include "core/timing.hpp"
#include "vbscf/structures/assembly/sources.hpp"
#include "vbscf/structures/orthogonal_ci/exterior_transform.hpp"
#include "vbscf/structures/overlap/representation.hpp"

namespace xmvb::vb {
namespace {

/** @brief A bounded single-spin factor and, when requested, its tangent. */
struct Factor {
  Eigen::MatrixXd value, tangent;
  Factor(int rows, int columns, bool directional)
      : value(Eigen::MatrixXd::Zero(rows, columns)) {
    if (directional) tangent = value;
  }
};

void transform(const ExteriorOrbitalTransform& map,
    const ExteriorTransformDirection* direction, Factor& factor, bool adjoint = false) {
  if (direction) {
    if (adjoint) map.apply_directional_adjoint_left(*direction, &factor.value, &factor.tangent);
    else map.apply_directional_left(*direction, &factor.value, &factor.tangent);
  } else {
    if (adjoint) map.apply_adjoint_left(&factor.value);
    else map.apply_left(&factor.value);
  }
}

void transform_right(const ExteriorOrbitalTransform& map,
    const ExteriorTransformDirection* direction,Factor& factor,bool adjoint=false) {
  if (direction) {
    if (adjoint) map.apply_directional_adjoint_right(*direction,&factor.value,&factor.tangent);
    else map.apply_directional_right(*direction,&factor.value,&factor.tangent);
  } else {
    if (adjoint) map.apply_adjoint_right(&factor.value);
    else map.apply_right(&factor.value);
  }
}

Factor physical_product(const StructureAction& action,const Eigen::Ref<const Eigen::VectorXd>& coefficients,
    const ExteriorOrbitalTransform& alpha,const ExteriorOrbitalTransform& beta,
    const ExteriorTransformDirection* da=nullptr,const ExteriorTransformDirection* db=nullptr) {
  Factor result(0,0,false);
  result.value=action.expand_structure_block(coefficients);
  if (da) result.tangent=Eigen::MatrixXd::Zero(result.value.rows(),result.value.cols());
  transform(alpha,da,result);
  transform_right(beta,db,result);
  return result;
}

bool same_chart(const ExteriorOrbitalTransform& a,const ExteriorOrbitalTransform& b) {
  return a.n_orbitals()==b.n_orbitals() &&
      (a.orbital_transform().array()==b.orbital_transform().array()).all();
}

int panel_width(std::size_t bytes, int na, int nb, int columns, int copies) {
  const std::size_t row_bytes = copies * sizeof(double) *
      (static_cast<std::size_t>(na) + nb);
  if (bytes < row_bytes) throw std::length_error("cross overlap budget cannot hold one factor column");
  return static_cast<int>(std::max<std::size_t>(1, std::min<std::size_t>(columns,
      std::min(bytes / row_bytes, static_cast<std::size_t>(effective_openmp_thread_count()) * 8))));
}

} // namespace

Eigen::MatrixXd StructureAction::cross_overlap(
    const Eigen::Ref<const Eigen::MatrixXd>& vectors,
    const ExteriorOrbitalTransform& la, const ExteriorOrbitalTransform& lb,
    const ExteriorOrbitalTransform& ra, const ExteriorOrbitalTransform& rb,
    const ExteriorTransformDirection* dla, const ExteriorTransformDirection* dlb,
    const ExteriorTransformDirection* dra, const ExteriorTransformDirection* drb,
    std::size_t workspace_bytes) const {
  const bool directional = dla || dlb || dra || drb;
  if (vectors.rows() != n_structures_ || vectors.cols() <= 0 ||
      la.dimension() != n_unique_alpha_ || ra.dimension() != n_unique_alpha_ ||
      lb.dimension() != n_unique_beta_ || rb.dimension() != n_unique_beta_ ||
      (directional && !(dla && dlb && dra && drb)))
    throw std::invalid_argument("invalid selected cross overlap");
  const auto& sources = *sources_;
  const std::size_t products=static_cast<std::size_t>(n_unique_alpha_)*n_unique_beta_;
  if (spin_products_.size()==products &&
      products<=workspace_bytes/((directional ? 2:1)*sizeof(double))) {
    Eigen::MatrixXd result(n_structures_,vectors.cols());
    for (Eigen::Index column=0; column<vectors.cols(); ++column) {
      auto image=physical_product(*this,vectors.col(column),ra,rb,dra,drb);
      transform_right(lb,dlb,image,true);
      transform(la,dla,image,true);
      result.col(column)=contract_spin_product_block(directional ? image.tangent:image.value);
    }
    return result;
  }
  // Column-singleton support (including covalent structures) has one alpha string
  // per touched beta. Its scalar coefficient commutes with both exterior
  // maps. Generate each common basis image once, not once per CI RHS.
  const bool singleton=std::all_of(sources.columns.begin(),sources.columns.end(),
      [&](int b) { return sources.beta_groups[b].size()==1; });
  if (singleton && vectors.cols()>1) {
    const std::size_t groups=spin_products_.size();
    const std::size_t factor_rows=(directional ? 4:2)*
        (static_cast<std::size_t>(n_unique_alpha_)+n_unique_beta_);
    // One factor panel, one sampled Gram panel and two RHS product blocks.
    // The same budget bounds all of them; there is no groups-squared cache.
    const std::size_t capacity=workspace_bytes/sizeof(double);
    if (capacity>=factor_rows+3*groups) {
      const int width=std::min<std::size_t>(vectors.cols(),
          (capacity-factor_rows-groups)/(2*groups));
      const int extent=std::min<std::size_t>(sources.columns.size(),
          std::min<std::size_t>(8*effective_openmp_thread_count(),
              (capacity-2*groups*width)/(factor_rows+groups)));
      Eigen::MatrixXd result=Eigen::MatrixXd::Zero(n_structures_,vectors.cols());
      for (int first=0; first<vectors.cols(); first+=width) {
        const int count_rhs=std::min<int>(width,vectors.cols()-first);
        Eigen::MatrixXd coefficients=Eigen::MatrixXd::Zero(groups,count_rhs);
#pragma omp parallel for schedule(static) if(groups>32)
        for (int g=0;g<static_cast<int>(groups);++g)
          for (std::size_t t=spin_term_offsets_[g];t<spin_term_offsets_[g+1];++t) {
            const auto& term=spin_terms_[t];
            coefficients.row(g).noalias()+=term.coefficient*vectors.row(term.structure).segment(first,count_rhs);
          }
        Eigen::MatrixXd image=Eigen::MatrixXd::Zero(groups,count_rhs);
        for (int begin=0;begin<static_cast<int>(groups);begin+=extent) {
          const int count=std::min<int>(extent,groups-begin);
          Factor a(n_unique_alpha_,count,directional),b(n_unique_beta_,count,directional);
          for (int j=0;j<count;++j) {
            a.value(spin_products_[begin+j]/n_unique_beta_,j)=1;
            b.value(spin_products_[begin+j]%n_unique_beta_,j)=1;
          }
          transform(ra,dra,a); transform(la,dla,a,true);
          transform(rb,drb,b); transform(lb,dlb,b,true);
          Eigen::MatrixXd gram(groups,count);
#pragma omp parallel for schedule(static) if(groups>32)
          for (int g=0;g<static_cast<int>(groups);++g) {
            const int ai=spin_products_[g]/n_unique_beta_,bi=spin_products_[g]%n_unique_beta_;
            if (directional)
              gram.row(g)=(a.tangent.row(ai).array()*b.value.row(bi).array()+
                  a.value.row(ai).array()*b.tangent.row(bi).array()).matrix();
            else gram.row(g)=(a.value.row(ai).array()*b.value.row(bi).array()).matrix();
          }
          image.noalias()+=gram*coefficients.middleRows(begin,count);
        }
#pragma omp parallel for schedule(static) if(n_structures_>32)
        for (int s=0;s<n_structures_;++s)
          for (auto t=structure_term_offsets_[s];t<structure_term_offsets_[s+1];++t) {
            const auto& term=spin_terms_[structure_term_indices_[t]];
            result.row(s).segment(first,count_rhs).noalias()+=term.coefficient*image.row(term.group);
          }
      }
      return result;
    }
  }
  const int extent = panel_width(workspace_bytes, n_unique_alpha_, n_unique_beta_,
      sources.columns.size(), directional ? 8 : 4);
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(n_structures_, vectors.cols());
  for (int column = 0; column < vectors.cols(); ++column) {
    Eigen::VectorXd coefficients = Eigen::VectorXd::Zero(spin_products_.size());
    Eigen::VectorXd products = Eigen::VectorXd::Zero(spin_products_.size());
#pragma omp parallel for schedule(static) if(spin_products_.size() > 32)
    for (int group = 0; group < static_cast<int>(spin_products_.size()); ++group)
      for (std::size_t term = spin_term_offsets_[group]; term < spin_term_offsets_[group + 1]; ++term)
        coefficients[group] += spin_terms_[term].coefficient * vectors(spin_terms_[term].structure, column);

    for (int begin = 0; begin < static_cast<int>(sources.columns.size()); begin += extent) {
      const int count = std::min<int>(extent, sources.columns.size() - begin);
      Factor a(n_unique_alpha_, count, directional), b(n_unique_beta_, count, directional);
      for (int j = 0; j < count; ++j) {
        const int beta = sources.columns[begin + j];
        b.value(beta, j) = 1.0;
        for (const int group : sources.beta_groups[beta])
          a.value(spin_products_[group] / n_unique_beta_, j) = coefficients[group];
      }
      // K_alpha X_J and K_beta E_J, with K_sigma = T_left^T T_right.
      // Associativity removes the transformed alpha-by-beta product entirely.
      transform(ra, dra, a); transform(la, dla, a, true);
      transform(rb, drb, b); transform(lb, dlb, b, true);
#pragma omp parallel for schedule(static) if(spin_products_.size() > 32)
      for (int group = 0; group < static_cast<int>(spin_products_.size()); ++group) {
        const int alpha = spin_products_[group] / n_unique_beta_;
        const int beta = spin_products_[group] % n_unique_beta_;
        products[group] += directional
            ? a.tangent.row(alpha).dot(b.value.row(beta)) +
                  a.value.row(alpha).dot(b.tangent.row(beta))
            : a.value.row(alpha).dot(b.value.row(beta));
      }
    }
#pragma omp parallel for schedule(static) if(n_structures_ > 32)
    for (int structure = 0; structure < n_structures_; ++structure)
      for (std::size_t offset = structure_term_offsets_[structure];
           offset < structure_term_offsets_[structure + 1]; ++offset) {
        const auto& term = spin_terms_[structure_term_indices_[offset]];
        result(structure, column) += term.coefficient * products[term.group];
      }
  }
  return result;
}

UpperPair StructureAction::cross_overlap_adjoint(
    const Eigen::VectorXd& left, const Eigen::VectorXd& right,
    const ExteriorOrbitalTransform& la, const ExteriorOrbitalTransform& lb,
    const ExteriorOrbitalTransform& ra, const ExteriorOrbitalTransform& rb,
    const Eigen::MatrixXd* delta_left, const Eigen::MatrixXd* delta_right,
    std::size_t workspace_bytes, OverlapSide side) const {
  const bool directional = delta_left != nullptr;
  if (left.size() != n_structures_ || right.size() != n_structures_ ||
      directional != (delta_right != nullptr))
    throw std::invalid_argument("invalid cross overlap adjoint");
  // Scalar cross-Gram symmetry turns a left-only pullback into the same
  // right-only contraction, including both chart directions when present.
  if (side == OverlapSide::Left) {
    auto swapped = cross_overlap_adjoint(right,left,ra,rb,la,lb,
        delta_right,delta_left,workspace_bytes,OverlapSide::Right);
    return {std::move(swapped.right),{}};
  }
  const bool left_requested = side == OverlapSide::Both;
  const std::size_t products=static_cast<std::size_t>(n_unique_alpha_)*n_unique_beta_;
  const bool dl=directional && !delta_left->isZero(0.0);
  const bool dr=directional && !delta_right->isZero(0.0);
  if (spin_products_.size()==products && same_chart(la,lb) && same_chart(ra,rb) &&
      products<=workspace_bytes/((2+static_cast<int>(dl)+static_cast<int>(dr))*sizeof(double))) {
    if (directional && !dl && !dr)
      return {left_requested ? Eigen::MatrixXd::Zero(la.n_orbitals(),la.n_orbitals()).eval() : Eigen::MatrixXd{},
              Eigen::MatrixXd::Zero(ra.n_orbitals(),ra.n_orbitals())};
    ExteriorTransformDirection dla,dlb,dra,drb;
    if (dl) { dla=la.direction(*delta_left); dlb=lb.direction(*delta_left); }
    if (dr) { dra=ra.direction(*delta_right); drb=rb.direction(*delta_right); }
    const auto l=physical_product(*this,left,la,lb,dl ? &dla:nullptr,dl ? &dlb:nullptr);
    const auto r=physical_product(*this,right,ra,rb,dr ? &dra:nullptr,dr ? &drb:nullptr);
    // D(R,L)_{pq}=<R|E_pq|L>; the opposite orientation is D(R,L)^T.
    const Eigen::MatrixXd moment=carrier_generator_adjoint(r.value,l.value);
    const Eigen::MatrixXd transpose=moment.transpose();
    if (!directional)
      return {left_requested ? la.moment_pullback(moment) : Eigen::MatrixXd{},
              ra.moment_pullback(transpose)};
    Eigen::MatrixXd delta=Eigen::MatrixXd::Zero(moment.rows(),moment.cols());
    if (dl) delta+=carrier_generator_adjoint(r.value,l.tangent);
    if (dr) delta+=carrier_generator_adjoint(r.tangent,l.value);
    const Eigen::MatrixXd delta_transpose=delta.transpose();
    return {left_requested ? la.moment_pullback(moment,&delta,delta_left) : Eigen::MatrixXd{},
            ra.moment_pullback(transpose,&delta_transpose,delta_right)};
  }
  const auto& sources = *sources_;
  const int extent = panel_width(workspace_bytes, n_unique_alpha_, n_unique_beta_,
      sources.columns.size(), directional ? 16 : 8);
  ExteriorTransformDirection dla, dlb, dra, drb;
  if (directional) {
    dla = la.direction(*delta_left); dlb = lb.direction(*delta_left);
    dra = ra.direction(*delta_right); drb = rb.direction(*delta_right);
  }
  Eigen::MatrixXd coefficients = Eigen::MatrixXd::Zero(spin_products_.size(), 2);
#pragma omp parallel for schedule(static) if(spin_products_.size() > 32)
  for (int group = 0; group < static_cast<int>(spin_products_.size()); ++group)
    for (std::size_t t = spin_term_offsets_[group]; t < spin_term_offsets_[group + 1]; ++t) {
      const auto& term = spin_terms_[t];
      coefficients(group, 0) += term.coefficient * left[term.structure];
      coefficients(group, 1) += term.coefficient * right[term.structure];
    }
  UpperPair result{left_requested ? Eigen::MatrixXd::Zero(la.n_orbitals(), la.n_orbitals()).eval() : Eigen::MatrixXd{},
                   Eigen::MatrixXd::Zero(ra.n_orbitals(), ra.n_orbitals())};
  for (int begin = 0; begin < static_cast<int>(sources.columns.size()); begin += extent) {
    const int count = std::min<int>(extent, sources.columns.size() - begin);
    Factor u(n_unique_alpha_, count, directional), v(n_unique_beta_, count, directional);
    for (int j = 0; j < count; ++j) {
      const int beta = sources.columns[begin + j];
      v.value(beta, j) = 1.0;
      for (const int group : sources.beta_groups[beta])
        u.value(spin_products_[group] / n_unique_beta_, j) = coefficients(group, 1);
    }
    transform(ra, directional ? &dra : nullptr, u);
    transform(rb, directional ? &drb : nullptr, v);
    Factor a = u, b = v;
    transform(la, directional ? &dla : nullptr, a, true);
    transform(lb, directional ? &dlb : nullptr, b, true);
    Factor za(n_unique_alpha_, count, directional), zb(n_unique_beta_, count, directional);
    // Seeds X_left B and X_left^T A use only the original touched products.
#pragma omp parallel for schedule(static) if(n_unique_alpha_ > 32)
    for (int alpha = 0; alpha < n_unique_alpha_; ++alpha)
      for (int group = sources.alpha_offsets[alpha]; group < sources.alpha_offsets[alpha + 1]; ++group) {
        const int beta = spin_products_[group] % n_unique_beta_;
        za.value.row(alpha).noalias() += coefficients(group, 0) * b.value.row(beta);
        if (directional) za.tangent.row(alpha).noalias() += coefficients(group, 0) * b.tangent.row(beta);
      }
#pragma omp parallel for schedule(static) if(n_unique_beta_ > 32)
    for (int beta = 0; beta < n_unique_beta_; ++beta)
      for (const int group : sources.beta_groups[beta]) {
        const int alpha = spin_products_[group] / n_unique_beta_;
        zb.value.row(beta).noalias() += coefficients(group, 0) * a.value.row(alpha);
        if (directional) zb.tangent.row(beta).noalias() += coefficients(group, 0) * a.tangent.row(alpha);
      }
    transform(la, directional ? &dla : nullptr, za);
    transform(lb, directional ? &dlb : nullptr, zb);
    const auto pullback = [&](const ExteriorOrbitalTransform& map,
        const Factor& seed, const Factor& image, const Eigen::MatrixXd* delta) {
      return map.generator_pullback(seed.value, image.value,
          directional ? &seed.tangent : nullptr, directional ? &image.tangent : nullptr, delta);
    };
    if (left_requested)
      result.left += pullback(la, u, za, delta_left) + pullback(lb, v, zb, delta_left);
    result.right += pullback(ra, za, u, delta_right) + pullback(rb, zb, v, delta_right);
  }
  return result;
}

OverlapTransportCurvature StructureAction::cross_overlap_transport(
    const Eigen::VectorXd& lambda, const Eigen::VectorXd& coefficients,
    const Eigen::VectorXd& coefficient_direction,
    const ExteriorOrbitalTransform& alpha, const ExteriorOrbitalTransform& beta,
    const Eigen::MatrixXd& direction, std::size_t workspace_bytes) const {
  profiling::ScopedPhase phase(profiling::Phase::CiOverlapTransport);
  const int n=alpha.n_orbitals();
  if (lambda.size()!=n_structures_ || coefficients.size()!=n_structures_ ||
      coefficient_direction.size()!=n_structures_ || alpha.dimension()!=n_unique_alpha_ ||
      beta.dimension()!=n_unique_beta_ || beta.n_orbitals()!=n ||
      direction.rows()!=n || direction.cols()!=n)
    throw std::invalid_argument("invalid overlap transport tile inputs");
  const auto da=alpha.direction(direction), db=beta.direction(direction);
  OverlapTransportCurvature result{Eigen::MatrixXd::Zero(n,n),
      Eigen::MatrixXd::Zero(n,n),Eigen::VectorXd::Zero(n_structures_)};
  const std::size_t products=static_cast<std::size_t>(n_unique_alpha_)*n_unique_beta_;
  // Broad support already has a cheaper one-body-moment contraction order.
  // Keep it within the existing 64 MiB overlap allowance; narrower selections
  // never expand here. Three shared images suffice for all consumers.
  if (spin_products_.size()==products && same_chart(alpha,beta) &&
      products<=workspace_bytes/(3*sizeof(double))) {
    auto accepted=physical_product(*this,lambda,alpha,beta);
    auto right=physical_product(*this,coefficients,alpha,beta,&da,&db);
    const Eigen::MatrixXd moment=carrier_generator_adjoint(accepted.value,right.value);
    const Eigen::MatrixXd delta=carrier_generator_adjoint(accepted.value,right.tangent);
    result.directional=alpha.moment_pullback(moment,&delta,&direction);
    right.value.resize(0,0);
    right.value=expand_structure_block(coefficient_direction);
    alpha.apply_left(&right.value); beta.apply_right(&right.value);
    result.mixed=alpha.moment_pullback(carrier_generator_adjoint(accepted.value,right.value));
    right.value.resize(0,0);
    right.tangent.setZero();
    // Reuse the dead directional image as the transpose-action workspace.
    beta.apply_directional_adjoint_right(db,&accepted.value,&right.tangent);
    alpha.apply_directional_adjoint_left(da,&accepted.value,&right.tangent);
    result.image=contract_spin_product_block(right.tangent).col(0);
    return result;
  }
  const auto& sources=*sources_;
  // Eight alpha and seven beta row panels coexist at most. One allowance
  // bounds all three consumers together, not three independent workspaces.
  const std::size_t minimum_bytes=8*sizeof(double)*
      (static_cast<std::size_t>(n_unique_alpha_)+n_unique_beta_);
  const auto panel_bytes=std::max(std::min<std::size_t>(workspace_bytes,8ULL*1024*1024),minimum_bytes);
  const int extent=panel_width(panel_bytes,n_unique_alpha_,n_unique_beta_,
      sources.columns.size(),8);
  const Eigen::VectorXd weights=expand_support(lambda);
  const auto fixed_adjoint=[](const ExteriorOrbitalTransform& map,Factor& factor) {
    map.apply_adjoint_left(&factor.value);
    if (factor.tangent.size()) map.apply_adjoint_left(&factor.tangent);
  };
  const auto fixed_forward=[](const ExteriorOrbitalTransform& map,Factor& factor) {
    map.apply_left(&factor.value);
    if (factor.tangent.size()) map.apply_left(&factor.tangent);
  };
  for (int begin=0; begin<static_cast<int>(sources.columns.size()); begin+=extent) {
    const int count=std::min<int>(extent,sources.columns.size()-begin);
    Factor u(n_unique_alpha_,count,true), v(n_unique_beta_,count,true);
    Factor mixed_u(n_unique_alpha_,count,false);
    for (int j=0; j<count; ++j) {
      const int b=sources.columns[begin+j];
      v.value(b,j)=1;
      for (const int g:sources.beta_groups[b]) {
        const int a=spin_products_[g]/n_unique_beta_;
        for (auto t=spin_term_offsets_[g]; t<spin_term_offsets_[g+1]; ++t) {
          const auto& term=spin_terms_[t];
          u.value(a,j)+=term.coefficient*coefficients[term.structure];
          mixed_u.value(a,j)+=term.coefficient*coefficient_direction[term.structure];
        }
      }
    }
    transform(alpha,&da,u); transform(beta,&db,v);
    transform(alpha,nullptr,mixed_u);
    Factor a=u,b=v, mixed_a=mixed_u;
    fixed_adjoint(alpha,a); fixed_adjoint(beta,b); fixed_adjoint(alpha,mixed_a);
    Factor za(n_unique_alpha_,count,true), zb(n_unique_beta_,count,true);
    Factor mixed_zb(n_unique_beta_,count,false);
#pragma omp parallel for schedule(static) if(n_unique_alpha_>32)
    for (int ai=0; ai<n_unique_alpha_; ++ai)
      for (int g=sources.alpha_offsets[ai]; g<sources.alpha_offsets[ai+1]; ++g) {
        const int bi=spin_products_[g]%n_unique_beta_;
        za.value.row(ai).noalias()+=weights[g]*b.value.row(bi);
        za.tangent.row(ai).noalias()+=weights[g]*b.tangent.row(bi);
      }
#pragma omp parallel for schedule(static) if(n_unique_beta_>32)
    for (int bi=0; bi<n_unique_beta_; ++bi)
      for (const int g:sources.beta_groups[bi]) {
        const int ai=spin_products_[g]/n_unique_beta_;
        zb.value.row(bi).noalias()+=weights[g]*a.value.row(ai);
        zb.tangent.row(bi).noalias()+=weights[g]*a.tangent.row(ai);
        mixed_zb.value.row(bi).noalias()+=weights[g]*mixed_a.value.row(ai);
      }
    fixed_forward(alpha,za); fixed_forward(beta,zb); fixed_forward(beta,mixed_zb);
    result.directional+=alpha.generator_pullback(za.value,u.value,
        &za.tangent,&u.tangent,&direction)+beta.generator_pullback(zb.value,v.value,
        &zb.tangent,&v.tangent,&direction);
    result.mixed+=alpha.generator_pullback(za.value,mixed_u.value)+
        beta.generator_pullback(mixed_zb.value,v.value);
    // The transpose directional image reuses za=T_a X_lambda T_b^T v:
    // d(T_a^T za)=dT_a^T za + T_a^T T_a X_lambda T_b^T dv.
    // Sampling only touched original products is exactly W^T, not closure.
    alpha.apply_directional_adjoint_left(da,&za.value,&za.tangent);
    for (int j=0; j<count; ++j)
      for (const int g:sources.beta_groups[sources.columns[begin+j]]) {
        const double value=za.tangent(spin_products_[g]/n_unique_beta_,j);
        for (auto t=spin_term_offsets_[g]; t<spin_term_offsets_[g+1]; ++t) {
          const auto& term=spin_terms_[t];
          result.image[term.structure]+=term.coefficient*value;
        }
      }
  }
  return result;
}

} // namespace xmvb::vb
