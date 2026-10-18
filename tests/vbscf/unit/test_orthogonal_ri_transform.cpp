#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include <Eigen/Cholesky>
#include <Eigen/Core>

#include "integrals/active/ri_transform.hpp"
#include "vbscf/structures/orthogonal_ci/integrals.hpp"

namespace {

using Matrix = Eigen::MatrixXd;

Matrix fixture(int rows, int cols, double phase) {
  Matrix result(rows, cols);
  for (int j=0; j<cols; ++j)
    for (int i=0; i<rows; ++i)
      result(i,j)=std::sin(phase+0.37*(i+1)*(j+2));
  return result;
}

void close(const Matrix& actual, const Matrix& expected, double tolerance=2e-10) {
  if (actual.rows()!=expected.rows() || actual.cols()!=expected.cols() ||
      !actual.allFinite() ||
      (actual-expected).norm()>tolerance*std::max(1.0,expected.norm()))
    throw std::runtime_error("fused orthogonal RI transform differs from independent reference");
}

std::vector<double> values(const Matrix& matrix) {
  return {matrix.data(),matrix.data()+matrix.size()};
}

std::vector<double> pack(const Matrix& matrix) {
  std::vector<double> result;
  for (int j=0; j<matrix.cols(); ++j)
    for (int i=0; i<=j; ++i) result.push_back(matrix(i,j));
  return result;
}

void check(int n) {
  const int n_bf=n+3, n_aux=7, pairs=n*(n+1)/2;
  const Matrix c=fixture(n_bf,n,0.1)+2*Matrix::Identity(n_bf,n);
  const Matrix dc=0.13*fixture(n_bf,n,0.7);
  const Matrix ao_h=fixture(n_bf,n_bf,0.2);
  const Matrix h=0.5*(ao_h+ao_h.transpose());
  std::vector<Matrix> l;
  for (int p=0;p<n_aux;++p) {
    const Matrix raw=fixture(n_bf,n_bf,0.4+p);
    l.push_back(0.5*(raw+raw.transpose()));
  }
  const auto factors=[&](const Matrix& coefficients) {
    Matrix result(n_aux,pairs);
    for (int p=0;p<n_aux;++p) {
      const Matrix transformed=coefficients.transpose()*l[p]*coefficients;
      for (int a=0;a<n;++a)
        for (int b=0;b<=a;++b) result(p,a*(a+1)/2+b)=transformed(a,b);
    }
    return result;
  };
  const auto prepare=[&](const Matrix& coefficients, bool ri) {
    xmvb::integrals::ActiveSpaceTwoElectronResult two;
    two.ri_active_pair_factors=factors(coefficients);
    two.n_auxiliary_functions=n_aux;
    two.representation=ri
        ? xmvb::integrals::ActiveSpaceTwoElectronRepresentation::ResolutionOfIdentity
        : xmvb::integrals::ActiveSpaceTwoElectronRepresentation::PackedExact;
    if (!ri) two.packed_active_two_electron_integrals=pack(
        two.ri_active_pair_factors.transpose()*two.ri_active_pair_factors);
    return xmvb::vb::orthogonalize_active_integrals(
        values(coefficients.transpose()*coefficients),
        coefficients.transpose()*h*coefficients,two,n);
  };
  const auto accepted=prepare(c,true), exact=prepare(c,false);
  close(accepted.pair_kernel,exact.pair_kernel);
  const Matrix overlap=c.transpose()*c;
  Eigen::LLT<Matrix> chol(overlap);
  const Matrix q=chol.matrixU().solve(Matrix::Identity(n,n));
  const Matrix independent=factors(c*q);
  close(accepted.ri_pair_factors,independent);
  close(accepted.pair_kernel,independent.transpose()*independent);

  const Matrix b=factors(c);
  Matrix db(n_aux,pairs);
  for (int p=0;p<n_aux;++p) {
    const Matrix direction=dc.transpose()*l[p]*c+c.transpose()*l[p]*dc;
    for (int a=0;a<n;++a)
      for (int k=0;k<=a;++k) db(p,a*(a+1)/2+k)=direction(a,k);
  }
  const Matrix ds=dc.transpose()*c+c.transpose()*dc;
  const Matrix dh=dc.transpose()*h*c+c.transpose()*h*dc;
  const Matrix dv=db.transpose()*b+b.transpose()*db;
  const auto direction=xmvb::vb::orthogonalize_active_integral_direction(
      accepted,values(ds),dh,{},n,&db);
  const auto reference=xmvb::vb::orthogonalize_active_integral_direction(
      exact,values(ds),dh,pack(dv),n);
  close(direction.pair_kernel,reference.pair_kernel);
  close(direction.one_electron,reference.one_electron);
  constexpr double step=1e-5;
  const auto plus=prepare(c+step*dc,true), minus=prepare(c-step*dc,true);
  close(direction.pair_kernel,(plus.pair_kernel-minus.pair_kernel)/(2*step),2e-8);

  const Matrix raw_g=fixture(pairs,pairs,0.9), raw_dg=fixture(pairs,pairs,1.3);
  const Matrix g=0.5*(raw_g+raw_g.transpose());
  const Matrix dg=0.5*(raw_dg+raw_dg.transpose());
  const Matrix raw_hg=fixture(n,n,0.9), raw_dhg=fixture(n,n,1.5);
  const Matrix hg=0.5*(raw_hg+raw_hg.transpose());
  const Matrix dhg=0.5*(raw_dhg+raw_dhg.transpose());
  const Matrix generator=fixture(n,n,0.3), dgenerator=fixture(n,n,0.6);
  const auto derivative=xmvb::vb::backpropagate_orthogonal_active_integral_adjoint_direction(
      accepted,direction,hg,g,generator,dhg,dg,dgenerator);
  const auto gp=xmvb::vb::backpropagate_orthogonal_active_integral_adjoint(
      plus,hg+step*dhg,g+step*dg,generator+step*dgenerator);
  const auto gm=xmvb::vb::backpropagate_orthogonal_active_integral_adjoint(
      minus,hg-step*dhg,g-step*dg,generator-step*dgenerator);
  close(derivative.overlap,(gp.overlap-gm.overlap)/(2*step),2e-8);
  close(derivative.one_electron,(gp.one_electron-gm.one_electron)/(2*step),2e-8);
  close(derivative.pair_kernel,(gp.pair_kernel-gm.pair_kernel)/(2*step),2e-8);
}

}  // namespace

int main() {
  for (int n : {1,3,5}) check(n);
}
