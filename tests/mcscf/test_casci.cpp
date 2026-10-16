#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <Eigen/Eigenvalues>
#include "casscf/driver/casci.hpp"
#include "mcscf/ci/sigma.hpp"
#include "mcscf/ci/spin.hpp"
#include "core/pair_index.hpp"
#include "integrals/ao/graph.hpp"

namespace {

void require(bool condition,const char* message) {
  if (!condition) throw std::runtime_error(message);
}

void close(double a,double b,const char* message) {
  require(std::abs(a-b)<2e-9,message);
}

xmvb::casscf::CasIntegrals hubbard() {
  xmvb::casscf::CasIntegrals result;
  result.active.one_electron=Eigen::Matrix2d::Zero();
  result.active.one_electron(0,1)=result.active.one_electron(1,0)=-1;
  result.active.pair_kernel=Eigen::Matrix3d::Zero();
  result.active.pair_kernel(0,0)=result.active.pair_kernel(2,2)=4;
  result.core_energy=3.25;
  return result;
}

void check_spin() {
  for (int n=2; n<=4; ++n) {
    xmvb::mcscf::SpinSpace space(n,2,1);
    Eigen::MatrixXd x(space.alpha().size(),space.beta().size());
    for (int i=0; i<x.size(); ++i) x.data()[i]=std::sin(0.37*(i+1));
    Eigen::MatrixXd sum=Eigen::MatrixXd::Zero(x.rows(),x.cols());
    for (int spin=1; spin<=std::min(3,2*n-3); spin+=2) {
      const auto p=space.project(x,spin);
      const auto pp=space.project(p,spin);
      require((pp-p).norm()<1e-12,"spin projection is not idempotent");
      require((space.apply_squared(p)-0.25*spin*(spin+2)*p).norm()<1e-12,
          "spin projector selects the wrong eigenvalue");
      sum+=p;
    }
    require((sum-x).norm()<1e-12,"spin projectors do not resolve the identity");
    xmvb::mcscf::SpinSpace reversed(n,1,2);
    require((reversed.project(x.transpose(),1)-space.project(x,1).transpose()).norm()<1e-12,
        "spin projection depends on alpha/beta imbalance orientation");
  }
}

void check_casci() {
  const auto input=hubbard();
  const auto singlet=xmvb::casscf::solve_casci(input,2,0);
  close(singlet.energy(),3.25+2-std::sqrt(8.0),"Hubbard singlet energy is wrong");
  close(singlet.spin_squared[0],0,"singlet has spin contamination");
  const auto triplet=xmvb::casscf::solve_casci(input,2,2);
  close(triplet.energy(),3.25,"Hubbard triplet energy is wrong");
  close(triplet.spin_squared[0],2,"triplet has spin contamination");
  xmvb::casscf::CasciOptions options;
  options.n_roots=3;
  const auto roots=xmvb::casscf::solve_casci(input,2,0,options);
  close(roots.energy(0),singlet.energy(),"multiroot ground state changed");
  close(roots.energy(1),3.25+4,"second singlet root is wrong");
  close(roots.energy(2),3.25+2+std::sqrt(8.0),"third singlet root is wrong");
  options.workspace_bytes=1024;
  bool rejected=false;
  try { (void)xmvb::casscf::solve_casci(input,2,0,options); }
  catch (const std::length_error&) { rejected=true; }
  require(rejected,"CASCI ignored its caller memory budget");
  xmvb::casscf::CasIntegrals large;
  large.active.one_electron=Eigen::MatrixXd::Zero(16,16);
  large.active.pair_kernel=Eigen::MatrixXd::Zero(136,136);
  rejected=false;
  try { (void)xmvb::casscf::solve_casci(large,16,0); }
  catch (const std::length_error&) { rejected=true; }
  require(rejected,"CAS(16,16) allocated a carrier before checking its budget");
}

void check_inactive_energy() {
  xmvb::mcscf::ActiveIntegrals mo;
  mo.one_electron=Eigen::Matrix3d::Zero();
  mo.one_electron.diagonal()<<-2,-0.5,0.7;
  mo.pair_kernel=Eigen::MatrixXd::Zero(6,6);
  using P=xmvb::core::PairIndex;
  mo.pair_kernel(P::packed_pair_index(0,0),P::packed_pair_index(0,0))=0.8;
  mo.pair_kernel(P::packed_pair_index(1,1),P::packed_pair_index(0,0))=0.4;
  mo.pair_kernel(P::packed_pair_index(0,0),P::packed_pair_index(1,1))=0.4;
  mo.pair_kernel(P::packed_pair_index(1,0),P::packed_pair_index(1,0))=0.1;
  const auto eri=[&](int p,int q,int r,int s) {
    return mo.pair_kernel(P::packed_pair_index(p,q),P::packed_pair_index(r,s));
  };
  const auto cas=xmvb::casscf::make_active_integrals(mo.one_electron,eri,1,2,1.2);
  close(cas.core_energy,-2,"inactive plus nuclear energy is wrong");
  close(cas.active.one_electron(0,0),0.2,"inactive potential is missing");
  const auto doublet=xmvb::casscf::solve_casci(cas,1,1);
  close(doublet.energy(),-1.8,"inactive CASCI total energy is wrong");

  xmvb::integrals::AoIntegralInput ao;
  ao.n_basis_functions=3;
  ao.ao_core_hamiltonian_matrix=mo.one_electron;
  std::vector<int> left,right;
  std::vector<double> values;
  for (int p=0; p<6; ++p)
    for (int q=0; q<=p; ++q) {
      left.push_back(p); right.push_back(q); values.push_back(mo.pair_kernel(p,q));
    }
  ao.pair_graph=xmvb::integrals::build_ao_pair_graph(left,right,values,3);
  Eigen::MatrixXd orbitals=Eigen::Matrix3d::Identity();
  orbitals(0,0)=orbitals(1,1)=std::cos(0.37);
  orbitals(0,1)=-std::sin(0.37);
  orbitals(1,0)=std::sin(0.37);
  const Eigen::MatrixXd h=orbitals.transpose()*mo.one_electron*orbitals;
  const auto rotated=[&](int p,int q,int r,int s) {
    double value=0;
    for (int i=0; i<3; ++i)
      for (int j=0; j<3; ++j)
        for (int k=0; k<3; ++k)
          for (int l=0; l<3; ++l)
            value+=orbitals(i,p)*orbitals(j,q)*orbitals(k,r)*orbitals(l,s)*eri(i,j,k,l);
    return value;
  };
  const auto reference=xmvb::casscf::make_active_integrals(h,rotated,1,2,1.2);
  const auto shared=xmvb::casscf::make_active_integrals(ao,orbitals,1,2,1.2);
  close(shared.core_energy,reference.core_energy,"AO preparation changed inactive energy");
  require((shared.active.one_electron-reference.active.one_electron).norm()<1e-12,
      "shared AO preparation changed the active one-electron operator");
  require((shared.active.pair_kernel-reference.active.pair_kernel).norm()<1e-12,
      "shared AO preparation changed active Coulomb integrals");
  close(xmvb::casscf::solve_casci(shared,1,1).energy(),
      xmvb::casscf::solve_casci(reference,1,1).energy(),"AO-to-CASCI energy changed");
}

}  // namespace

int main(int argc,char** argv) {
  try {
    check_spin();
    check_casci();
    check_inactive_energy();
    if (argc==2) {
      // Independent molecular reference: active h, packed-pair ERIs, and
      // the inactive+nuclear constant exported by the reference program.
      std::ifstream stream(argv[1]);
      int n,electrons,spin;
      double core,reference;
      if (!(stream>>n>>electrons>>spin>>core>>reference) || n<1 || n>63)
        throw std::runtime_error("invalid molecular CASCI reference header");
      xmvb::casscf::CasIntegrals input;
      input.core_energy=core;
      input.active.one_electron.resize(n,n);
      input.active.pair_kernel.resize(n*(n+1)/2,n*(n+1)/2);
      for (int i=0; i<n; ++i)
        for (int j=0; j<n; ++j)
          if (!(stream>>input.active.one_electron(i,j)))
            throw std::runtime_error("truncated molecular one-electron integrals");
      for (int i=0; i<input.active.pair_kernel.rows(); ++i)
        for (int j=0; j<input.active.pair_kernel.cols(); ++j)
          if (!(stream>>input.active.pair_kernel(i,j)))
            throw std::runtime_error("truncated molecular pair integrals");
      const auto result=xmvb::casscf::solve_casci(input,electrons,spin);
      close(result.energy(),reference,"molecular CASCI differs from independent reference");
      std::cout<<std::setprecision(16)<<"CASCI energy="<<result.energy()
          <<" reference="<<reference<<" S2="<<result.spin_squared[0]
          <<" Davidson="<<result.active_ci.iterations<<'\n';
    }
    std::cout<<"shared direct-CI, pure-spin CASCI and inactive energy: PASS\n";
  } catch (const std::exception& error) {
    std::cerr<<error.what()<<'\n';
    return 1;
  }
}
