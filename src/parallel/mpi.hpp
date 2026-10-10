#pragma once

#include "vbscf/integrals/ao/contracts/service.hpp"
#include "vbscf/integrals/ao/one_electron/direct_operator.hpp"

namespace xmvb::parallel {

/** @brief One integral rank per node, with a root-only SCF driver. */
class MpiSession final : public vb::IntegralService {
public:
  MpiSession(int& argc, char**& argv);
  ~MpiSession() override;
  bool run_workers();
  vb::AoIntegralInput build(const vb::LibcintInput&) override;
  bool dispatching() const noexcept override { return dispatching_; }
  std::vector<double> h1e(const vb::AoIntegralInput&, const double*, bool) override;
  void fused(const vb::AoIntegralInput&, const double*, const double*,
             std::vector<double>*, std::vector<double>*) override;
  void pair(const vb::AoPairGraph&, int n_bf, const Eigen::Ref<const vb::PairMatrix>&,
            vb::PairMatrix*) override;
  void generated(const vb::AoIntegralInput&, const Eigen::Ref<const Eigen::MatrixXd>&,
                 const Eigen::MatrixXd*, Eigen::Index, Eigen::Index,
                 vb::PairMatrix*, vb::PairMatrix*) override;

private:
  void execute(int command, const vb::AoIntegralInput&, int rows, int columns);
  void execute_pair(const vb::AoPairGraph&, int n_bf, int rows, int columns);
  void execute_generated(const vb::AoIntegralInput&, int n_bf, int n_active);
  void send(int command, int rows, int columns);
  void receive(int command, int rows, int columns);
  vb::AoIntegralInput generate(vb::LibcintInput);
  int rank_ = 0, size_ = 1;
  bool dispatching_ = false;
  std::size_t commands_ = 0, peak_tile_bytes_ = 0, owned_integrals_ = 0;
  std::size_t generated_commands_ = 0, generated_direction_commands_ = 0;
  double action_seconds_ = 0;
  vb::AoIntegralInput worker_;
  std::vector<double> source_, adjoint_, forward_, transpose_;
  vb::PairMatrix pair_image_;
  vb::PairMatrix directional_image_;
  Eigen::MatrixXd direction_factors_;
  int generated_args_[3]{};
  vb::AoH1eFusedWorkspace h1e_workspace_;
};

}  // namespace xmvb::parallel
