#include <exception>
#include <iostream>
#include <utility>

#include "cli/options.hpp"
#include "cli/run.hpp"
#ifdef XMVB_ENABLE_MPI
#include "parallel/mpi.hpp"
#endif

int main(int argc, char** argv) {
  try {
#ifdef XMVB_ENABLE_MPI
    xmvb::parallel::MpiSession session(argc, argv);
    if (session.run_workers()) return 0;
#endif
    auto options = xmvb::cli::parse_options(argc, argv);
    if (!options.has_value()) {
      return 1;
    }
    return xmvb::cli::run(std::move(*options));
  } catch (const std::exception& error) {
    std::cerr << "XMVB-CPP error: " << error.what() << '\n';
    return 2;
  }
}
