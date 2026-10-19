#include "dft/xc/functional.hpp"

#include <xc.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace xmvb::dft {
namespace {

std::string trim(const std::string& text) {
  const auto first = text.find_first_not_of(" \t\n\r");
  if (first == std::string::npos) return {};
  const auto last = text.find_last_not_of(" \t\n\r");
  return text.substr(first, last - first + 1);
}

struct Component {
  xc_func_type function{};
  int family = 0;
  bool initialized = false;

  Component(const std::string& name, int identifier, int polarization) {
    if (xc_func_init(&function, identifier, polarization) != 0) {
      throw std::invalid_argument("Libxc could not initialize " + name);
    }
    initialized = true;
    const auto* info = xc_func_get_info(&function);
    const int flags = xc_func_info_get_flags(info);
    family = xc_func_info_get_family(info);
    if ((family != XC_FAMILY_LDA && family != XC_FAMILY_GGA) ||
        xc_func_info_get_kind(info) == XC_KINETIC ||
        (flags & (XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC)) !=
            (XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC)) {
      xc_func_end(&function);
      initialized = false;
      throw std::invalid_argument(
          "XC component must be pure LDA/GGA with EXC and VXC: " + name);
    }
  }

  ~Component() {
    if (initialized) xc_func_end(&function);
  }
  Component(const Component&) = delete;
  Component& operator=(const Component&) = delete;
};

}  // namespace

struct Functional::Impl {
  int n_spin = 1;
  bool gga = false;
  std::vector<std::unique_ptr<Component>> components;
};

Functional::Functional(const std::string& name, int n_spin)
    : impl_(std::make_unique<Impl>()) {
  if (n_spin != 1 && n_spin != 2) {
    throw std::invalid_argument("Libxc n_spin must be 1 or 2");
  }
  impl_->n_spin = n_spin;
  const int polarization = n_spin == 1 ? XC_UNPOLARIZED : XC_POLARIZED;
  std::size_t first = 0;
  for (;;) {
    const auto comma = name.find(',', first);
    const std::string token = trim(name.substr(first, comma - first));
    if (token.empty()) {
      throw std::invalid_argument("empty Libxc XC component");
    }
    const int identifier = xc_functional_get_number(token.c_str());
    if (identifier <= 0) {
      throw std::invalid_argument("unknown Libxc functional: " + token);
    }
    auto component = std::make_unique<Component>(token, identifier, polarization);
    impl_->gga = impl_->gga || component->family == XC_FAMILY_GGA;
    impl_->components.push_back(std::move(component));
    if (comma == std::string::npos) break;
    first = comma + 1;
  }
}

Functional::~Functional() = default;
Functional::Functional(Functional&&) noexcept = default;
Functional& Functional::operator=(Functional&&) noexcept = default;

std::string Functional::library_version() {
  int major = 0, minor = 0, micro = 0;
  xc_version(&major, &minor, &micro);
  return std::to_string(major) + "." + std::to_string(minor) + "." +
         std::to_string(micro);
}

int Functional::n_spin() const {
  if (!impl_) throw std::invalid_argument("moved-from XC functional");
  return impl_->n_spin;
}

bool Functional::is_gga() const {
  if (!impl_) throw std::invalid_argument("moved-from XC functional");
  return impl_->gga;
}

void Functional::evaluate(const Eigen::Ref<const Eigen::MatrixXd>& rho,
                          const Eigen::Ref<const Eigen::MatrixXd>& sigma,
                          XcValues* out) const {
  if (!impl_ || !out) {
    throw std::invalid_argument("invalid XC functional or output");
  }
  const int n_spin = impl_->n_spin;
  const int n_sigma = n_spin == 1 ? 1 : 3;
  const Eigen::Index n_points = rho.rows();
  if (rho.cols() != n_spin || sigma.rows() != n_points ||
      (impl_->gga ? sigma.cols() != n_sigma
                  : (sigma.cols() != 0 && sigma.cols() != n_sigma))) {
    throw std::invalid_argument("XC rho/sigma matrix shape mismatch");
  }

  double density_scale = 0.0;
  for (Eigen::Index i = 0; i < n_points; ++i) {
    for (int spin = 0; spin < n_spin; ++spin) {
      const double value = rho(i, spin);
      if (!std::isfinite(value)) {
        throw std::invalid_argument("XC density must be finite");
      }
      density_scale = std::max(density_scale, value);
    }
  }
  const double density_noise =
      32.0 * std::numeric_limits<double>::epsilon() * density_scale;
  std::vector<double> packed_rho(static_cast<std::size_t>(n_points) * n_spin);
  for (Eigen::Index i = 0; i < n_points; ++i) {
    for (int spin = 0; spin < n_spin; ++spin) {
      const double value = rho(i, spin);
      if (value < -density_noise) {
        throw std::invalid_argument("XC density is materially negative");
      }
      packed_rho[static_cast<std::size_t>(i) * n_spin + spin] =
          std::max(0.0, value);
    }
  }

  std::vector<double> packed_sigma;
  if (impl_->gga) {
    double sigma_scale = 0.0;
    for (Eigen::Index i = 0; i < n_points; ++i) {
      for (int channel = 0; channel < n_sigma; ++channel) {
        const double value = sigma(i, channel);
        if (!std::isfinite(value)) {
          throw std::invalid_argument("XC sigma must be finite");
        }
        if (n_spin == 1 || channel != 1) {
          sigma_scale = std::max(sigma_scale, value);
        }
      }
    }
    const double sigma_noise =
        32.0 * std::numeric_limits<double>::epsilon() * sigma_scale;
    packed_sigma.resize(static_cast<std::size_t>(n_points) * n_sigma);
    for (Eigen::Index i = 0; i < n_points; ++i) {
      for (int channel = 0; channel < n_sigma; ++channel) {
        const double value = sigma(i, channel);
        if ((n_spin == 1 || channel != 1) && value < -sigma_noise) {
          throw std::invalid_argument("XC diagonal sigma is materially negative");
        }
        packed_sigma[static_cast<std::size_t>(i) * n_sigma + channel] =
            (n_spin == 2 && channel == 1) ? value : std::max(0.0, value);
      }
    }
  }

  out->exc.setZero(n_points);
  out->vrho.setZero(n_points, n_spin);
  out->vsigma.setZero(n_points, impl_->gga ? n_sigma : 0);
  if (n_points == 0) return;

  std::vector<double> component_exc(static_cast<std::size_t>(n_points));
  std::vector<double> component_vrho(static_cast<std::size_t>(n_points) * n_spin);
  std::vector<double> component_vsigma;
  if (impl_->gga) {
    component_vsigma.resize(static_cast<std::size_t>(n_points) * n_sigma);
  }
  for (const auto& component : impl_->components) {
    if (component->family == XC_FAMILY_GGA) {
      xc_gga_exc_vxc(&component->function, static_cast<std::size_t>(n_points),
                     packed_rho.data(), packed_sigma.data(),
                     component_exc.data(), component_vrho.data(),
                     component_vsigma.data());
    } else {
      xc_lda_exc_vxc(&component->function, static_cast<std::size_t>(n_points),
                     packed_rho.data(), component_exc.data(),
                     component_vrho.data());
    }
    for (Eigen::Index i = 0; i < n_points; ++i) {
      const std::size_t point = static_cast<std::size_t>(i);
      if (!std::isfinite(component_exc[point])) {
        throw std::runtime_error("Libxc returned nonfinite XC energy");
      }
      out->exc[i] += component_exc[point];
      for (int spin = 0; spin < n_spin; ++spin) {
        const double value = component_vrho[point * n_spin + spin];
        if (!std::isfinite(value)) {
          throw std::runtime_error("Libxc returned nonfinite XC vrho");
        }
        out->vrho(i, spin) += value;
      }
      if (component->family == XC_FAMILY_GGA) {
        for (int channel = 0; channel < n_sigma; ++channel) {
          const double value = component_vsigma[point * n_sigma + channel];
          if (!std::isfinite(value)) {
            throw std::runtime_error("Libxc returned nonfinite XC vsigma");
          }
          out->vsigma(i, channel) += value;
        }
      }
    }
  }
}

}  // namespace xmvb::dft
