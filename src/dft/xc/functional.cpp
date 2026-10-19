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
  bool initialized = false;

  Component(const std::string& name, int identifier) {
    if (xc_func_init(&function, identifier, XC_UNPOLARIZED) != 0) {
      throw std::invalid_argument("Libxc could not initialize " + name);
    }
    initialized = true;
    const auto* info = xc_func_get_info(&function);
    const int flags = xc_func_info_get_flags(info);
    const int kind = xc_func_info_get_kind(info);
    if (xc_func_info_get_family(info) != XC_FAMILY_LDA ||
        kind == XC_KINETIC ||
        (flags & (XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC)) !=
            (XC_FLAGS_HAVE_EXC | XC_FLAGS_HAVE_VXC)) {
      xc_func_end(&function);
      initialized = false;
      throw std::invalid_argument(
          "XC component must be a nonkinetic LDA with EXC and VXC: " + name);
    }
  }

  ~Component() {
    if (initialized) xc_func_end(&function);
  }
  Component(const Component&) = delete;
  Component& operator=(const Component&) = delete;
};

}  // namespace

struct LdaFunctional::Impl {
  std::vector<std::unique_ptr<Component>> components;
};

LdaFunctional::LdaFunctional(const std::string& name)
    : impl_(std::make_unique<Impl>()) {
  std::size_t first = 0;
  for (;;) {
    const auto comma = name.find(',', first);
    const std::string token = trim(name.substr(first, comma - first));
    if (token.empty()) {
      throw std::invalid_argument("empty Libxc LDA component");
    }
    const int identifier = xc_functional_get_number(token.c_str());
    if (identifier <= 0) {
      throw std::invalid_argument("unknown Libxc functional: " + token);
    }
    impl_->components.push_back(std::make_unique<Component>(token, identifier));
    if (comma == std::string::npos) break;
    first = comma + 1;
  }
}

LdaFunctional::~LdaFunctional() = default;
LdaFunctional::LdaFunctional(LdaFunctional&&) noexcept = default;
LdaFunctional& LdaFunctional::operator=(LdaFunctional&&) noexcept = default;

std::string LdaFunctional::library_version() {
  int major = 0, minor = 0, micro = 0;
  xc_version(&major, &minor, &micro);
  return std::to_string(major) + "." + std::to_string(minor) + "." +
         std::to_string(micro);
}

void LdaFunctional::evaluate(const Eigen::Ref<const Eigen::VectorXd>& density,
                             Eigen::VectorXd* exc,
                             Eigen::VectorXd* vrho) const {
  if (!impl_ || !exc || !vrho || exc == vrho) {
    throw std::invalid_argument("invalid LDA functional or output vectors");
  }
  double positive_scale = 0.0;
  for (Eigen::Index i = 0; i < density.size(); ++i) {
    const double value = density[i];
    if (!std::isfinite(value)) {
      throw std::invalid_argument("LDA density must be finite");
    }
    positive_scale = std::max(positive_scale, value);
  }
  const double noise = 32.0 * std::numeric_limits<double>::epsilon() *
                       positive_scale;
  Eigen::VectorXd nonnegative(density.size());
  for (Eigen::Index i = 0; i < density.size(); ++i) {
    if (density[i] < -noise) {
      throw std::invalid_argument("LDA density is materially negative");
    }
    nonnegative[i] = std::max(0.0, density[i]);
  }

  exc->setZero(density.size());
  vrho->setZero(density.size());
  if (density.size() == 0) return;
  Eigen::VectorXd component_exc(density.size());
  Eigen::VectorXd component_vrho(density.size());
  for (const auto& component : impl_->components) {
    xc_lda_exc_vxc(&component->function,
                   static_cast<std::size_t>(density.size()),
                   nonnegative.data(), component_exc.data(),
                   component_vrho.data());
    if (!component_exc.allFinite() || !component_vrho.allFinite()) {
      throw std::runtime_error("Libxc returned nonfinite LDA results");
    }
    *exc += component_exc;
    *vrho += component_vrho;
  }
}

}  // namespace xmvb::dft
