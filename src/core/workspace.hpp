#pragma once

#include <atomic>
#include <cstddef>
#include <memory>
#include <limits>
#include <stdexcept>
#include <utility>

namespace xmvb::workspace {

/** @brief Shared CI/RI payload allowance; not a process RSS limit.
 * Retained panels and accepted RI images are charged once for their lifetime.
 * Each sequential kernel plans temporary panels from available(). Static AO
 * integrals, immutable excitation links and eigensolver vectors are separate.
 */
inline constexpr std::size_t default_bytes = 4ULL * 1024 * 1024 * 1024;
namespace detail {
inline std::atomic<std::size_t> limit{default_bytes};
inline std::atomic<std::size_t> retained{0};
}

inline std::size_t limit() noexcept { return detail::limit.load(); }
inline std::size_t retained() noexcept { return detail::retained.load(); }
inline std::size_t available() noexcept {
  const auto cap = limit(), used = retained();
  return used < cap ? cap - used : 0;
}

/** @brief Configure before work; existing reservations remain charged. */
inline void set_limit(std::size_t bytes) {
  if (!bytes) throw std::invalid_argument("workspace must be positive");
  detail::limit.store(bytes);
}

/** @brief Move-only lifetime accounting for retained numerical payload. */
class Reservation {
 public:
  Reservation() = default;
  Reservation(const Reservation&) = delete;
  Reservation& operator=(const Reservation&) = delete;
  Reservation(Reservation&& other) noexcept : bytes_(std::exchange(other.bytes_, 0)) {}
  Reservation& operator=(Reservation&& other) noexcept {
    if (this != &other) {
      try_resize(0);
      bytes_ = std::exchange(other.bytes_, 0);
    }
    return *this;
  }
  ~Reservation() { try_resize(0); }
  std::size_t bytes() const noexcept { return bytes_; }
  bool try_resize(std::size_t bytes) noexcept {
    if (bytes <= bytes_) {
      detail::retained.fetch_sub(bytes_ - bytes);
      bytes_ = bytes;
      return true;
    }
    const auto increment = bytes - bytes_;
    auto used = retained();
    for (;;) {
      const auto cap = limit();
      if (used > cap || increment > cap - used) return false;
      if (detail::retained.compare_exchange_weak(used, used + increment)) {
        bytes_ = bytes;
        return true;
      }
    }
  }
 private:
  std::size_t bytes_ = 0;
};

/** @brief Allocate a shared retained matrix, leaving the caller's minimum
 * scratch available. Aliasing owners/copies are charged only once.
 * An empty pointer requests mathematically identical tiled regeneration.
 */
template<class Matrix>
std::shared_ptr<Matrix> make_matrix(std::size_t rows, std::size_t cols,
                                  std::size_t scratch = 0) {
  constexpr auto scalar_bytes=sizeof(typename Matrix::Scalar);
  if (rows && cols>std::numeric_limits<std::size_t>::max()/scalar_bytes/rows)
    throw std::length_error("workspace matrix size overflow");
  const auto bytes = rows * cols * scalar_bytes;
  const auto remaining=available();
  if (scratch > remaining || bytes > remaining - scratch) return {};
  struct Storage {
    Reservation reservation;
    Matrix matrix;
  };
  auto storage = std::make_shared<Storage>();
  if (!storage->reservation.try_resize(bytes)) return {};
  storage->matrix.resize(rows, cols);
  auto* matrix = &storage->matrix;
  return std::shared_ptr<Matrix>(std::move(storage), matrix);
}

} // namespace xmvb::workspace
