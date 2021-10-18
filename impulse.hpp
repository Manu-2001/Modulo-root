#ifndef IMPULSE_HPP
#define IMPULSE_HPP

#include <cmath>
#include <type_traits>

template <typename T>
struct Impulse {
  static_assert(std::is_floating_point<T>::value,
                "Impulse: invalid template argument, floating point required");

  T x{}, y{}, z{};

  T Norm() const { return std::hypot(this->x, this->y, this->z); }

  Impulse<T>& operator+=(Impulse<T> const& other) noexcept {
    this->x += other.x;
    this->y += other.y;
    this->z += other.z;

    return *this;
  }

  Impulse<T>& operator-=(Impulse<T> const& other) noexcept {
    this->x -= other.x;
    this->y -= other.y;
    this->z -= other.z;

    return *this;
  }

  template <typename K>
  Impulse<T>& operator*=(K const& scalar) noexcept {
    static_assert(
        std::is_arithmetic<K>::value,
        "Impulse::operator*= invalid template argument, arithmetic required");

    T const Tscalar{static_cast<T>(scalar)};

    this->x *= Tscalar;
    this->y *= Tscalar;
    this->z *= Tscalar;

    return *this;
  }
};

template <class from, class to>
Impulse<to> Convert_to(Impulse<to> const& other) {
  return Impulse<to>{static_cast<to>(other.x), static_cast<to>(other.y),
                     static_cast<to>(other.z)};
}

template <typename T>
bool operator==(Impulse<T> const& lhs, Impulse<T> const& rhs) noexcept {
  return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

template <typename T>
bool operator!=(Impulse<T> const& lhs, Impulse<T> const& rhs) noexcept {
  return !(lhs == rhs);
}

template <typename T>
Impulse<T> operator+(Impulse<T> const& lhs, Impulse<T> const& rhs) noexcept {
  return Impulse<T>{lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

template <typename T>
Impulse<T> operator-(Impulse<T> const& lhs, Impulse<T> const& rhs) noexcept {
  return Impulse<T>{lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

template <typename T>
Impulse<T> operator-(Impulse<T> const& lhs) noexcept {
  return Impulse<T>{-lhs.x, -lhs.y, -lhs.z};
}

template <typename T, typename K>
Impulse<T> operator*(Impulse<T> const& lhs, K const& scalar) noexcept {
  static_assert(
      std::is_arithmetic<K>::value,
      "Impulse::operator*= invalid template argument, arithmetic required");

  T const Tsclare{static_cast<T>(scalar)};

  return Impulse<T>{lhs.x * Tsclare, lhs.y * Tsclare, lhs.z * Tsclare};
}

#endif