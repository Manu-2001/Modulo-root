#ifndef POINT_HPP
#define POINT_HPP

#include <cmath>
#include <type_traits>

template <typename T>
struct Point {
  static_assert(std::is_floating_point<T>::value,
                "Point: invalid template argument, floating point required");

  T x, y, z;

  Point<T>(T x_ = T{}, T y_ = T{}, T z_ = T{}) : x{x_}, y{y_}, z{z_} {}

  T Norm() const { return sqrt(x * x + y * y + z * z); }

  T Norm2() const { return x * x + y * y + z * z; }

  Point<T>& operator+=(Point<T> const& other) noexcept {
    x += other.x;
    y += other.y;
    z += other.z;

    return *this;
  }

  Point<T>& operator-=(Point<T> const& other) noexcept {
    x -= other.x;
    y -= other.y;
    z -= other.z;

    return *this;
  }

  template <typename K>
  Point<T>& operator*=(K scalar) noexcept {
    static_assert(
        std::is_arithmetic<K>::value,
        "Point::operator*= invalid template argument, arithmetic required");

    T const Tscalar{static_cast<T>(scalar)};

    x *= Tscalar;
    y *= Tscalar;
    z *= Tscalar;

    return *this;
  }
};

template <class from, class to>
Point<to> Convert_to(Point<to> const& other) {
  return Point<to>{static_cast<to>(other.x), static_cast<to>(other.y),
                   static_cast<to>(other.z)};
}

template <typename T>
bool operator==(Point<T> const& lhs, Point<T> const& rhs) noexcept {
  return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

template <typename T>
bool operator!=(Point<T> const& lhs, Point<T> const& rhs) noexcept {
  return !(lhs == rhs);
}

template <typename T>
Point<T> operator+(Point<T> const& lhs, Point<T> const& rhs) noexcept {
  return Point<T>{lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

template <typename T>
Point<T> operator-(Point<T> const& lhs, Point<T> const& rhs) noexcept {
  return Point<T>{lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

template <typename T>
Point<T> operator-(Point<T> const& lhs) noexcept {
  return Point<T>{-lhs.x, -lhs.y, -lhs.z};
}

template <typename T>
T operator*(Point<T> const& lhs, Point<T> const& rhs) noexcept {
  return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

template <typename T, typename K>
Point<T> operator*(Point<T> const& lhs, K scalar) noexcept {
  static_assert(
      std::is_arithmetic<K>::value,
      "Point::operator*= invalid template argument, arithmetic required");

  T const Tscalar{static_cast<T>(scalar)};

  return Point<T>{lhs.x * Tscalar, lhs.y * Tscalar, lhs.z * Tscalar};
}

template <typename T, typename K>
Point<T> operator/(Point<T> const& lhs, K scalar) noexcept {
  return lhs * (1 / scalar);
}

#endif