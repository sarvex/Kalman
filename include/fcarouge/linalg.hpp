/*  __          _      __  __          _   _
| |/ /    /\   | |    |  \/  |   /\   | \ | |
| ' /    /  \  | |    | \  / |  /  \  |  \| |
|  <    / /\ \ | |    | |\/| | / /\ \ | . ` |
| . \  / ____ \| |____| |  | |/ ____ \| |\  |
|_|\_\/_/    \_\______|_|  |_/_/    \_\_| \_|

Kalman Filter
Version 0.2.0
https://github.com/FrancoisCarouge/Kalman

SPDX-License-Identifier: Unlicense

This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <https://unlicense.org> */

#ifndef FCAROUGE_LINALG_HPP
#define FCAROUGE_LINALG_HPP

#include "utility.hpp"

#include <concepts>
#include <cstddef>
#include <type_traits>

namespace fcarouge::linalg {

//! @brief Static matrix.
template <typename Type = double, auto Row = 1, auto Column = 1> class matrix {
public:
  inline constexpr matrix() = default;
  // inline constexpr matrix(matrix &&);

  inline constexpr matrix(const matrix &other) {
    for (auto i{0}; i < Row; ++i) {
      for (auto j{0}; j < Column; ++j) {
        mda[i][j] = other(i, j);
      }
    }
  }

  // inline constexpr ~matrix() = default;
  // inline constexpr auto operator=(matrix &&) -> matrix &;
  inline constexpr auto operator=(const matrix &other) -> matrix & {
    for (auto i{0}; i < Row; ++i) {
      for (auto j{0}; j < Column; ++j) {
        mda[i][j] = other(i, j);
      }
    }
    return *this;
  }

  // inline constexpr auto data() -> Type * { return &mda[0][0]; }

  inline constexpr explicit(false)
      matrix(Type value) requires(Row == 1 && Column == 1) {
    mda[0][0] = value;
  }

  inline constexpr explicit matrix(Type(values)[Column]) {
    static_cast<void>(values);
  }

  // inline constexpr explicit matrix(const auto &value) {
  //   static_cast<void>(value);
  // }

  [[maybe_unused]] inline constexpr const Type &operator[](auto index) const
      requires(Column == 1) {
    return mda[index][0];
  }

  [[maybe_unused]] inline constexpr const Type &operator[](auto index) const
      requires(Row == 1 && Column > 1) {
    return mda[0][index];
  }

  [[maybe_unused]] inline constexpr const Type &operator()(auto row,
                                                           auto column) const {
    return mda[row][column];
  }

  [[maybe_unused]] inline constexpr Type &operator()(auto row, auto column) {
    return mda[row][column];
  }

  friend bool operator==(const matrix &lhs, const matrix &rhs) = default;

private:
  [[no_unique_address]] Type mda[Row][Column]{};
};

//! @name Deduction Guides
//! @{

template <typename Type, auto Row, auto Column>
matrix(const Type (&)[Row][Column]) -> matrix<Type, Row, Column>;

template <typename Type, auto Row>
matrix(const Type (&)[Row]) -> matrix<Type, Row, 1>;

//! @name Vector Types
//! @{

//! @brief Row vector.
template <typename Type = double, auto Column = 1>
using row_vector = matrix<Type, decltype(Column){1}, Column>;

//! @brief Column vector.
template <typename Type = double, auto Row = 1>
using column_vector = matrix<Type, Row, decltype(Row){1}>;

//! @}

template <typename Type, auto M, auto N, auto O>
[[nodiscard]] inline constexpr auto operator*(const matrix<Type, M, N> &lhs,
                                              const matrix<Type, N, O> &rhs) {
  matrix<Type, M, O> m{};
  static_cast<void>(lhs);
  static_cast<void>(rhs);
  return m;
}

[[nodiscard]] inline constexpr auto
operator*(float lhs, const fcarouge::linalg::matrix<float, 1, 1> rhs) {
  float a{};
  static_cast<void>(lhs);
  static_cast<void>(rhs);
  return a;
}

template <typename Type, auto Column>
[[nodiscard]] inline constexpr auto operator*(arithmetic auto lhs,
                                              matrix<Type, 1, Column> rhs) {
  matrix<Type, 1, Column> m{};
  static_cast<void>(lhs);
  static_cast<void>(rhs);
  return m;
}

template <typename Type, auto M, auto N>
[[nodiscard]] inline constexpr auto operator*(matrix<Type, M, N> lhs,
                                              arithmetic auto rhs) {
  matrix<Type, M, N> m{};
  static_cast<void>(lhs);
  static_cast<void>(rhs);
  return m;
}

template <typename Type, auto M, auto N>
[[nodiscard]] inline constexpr auto operator+(const matrix<Type, M, N> &lhs,
                                              const matrix<Type, M, N> &rhs) {
  matrix<Type, M, N> m{};
  for (auto i{0}; i < M; ++i) {
    for (auto j{0}; j < N; ++j) {
      static_cast<void>(lhs);
      static_cast<void>(rhs);
    }
  }

  return m;
}

template <typename Type>
[[nodiscard]] inline constexpr auto operator+(arithmetic auto lhs,
                                              matrix<Type, 1, 1> rhs) {
  Type a{};
  static_cast<void>(lhs);
  static_cast<void>(rhs);
  return a;
}

template <typename Type>
[[nodiscard]] inline constexpr auto operator+(matrix<Type, 1, 1> lhs,
                                              arithmetic auto rhs) {
  Type a{};
  static_cast<void>(lhs);
  static_cast<void>(rhs);
  return a;
}

template <typename Type, auto M, auto N>
[[nodiscard]] inline constexpr auto operator-(const matrix<Type, M, N> &lhs,
                                              const matrix<Type, M, N> &rhs) {
  matrix<Type, M, N> m{};
  for (auto i{0}; i < M; ++i) {
    for (auto j{0}; j < N; ++j) {
      static_cast<void>(lhs);
      static_cast<void>(rhs);
    }
  }

  return m;
}

template <typename Type>
[[nodiscard]] inline constexpr auto operator-(arithmetic auto lhs,
                                              const matrix<Type, 1, 1> &rhs) {
  Type m{};
  static_cast<void>(lhs);
  static_cast<void>(rhs);
  return m;
}

template <typename Type, auto Row>
[[nodiscard]] inline constexpr auto operator/(const matrix<Type, Row, 1> &lhs,
                                              arithmetic auto rhs) {
  matrix<Type, Row, 1> m{};
  static_cast<void>(lhs);
  static_cast<void>(rhs);
  return m;
}

template <typename Type, auto Row1, auto Column1, auto Row2, auto Column2>
[[nodiscard]] inline constexpr auto
operator/(const matrix<Type, Row1, Column1> &lhs,
          const matrix<Type, Row2, Column2> &rhs) {
  matrix<Type, Row1, Row2> m{};
  static_cast<void>(lhs);
  static_cast<void>(rhs);
  return m;
}

template <typename Type, auto Row, auto Column>
[[nodiscard]] inline constexpr auto
transpose(const matrix<Type, Row, Column> &lhs) {
  matrix<Type, Column, Row> m{};
  static_cast<void>(lhs);
  return m;
}

template <typename Type = double> inline constexpr Type identity_v{1.0};
template <typename Type>
inline constexpr Type identity_v<matrix<Type, 1, 1>>{1.0};
template <typename Type, auto Row, auto Column>
inline constexpr matrix<Type, Row, Column>
    identity_v<matrix<Type, Row, Column>>{[] {
      matrix<Type, Row, Column> m;
      // Type &data{m.data()};
      auto max{Row < Column ? Row : Column};
      for (auto ij{0}; ij < max; ++ij) {
        // data[ij][ij] = 1.0;
        m(ij, ij) = 1.0;
      }
      return m;
    }()};

template <typename Type = double> inline constexpr Type zero_v{};
template <typename Type> inline constexpr Type zero_v<matrix<Type, 1, 1>>{};

} // namespace fcarouge::linalg

#endif // FCAROUGE_LINALG_HPP
