/*  __          _      __  __          _   _
| |/ /    /\   | |    |  \/  |   /\   | \ | |
| ' /    /  \  | |    | \  / |  /  \  |  \| |
|  <    / /\ \ | |    | |\/| | / /\ \ | . ` |
| . \  / ____ \| |____| |  | |/ ____ \| |\  |
|_|\_\/_/    \_\______|_|  |_/_/    \_\_| \_|

Kalman Filter for C++
Version 0.1.0
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

#include "fcarouge/eigen/kalman.hpp"

#include <Eigen/Eigen>

#include <cassert>

namespace fcarouge::eigen::test {
namespace {

template <typename Type, auto Size> using vector = Eigen::Vector<Type, Size>;

template <typename Type, auto RowSize, auto ColumnSize>
using matrix = Eigen::Matrix<Type, RowSize, ColumnSize>;

//! @test Verifies the state transition matrix F management overloads for
//! the Eigen filter type.
[[maybe_unused]] auto f_5x4x3{[] {
  using kalman = kalman<vector<double, 5>, vector<double, 4>, vector<double, 3>,
                        std::tuple<double, float, int, char>,
                        std::tuple<char, int, float, double>>;

  kalman filter;

  const auto i5x5{matrix<double, 5, 5>::Identity()};
  const auto z5x5{matrix<double, 5, 5>::Zero()};
  const vector<double, 3> z3{vector<double, 3>::Zero()};

  assert(filter.f() == i5x5);

  {
    const auto f{i5x5};
    filter.f(f);
    assert(filter.f() == i5x5);
  }

  {
    const auto f{z5x5};
    filter.f(std::move(f));
    assert(filter.f() == z5x5);
  }

  {
    const auto f{i5x5};
    filter.f(f);
    assert(filter.f() == i5x5);
  }

  {
    const auto f{z5x5};
    filter.f(std::move(f));
    assert(filter.f() == z5x5);
  }

  {
    const auto f{[](const kalman::state &x, const kalman::input &u,
                    const char &c, const int &i, const float &fp,
                    const double &d) -> kalman::state_transition {
      static_cast<void>(x);
      static_cast<void>(d);
      static_cast<void>(fp);
      static_cast<void>(i);
      static_cast<void>(c);
      static_cast<void>(u);
      return matrix<double, 5, 5>::Identity();
    }};
    filter.f(f);
    assert(filter.f() == z5x5);
    filter.predict(char(0), 0, 0.f, 0., z3);
    assert(filter.f() == i5x5);
  }

  {
    const auto f{[](const kalman::state &x, const kalman::input &u,
                    const char &c, const int &i, const float &fp,
                    const double &d) -> kalman::state_transition {
      static_cast<void>(x);
      static_cast<void>(d);
      static_cast<void>(fp);
      static_cast<void>(i);
      static_cast<void>(c);
      static_cast<void>(u);
      return matrix<double, 5, 5>::Zero();
    }};
    filter.f(std::move(f));
    assert(filter.f() == i5x5);
    filter.predict(0, 0, 0.f, 0., z3);
    assert(filter.f() == z5x5);
  }

  return 0;
}()};

} // namespace
} // namespace fcarouge::eigen::test