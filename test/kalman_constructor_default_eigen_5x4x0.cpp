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

#include "fcarouge/kalman.hpp"

#include <Eigen/Eigen>

#include <cassert>

namespace fcarouge::test {
namespace {

template <auto Size> using vector = Eigen::Vector<double, Size>;

template <auto Row, auto Column>
using matrix = Eigen::Matrix<double, Row, Column>;

//! @test Verifies default values are initialized for multi-dimension filters,
//! no input.
[[maybe_unused]] auto test{[] {
  using kalman = kalman<vector<5>, vector<4>>;
  kalman filter;

  const auto i4x4{matrix<4, 4>::Identity()};
  const auto i4x5{matrix<4, 5>::Identity()};
  const auto i5x4{matrix<5, 4>::Identity()};
  const auto i5x5{matrix<5, 5>::Identity()};
  const auto z4x1{vector<4>::Zero()};
  const auto z4x4{matrix<4, 4>::Zero()};
  const auto z5x1{vector<5>::Zero()};
  const auto z5x5{matrix<5, 5>::Zero()};

  assert(filter.f() == i5x5);
  assert(filter.h() == i4x5);
  assert(filter.k() == i5x4);
  assert(filter.p() == i5x5);
  assert(filter.q() == z5x5 && "No process noise by default.");
  assert(filter.r() == z4x4 && "No observation noise by default.");
  assert(filter.s() == i4x4);
  assert(filter.x() == z5x1 && "Origin state.");
  assert(filter.y() == z4x1);
  assert(filter.z() == z4x1);

  return 0;
}()};

} // namespace
} // namespace fcarouge::test
