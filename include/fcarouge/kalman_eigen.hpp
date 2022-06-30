/*_  __          _      __  __          _   _
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

#ifndef FCAROUGE_KALMAN_EIGEN_HPP
#define FCAROUGE_KALMAN_EIGEN_HPP

//! @file
//! @brief Kalman operation for Eigen 3 types.

#include "kalman.hpp"
#include "internal/kalman_eigen_operator.hpp"

#include <Eigen/Eigen>

#include <cstddef>
#include <functional>
#include <tuple>

namespace fcarouge::eigen
{
//! @brief Eigen-based Kalman filter.
//!
//! @details Implemented with the Eigen linear algebra library matrices with
//! sizes fixed at compile-time.
//!
//! @tparam Type The type template parameter of the matrices data.
//! @tparam State The non-type template parameter size of the state vector x.
//! @tparam Output The non-type template parameter size of the measurement
//! vector z.
//! @tparam Input The non-type template parameter size of the control u.
//! @tparam PredictionArguments The variadic type template parameter for
//! additional prediction function parameters. Time, or a delta thereof, is
//! often a prediction parameter.
template <typename Type = double, std::size_t State = 1, std::size_t Output = 1,
          std::size_t Input = 0, typename UpdateArguments = std::tuple<>,
          typename PredictionArguments = std::tuple<>>
using kalman = fcarouge::kalman<
    Type, Eigen::Vector<Type, State>, Eigen::Vector<Type, Output>,
    Eigen::Vector<Type, Input>, internal::transpose, internal::symmetrize,
    internal::divide, internal::identity, UpdateArguments, PredictionArguments>;

} // namespace fcarouge::eigen

#endif // FCAROUGE_KALMAN_EIGEN_HPP