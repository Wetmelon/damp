// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file pade.hpp
 * @brief Padé approximations of a pure time delay, continuous and Tustin-discretized
 */

#include "damp/design/model_result.hpp"
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/transfer_function.hpp"
#include "damp/systems/zpk.hpp"

namespace damp {
namespace design {

/**
 * @brief First-order Pade approximation of time delay
 *
 * Approximates e^(-sT) with a first-order rational transfer function:
 * H(s) = (1 - sT/2) / (1 + sT/2)
 *
 * @param T_delay Time delay [s] (> 0, finite)
 * @return TransferFunctionResult<2, 2, T> with the delay approximation; refused (invalid_parameter) unless 0 < T_delay < ∞
 */
template<typename T = double>
[[nodiscard]] constexpr TransferFunctionResult<2, 2, T> pade_delay_1st(T T_delay) {
    if (!(T_delay > T{0}) || !damp::isfinite(T_delay)) {
        return ::damp::detail::refused<TransferFunctionResult<2, 2, T>>(DesignStatus::invalid_parameter);
    }
    const T half_T = T_delay / T{2};
    return detail::designed(TransferFunction<2, 2, T>{
        {-half_T, T{1}}, // num: −(T/2)s + 1
        {half_T, T{1}}   // den: (T/2)s + 1
    });
}

/**
 * @brief First-order Pade approximation of time delay (discrete-time)
 *
 * Discrete-time version of the first-order Pade delay approximation.
 *
 * @param T_delay Time delay [s]
 * @param Ts Sample time [s]
 * @return TransferFunctionResult<2, 2, T> with the discrete-time delay approximation; refused with
 *         invalid_sample_time if Ts ≤ 0, invalid_parameter unless 0 < T_delay < ∞, singular if the
 *         realization fails
 */
template<typename T = float>
[[nodiscard]] constexpr TransferFunctionResult<2, 2, T> pade_delay_1st(T T_delay, T Ts) {
    if (!(Ts > T{0}) || !damp::isfinite(Ts)) {
        return ::damp::detail::refused<TransferFunctionResult<2, 2, T>>(DesignStatus::invalid_sample_time);
    }
    const auto tf_c = pade_delay_1st<T>(T_delay);
    if (!tf_c.success) {
        return tf_c;
    }
    const auto sys_c_opt = tf_c.tf.to_state_space();
    if (!sys_c_opt) {
        return ::damp::detail::refused<TransferFunctionResult<2, 2, T>>(DesignStatus::singular);
    }
    const auto sys_d_opt = discretize(*sys_c_opt, Ts, DiscretizationMethod::Tustin);
    if (!sys_d_opt) {
        return ::damp::detail::refused<TransferFunctionResult<2, 2, T>>(DesignStatus::singular);
    }
    return detail::designed(::damp::to_transfer_function(*sys_d_opt));
}

/**
 * @brief Second-order Pade approximation of time delay
 *
 * Approximates e^(-sT) with a second-order rational transfer function:
 * H(s) = (1 - sT/2 + (sT)²/12) / (1 + sT/2 + (sT)²/12)
 *
 * @param T_delay Time delay [s] (> 0, finite)
 * @return TransferFunctionResult<3, 3, T> with the delay approximation; refused (invalid_parameter) unless 0 < T_delay < ∞
 */
template<typename T = double>
[[nodiscard]] constexpr TransferFunctionResult<3, 3, T> pade_delay_2nd(T T_delay) {
    if (!(T_delay > T{0}) || !damp::isfinite(T_delay)) {
        return ::damp::detail::refused<TransferFunctionResult<3, 3, T>>(DesignStatus::invalid_parameter);
    }
    const T half_T = T_delay / T{2};
    const T T_sq_12 = T_delay * T_delay / T{12};
    return detail::designed(TransferFunction<3, 3, T>{
        {T_sq_12, -half_T, T{1}}, // num: (T²/12)s² − (T/2)s + 1
        {T_sq_12, half_T, T{1}}   // den: (T²/12)s² + (T/2)s + 1
    });
}

/**
 * @brief Second-order Pade approximation of time delay (discrete-time)
 *
 * Discrete-time version of the second-order Pade delay approximation.
 *
 * @param T_delay Time delay [s]
 * @param Ts Sample time [s]
 * @return TransferFunctionResult<3, 3, T> with the discrete-time delay approximation; refused with
 *         invalid_sample_time if Ts ≤ 0, invalid_parameter unless 0 < T_delay < ∞, singular if the
 *         realization fails
 */
template<typename T = float>
[[nodiscard]] constexpr TransferFunctionResult<3, 3, T> pade_delay_2nd(T T_delay, T Ts) {
    if (!(Ts > T{0}) || !damp::isfinite(Ts)) {
        return ::damp::detail::refused<TransferFunctionResult<3, 3, T>>(DesignStatus::invalid_sample_time);
    }
    const auto tf_c = pade_delay_2nd<T>(T_delay);
    if (!tf_c.success) {
        return tf_c;
    }
    const auto sys_c_opt = tf_c.tf.to_state_space();
    if (!sys_c_opt) {
        return ::damp::detail::refused<TransferFunctionResult<3, 3, T>>(DesignStatus::singular);
    }
    const auto sys_d_opt = discretize(*sys_c_opt, Ts, DiscretizationMethod::Tustin);
    if (!sys_d_opt) {
        return ::damp::detail::refused<TransferFunctionResult<3, 3, T>>(DesignStatus::singular);
    }
    return detail::designed(::damp::to_transfer_function(*sys_d_opt));
}

} // namespace design
} // namespace damp
