// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file adrc.hpp
 * @brief Active Disturbance Rejection Control design and runtime
 */

#include <cstddef>

#include "damp/backend.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"

namespace damp {

namespace design {

/**
 * @struct ADRCResult
 * @brief Active Disturbance Rejection Control design result
 *
 * Contains the computed observer and controller gains for ADRC.
 * Use .as\<U\>() to convert for type conversion (e.g., double to float).
 *
 * @see "From PID to Active Disturbance Rejection Control" (Han, 1998)
 */
template<size_t NX, typename T = double>
struct ADRCResult {
    static_assert(NX == 1 || NX == 2, "ADRC supports 1st- and 2nd-order plants only");

    T wc{}; ///< Controller bandwidth
    T wo{}; ///< Luenberger bandwidth
    T b0{}; ///< Plant gain

    damp::array<T, NX + 1> beta{}; ///< ESO gains

    T Kp{}; ///< Proportional gain
    T Kd{}; ///< Derivative gain

    bool success{false}; ///< true if wc, wo, b0 are finite and strictly positive

    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        return ADRCResult<NX, U>{
            static_cast<U>(wc),
            static_cast<U>(wo),
            static_cast<U>(b0),
            array_as<U>(beta),
            static_cast<U>(Kp),
            static_cast<U>(Kd),
            success,
        };
    }
};

/**
 * @brief Active Disturbance Rejection Control design
 *
 * Designs ESO gains using pole placement for observer poles at -wo.
 * Places all poles at -wo for the extended system.
 * Controller gains are computed based on system order:
 * - 1st order: Kp = wc/b0, Kd = 0
 * - 2nd order: Kp = wc²/b0, Kd = 2*wc/b0
 *
 * Only plant orders NX = 1 and NX = 2 are supported (same surface as typical
 * linear ADRC toolboxes / MATLAB® discrete ADRC). Higher-order ADRC is out of
 * scope.
 *
 * Requires wc > 0, wo > 0, b0 > 0 (all finite). On failure returns success=false
 * and zero gains so the runtime controller cannot divide by an invalid b0.
 *
 * @param wc Controller bandwidth
 * @param wo Luenberger bandwidth
 * @param b0 Plant gain
 *
 * @return ADRCResult with computed gains
 *
 * @see "From PID to Active Disturbance Rejection Control" (Han, 1998)
 */
template<size_t NX, typename T = double>
[[nodiscard]] constexpr ADRCResult<NX, T> adrc(T wc, T wo, T b0) {
    static_assert(NX == 1 || NX == 2, "ADRC supports 1st- and 2nd-order plants only");

    const bool ok = damp::isfinite(wc) && damp::isfinite(wo) && damp::isfinite(b0) && (wc > T{0})
                 && (wo > T{0}) && (b0 > T{0});
    if (!ok) {
        return ADRCResult<NX, T>{};
    }

    damp::array<T, NX + 1> beta{};
    constexpr size_t       n = NX + 1;
    for (size_t i = 1; i <= n; ++i) {
        T binom = T{1};
        for (size_t k = 1; k <= i; ++k) {
            binom *= static_cast<T>(n - k + 1);
            binom /= static_cast<T>(k);
        }
        beta[i - 1] = binom * damp::pow(wo, static_cast<int>(i));
    }

    if constexpr (NX == 1) {
        return ADRCResult<NX, T>{
            .wc = wc,
            .wo = wo,
            .b0 = b0,
            .beta = beta,
            .Kp = wc / b0,
            .Kd = T{0},
            .success = true,
        };
    } else {
        // NX == 2 (enforced by static_assert)
        return ADRCResult<NX, T>{
            .wc = wc,
            .wo = wo,
            .b0 = b0,
            .beta = beta,
            .Kp = (wc * wc) / b0,
            .Kd = (T{2} * wc) / b0,
            .success = true,
        };
    }
}

} // namespace design

/**
 * @ingroup discrete_controllers
 * @brief Active Disturbance Rejection Control (ADRC)
 *
 * Discrete ADRC with Extended State Luenberger (ESO) for SISO systems.
 * Supports both 1st-order and 2nd-order plant models.
 * Based on the work of Jingqing Han and Zhiqiang Gao.
 *
 * References:
 * - Han, J. (1998). From PID to Active Disturbance Rejection Control. IEEE Transactions on Industrial Electronics, 45(5), 900-906.
 * - Gao, Z. (2006). Active Disturbance Rejection Control: A Paradigm Shift in Feedback Control System Design. American Control Conference.
 *
 * @tparam NX Plant order: 1 (1st-order) or 2 (2nd-order) only
 * @tparam T Scalar type (default: float)
 */
template<size_t NX, typename T = float>
class ADRCController {
    static_assert(NX == 1 || NX == 2, "ADRC supports 1st- and 2nd-order plants only");

    template<size_t, typename>
    friend class ADRCController; // cross-precision converting ctor

    T b0{1.0f}; ///< input gain
    T Kp{1.0f}; ///< proportional gain
    T Kd{1.0f}; ///< derivative gain

    damp::array<T, NX + 1> beta{};        ///< ESO gains [β1, β2, ..., β_{NX+1}]
    ColVec<NX + 1, T>      z{};           ///< ESO state: [ŷ, ŷ̇, ..., f̂]
    T                      u_prev_{T{0}}; ///< Last applied command (post-saturation)
    T                      Ts_{T{0}};     ///< Stored sample time for the 2-arg control(r, y) form [s]
    bool                   valid_{false}; ///< From the design's success flag; gates control()

public:
    constexpr ADRCController() = default;
    constexpr ADRCController(const design::ADRCResult<NX, T>& result, T Ts = T{0})
        : b0(result.b0), Kp(result.Kp), Kd(result.Kd), beta(result.beta), Ts_(Ts), valid_(result.success && result.b0 > T{0}) {}

    template<typename U>
    constexpr ADRCController(const ADRCController<NX, U>& other)
        : b0(static_cast<T>(other.b0)),
          Kp(static_cast<T>(other.Kp)),
          Kd(static_cast<T>(other.Kd)),
          beta(array_as<T>(other.beta)),
          z(other.z.template as<T>()),
          u_prev_(static_cast<T>(other.u_prev_)),
          Ts_(static_cast<T>(other.Ts_)),
          valid_(other.valid_) {}

    /**
     * @brief Reference-tracking control step (Gao's linear ADRC).
     *
     * Runs the linear extended-state observer one Euler step using the
     * previously-applied command, then forms the control law
     *
     *     u = (Kp·(r − ẑ₁) − Kd·ẑ₂ − f̂) / b0
     *
     * where `f̂ = z_{NX+1}` is the ESO's total-disturbance estimate. Cancelling
     * `f̂` is what gives ADRC its disturbance rejection and (integral-free) zero
     * steady-state error. The two design knobs are the controller and observer
     * bandwidths `wc`/`wo` (see @ref design::adrc); there is no separate
     * integrator gain.
     *
     * The ESO uses the previously-applied command (stored in `u_prev_`). After
     * computing a new command this overload optimistically stores the
     * unsaturated value in `u_prev_`; if a downstream stage clamps the command,
     * follow up with `back_calculate(u_unsat, u_sat)` so the next ESO tick uses
     * what was actually applied.
     *
     * Failed designs (success=false or b0 ≤ 0) emit zero and hold the ESO state
     * so the control law never divides by an invalid b0.
     *
     * @param r  Reference (setpoint).
     * @param y  Measurement (plant output).
     * @param Ts Sample time [s] for the explicit-Euler ESO update.
     * @return Control command `u`.
     */
    [[nodiscard]] constexpr T control(T r, T y, T Ts) {
        // Invalid design or non-positive b0 would divide-by-zero in the control law.
        // Hold the ESO and emit no correction (same gate pattern as SMC).
        if (!valid_ || !(b0 > T{0})) {
            return T{0};
        }

        // --- Linear ESO, explicit-Euler update using the last applied command.
        // Plant model: y^(NX) = f + b0·u, with z = [ŷ, ŷ̇, …, ŷ^(NX-1), f̂].
        const T e = z[0] - y; // estimation error (ŷ − y)

        damp::array<T, NX + 1> dz{};
        for (size_t i = 0; i < NX; ++i) {
            dz[i] = z[i + 1] - (beta[i] * e); // ż_i = ẑ_{i+1} − β_i·e
        }
        dz[NX - 1] += b0 * u_prev_; // b0·u enters the highest derivative state
        dz[NX] = -(beta[NX] * e);   // ḟ̂ = −β_{NX+1}·e

        for (size_t i = 0; i <= NX; ++i) {
            z[i] += Ts * dz[i];
        }

        // --- Control law: PD on the estimated state, minus the disturbance.
        const T u0 = (Kp * (r - z[0])) - (Kd * z[1]); // Kd == 0 for NX == 1
        const T u = (u0 - z[NX]) / b0;                // cancel f̂ = z_{NX+1}

        u_prev_ = u;
        return u;
    }

    /**
     * @brief 2-arg fixed-rate form using the stored sample time.
     *
     * The SisoController concept shape: construct with Ts and call
     * control(r, y) each tick. With no stored Ts the ESO holds (no state
     * advance), exactly like the 3-arg form with Ts = 0.
     */
    [[nodiscard]] constexpr T control(T r, T y) { return control(r, y, Ts_); }

    /**
     * @brief Anti-windup hook for cascade-level saturation propagation.
     *
     * Overrides the stored `u_prev_` with the actually-realized command so
     * the next ESO tick reflects what the plant received, not what the
     * controller wanted. This is the natural ADRC anti-windup behavior:
     * the ESO's disturbance estimate self-corrects once it sees the true
     * input, without an explicit integrator unwind step.
     */
    constexpr void back_calculate(T u_unsat, T u_sat) {
        (void)u_unsat;
        u_prev_ = u_sat;
    }

    [[nodiscard]] constexpr bool valid() const { return valid_; }

    constexpr void reset() {
        z = ColVec<NX + 1, T>{};
        u_prev_ = T{0};
    }
};
} // namespace damp
