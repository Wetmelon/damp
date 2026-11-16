// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file state_feedback_types.hpp
 * @brief The types of an LQ design: what it is given (LQRConfig, LQIConfig, IntegralAction) and what it hands
 *        its runtime (StateFeedbackResult, NoEstimator)
 *
 * The designs are in lqr.hpp, lqi.hpp, lqg.hpp, lqgi.hpp and oflqg.hpp (detail/state_feedback_design.hpp
 * behind them), the runtime they share in state_feedback_core.hpp.
 */

#include "damp/backend.hpp"
#include "damp/controllers/law_detail.hpp"
#include "damp/design/status.hpp"
#include "damp/design/steady_state_target.hpp"
#include "damp/estimation/fixed_gain_observer.hpp"
#include "damp/math/complex.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/toolbox/bounds.hpp"

namespace damp {
namespace design {

/** @brief The filter half of an estimated design (defined in kalman.hpp). */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t ND>
struct KalmanResult;

} // namespace design

/**
 * @brief The estimator of a law that measures its state: none, an empty type, so a
 *        `[[no_unique_address]]` member of it takes no storage.
 */
struct NoEstimator {
    /** @brief Any scalar: still nothing. */
    template<typename U>
    [[nodiscard]] constexpr NoEstimator as() const {
        return {};
    }
};

namespace detail {

/** @brief What an estimator type carries: none (NoEstimator). */
template<typename E>
struct estimator_traits {
    static constexpr bool   present = false; ///< A filter runs
    static constexpr size_t NXE = 0;         ///< Filter states before the disturbance
    static constexpr size_t ND = 0;          ///< Disturbance states
};

/** @brief A FixedGainObserver's sizes. */
template<size_t NXK, size_t NU, size_t NY, typename T, size_t ND_, typename Store>
struct estimator_traits<FixedGainObserver<NXK, NU, NY, T, ND_, Store>> {
    static constexpr bool   present = true;  ///< A filter runs
    static constexpr size_t NXE = NXK - ND_; ///< Filter states before the disturbance
    static constexpr size_t ND = ND_;        ///< Disturbance states
};

/** @brief The runtime of a design's filter half: none for a measured state. */
template<typename Kalman, typename T>
struct runtime_estimator {
    using type = NoEstimator; ///< NoEstimator
};

/** @brief A KalmanResult runs as the FixedGainObserver of its sizes (its fixed gain M), in the runtime's scalar @p T. */
template<size_t NXK, size_t NU, size_t NY, typename U, size_t NW, size_t NV, size_t ND, typename T>
struct runtime_estimator<design::KalmanResult<NXK, NU, NY, U, NW, NV, ND>, T> {
    using type = FixedGainObserver<NXK, NU, NY, T, ND>; ///< The observer and the disturbance
};

/** @brief The runtime of a design's filter half (@p Kalman: KalmanResult or NoEstimator), in the scalar @p T. */
template<typename Kalman, typename T>
using runtime_estimator_t = runtime_estimator<Kalman, T>::type;

} // namespace detail

namespace design {

/**
 * @brief Offset-free by integral action: integrators on r − y, every measured output
 *        (`IntegralAction{}`) or the listed ones (`IntegralAction<2>{{0, 3}}`, in that order:
 *        the references r follow it).
 * @tparam NI Outputs integrated, listed in @ref outputs (0: every measured output)
 */
template<size_t NI = 0>
struct IntegralAction {
    damp::array<size_t, NI> outputs{}; ///< Indices of the integrated outputs, each once (none listed: all)
};

/**
 * @struct LQRConfig
 * @brief What a regulator design (lqr, lqg, oflqg) is given with the plant: the cost, what to
 *        penalize and how much, and the limit its command runs within
 *
 * The cost is the weights Q, R and N, in the plant's own time domain: a cost integral on a
 * continuous plant (with a Ts argument it is discretized exactly), per step on a discrete one.
 * Braces size it from the call: `{.Q = Q, .R = R}`. A law with integral action takes LQIConfig.
 *
 * | Field | Penalizes | Size |
 * | --- | --- | --- |
 * | Q | the plant states x | NX × NX, positive semidefinite |
 * | R | the inputs u | NU × NU, positive definite |
 * | N | the cross term 2xᵀNu (default zero) | NX × NU |
 *
 * An R that is not positive definite is refused with DesignStatus::bad_weights.
 *
 * u_lim is the runtime's, not the cost's: the gain is the same with or without it.
 * `{.Q = Q, .R = R, .u_lim = {.lower = -12.0, .upper = 12.0}}` clamps the command control
 * returns to ±12, per input (default unbounded); back_calculate stays for a stage that clamps
 * after the law. It is a member of the runtime, `law.u_lim`, so a rail that moves is set at run
 * time. A limit that is NaN, or a lower limit above its upper one, is refused with
 * DesignStatus::invalid_parameter.
 */
template<size_t NX, size_t NU, typename T = double>
struct LQRConfig {
    Matrix<NX, NX, T> Q{};     ///< On the plant states x (positive semidefinite)
    Matrix<NU, NU, T> R{};     ///< On the inputs u (positive definite)
    Matrix<NX, NU, T> N{};     ///< Cross term between x and u (default zero)
    Bounds<NU, T>     u_lim{}; ///< Limits on the command u, per input (default unbounded)
};

/** @brief Sizes from the matrices: `LQRConfig{Q, R}`. */
template<size_t NX, size_t NU, typename T>
LQRConfig(Matrix<NX, NX, T>, Matrix<NU, NU, T>) -> LQRConfig<NX, NU, T>;

/** @brief Sizes from the matrices: `LQRConfig{Q, R, N}`. */
template<size_t NX, size_t NU, typename T>
LQRConfig(Matrix<NX, NX, T>, Matrix<NU, NU, T>, Matrix<NX, NU, T>) -> LQRConfig<NX, NU, T>;

/**
 * @struct LQIConfig
 * @brief What a design with integral action (lqi, lqgi) is given with the plant: LQRConfig's
 *        cost and limit, and the weight and the limits on the NI integrators
 *
 * Q, R, N and u_lim are LQRConfig's. Braces size it from the call:
 * `{.Q = Q, .R = R, .Q_integral = Qi}`. Template order is NX, NU, NI, T: `LQIConfig<2, 1, 1>`
 * is one integrator on double.
 *
 * | Field | | Size |
 * | --- | --- | --- |
 * | Q_integral | penalizes the discrete integral state ξ | NI × NI, positive definite |
 * | integral_lim | limits ξ, per integrator | NI, default unbounded |
 *
 * Q_integral is the cost on the Euler-integrated error ξ ≈ ∫e dt, not a block of a generic
 * augmented Q. Q, R and N are the plant cost (Van Loan on a continuous plant); Q_integral keeps
 * that meaning on both the continuous and discrete routes. The stacked DARE cost is
 * block-diagonal in Q and Q_integral: x–ξ cross terms are not a public weight and are
 * refused. The closed loop still couples x and ξ, so raising Q_integral
 * speeds the integral and moves the regulator poles too. An integrator Q_integral does not
 * weight is a mode at z = 1 the cost does not see; a Q_integral that is not positive definite
 * (left out of the braces, for one) is refused with DesignStatus::bad_weights.
 *
 * The law holds its integrators against the clamp of u_lim, with no back_calculate.
 * integral_lim is on ξ itself, in the units of ∫e dt (the error's times seconds): through the
 * gain, ξ moves input j by K(j, NX + i)·ξᵢ. It is a member of the runtime, `law.integral.lim`,
 * set at run time like `law.u_lim`, and refused as u_lim is.
 */
template<size_t NX, size_t NU, size_t NI, typename T = double>
struct LQIConfig {
    Matrix<NX, NX, T> Q{};            ///< On the plant states x (positive semidefinite)
    Matrix<NU, NU, T> R{};            ///< On the inputs u (positive definite)
    Matrix<NX, NU, T> N{};            ///< Cross term between x and u (default zero)
    Matrix<NI, NI, T> Q_integral{};   ///< On the discrete integral state ξ ≈ ∫e dt (positive definite)
    Bounds<NU, T>     u_lim{};        ///< Limits on the command u, per input (default unbounded)
    Bounds<NI, T>     integral_lim{}; ///< Limits on the integral state ξ, per integrator (default unbounded)
};

/** @brief Sizes from the matrices: `LQIConfig{Q, R, N, Q_integral}`. */
template<size_t NX, size_t NU, typename T, size_t NI>
LQIConfig(Matrix<NX, NX, T>, Matrix<NU, NU, T>, Matrix<NX, NU, T>, Matrix<NI, NI, T>) -> LQIConfig<NX, NU, NI, T>;

namespace detail {

/** @brief The config of a law with @p NI integrators: LQRConfig for none, LQIConfig with them. */
template<size_t NX, size_t NU, size_t NI, typename T>
using lq_config = damp::conditional_t<(NI == 0), LQRConfig<NX, NU, T>, LQIConfig<NX, NU, NI, T>>;

} // namespace detail

/**
 * @brief A limit on what a law's outer terms ask of one of its states: a cascade's outer limit,
 *        inside one gain
 *
 * Read the law u = u* − K·z, z = [x − x*; ξ], as a cascade on the state s: u = −K_s·(z_s − z_c) − (the
 * other terms), where the outer terms (the plant states and integrators listed in @ref outer) ask for
 * z_c = −(Σ K_j·z_j)/K_s. The law clamps z_c into @ref lim, per input that acts on s (K_s ≠ 0). A
 * balancing robot's position, speed and integral ask for a lean; this bounds the lean they may ask
 * for, not the lean itself, which no law without a prediction can bound. The integrators hold against
 * the clamp as against u_lim's. At a steady state the outer terms carry the input that holds it (none
 * for a robot balanced upright; a spring's force through its position): a limit below that holds the
 * plant short of its reference.
 *
 * Off by default (state = none). Set it on the runtime: `law.demand_lim = {.state = 2, .outer = {true,
 * true, false, false, true}, .lim = {.lower = -0.26, .upper = 0.26}}`.
 *
 * @tparam N The gain's plant states and integrators, NX + NI (the input registers are never outer)
 */
template<size_t N, typename T = double>
struct DemandLimit {
    static constexpr size_t none = static_cast<size_t>(-1); ///< No state: the limit is off

    size_t               state{none}; ///< s: the state asked for, an index of the plant state x (none: off)
    damp::array<bool, N> outer{};     ///< The terms that ask: the plant states, then the integrators
    Bounds<1, T>         lim{};       ///< The most they may ask of x_s, as a deviation from the reference point (the target with one)
};

/** @brief The integral action of a design: the integrated outputs, their maps of the state and the input, the per-sample step and the limits on ξ. */
template<size_t NI, size_t NX, size_t NU, size_t NY, typename T = double>
struct IntegralDesign {
    Matrix<NI, NY, T> Ci{};    ///< The integrated outputs: e = r − Ci·y, y = C·x + D·u
    Matrix<NI, NX, T> Ci_x{};  ///< Ci·C: the integrated outputs as a map of the plant state (what a law that measures its state integrates)
    Matrix<NI, NU, T> Ci_u{};  ///< Ci·D: the integrated outputs' feedthrough of the input applied over the period
    T                 step{1}; ///< Per-sample integrator gain: the design's Ts (xi += Ts·e ≈ ∫e dt)
    Bounds<NI, T>     lim{};   ///< Limits on the integral state ξ (LQIConfig::integral_lim; default unbounded)

    /** @brief The same integral action in the scalar U. */
    template<typename U>
    [[nodiscard]] constexpr IntegralDesign<NI, NX, NU, NY, U> as() const {
        return {Ci.template as<U>(), Ci_x.template as<U>(), Ci_u.template as<U>(), static_cast<U>(step), lim.template as<U>()};
    }
};

/**
 * @struct StateFeedbackResult
 * @brief A state-feedback design: the gain, and the parts the law's options add
 *
 * Its template arguments are the runtime's: LQR, LQI, LQG, LQGI or OFLQG, the one its design names, deduces them from it. The
 * gain is on the law's state [x; xi]: the plant, then the integrators (NI > 0).
 *
 * | Member | What it is | Present |
 * | --- | --- | --- |
 * | K, S, e | the gain, the Riccati solution, the closed-loop poles | always |
 * | success, status, Ts | ok or why not; the law's sample time | always |
 * | u_lim | the limits on the command (LQRConfig::u_lim) | always, default unbounded |
 * | kalman | the Kalman filter (design::kalman's result) | lqg, lqgi, oflqg |
 * | integral | the integrated outputs Ci, Ci·C, Ci·D, the step and the limits on ξ | NI > 0 |
 * | target | the steady-state target (design::steady_state_target's result) | oflqg; lqr and lqg with TrackedOutputs (NZ > 0) |
 *
 * A design on a plant whose input arrives late (an InputDelay) is on the law's model
 * [x; v(k−1); …; v(k−NRU)]: K, S, e, the integrated outputs, the filter's model and the target
 * cover the registers. A measured-state runtime holds them, a filter carries them.
 *
 * @tparam NX     Plant states
 * @tparam NI     Integrators (0: none)
 * @tparam NY     Measured outputs (NX for a law without integrators or estimator)
 * @tparam Kalman The filter's design: KalmanResult, or NoEstimator for a measured state
 * @tparam NZ     Controlled outputs of the target (0: none)
 * @tparam NRU    Input registers the law holds, NRU·NU states after the plant's (0: none)
 *
 * @see "Optimal Control" (Anderson & Moore, 1990), §4 and §8
 */
template<size_t NX, size_t NU, typename T = double, size_t NI = 0, size_t NY = NX, typename Kalman = NoEstimator, size_t NZ = 0, size_t NRU = 0>
struct [[nodiscard]] StateFeedbackResult {
    using Filter = ::damp::detail::estimator_traits<::damp::detail::runtime_estimator_t<Kalman, T>>;
    static constexpr bool   kEstimated = Filter::present; ///< A Kalman filter estimates the state
    static constexpr bool   kIntegral = NI > 0;           ///< The law has integral action
    static constexpr size_t NXL = NX + (NRU * NU);        ///< The law's model: the plant, then its input registers
    static constexpr size_t NK = NXL + NI;                ///< The gain's state: the law's model and the integrators
    static constexpr size_t ND = Filter::ND;              ///< Disturbance states of the filter
    using Integral = ::damp::detail::maybe<(NI > 0), IntegralDesign<NI, NXL, NU, NY, T>, ::damp::detail::OffSlot::design_integral>;
    using Target = ::damp::detail::maybe<(NZ > 0), SteadyStateTargetResult<NXL, NU, ND, NZ, T>, ::damp::detail::OffSlot::design_target>;

    Matrix<NU, NK, T>                 K{};                                 ///< Optimal gain: u = −K·[x; xi]
    Matrix<NK, NK, T>                 S{};                                 ///< Riccati solution (positive semidefinite)
    damp::array<damp::complex<T>, NK> e{};                                 ///< Closed-loop poles of the law's model
    bool                              success{false};                      ///< status == DesignStatus::ok
    DesignStatus                      status{DesignStatus::not_converged}; ///< ok, or the first failure: the gain's, the filter's, the target's
    T                                 Ts{0};                               ///< Sample time of the law: 0 continuous (e in the s-plane), > 0 discrete
    [[no_unique_address]] Kalman      kalman{};                            ///< The Kalman filter (lqg, lqgi, oflqg)
    [[no_unique_address]] Integral    integral{};                          ///< The integral action (NI > 0)
    [[no_unique_address]] Target      target{};                            ///< The steady-state target (NZ > 0)
    Bounds<NU, T>                     u_lim{};                             ///< Limits on the command the runtime returns (default unbounded)

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /** @brief The same design in the scalar @p U: `res.as<float>()` for the embedded runtime. */
    template<typename U>
    [[nodiscard]] constexpr auto as() const {
        using Out = StateFeedbackResult<NX, NU, U, NI, NY, damp::remove_cvref_t<decltype(kalman.template as<U>())>, NZ, NRU>;
        Out out{K.template as<U>(), S.template as<U>(), array_as<damp::complex<U>>(e), success, status, static_cast<U>(Ts), kalman.template as<U>()};
        if constexpr (NI > 0) {
            out.integral = integral.template as<U>();
        }
        if constexpr (NZ > 0) {
            out.target = target.template as<U>();
        }
        out.u_lim = u_lim.template as<U>();
        return out;
    }

    /**
     * @brief Check if the closed loop is stable: every pole strictly inside the unit circle
     *        (Ts > 0), or strictly in the left half plane (Ts = 0). False for a failed design.
     */
    [[nodiscard]] constexpr bool is_stable() const {
        if (!success) {
            return false; // a failed design has no closed loop (e is zero)
        }
        for (size_t i = 0; i < NK; ++i) {
            if (Ts > T{0} ? (e[i].abs() >= T{1}) : (e[i].real() >= T{0})) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief The law as a discrete state-space block, for closing the loop in analysis
     *
     * The realization of the runtime's tick, by the law's features:
     *
     * | Law | State | Inputs | Output |
     * | --- | --- | --- | --- |
     * | measured state | none | x | u = −Kx |
     * | measured state, integrators | ξ | [r; x] | u |
     * | Kalman filter | x̂(k−1\|k−1) | y | u |
     * | Kalman filter, integrators | [x̂(k\|k−1); ξ; u_prev] | [r; y] | u |
     *
     * With input registers, x is the law's model [x; v(k−1); …], the state of the delayed plant the
     * loop closes on.
     *
     * With integrators, @f$ u = -K_x x - K_i \xi @f$ reads the current ξ, then
     * @f$ \xi \leftarrow \xi + h(r - C_i (C x + D u)) @f$ (h = integral.step), the measurement
     * less the feedthrough of the input it was taken under standing in for C x with a filter. With a filter it is the
     * current estimator: predict with the last command, the measurement update with the
     * innovation gain M (KalmanResult::M), then the law,
     * @f[
     *   A_c = (I-MC)(A-BK) + MDK,\quad B_c = M,\quad C_c = -K A_c,\quad D_c = -KM.
     * @f]
     * A law with a disturbance model holds more than these
     * states, so it has no such block. The steady-state target of TrackedOutputs is a
     * feedforward on r outside this block: the block is the regulator (r = 0). The block is the
     * law within its limits: u_lim, the integrators' limits and demand_lim are not part of it.
     */
    [[nodiscard]] constexpr auto to_ss() const
        requires(ND == 0)
    {
        if constexpr (!kEstimated && NI == 0) {
            return StateSpace<0, NXL, NU, T>{.A = Matrix<0, 0, T>{}, .B = Matrix<0, NXL, T>{}, .C = Matrix<NU, 0, T>{}, .D = -K, .Ts = Ts};
        } else if constexpr (!kEstimated) {
            return integral_ss();
        } else if constexpr (NI == 0) {
            return estimator_ss();
        } else {
            return integral_estimator_ss();
        }
    }

private:
    /**
     * @brief The law with integrators on a measured state: u = C_c ξ + D_c [r; x], C_c = −K_i,
     *        D_c = [0, −K_x], and ξ ← ξ + h·(r − Ci·C·x − Ci·D·u): A_c = I − h·Ci·D·C_c,
     *        B_c = [hI, −h·Ci·C] − h·Ci·D·D_c.
     */
    [[nodiscard]] constexpr StateSpace<NI, NI + NXL, NU, T> integral_ss() const {
        const T                 h = integral.step;
        const Matrix<NU, NI, T> Cc = -K.template block<NU, NI>(0, NXL);
        Matrix<NU, NI + NXL, T> Dc{};
        Dc.template block<NU, NXL>(0, NI) = -K.template block<NU, NXL>(0, 0);
        Matrix<NI, NI + NXL, T> Bc{};
        Bc.template block<NI, NI>(0, 0) = Matrix<NI, NI, T>::identity() * h;
        Bc.template block<NI, NXL>(0, NI) = -integral.Ci_x * h;
        const Matrix<NI, NU, T> hCiD = integral.Ci_u * h;
        const Matrix<NI, NI, T> Ac = Matrix<NI, NI, T>::identity() - (hCiD * Cc);
        return {.A = Ac, .B = Matrix<NI, NI + NXL, T>(Bc - (hCiD * Dc)), .C = Cc, .D = Dc, .Ts = Ts};
    }

    /** @brief The current-estimator law: the filter's estimate x̂(k−1|k−1) is the state. */
    [[nodiscard]] constexpr StateSpace<NXL, NY, NU, T> estimator_ss() const {
        const auto& sys = kalman.sys;
        const auto  ImMC = Matrix<NXL, NXL, T>::identity() - (kalman.M * sys.C);
        const auto  Ac = (ImMC * (sys.A - (sys.B * K))) + ((kalman.M * sys.D) * K);
        return {.A = Ac, .B = kalman.M, .C = -(K * Ac), .D = -(K * kalman.M), .Ts = Ts};
    }

    /** @brief The current-estimator law with integrators: the state is [x̂(k|k−1); ξ; u_prev], the inputs [r; y] (to_ss()). */
    [[nodiscard]] constexpr StateSpace<NXL + NI + NU, NI + NY, NU, T> integral_estimator_ss() const {
        constexpr size_t         NS = NXL + NI + NU;
        const auto&              sys = kalman.sys;
        const auto&              M = kalman.M;
        const Matrix<NU, NXL, T> Kx = K.template block<NU, NXL>(0, 0);
        const Matrix<NU, NI, T>  Ki = K.template block<NU, NI>(0, NXL);
        const auto               ImMC = Matrix<NXL, NXL, T>::identity() - (M * sys.C);
        const auto               MD = M * sys.D;
        const auto               KxImMC = Kx * ImMC;
        const auto               AmBKx = sys.A - (sys.B * Kx);
        const auto               KxM = Kx * M;
        const T                  h = integral.step;

        Matrix<NS, NS, T> Ac{};
        Ac.template block<NXL, NXL>(0, 0) = AmBKx * ImMC;
        Ac.template block<NXL, NI>(0, NXL) = -(sys.B * Ki);
        Ac.template block<NXL, NU>(0, NXL + NI) = -(AmBKx * MD);
        Ac.template block<NI, NI>(NXL, NXL) = Matrix<NI, NI, T>::identity();
        Ac.template block<NU, NXL>(NXL + NI, 0) = -KxImMC;
        Ac.template block<NU, NI>(NXL + NI, NXL) = -Ki;
        Ac.template block<NU, NU>(NXL + NI, NXL + NI) = Kx * MD;
        Matrix<NS, NI + NY, T> Bc{};
        Bc.template block<NI, NI>(NXL, 0) = Matrix<NI, NI, T>::identity() * h;
        Bc.template block<NXL, NY>(0, NI) = AmBKx * M;
        Bc.template block<NI, NY>(NXL, NI) = -integral.Ci * h;
        Bc.template block<NU, NY>(NXL + NI, NI) = -KxM;
        Matrix<NU, NS, T> Cc{};
        Cc.template block<NU, NXL>(0, 0) = -KxImMC;
        Cc.template block<NU, NI>(0, NXL) = -Ki;
        Cc.template block<NU, NU>(0, NXL + NI) = Kx * MD;
        Matrix<NU, NI + NY, T> Dc{};
        Dc.template block<NU, NY>(0, NI) = -KxM;
        //! ξ ← ξ + h·(r − Ci·y + Ci·D·u_prev − Ci·D·u): the measurement less the feedthrough of the
        //! input it was taken under, then the feedthrough of the input the law applies, u = C_c s + D_c [r; y].
        const Matrix<NI, NU, T> hCiD = integral.Ci_u * h;
        Ac.template block<NI, NU>(NXL, NXL + NI) = hCiD;
        const Matrix<NI, NS, T>      fed_a = hCiD * Cc;
        const Matrix<NI, NI + NY, T> fed_b = hCiD * Dc;
        for (size_t i = 0; i < NI; ++i) {
            for (size_t j = 0; j < NS; ++j) {
                Ac(NXL + i, j) -= fed_a(i, j);
            }
            for (size_t j = 0; j < NI + NY; ++j) {
                Bc(NXL + i, j) -= fed_b(i, j);
            }
        }
        return {.A = Ac, .B = Bc, .C = Cc, .D = Dc, .Ts = Ts};
    }
};

} // namespace design

namespace detail {

/** @brief Not a StateFeedbackResult. */
template<typename R>
struct is_state_feedback_result : damp::false_type {};

/** @brief A StateFeedbackResult. */
template<size_t NX, size_t NU, typename T, size_t NI, size_t NY, typename Kalman, size_t NZ, size_t NRU>
struct is_state_feedback_result<design::StateFeedbackResult<NX, NU, T, NI, NY, Kalman, NZ, NRU>> : damp::true_type {};

/**
 * @brief A StateFeedbackResult of any sizes. A runtime's constructor takes one and reads its kind
 *        from the result's own constants (kIntegral, kEstimated, ND), so a parameter added to the
 *        result does not reach the runtimes' constructors.
 */
template<typename R>
concept state_feedback_result = is_state_feedback_result<R>::value;

} // namespace detail

} // namespace damp
