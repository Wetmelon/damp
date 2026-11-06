// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file matlab.hpp
 * @brief MATLAB®-style free-function aliases for Damp APIs (host-oriented)
 *
 * Thin spellings over Damp design/analysis/systems for Control System Toolbox
 * muscle memory. They are MATLAB-inspired names, not drop-in MATLAB: argument
 * order, defaults, and return types follow the Damp primary.
 *
 * | MATLAB® | Damp |
 * | ------- | ---- |
 * | `ss`, `tf`, `zpk` | `matlab::ss` / `tf` / `zpk` → `StateSpace`, `TransferFunction`, `ZPK` |
 * | `c2d`, `d2d` | `matlab::c2d` / `d2d` → `discretize` |
 * | `lqr`, `dlqr`, `lqrd`, `lqi`, `lqg` | `matlab::lqr` … → `design::lqr` / `lqi` / `lqg` |
 * | `place`, `acker` | `matlab::place` / `acker` → `design::place` |
 * | `pidtune` | `matlab::pidtune` → `design::pidtune` |
 * | `kalmd`, `dlqe` | `matlab::kalmd` / `dlqe` → `design::kalman` |
 * | `bode`, `nyquist`, `margin` | `matlab::bode` … → `analysis::` |
 *
 * Include via `damp/workbench.hpp`. Each alias returns what its primary returns.
 * Constructors (`ss`, `tf`, `zpk`) return system types. The design and conversion aliases
 * (@c c2d, @c place, @c acker, @c pidtune, @c lqr / @c dlqr / @c lqrd / @c lqi / @c lqg /
 * @c kalmd, @c canon, @c minreal, @c balreal, @c balred) return the primary's Result,
 * whose @c success and @c status say whether and why it failed; @c gram, @c lyap,
 * @c dlyap, @c lyapchol, @c norm, @c expm / @c sqrtm / @c logm, @c eig, @c pole,
 * @c hankelsv, @c modred and @c thiran return @c damp::optional. Frequency aliases
 * (@c bode, @c nyquist, @c sigma, @c rlocus, …) skip contour poles honestly (shorter
 * grids) via @c analysis::.
 */

#include <cstddef>
#include <limits>
#include <string_view>
#include <type_traits>
#include <vector>

#include "damp/analysis/analysis.hpp"
#include "damp/backend.hpp"
#include "damp/controllers/lqg.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/design/linearization.hpp"
#include "damp/design/minreal.hpp"
#include "damp/design/model_reduction.hpp"
#include "damp/design/model_result.hpp"
#include "damp/design/pid_design.hpp"
#include "damp/design/pole_placement.hpp"
#include "damp/design/qp.hpp"
#include "damp/design/riccati.hpp"
#include "damp/design/status.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/filters/biquad.hpp"
#include "damp/filters/pade.hpp"
#include "damp/math/complex.hpp"
#include "damp/matrix/eigen.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/svd.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/realization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"
#include "damp/systems/zeros.hpp"
#include "damp/systems/zpk.hpp"
#include "design/stability.hpp"

namespace damp {
// MATLAB®-style matrix functions
namespace matlab {

/**
 * @brief MATLAB®-style transfer function constructor
 *
 * Coefficients are highest power of s or z first, as MATLAB®'s: MATLAB® tf([1 2], [1 3 2]) is
 * matlab::tf({1.0, 2.0}, {1.0, 3.0, 2.0}), (s + 2)/(s² + 3s + 2).
 */
template<size_t Nnum, size_t Nden, typename T = double>
[[nodiscard]] constexpr TransferFunction<Nnum, Nden, T>
tf(const damp::array<T, Nnum>& num, const damp::array<T, Nden>& den) noexcept {
    return TransferFunction<Nnum, Nden, T>{.num = num, .den = den};
}

/**
 * @brief MATLAB®-style transfer function constructor from braced coefficient lists
 *
 * Enables calls like: tf({1.0, 2.0}, {1.0, 3.0, 2.0})
 * Coefficients are highest power of s or z first, as MATLAB®'s.
 */
template<typename TNum, size_t Nnum, typename TDen, size_t Nden>
[[nodiscard]] constexpr TransferFunction<Nnum, Nden, std::common_type_t<TNum, TDen>>
tf(const TNum (&num)[Nnum], const TDen (&den)[Nden]) noexcept {
    using T = std::common_type_t<TNum, TDen>;
    TransferFunction<Nnum, Nden, T> result{};

    for (size_t i = 0; i < Nnum; ++i) {
        result.num[i] = static_cast<T>(num[i]);
    }
    for (size_t i = 0; i < Nden; ++i) {
        result.den[i] = static_cast<T>(den[i]);
    }

    return result;
}

/**
 * @brief MATLAB®-style state-space model constructor.
 *
 * Builds a StateSpace from (A, B, C, D) with an optional sampling period
 * (Ts = 0 for continuous, > 0 for discrete). Process/measurement noise
 * channels are omitted (NW = NV = 0); construct StateSpace directly for those.
 *
 * @note Compare with MATLAB®'s sys = ss(A, B, C, D) / ss(A, B, C, D, Ts).
 */
template<size_t NX, size_t NU, size_t NY, typename T = double>
[[nodiscard]] constexpr StateSpace<NX, NU, NY, T> ss(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NY, NX, T>& C,
    const Matrix<NY, NU, T>& D = Matrix<NY, NU, T>::zeros(),
    T                        Ts = T{0}
) noexcept {
    return StateSpace<NX, NU, NY, T>{
        .A = A,
        .B = B,
        .C = C,
        .D = D,
        .Ts = Ts,
    };
}

/**
 * @brief MATLAB®-style zero-pole-gain model constructor.
 *
 * @param z Zeros, @param p Poles, @param k Scalar gain.
 * @note Compare with MATLAB®'s sys = zpk(z, p, k).
 */
template<size_t Nz, size_t Np, typename T = double>
[[nodiscard]] constexpr ZPK<Nz, Np, T> zpk(
    const damp::array<damp::complex<T>, Nz>& z,
    const damp::array<damp::complex<T>, Np>& p,
    T                                        k
) noexcept {
    return ZPK<Nz, Np, T>{z, p, k};
}

/**
 * @brief MATLAB® short alias for to_zpk (transfer function → zero-pole-gain)
 * @note Compare with MATLAB®'s [z, p, k] = tf2zpk(b, a) / zpk(tf).
 */
template<size_t Nnum, size_t Nden, typename T = double>
[[nodiscard]] constexpr auto tf2zpk(const TransferFunction<Nnum, Nden, T>& tf) {
    return to_zpk(tf);
}

/**
 * @brief MATLAB® short alias for ZPK::to_transfer_function
 * @note Compare with MATLAB®'s [b, a] = zpk2tf(z, p, k) / tf(zpk).
 */
template<size_t Nz, size_t Np, typename T = double>
[[nodiscard]] constexpr TransferFunction<Nz + 1, Np + 1, T> zpk2tf(const ZPK<Nz, Np, T>& sys) {
    return sys.to_transfer_function();
}

/**
 * @brief MATLAB®-style transfer function → second-order sections, from MATLAB®'s z⁻¹ vectors
 *
 * @p b and @p a are as MATLAB® takes them, ascending powers of z⁻¹ (b(1) + b(2)z⁻¹ + …), the
 * shorter one padded with trailing zeros: MATLAB® tf2sos([1 3 3 1]/6, [3 0 1 0]/3) is
 * matlab::tf2sos({1.0 / 6, 0.5, 0.5, 1.0 / 6}, {1.0, 0.0, 1.0 / 3, 0.0}). Run the result with
 * BiquadCascade.
 *
 * @note Compare with MATLAB®'s sos = tf2sos(b, a): order 'up', the gain in the first section.
 * @see design::second_order_sections
 */
template<size_t NB, size_t NA, typename T = double>
[[nodiscard]] constexpr auto tf2sos(const damp::array<T, NB>& b, const damp::array<T, NA>& a) {
    constexpr size_t          N = (NB > NA) ? NB : NA;
    TransferFunction<N, N, T> sys{};
    for (size_t i = 0; i < NB; ++i) {
        sys.num[i] = b[i];
    }
    for (size_t i = 0; i < NA; ++i) {
        sys.den[i] = a[i];
    }
    return design::second_order_sections(sys);
}

/**
 * @brief MATLAB®-style tf2sos from braced coefficient lists, ascending powers of z⁻¹.
 * @note Compare with MATLAB®'s sos = tf2sos(b, a).
 */
template<typename TB, size_t NB, typename TA, size_t NA>
[[nodiscard]] constexpr auto tf2sos(const TB (&b)[NB], const TA (&a)[NA]) {
    using T = std::common_type_t<TB, TA>;
    damp::array<T, NB> bb{};
    damp::array<T, NA> aa{};
    for (size_t i = 0; i < NB; ++i) {
        bb[i] = static_cast<T>(b[i]);
    }
    for (size_t i = 0; i < NA; ++i) {
        aa[i] = static_cast<T>(a[i]);
    }
    return tf2sos(bb, aa);
}

/**
 * @brief MATLAB®-style zero-pole-gain → second-order sections, H(z) = k·Π(z − zᵢ)/Π(z − pⱼ).
 * @note Compare with MATLAB®'s sos = zp2sos(z, p, k).
 * @see design::second_order_sections
 */
template<size_t Nz, size_t Np, typename T = double>
[[nodiscard]] constexpr auto zp2sos(const damp::array<damp::complex<T>, Nz>& z, const damp::array<damp::complex<T>, Np>& p, T k) {
    return design::second_order_sections(ZPK<Nz, Np, T>{z, p, k});
}

/**
 * @brief MATLAB® short alias for SISO to_zpk (state space → zero-pole-gain)
 * @note Compare with MATLAB®'s [z, p, k] = ss2zpk(A, B, C, D) / zpk(ss).
 */
template<size_t NX, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr StateSpaceZPKResult<NX, T> ss2zpk(const StateSpace<NX, 1, 1, T, NW, NV>& sys) {
    return to_zpk(sys);
}

/**
 * @brief MATLAB®-style SISO SS → TF.
 * @note Compare with MATLAB®'s ss2tf.
 */
using damp::ss2tf;

/**
 * @brief MATLAB®-style parallel-form continuous PID constructor.
 *
 * Returns the continuous design::PIDGains, as the rest of the PID family (pidstd,
 * pidstd2, make1DOF, make2DOF, pidtune) does; @p Tf is the derivative filter time
 * constant (0 = unfiltered). Deploy fixed-rate with `PIDController{matlab::pid(...).discretize(Ts)}`,
 * or hand the gains to `ContinuousPID{matlab::pid(...)}` when the period is measured each tick.
 *
 * @note Compare with MATLAB®'s C = pid(Kp, Ki, Kd, Tf). MATLAB®'s filter
 *       coefficient N maps to Tf = 1/N.
 */
template<typename T = double>
[[nodiscard]] constexpr design::PIDGains<T> pid(T Kp, T Ki = T{0}, T Kd = T{0}, T Tf = T{0}) noexcept {
    return design::PIDGains<T>{.Kp = Kp, .Ki = Ki, .Kd = Kd, .Tf = Tf};
}

/**
 * @brief Standard-form continuous PID constructor (1-DOF).
 *
 * Maps industrial standard form
 * @f[
 *   C(s) = K_p\Bigl(1 + \frac{1}{T_i s} + \frac{T_d s}{1+(T_d/N)s}\Bigr)
 * @f]
 * to the shipped parallel @ref design::PIDGains with
 * @f$K_i = K_p/T_i@f$, @f$K_d = K_p T_d@f$, @f$T_f = T_d/N@f$ when @p N > 0
 * (else unfiltered). @p N = 0 (the default) means no derivative filter, matching
 * MATLAB®'s default N = ∞. @p Ti <= 0 yields @f$K_i = 0@f$ (P/PD).
 *
 * Deploy with `.discretize(Ts)` → PIDController, same as other PID designs.
 *
 * @note Compare with MATLAB®'s C = pidstd(Kp, Ti, Td, N) (default N = ∞, spelled 0 here).
 * @see design::pid for the parallel-form constructor
 * @see pidstd2 for the 2-DOF standard-form variant
 */
template<typename T = double>
[[nodiscard]] constexpr design::PIDGains<T> pidstd(T Kp, T Ti, T Td = T{0}, T N = T{0}) noexcept {
    const T Ki = (Ti > T{0}) ? (Kp / Ti) : T{0};
    const T Kd = Kp * Td;
    const T Tf = (N > T{0} && Td > T{0}) ? (Td / N) : T{0};
    return design::PIDGains<T>{.Kp = Kp, .Ki = Ki, .Kd = Kd, .Tf = Tf};
}

/**
 * @brief Standard-form continuous 2-DOF PID constructor.
 *
 * Same Ti/Td → Ki/Kd map as @ref pidstd, with proportional/derivative setpoint
 * weights @p b and @p c (see @ref design::PIDGains). Both default to 1, as MATLAB®'s;
 * c = 0 puts the derivative on the measurement (PI-D).
 *
 * @note Compare with MATLAB®'s C = pidstd2(Kp, Ti, Td, N, b, c).
 * @see pidstd, make1DOF, make2DOF
 */
template<typename T = double>
[[nodiscard]] constexpr design::PIDGains<T>
pidstd2(T Kp, T Ti, T Td = T{0}, T N = T{0}, T b = T{1}, T c = T{1}) noexcept {
    auto r = pidstd(Kp, Ti, Td, N);
    r.b = b;
    r.c = c;
    return r;
}

/**
 * @brief Force 1-DOF setpoint weights on a PID design result (b=c=1).
 *
 * @note Compare with MATLAB®'s make1DOF (Control System Toolbox).
 */
template<typename T = double>
[[nodiscard]] constexpr design::PIDGains<T> make1DOF(design::PIDGains<T> r) noexcept {
    r.b = T{1};
    r.c = T{1};
    return r;
}

/**
 * @brief Apply 2-DOF setpoint weights on a PID design result.
 *
 * @p b and @p c default to 1, as MATLAB®'s, so make2DOF(r) keeps the response of r;
 * c = 0 puts the derivative on the measurement (PI-D).
 *
 * @note Compare with MATLAB®'s make2DOF(C1, b, c) (Control System Toolbox).
 */
template<typename T = double>
[[nodiscard]] constexpr design::PIDGains<T>
make2DOF(design::PIDGains<T> r, T b = T{1}, T c = T{1}) noexcept {
    r.b = b;
    r.c = c;
    return r;
}

/**
 * @brief MATLAB® short alias for design::minreal, with its default tolerance √ε
 *
 * @note Compare with MATLAB®'s minreal(sys, tol).
 * @see design::minreal()
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr design::MinimalRealizationResult<NX, NU, NY, T, NW, NV>
minreal(const StateSpace<NX, NU, NY, T, NW, NV>& sys, T tol = damp::sqrt(eps<T>())) {
    return design::minreal(sys, tol);
}

/**
 * @brief MATLAB® short alias for design::balanced_realization
 * @note Compare with MATLAB®'s [sysb, g] = balreal(sys).
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr design::BalancedRealizationResult<NX, NU, NY, T, NW, NV> balreal(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    return design::balanced_realization(sys);
}

/**
 * @brief MATLAB® short alias for design::hankel_singular_values
 * @note Compare with MATLAB®'s hsv = hankelsv(sys).
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr damp::optional<ColVec<NX, T>> hankelsv(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    return design::hankel_singular_values(sys);
}

/**
 * @brief MATLAB® short alias for design::balanced_reduction to NR states
 * @note Compare with MATLAB®'s balred(sys, NR, balredOptions('StateElimMethod', …)); the
 *       default here is Truncate, MATLAB®'s is MatchDC.
 */
template<size_t NR, size_t NX, size_t NU, size_t NY, size_t NW = 0, size_t NV = 0, typename T = double>
[[nodiscard]] constexpr design::BalancedReductionResult<NR, NU, NY, T, NW, NV> balred(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    design::ModelReductionMethod             method = design::ModelReductionMethod::Truncate
) {
    return design::balanced_reduction<NR>(sys, method);
}

/**
 * @brief MATLAB® short alias for design::eliminate_states (0-based indices)
 * @note Compare with MATLAB®'s modred(sys, elim, 'MatchDC'|'Truncate'); MATLAB® indices are 1-based.
 */
template<size_t... Elim, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
    requires(design::detail::valid_elimination<NX, Elim...>())
[[nodiscard]] constexpr damp::optional<StateSpace<NX - sizeof...(Elim), NU, NY, T, NW, NV>> modred(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    design::ModelReductionMethod             method = design::ModelReductionMethod::MatchDC
) {
    return design::eliminate_states<Elim...>(sys, method);
}

/**
 * @brief MATLAB® short alias for controllability_matrix
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr Matrix<NX, NX * NU, T> ctrb(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B) noexcept {
    return stability::controllability_matrix(A, B);
}

/**
 * @brief MATLAB® short alias for controllability_matrix taking StateSpace
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr Matrix<NX, NX * NU, T> ctrb(const StateSpace<NX, NU, NY, T, NW, NV>& sys) noexcept {
    return stability::controllability_matrix(sys.A, sys.B);
}

/**
 * @brief MATLAB® short alias for observability_matrix
 */
template<size_t NX, size_t NY, typename T = double>
[[nodiscard]] constexpr Matrix<NX * NY, NX, T> obsv(const Matrix<NX, NX, T>& A, const Matrix<NY, NX, T>& C) noexcept {
    return stability::observability_matrix(A, C);
}

/**
 * @brief MATLAB® short alias for observability_matrix
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr Matrix<NX * NY, NX, T> obsv(const StateSpace<NX, NU, NY, T, NW, NV>& sys) noexcept {
    return stability::observability_matrix(sys.A, sys.C);
}

/**
 * @brief MATLAB® interface function c2d to discretize a continuous-time state-space system
 *
 * @param sys               Continuous-time state-space model (Ts should be 0)
 * @param sampling_time     Desired sampling period for discrete system
 * @param method            Discretization method (ZOH or Tustin)
 * @param prewarp_frequency [rad/s] Tustin pre-warp frequency (c2dOptions 'PrewarpFrequency'); 0 for none
 *
 * @return Discrete StateSpace and status, as damp::discretize (e.g. DesignStatus::singular for a singular Tustin factor)
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr DiscretizationResult<NX, NU, NY, T, NW, NV> c2d(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        sampling_time,
    DiscretizationMethod                     method = DiscretizationMethod::ZOH,
    T                                        prewarp_frequency = T{0}
) {
    return discretize(sys, sampling_time, method, prewarp_frequency);
}

/**
 * @brief MATLAB®-style c2d for SISO transfer functions
 *
 * Realizes the TF then discretizes. Refused with DesignStatus::singular if the realization
 * fails (den[0] = 0), otherwise with damp::discretize's status.
 */
template<size_t Nnum, size_t Nden, typename T = double>
[[nodiscard]] constexpr DiscretizationResult<Nden - 1, 1, 1, T> c2d(
    const TransferFunction<Nnum, Nden, T>& tf_sys,
    T                                      sampling_time,
    DiscretizationMethod                   method = DiscretizationMethod::ZOH
) {
    const auto ss = tf_sys.to_state_space();
    if (!ss) {
        return ::damp::detail::refused<DiscretizationResult<Nden - 1, 1, 1, T>>(design::DesignStatus::singular);
    }
    return discretize(*ss, sampling_time, method);
}

/**
 * @brief MATLAB® d2d: resample a ZOH-discrete model to a new sample time
 * @note Forwards to damp::resample. Compare with MATLAB®'s d2d(sys, Ts) (ZOH method).
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr DiscretizationResult<NX, NU, NY, T, NW, NV> d2d(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    std::type_identity_t<T>                  sampling_time
) {
    return resample(sys, sampling_time);
}

/**
 * @brief MATLAB® upsample: sys(z) → sys(zᴸ) at Ts / L
 * @note Forwards to damp::upsample. Compare with MATLAB®'s upsample(sys, L); L is a template argument here.
 */
template<size_t L, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
    requires(L >= 1)
[[nodiscard]] constexpr StateSpace<NX * L, NU, NY, T, NW, NV> upsample(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    return damp::upsample<L>(sys);
}

/**
 * @brief MATLAB® prescale: diagonal state scaling that balances A, B and C
 * @note Forwards to damp::balance_states, the ssbal-style balance. MATLAB®'s prescale(sys)
 *       chooses its factors over a frequency band; both return a scaled model and the scaling.
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr StateScaling<NX, NU, NY, T, NW, NV> prescale(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    return balance_states(sys);
}

/**
 * @brief MATLAB® ssbal: diagonal state balancing of [A B; C 0]
 * @note Forwards to damp::balance_states. Compare with MATLAB®'s [sysb, T] = ssbal(sys): T = diag(1 / scale()).
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr StateScaling<NX, NU, NY, T, NW, NV> ssbal(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    return balance_states(sys);
}

/**
 * @brief MATLAB® canon(sys) / canon(sys, 'modal'): real modal canonical form and its transform
 * @note Forwards to damp::modal_form; MATLAB®'s default form is 'modal'.
 *       Compare with MATLAB®'s [csys, T] = canon(sys, 'modal'): csys = result.sys, T = result.transform.
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr StateTransformResult<NX, NU, NY, T, NW, NV> canon(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    return modal_form(sys);
}

/**
 * @brief MATLAB® canon(sys, type): 'modal' or 'companion' canonical form and its transform
 * @note Forwards to damp::modal_form or damp::companion_form (SISO only); any other type, or
 *       'companion' on a MIMO model, is refused with DesignStatus::bad_weights.
 *       Compare with MATLAB®'s [csys, T] = canon(sys, type).
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr StateTransformResult<NX, NU, NY, T, NW, NV>
canon(const StateSpace<NX, NU, NY, T, NW, NV>& sys, std::string_view type) {
    if (type == "modal") {
        return modal_form(sys);
    }
    if constexpr (NU == 1 && NY == 1) {
        if (type == "companion") {
            return companion_form(sys);
        }
    }
    StateTransformResult<NX, NU, NY, T, NW, NV> refused{};
    return damp::detail::set_status(refused, design::DesignStatus::bad_weights);
}

/**
 * @brief MATLAB® ss2ss(sys, T): the model in the state z = T·x
 * @note Forwards to damp::similarity_transform. Compare with MATLAB®'s ss2ss(sys, T).
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr StateTransformResult<NX, NU, NY, T, NW, NV>
ss2ss(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const Matrix<NX, NX, T>& transform) {
    return similarity_transform(sys, transform);
}

/**
 * @brief MATLAB® ctrbf(A, B, C): controllability staircase form, its orthogonal transform and rank
 * @note Forwards to damp::controllability_staircase. Compare with MATLAB®'s
 *       [Abar, Bbar, Cbar, T, k] = ctrbf(A, B, C): result.sys, result.transform, result.rank = sum(k).
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr StaircaseFormResult<NX, NU, NY, T, NW, NV> ctrbf(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    return controllability_staircase(sys);
}

/**
 * @brief MATLAB® obsvf(A, B, C): observability staircase form, its orthogonal transform and rank
 * @note Forwards to damp::observability_staircase. Compare with MATLAB®'s
 *       [Abar, Bbar, Cbar, T, k] = obsvf(A, B, C): result.sys, result.transform, result.rank = sum(k).
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr StaircaseFormResult<NX, NU, NY, T, NW, NV> obsvf(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    return observability_staircase(sys);
}

/**
 * @brief MATLAB® tzero(sys): the finite invariant zeros of a model with any number of inputs and outputs
 * @note Forwards to damp::transmission_zeros. Compare with MATLAB®'s z = tzero(sys):
 *       result.zeros[0 .. result.n_zeros).
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr TransmissionZerosResult<NX, T> tzero(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    return transmission_zeros(sys);
}

/**
 * @brief MATLAB® thiran: fractional-delay all-pass of order N for a delay tau at Ts
 * @note Forwards to damp::fractional_delay. Compare with MATLAB®'s thiran(tau, Ts), which
 *       chooses N = ⌈tau/Ts⌉; here N is a template argument.
 */
template<size_t N, typename T>
    requires(N >= 1)
[[nodiscard]] constexpr damp::optional<TransferFunction<N + 1, N + 1, T>> thiran(T tau, std::type_identity_t<T> Ts) {
    return fractional_delay<N>(tau, Ts);
}

/**
 * @brief MATLAB®-style nonlinear linearization about an operating point
 *
 * Forwards to design::linearize: the Jacobians of the Dynamics ẋ = f(x, u) and the
 * Measurement y = h(x, u) about (x_op, u_op), analytic when the model provides them,
 * else by central differences.
 *
 * @return design::LinearizationResult with Jacobians A, B, C, D and the operating point
 * @note Compare with MATLAB®'s [A,B,C,D] = linmod(model, x, u). MATLAB®'s several
 *       outputs are one design::LinearizationResult here.
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, typename F, typename H>
[[nodiscard]] constexpr design::LinearizationResult<NX, NU, NY, T> linmod(
    const F&             dynamics,
    const H&             output,
    const ColVec<NX, T>& x_op,
    const ColVec<NU, T>& u_op,
    T                    epsilon = T{0}
) {
    return design::linearize<NX, NU, NY, T>(dynamics, output, x_op, u_op, epsilon);
}

/** @brief MATLAB® spellings of core size-deducing vocabulary (see matrix/concat.hpp). */
using damp::blkdiag;
using damp::diag;

/**
 * @brief Create an identity matrix of size n x n
 *
 * @tparam N Size of the identity matrix
 * @return constexpr Matrix<N, N, T>
 */
template<size_t N, typename T = double>
[[nodiscard]] constexpr auto eye() noexcept {
    return Matrix<N, N, T>::identity();
}

/**
 * @brief N×N matrix of ones
 * @note Compare with MATLAB®'s ones(N).
 */
template<size_t N, typename T = double>
[[nodiscard]] constexpr auto ones() noexcept {
    return Matrix<N, N, T>::constant(T{1});
}

/**
 * @brief MATLAB® short alias for the singular value decomposition.
 *
 * Returns a mat::SVDResult{singular_U, singular_values, singular_V} rather than
 * the MATLAB® [U, S, V] tuple; the singular values are a descending array, not a
 * diagonal matrix.
 *
 * @note Compare with MATLAB®'s [U, S, V] = svd(A).
 */
template<size_t M, size_t N, typename T>
[[nodiscard]] constexpr auto svd(const Matrix<M, N, T>& A) {
    return mat::svd(A);
}

/**
 * @brief MATLAB® short alias for the Moore–Penrose pseudoinverse.
 * @note Compare with MATLAB®'s pinv(A).
 */
template<size_t M, size_t N, typename T>
[[nodiscard]] constexpr Matrix<N, M, T> pinv(const Matrix<M, N, T>& A) {
    return mat::pseudo_inverse(A);
}

/**
 * @brief MATLAB® short alias for an orthonormal null-space basis.
 *
 * Returns a mat::NullSpace whose trailing `dim` columns are the kernel basis
 * (MATLAB®'s null returns those columns directly as an N×dim matrix).
 *
 * @note Compare with MATLAB®'s null(A).
 */
template<size_t M, size_t N, typename T>
[[nodiscard]] constexpr mat::NullSpace<N, T> null(const Matrix<M, N, T>& A) {
    return mat::null_space(A);
}

/**
 * @brief MATLAB® short alias for the eigenvalues of a square matrix.
 *
 * Returns the eigenvalues, from balancing, Hessenberg reduction and Francis double-shift QR,
 * or nullopt when the iteration does not converge or A has a non-finite entry. For the real
 * Schur basis, call mat::compute_eigenvalues() directly; neither returns eigenvectors of a
 * general A.
 *
 * @note Compare with MATLAB®'s e = eig(A).
 */
template<size_t N, typename T>
[[nodiscard]] constexpr damp::optional<damp::array<damp::complex<T>, N>> eig(const Matrix<N, N, T>& A) {
    return analysis::poles(A);
}

/**
 * @brief MATLAB® short alias for the matrix exponential (mat::expm).
 *
 * Returns nullopt for a non-finite input or a failed Padé solve.
 *
 * @note Compare with MATLAB®'s expm(A).
 */
template<size_t N, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<N, N, T>> expm(const Matrix<N, N, T>& A) {
    return mat::expm(A);
}

/**
 * @brief MATLAB® short alias for the principal matrix square root (mat::sqrt).
 *
 * Returns nullopt when mat::sqrt does: a non-finite entry, no real principal root, or a result
 * that fails its residual check.
 *
 * @note Compare with MATLAB®'s sqrtm(A).
 */
template<size_t N, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<N, N, T>> sqrtm(const Matrix<N, N, T>& A) {
    return mat::sqrt(A);
}

/**
 * @brief MATLAB® short alias for the principal matrix logarithm (mat::log).
 *
 * Returns nullopt when mat::log does: an eigenvalue on the closed negative real axis, or a square
 * root step that fails.
 *
 * @note Compare with MATLAB®'s logm(A).
 */
template<size_t N, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<N, N, T>> logm(const Matrix<N, N, T>& A) {
    return mat::log(A);
}

/**
 * @brief Form the state estimator from a plant and an estimator gain
 *
 * est = estim(sys, L) is the full-order estimator of gain L driven by the measured outputs y,
 * with the state estimate as its output:
 *
 * @f[
 *   \dot{\hat x} = (A - LC)\,\hat x + L y, \qquad \hat y = \hat x .
 * @f]
 *
 * The plant inputs u are not estimator inputs (as in matlab::reg), and the noise channels
 * G, H do not carry over. @p sys.Ts carries over, so a discrete plant gives a discrete estimator.
 *
 * @note Compare with MATLAB®'s est = estim(sys, L) with every output measured and no known input.
 *
 * @param sys State-space plant (discrete or continuous)
 * @param L   Estimator gain
 * @return Estimator StateSpace<NX, NY, NX> (input y, output x̂), with sys.Ts
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr StateSpace<NX, NY, NX, T> estim(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const Matrix<NX, NY, T>& L) noexcept {
    return StateSpace<NX, NY, NX, T>{
        .A = sys.A - (L * sys.C),
        .B = L,
        .C = Matrix<NX, NX, T>::identity(),
        .Ts = sys.Ts,
    };
}

/**
 * @brief Form the observer-based compensator y → u from a plant, K and L
 *
 * rsys = reg(sys, K, L) connects the state-feedback law u = −K x̂ with the
 * full-order estimator of gain L (see estim). All plant inputs are controls and
 * all outputs are measured, so the regulator's input is y and its output is u:
 *
 * @f[
 *   \dot{\hat x} = (A - BK - LC + LDK)\,\hat x + L y, \qquad u = -K \hat x .
 * @f]
 *
 * Closing the loop with the plant places the poles at eig(A − BK) ∪ eig(A − LC)
 * (separation principle), in continuous or discrete time alike.
 *
 * @note Compare with MATLAB®'s rsys = reg(sys, K, L).
 *
 * @param sys State-space plant (discrete or continuous)
 * @param K   State-feedback gain (u = −K x̂)
 * @param L   Estimator gain
 * @return Compensator StateSpace<NX, NY, NU> (input y, output u), with sys.Ts
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr StateSpace<NX, NY, NU, T> reg(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NU, NX, T>&                 K,
    const Matrix<NX, NY, T>&                 L
) noexcept {
    return StateSpace<NX, NY, NU, T>{
        .A = sys.A - (sys.B * K) - (L * sys.C) + (L * sys.D * K),
        .B = L,
        .C = -K,
        .D = Matrix<NU, NY, T>::zeros(),
        .Ts = sys.Ts,
    };
}

namespace detail {

/** @brief MATLAB®-shaped (A, B, p) as the design::place arguments: a continuous plant and complex poles. */
template<size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr StateSpace<NX, NU, NX, T> place_plant(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B) {
    return StateSpace<NX, NU, NX, T>{.A = A, .B = B, .C = Matrix<NX, NX, T>::identity()};
}

/** @brief Poles from any tuple-like of damp::complex / std::complex / real values. */
template<size_t NX, typename T>
[[nodiscard]] constexpr damp::array<damp::complex<T>, NX> place_poles(const auto& p) {
    const auto poles = std::apply(
        [](const auto&... elems) { return damp::array<damp::complex<T>, sizeof...(elems)>{damp::complex<T>(elems)...}; }, p
    );
    static_assert(poles.size() == NX, "place / acker need one pole per state");
    return poles;
}

/**
 * @brief MATLAB®'s Q on [x; xi] as the design's Q and Q_integral; nullopt when its x–xi blocks
 *        are not zero (design::lqi weights the two blocks apart).
 */
template<size_t NX, size_t NI, size_t NU, typename T>
[[nodiscard]] constexpr damp::optional<design::LQRWeights<NX, NU, T, NI>> split_integral_weight(const Matrix<NX + NI, NX + NI, T>& Q_aug, const Matrix<NU, NU, T>& R) {
    for (size_t i = 0; i < NX; ++i) {
        for (size_t j = 0; j < NI; ++j) {
            if (Q_aug(i, NX + j) != T{0} || Q_aug(NX + j, i) != T{0}) {
                return damp::nullopt;
            }
        }
    }
    return design::LQRWeights<NX, NU, T, NI>{.Q = Q_aug.template block<NX, NX>(0, 0), .R = R, .Q_integral = Q_aug.template block<NI, NI>(NX, NX)};
}

/** @brief An index list as the design's IntegralAction. */
template<size_t NI>
[[nodiscard]] constexpr design::IntegralAction<NI> integrated_list(const size_t (&outputs)[NI]) {
    design::IntegralAction<NI> sel{};
    for (size_t i = 0; i < NI; ++i) {
        sel.outputs[i] = outputs[i];
    }
    return sel;
}

/** @brief The LQ design with integral action on a discrete plant, with MATLAB®'s Q on [x; xi]: bad_weights for x–xi cross weights. */
template<size_t NI, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, typename Options>
[[nodiscard]] constexpr auto integral_design(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const Matrix<NX + NI, NX + NI, T>& Q_aug, const Matrix<NU, NU, T>& R, const Options& options) {
    const auto w = split_integral_weight<NX, NI>(Q_aug, R);
    auto       res = design::detail::state_feedback(sys, w ? *w : design::LQRWeights<NX, NU, T, NI>{}, options);
    if (!w) {
        res = ::damp::detail::refused<decltype(res)>(design::DesignStatus::bad_weights); // x–xi cross weights
    }
    return res;
}

} // namespace detail

/**
 * @brief Single-input pole placement (MATLAB®'s acker)
 *
 * Forwards to design::place, which gives the unique single-input gain and assigns
 * repeated poles as one Jordan chain, the case Ackermann's formula covers. The
 * eigenvalues are read as continuous (MATLAB® place / acker carry no time domain);
 * for a discrete law call design::place on the discrete StateSpace.
 *
 * @param A State matrix
 * @param B Input matrix (one column)
 * @param p Desired poles (damp::array / std::array of damp::complex or std::complex)
 * @return design::PlaceResult with K, the achieved poles e and success / status
 *
 * @note Compare with MATLAB®'s K = acker(A, B, p).
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr design::PlaceResult<NX, NU, T> acker(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const auto&              p
) {
    static_assert(NU == 1, "acker is single-input only; use place for multi-input systems");
    return design::place(detail::place_plant(A, B), detail::place_poles<NX, T>(p));
}

/**
 * @brief Robust multi-input pole placement (MATLAB®'s place)
 *
 * Forwards to design::place: Kautsky–Nichols–Van Dooren robust eigenvalue assignment,
 * spending the multi-input freedom to minimize eigenvector conditioning. The
 * eigenvalues are read as continuous, as for acker.
 *
 * @param A State matrix
 * @param B Input matrix
 * @param p Desired poles (damp::array / std::array of damp::complex or std::complex)
 * @return design::PlaceResult with K, the achieved poles e and success / status
 *
 * @note Compare with MATLAB®'s K = place(A, B, p).
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr design::PlaceResult<NX, NU, T> place(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const auto&              p
) {
    return design::place(detail::place_plant(A, B), detail::place_poles<NX, T>(p));
}

/**
 * @brief Continuous-time LQR design (MATLAB®'s lqr)
 *
 * Thin alias for design::continuous_lqr — the optimal gain K for u = −Kx
 * minimizing ∫(xᵀQx + uᵀRu + 2xᵀNu) dt via the continuous ARE (care()).
 *
 * @return design::LQRResult with K, the CARE solution S, closed-loop poles e, and success / status
 * @note Compare with MATLAB®'s [K,S,e] = lqr(A, B, Q, R, N). MATLAB®'s several
 *       outputs are one design::LQRResult here.
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr design::LQRResult<NX, NU, T> lqr(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    const Matrix<NX, NU, T>& N = Matrix<NX, NU, T>{}
) {
    return design::continuous_lqr(A, B, Q, R, N);
}

/**
 * @brief Discrete-time Linear-Quadratic Regulator design (MATLAB®'s dlqr)
 * @note Forwards to design::discrete_lqr. Compare with MATLAB®'s [K,S,e] = dlqr(A, B, Q, R, N).
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr design::LQRResult<NX, NU, T> dlqr(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    const Matrix<NX, NU, T>& N = Matrix<NX, NU, T>{}
) {
    return design::discrete_lqr(A, B, Q, R, N);
}

/**
 * @brief Sampled-data LQR from a continuous plant: ZOH dynamics + Van Loan cost (MATLAB®'s lqrd)
 * @note Forwards to design::lqr(sys, Ts, w). Compare with MATLAB®'s lqrd(A, B, Q, R, Ts, N).
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr auto lqrd(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    T                        Ts,
    const Matrix<NX, NU, T>& N = Matrix<NX, NU, T>{}
) {
    return design::lqr(StateSpace<NX, NU, 1, T>{.A = A, .B = B}, Ts, {Q, R, N});
}

/**
 * @brief Sampled-data LQR from a continuous state-space plant (MATLAB®'s lqrd)
 * @note Forwards to design::lqr(sys, Ts, w). Compare with MATLAB®'s lqrd(sys, Q, R, Ts).
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr auto lqrd(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NX, NX, T>&                 Q,
    const Matrix<NU, NU, T>&                 R,
    T                                        Ts,
    const Matrix<NX, NU, T>&                 N = Matrix<NX, NU, T>{}
) {
    return design::lqr(sys, Ts, {Q, R, N});
}

/**
 * @brief Discrete linear-quadratic estimator from matrices (MATLAB®'s dlqe)
 *
 * Forwards to design::kalman on a carrier plant: process noise through @p G, H = I,
 * B = D = 0 (they do not enter the gain) and Ts = 1 (the matrices are already discrete).
 * The result's @c sys is that carrier: deploy the observer from the real plant and this M.
 *
 * @note Compare with MATLAB®'s [M,P,Z,E] = dlqe(A, G, C, Q, R): M (the innovation gain) and P
 *       are on the result, with the predictor gain L = A M.
 */
template<size_t NX, size_t NW, size_t NY, typename T = double>
[[nodiscard]] constexpr design::KalmanResult<NX, 1, NY, T, NW, NY> dlqe(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NW, T>& G,
    const Matrix<NY, NX, T>& C,
    const Matrix<NW, NW, T>& Q,
    const Matrix<NY, NY, T>& R
) {
    return design::kalman(StateSpace<NX, 1, NY, T, NW, NY>{.A = A, .C = C, .G = G, .H = Matrix<NY, NY, T>::identity(), .Ts = T{1}}, {.process = Q, .measurement = R});
}

/**
 * @brief Discrete Kalman estimator for a continuous plant and continuous noise (MATLAB®'s kalmd)
 *
 * Keeps MATLAB's convention: @p Rn is a continuous intensity, and the measurement variance
 * per reading is Rd = Rn/Ts (an integrate-and-dump sampler). Converts, then forwards to
 * design::kalman(sys, Ts, {.process = Qn, .measurement = Rn/Ts}): ZOH plant, Van Loan process
 * noise. Qn enters through the plant's G and Rn through its H. design::kalman itself takes
 * the per-reading variance directly (see design::KalmanNoise).
 *
 * @note Compare with MATLAB®'s [kest,L,P,M] = kalmd(sys, Qn, Rn, Ts): L (the predictor gain),
 *       P and M (the innovation gain) are on the result.
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr design::KalmanResult<NX, NU, NY, T, NX, NV> kalmd(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NW, NW, T>&                 Qn,
    const Matrix<NV, NV, T>&                 Rn,
    std::type_identity_t<T>                  Ts
) {
    const T per_reading = damp::finite_positive(Ts) ? T{1} / Ts : T{1}; // a bad Ts is refused by design::kalman
    return design::kalman(sys, Ts, {.process = Qn, .measurement = Matrix<NV, NV, T>(Rn * per_reading)});
}

/**
 * @brief Output-weighted continuous LQR (state cost Q = Cᵀ Q_y C).
 *
 * Builds the equivalent state-feedback LQR problem from output weights:
 * @f[
 *   Q_x = C^\top Q_y C,\quad
 *   N = C^\top Q_y D,\quad
 *   R_u = R + D^\top Q_y D
 * @f]
 * then calls @ref design::continuous_lqr. When @p D is zero this reduces to
 * pure output weighting @f$Q_x = C^\top Q_y C@f$.
 *
 * @note Compare with MATLAB®'s [K,S,e] = lqry(A, B, C, D, Qy, R).
 * @see design::continuous_lqr, dlqr, lqrd
 */
template<size_t NX, size_t NU, size_t NY, typename T = double>
[[nodiscard]] constexpr design::LQRResult<NX, NU, T> lqry(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NY, NX, T>& C,
    const Matrix<NY, NU, T>& D,
    const Matrix<NY, NY, T>& Qy,
    const Matrix<NU, NU, T>& R
) {
    const Matrix<NX, NX, T> Qx = C.transpose() * Qy * C;
    const Matrix<NX, NU, T> N = C.transpose() * Qy * D;
    const Matrix<NU, NU, T> Ru = R + (D.transpose() * Qy * D);
    return design::continuous_lqr(A, B, Qx, Ru, N);
}

/**
 * @brief Output-weighted LQR for a state-space plant.
 *
 * Follows the plant's time domain through @ref design::lqr(): CARE when
 * @p sys.Ts == 0, DARE when @p sys.Ts > 0. Same @f$Q_x/N/R_u@f$ map as the
 * matrix overload.
 *
 * @note Compare with MATLAB®'s [K,S,e] = lqry(sys, Qy, R).
 * @see lqry(A,B,C,D,Qy,R), design::lqr
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr auto lqry(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NY, NY, T>&                 Qy,
    const Matrix<NU, NU, T>&                 R
) {
    const Matrix<NX, NX, T> Qx = sys.C.transpose() * Qy * sys.C;
    const Matrix<NX, NU, T> N = sys.C.transpose() * Qy * sys.D;
    const Matrix<NU, NU, T> Ru = R + (sys.D.transpose() * Qy * sys.D);
    return design::lqr(sys, {Qx, Ru, N});
}

/**
 * @brief Linear-Quadratic Integral design for tracking
 *
 * Discrete plants only (@p sys.Ts > 0): a continuous @p sys is refused with
 * DesignStatus::invalid_sample_time rather than a DARE solution of the wrong problem.
 *
 * Q is MATLAB®'s, on [x; xi], (NX+NI)×(NX+NI): its diagonal blocks are the design's Q and
 * Q_integral. Cross weights between x and xi are refused with DesignStatus::bad_weights.
 * Integrates every output. To integrate a subset, add an index list, `lqi(sys, Q, R, {0, 3})`,
 * or an NI×NY selector Ci.
 *
 * @note Alias for design::lqi. Compare with MATLAB®'s lqi(sys, Q, R).
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr auto lqi(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NX + NY, NX + NY, T>&       Q,
    const Matrix<NU, NU, T>&                 R
) {
    return detail::integral_design<NY>(sys, Q, R, design::detail::StateFeedbackOptions{.offset_free = design::IntegralAction{}});
}

/** @brief LQI integrating the outputs @p Ci selects (NI×NY). Discrete plants only. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NI>
[[nodiscard]] constexpr auto lqi(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NX + NI, NX + NI, T>&       Q,
    const Matrix<NU, NU, T>&                 R,
    const Matrix<NI, NY, T>&                 Ci
) {
    return detail::integral_design<NI>(sys, Q, R, design::detail::StateFeedbackOptions{.offset_free = Ci});
}

/** @brief LQI integrating the listed outputs, e.g. `lqi(sys, Q, R, {0, 3})`. Discrete plants only. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NI>
[[nodiscard]] constexpr auto lqi(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NX + NI, NX + NI, T>&       Q,
    const Matrix<NU, NU, T>&                 R,
    const size_t (&outputs)[NI]
) {
    return detail::integral_design<NI>(sys, Q, R, design::detail::StateFeedbackOptions{.offset_free = detail::integrated_list(outputs)});
}

/**
 * @brief Linear-Quadratic-Gaussian regulator design
 *
 * Discrete plants only (@p sys.Ts > 0); a continuous @p sys is refused with
 * DesignStatus::invalid_sample_time.
 *
 * @note Alias for design::lqg. Compare with MATLAB®'s lqg(sys, ...).
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr auto lqg(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NX, NX, T>&                 Q_lqr,
    const Matrix<NU, NU, T>&                 R_lqr,
    const Matrix<NW, NW, T>&                 Q_kf,
    const Matrix<NV, NV, T>&                 R_kf,
    const Matrix<NX, NU, T>&                 N = Matrix<NX, NU, T>{}
) {
    return design::lqg(sys, {Q_lqr, R_lqr, N}, {Q_kf, R_kf});
}

/**
 * @brief Combine separate Kalman filter and LQR designs into an LQG controller
 * @note Alias for design::lqg_from_parts. Compare with MATLAB®'s lqgreg(kest, k).
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = NX, size_t NV = NY, typename Reg>
[[nodiscard]] constexpr auto lqgreg(
    const design::KalmanResult<NX, NU, NY, T, NW, NV>& kest,
    const Reg&                                         reg // design::lqr's result, or dlqr's
) {
    return design::lqg_from_parts(kest, reg);
}

/**
 * @brief Linear-Quadratic-Gaussian design with integral action for tracking
 *
 * Discrete plants only (@p sys.Ts > 0); a continuous @p sys is refused with
 * DesignStatus::invalid_sample_time.
 *
 * Q_aug is MATLAB®'s, on [x; xi], as for lqi: cross weights between x and xi are refused with
 * DesignStatus::bad_weights. Integrates every output. To integrate a subset, add an index
 * list, `lqgtrack(sys, Q_aug, R, Q_kf, R_kf, {0, 3})`, or an NI×NY selector Ci.
 *
 * @note Alias for design::lqgi. Compare with
 *       MATLAB®'s lqgtrack(...).
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr auto lqgtrack(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NX + NY, NX + NY, T>&       Q_aug,
    const Matrix<NU, NU, T>&                 R,
    const Matrix<NW, NW, T>&                 Q_kf,
    const Matrix<NV, NV, T>&                 R_kf
) {
    return detail::integral_design<NY>(sys, Q_aug, R, design::detail::StateFeedbackOptions{.noise = design::KalmanNoise{Q_kf, R_kf}, .offset_free = design::IntegralAction{}});
}

/** @brief LQG servo integrating the outputs @p Ci selects (NI×NY). Discrete plants only. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NI>
[[nodiscard]] constexpr auto lqgtrack(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NX + NI, NX + NI, T>&       Q_aug,
    const Matrix<NU, NU, T>&                 R,
    const Matrix<NW, NW, T>&                 Q_kf,
    const Matrix<NV, NV, T>&                 R_kf,
    const Matrix<NI, NY, T>&                 Ci
) {
    return detail::integral_design<NI>(sys, Q_aug, R, design::detail::StateFeedbackOptions{.noise = design::KalmanNoise{Q_kf, R_kf}, .offset_free = Ci});
}

/**
 * @brief LQG servo integrating the listed outputs, e.g. `lqgtrack(sys, Q_aug, R, Q_kf, R_kf, {0, 3})`.
 * Discrete plants only.
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NI>
[[nodiscard]] constexpr auto lqgtrack(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NX + NI, NX + NI, T>&       Q_aug,
    const Matrix<NU, NU, T>&                 R,
    const Matrix<NW, NW, T>&                 Q_kf,
    const Matrix<NV, NV, T>&                 R_kf,
    const size_t (&outputs)[NI]
) {
    return detail::integral_design<NI>(sys, Q_aug, R, design::detail::StateFeedbackOptions{.noise = design::KalmanNoise{Q_kf, R_kf}, .offset_free = detail::integrated_list(outputs)});
}

/**
 * @brief Plant-aware PID tune at crossover ωc (MATLAB® pidtune spelling): 60° PID
 *
 * Forwards to @ref design::pidtune with `{.wc = wc}`; pass a design::PIDTuneSpec to
 * design::pidtune for the phase margin, law, Tf and 2-DOF weights.
 *
 * @note Compare with MATLAB®'s pidtune(sys, 'PID', wc).
 */
template<size_t NX, typename T = double>
[[nodiscard]] constexpr design::PIDGains<T> pidtune(const StateSpace<NX, 1, 1, T>& sys, T wc) noexcept {
    return design::pidtune(sys, design::PIDTuneSpec<T>{.wc = wc});
}

// ===========================================================================
// Frequency-domain analysis — MATLAB® spellings over damp::analysis
// ===========================================================================

// These already carry their MATLAB® names in analysis::; surface them under
// matlab:: too so a MATLAB®-style call site finds them in one namespace.
using analysis::bode;
using analysis::damp;
using analysis::dcgain;
using analysis::impulse;
using analysis::initial;
using analysis::linspace;
using analysis::lsim;
using analysis::lsiminfo;
using analysis::nichols;
using analysis::nyquist;
using analysis::pzmap;
using analysis::rlocus;
using analysis::sigma;
using analysis::step;
using analysis::stepinfo;

/**
 * @brief @p n points from 10ᵃ to 10ᵇ, logarithmically spaced (decade exponents, as MATLAB®'s)
 *
 * analysis::geomspace(10ᵃ, 10ᵇ, n): the endpoints are exact powers of ten.
 *
 * @note Compare with MATLAB®'s y = logspace(a, b, n).
 */
template<typename T = double>
    requires std::is_floating_point_v<T>
[[nodiscard]] std::vector<T> logspace(T a, T b, size_t n) {
    return analysis::geomspace(damp::pow(T{10}, a), damp::pow(T{10}, b), n);
}

/**
 * @brief MATLAB® short alias for the open-loop poles of a system.
 *
 * nullopt when the eigensolver does not converge or A has a non-finite entry.
 *
 * @note Compare with MATLAB®'s p = pole(sys). (analysis names it `poles`.)
 */
template<size_t NX, typename T = double>
[[nodiscard]] constexpr damp::optional<damp::array<damp::complex<T>, NX>> pole(const Matrix<NX, NX, T>& A) {
    return analysis::poles(A);
}

/**
 * @brief Gain/phase margins and their crossover frequencies.
 *
 * @note Compare with MATLAB®'s [Gm, Pm, Wcg, Wcp] = margin(...). Gm is a linear
 *       ratio; numeric_limits<T>::max() (the
 *       library's finite stand-in for +inf under -ffast-math) marks a missing crossover.
 *       With several crossovers, the smallest margins, as MATLAB®'s.
 *
 * @c Wcg is the phase-crossover frequency (GM site); @c Wcp is the
 * gain-crossover frequency (PM site).
 */
template<typename T = double>
struct MarginResult {
    T Gm{};  ///< Gain margin (linear ratio); +inf if phase never crosses -180°
    T Pm{};  ///< Phase margin (degrees); +inf if gain never crosses 0 dB
    T Wcg{}; ///< Phase-crossover frequency (rad/s, where phase = -180°) — GM site
    T Wcp{}; ///< Gain-crossover frequency (rad/s, where |G| = 1) — PM site
};

/**
 * @brief Gain and phase margins of a SISO loop over a frequency grid.
 *
 * Thin composition of analysis::bode and BodeResult::gain_margin /
 * phase_margin. Unlike MATLAB®'s margin(sys), the grid is explicit — pass
 * e.g. matlab::logspace(...) — since damp does no automatic frequency gridding.
 *
 * @param sys   SISO state-space loop (continuous or discrete)
 * @param omega Frequency grid (rad/s)
 * @return MarginResult{Gm, Pm, Wcg, Wcp}
 *
 * @note Compare with MATLAB®'s [Gm, Pm, Wcg, Wcp] = margin(sys, w).
 */
template<size_t NX, size_t NW, size_t NV, typename T>
[[nodiscard]] MarginResult<T> margin(const StateSpace<NX, 1, 1, T, NW, NV>& sys, const std::vector<T>& omega) {
    const auto      bode_data = analysis::bode(sys, omega);
    constexpr T     inf = std::numeric_limits<T>::max();
    MarginResult<T> m{inf, inf, T{0}, T{0}};
    if (const auto gm = bode_data.gain_margin()) {
        m.Gm = damp::pow(T{10}, gm->first / T{20}); // dB -> linear ratio
        m.Wcg = gm->second;
    }
    if (const auto pm = bode_data.phase_margin()) {
        m.Pm = pm->first;
        m.Wcp = pm->second;
    }
    return m;
}

/**
 * @brief All classical margins including delay margin (superset of @ref margin).
 *
 * Delay margin @f$D_m = P_m\mathrm{(rad)} / \omega_{gc}@f$ where @f$\omega_{gc}@f$
 * is the gain-crossover frequency stored in @ref MarginResult::Wcp. Missing
 * crossovers keep the corresponding margin at +inf (same sentinel as @ref margin).
 *
 * @note Compare with MATLAB®'s S = allmargin(sys) (subset: the smallest GM and PM and the
 *       DM at the PM's crossover; multi-crossing arrays are not collected).
 * @see margin
 */
template<typename T = double>
struct AllMarginResult {
    T Gm{};  ///< Gain margin (linear ratio); +inf if no phase crossover
    T Pm{};  ///< Phase margin (degrees); +inf if no gain crossover
    T Dm{};  ///< Delay margin (seconds); +inf if no gain crossover or ω = 0
    T Wcg{}; ///< Phase-crossover frequency (rad/s) — GM site
    T Wcp{}; ///< Gain-crossover frequency (rad/s) — PM / DM site
};

/**
 * @brief Gain, phase, and delay margins of a SISO loop over a frequency grid.
 *
 * Thin extension of margin: delay margin is the classical conversion of
 * phase margin at the gain crossover. Grid is explicit (pass @ref logspace).
 *
 * @param sys   SISO state-space loop (continuous or discrete)
 * @param omega Frequency grid (rad/s)
 * @return AllMarginResult{Gm, Pm, Dm, Wcg, Wcp}
 *
 * @note Compare with MATLAB®'s S = allmargin(sys).
 */
template<size_t NX, size_t NW, size_t NV, typename T>
[[nodiscard]] AllMarginResult<T> allmargin(const StateSpace<NX, 1, 1, T, NW, NV>& sys, const std::vector<T>& omega) {
    const MarginResult<T> m = margin(sys, omega);
    constexpr T           inf = std::numeric_limits<T>::max();
    AllMarginResult<T>    a{m.Gm, m.Pm, inf, m.Wcg, m.Wcp};
    if (m.Pm < inf && m.Wcp > T{0}) {
        a.Dm = (m.Pm * damp::numbers::pi_v<T> / T{180}) / m.Wcp;
    }
    return a;
}

/**
 * @brief Disk margin of a SISO loop.
 * @note Forwards to analysis::disk_margin. Compare with MATLAB®'s DM = diskmargin(L, sigma):
 *       DiskMargin = alpha, GainMargin = [gain_low gain_high], PhaseMargin = ±phase_margin_deg,
 *       Frequency = frequency.
 */
template<size_t NX, size_t NW, size_t NV, typename T>
[[nodiscard]] constexpr analysis::DiskMarginResult<T> diskmargin(const StateSpace<NX, 1, 1, T, NW, NV>& L, std::type_identity_t<T> sigma = T{0}) {
    return analysis::disk_margin(L, sigma);
}

/**
 * @brief Loop-at-a-time disk margins of a square MIMO loop.
 * @note Forwards to analysis::loop_at_a_time_disk_margin. Compare with MATLAB®'s
 *       DM = diskmargin(L, sigma) for a MIMO L: DM(i) is channel[i]. The multiloop margin MM
 *       needs μ and is not computed.
 */
template<size_t NX, size_t N, size_t NW, size_t NV, typename T>
    requires(N > 1)
[[nodiscard]] constexpr analysis::LoopDiskMarginResult<N, T> diskmargin(const StateSpace<NX, N, N, T, NW, NV>& L, std::type_identity_t<T> sigma = T{0}) {
    return analysis::loop_at_a_time_disk_margin(L, sigma);
}

/**
 * @brief Stability predicate for a state-space model.
 *
 * Continuous (@p sys.Ts == 0): left half-plane. Discrete (@p sys.Ts > 0): unit disk.
 *
 * @note Compare with MATLAB®'s tf = isstable(sys).
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr bool isstable(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    return stability::is_stable(sys);
}

/**
 * @brief Stability predicate for a continuous transfer function.
 *
 * Realizes @p tf as companion state-space and tests the poles.
 *
 * @note Compare with MATLAB®'s isstable(sys).
 */
template<size_t Nnum, size_t Nden, typename T = double>
[[nodiscard]] constexpr bool isstable(const TransferFunction<Nnum, Nden, T>& tf) {
    const auto ss = tf.to_state_space();
    return ss ? isstable(*ss) : false;
}

/**
 * @brief Padé approximation of order N (1 or 2) of the pure delay e^{−sT}
 *
 * Forwards to design::pade_delay_1st / design::pade_delay_2nd. The order is a template
 * argument because the transfer function's size depends on it; a sampled delay is
 * design::pade_delay_1st(T, Ts) / design::pade_delay_2nd(T, Ts), or matlab::c2d.
 *
 * @note Compare with MATLAB®'s sys = pade(T, N), spelled matlab::pade<N>(T).
 * @see design::pade_delay_1st, design::pade_delay_2nd
 */
template<size_t N, typename T = double>
    requires(N == 1 || N == 2)
[[nodiscard]] constexpr auto pade(T T_delay) {
    if constexpr (N == 1) {
        return design::pade_delay_1st(T_delay);
    } else {
        return design::pade_delay_2nd(T_delay);
    }
}

/**
 * @brief -3 dB bandwidth of a SISO system over a frequency grid.
 *
 * @param sys   SISO state-space system
 * @param omega Frequency grid (rad/s), ascending from ~DC
 * @return Bandwidth (rad/s), or nullopt if the response never drops 3 dB.
 *
 * @note Compare with MATLAB®'s fb = bandwidth(sys). Grid is explicit here.
 */
template<size_t NX, size_t NW, size_t NV, typename T>
[[nodiscard]] damp::optional<T> bandwidth(const StateSpace<NX, 1, 1, T, NW, NV>& sys, const std::vector<T>& omega) {
    return analysis::bode(sys, omega).bandwidth();
}

/**
 * @brief MATLAB® alias for the continuous Lyapunov solve AX+XAᵀ+Q=0.
 * @note Compare with MATLAB®'s @c X=lyap(A,Q). Returns nullopt if no unique solution.
 */
template<size_t NX, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NX, T>>
lyap(const Matrix<NX, NX, T>& A, const Matrix<NX, NX, T>& Q) {
    return damp::lyap(A, Q);
}

/**
 * @brief MATLAB® alias for the discrete Lyapunov solve AXAᵀ−X+Q=0.
 * @note Compare with MATLAB®'s @c X=dlyap(A,Q). Returns nullopt if no unique solution.
 */
template<size_t NX, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NX, T>>
dlyap(const Matrix<NX, NX, T>& A, const Matrix<NX, NX, T>& Q) {
    return damp::dlyap(A, Q);
}

/**
 * @brief MATLAB® alias for the Cholesky factor of the continuous Lyapunov solution (X = RᵀR).
 * @note Compare with MATLAB®'s @c R=lyapchol(A,B). Returns nullopt unless A is Hurwitz.
 */
template<size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NX, T>>
lyapchol(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B) {
    return damp::lyapunov_cholesky(A, B);
}

/**
 * @brief MATLAB® alias for the Cholesky factor of the discrete Lyapunov solution (X = RᵀR).
 * @note Compare with MATLAB®'s @c R=dlyapchol(A,B). Returns nullopt unless A is Schur stable.
 */
template<size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NX, T>>
dlyapchol(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B) {
    return damp::discrete_lyapunov_cholesky(A, B);
}

/**
 * @brief MATLAB® alias for the dense inequality-constrained QP solve.
 * @note Compare with MATLAB®'s @c x=quadprog(H,f,A,b). Returns only the
 *       minimizer (@c nullopt unless the solve succeeded); for multipliers,
 *       objective, and status call design::solve_qp directly.
 */
template<size_t NV, size_t NI, typename T = double>
[[nodiscard]] constexpr damp::optional<ColVec<NV, T>> quadprog(
    const Matrix<NV, NV, T>& H,
    const ColVec<NV, T>&     f,
    const Matrix<NI, NV, T>& A,
    const ColVec<NI, T>&     b
) {
    const auto res = design::solve_qp(H, f, A, b);
    if (!res.success) {
        return damp::nullopt;
    }
    return res.x;
}

/**
 * @brief MATLAB® alias for the equality-constrained QP solve.
 * @note Compare with MATLAB®'s @c x=quadprog(H,f,A,b,Aeq,beq). Returns only the
 *       minimizer (@c nullopt unless the solve succeeded); call design::solve_qp
 *       directly for multipliers and status.
 */
template<size_t NV, size_t NI, size_t NE, typename T = double>
[[nodiscard]] constexpr damp::optional<ColVec<NV, T>> quadprog(
    const Matrix<NV, NV, T>& H,
    const ColVec<NV, T>&     f,
    const Matrix<NI, NV, T>& A,
    const ColVec<NI, T>&     b,
    const Matrix<NE, NV, T>& Aeq,
    const ColVec<NE, T>&     beq
) {
    const auto res = design::solve_qp(H, f, A, b, Aeq, beq);
    if (!res.success) {
        return damp::nullopt;
    }
    return res.x;
}

/**
 * @brief MATLAB® alias for the controllability/observability Gramian of a system.
 *
 * @param sys  State-space system (continuous or discrete, dispatched on Ts).
 * @param kind @c 'c' for the controllability Gramian, @c 'o' for observability.
 * @return Gramian (NX × NX), or nullopt if the Lyapunov solve fails.
 * @note Compare with MATLAB®'s @c gram(sys,'c') / @c gram(sys,'o').
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NX, T>>
gram(const StateSpace<NX, NU, NY, T, NW, NV>& sys, char kind = 'c') {
    if (kind == 'o' || kind == 'O') {
        return stability::observability_gramian(sys.A, sys.C, sys.is_discrete());
    }
    return stability::controllability_gramian(sys.A, sys.B, sys.is_discrete());
}

/**
 * @brief MATLAB® alias for the H2 system norm @c norm(sys,2).
 * @note Returns nullopt if infinite (non-strictly-proper continuous system) or the
 *       Gramian solve fails. @see analysis::norm_h2.
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr damp::optional<T> norm(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    return analysis::norm_h2(sys);
}

/**
 * @brief MATLAB® alias for the H∞ system norm @c [ninf, fpeak] = hinfnorm(sys, tol).
 * @note Forwards to analysis::norm_hinf: ninf = gain, fpeak = frequency [rad/s]; an unstable
 *       sys returns success = false. The default tol is 1e-12 (1e-6 for float), not MATLAB®'s 0.01.
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr analysis::HinfNormResult<T> hinfnorm(const StateSpace<NX, NU, NY, T, NW, NV>& sys, std::type_identity_t<T> tol = default_tol<T>()) {
    return analysis::norm_hinf(sys, tol);
}

} // namespace matlab

} // namespace damp