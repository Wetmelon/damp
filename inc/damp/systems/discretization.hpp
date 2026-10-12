// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file discretization.hpp
 * @brief Continuous-to-discrete state-space discretization
 */

#include <cstddef>

#include "damp/backend.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "state_space.hpp"

namespace damp {

/**
 * @brief Discretization methods for continuous-time state-space systems
 *
 * @see "Feedback Control of Dynamic Systems" (Franklin et al., 2015), Chapter 8
 */
enum class DiscretizationMethod {
    ForwardEuler, ///< Explicit Euler: Ad = I + ATs (first-order, simple, low overhead)
    ZOH,          ///< Zero-Order Hold: Ad = e^(ATs), exact for piecewise-constant inputs
    Tustin,       ///< Bilinear transform: s → (2/Ts)(z−1)/(z+1), preserves stability; pre-warp with discretize's last argument
};

namespace detail {

/**
 * @brief Discretize using Forward Euler (explicit Euler)
 *
 *     A_d = I + ATs,  B_d = BTs,  C_d = C,  D_d = D
 *
 * First-order approximation of the matrix exponential. Simple and fast but
 * only accurate when ‖A‖·Ts ≪ 1. Does not preserve stability for stiff systems.
 *
 * @see "Feedback Control of Dynamic Systems" (Franklin et al., 2015), §8.3
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr StateSpace<NX, NU, NY, T, NW, NV> discretize_forward_euler_impl(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        sampling_time
) {
    // A_d = I + A·Ts, B_d = B·Ts; C/D unchanged. Scale G like B when NW > 0.
    return StateSpace<NX, NU, NY, T, NW, NV>{
        .A = Matrix<NX, NX, T>::identity() + (sys.A * sampling_time),
        .B = sys.B * sampling_time,
        .C = sys.C,
        .D = sys.D,
        .G = sys.G * sampling_time,
        .H = sys.H,
        .Ts = sampling_time,
    };
}

/// Copy the NX × NC block of M starting at column c0 (rows 0 … NX−1).
template<size_t NX, size_t NC, size_t N, typename T>
[[nodiscard]] constexpr Matrix<NX, NC, T> leading_rows_block(const Matrix<N, N, T>& M, size_t c0) {
    Matrix<NX, NC, T> out{};
    for (size_t i = 0; i < NX; ++i) {
        for (size_t j = 0; j < NC; ++j) {
            out(i, j) = M(i, c0 + j);
        }
    }
    return out;
}

/**
 * @brief Discretize using Zero-Order Hold (ZOH)
 *
 *     A_d = e^(ATs)
 *     B_d = ∫₀^{Ts} e^{Aτ} B dτ
 *     G_d = ∫₀^{Ts} e^{Aτ} G dτ   (same map as B; process-noise input)
 *
 * Exact discretization assuming piecewise-constant input between samples. All
 * three come from one exponential of the augmented matrix (Van Loan, 1978):
 * @f[
 *   \exp\left(\begin{bmatrix} A & B & G \\ 0 & 0 & 0 \end{bmatrix} T_s\right)
 *   = \begin{bmatrix} A_d & B_d & G_d \\ 0 & I & 0 \\ 0 & 0 & I \end{bmatrix}.
 * @f]
 * Unlike @f$ A^{-1}(e^{AT_s} - I)B @f$ this has no cancellation when @f$ A @f$ is
 * singular or nearly so (integrators, very slow modes).
 *
 * @note Compare with MATLAB®'s c2d(sys, Ts, 'zoh').
 * @see "Computing integrals involving the matrix exponential" (Van Loan, 1978)
 * @see "Feedback Control of Dynamic Systems" (Franklin et al., 2015), §8.3
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr StateSpace<NX, NU, NY, T, NW, NV> discretize_zoh_impl(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        sampling_time
) {
    constexpr size_t N = NX + NU + NW;
    Matrix<N, N, T>  M{};
    for (size_t i = 0; i < NX; ++i) {
        for (size_t j = 0; j < NX; ++j) {
            M(i, j) = sys.A(i, j) * sampling_time;
        }
        for (size_t j = 0; j < NU; ++j) {
            M(i, NX + j) = sys.B(i, j) * sampling_time;
        }
        for (size_t j = 0; j < NW; ++j) {
            M(i, NX + NU + j) = sys.G(i, j) * sampling_time;
        }
    }
    const Matrix<N, N, T> E = mat::expm(M);

    return StateSpace<NX, NU, NY, T, NW, NV>{
        .A = leading_rows_block<NX, NX>(E, 0),
        .B = leading_rows_block<NX, NU>(E, NX),
        .C = sys.C,
        .D = sys.D,
        .G = leading_rows_block<NX, NW>(E, NX + NU),
        .H = sys.H,
        .Ts = sampling_time,
    };
}

/**
 * @brief Discretize a continuous-time state-space system using Tustin method
 *
 * Tustin (Bilinear Transform) implementation
 *    Maps: s → (2/Ts) · (z − 1) / (z + 1)
 *
 *    The bilinear map compresses the frequency axis: the discrete response at ω is
 *    the continuous one at (2/Ts)·tan(ωTs/2). With @p prewarp_frequency ω₀ > 0 the map
 *    is s → (ω₀ / tan(ω₀Ts/2)) · (z − 1)/(z + 1) instead, which makes the two responses
 *    agree exactly at ω₀ (a notch, a resonance, a crossover) at the cost of a little
 *    more warping elsewhere. That is the plain map with 2/Ts replaced, so it runs the
 *    same formulas with the effective period Ts' = 2·tan(ω₀Ts/2)/ω₀; the result still
 *    carries the true Ts. ω₀Ts must be below π (the Nyquist rate).
 *
 *    Let L = (I − A·Ts/2). All results are computed via LU solve
 *    against L rather than forming L⁻¹ explicitly:
 *      A_d:  solve L · A_d = (I + A·Ts/2)
 *      B_d:  solve L · X = B,  then B_d = Ts · X
 *      C_d:  solve Lᵀ · X = Cᵀ, then C_d = Xᵀ
 *      D_d = D + (Ts/2) · C_d · B
 *
 * @note Compare with MATLAB®'s c2d(sys, Ts, c2dOptions('Method', 'tustin', 'PrewarpFrequency', w0)).
 * @see "Feedback Control of Dynamic Systems" (Franklin et al., 2015), §8.6
 * @see "Discrete-Time Signal Processing" (Oppenheim & Schafer, 2010), §7.3 (bilinear transform, frequency warping)
 *
 * @param sys               Continuous-time state-space model
 * @param sampling_time     Desired sampling period for discrete system
 * @param prewarp_frequency [rad/s] frequency at which the discrete response must match; 0 for none
 * @return Discretized state-space system, nullopt if L is singular or ω₀Ts ≥ π
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr damp::optional<StateSpace<NX, NU, NY, T, NW, NV>> discretize_tustin_impl(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        sampling_time,
    T                                        prewarp_frequency = T{0}
) {
    T ts_eff = sampling_time;
    if (prewarp_frequency > T{0}) {
        const T half_angle = prewarp_frequency * sampling_time / T{2};
        if (!(half_angle < damp::numbers::pi_v<T> / T{2})) {
            return damp::nullopt; // at or past Nyquist: tan has no finite value
        }
        ts_eff = T{2} * damp::tan(half_angle) / prewarp_frequency;
    }
    const T      ts_half = ts_eff / T{2};
    const Matrix I = Matrix<NX, NX, T>::identity();
    const Matrix L = I - sys.A * ts_half;

    // A_d: solve L · A_d = (I + A·Ts/2)
    const auto A_d_opt = mat::lu_solve(L, I + sys.A * ts_half);
    if (!A_d_opt) {
        return damp::nullopt;
    }
    const Matrix A_d = A_d_opt.value();

    // B_d: solve L · X = B, then B_d = Ts · X
    const auto B_d_opt = mat::lu_solve(L, sys.B);
    if (!B_d_opt) {
        return damp::nullopt;
    }
    const Matrix B_d = B_d_opt.value() * ts_eff;

    // C_d = C · L⁻¹ → solve Lᵀ · X = Cᵀ, then C_d = Xᵀ
    const auto C_d_t_opt = mat::lu_solve(L.transpose(), sys.C.transpose());
    if (!C_d_t_opt) {
        return damp::nullopt;
    }
    const Matrix C_d = C_d_t_opt.value().transpose();

    // D_d = D + (Ts/2) · C_d · B
    const Matrix D_d = sys.D + C_d * sys.B * ts_half;

    // G_d: same map as B
    const auto G_d_opt = mat::lu_solve(L, sys.G);
    if (!G_d_opt) {
        return damp::nullopt;
    }
    const Matrix G_d = G_d_opt.value() * ts_eff;
    const Matrix H_d = sys.H;

    return StateSpace{A_d, B_d, C_d, D_d, G_d, H_d, sampling_time};
}

} // namespace detail

/**
 * @brief Discretize a continuous-time state-space system
 *
 * Converts @f$ \dot x = Ax + Bu @f$ to @f$ x[k+1] = A_d x[k] + B_d u[k] @f$.
 * Already-discrete models and non-positive @p sampling_time return @p sys unchanged.
 * Tustin returns @c nullopt if the bilinear factor @f$ I - A T_s/2 @f$ is singular
 * (no silent fallback to ZOH). Tustin warps the frequency axis; @p prewarp_frequency
 * ω₀ > 0 makes the discrete response exact at ω₀ (see detail::discretize_tustin_impl)
 * and is ignored by the other methods.
 *
 * @note Compare with MATLAB®'s c2d(sys, Ts, method) and c2dOptions('PrewarpFrequency', w0).
 * @see "Feedback Control of Dynamic Systems" (Franklin et al., 2015), Chapter 8
 *
 * @param sys               Continuous-time state-space model (@c Ts = 0)
 * @param sampling_time     Desired sampling period [s]
 * @param method            Discretization method (ForwardEuler, ZOH, or Tustin)
 * @param prewarp_frequency [rad/s] Tustin pre-warp frequency; 0 (default) for the plain bilinear map
 * @return Discrete system, or @c nullopt if the method cannot be applied
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr damp::optional<StateSpace<NX, NU, NY, T, NW, NV>> discretize(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        sampling_time,
    DiscretizationMethod                     method = DiscretizationMethod::ZOH,
    T                                        prewarp_frequency = T{0}
) {
    if (sampling_time <= T{0}) {
        return sys;
    }

    if (sys.Ts > T{0}) {
        return sys;
    }

    switch (method) {
        case DiscretizationMethod::ForwardEuler:
            return detail::discretize_forward_euler_impl(sys, sampling_time);
        case DiscretizationMethod::ZOH:
            return detail::discretize_zoh_impl(sys, sampling_time);
        case DiscretizationMethod::Tustin:
            return detail::discretize_tustin_impl(sys, sampling_time, prewarp_frequency);
        default:
            return detail::discretize_forward_euler_impl(sys, sampling_time);
    }
}

// ---- Input delay ------------------------------------------------------------

namespace detail {

/**
 * @brief Assemble a delay-line model: the plant state followed by D blocks of past inputs.
 *
 * State [x; u(k−1); …; u(k−D)]. @p Gamma0 multiplies u(k−lag0) and @p Gamma1 multiplies
 * u(k−lag0−1) (lag 0 is the current input, through B). The output feedthrough Dy reads the
 * input applied at the sample instant, u(k−lag_y).
 */
template<size_t D, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr StateSpace<NX + (NU * D), NU, NY, T, NW, NV> delay_line(
    const StateSpace<NX, NU, NY, T, NW, NV>& base, // A, C, G, H, Ts of the plant
    const Matrix<NX, NU, T>&                 Gamma0,
    size_t                                   lag0,
    const Matrix<NX, NU, T>&                 Gamma1,
    const Matrix<NY, NU, T>&                 Dy,
    size_t                                   lag_y
) {
    constexpr size_t                  NA = NX + (NU * D);
    StateSpace<NA, NU, NY, T, NW, NV> out{};
    const auto                        place_gamma = [&](const Matrix<NX, NU, T>& g, size_t lag) {
        if (lag == 0) {
            out.B.template block<NX, NU>(0, 0) = g;
        } else if (lag <= D) {
            out.A.template block<NX, NU>(0, NX + (NU * (lag - 1))) = g;
        }
    };
    out.A.template block<NX, NX>(0, 0) = base.A;
    place_gamma(Gamma0, lag0);
    place_gamma(Gamma1, lag0 + 1);
    if constexpr (D > 0) {
        out.B.template block<NU, NU>(NX, 0) = Matrix<NU, NU, T>::identity(); // u(k) enters the line
        for (size_t j = 1; j < D; ++j) {                                     // and shifts along it
            out.A.template block<NU, NU>(NX + (NU * j), NX + (NU * (j - 1))) = Matrix<NU, NU, T>::identity();
        }
    }
    out.C.template block<NY, NX>(0, 0) = base.C;
    if (lag_y == 0) {
        out.D = Dy;
    } else {
        out.C.template block<NY, NU>(0, NX + (NU * (lag_y - 1))) = Dy;
    }
    out.G.template block<NX, NW>(0, 0) = base.G; // the past inputs are known: no process noise
    out.H = base.H;
    out.Ts = base.Ts;
    return out;
}

} // namespace detail

/**
 * @brief Number of past-input blocks a delay needs: ⌈delay / Ts⌉ (0 for no delay).
 *
 * Sizes an InputDelay at compile time: `InputDelay<input_delay_steps(kDelay, kTs)>{kDelay}`.
 */
template<typename T>
[[nodiscard]] constexpr size_t input_delay_steps(T delay, T sampling_time) {
    if (!(delay > T{0}) || !(sampling_time > T{0})) {
        return 0;
    }
    const T      ratio = delay / sampling_time;
    const T      tol = T{16} * eps<T>() * damp::max(T{1}, ratio);
    const size_t whole = static_cast<size_t>(ratio + tol); // floor, tolerant of a delay that is a whole number of samples
    return (ratio - static_cast<T>(whole) > tol) ? whole + 1 : whole;
}

/**
 * @brief An input delay for discretize. @p D past-input blocks become states, so D is part of
 *        the type. `InputDelay<2>{}` is exactly two samples; give @ref seconds only for a
 *        fractional delay, with D = input_delay_steps(seconds, Ts).
 */
template<size_t D, typename T = double>
struct InputDelay {
    T seconds{T{0}}; ///< [s] delay between a command and the plant receiving it; 0: exactly D samples
};

/**
 * @brief Discretize with an input delay, exactly; the past inputs become states.
 *
 * The returned model carries @f$ x_a = [x;\ u_{k-1};\ \dots;\ u_{k-D}] @f$, so any controller
 * or estimator designed on it plans around the delay. Those states are commands already
 * sent: they are known, and get no process noise (G is zero on them).
 *
 * G is not delayed: put a physical disturbance (a load force, friction) in G of the
 * continuous plant and the discrete G is its exact ZOH entry on the plant states, ready to be
 * an OffsetFreeModel's Bd; padding B with it would delay it with the command.
 *
 * The plant states come first, and every design and runtime reads a weight, covariance or
 * state reference sized to the plant as that leading block: the delay states get no cost and
 * no process noise, and their reference is filled in (the past u*, the outputs the reference
 * produced). A state-feedback law still needs the whole state: the plant's plus the commands
 * already sent.
 *
 * Write the delay as @f$ d\,T_s + \tau @f$, @f$ 0 \le \tau < T_s @f$. A whole delay (τ = 0)
 * is a shift of the discretized input and works with any @p method, and on a model that is
 * already discrete (same Ts). A fractional one needs ZOH: within each sample the plant holds
 * u(k−d−1) until τ, then u(k−d), so
 * @f[
 *   x_{k+1} = \Phi x_k + \Gamma_0 u_{k-d} + \Gamma_1 u_{k-d-1},\qquad
 *   \Gamma_0 = \int_0^{T_s-\tau} e^{A\eta}\,d\eta\,B,\quad
 *   \Gamma_1 = e^{A(T_s-\tau)} \int_0^{\tau} e^{A\eta}\,d\eta\,B,
 * @f]
 * exact for a held input; no Padé approximation.
 *
 * @code
 * auto a = discretize(plant, 1e-3, InputDelay<2>{});       // exactly two samples
 * auto b = discretize(plant, 1e-3, InputDelay<2>{1.4e-3}); // 1.4 ms CAN round trip at 1 kHz
 * @endcode
 *
 * @note Compare with MATLAB®'s c2d on a model with InputDelay, then absorbDelay.
 * @see "Digital Control of Dynamic Systems" (Franklin, Powell, Workman, 1998), §4.3.3
 *
 * @return The augmented discrete model, or nullopt if D does not match the delay, the
 *         delay is negative, a fractional delay is asked of a non-ZOH method or of a
 *         discrete model, or a discrete model's Ts differs.
 */
template<size_t D, size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr damp::optional<StateSpace<NX + (NU * D), NU, NY, T, NW, NV>> discretize(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        sampling_time,
    InputDelay<D, T>                         delay,
    DiscretizationMethod                     method = DiscretizationMethod::ZOH
) {
    const T Ts = sampling_time;
    if (delay.seconds == T{0}) {
        delay.seconds = static_cast<T>(D) * Ts; // unspecified: exactly D samples
    }
    if (!(Ts > T{0}) || !(delay.seconds >= T{0}) || input_delay_steps(delay.seconds, Ts) != D) { // NaN refused too
        return damp::nullopt;
    }
    const T    ratio = delay.seconds / Ts;
    const T    tol = T{16} * eps<T>() * damp::max(T{1}, ratio);
    const auto d = static_cast<size_t>(ratio + tol);
    const T    tau = (ratio - static_cast<T>(d) > tol) ? (delay.seconds - (static_cast<T>(d) * Ts)) : T{0};

    if (sys.Ts > T{0} || tau == T{0}) { // whole samples: shift the discretized input
        if (sys.Ts > T{0} && (tau != T{0} || damp::abs(sys.Ts - Ts) > tol * Ts)) {
            return damp::nullopt;
        }
        const auto base = discretize(sys, Ts, method);
        if (!base) {
            return damp::nullopt;
        }
        return detail::delay_line<D>(*base, base->B, d, Matrix<NX, NU, T>{}, base->D, d);
    }
    if (method != DiscretizationMethod::ZOH) {
        return damp::nullopt;
    }
    const auto              full = detail::discretize_zoh_impl(sys, Ts);
    const auto              head = detail::discretize_zoh_impl(sys, tau);      // ∫₀^τ e^{Aη}dη B
    const auto              tail = detail::discretize_zoh_impl(sys, Ts - tau); // e^{A(Ts−τ)}, ∫₀^{Ts−τ} e^{Aη}dη B
    const Matrix<NX, NU, T> gamma1 = tail.A * head.B;
    return detail::delay_line<D>(full, tail.B, d, gamma1, sys.D, d + 1);
}

// ---- Output delay -----------------------------------------------------------

/**
 * @brief A sensor latency for discretize: the measurement read at sample k is the plant
 *        output of sample k − D. Whole samples; D is part of the type.
 *
 * By default every output is late. When only some sensors are (an encoder read over a field
 * bus next to a local scale), list them: `OutputDelay<2, 2>{{2, 3}}` delays outputs 2 and 3
 * and leaves the others current. Only the late outputs get register states.
 *
 * @tparam D  Delay in samples.
 * @tparam NS Number of late outputs listed in @ref rows; 0 (the default): every output.
 */
template<size_t D, size_t NS = 0>
struct OutputDelay {
    damp::array<size_t, NS> rows{}; ///< indices of the late outputs into y, each listed once
};

namespace detail {

/// Number of late outputs: the listed ones, or all NY when none are listed.
template<size_t NS, size_t NY>
inline constexpr size_t late_outputs = (NS == 0) ? NY : NS;

/// The late output indices: the listed ones, checked (in range, each once), or 0 … NY−1.
template<size_t NY, size_t D, size_t NS>
[[nodiscard]] constexpr damp::optional<damp::array<size_t, late_outputs<NS, NY>>> late_rows(const OutputDelay<D, NS>& delay) {
    static_cast<void>(delay); // only read when outputs are listed
    damp::array<size_t, late_outputs<NS, NY>> rows{};
    for (size_t i = 0; i < late_outputs<NS, NY>; ++i) {
        if constexpr (NS == 0) {
            rows[i] = i;
        } else {
            rows[i] = delay.rows[i];
        }
        if (rows[i] >= NY) {
            return damp::nullopt;
        }
        for (size_t j = 0; j < i; ++j) {
            if (rows[j] == rows[i]) {
                return damp::nullopt;
            }
        }
    }
    return rows;
}

/**
 * @brief Append a D-deep shift register of the late outputs to a discrete model.
 *
 * State [x; y_L(k−1); …; y_L(k−D)] with y_L the NL late outputs of y = Cx + Du; each late
 * output now reads its register tail, the others read the plant as before. The register is
 * deterministic (no process noise); measurement noise stays on the reading.
 */
template<size_t D, size_t NL, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr StateSpace<NX + (NL * D), NU, NY, T, NW, NV> output_delay_line(
    const StateSpace<NX, NU, NY, T, NW, NV>& base,
    const damp::array<size_t, NL>&           rows
) {
    if constexpr (D == 0 || NL == 0) {
        return base;
    } else {
        StateSpace<NX + (NL * D), NU, NY, T, NW, NV> out{};
        out.A.template block<NX, NX>(0, 0) = base.A;
        out.B.template block<NX, NU>(0, 0) = base.B;
        out.C.template block<NY, NX>(0, 0) = base.C; // current outputs read the plant
        out.D = base.D;
        constexpr size_t tail = NX + (NL * (D - 1));
        for (size_t i = 0; i < NL; ++i) {
            const size_t r = rows[i];
            for (size_t j = 0; j < NX; ++j) { // y_r(k) = C_r x + D_r u enters the register
                out.A(NX + i, j) = base.C(r, j);
                out.C(r, j) = T{0};
            }
            for (size_t j = 0; j < NU; ++j) {
                out.B(NX + i, j) = base.D(r, j);
                out.D(r, j) = T{0};
            }
            out.C(r, tail + i) = T{1}; // and is read D samples later
        }
        for (size_t j = 1; j < D; ++j) { // shifts along it
            out.A.template block<NL, NL>(NX + (NL * j), NX + (NL * (j - 1))) = Matrix<NL, NL, T>::identity();
        }
        out.G.template block<NX, NW>(0, 0) = base.G;
        out.H = base.H;
        out.Ts = base.Ts;
        return out;
    }
}

} // namespace detail

/**
 * @brief Discretize with a sensor latency, exactly; the past outputs become states.
 *
 * The returned model carries @f$ x_a = [x;\ y_{L,k-1};\ \dots;\ y_{L,k-D}] @f$, where
 * @f$ y_L @f$ are the late outputs (all of them by default), and measures them D samples
 * late while the other outputs stay current:
 * @f[
 *   x_{k+1} = A x_k + B u_k,\qquad
 *   y_{L,k-j+1} \leftarrow y_{L,k-j},\qquad
 *   y_{L,k} = C_L x_k + D_L u_k \text{ enters the register.}
 * @f]
 * An estimator designed on it reconstructs the current state from stale readings, and a
 * controller plans on that estimate. The register is deterministic (no process noise);
 * measurement noise stays on the reading. Works with any @p method and on a model that is
 * already discrete (same Ts). Every design family reads the result the same way: weights,
 * covariances and state references sized to the plant are its leading block.
 *
 * @code
 * auto a = discretize(plant, 1e-3, OutputDelay<1>{});                        // every sensor one sample late
 * auto b = discretize(plant, 1e-3, OutputDelay<2, 2>{{2, 3}});               // outputs 2 and 3 over a field bus
 * auto c = discretize(plant, 1e-3, InputDelay<2>{1.4e-3}, OutputDelay<1>{}); // and a CAN command link
 * @endcode
 *
 * @note Compare with MATLAB®'s c2d on a model with OutputDelay, then absorbDelay.
 * @see "Digital Control of Dynamic Systems" (Franklin, Powell, Workman, 1998), §4.3.3
 *
 * @return The augmented discrete model, or nullopt if Ts ≤ 0, @p method fails, a discrete
 *         model's Ts differs, or a listed output is out of range or listed twice.
 */
template<size_t D, size_t NS, size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr damp::optional<StateSpace<NX + (detail::late_outputs<NS, NY> * D), NU, NY, T, NW, NV>> discretize(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        sampling_time,
    const OutputDelay<D, NS>&                delay,
    DiscretizationMethod                     method = DiscretizationMethod::ZOH
) {
    const T    Ts = sampling_time;
    const auto rows = detail::late_rows<NY>(delay);
    if (!rows || !(Ts > T{0}) || (sys.Ts > T{0} && damp::abs(sys.Ts - Ts) > T{16} * eps<T>() * Ts)) {
        return damp::nullopt;
    }
    const auto base = discretize(sys, Ts, method);
    if (!base) {
        return damp::nullopt;
    }
    return detail::output_delay_line<D>(*base, *rows);
}

/**
 * @brief Discretize with both an input delay and a sensor latency, exactly.
 *
 * State @f$ [x;\ u_{k-1};\ \dots;\ u_{k-D_i};\ y_{L,k-1};\ \dots;\ y_{L,k-D_o}] @f$: the
 * input delay as in discretize(sys, Ts, InputDelay), then the output register on that model
 * (all outputs, or the ones @p output_delay lists).
 *
 * @code
 * auto m = discretize(plant, 1e-3, InputDelay<2>{1.4e-3}, OutputDelay<1>{});
 * auto n = discretize(plant, 1e-3, InputDelay<2>{2e-3}, OutputDelay<2, 2>{{2, 3}}); // a CAN drive: torque out, encoder back
 * @endcode
 *
 * @return The augmented discrete model, or nullopt where either delay alone would be refused.
 */
template<size_t DI, size_t DO, size_t NS, size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr damp::optional<StateSpace<NX + (NU * DI) + (detail::late_outputs<NS, NY> * DO), NU, NY, T, NW, NV>> discretize(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        sampling_time,
    InputDelay<DI, T>                        input_delay,
    const OutputDelay<DO, NS>&               output_delay,
    DiscretizationMethod                     method = DiscretizationMethod::ZOH
) {
    const auto rows = detail::late_rows<NY>(output_delay);
    const auto delayed = discretize(sys, sampling_time, input_delay, method);
    if (!rows || !delayed) {
        return damp::nullopt;
    }
    return detail::output_delay_line<DO>(*delayed, *rows);
}

} // namespace damp
