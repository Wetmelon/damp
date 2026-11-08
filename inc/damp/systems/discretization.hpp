// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file discretization.hpp
 * @brief c2d / d2d: discretize a plant; InputDelay in N registers; Tustin pre-warp
 */

#include <cstddef>
#include <limits>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "state_space.hpp"
#include "transfer_function.hpp"

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

/**
 * @brief A discretized or resampled model with its success / status, read like an optional
 *
 * `*r`, `r->A`, `r.value()`, `r.has_value()` and `if (r)` read it as they read a
 * @c damp::optional, so code that needs only the model is unchanged; @c status names the cause
 * when the model is refused. A refused result holds the zero model: `*r` returns it, and
 * @c value() does not compile when constant-evaluated and traps at run time (as
 * @c std::optional::value() does without exceptions).
 *
 * @tparam NX Number of states
 * @tparam NU Number of inputs
 * @tparam NY Number of outputs
 * @tparam T  Scalar type
 * @tparam NW Number of process-noise inputs
 * @tparam NV Number of measurement-noise inputs
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
struct [[nodiscard]] DiscretizationResult {
    using value_type = StateSpace<NX, NU, NY, T, NW, NV>; ///< The model type, as std::optional::value_type

    value_type           sys{};                                       ///< The discrete model (zero when refused)
    bool                 success{false};                              ///< status == DesignStatus::ok
    design::DesignStatus status{design::DesignStatus::not_converged}; ///< ok, or why the model was refused

    template<typename U>
    [[nodiscard]] constexpr DiscretizationResult<NX, NU, NY, U, NW, NV> as() const {
        return DiscretizationResult<NX, NU, NY, U, NW, NV>{sys.template as<U>(), success, status};
    }

    /** @brief True when the model was produced (success). */
    [[nodiscard]] constexpr bool has_value() const { return success; }

    /** @brief Contextual test, as std::optional: `if (r)`, `!r`, `r && …`. */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /** @brief The model (the zero model when refused). */
    [[nodiscard]] constexpr const value_type& operator*() const { return sys; }

    /** @brief The model (the zero model when refused). */
    [[nodiscard]] constexpr value_type& operator*() { return sys; }

    /** @brief Member access to the model: `r->A`. */
    [[nodiscard]] constexpr const value_type* operator->() const { return &sys; }

    /** @brief Member access to the model: `r->A`. */
    [[nodiscard]] constexpr value_type* operator->() { return &sys; }

    /** @brief The model; a failed result does not compile when constant-evaluated and traps at run time. */
    [[nodiscard]] constexpr const value_type& value() const {
        if (!success) {
            ::damp::detail::refuse_failed_design(status); // constant-evaluated: does not compile, naming the status; run time: traps
        }
        return sys;
    }

    /** @brief The model; a failed result does not compile when constant-evaluated and traps at run time. */
    [[nodiscard]] constexpr value_type& value() {
        if (!success) {
            ::damp::detail::refuse_failed_design(status); // constant-evaluated: does not compile, naming the status; run time: traps
        }
        return sys;
    }
};

namespace detail {

/** @brief A successful DiscretizationResult holding @p sys. */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr DiscretizationResult<NX, NU, NY, T, NW, NV> discretized(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    return DiscretizationResult<NX, NU, NY, T, NW, NV>{sys, true, design::DesignStatus::ok};
}

/**
 * @brief True when every entry of A, B and G is finite and below the overflow bound.
 *
 * damp::exp saturates at the largest finite value instead of returning ∞, so an exponential
 * that overflowed holds that value; it counts as not finite here.
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr bool finite_dynamics(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    const auto finite = [](T v) { return damp::isfinite(v) && damp::abs(v) < std::numeric_limits<T>::max(); };
    return mat::all_of(sys.A, finite) && mat::all_of(sys.B, finite) && mat::all_of(sys.G, finite);
}

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

/** @brief Copy the NX × NC block of M starting at column c0 (rows 0 … NX−1). */
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
[[nodiscard]] constexpr damp::optional<StateSpace<NX, NU, NY, T, NW, NV>> discretize_zoh_impl(
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
    const auto E = mat::expm(M);
    if (!E) {
        return damp::nullopt;
    }

    return StateSpace<NX, NU, NY, T, NW, NV>{
        .A = leading_rows_block<NX, NX>(*E, 0),
        .B = leading_rows_block<NX, NU>(*E, NX),
        .C = sys.C,
        .D = sys.D,
        .G = leading_rows_block<NX, NW>(*E, NX + NU),
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
 *    The bilinear image also has a w → y feedthrough (Ts/2)·C·L⁻¹·G; the (G, H) noise split
 *    has no slot for it, so it is dropped.
 *
 * @note Compare with MATLAB®'s c2d(sys, Ts, c2dOptions('Method', 'tustin', 'PrewarpFrequency', w0)).
 * @see "Feedback Control of Dynamic Systems" (Franklin et al., 2015), §8.6
 * @see "Discrete-Time Signal Processing" (Oppenheim & Schafer, 2010), §7.3 (bilinear transform, frequency warping)
 *
 * @param sys               Continuous-time state-space model
 * @param sampling_time     Desired sampling period for discrete system
 * @param prewarp_frequency [rad/s] frequency at which the discrete response must match; 0 for none
 * @return Discretized state-space system; refused with DesignStatus::singular if L is singular,
 *         DesignStatus::invalid_parameter if ω₀ is negative or not finite, or ω₀Ts ≥ π
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr DiscretizationResult<NX, NU, NY, T, NW, NV> discretize_tustin_impl(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        sampling_time,
    T                                        prewarp_frequency = T{0}
) {
    using Result = DiscretizationResult<NX, NU, NY, T, NW, NV>;
    T ts_eff = sampling_time;
    if (!damp::isfinite(prewarp_frequency) || prewarp_frequency < T{0}) {
        return refused<Result>(design::DesignStatus::invalid_parameter);
    }
    if (prewarp_frequency > T{0}) {
        const T half_angle = prewarp_frequency * sampling_time / T{2};
        if (!(half_angle < damp::numbers::pi_v<T> / T{2})) {
            return refused<Result>(design::DesignStatus::invalid_parameter); // at or past Nyquist: tan has no finite value
        }
        ts_eff = T{2} * damp::tan(half_angle) / prewarp_frequency;
    }
    const T      ts_half = ts_eff / T{2};
    const Matrix I = Matrix<NX, NX, T>::identity();
    const Matrix L = I - sys.A * ts_half;

    // A_d: solve L · A_d = (I + A·Ts/2)
    const auto A_d_opt = mat::lu_solve(L, I + sys.A * ts_half);
    if (!A_d_opt) {
        return refused<Result>(design::DesignStatus::singular);
    }
    const Matrix A_d = A_d_opt.value();

    // B_d: solve L · X = B, then B_d = Ts · X
    const auto B_d_opt = mat::lu_solve(L, sys.B);
    if (!B_d_opt) {
        return refused<Result>(design::DesignStatus::singular);
    }
    const Matrix B_d = B_d_opt.value() * ts_eff;

    // C_d = C · L⁻¹ → solve Lᵀ · X = Cᵀ, then C_d = Xᵀ
    const auto C_d_t_opt = mat::lu_solve(L.transpose(), sys.C.transpose());
    if (!C_d_t_opt) {
        return refused<Result>(design::DesignStatus::singular);
    }
    const Matrix C_d = C_d_t_opt.value().transpose();

    // D_d = D + (Ts/2) · C_d · B
    const Matrix D_d = sys.D + C_d * sys.B * ts_half;

    // G_d: same map as B
    const auto G_d_opt = mat::lu_solve(L, sys.G);
    if (!G_d_opt) {
        return refused<Result>(design::DesignStatus::singular);
    }
    const Matrix G_d = G_d_opt.value() * ts_eff;
    const Matrix H_d = sys.H;

    return discretized(StateSpace{A_d, B_d, C_d, D_d, G_d, H_d, sampling_time});
}

} // namespace detail

/**
 * @brief Discretize a continuous-time state-space system
 *
 * Converts @f$ \dot x = Ax + Bu @f$ to @f$ x[k+1] = A_d x[k] + B_d u[k] @f$.
 * A model already discrete at @p sampling_time is returned unchanged; resample changes the
 * rate of a discrete model.
 * Tustin warps the frequency axis; @p prewarp_frequency ω₀ > 0 makes the discrete response
 * exact at ω₀ (see detail::discretize_tustin_impl) and is ignored by the other methods.
 *
 * Statuses: DesignStatus::invalid_sample_time when @p sampling_time is not finite and
 * positive, or @p sys is discrete at another sample time (or its Ts is negative or not
 * finite); DesignStatus::bad_weights when A, B or G has a non-finite entry;
 * DesignStatus::singular when Tustin's bilinear factor @f$ I - A T_s/2 @f$ is singular (no
 * silent fallback to ZOH); DesignStatus::invalid_parameter when ω₀ is negative or not finite,
 * or ω₀Ts ≥ π; DesignStatus::not_converged when the discrete model overflows (an unstable
 * mode grows past the range of T within one period).
 *
 * @note Compare with MATLAB®'s c2d(sys, Ts, method) and c2dOptions('PrewarpFrequency', w0).
 * @see "Feedback Control of Dynamic Systems" (Franklin et al., 2015), Chapter 8
 *
 * @param sys               Continuous-time state-space model (@c Ts = 0)
 * @param sampling_time     Desired sampling period [s], finite and positive
 * @param method            Discretization method (ForwardEuler, ZOH, or Tustin)
 * @param prewarp_frequency [rad/s] Tustin pre-warp frequency; 0 (default) for the plain bilinear map
 * @return Discrete system and status; reads like an optional (`*r`, `r->A`, `if (r)`)
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr DiscretizationResult<NX, NU, NY, T, NW, NV> discretize(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        sampling_time,
    DiscretizationMethod                     method = DiscretizationMethod::ZOH,
    T                                        prewarp_frequency = T{0}
) {
    using Result = DiscretizationResult<NX, NU, NY, T, NW, NV>;
    using design::DesignStatus;
    if (!damp::finite_positive(sampling_time)) {
        return detail::refused<Result>(DesignStatus::invalid_sample_time);
    }
    if (sys.Ts != T{0} || !damp::isfinite(sys.Ts)) {
        if (detail::same_sample_time(sys.Ts, sampling_time)) {
            return detail::discretized(sys); // already discrete at this rate
        }
        return detail::refused<Result>(DesignStatus::invalid_sample_time); // another rate: resample, not discretize
    }
    if (!detail::finite_dynamics(sys)) {
        return detail::refused<Result>(DesignStatus::bad_weights);
    }

    const auto checked = [](const Result& r) {
        if (r && !detail::finite_dynamics(*r)) {
            return detail::refused<Result>(DesignStatus::not_converged); // overflowed within one period
        }
        return r;
    };
    switch (method) {
        case DiscretizationMethod::ForwardEuler:
            return checked(detail::discretized(detail::discretize_forward_euler_impl(sys, sampling_time)));
        case DiscretizationMethod::ZOH: {
            const auto d = detail::discretize_zoh_impl(sys, sampling_time);
            if (!d) {
                return detail::refused<Result>(DesignStatus::not_converged); // the exponential overflowed
            }
            return checked(detail::discretized(*d));
        }
        case DiscretizationMethod::Tustin:
            return checked(detail::discretize_tustin_impl(sys, sampling_time, prewarp_frequency));
    }
    return detail::refused<Result>(DesignStatus::invalid_parameter); // not a DiscretizationMethod
}

// ---- Resampling a discrete model -------------------------------------------

/**
 * @brief Resample a ZOH-discrete model to a new sample time
 *
 * A model discretized by ZOH at @f$ T_s @f$ is the exponential of its continuous parent,
 * @f$ M = \exp\!\left(\left[\begin{smallmatrix} A_c & B_c \\ 0 & 0 \end{smallmatrix}\right] T_s\right)
 *   = \left[\begin{smallmatrix} A & B \\ 0 & I \end{smallmatrix}\right] @f$,
 * so the same plant held at @f$ T_s' @f$ is @f$ M^{T_s'/T_s} @f$. An integer ratio (slower
 * rate) is a plain matrix power and always exists. Any other ratio goes through the principal
 * logarithm (the d2c step), which needs A free of eigenvalues on the closed negative real
 * axis: a mode at z ≤ 0 has no continuous parent, and is refused. The noise input G is
 * resampled with B (held, as in discretize); C, D and H do not depend on the rate.
 *
 * @note Compare with MATLAB®'s d2d(sys, Ts) (ZOH method); matlab::d2d is the alias.
 * @see "Computing Integrals Involving the Matrix Exponential" (Van Loan, 1978)
 * @see "Functions of Matrices: Theory and Computation" (Higham, 2008), §11 (logarithm)
 *
 * Statuses: DesignStatus::invalid_sample_time if @p sys is not discrete or the new rate is
 * not finite and positive; DesignStatus::not_converged if a fractional ratio has no real
 * logarithm (the matrix logarithm fails or does not reach a finite real parent).
 *
 * @param sys            Discrete model (Ts > 0), discretized by ZOH
 * @param sampling_time  New sample time [s], finite and positive
 * @return The model at the new rate and status; reads like an optional (`*r`, `r->A`, `if (r)`)
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr DiscretizationResult<NX, NU, NY, T, NW, NV> resample(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    std::type_identity_t<T>                  sampling_time
) {
    using Result = DiscretizationResult<NX, NU, NY, T, NW, NV>;
    if (!finite_positive(sys.Ts) || !finite_positive(sampling_time)) {
        return detail::refused<Result>(design::DesignStatus::invalid_sample_time);
    }
    constexpr size_t N = NX + NU + NW;
    Matrix<N, N, T>  M = Matrix<N, N, T>::identity();
    M.template block<NX, NX>(0, 0) = sys.A;
    M.template block<NX, NU>(0, NX) = sys.B;
    M.template block<NX, NW>(0, NX + NU) = sys.G;

    //! An integer ratio is a power of M; anything else needs log(M).
    const T    ratio = sampling_time / sys.Ts;
    const T    whole = damp::floor(ratio + static_cast<T>(0.5));
    const T    tol = T{16} * eps<T>() * damp::max(T{1}, ratio);
    const bool integer = whole >= T{1} && whole < static_cast<T>(1 << 20) && damp::abs(ratio - whole) <= tol;
    const auto E = integer ? mat::pow(M, static_cast<int>(whole)) : mat::pow(M, ratio);
    if (!E) {
        return detail::refused<Result>(design::DesignStatus::not_converged);
    }
    for (size_t i = 0; i < NX; ++i) {
        for (size_t j = 0; j < N; ++j) {
            if (!damp::isfinite((*E)(i, j))) {
                return detail::refused<Result>(design::DesignStatus::not_converged); // the logarithm did not converge to a real parent
            }
        }
    }
    return detail::discretized(StateSpace<NX, NU, NY, T, NW, NV>{
        .A = detail::leading_rows_block<NX, NX>(*E, 0),
        .B = detail::leading_rows_block<NX, NU>(*E, NX),
        .C = sys.C,
        .D = sys.D,
        .G = detail::leading_rows_block<NX, NW>(*E, NX + NU),
        .H = sys.H,
        .Ts = sampling_time,
    });
}

/**
 * @brief Upsample a discrete model by an integer factor: sys(z) → sys(zᴸ) at Ts / L
 *
 * The result runs @p L times faster and holds the original response between its samples by
 * substituting @f$ z^L @f$ for @f$ z @f$: its impulse response is the original one with
 * L − 1 zeros between samples. The realization chains L copies of the state,
 * @f$ \xi_1^+ = A \xi_L + B u @f$, @f$ \xi_{i+1}^+ = \xi_i @f$, @f$ y = C \xi_L + D u @f$,
 * so the order is L·NX; G enters with B, on the first copy only. Exact, with no matrix function
 * and no failure mode. @p sys must be discrete: a continuous model (Ts = 0) gives a result that
 * models nothing.
 *
 * @note Compare with MATLAB®'s upsample(sys, L); matlab::upsample is the alias.
 *
 * @tparam L  Upsampling factor (≥ 1)
 * @param sys Discrete model (Ts > 0)
 * @return Model of order L·NX at sample time sys.Ts / L
 */
template<size_t L, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
    requires(L >= 1)
[[nodiscard]] constexpr StateSpace<NX * L, NU, NY, T, NW, NV> upsample(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    StateSpace<NX * L, NU, NY, T, NW, NV> out{};
    out.G = Matrix<NX * L, NW, T>{}; // not the square default identity: no noise on the chained copies
    out.A.template block<NX, NX>(0, NX * (L - 1)) = sys.A;
    for (size_t i = 1; i < L; ++i) {
        for (size_t k = 0; k < NX; ++k) {
            out.A((i * NX) + k, ((i - 1) * NX) + k) = T{1};
        }
    }
    out.B.template block<NX, NU>(0, 0) = sys.B;
    out.G.template block<NX, NW>(0, 0) = sys.G;
    out.C.template block<NY, NX>(0, NX * (L - 1)) = sys.C;
    out.D = sys.D;
    out.H = sys.H;
    out.Ts = sys.Ts / static_cast<T>(L);
    return out;
}

// ---- Input delay ------------------------------------------------------------

/**
 * @brief Fractional-delay all-pass filter of order N (Thiran): z⁻ᴰ for a real D = delay / Ts
 *
 * The all-pass
 * @f[
 *   H(z) = \frac{a_N + a_{N-1} z^{-1} + \dots + z^{-N}}{1 + a_1 z^{-1} + \dots + a_N z^{-N}},\qquad
 *   a_k = (-1)^k \binom{N}{k} \prod_{n=0}^{N} \frac{D - N + n}{D - N + k + n}
 * @f]
 * has unit gain at every frequency and a group delay that is maximally flat at DC, where it
 * equals D exactly. It is stable for D > N − 1; the usual choice is N = ⌈D⌉, so the filter
 * spends its N samples on the delay itself. A whole D = N gives the pure delay z⁻ᴺ.
 *
 * The result is a transfer function in z, highest power first as TransferFunction stores it
 * (den = {1, a_1, …, a_N}, MATLAB® thiran's den); realize it with to_state_space and set Ts for a
 * model.
 *
 * @note Compare with MATLAB®'s thiran(tau, Ts), which picks N = ⌈tau/Ts⌉; matlab::thiran is
 *       the alias, with N a template argument.
 * @see "Recursive digital filters with maximally flat group delay" (Thiran, 1971),
 *      IEEE Trans. Circuit Theory 18(6), https://doi.org/10.1109/TCT.1971.1083363
 * @see "Splitting the Unit Delay" (Laakso, Välimäki, Karjalainen & Laine, 1996),
 *      IEEE Signal Processing Magazine 13(1) (all-pass fractional delay)
 *
 * @tparam N              Filter order (≥ 1)
 * @param delay           Delay to approximate [s], ≥ 0
 * @param sampling_time   Sample time [s], finite and positive
 * @return The all-pass H(z), or @c nullopt if D = delay / Ts is not finite or D ≤ N − 1
 *         (an unstable filter)
 */
template<size_t N, typename T>
    requires(N >= 1)
[[nodiscard]] constexpr damp::optional<TransferFunction<N + 1, N + 1, T>> fractional_delay(T delay, std::type_identity_t<T> sampling_time) {
    if (!finite_positive(sampling_time) || !damp::isfinite(delay)) {
        return damp::nullopt;
    }
    const T D = delay / sampling_time;
    const T order = static_cast<T>(N);
    if (!damp::isfinite(D) || !(D > order - T{1})) {
        return damp::nullopt; // D ≤ N − 1: a pole on or outside the unit circle
    }
    damp::array<T, N + 1> a{};
    a[0] = T{1};       // the product is 1 at k = 0 (0/0 term-wise when D is whole)
    T binomial = T{1}; // C(N, k), built up term by term
    for (size_t k = 1; k <= N; ++k) {
        binomial = binomial * static_cast<T>(N - k + 1) / static_cast<T>(k);
        T product = T{1};
        for (size_t n = 0; n <= N; ++n) {
            product *= (D - order + static_cast<T>(n)) / (D - order + static_cast<T>(k + n)); // denominator ≥ D − N + 1 > 0 for k ≥ 1
        }
        a[k] = ((k % 2 == 0) ? binomial : -binomial) * product;
    }
    TransferFunction<N + 1, N + 1, T> H{};
    for (size_t k = 0; k <= N; ++k) {
        H.num[k] = a[N - k]; // z^N · (a_N + … + z⁻ᴺ): coefficient of z^(N−k) is a_(N−k)
        H.den[k] = a[k];     // z^N · (1 + a_1 z⁻¹ + … + a_N z⁻ᴺ): coefficient of z^(N−k) is a_k
    }
    return H;
}

/**
 * @brief An input delay of θ = @ref seconds, in up to N registers (0 ≤ θ ≤ N·Ts): each register
 *        state holds an input already sent, v(k−1) … v(k−N). `discretize(sys, Ts, InputDelay<N>{theta})`
 *        returns that model; design on it like any other. `InputDelay{theta}` is N = 1.
 *
 * N sizes the model and θ is a value, so the model's type is the same for any delay up to N·Ts,
 * one measured at init included. Each register is a state of every matrix the law runs on, at a
 * cost of (NX + N·NU)² per tick for a filter: for many samples of dead time, sample slower or,
 * for a single loop, use SmithPredictor.
 *
 * @tparam N Registers, the longest delay the model can hold in samples
 * @tparam T Scalar
 */
template<size_t N = 1, typename T = double>
struct InputDelay {
    static_assert(N >= 1, "InputDelay holds at least one register: a plant without delay is discretize(sys, Ts).");
    T seconds{T{0}}; ///< [s] delay between a command and the plant receiving it, 0 ≤ θ ≤ N·Ts
};

/** @brief One register, the scalar from the seconds: `InputDelay{theta}`. */
template<typename T>
InputDelay(T) -> InputDelay<1, T>;

namespace detail {

/** @brief A delay written as d·Ts + τ. */
template<typename T>
struct RegisterDelay {
    size_t whole{0};       ///< d, the whole periods
    T      fraction{T{0}}; ///< [s] τ, 0 ≤ τ < Ts: when the held input switches within the period
};

/**
 * @brief A delay of @p seconds as whole periods and a remainder, the remainder snapped to 0 or Ts
 *        within 16·eps·max(θ, Ts).
 * @return The split, or nullopt if the delay is negative, not finite or longer than @p capacity
 *         samples, or Ts is not finite and positive.
 */
template<typename T>
[[nodiscard]] constexpr damp::optional<RegisterDelay<T>> register_delay(T seconds, T sampling_time, size_t capacity) {
    if (!damp::isfinite(seconds) || !(seconds >= T{0}) || !damp::finite_positive(sampling_time)) {
        return damp::nullopt;
    }
    const T tol = T{16} * eps<T>() * damp::max(seconds, sampling_time);
    if (seconds > (static_cast<T>(capacity) * sampling_time) + tol) {
        return damp::nullopt;
    }
    RegisterDelay<T> out{static_cast<size_t>(seconds / sampling_time), T{0}}; // at most capacity: checked above
    const T          rest = seconds - (static_cast<T>(out.whole) * sampling_time);
    if (rest > sampling_time - tol) {
        ++out.whole;
    } else if (rest > tol && out.whole < capacity) {
        out.fraction = rest;
    }
    return out;
}

/**
 * @brief The plant state followed by N registers holding the inputs of the last N periods.
 *
 * State [x; v(k−1); …; v(k−N)], input v(k): x(k+1) = Φx + Γ₀v(k−d) + Γ₁v(k−d−1), the first
 * register loads v(k) and each of the others the one before it. The output feedthrough @p Dy
 * reads the input held at the sample instant: v(k−d−1) when @p fraction (the switch falls inside
 * the period), v(k−d) otherwise. The registers are commands already sent: no process noise.
 */
template<size_t N, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr StateSpace<NX + (N * NU), NU, NY, T, NW, NV> input_registers(
    const StateSpace<NX, NU, NY, T, NW, NV>& base, // A, C, G, H, Ts of the plant
    const Matrix<NX, NU, T>&                 Gamma0,
    const Matrix<NX, NU, T>&                 Gamma1,
    const Matrix<NY, NU, T>&                 Dy,
    size_t                                   whole, // d
    bool                                     fraction
) {
    constexpr size_t                  NS = NX + (N * NU);
    StateSpace<NS, NU, NY, T, NW, NV> out{};
    const auto                        Iu = Matrix<NU, NU, T>::identity();
    out.G = Matrix<NS, NW, T>{}; // not the square default identity: no noise on the registers
    out.A.template block<NX, NX>(0, 0) = base.A;
    out.B.template block<NU, NU>(NX, 0) = Iu;
    for (size_t j = 1; j < N; ++j) {
        out.A.template block<NU, NU>(NX + (j * NU), NX + ((j - 1) * NU)) = Iu;
    }
    if (whole == 0) {
        out.B.template block<NX, NU>(0, 0) = Gamma0;
    } else {
        out.A.template block<NX, NU>(0, NX + ((whole - 1) * NU)) = Gamma0;
    }
    if (fraction) {
        out.A.template block<NX, NU>(0, NX + (whole * NU)) = Gamma1;
    }
    out.C.template block<NY, NX>(0, 0) = base.C;
    const size_t lag = fraction ? whole + 1 : whole;
    if (lag == 0) {
        out.D = Dy;
    } else {
        out.C.template block<NY, NU>(0, NX + ((lag - 1) * NU)) = Dy;
    }
    out.G.template block<NX, NW>(0, 0) = base.G;
    out.H = base.H;
    out.Ts = base.Ts;
    return out;
}

} // namespace detail

/**
 * @brief Discretize with an input delay of up to N samples: NX + N·NU states, whatever the delay.
 *
 * The model carries the delay θ = d·Ts + τ (0 ≤ τ < Ts) in N registers: within each period the
 * plant holds v(k−d−1) until τ, then v(k−d), so
 * @f[
 *   x_{k+1} = \Phi x_k + \Gamma_0 v_{k-d} + \Gamma_1 v_{k-d-1},\qquad
 *   \Gamma_0 = \int_0^{T_s-\tau} e^{A\eta}\,d\eta\,B,\quad
 *   \Gamma_1 = e^{A(T_s-\tau)} \int_0^{\tau} e^{A\eta}\,d\eta\,B,
 * @f]
 * exact for a held input (no Padé approximation): @f$ \Gamma_1 = 0 @f$ and
 * @f$ \Gamma_0 = \Gamma @f$ for a whole number of samples. A fraction of a sample needs ZOH and
 * a continuous plant; a whole number works with any @p method and on a model that is already
 * discrete (same Ts). A shorter delay than N samples leaves the last registers unread.
 *
 * G is not delayed: a physical disturbance in G of the continuous plant (a load force,
 * friction) gets its exact ZOH entry on the plant states, the B_d that design::PlantDisturbance
 * takes from it.
 *
 * @code
 * constexpr auto model = discretize(plant, 2.0, InputDelay{1.4});     // NX + NU states
 * constexpr auto longer = discretize(plant, 2.0, InputDelay<3>{5.4}); // NX + 3·NU states
 * @endcode
 * The model has the same type for any delay up to N samples; design on it like any plant:
 * `design::mpc<NP, NC>(model, weights, constraints)`, with NP long enough for a command to land
 * (design::mpc refuses one too short). A state weight covers the registers, none on them:
 * `{.Q = blkdiag(Q, Matrix<N * NU, N * NU>{}), .R = R}`. A law that measures its state reads
 * [x; u(k−1); …; u(k−N)], the inputs it applied; a filter estimates the registers like any other
 * state. design::luenberger needs every register read, and A pins NU of the register poles at 0,
 * one per input. It places a pole at most NY times, so with more registers than that put the
 * other register poles near 0 but apart (1e-3, 2e-3, …), or use design::kalman.
 *
 * @note Compare with MATLAB®'s c2d on a model with InputDelay, then absorbDelay.
 * @see "Digital Control of Dynamic Systems" (Franklin, Powell, Workman, 1998), §4.3.3
 *
 * Statuses: DesignStatus::invalid_sample_time if Ts is not finite and positive, a discrete
 * model's Ts differs, or a fraction of a sample is asked of a discrete model;
 * DesignStatus::invalid_parameter if the delay is negative, not finite or longer than N samples,
 * or a fraction of a sample is asked of a non-ZOH method; otherwise the status of the plain
 * discretize step (run at Ts for a fraction of a sample too).
 *
 * @return The model and status; reads like an optional (`*r`, `r->A`, `if (r)`)
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0, size_t N = 1, typename TD = T>
[[nodiscard]] constexpr DiscretizationResult<NX + (N * NU), NU, NY, T, NW, NV> discretize(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        sampling_time,
    InputDelay<N, TD>                        delay, // any scalar: InputDelay<N>{θ} is double, on a float plant too
    DiscretizationMethod                     method = DiscretizationMethod::ZOH
) {
    using Result = DiscretizationResult<NX + (N * NU), NU, NY, T, NW, NV>;
    using design::DesignStatus;
    if (!damp::finite_positive(sampling_time)) {
        return detail::refused<Result>(DesignStatus::invalid_sample_time);
    }
    const auto split = detail::register_delay(static_cast<T>(delay.seconds), sampling_time, N);
    if (!split) {
        return detail::refused<Result>(DesignStatus::invalid_parameter);
    }
    const T tau = split->fraction;
    if (tau == T{0}) { // a whole number of samples: Γ₀ = Γ, Γ₁ = 0
        const auto base = discretize(sys, sampling_time, method);
        if (!base) {
            return detail::refused<Result>(base.status);
        }
        return detail::discretized(detail::input_registers<N>(*base, base->B, Matrix<NX, NU, T>{}, base->D, split->whole, false));
    }
    if (sys.Ts != T{0}) {
        return detail::refused<Result>(DesignStatus::invalid_sample_time); // a fraction of a sample needs the continuous plant
    }
    if (method != DiscretizationMethod::ZOH) {
        return detail::refused<Result>(DesignStatus::invalid_parameter); // a fraction of a sample needs ZOH
    }
    const auto full = discretize(sys, sampling_time, DiscretizationMethod::ZOH); // checks A, B, G and the result
    if (!full) {
        return detail::refused<Result>(full.status);
    }
    const auto head = detail::discretize_zoh_impl(sys, tau);                 // ∫₀^τ e^{Aη}dη B
    const auto tail = detail::discretize_zoh_impl(sys, sampling_time - tau); // e^{A(Ts−τ)}, ∫₀^{Ts−τ} e^{Aη}dη B
    if (!head || !tail) {
        return detail::refused<Result>(DesignStatus::not_converged); // the exponential overflowed
    }
    const Matrix<NX, NU, T> gamma1 = tail->A * head->B;
    return detail::discretized(detail::input_registers<N>(*full, tail->B, gamma1, sys.D, split->whole, true));
}

} // namespace damp
