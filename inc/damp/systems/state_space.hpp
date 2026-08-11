// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file state_space.hpp
 * @brief State-space representation for linear time-invariant systems
 *
 * Shared LTI type: continuous (@c Ts == 0) or discrete (@c Ts > 0), with
 * optional process/measurement noise channels @c G/@c H preserved through
 * interconnections.
 *
 * Interconnection operators (MATLAB®-compatible):
 * | damp | MATLAB® | Signal path |
 * | --- | ------- | ----------- |
 * | @c series(sys1, sys2) | @c series(sys1,sys2) | @f$ u \to \mathrm{sys1} \to \mathrm{sys2} \to y @f$ |
 * | @c sys2 * sys1 | @c sys2*sys1 | same as @c series(sys1,sys2) |
 * | @c parallel / @c + | @c parallel / @c + | shared @f$u@f$, sum outputs |
 * | @c feedback / @c / | @c feedback | negative feedback @f$ y = \mathrm{sys1}(r - \mathrm{sys2}(y)) @f$ |
 *
 * Mismatched sample times or a singular algebraic loop return
 * @c damp::nullopt (MATLAB® errors in those cases).
 *
 * @code
 * using namespace damp;
 * constexpr StateSpace<2, 1, 1> sys{
 *     .A = Matrix<2,2>{{0.0, 1.0}, {0.0, -0.1}},
 *     .B = Matrix<2,1>{{0.0}, {1.0}},
 *     .C = Matrix<1,2>{{1.0, 0.0}},
 *     .Ts = 0.0, // continuous
 * };
 * static_assert(sys.is_continuous());
 * @endcode
 *
 * @see discretize() in systems/discretization.hpp
 * @see TransferFunction::to_state_space()
 * @see "Feedback Control of Dynamic Systems" (Franklin, Powell & Emami-Naeini)
 * @note Compare with MATLAB®'s ss, series, parallel, feedback, and * / + operators.
 */

#include <cstddef>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/matrix_traits.hpp"
#include "damp/matrix/solve.hpp"

namespace damp {

namespace detail {

/**
 * @brief Default noise matrices: identity when square, zero otherwise
 *
 * When @p NW==@p NX (resp. @p NV==@p NY), process (measurement) noise enters
 * every state (output) channel one-to-one — the usual Kalman-shaped default.
 * Non-square noise sizes start at zero and must be filled by the user.
 */
template<size_t R, size_t C, typename T>
constexpr Matrix<R, C, T> default_noise_matrix() {
    if constexpr (R == C && R > 0) {
        return Matrix<R, C, T>::identity();
    } else {
        return Matrix<R, C, T>{};
    }
}

/**
 * @brief True if both systems are continuous or share the same sample period
 *
 * Matches MATLAB®: series/parallel/feedback require compatible time domains.
 */
template<typename T>
[[nodiscard]] constexpr bool same_sample_time(T ts1, T ts2) {
    if (ts1 == T{0} && ts2 == T{0}) {
        return true;
    }
    if (ts1 > T{0} && ts2 > T{0}) {
        const T scale = damp::max(T{1}, damp::max(ts1, ts2));
        return damp::abs(ts1 - ts2) <= default_tol<T>() * scale;
    }
    return false;
}

} // namespace detail

/**
 * @brief State-space representation for linear time-invariant systems (discrete or continuous)
 *
 * Fixed-size, stack-allocated LTI container. @c Ts == 0 is continuous; @c Ts > 0 is discrete.
 * Process/measurement noise channels @c G (NX×NW) and @c H (NY×NV) are preserved through
 * interconnection operators (* series, + parallel, − differencing, / feedback) when the
 * composition succeeds.
 *
 * Discrete-time (@c Ts > 0):
 * @f[
 *   x_{k+1} = A x_k + B u_k + G w_k,\qquad
 *   y_k = C x_k + D u_k + H v_k
 * @f]
 *
 * Continuous-time (@c Ts = 0):
 * @f[
 *   \dot x = A x + B u + G w,\qquad
 *   y = C x + D u + H v
 * @f]
 *
 * @see discretize() in systems/discretization.hpp for continuous→discrete maps
 * @see TransferFunction::to_state_space() for companion-form realization
 *
 * @tparam NX Number of states (rows/cols of A)
 * @tparam NU Number of control inputs (cols of B and D)
 * @tparam NY Number of outputs (rows of C and D)
 * @tparam T  Scalar type (default double); requires floating-point
 * @tparam NW Number of process noise inputs (cols of G); default 0
 * @tparam NV Number of measurement noise inputs (cols of H); default 0
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
    requires std::is_floating_point_v<T>
struct StateSpace {
    Matrix<NX, NX, T> A{}; ///< State dynamics matrix
    Matrix<NX, NU, T> B{}; ///< Control input matrix
    Matrix<NY, NX, T> C{}; ///< Output matrix
    Matrix<NY, NU, T> D{}; ///< Direct feedthrough matrix

    Matrix<NX, NW, T> G = detail::default_noise_matrix<NX, NW, T>(); ///< Process noise input matrix
    Matrix<NY, NV, T> H = detail::default_noise_matrix<NY, NV, T>(); ///< Measurement noise input matrix

    T Ts = T{0}; ///< Sampling period (0 for continuous, > 0 for discrete)

    template<typename U>
    [[nodiscard]] constexpr StateSpace<NX, NU, NY, U, NW, NV> as() const {
        return StateSpace<NX, NU, NY, U, NW, NV>{
            A.template as<U>(),
            B.template as<U>(),
            C.template as<U>(),
            D.template as<U>(),
            G.template as<U>(),
            H.template as<U>(),
            static_cast<U>(Ts)
        };
    }

    [[nodiscard]] constexpr bool is_discrete() const { return Ts > T{0}; }
    [[nodiscard]] constexpr bool is_continuous() const { return Ts == T{0}; }
};

/**
 * @brief Evaluate frequency response of a state-space system
 *
 * Continuous: @f$ G(s) = C(sI-A)^{-1}B + D @f$.
 * Discrete: @f$ G(z) = C(zI-A)^{-1}B + D @f$.
 *
 * Uses mat::solve (no explicit inverse). Returns @c nullopt if
 * @f$ sI-A @f$ (or @f$ zI-A @f$) is singular — e.g. evaluating exactly at a pole.
 *
 * @note Compare with MATLAB®'s evalfr / freqresp at a single complex point.
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr damp::optional<Matrix<NY, NU, damp::complex<T>>>
eval_frf(const StateSpace<NX, NU, NY, T, NW, NV>& sys, damp::complex<T> s) {
    using Cplx = damp::complex<T>;
    auto I = Matrix<NX, NX, Cplx>::identity();
    auto sI_minus_A = s * I - sys.A.template as<Cplx>();
    auto X_opt = mat::solve(sI_minus_A, sys.B.template as<Cplx>());
    if (!X_opt) {
        return damp::nullopt;
    }
    auto temp = sys.C.template as<Cplx>() * (*X_opt);
    return temp + sys.D.template as<Cplx>();
}

/**
 * @brief Series connection: @p sys2 follows @p sys1 (u → sys1 → sys2 → y)
 *
 * @f[
 *   A = \begin{bmatrix} A_1 & 0 \\ B_2 C_1 & A_2 \end{bmatrix},\;
 *   B = \begin{bmatrix} B_1 \\ B_2 D_1 \end{bmatrix},\;
 *   C = \begin{bmatrix} D_2 C_1 & C_2 \end{bmatrix},\;
 *   D = D_2 D_1
 * @f]
 *
 * @return Composed system, or @c nullopt if sample times differ
 * @note Compare with MATLAB®'s series(sys1,sys2) (≡ @c sys2*sys1).
 */
template<
    size_t NX1, size_t NU1, size_t NY1, size_t NW1, size_t NV1,
    size_t NX2, size_t NU2, size_t NY2, size_t NW2, size_t NV2,
    typename T>
    requires(NY1 == NU2)
[[nodiscard]] constexpr auto series(
    const StateSpace<NX1, NU1, NY1, T, NW1, NV1>& sys1,
    const StateSpace<NX2, NU2, NY2, T, NW2, NV2>& sys2
) -> damp::optional<StateSpace<NX1 + NX2, NU1, NY2, T, NW1 + NW2, NV1 + NV2>> {
    if (!detail::same_sample_time(sys1.Ts, sys2.Ts)) {
        return damp::nullopt;
    }

    constexpr size_t n1 = NX1;
    constexpr size_t n2 = NX2;
    constexpr size_t m = NU1;
    constexpr size_t p = NY2;
    constexpr size_t nw1 = NW1;
    constexpr size_t nw2 = NW2;
    constexpr size_t nv1 = NV1;
    constexpr size_t nv2 = NV2;

    Matrix<n1 + n2, n1 + n2, T>   A{};
    Matrix<n1 + n2, m, T>         B{};
    Matrix<p, n1 + n2, T>         C{};
    Matrix<p, m, T>               D{};
    Matrix<n1 + n2, nw1 + nw2, T> G{};
    Matrix<p, nv1 + nv2, T>       H{};

    A.template block<n1, n1>(0, 0) = sys1.A;
    A.template block<n2, n1>(n1, 0) = sys2.B * sys1.C;
    A.template block<n2, n2>(n1, n1) = sys2.A;

    B.template block<n1, m>(0, 0) = sys1.B;
    B.template block<n2, m>(n1, 0) = sys2.B * sys1.D;

    C.template block<p, n1>(0, 0) = sys2.D * sys1.C;
    C.template block<p, n2>(0, n1) = sys2.C;

    D = sys2.D * sys1.D;

    G.template block<n1, nw1>(0, 0) = sys1.G;
    G.template block<n2, nw1>(n1, 0) = sys2.B * sys1.H;
    G.template block<n2, nw2>(n1, nw1) = sys2.G;

    H.template block<p, nv1>(0, 0) = sys2.D * sys1.H;
    H.template block<p, nv2>(0, nv1) = sys2.H;

    return StateSpace<n1 + n2, m, p, T, nw1 + nw2, nv1 + nv2>{A, B, C, D, G, H, sys1.Ts};
}

/**
 * @brief Parallel connection (shared input, summed outputs)
 *
 * @return Composed system, or @c nullopt if sample times differ
 * @note Compare with MATLAB®'s parallel(sys1,sys2) / @c sys1+sys2.
 */
template<
    size_t NX1, size_t NU1, size_t NY1, size_t NW1, size_t NV1,
    size_t NX2, size_t NU2, size_t NY2, size_t NW2, size_t NV2,
    typename T>
    requires(NU1 == NU2) && (NY1 == NY2)
[[nodiscard]] constexpr auto parallel(
    const StateSpace<NX1, NU1, NY1, T, NW1, NV1>& sys1,
    const StateSpace<NX2, NU2, NY2, T, NW2, NV2>& sys2
) -> damp::optional<StateSpace<NX1 + NX2, NU1, NY1, T, NW1 + NW2, NV1 + NV2>> {
    if (!detail::same_sample_time(sys1.Ts, sys2.Ts)) {
        return damp::nullopt;
    }

    constexpr size_t n1 = NX1;
    constexpr size_t n2 = NX2;
    constexpr size_t m = NU1;
    constexpr size_t p = NY1;
    constexpr size_t nw1 = NW1;
    constexpr size_t nw2 = NW2;
    constexpr size_t nv1 = NV1;
    constexpr size_t nv2 = NV2;

    Matrix<n1 + n2, n1 + n2, T>   A{};
    Matrix<n1 + n2, m, T>         B{};
    Matrix<p, n1 + n2, T>         C{};
    Matrix<p, m, T>               D{};
    Matrix<n1 + n2, nw1 + nw2, T> G{};
    Matrix<p, nv1 + nv2, T>       H{};

    A.template block<n1, n1>(0, 0) = sys1.A;
    A.template block<n2, n2>(n1, n1) = sys2.A;

    B.template block<n1, m>(0, 0) = sys1.B;
    B.template block<n2, m>(n1, 0) = sys2.B;

    C.template block<p, n1>(0, 0) = sys1.C;
    C.template block<p, n2>(0, n1) = sys2.C;

    D = sys1.D + sys2.D;

    G.template block<n1, nw1>(0, 0) = sys1.G;
    G.template block<n2, nw2>(n1, nw1) = sys2.G;

    H.template block<p, nv1>(0, 0) = sys1.H;
    H.template block<p, nv2>(0, nv1) = sys2.H;

    return StateSpace<n1 + n2, m, p, T, nw1 + nw2, nv1 + nv2>{A, B, C, D, G, H, sys1.Ts};
}

/**
 * @brief Negative feedback: y = sys1(u − sys2(y))
 *
 * With @f$ V=(I+D_2 D_1)^{-1} @f$, @f$ W=(I+D_1 D_2)^{-1} @f$ (solved via LU, not explicit inverse).
 * Returns @c nullopt if sample times differ or either matrix is singular (ill-posed algebraic loop).
 *
 * @note Process/measurement noise (G, H) use the D=0 structural placement through the loop.
 * @note Compare with MATLAB®'s feedback(sys1,sys2) (negative feedback).
 */
template<
    size_t NX1, size_t NU1, size_t NY1, size_t NW1, size_t NV1,
    size_t NX2, size_t NU2, size_t NY2, size_t NW2, size_t NV2,
    typename T>
    requires(NY1 == NU2) && (NU1 == NY2)
[[nodiscard]] constexpr auto feedback(
    const StateSpace<NX1, NU1, NY1, T, NW1, NV1>& sys1,
    const StateSpace<NX2, NU2, NY2, T, NW2, NV2>& sys2
) -> damp::optional<StateSpace<NX1 + NX2, NU1, NY1, T, NW1 + NW2, NV1 + NV2>> {
    if (!detail::same_sample_time(sys1.Ts, sys2.Ts)) {
        return damp::nullopt;
    }

    constexpr size_t n1 = NX1;
    constexpr size_t n2 = NX2;
    constexpr size_t m = NU1;
    constexpr size_t p = NY1;
    constexpr size_t nw1 = NW1;
    constexpr size_t nw2 = NW2;
    constexpr size_t nv1 = NV1;
    constexpr size_t nv2 = NV2;

    Matrix<n1 + n2, n1 + n2, T>   A{};
    Matrix<n1 + n2, m, T>         B{};
    Matrix<p, n1 + n2, T>         C{};
    Matrix<p, m, T>               D{};
    Matrix<n1 + n2, nw1 + nw2, T> G{};
    Matrix<p, nv1 + nv2, T>       H{};

    const Matrix<m, m, T> I_m = Matrix<m, m, T>::identity();
    const Matrix<p, p, T> I_p = Matrix<p, p, T>::identity();

    // V = (I_m + D2 D1)^{-1}, W = (I_p + D1 D2)^{-1} via multi-RHS LU solve
    const auto V_opt = mat::lu_solve(I_m + sys2.D * sys1.D, I_m);
    const auto W_opt = mat::lu_solve(I_p + sys1.D * sys2.D, I_p);
    if (!V_opt || !W_opt) {
        return damp::nullopt;
    }
    const Matrix<m, m, T> V = *V_opt;
    const Matrix<p, p, T> W = *W_opt;

    const Matrix B1_V = sys1.B * V;
    const Matrix W_C1 = W * sys1.C;
    const Matrix W_D1 = W * sys1.D;

    A.template block<n1, n1>(0, 0) = sys1.A - B1_V * sys2.D * sys1.C;
    if constexpr (n2 > 0) {
        A.template block<n1, n2>(0, n1) = -B1_V * sys2.C;
        A.template block<n2, n1>(n1, 0) = sys2.B * W_C1;
        A.template block<n2, n2>(n1, n1) = sys2.A - sys2.B * W_D1 * sys2.C;
    }

    B.template block<n1, m>(0, 0) = B1_V;
    if constexpr (n2 > 0) {
        B.template block<n2, m>(n1, 0) = sys2.B * W_D1;
    }

    C.template block<p, n1>(0, 0) = W_C1;
    if constexpr (n2 > 0) {
        C.template block<p, n2>(0, n1) = -W_D1 * sys2.C;
    }

    D = W_D1;

    G.template block<n1, nw1>(0, 0) = sys1.G;
    if constexpr (n2 > 0) {
        G.template block<n2, nw1>(n1, 0) = sys2.B * sys1.H;
        G.template block<n2, nw2>(n1, nw1) = sys2.G;
    }

    H.template block<p, nv1>(0, 0) = sys1.H;
    if constexpr (nv2 > 0) {
        H.template block<p, nv2>(0, nv1) = sys2.H;
    }

    return StateSpace<n1 + n2, m, p, T, nw1 + nw2, nv1 + nv2>{A, B, C, D, G, H, sys1.Ts};
}

/**
 * @brief Differencing connection: outputs y = y₁ − y₂
 *
 * @return Composed system, or @c nullopt if sample times differ
 */
template<
    size_t NX1, size_t NU1, size_t NY1, size_t NW1, size_t NV1,
    size_t NX2, size_t NU2, size_t NY2, size_t NW2, size_t NV2,
    typename T>
    requires(NU1 == NU2) && (NY1 == NY2)
[[nodiscard]] constexpr auto subtract(
    const StateSpace<NX1, NU1, NY1, T, NW1, NV1>& sys1,
    const StateSpace<NX2, NU2, NY2, T, NW2, NV2>& sys2
) -> damp::optional<StateSpace<NX1 + NX2, NU1, NY1, T, NW1 + NW2, NV1 + NV2>> {
    if (!detail::same_sample_time(sys1.Ts, sys2.Ts)) {
        return damp::nullopt;
    }

    constexpr size_t n1 = NX1;
    constexpr size_t n2 = NX2;
    constexpr size_t m = NU1;
    constexpr size_t p = NY1;
    constexpr size_t nw1 = NW1;
    constexpr size_t nw2 = NW2;
    constexpr size_t nv1 = NV1;
    constexpr size_t nv2 = NV2;

    Matrix<n1 + n2, n1 + n2, T>   A{};
    Matrix<n1 + n2, m, T>         B{};
    Matrix<p, n1 + n2, T>         C{};
    Matrix<p, m, T>               D{};
    Matrix<n1 + n2, nw1 + nw2, T> G{};
    Matrix<p, nv1 + nv2, T>       H{};

    A.template block<n1, n1>(0, 0) = sys1.A;
    A.template block<n2, n2>(n1, n1) = sys2.A;

    B.template block<n1, m>(0, 0) = sys1.B;
    B.template block<n2, m>(n1, 0) = sys2.B;

    C.template block<p, n1>(0, 0) = sys1.C;
    C.template block<p, n2>(0, n1) = -sys2.C;

    D = sys1.D - sys2.D;

    G.template block<n1, nw1>(0, 0) = sys1.G;
    G.template block<n2, nw2>(n1, nw1) = sys2.G;

    H.template block<p, nv1>(0, 0) = sys1.H;
    H.template block<p, nv2>(0, nv1) = sys2.H;

    return StateSpace<n1 + n2, m, p, T, nw1 + nw2, nv1 + nv2>{A, B, C, D, G, H, sys1.Ts};
}

// --- Operators (MATLAB®-compatible; return optional like the free functions) ---

template<
    size_t NX1, size_t NU1, size_t NY1, size_t NW1, size_t NV1,
    size_t NX2, size_t NU2, size_t NY2, size_t NW2, size_t NV2,
    typename T>
[[nodiscard]] constexpr auto operator+(
    const StateSpace<NX1, NU1, NY1, T, NW1, NV1>& sys1,
    const StateSpace<NX2, NU2, NY2, T, NW2, NV2>& sys2
) {
    return parallel(sys1, sys2);
}

template<
    size_t NX1, size_t NU1, size_t NY1, size_t NW1, size_t NV1,
    size_t NX2, size_t NU2, size_t NY2, size_t NW2, size_t NV2,
    typename T>
[[nodiscard]] constexpr auto operator-(
    const StateSpace<NX1, NU1, NY1, T, NW1, NV1>& sys1,
    const StateSpace<NX2, NU2, NY2, T, NW2, NV2>& sys2
) {
    return subtract(sys1, sys2);
}

/**
 * @brief Series product matching MATLAB®: @c sys1 * sys2 ≡ @c series(sys2, sys1)
 *
 * Signal path: @f$ u \to \mathrm{sys2} \to \mathrm{sys1} \to y @f$
 * (right factor acts first). Prefer @ref series when reading left-to-right.
 */
template<
    size_t NX1, size_t NU1, size_t NY1, size_t NW1, size_t NV1,
    size_t NX2, size_t NU2, size_t NY2, size_t NW2, size_t NV2,
    typename T>
[[nodiscard]] constexpr auto operator*(
    const StateSpace<NX1, NU1, NY1, T, NW1, NV1>& sys1,
    const StateSpace<NX2, NU2, NY2, T, NW2, NV2>& sys2
) {
    return series(sys2, sys1);
}

template<
    size_t NX1, size_t NU1, size_t NY1, size_t NW1, size_t NV1,
    size_t NX2, size_t NU2, size_t NY2, size_t NW2, size_t NV2,
    typename T>
[[nodiscard]] constexpr auto operator/(
    const StateSpace<NX1, NU1, NY1, T, NW1, NV1>& sys1,
    const StateSpace<NX2, NU2, NY2, T, NW2, NV2>& sys2
) {
    return feedback(sys1, sys2);
}

} // namespace damp
