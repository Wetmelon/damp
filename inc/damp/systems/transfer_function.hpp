// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file transfer_function.hpp
 * @brief SISO polynomial transfer function G(s) or G(z) = num/den
 *
 * Coefficients are highest power first, as MATLAB® writes them and as a characteristic
 * equation reads: num = {b_m, …, b_1, b_0} is b_m s^m + … + b_1 s + b_0.
 * There is no @c Ts field: continuous vs discrete is the caller's convention
 * (same as ZPK). Realize with to_state_space then discretize for c2d.
 *
 * Algebra (@c *, @c +, @c feedback) expands polynomials without cancellation —
 * use design::minreal or ZPK @c cancel_matching when degree growth matters.
 *
 * @note Compare with MATLAB®'s tf.
 * @see TransferFunction::to_state_space, systems/state_space.hpp, systems/zpk.hpp
 */


#include "damp/backend.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix_traits.hpp"
#include "state_space.hpp"

namespace damp {

namespace tf_detail {

/** @brief dst += scale·src, aligned on the constant term (the last entry: coefficients are highest power first). */
template<size_t Nres, size_t Nsrc, typename T>
constexpr void accumulate_poly(damp::array<T, Nres>& dst, const damp::array<T, Nsrc>& src, T scale = T{1}) {
    constexpr size_t N = (Nsrc < Nres) ? Nsrc : Nres;
    for (size_t i = 0; i < N; ++i) {
        dst[Nres - 1 - i] += scale * src[Nsrc - 1 - i];
    }
}

template<size_t Na, size_t Nb, typename T>
[[nodiscard]] constexpr damp::array<T, Na + Nb - 1> convolve_poly(const damp::array<T, Na>& a, const damp::array<T, Nb>& b) {
    damp::array<T, Na + Nb - 1> result{};
    for (size_t i = 0; i < Na; ++i) {
        for (size_t j = 0; j < Nb; ++j) {
            result[i + j] += a[i] * b[j];
        }
    }
    return result;
}

/**
 * @brief True if every root of z^(N−1)·p(z⁻¹) lies strictly inside the unit circle (Schur–Cohn).
 *
 * Step-down recursion on the reflection coefficients kₙ = aₙ/a₀: the polynomial is Schur-stable
 * iff every |kₙ| < 1. A root on the unit circle (|kₙ| = 1 to rounding) is refused. kₙ is
 * dimensionless, so the margin N·ε of T does not depend on the coefficients' units.
 *
 * Float coefficients are judged in double, where they are exact: the recursion's own rounding
 * in float can pass a polynomial with a root outside the circle (a fourth-order cluster near
 * z = −1, depending on whether the compiler fuses its multiply-adds).
 */
template<size_t N, typename T>
[[nodiscard]] constexpr bool schur_stable(const damp::array<T, N>& p) {
    using W = damp::conditional_t<damp::is_same_v<T, float>, double, T>;
    using R = scalar_type_t<W>;
    damp::array<W, N> a{};
    for (size_t i = 0; i < N; ++i) {
        a[i] = static_cast<W>(p[i]);
    }
    const R margin = static_cast<R>(N) * static_cast<R>(eps<T>());
    for (size_t n = N; n-- > 1;) {
        const W k = a[n] / a[0];
        if (!(damp::abs(k) < R{1} - margin)) {
            return false;
        }
        const damp::array<W, N> prev = a;
        for (size_t i = 0; i < n; ++i) {
            a[i] = prev[i] - (k * prev[n - i]);
        }
    }
    return true;
}

} // namespace tf_detail

template<typename TNum, typename TDen>
using transfer_function_scalar_t = damp::conditional_t<
    damp::is_floating_point_v<damp::common_type_t<TNum, TDen>>,
    damp::common_type_t<TNum, TDen>,
    double>;

/**
 * @brief SISO polynomial transfer function G(s) = num(s)/den(s)
 *
 * Coefficients are highest power first, as in MATLAB®'s tf: den[0] multiplies the highest
 * power, den[Nden-1] is the constant term. Same convention for num, so
 * `{.num = {1, 2}, .den = {1, 3, 2}}` is (s + 2)/(s² + 3s + 2).
 *
 * @tparam Nnum Number of numerator coefficients
 * @tparam Nden Number of denominator coefficients
 * @tparam T    Scalar type (default double)
 */
template<size_t Nnum, size_t Nden, typename T = double>
    requires damp::is_floating_point_v<T>
struct TransferFunction {
    damp::array<T, Nnum> num{}; ///< Numerator coefficients, highest power of s first
    damp::array<T, Nden> den{}; ///< Denominator coefficients, highest power of s first

    template<typename U>
    [[nodiscard]] constexpr TransferFunction<Nnum, Nden, U> as() const {
        TransferFunction<Nnum, Nden, U> result{};
        for (size_t i = 0; i < Nnum; ++i) {
            result.num[i] = static_cast<U>(num[i]);
        }
        for (size_t i = 0; i < Nden; ++i) {
            result.den[i] = static_cast<U>(den[i]);
        }
        return result;
    }

    /**
     * @brief Controllable companion-form state-space realization
     *
     * Realizes @f$ G = \mathrm{num}/\mathrm{den} @f$ as SISO StateSpace with
     * @c Nden−1 states. Requires a proper TF (@c Nnum ≤ Nden) and nonzero leading
     * denominator coefficient. Biproper systems set @f$ D = \mathrm{num}[0]/\mathrm{den}[0] @f$ and reduce C.
     *
     * Only an exactly zero leading coefficient is rejected: polynomial
     * coefficients scale as powers of the frequency unit, so no magnitude
     * threshold on them is unit-invariant.
     *
     * @return Continuous realization (@c Ts = 0), or @c nullopt if
     *         @f$ \mathrm{den}[0] = 0 @f$ (padded / improper array)
     * @note Compare with MATLAB®'s ss(tf(...)) / tf2ss (companion form).
     * @see to_observable_canonical_form for the dual realization
     */
    [[nodiscard]] constexpr damp::optional<StateSpace<Nden - 1, 1, 1, T>> to_state_space() const {
        static_assert(Nnum <= Nden, "TransferFunction must be proper (deg num <= deg den)");

        StateSpace<Nden - 1, 1, 1, T> sys{};

        if constexpr (Nden == 1) {
            if constexpr (Nnum == 1) {
                if (den[0] == T{0}) {
                    return damp::nullopt;
                }
                sys.D(0, 0) = num[0] / den[0];
            }
            return sys;
        }

        const T a_lead = den[0];
        if (a_lead == T{0}) {
            return damp::nullopt;
        }

        // The companion form indexes by power: a[k] and b[k] multiply s^k (num aligned on its
        // constant term, zero above its degree).
        damp::array<T, Nden> a{};
        damp::array<T, Nden> b{};
        for (size_t k = 0; k < Nden; ++k) {
            a[k] = den[Nden - 1 - k];
        }
        for (size_t k = 0; k < Nnum; ++k) {
            b[k] = num[Nnum - 1 - k];
        }

        for (size_t i = 0; i < Nden - 2; ++i) {
            sys.A(i, i + 1) = T{1};
        }
        for (size_t j = 0; j < Nden - 1; ++j) {
            sys.A(Nden - 2, j) = -a[j] / a_lead;
        }

        sys.B(Nden - 2, 0) = T{1};

        const T d_val = (Nnum == Nden) ? (b[Nden - 1] / a_lead) : T{0};
        sys.D(0, 0) = d_val;

        for (size_t k = 0; k < Nden - 1; ++k) {
            sys.C(0, k) = (b[k] - d_val * a[k]) / a_lead;
        }

        return sys;
    }

    /**
     * @copybrief to_state_space
     * The controllable companion form by its canonical name; see @ref to_state_space.
     */
    [[nodiscard]] constexpr damp::optional<StateSpace<Nden - 1, 1, 1, T>> to_controllable_canonical_form() const {
        return to_state_space();
    }

    /**
     * @brief Observable companion-form state-space realization
     *
     * The dual of @ref to_state_space(): @f$ A_o = A_c^\top,\; B_o = C_c^\top,\; C_o = B_c^\top @f$,
     * @f$ D @f$ unchanged, so the denominator coefficients sit in the last column of A and
     * the numerator coefficients in B. Same transfer function and poles; (A, C) is observable
     * by construction, which suits an observer or Kalman design the way the controllable
     * form suits state feedback.
     *
     * @return Continuous realization (@c Ts = 0), or @c nullopt as @ref to_state_space
     * @note Compare with MATLAB®'s canon(sys, 'companion') on the transposed system.
     * @see "Linear System Theory and Design" (Chen, 1999), §7.2
     */
    [[nodiscard]] constexpr damp::optional<StateSpace<Nden - 1, 1, 1, T>> to_observable_canonical_form() const {
        const auto ccf = to_state_space();
        if (!ccf) {
            return damp::nullopt;
        }
        StateSpace<Nden - 1, 1, 1, T> sys{};
        sys.A = ccf->A.transpose();
        sys.B = ccf->C.transpose();
        sys.C = ccf->B.transpose();
        sys.D = ccf->D;
        return sys;
    }
};

template<typename TNum, size_t Nnum, typename TDen, size_t Nden>
TransferFunction(const damp::array<TNum, Nnum>&, const damp::array<TDen, Nden>&)
    -> TransferFunction<Nnum, Nden, transfer_function_scalar_t<TNum, TDen>>;

template<typename TNum, size_t Nnum, typename TDen, size_t Nden>
TransferFunction(const TNum (&)[Nnum], const TDen (&)[Nden])
    -> TransferFunction<Nnum, Nden, transfer_function_scalar_t<TNum, TDen>>;

template<size_t Nnum1, size_t Nden1, size_t Nnum2, size_t Nden2, typename T>
[[nodiscard]] constexpr auto operator*(
    const TransferFunction<Nnum1, Nden1, T>& tf1,
    const TransferFunction<Nnum2, Nden2, T>& tf2
) {
    constexpr size_t Nnum_res = Nnum1 + Nnum2 - 1;
    constexpr size_t Nden_res = Nden1 + Nden2 - 1;

    return TransferFunction<Nnum_res, Nden_res, T>{
        tf_detail::convolve_poly(tf1.num, tf2.num),
        tf_detail::convolve_poly(tf1.den, tf2.den),
    };
}

template<size_t Nnum1, size_t Nden1, size_t Nnum2, size_t Nden2, typename T>
[[nodiscard]] constexpr auto operator+(
    const TransferFunction<Nnum1, Nden1, T>& tf1,
    const TransferFunction<Nnum2, Nden2, T>& tf2
) {
    constexpr size_t Nnum_l = Nnum1 + Nden2 - 1;
    constexpr size_t Nnum_r = Nnum2 + Nden1 - 1;
    constexpr size_t Nnum_res = (Nnum_l > Nnum_r) ? Nnum_l : Nnum_r;
    constexpr size_t Nden_res = Nden1 + Nden2 - 1;

    const auto left_num = tf_detail::convolve_poly(tf1.num, tf2.den);
    const auto right_num = tf_detail::convolve_poly(tf2.num, tf1.den);
    const auto den = tf_detail::convolve_poly(tf1.den, tf2.den);

    damp::array<T, Nnum_res> num{};
    tf_detail::accumulate_poly(num, left_num);
    tf_detail::accumulate_poly(num, right_num);

    return TransferFunction<Nnum_res, Nden_res, T>{num, den};
}

template<size_t Nnum1, size_t Nden1, size_t Nnum2, size_t Nden2, typename T>
[[nodiscard]] constexpr auto operator-(
    const TransferFunction<Nnum1, Nden1, T>& tf1,
    const TransferFunction<Nnum2, Nden2, T>& tf2
) {
    constexpr size_t Nnum_l = Nnum1 + Nden2 - 1;
    constexpr size_t Nnum_r = Nnum2 + Nden1 - 1;
    constexpr size_t Nnum_res = (Nnum_l > Nnum_r) ? Nnum_l : Nnum_r;
    constexpr size_t Nden_res = Nden1 + Nden2 - 1;

    const auto left_num = tf_detail::convolve_poly(tf1.num, tf2.den);
    const auto right_num = tf_detail::convolve_poly(tf2.num, tf1.den);
    const auto den = tf_detail::convolve_poly(tf1.den, tf2.den);

    damp::array<T, Nnum_res> num{};
    tf_detail::accumulate_poly(num, left_num);
    tf_detail::accumulate_poly(num, right_num, T{-1});

    return TransferFunction<Nnum_res, Nden_res, T>{num, den};
}

template<size_t Nnum1, size_t Nden1, size_t Nnum2, size_t Nden2, typename T>
[[nodiscard]] constexpr auto series(
    const TransferFunction<Nnum1, Nden1, T>& tf1,
    const TransferFunction<Nnum2, Nden2, T>& tf2
) {
    return tf1 * tf2;
}

template<size_t Nnum1, size_t Nden1, size_t Nnum2, size_t Nden2, typename T>
[[nodiscard]] constexpr auto parallel(
    const TransferFunction<Nnum1, Nden1, T>& tf1,
    const TransferFunction<Nnum2, Nden2, T>& tf2
) {
    return tf1 + tf2;
}

template<size_t Nnum1, size_t Nden1, size_t Nnum2, size_t Nden2, typename T>
[[nodiscard]] constexpr auto subtract(
    const TransferFunction<Nnum1, Nden1, T>& tf1,
    const TransferFunction<Nnum2, Nden2, T>& tf2
) {
    return tf1 - tf2;
}

/**
 * @brief Negative feedback of two TFs: G/(1+GH) (polynomial form)
 *
 * Expands without cancellation; pair with ZPK cancel / design::minreal after
 * realization when a minimal order is needed.
 *
 * @note Compare with MATLAB®'s feedback(sys1,sys2).
 */
template<size_t Nnum1, size_t Nden1, size_t Nnum2, size_t Nden2, typename T>
[[nodiscard]] constexpr auto feedback(
    const TransferFunction<Nnum1, Nden1, T>& sys1,
    const TransferFunction<Nnum2, Nden2, T>& sys2
) {
    // Negative feedback: sys1 / (1 + sys1*sys2)
    constexpr size_t Nnum_res = Nnum1 + Nden2 - 1;
    constexpr size_t Nden_l = Nden1 + Nden2 - 1;
    constexpr size_t Nden_r = Nnum1 + Nnum2 - 1;
    constexpr size_t Nden_res = (Nden_l > Nden_r) ? Nden_l : Nden_r;

    const auto num = tf_detail::convolve_poly(sys1.num, sys2.den);
    const auto den_l = tf_detail::convolve_poly(sys1.den, sys2.den);
    const auto den_r = tf_detail::convolve_poly(sys1.num, sys2.num);

    damp::array<T, Nden_res> den{};
    tf_detail::accumulate_poly(den, den_l);
    tf_detail::accumulate_poly(den, den_r);

    return TransferFunction<Nnum_res, Nden_res, T>{num, den};
}

} // namespace damp
