// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file transfer_function.hpp
 * @brief SISO polynomial transfer function G(s) or G(z) = num/den
 *
 * Coefficients are in ascending powers (constant first), matching the rest of Damp.
 * There is no @c Ts field: continuous vs discrete is the caller's convention
 * (same as ZPK). Realize with to_state_space then discretize for c2d.
 *
 * Algebra (@c *, @c +, @c feedback) expands polynomials without cancellation —
 * use design::minreal or ZPK @c cancel_matching when degree growth matters.
 *
 * @note Compare with MATLAB®'s tf (note: MATLAB® uses descending powers by default).
 * @see TransferFunction::to_state_space, systems/state_space.hpp, systems/zpk.hpp
 */

#include <cstddef>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix_traits.hpp"
#include "state_space.hpp"

namespace damp {

namespace tf_detail {

template<size_t Nres, size_t Nsrc, typename T>
constexpr void accumulate_poly(damp::array<T, Nres>& dst, const damp::array<T, Nsrc>& src, T scale = T{1}) {
    constexpr size_t N = (Nsrc < Nres) ? Nsrc : Nres;
    for (size_t i = 0; i < N; ++i) {
        dst[i] += scale * src[i];
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

} // namespace tf_detail

template<typename TNum, typename TDen>
using transfer_function_scalar_t = std::conditional_t<
    std::is_floating_point_v<std::common_type_t<TNum, TDen>>,
    std::common_type_t<TNum, TDen>,
    double>;

/**
 * @brief SISO polynomial transfer function G(s) = num(s)/den(s)
 *
 * Coefficients are in ascending powers of s (or z): den[0] is the constant
 * term, den[Nden-1] the highest power. Same convention for num.
 *
 * @tparam Nnum Number of numerator coefficients
 * @tparam Nden Number of denominator coefficients
 * @tparam T    Scalar type (default double)
 */
template<size_t Nnum, size_t Nden, typename T = double>
    requires std::is_floating_point_v<T>
struct TransferFunction {
    damp::array<T, Nnum> num{}; ///< Numerator coefficients, ascending powers of s
    damp::array<T, Nden> den{}; ///< Denominator coefficients, ascending powers of s

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
     * denominator coefficient. Biproper systems set @f$ D = n_n/d_n @f$ and reduce C.
     *
     * @return Continuous realization (@c Ts = 0), or @c nullopt if
     *         @f$ |\mathrm{den}[\mathrm{end}]| @f$ is below default_tol
     * @note Compare with MATLAB®'s ss(tf(...)) / tf2ss (companion form).
     */
    [[nodiscard]] constexpr damp::optional<StateSpace<Nden - 1, 1, 1, T>> to_state_space() const {
        static_assert(Nnum <= Nden, "TransferFunction must be proper (deg num <= deg den)");

        StateSpace<Nden - 1, 1, 1, T> sys{};

        if constexpr (Nden == 1) {
            if constexpr (Nnum == 1) {
                if (!(damp::abs(den[0]) > default_tol<T>())) {
                    return damp::nullopt;
                }
                sys.D(0, 0) = num[0] / den[0];
            }
            return sys;
        }

        const T a_lead = den[Nden - 1];
        if (!(damp::abs(a_lead) > default_tol<T>())) {
            return damp::nullopt;
        }

        for (size_t i = 0; i < Nden - 2; ++i) {
            sys.A(i, i + 1) = T{1};
        }
        for (size_t j = 0; j < Nden - 1; ++j) {
            sys.A(Nden - 2, j) = -den[j] / a_lead;
        }

        sys.B(Nden - 2, 0) = T{1};

        const T d_val = (Nnum == Nden) ? (num[Nden - 1] / a_lead) : T{0};
        sys.D(0, 0) = d_val;

        for (size_t k = 0; k < Nden - 1; ++k) {
            const T n_k = (k < Nnum) ? num[k] : T{0};
            sys.C(0, k) = (n_k - d_val * den[k]) / a_lead;
        }

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
 * @brief Negative feedback of two TFs: @f$ G/(1+GH) @f$ (polynomial form)
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

template<size_t Nnum1, size_t Nden1, size_t Nnum2, size_t Nden2, typename T>
[[nodiscard]] constexpr auto operator/(
    const TransferFunction<Nnum1, Nden1, T>& sys1,
    const TransferFunction<Nnum2, Nden2, T>& sys2
) {
    return feedback(sys1, sys2);
}

} // namespace damp
