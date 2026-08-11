// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file zpk.hpp
 * @brief Zero-pole-gain (ZPK) SISO LTI representation
 *
 * # Upgrade path (TF / SS / ZPK)
 *
 * | From → To | API | Notes |
 * | --------- | --- | ----- |
 * | ZPK → TF  | ZPK::to_transfer_function / zpk2tf | Always succeeds; expands factored form |
 * | ZPK → SS  | ZPK::to_state_space | Via TF companion realization |
 * | TF → ZPK  | to_zpk / tf2zpk | Poly roots via companion + eigen; @c success if leads ≠ 0 and QR converges |
 * | SS → ZPK  | to_zpk (SISO) / ss2zpk | Biproper: zeros = eig(@f$A-BD^{-1}C@f$), @f$k=D@f$. Strict: Leverrier num roots |
 * | SS → TF   | to_transfer_function / ss2tf | Faddeev–LeVerrier (SISO) |
 *
 * Domain (continuous vs discrete). Like TransferFunction, @c ZPK does not
 * store @c Ts. Factors are in @f$s@f$ or @f$z@f$ per the caller's convention. Use
 * ZPK::dcgain for continuous @f$H(0)@f$ and ZPK::dcgain_discrete for
 * discrete @f$H(1)@f$. StateSpace carries @c Ts; SS→ZPK drops it (re-attach
 * on a realized SS if needed).
 *
 * Limits. Real polynomial roots use the companion-matrix eigenproblem
 * (mat::compute_eigenvalues): closed form for degree ≤ 2, Francis QR for
 * higher degree (design-time; non-convergence → @c success=false). Leading
 * numerator/denominator coefficients must be nonzero for TF→ZPK. Strictly
 * proper SS→ZPK recovers finite zeros into StateSpaceZPKResult (runtime
 * @c n_zeros); call StateSpaceZPKResult::as_zpk when the count is known at
 * compile time. High-order root finding (degree ≫ 8–12) is the same QR path as
 * the rest of the library — solid for typical plant orders, not a multiprecision
 * polynomial solver.
 *
 * Example: round-trip a lead compensator
 * @code
 * #include "damp/systems/zpk.hpp"
 * using namespace damp;
 *
 * constexpr TransferFunction<2, 2> tf{.num = {2.0, 1.0}, .den = {10.0, 1.0}}; // (s+2)/(s+10)
 * constexpr auto art = to_zpk(tf);
 * static_assert(art.success);
 * constexpr auto tf2 = art.zpk.to_transfer_function();
 * @endcode
 *
 * @note Compare with MATLAB®'s zpk, zpkdata, tf2zpk, zpk2tf, ss2tf, minreal (factored form).
 */

#include <cstddef>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/eigen.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/matrix_traits.hpp"
#include "state_space.hpp"
#include "transfer_function.hpp"

namespace damp {

// ============================================================================
// Fixed-size polynomial roots (embeddable)
// ============================================================================

/**
 * @brief Roots of a real polynomial given in ascending powers
 *
 * For coefficients @f$ c_0 + c_1 x + \cdots + c_{N-1} x^{N-1} @f$, returns the
 * @f$ N-1 @f$ roots as eigenvalues of the companion matrix. Requires
 * @f$ |c_{N-1}| @f$ above @p tol (nonzero leading coefficient).
 *
 * - Degree 0 (@c N == 1): empty root set, @c success = true
 * - Degree 1–2: closed form
 * - Higher degree: Francis QR via mat::compute_eigenvalues; @c success
 *   mirrors the eigen convergence flag
 *
 * @note Compare with MATLAB®'s roots(flip(c)) for ascending @p c.
 * @see to_zpk() for the TF factorization path that uses this
 *
 * @tparam N Number of coefficients (degree N−1 when N ≥ 1)
 * @tparam T Floating-point scalar type
 */
template<size_t N, typename T = double>
struct PolyRootsResult {
    damp::array<damp::complex<T>, (N > 0 ? N - 1 : 0)> roots{};        ///< Roots (empty when N ≤ 1)
    bool                                               success{false}; ///< true if leading coeff OK and eigen converged

    template<typename U>
    [[nodiscard]] constexpr PolyRootsResult<N, U> as() const {
        PolyRootsResult<N, U> out{};
        out.success = success;
        if constexpr (N > 1) {
            for (size_t i = 0; i < N - 1; ++i) {
                out.roots[i] = damp::complex<U>{
                    static_cast<U>(roots[i].real()),
                    static_cast<U>(roots[i].imag())
                };
            }
        }
        return out;
    }
};

/**
 * @brief Compute roots of an ascending-power real polynomial
 *
 * @param c   Coefficients c₀ + c₁ x + ·s + c_N-1 x^N-1
 * @param tol Leading-coefficient magnitude floor (default: default_tol)
 * @return @ref PolyRootsResult
 */
template<size_t N, typename T = double>
    requires std::is_floating_point_v<T>
[[nodiscard]] constexpr PolyRootsResult<N, T> poly_roots(
    const damp::array<T, N>& c,
    T                        tol = default_tol<T>()
) {
    if constexpr (N <= 1) {
        return PolyRootsResult<N, T>{.success = true};
    } else if constexpr (N == 2) {
        if (!(damp::abs(c[1]) > tol)) {
            return PolyRootsResult<N, T>{};
        }
        return PolyRootsResult<N, T>{
            .roots = {damp::complex<T>{-c[0] / c[1], T{0}}},
            .success = true,
        };
    } else {
        if (!(damp::abs(c[N - 1]) > tol)) {
            return PolyRootsResult<N, T>{};
        }
        constexpr size_t M = N - 1;
        Matrix<M, M, T>  companion{};
        for (size_t i = 0; i + 1 < M; ++i) {
            companion(i, i + 1) = T{1};
        }
        const T inv_lead = T{1} / c[N - 1];
        for (size_t j = 0; j < M; ++j) {
            companion(M - 1, j) = -c[j] * inv_lead;
        }
        const auto            eig = mat::compute_eigenvalues(companion);
        PolyRootsResult<N, T> result{};
        for (size_t i = 0; i < M; ++i) {
            result.roots[i] = eig.values[i];
        }
        result.success = eig.converged;
        return result;
    }
}

// ============================================================================
// ZPK type
// ============================================================================

/**
 * @brief Zero-pole-gain (ZPK) representation of a SISO LTI system
 *
 * Represents the factored transfer function
 * @f[
 *   H(s) = k \frac{\prod_{i} (s - z_i)}{\prod_{j} (s - p_j)}
 * @f]
 * with zeros @f$z_i@f$, poles @f$p_j@f$, and scalar gain @f$k@f$. The same form
 * holds in @f$z@f$ for discrete models (caller tracks the domain; no @c Ts field).
 *
 * MATLAB® equivalent: `zpk(z, p, k)`.
 *
 * @tparam Nz   Number of zeros
 * @tparam Np   Number of poles
 * @tparam T    Floating-point scalar type
 */
template<size_t Nz, size_t Np, typename T = double>
    requires std::is_floating_point_v<T>
struct ZPK {
    damp::array<damp::complex<T>, Nz> zeros{}; ///< Zeros of the transfer function
    damp::array<damp::complex<T>, Np> poles{}; ///< Poles of the transfer function
    T                                 gain{1}; ///< Gain k of the monic factored form

    /**
     * @brief Result of @ref cancel_matching (factored minreal)
     */
    struct CancelResult {
        ZPK    zpk{};          ///< Model with cancelled pairs re-appended as identical z=p
        size_t n_cancelled{0}; ///< Number of pole-zero pairs cancelled
    };

    template<typename U>
    [[nodiscard]] constexpr ZPK<Nz, Np, U> as() const {
        ZPK<Nz, Np, U> result{};
        for (size_t i = 0; i < Nz; ++i) {
            result.zeros[i] = damp::complex<U>{static_cast<U>(zeros[i].real()), static_cast<U>(zeros[i].imag())};
        }
        for (size_t j = 0; j < Np; ++j) {
            result.poles[j] = damp::complex<U>{static_cast<U>(poles[j].real()), static_cast<U>(poles[j].imag())};
        }
        result.gain = static_cast<U>(gain);
        return result;
    }

    /**
     * @brief Convert ZPK representation to transfer function form
     *
     * Expands @f$ k \prod(s - z_i) / \prod(s - p_j) @f$ into ascending-power real
     * polynomials. Intermediate products use complex coefficients so conjugate
     * pairs cancel imaginary parts; the returned TF stores the real parts.
     *
     * @note Compare with MATLAB®'s tf(zpk_sys) / zpk2tf.
     * @return Transfer function with Nnum = Nz+1, Nden = Np+1
     */
    [[nodiscard]] constexpr TransferFunction<Nz + 1, Np + 1, T> to_transfer_function() const {
        TransferFunction<Nz + 1, Np + 1, T> tf{};

        auto mul_linear = [](auto& poly, size_t degree, damp::complex<T> root) {
            poly[degree + 1] = poly[degree];
            for (size_t j = degree; j > 0; --j) {
                poly[j] = poly[j - 1] - root * poly[j];
            }
            poly[0] = -root * poly[0];
        };

        damp::array<damp::complex<T>, Nz + 1> num_c{};
        num_c[0] = damp::complex<T>{gain, T{0}};
        for (size_t i = 0; i < Nz; ++i) {
            mul_linear(num_c, i, zeros[i]);
        }
        for (size_t k = 0; k < Nz + 1; ++k) {
            tf.num[k] = num_c[k].real();
        }

        damp::array<damp::complex<T>, Np + 1> den_c{};
        den_c[0] = damp::complex<T>{T{1}, T{0}};
        for (size_t i = 0; i < Np; ++i) {
            mul_linear(den_c, i, poles[i]);
        }
        for (size_t k = 0; k < Np + 1; ++k) {
            tf.den[k] = den_c[k].real();
        }

        return tf;
    }

    /**
     * @brief Convert ZPK to state-space (companion realization of the TF)
     *
     * @note Compare with MATLAB®'s ss(zpk_sys).
     * @return Continuous StateSpace (Ts = 0) with NX = Np
     */
    /**
     * @brief Companion realization via TF (leading den is monic after expansion)
     */
    [[nodiscard]] constexpr StateSpace<Np, 1, 1, T> to_state_space() const {
        auto ss = to_transfer_function().to_state_space();
        // Expanded monic den has leading 1; companion always succeeds for valid ZPK.
        return ss ? *ss : StateSpace<Np, 1, 1, T>{};
    }

    /**
     * @brief Evaluate H(s) = k ∏(s−zᵢ)/∏(s−pⱼ)
     *
     * @param s   Complex frequency (s or z)
     * @param tol Pole-coincidence tolerance
     * @return H(s), or nullopt if any pole lies within @p tol of @p s
     */
    [[nodiscard]] constexpr damp::optional<damp::complex<T>> eval(damp::complex<T> s, T tol = default_tol<T>()) const {
        damp::complex<T> num{gain, T{0}};
        damp::complex<T> den{T{1}, T{0}};
        for (size_t i = 0; i < Nz; ++i) {
            num = num * (s - zeros[i]);
        }
        for (size_t j = 0; j < Np; ++j) {
            const damp::complex<T> d = s - poles[j];
            if (damp::abs(d) <= tol) {
                return damp::nullopt;
            }
            den = den * d;
        }
        return num / den;
    }

    /**
     * @brief Continuous-time DC gain H(0)
     *
     * @return Real DC gain, or nullopt if a pole is at the origin
     * @note Compare with MATLAB®'s dcgain for continuous zpk models.
     */
    [[nodiscard]] constexpr damp::optional<T> dcgain(T tol = default_tol<T>()) const {
        const auto h = eval(damp::complex<T>{T{0}, T{0}}, tol);
        if (!h) {
            return damp::nullopt;
        }
        return h->real();
    }

    /**
     * @brief Discrete-time DC gain H(1)
     *
     * @return Real DC gain, or nullopt if a pole is at z = 1
     * @note Compare with MATLAB®'s dcgain for discrete zpk models.
     */
    [[nodiscard]] constexpr damp::optional<T> dcgain_discrete(T tol = default_tol<T>()) const {
        const auto h = eval(damp::complex<T>{T{1}, T{0}}, tol);
        if (!h) {
            return damp::nullopt;
        }
        return h->real();
    }

    /**
     * @brief Cancel matching pole-zero pairs within tolerance (SISO minreal on factors)
     *
     * Greedy one-to-one matching: each zero cancels at most one pole with
     * @f$ |z-p| \le tol @f$. Survivors are compacted to the front of the arrays;
     * each cancelled value is re-appended as an identical @f$ z = p @f$ pair so
     * to_transfer_function remains algebraically equivalent (common factors),
     * while @c n_cancelled reports how many pairs were cancelled.
     *
     * Gain is unchanged (cancelled monic factors contribute 1).
     *
     * @note Compare with MATLAB®'s minreal on a zpk model (factored form only;
     *       does not perform SS uncontrollable/unobservable reduction).
     * @param tol Match tolerance on |z - p| (default: default_tol)
     * @return @ref CancelResult
     */
    [[nodiscard]] constexpr CancelResult cancel_matching(T tol = default_tol<T>()) const {
        CancelResult out{};
        out.zpk.gain = gain;

        damp::array<bool, Nz>                              z_gone{};
        damp::array<bool, Np>                              p_gone{};
        damp::array<damp::complex<T>, (Nz < Np ? Nz : Np)> cancelled{};
        size_t                                             n_cancel = 0;

        for (size_t i = 0; i < Nz; ++i) {
            for (size_t j = 0; j < Np; ++j) {
                if (z_gone[i] || p_gone[j]) {
                    continue;
                }
                if (damp::abs(zeros[i] - poles[j]) <= tol) {
                    z_gone[i] = true;
                    p_gone[j] = true;
                    cancelled[n_cancel] = zeros[i];
                    ++n_cancel;
                    break;
                }
            }
        }

        size_t zi = 0;
        size_t pi = 0;
        for (size_t i = 0; i < Nz; ++i) {
            if (!z_gone[i]) {
                out.zpk.zeros[zi++] = zeros[i];
            }
        }
        for (size_t j = 0; j < Np; ++j) {
            if (!p_gone[j]) {
                out.zpk.poles[pi++] = poles[j];
            }
        }
        for (size_t c = 0; c < n_cancel; ++c) {
            if (zi < Nz) {
                out.zpk.zeros[zi++] = cancelled[c];
            }
            if (pi < Np) {
                out.zpk.poles[pi++] = cancelled[c];
            }
        }

        out.n_cancelled = n_cancel;
        return out;
    }

    /**
     * @brief Factored minreal: same as @ref cancel_matching but returns only the ZPK
     * @note Compare with MATLAB®'s minreal(zpk_sys).
     */
    [[nodiscard]] constexpr ZPK minreal(T tol = default_tol<T>()) const {
        return cancel_matching(tol).zpk;
    }
};

// ============================================================================
// TF → ZPK
// ============================================================================

/**
 * @brief Result of converting a transfer function to zero-pole-gain form
 *
 * @tparam Nz Number of zeros (= Nnum − 1)
 * @tparam Np Number of poles (= Nden − 1)
 * @tparam T  Scalar type
 *
 * @note Compare with MATLAB®'s [z,p,k] = zpkdata(tf) / tf2zpk.
 */
template<size_t Nz, size_t Np, typename T = double>
struct ZPKResult {
    ZPK<Nz, Np, T> zpk{};          ///< Factored model (valid when success)
    bool           success{false}; ///< true if both poly root solves succeeded

    template<typename U>
    [[nodiscard]] constexpr ZPKResult<Nz, Np, U> as() const {
        return ZPKResult<Nz, Np, U>{zpk.template as<U>(), success};
    }
};

/**
 * @brief Convert a SISO transfer function to zero-pole-gain form
 *
 * Factors @f$ G = \mathrm{num}/\mathrm{den} @f$ by rooting both real polynomials.
 * The ZPK gain is @f$ k = a_{\mathrm{num,lead}} / a_{\mathrm{den,lead}} @f$ so
 * @f$ G(s) = k \prod(s-z)/\prod(s-p) @f$ with monic factors.
 *
 * Requires nonzero leading coefficients. Degree ≥ 3 uses iterative QR (design-time).
 *
 * @note Compare with MATLAB®'s tf2zpk / zpk(tf_sys).
 * @see ZPK::to_transfer_function for the inverse map
 * @see poly_roots
 *
 * @param tf  Ascending-power transfer function
 * @param tol Leading-coefficient floor and root-solver tolerance
 * @return ZPKResult with success flag
 */
template<size_t Nnum, size_t Nden, typename T = double>
    requires std::is_floating_point_v<T>
[[nodiscard]] constexpr ZPKResult<(Nnum > 0 ? Nnum - 1 : 0), (Nden > 0 ? Nden - 1 : 0), T>
to_zpk(const TransferFunction<Nnum, Nden, T>& tf, T tol = default_tol<T>()) {
    constexpr size_t     Nz = (Nnum > 0 ? Nnum - 1 : 0);
    constexpr size_t     Np = (Nden > 0 ? Nden - 1 : 0);
    ZPKResult<Nz, Np, T> art{};

    if constexpr (Nden == 0 || Nnum == 0) {
        return art;
    }

    if (!(damp::abs(tf.den[Nden - 1]) > tol)) {
        return art;
    }
    if (!(damp::abs(tf.num[Nnum - 1]) > tol)) {
        return art;
    }

    const auto zr = poly_roots(tf.num, tol);
    const auto pr = poly_roots(tf.den, tol);
    if (!zr.success || !pr.success) {
        return art;
    }

    art.zpk.gain = tf.num[Nnum - 1] / tf.den[Nden - 1];
    if constexpr (Nz > 0) {
        for (size_t i = 0; i < Nz; ++i) {
            art.zpk.zeros[i] = zr.roots[i];
        }
    }
    if constexpr (Np > 0) {
        for (size_t j = 0; j < Np; ++j) {
            art.zpk.poles[j] = pr.roots[j];
        }
    }
    art.success = true;
    return art;
}

/**
 * @brief MATLAB®-style alias for to_zpk (transfer function)
 * @note Compare with MATLAB®'s tf2zpk.
 */
template<size_t Nnum, size_t Nden, typename T = double>
[[nodiscard]] constexpr auto tf2zpk(const TransferFunction<Nnum, Nden, T>& tf, T tol = default_tol<T>()) {
    return to_zpk(tf, tol);
}

/**
 * @brief MATLAB®-style alias for ZPK::to_transfer_function
 * @note Compare with MATLAB®'s zpk2tf / tf(zpk_sys).
 */
template<size_t Nz, size_t Np, typename T = double>
[[nodiscard]] constexpr TransferFunction<Nz + 1, Np + 1, T> zpk2tf(const ZPK<Nz, Np, T>& sys) {
    return sys.to_transfer_function();
}

// ============================================================================
// SS → TF / ZPK (SISO)
// ============================================================================

namespace zpk_detail {

/**
 * @brief Faddeev–LeVerrier characteristic polynomial and adjoint coefficient matrices
 *
 * Builds monic @f$ \det(sI-A) @f$ (ascending in @c den) and @f$ M_1,\ldots,M_n @f$
 * with @f$ M_1 = I @f$ such that
 * @f$ \mathrm{adj}(sI-A) = \sum_{j=0}^{n-1} s^{n-1-j} M_{j+1} @f$.
 *
 * @see Faddeev–LeVerrier; Golub & Van Loan, Matrix Computations
 */
template<size_t NX, typename T>
struct LeverrierResult {
    damp::array<T, NX + 1>             den{}; ///< Ascending monic char poly
    damp::array<Matrix<NX, NX, T>, NX> M{};   ///< M_1 .. M_NX
};

template<size_t NX, typename T>
[[nodiscard]] constexpr LeverrierResult<NX, T> leverrier(const Matrix<NX, NX, T>& A) {
    if constexpr (NX == 0) {
        return LeverrierResult<NX, T>{.den = {T{1}}};
    } else {
        LeverrierResult<NX, T> r{};

        const auto        I = Matrix<NX, NX, T>::identity();
        Matrix<NX, NX, T> Mk = Matrix<NX, NX, T>::zeros();
        T                 c_prev = T{1};
        r.den[NX] = T{1};

        for (size_t k = 1; k <= NX; ++k) {
            Mk = A * Mk + c_prev * I;
            r.M[k - 1] = Mk;
            const T ck = -(A * Mk).trace() / static_cast<T>(k);
            r.den[NX - k] = ck;
            c_prev = ck;
        }
        return r;
    }
}

template<size_t NX, typename T>
[[nodiscard]] constexpr TransferFunction<NX + 1, NX + 1, T>
ss_to_tf_full(const StateSpace<NX, 1, 1, T>& sys) {
    if constexpr (NX == 0) {
        return TransferFunction<NX + 1, NX + 1, T>{
            {sys.D(0, 0)},
            {T{1}},
        };
    } else {
        TransferFunction<NX + 1, NX + 1, T> tf{};
        const auto                          lev = leverrier(sys.A);
        tf.den = lev.den;

        const T D0 = sys.D(0, 0);
        for (size_t k = 0; k <= NX; ++k) {
            tf.num[k] = D0 * tf.den[k];
        }

        const ColVec<NX, T> b = sys.B;
        for (size_t j = 0; j < NX; ++j) {
            const ColVec<NX, T> Mb = lev.M[j] * b;
            T                   cMb = T{0};
            for (size_t r = 0; r < NX; ++r) {
                cMb += sys.C(0, r) * Mb[r];
            }
            tf.num[NX - 1 - j] += cMb;
        }
        return tf;
    }
}

} // namespace zpk_detail

/**
 * @brief SISO state-space → ZPK conversion result with runtime zero count
 *
 * Poles fill @c poles[0..NX). Finite zeros fill @c zeros[0..n_zeros).
 * Use @ref as_zpk when @c n_zeros is known as a compile-time @c Nz.
 *
 * @note Compare with MATLAB®'s zpk(ss_sys) / ss2zpk.
 */
template<size_t NX, typename T = double>
struct StateSpaceZPKResult {
    damp::array<damp::complex<T>, NX> poles{};
    damp::array<damp::complex<T>, NX> zeros{}; ///< First @ref n_zeros entries valid
    size_t                            n_zeros{0};
    T                                 gain{T{1}};
    bool                              success{false};

    template<typename U>
    [[nodiscard]] constexpr StateSpaceZPKResult<NX, U> as() const {
        StateSpaceZPKResult<NX, U> o{};
        o.n_zeros = n_zeros;
        o.gain = static_cast<U>(gain);
        o.success = success;
        for (size_t i = 0; i < NX; ++i) {
            o.poles[i] = damp::complex<U>{static_cast<U>(poles[i].real()), static_cast<U>(poles[i].imag())};
            o.zeros[i] = damp::complex<U>{static_cast<U>(zeros[i].real()), static_cast<U>(zeros[i].imag())};
        }
        return o;
    }

    /**
     * @brief Materialize a fixed-size ZPK when @p Nz equals @ref n_zeros
     */
    template<size_t Nz>
    [[nodiscard]] constexpr damp::optional<ZPK<Nz, NX, T>> as_zpk() const {
        if (!success || Nz != n_zeros) {
            return damp::nullopt;
        }
        ZPK<Nz, NX, T> z{};
        z.gain = gain;
        for (size_t i = 0; i < Nz; ++i) {
            z.zeros[i] = zeros[i];
        }
        for (size_t j = 0; j < NX; ++j) {
            z.poles[j] = poles[j];
        }
        return z;
    }
};

namespace zpk_detail {

template<size_t Deg, size_t NX, typename T>
[[nodiscard]] constexpr bool fill_zeros_from_num(
    const damp::array<T, NX + 1>& num,
    StateSpaceZPKResult<NX, T>&   art,
    T                             tol
) {
    static_assert(Deg >= 1 && Deg <= NX);
    damp::array<T, Deg + 1> coeffs{};
    for (size_t i = 0; i <= Deg; ++i) {
        coeffs[i] = num[i];
    }
    const auto zr = poly_roots(coeffs, tol);
    if (!zr.success) {
        return false;
    }
    for (size_t i = 0; i < Deg; ++i) {
        art.zeros[i] = zr.roots[i];
    }
    return true;
}

template<size_t NX, typename T>
[[nodiscard]] constexpr bool fill_zeros_dispatch(
    size_t                        deg,
    const damp::array<T, NX + 1>& num,
    StateSpaceZPKResult<NX, T>&   art,
    T                             tol
) {
    if (deg == 0) {
        return true;
    }
    // Explicit ladder keeps constexpr-friendly fixed sizes (typical plant orders).
    if constexpr (NX >= 1) {
        if (deg == 1) {
            return fill_zeros_from_num<1, NX, T>(num, art, tol);
        }
    }
    if constexpr (NX >= 2) {
        if (deg == 2) {
            return fill_zeros_from_num<2, NX, T>(num, art, tol);
        }
    }
    if constexpr (NX >= 3) {
        if (deg == 3) {
            return fill_zeros_from_num<3, NX, T>(num, art, tol);
        }
    }
    if constexpr (NX >= 4) {
        if (deg == 4) {
            return fill_zeros_from_num<4, NX, T>(num, art, tol);
        }
    }
    if constexpr (NX >= 5) {
        if (deg == 5) {
            return fill_zeros_from_num<5, NX, T>(num, art, tol);
        }
    }
    if constexpr (NX >= 6) {
        if (deg == 6) {
            return fill_zeros_from_num<6, NX, T>(num, art, tol);
        }
    }
    if constexpr (NX >= 7) {
        if (deg == 7) {
            return fill_zeros_from_num<7, NX, T>(num, art, tol);
        }
    }
    if constexpr (NX >= 8) {
        if (deg == 8) {
            return fill_zeros_from_num<8, NX, T>(num, art, tol);
        }
    }
    if constexpr (NX >= 9) {
        if (deg == 9) {
            return fill_zeros_from_num<9, NX, T>(num, art, tol);
        }
    }
    if constexpr (NX >= 10) {
        if (deg == 10) {
            return fill_zeros_from_num<10, NX, T>(num, art, tol);
        }
    }
    if constexpr (NX >= 11) {
        if (deg == 11) {
            return fill_zeros_from_num<11, NX, T>(num, art, tol);
        }
    }
    if constexpr (NX >= 12) {
        if (deg == 12) {
            return fill_zeros_from_num<12, NX, T>(num, art, tol);
        }
    }
    // deg == NX (any size): use full-order poly_roots on NX+1 coeffs
    if (deg == NX) {
        return fill_zeros_from_num<NX, NX, T>(num, art, tol);
    }
    return false;
}

} // namespace zpk_detail

/**
 * @brief Convert a SISO state-space model to zero-pole-gain form
 *
 * - Poles: eigenvalues of @f$ A @f$
 * - Biproper (@f$ |D| > tol @f$): @f$ k = D @f$, zeros = eig(@f$ A - B D^{-1} C @f$)
 * - Strictly proper: zeros and @f$ k @f$ from the Leverrier numerator (finite
 *   zeros only; @c n_zeros ≤ NX−1). Effective zero degrees above 12 require
 *   @c n_zeros == NX (full-order path) or use the TF→ZPK path with an explicit
 *   TransferFunction of the correct size
 *
 * @note Compare with MATLAB®'s zpk(sys) for SISO ss.
 * @param sys SISO plant (NU = NY = 1)
 * @param tol Magnitude floor for D and poly leads
 * @return StateSpaceZPKResult
 */
template<size_t NX, typename T = double, size_t NW = 0, size_t NV = 0>
    requires std::is_floating_point_v<T>
[[nodiscard]] constexpr StateSpaceZPKResult<NX, T>
to_zpk(const StateSpace<NX, 1, 1, T, NW, NV>& sys, T tol = default_tol<T>()) {
    StateSpaceZPKResult<NX, T> art{};

    if constexpr (NX == 0) {
        art.gain = sys.D(0, 0);
        art.n_zeros = 0;
        art.success = true;
        return art;
    }

    // Drop noise dimensions: work on the core A,B,C,D
    StateSpace<NX, 1, 1, T> core{
        .A = sys.A,
        .B = sys.B,
        .C = sys.C,
        .D = sys.D,
        .Ts = sys.Ts
    };

    const auto pe = mat::compute_eigenvalues(core.A);
    if (!pe.converged) {
        return art;
    }
    for (size_t i = 0; i < NX; ++i) {
        art.poles[i] = pe.values[i];
    }

    const T D0 = core.D(0, 0);
    if (damp::abs(D0) > tol) {
        Matrix<NX, NX, T> Abar = core.A;
        const T           invD = T{1} / D0;
        for (size_t i = 0; i < NX; ++i) {
            for (size_t j = 0; j < NX; ++j) {
                Abar(i, j) -= core.B(i, 0) * invD * core.C(0, j);
            }
        }
        const auto ze = mat::compute_eigenvalues(Abar);
        if (!ze.converged) {
            return art;
        }
        for (size_t i = 0; i < NX; ++i) {
            art.zeros[i] = ze.values[i];
        }
        art.n_zeros = NX;
        art.gain = D0;
        art.success = true;
        return art;
    }

    const auto tf = zpk_detail::ss_to_tf_full(core);
    if (!(damp::abs(tf.den[NX]) > tol)) {
        return art;
    }

    size_t deg_num = 0;
    bool   found = false;
    for (size_t k = NX; k > 0; --k) {
        if (damp::abs(tf.num[k]) > tol) {
            deg_num = k;
            found = true;
            break;
        }
    }
    if (!found) {
        art.n_zeros = 0;
        art.gain = tf.num[0] / tf.den[NX];
        art.success = true;
        return art;
    }

    art.gain = tf.num[deg_num] / tf.den[NX];
    art.n_zeros = deg_num;
    art.success = zpk_detail::fill_zeros_dispatch(deg_num, tf.num, art, tol);
    return art;
}

/**
 * @brief MATLAB®-style alias for SISO to_zpk (state-space)
 * @note Compare with MATLAB®'s ss2zpk / zpk(ss).
 */
template<size_t NX, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr StateSpaceZPKResult<NX, T>
ss2zpk(const StateSpace<NX, 1, 1, T, NW, NV>& sys, T tol = default_tol<T>()) {
    return to_zpk(sys, tol);
}

/**
 * @brief SISO state-space → transfer function (Leverrier), size NX+1 / NX+1
 *
 * Monic denominator. Leading numerator coefficient is @c D (zero when strictly proper).
 *
 * @note Compare with MATLAB®'s tf(ss_sys) / ss2tf.
 */
template<size_t NX, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr TransferFunction<NX + 1, NX + 1, T>
to_transfer_function(const StateSpace<NX, 1, 1, T, NW, NV>& sys) {
    StateSpace<NX, 1, 1, T> core{
        .A = sys.A,
        .B = sys.B,
        .C = sys.C,
        .D = sys.D,
        .Ts = sys.Ts
    };
    return zpk_detail::ss_to_tf_full(core);
}

/**
 * @brief MATLAB®-style alias for SISO to_transfer_function
 * @note Compare with MATLAB®'s ss2tf.
 */
template<size_t NX, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr TransferFunction<NX + 1, NX + 1, T>
ss2tf(const StateSpace<NX, 1, 1, T, NW, NV>& sys) {
    return to_transfer_function(sys);
}

// ============================================================================
// design:: thin aliases (design-time conversion / minreal)
// ============================================================================

namespace design {

using damp::poly_roots;
using damp::ss2tf;
using damp::ss2zpk;
using damp::tf2zpk;
using damp::to_transfer_function;
using damp::to_zpk;
using damp::zpk2tf;

/**
 * @brief Cancel matching pole-zero pairs on a ZPK model
 * @note Compare with MATLAB®'s minreal(zpk_sys) in factored form.
 */
template<size_t Nz, size_t Np, typename T = double>
[[nodiscard]] constexpr typename ZPK<Nz, Np, T>::CancelResult
minreal_zpk(const ZPK<Nz, Np, T>& sys, T tol = default_tol<T>()) {
    return sys.cancel_matching(tol);
}

} // namespace design

} // namespace damp
