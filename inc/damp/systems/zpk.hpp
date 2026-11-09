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
 * | ZPK → TF  | ZPK::to_transfer_function | Always succeeds; expands factored form |
 * | ZPK → SS  | ZPK::to_state_space | Via TF companion realization |
 * | TF → ZPK  | to_zpk | Poly roots via companion + eigen; @c success if leads ≠ 0 and QR converges |
 * | SS → ZPK  | to_zpk (SISO) | Zeros: transmission_zeros (QZ on the reduced system pencil). @f$k@f$: the Markov parameter at the relative degree |
 * | SS → TF   | to_transfer_function / ss2tf | Expands the SS → ZPK factors (SISO), numerator checked against the Markov parameters; Faddeev–LeVerrier only if QR / QZ does not converge |
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
 * numerator/denominator coefficients must be nonzero for TF→ZPK. SS→ZPK
 * returns the finite zeros in StateSpaceZPKResult (runtime @c n_zeros); call
 * StateSpaceZPKResult::as_zpk when the count is known at compile time. Its zeros
 * come from the system pencil, not from the transfer-function numerator, so they
 * do not carry the error of polynomial coefficients. High-order root finding
 * (degree ≫ 8–12) for TF→ZPK is the same QR path as the rest of the library:
 * solid for typical plant orders, not a multiprecision polynomial solver.
 *
 * Example: round-trip a lead compensator
 * @code
 * #include "damp/systems/zpk.hpp"
 * using namespace damp;
 *
 * constexpr TransferFunction<2, 2> tf{.num = {1.0, 2.0}, .den = {1.0, 10.0}}; // (s+2)/(s+10)
 * constexpr auto art = to_zpk(tf);
 * static_assert(art.success);
 * constexpr auto tf2 = art.zpk.to_transfer_function();
 * @endcode
 *
 * @note Compare with MATLAB®'s zpk, zpkdata, tf2zpk, zpk2tf, ss2tf, minreal (factored form).
 */


#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/eigen.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/matrix_traits.hpp"
#include "state_space.hpp"
#include "transfer_function.hpp"
#include "zeros.hpp"

namespace damp {

// ============================================================================
// Fixed-size polynomial roots (embeddable)
// ============================================================================

/**
 * @brief Roots of a real polynomial given highest power first
 *
 * For coefficients @f$ c_0 x^{N-1} + c_1 x^{N-2} + \cdots + c_{N-1} @f$, returns the
 * @f$ N-1 @f$ roots as eigenvalues of the companion matrix. Requires a nonzero
 * leading coefficient @f$ c_0 @f$ (only an exact zero is rejected: the
 * coefficients scale as powers of the variable's unit, so no magnitude floor
 * on them is unit-invariant).
 *
 * - Degree 0 (@c N == 1): empty root set, @c success = true
 * - Degree 1–2: closed form
 * - Higher degree: Francis QR via mat::compute_eigenvalues; @c success
 *   mirrors the eigen convergence flag
 *
 * @note Compare with MATLAB®'s roots(c).
 * @see to_zpk() for the TF factorization path that uses this
 *
 * @tparam N Number of coefficients (degree N−1 when N ≥ 1)
 * @tparam T Floating-point scalar type
 */
template<size_t N, typename T = double>
struct [[nodiscard]] PolyRootsResult {
    damp::array<damp::complex<T>, (N > 0 ? N - 1 : 0)> roots{};        ///< Roots (empty when N ≤ 1)
    bool                                               success{false}; ///< true if leading coeff OK and eigen converged

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

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
 * @brief Compute roots of a real polynomial given highest power first
 *
 * @param c Coefficients c₀ x^{N−1} + c₁ x^{N−2} + … + c_{N−1}
 * @return @ref PolyRootsResult; @c success = false if c₀ = 0
 */
template<size_t N, typename T = double>
    requires damp::is_floating_point_v<T>
[[nodiscard]] constexpr PolyRootsResult<N, T> poly_roots(const damp::array<T, N>& c) {
    if constexpr (N <= 1) {
        return PolyRootsResult<N, T>{.success = true};
    } else if constexpr (N == 2) {
        if (c[0] == T{0}) {
            return PolyRootsResult<N, T>{};
        }
        return PolyRootsResult<N, T>{
            .roots = {damp::complex<T>{-c[1] / c[0], T{0}}},
            .success = true,
        };
    } else {
        if (c[0] == T{0}) {
            return PolyRootsResult<N, T>{};
        }
        constexpr size_t M = N - 1;
        Matrix<M, M, T>  companion{};
        for (size_t i = 0; i + 1 < M; ++i) {
            companion(i, i + 1) = T{1};
        }
        // Last row −[c_N−1 … c_1]/c₀: column j multiplies x^j.
        const T inv_lead = T{1} / c[0];
        for (size_t j = 0; j < M; ++j) {
            companion(M - 1, j) = -c[N - 1 - j] * inv_lead;
        }
        const auto eig = mat::compute_eigenvalues(companion);
        return PolyRootsResult<N, T>{.roots = eig.values, .success = eig.success};
    }
}

namespace zpk_detail {

/**
 * @brief Multiply a polynomial of @p degree by (s − @p root) in place; poly[k] multiplies sᵏ.
 */
template<typename Poly, typename C>
constexpr void mul_linear(Poly& poly, size_t degree, C root) {
    poly[degree + 1] = poly[degree];
    for (size_t j = degree; j > 0; --j) {
        poly[j] = poly[j - 1] - root * poly[j];
    }
    poly[0] = -root * poly[0];
}

} // namespace zpk_detail

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
 * A real model has its zeros and its poles each closed under conjugation; the conversions to
 * TF and SS assume it.
 *
 * MATLAB® equivalent: `zpk(z, p, k)`.
 *
 * @tparam Nz   Number of zeros
 * @tparam Np   Number of poles
 * @tparam T    Floating-point scalar type
 */
template<size_t Nz, size_t Np, typename T = double>
    requires damp::is_floating_point_v<T>
struct ZPK {
    damp::array<damp::complex<T>, Nz> zeros{}; ///< Zeros of the transfer function, closed under conjugation
    damp::array<damp::complex<T>, Np> poles{}; ///< Poles of the transfer function, closed under conjugation
    T                                 gain{1}; ///< Gain k of the monic factored form

    /**
     * @brief Result of @ref cancel_matching (factored minreal)
     */
    struct [[nodiscard]] CancelResult {
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
     * Expands @f$ k \prod(s - z_i) / \prod(s - p_j) @f$ into real polynomials,
     * highest power first. Intermediate products use complex coefficients so conjugate
     * pairs cancel imaginary parts; the returned TF stores the real parts. The zeros and the
     * poles must each be closed under conjugation: a complex root without its conjugate has
     * no real polynomial, and the imaginary parts it leaves are dropped.
     *
     * @note Compare with MATLAB®'s tf(zpk_sys) / zpk2tf.
     * @return Transfer function with Nnum = Nz+1, Nden = Np+1
     */
    [[nodiscard]] constexpr TransferFunction<Nz + 1, Np + 1, T> to_transfer_function() const {
        TransferFunction<Nz + 1, Np + 1, T> tf{};

        damp::array<damp::complex<T>, Nz + 1> num_c{};
        num_c[0] = damp::complex<T>{gain, T{0}};
        for (size_t i = 0; i < Nz; ++i) {
            zpk_detail::mul_linear(num_c, i, zeros[i]);
        }
        for (size_t k = 0; k < Nz + 1; ++k) {
            tf.num[k] = num_c[Nz - k].real(); // num_c[k] multiplies s^k
        }

        damp::array<damp::complex<T>, Np + 1> den_c{};
        den_c[0] = damp::complex<T>{T{1}, T{0}};
        for (size_t i = 0; i < Np; ++i) {
            zpk_detail::mul_linear(den_c, i, poles[i]);
        }
        for (size_t k = 0; k < Np + 1; ++k) {
            tf.den[k] = den_c[Np - k].real();
        }

        return tf;
    }

    /**
     * @brief Convert ZPK to state-space (companion realization of the TF)
     *
     * The expanded denominator is monic, so the companion realization always exists.
     *
     * @note Compare with MATLAB®'s ss(zpk_sys).
     * @return Continuous StateSpace (Ts = 0) with NX = Np
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
     * @param tol Pole-coincidence distance |s − pⱼ| (default 0: only an exact hit; an
     *            absolute default would refuse slow-but-nonzero poles)
     * @return H(s), or nullopt if any pole lies within @p tol of @p s
     */
    [[nodiscard]] constexpr damp::optional<damp::complex<T>> eval(damp::complex<T> s, T tol = T{0}) const {
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
    [[nodiscard]] constexpr damp::optional<T> dcgain(T tol = T{0}) const {
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
    [[nodiscard]] constexpr damp::optional<T> dcgain_discrete(T tol = T{0}) const {
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
     * @f$ |z-p| \le \mathrm{tol} \cdot \max(|z|, |p|) + 4 (N_z + N_p)\,\epsilon R @f$,
     * R the largest root magnitude of the model. Both terms scale with the frequency unit, so
     * the test is invariant to it. The first is relative to the pair, so two distinct slow
     * roots are not cancelled just for being small; the second is the roundoff of roots
     * computed from a model of size R, so a computed common factor at the origin (a zero at
     * 1e-17 against a pole at 0) cancels. Survivors are compacted to the front of the arrays;
     * each cancelled value is re-appended as an identical @f$ z = p @f$ pair so
     * to_transfer_function remains algebraically equivalent (common factors),
     * while @c n_cancelled reports how many pairs were cancelled.
     *
     * Gain is unchanged (cancelled monic factors contribute 1).
     *
     * @note Compare with MATLAB®'s minreal on a zpk model (factored form only;
     *       does not perform SS uncontrollable/unobservable reduction).
     * @param tol Relative match tolerance (default √eps). The computed roots of a simple common
     *            factor agree far closer than that; those of a k-fold factor spread by about
     *            ε^(1/k) of their magnitude, so a double or triple factor needs a tol above
     *            that spread (∛eps or more), at the risk of cancelling distinct roots as close)
     * @return @ref CancelResult
     */
    [[nodiscard]] constexpr CancelResult cancel_matching(T tol = damp::sqrt(eps<T>())) const {
        CancelResult out{};
        out.zpk.gain = gain;

        damp::array<bool, Nz>                              z_gone{};
        damp::array<bool, Np>                              p_gone{};
        damp::array<damp::complex<T>, (Nz < Np ? Nz : Np)> cancelled{};
        size_t                                             n_cancel = 0;

        T largest = T{0}; // R: the model's scale, for the roundoff floor of computed roots
        for (size_t i = 0; i < Nz; ++i) {
            largest = damp::max(largest, damp::abs(zeros[i]));
        }
        for (size_t j = 0; j < Np; ++j) {
            largest = damp::max(largest, damp::abs(poles[j]));
        }
        const T roundoff = static_cast<T>(Nz + Np) * T{4} * eps<T>() * largest;

        for (size_t i = 0; i < Nz; ++i) {
            for (size_t j = 0; j < Np; ++j) {
                if (z_gone[i] || p_gone[j]) {
                    continue;
                }
                const T zm = damp::abs(zeros[i]);
                const T pm = damp::abs(poles[j]);
                if (damp::abs(zeros[i] - poles[j]) <= (tol * ((zm > pm) ? zm : pm)) + roundoff) {
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
    [[nodiscard]] constexpr ZPK minreal(T tol = damp::sqrt(eps<T>())) const {
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
struct [[nodiscard]] ZPKResult {
    ZPK<Nz, Np, T> zpk{};          ///< Factored model (valid when success)
    bool           success{false}; ///< true if both poly root solves succeeded

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

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
 * Requires nonzero leading coefficients (a zero-padded array has fewer roots
 * than its fixed size and is rejected; a constant denominator of 0 is refused too). That
 * includes the numerator to_transfer_function(StateSpace) returns for a strictly proper model,
 * padded with leading zeros: to_zpk(StateSpace) is the direct route there. Degree ≥ 3 uses
 * iterative QR (design-time).
 *
 * @note Compare with MATLAB®'s tf2zpk / zpk(tf_sys).
 * @see ZPK::to_transfer_function for the inverse map
 * @see poly_roots
 *
 * @param tf Transfer function (coefficients highest power first)
 * @return ZPKResult with success flag
 */
template<size_t Nnum, size_t Nden, typename T = double>
    requires damp::is_floating_point_v<T>
[[nodiscard]] constexpr ZPKResult<(Nnum > 0 ? Nnum - 1 : 0), (Nden > 0 ? Nden - 1 : 0), T>
to_zpk(const TransferFunction<Nnum, Nden, T>& tf) {
    constexpr size_t     Nz = (Nnum > 0 ? Nnum - 1 : 0);
    constexpr size_t     Np = (Nden > 0 ? Nden - 1 : 0);
    ZPKResult<Nz, Np, T> art{};

    if constexpr (Nden == 0 || Nnum == 0) {
        return art;
    }
    if (tf.den[0] == T{0}) {
        return art; // no leading denominator coefficient (a constant 0 has no gain)
    }

    const auto zr = poly_roots(tf.num);
    const auto pr = poly_roots(tf.den);
    if (!zr.success || !pr.success) {
        return art;
    }

    art.zpk.gain = tf.num[0] / tf.den[0];
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

// ============================================================================
// SS → TF / ZPK (SISO)
// ============================================================================

namespace zpk_detail {

/**
 * @brief Faddeev–LeVerrier characteristic polynomial and adjoint coefficient matrices
 *
 * Builds monic @f$ \det(sI-A) @f$ (@c den, highest power first) and @f$ M_1,\ldots,M_n @f$
 * with @f$ M_1 = I @f$ such that
 * @f$ \mathrm{adj}(sI-A) = \sum_{j=0}^{n-1} s^{n-1-j} M_{j+1} @f$.
 *
 * @see Faddeev–LeVerrier; Golub & Van Loan, Matrix Computations
 */
template<size_t NX, typename T>
struct [[nodiscard]] LeverrierResult {
    damp::array<T, NX + 1>             den{}; ///< Monic char poly, highest power first
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
        r.den[0] = T{1};

        for (size_t k = 1; k <= NX; ++k) {
            Mk = A * Mk + c_prev * I;
            r.M[k - 1] = Mk;
            const T ck = -(A * Mk).trace() / static_cast<T>(k);
            r.den[k] = ck;
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
            tf.num[j + 1] += cMb; // M_{j+1} multiplies s^{NX−1−j}
        }
        return tf;
    }
}

/**
 * @brief Magnitude scale of each strictly proper Leverrier numerator coefficient
 *
 * scale[j+1] = Σᵣ |C|ᵣ (|M_{j+1}| |B|)ᵣ — the size of the terms whose sum forms
 * the coefficient, so a coefficient far below it is cancellation residue.
 */
template<size_t NX, typename T>
[[nodiscard]] constexpr damp::array<T, NX + 1> leverrier_num_scale(const StateSpace<NX, 1, 1, T>& sys) {
    damp::array<T, NX + 1> scale{};
    if constexpr (NX > 0) {
        const auto lev = leverrier(sys.A);
        for (size_t j = 0; j < NX; ++j) {
            T s = T{0};
            for (size_t r = 0; r < NX; ++r) {
                T row = T{0};
                for (size_t c = 0; c < NX; ++c) {
                    row += damp::abs(lev.M[j](r, c)) * damp::abs(sys.B(c, 0));
                }
                s += damp::abs(sys.C(0, r)) * row;
            }
            scale[j + 1] = s;
        }
    }
    return scale;
}

/**
 * @brief @p num, or the numerator the Markov parameters give over @p den when @p num misses it
 *
 * Any realization satisfies num[k] = Σ_{j≤k} den[j]·h_{k−j}, h₀ = D, hᵢ = C·A^{i−1}·B. The sum's
 * own roundoff is about (NX+1)·eps·sₖ, with sₖ = Σ_{j≤k} |den[j]|·|h_{k−j}|ₛ the size of its terms
 * (|hᵢ|ₛ = |C|·|A|^{i−1}·|B| entrywise, |h₀|ₛ = |D|); a @p num that misses it by more than
 * 8·(NX+1)·eps·sₖ is the less accurate of the two, and the sum replaces it. Zeros of a badly scaled pencil (a Tustin realization in float) miss it; the
 * factors of a model whose poles spread over decades pass, as sₖ is large there. Markov
 * parameters past the range of T keep @p num.
 */
template<size_t NX, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr damp::array<T, NX + 1> markov_checked_numerator(
    const StateSpace<NX, 1, 1, T, NW, NV>& sys,
    const damp::array<T, NX + 1>&          den,
    const damp::array<T, NX + 1>&          num
) {
    const auto             mag = [](T v) { return damp::abs(v); };
    const auto             A_abs = mat::transform(sys.A, mag);
    const auto             C_abs = mat::transform(sys.C, mag);
    Matrix<NX, 1, T>       v = sys.B;                           // A^{i−1}·B
    Matrix<NX, 1, T>       v_size = mat::transform(sys.B, mag); // |A|^{i−1}·|B|
    damp::array<T, NX + 1> h{sys.D(0, 0)};
    damp::array<T, NX + 1> h_size{damp::abs(sys.D(0, 0))};
    for (size_t i = 1; i <= NX; ++i) {
        h[i] = (sys.C * v)(0, 0);
        h_size[i] = (C_abs * v_size)(0, 0);
        v = sys.A * v;
        v_size = A_abs * v_size;
    }
    damp::array<T, NX + 1> markov{};
    bool                   missed = false;
    for (size_t k = 0; k <= NX; ++k) {
        T size = T{0};
        for (size_t j = 0; j <= k; ++j) {
            markov[k] += den[j] * h[k - j];
            size += damp::abs(den[j]) * h_size[k - j];
        }
        const T bound = static_cast<T>(NX + 1) * T{8} * eps<T>() * size;
        if (!damp::isfinite(bound) || !damp::isfinite(markov[k])) {
            return num;
        }
        if (damp::abs(num[k] - markov[k]) > bound) {
            missed = true;
        }
    }
    return missed ? markov : num;
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
struct [[nodiscard]] StateSpaceZPKResult {
    damp::array<damp::complex<T>, NX> poles{};        ///< Eigenvalues of A
    damp::array<damp::complex<T>, NX> zeros{};        ///< First @ref n_zeros entries valid
    size_t                            n_zeros{0};     ///< Number of finite zeros
    T                                 gain{T{1}};     ///< Gain k of the monic factored form
    bool                              success{false}; ///< true if the eigen/root solves converged

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

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

/**
 * @brief Convert a SISO state-space model to zero-pole-gain form
 *
 * - Poles: eigenvalues of @f$ A @f$
 * - Zeros: transmission_zeros, the finite zeros of the system pencil (QZ on the pencil reduced
 *   by orthogonal transformations; no polynomial is formed and D is not inverted). A D at
 *   roundoff level counts as zero, and the zeros it would put far out as infinite
 * - Gain: the leading numerator coefficient, the Markov parameter at the relative degree
 *   r = NX − n_zeros: @f$ k = D @f$ for r = 0, @f$ k = C A^{r-1} B @f$ otherwise
 *
 * The rank decisions are relative to the model scaled by powers of two, so @c n_zeros does not
 * depend on the units of time, of the input, of the output or of the states.
 *
 * @note Compare with MATLAB®'s zpk(sys) for SISO ss.
 * @param sys SISO plant (NU = NY = 1)
 * @return StateSpaceZPKResult
 */
template<size_t NX, typename T = double, size_t NW = 0, size_t NV = 0>
    requires damp::is_floating_point_v<T>
[[nodiscard]] constexpr StateSpaceZPKResult<NX, T> to_zpk(const StateSpace<NX, 1, 1, T, NW, NV>& sys) {
    StateSpaceZPKResult<NX, T> art{};

    if constexpr (NX == 0) {
        art.gain = sys.D(0, 0);
        art.n_zeros = 0;
        art.success = true;
        return art;
    } else {
        const auto pe = mat::compute_eigenvalues(sys.A);
        if (!pe.success) {
            return art;
        }
        art.poles = pe.values;
        const auto tz = transmission_zeros(sys);
        if (!tz.success) {
            return art;
        }
        art.zeros = tz.zeros;
        art.n_zeros = tz.n_zeros;
        const size_t relative_degree = NX - tz.n_zeros;
        if (relative_degree == 0) {
            art.gain = sys.D(0, 0);
        } else {
            Matrix<1, NX, T> row = sys.C; // C·A^(r−1)
            for (size_t k = 1; k < relative_degree; ++k) {
                row = row * sys.A;
            }
            art.gain = (row * sys.B)(0, 0);
        }
        art.success = true;
        return art;
    }
}

/**
 * @brief SISO state-space → transfer function, size NX+1 / NX+1
 *
 * Monic denominator, highest power first. The leading numerator coefficient, num[0], is @c D (zero
 * when strictly proper, so the numerator is padded with leading zeros, which to_zpk(TransferFunction)
 * refuses: to_zpk(StateSpace) is the direct route to the factors).
 *
 * The polynomials are the expanded factors of to_zpk(StateSpace): den = ∏(s − pᵢ) over the
 * eigenvalues of A, num = k·∏(s − zᵢ) over the transmission zeros, k the Markov parameter at
 * the relative degree. Expanding factors loses no digits to cancellation when the poles spread
 * over decades, which the Faddeev–LeVerrier trace recursion does (its coefficients come from
 * cancelling terms that grow like |λ_max|ᵏ). The numerator is checked against the Markov
 * parameters, num[k] = Σ_{j≤k} den[j]·h_{k−j} with h₀ = D and hᵢ = C·A^{i−1}·B: transmission
 * zeros of a badly scaled pencil (a Tustin realization in float, B ~ 1e−9 against C ~ 1e6) can
 * move by percents, and when the expanded factors miss that sum by more than its roundoff the
 * sum is the numerator. Faddeev–LeVerrier is the fallback only when the eigenvalue or zero
 * computation does not converge.
 *
 * @note Compare with MATLAB®'s tf(ss_sys) / ss2tf.
 * @see "Computation of zeros of linear multivariable systems" (Emami-Naeini & Van Dooren, 1982)
 */
template<size_t NX, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr TransferFunction<NX + 1, NX + 1, T>
to_transfer_function(const StateSpace<NX, 1, 1, T, NW, NV>& sys) {
    const StateSpaceZPKResult<NX, T> factors = to_zpk(sys);
    if (!factors.success) {
        return zpk_detail::ss_to_tf_full(StateSpace<NX, 1, 1, T>{.A = sys.A, .B = sys.B, .C = sys.C, .D = sys.D, .Ts = sys.Ts});
    }
    damp::array<damp::complex<T>, NX + 1> num_c{}; // num_c[k] multiplies sᵏ
    damp::array<damp::complex<T>, NX + 1> den_c{};
    num_c[0] = damp::complex<T>{factors.gain, T{0}};
    den_c[0] = damp::complex<T>{T{1}, T{0}};
    for (size_t i = 0; i < factors.n_zeros; ++i) {
        zpk_detail::mul_linear(num_c, i, factors.zeros[i]);
    }
    for (size_t i = 0; i < NX; ++i) {
        zpk_detail::mul_linear(den_c, i, factors.poles[i]);
    }
    TransferFunction<NX + 1, NX + 1, T> tf{};
    for (size_t k = 0; k <= NX; ++k) {
        tf.den[k] = den_c[NX - k].real();
        tf.num[k] = (NX - k <= factors.n_zeros) ? num_c[NX - k].real() : T{0}; // right-aligned on the constant term
    }
    tf.num = zpk_detail::markov_checked_numerator(sys, tf.den, tf.num);
    return tf;
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
using damp::to_transfer_function;
using damp::to_zpk;

/**
 * @brief Cancel matching pole-zero pairs on a ZPK model
 * @note Compare with MATLAB®'s minreal(zpk_sys) in factored form.
 */
template<size_t Nz, size_t Np, typename T = double>
[[nodiscard]] constexpr typename ZPK<Nz, Np, T>::CancelResult
minreal_zpk(const ZPK<Nz, Np, T>& sys, T tol = damp::sqrt(eps<T>())) {
    return sys.cancel_matching(tol);
}

} // namespace design

} // namespace damp
