// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file functions.hpp
 * @brief Matrix norms, determinant and rank, and matrix functions (expm, sqrt, log,
 *        pow, sin/cos, sinh/cosh)
 */


#include "core.hpp"
#include "damp/backend.hpp"
#include "damp/matrix/matrix_traits.hpp"
#include "eigen.hpp" // sqrt() uses the real Schur form
#include "solve.hpp" // expm() and log() use mat::solve
#include "svd.hpp"   // two_norm() and rank() use the SVD

namespace damp {

namespace mat {
/**
 * @brief Infinity norm ‖A‖∞: maximum absolute row sum
 *
 * Always returns a real scalar (@c scalar_type_t&lt;T&gt;), including for complex
 * element types where each entry contributes @f$ |a_{ij}| @f$.
 *
 * @note Compare with MATLAB®'s norm(A, inf).
 *
 * @tparam T Element type (real or damp::complex)
 * @tparam N Matrix dimension
 * @param A Square matrix
 * @return Infinity norm of the matrix
 */
template<typename T, size_t N>
[[nodiscard]] constexpr scalar_type_t<T> infinity_norm(const Matrix<N, N, T>& A) {
    using real_t = scalar_type_t<T>;
    real_t norm = real_t{0};
    for (size_t i = 0; i < N; ++i) {
        real_t row_sum = real_t{0};
        for (size_t j = 0; j < N; ++j) {
            row_sum += damp::abs(A(i, j));
        }
        if (row_sum > norm) {
            norm = row_sum;
        }
    }
    return norm;
}

/**
 * @brief One-norm ‖A‖₁: maximum absolute column sum
 *
 * Always returns a real scalar (@c scalar_type_t&lt;T&gt;), including for complex
 * element types where each entry contributes @f$ |a_{ij}| @f$.
 *
 * @note Compare with MATLAB®'s norm(A, 1).
 *
 * @tparam T Element type (real or damp::complex)
 * @tparam N Matrix dimension
 * @param A Square matrix
 * @return One-norm of the matrix
 */
template<typename T, size_t N>
[[nodiscard]] constexpr scalar_type_t<T> one_norm(const Matrix<N, N, T>& A) {
    using real_t = scalar_type_t<T>;
    real_t norm = real_t{0};
    for (size_t j = 0; j < N; ++j) {
        real_t col_sum = real_t{0};
        for (size_t i = 0; i < N; ++i) {
            col_sum += damp::abs(A(i, j));
        }
        if (col_sum > norm) {
            norm = col_sum;
        }
    }
    return norm;
}

/**
 * @brief Spectral norm ‖A‖₂ = σₘₐₓ(A)
 *
 * @f[
 *   \|A\|_2 = \sigma_{\max}(A) = \sqrt{\lambda_{\max}(A^{\mathrm{H}} A)}
 * @f]
 * where @f$ A^{\mathrm{H}} @f$ is the conjugate transpose.
 *
 * Computed as the largest singular value from mat::svd. Power iteration from a
 * fixed start vector fails when that vector is orthogonal to the dominant
 * singular vector (e.g. returns 0 for [[1, −1], [0, 0]]). Always returns a real
 * scalar.
 *
 * @note Compare with MATLAB®'s norm(A, 2).
 *
 * @tparam T Element type (real or damp::complex)
 * @tparam N Matrix dimension
 * @param A Square matrix
 * @return Spectral norm
 */
template<typename T, size_t N>
[[nodiscard]] constexpr scalar_type_t<T> two_norm(const Matrix<N, N, T>& A) {
    const auto s = svd(A);
    return s.success ? s.singular_values[0] : damp::numeric_limits<scalar_type_t<T>>::quiet_NaN();
}

namespace detail {

/**
 * @brief Gaussian elimination with partial pivoting on the leading n×n block of @p M, carrying @p b.
 *
 * On return @p M holds U in its leading block and @p b holds L⁻¹Pb, so U x = b solves M x = b.
 * Returns det = sign(P)·Π uᵢᵢ, or 0 when a pivot column is exactly zero. Unlike
 * lu_decomposition() there is no cancellation floor: an ill-conditioned but nonsingular
 * matrix keeps its (small) determinant.
 */
template<typename T, size_t K>
[[nodiscard]] constexpr T pivot_eliminate(Matrix<K, K, T>& M, ColVec<K, T>& b, size_t n) {
    T d = T{1};
    for (size_t i = 0; i < n; ++i) {
        size_t p = i;
        for (size_t r = i + 1; r < n; ++r) {
            if (damp::abs(M(r, i)) > damp::abs(M(p, i))) {
                p = r;
            }
        }
        if (!(damp::abs(M(p, i)) > scalar_type_t<T>{0})) {
            return T{0};
        }
        if (p != i) {
            for (size_t c = 0; c < n; ++c) {
                damp::swap(M(i, c), M(p, c));
            }
            damp::swap(b(i), b(p));
            d = -d;
        }
        d *= M(i, i);
        for (size_t r = i + 1; r < n; ++r) {
            const T f = M(r, i) / M(i, i);
            M(r, i) = T{0};
            for (size_t c = i + 1; c < n; ++c) {
                M(r, c) -= f * M(i, c);
            }
            b(r) -= f * b(i);
        }
    }
    return d;
}

} // namespace detail

/**
 * @brief Matrix determinant det(A)
 *
 * Cofactor expansion for @f$ N \le 4 @f$; Gaussian elimination with partial pivoting for
 * larger @f$ N @f$ (@f$ \det(A) = \mathrm{sign}(P)\prod_i U_{ii} @f$), which returns 0 only on
 * an exactly zero pivot column, so a nearly singular A keeps its small determinant.
 *
 * @note Compare with MATLAB®'s det(A).
 *
 * @tparam T Element type
 * @tparam N Matrix dimension
 * @param A Square matrix
 * @return Determinant of the matrix
 */
template<typename T, size_t N>
[[nodiscard]] constexpr T det(const Matrix<N, N, T>& A) {
    if constexpr (N == 1) {
        return A(0, 0);
    } else if constexpr (N == 2) {
        return A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
    } else if constexpr (N == 3) {
        return A(0, 0) * (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1))
             - A(0, 1) * (A(1, 0) * A(2, 2) - A(1, 2) * A(2, 0))
             + A(0, 2) * (A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0));
    } else if constexpr (N == 4) {
        T d = A(0, 0) * det(Matrix<3, 3, T>{
                  {A(1, 1), A(1, 2), A(1, 3)},
                  {A(2, 1), A(2, 2), A(2, 3)},
                  {A(3, 1), A(3, 2), A(3, 3)},
              });

        d -= A(0, 1) * det(Matrix<3, 3, T>{
                 {A(1, 0), A(1, 2), A(1, 3)},
                 {A(2, 0), A(2, 2), A(2, 3)},
                 {A(3, 0), A(3, 2), A(3, 3)},
             });

        d += A(0, 2) * det(Matrix<3, 3, T>{
                 {A(1, 0), A(1, 1), A(1, 3)},
                 {A(2, 0), A(2, 1), A(2, 3)},
                 {A(3, 0), A(3, 1), A(3, 3)},
             });

        d -= A(0, 3) * det(Matrix<3, 3, T>{
                 {A(1, 0), A(1, 1), A(1, 2)},
                 {A(2, 0), A(2, 1), A(2, 2)},
                 {A(3, 0), A(3, 1), A(3, 2)},
             });

        return d;
    } else {
        Matrix<N, N, T> U = A;
        ColVec<N, T>    unused{};
        return detail::pivot_eliminate(U, unused, N);
    }
}

/**
 * @brief Numerical rank: the number of singular values above rel_tol·σₘₐₓ
 *
 * Uses the SVD, the only rank-revealing factorization here; Gaussian
 * elimination with partial pivoting can miss a near-dependency. The threshold
 * is relative to @f$ \sigma_{\max} @f$, so scaling @p M never changes the result.
 * The default @f$ \max(R, C)\,\varepsilon @f$ matches MATLAB®. Works for
 * rectangular and complex matrices. Single source for the structural rank
 * tests in stability.hpp and riccati.hpp.
 *
 * @note Compare with MATLAB®'s rank(A).
 *
 * @tparam R   Rows
 * @tparam C   Columns
 * @tparam T   Element type (float/double or damp::complex thereof)
 * @param M        Input matrix
 * @param rel_tol  Singular values ≤ rel_tol·σₘₐₓ count as zero (default max(R, C)·eps)
 * @return Rank of the matrix (number of linearly independent rows/columns)
 */
template<size_t R, size_t C, typename T>
[[nodiscard]] constexpr size_t rank(
    const Matrix<R, C, T>& M,
    scalar_type_t<T>       rel_tol = detail::svd_default_rel_tol<R, C, T>()
) {
    return rank_from_svd(svd(M), rel_tol);
}

namespace detail {

// Al-Mohy & Higham 2009 scaling-and-squaring. Degree and the power-of-two
// scale come from ||A^k||_1^{1/k}. Triangular arguments rewrite the diagonal
// (and the first off-diagonal) from the scalar exponential while squaring
// (Higham 2008, §10.3, eq. 10.42; paper Code Fragment 2.1).

enum class ExpmTriangle { none,
                          upper,
                          lower,
                          diagonal };

struct ExpmScale {
    int    degree{0};
    size_t s{0};
};

template<typename T>
[[nodiscard]] constexpr T expm_log2(T x) {
    return damp::log(x) * damp::numbers::log2e_v<T>;
}

// Unit roundoff u = 2^{-digits}: 2^{-53} for double, 2^{-24} for float.
template<typename T>
[[nodiscard]] constexpr T expm_log2_u() {
    return -static_cast<T>(damp::numeric_limits<T>::digits);
}

// log2(c_m) for the coefficients in Higham 2005, (2.2) and (2.6). The 2009
// algorithm's ell() bound divides by these.
template<typename T>
[[nodiscard]] constexpr T expm_log2_coeff(int m) {
    if (m == 3) {
        return static_cast<T>(16.621136113274641);
    }
    if (m == 5) {
        return static_cast<T>(33.227772656854164);
    }
    if (m == 7) {
        return static_cast<T>(51.994974307382165);
    }
    if (m == 9) {
        return static_cast<T>(72.324718098934952);
    }
    return static_cast<T>(116.447004251763005);
}

template<typename T>
[[nodiscard]] constexpr T expm_entry_exp(T x) {
    if constexpr (is_complex_v<T>) {
        const auto e = damp::exp(x.real());
        return T{e * damp::cos(x.imag()), e * damp::sin(x.imag())};
    } else {
        return damp::exp(x);
    }
}

template<typename T, size_t N>
[[nodiscard]] constexpr Matrix<N, N, T> expm_div_pow2(Matrix<N, N, T> A, size_t s) {
    const T half = static_cast<T>(0.5);
    for (size_t i = 0; i < s; ++i) {
        A = A * half;
    }
    return A;
}

template<typename T>
[[nodiscard]] constexpr T expm_pow2(size_t e) {
    T v = T{1};
    for (size_t i = 0; i < e; ++i) {
        if (!damp::isfinite(v * T{2})) {
            return damp::numeric_limits<T>::max();
        }
        v *= T{2};
    }
    return v;
}

template<typename T, size_t N>
[[nodiscard]] constexpr Matrix<N, N, T> expm_pow(Matrix<N, N, T> base, int exponent) {
    Matrix<N, N, T> result = Matrix<N, N, T>::identity();
    unsigned        expn = static_cast<unsigned>(exponent);
    while (expn > 0U) {
        if ((expn & 1U) != 0U) {
            result = result * base;
        }
        base = base * base;
        expn >>= 1U;
    }
    return result;
}

template<typename T, size_t N>
[[nodiscard]] constexpr Matrix<N, N, scalar_type_t<T>> expm_abs(const Matrix<N, N, T>& A) {
    using real_t = scalar_type_t<T>;
    Matrix<N, N, real_t> out = Matrix<N, N, real_t>::zeros();
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            out(i, j) = damp::abs(A(i, j));
        }
    }
    return out;
}

// || |A|^{2m+1} ||_1 is formed on |A| / 2^e so the product cannot overflow.
// The power of two is added back in the log2 domain.
template<typename T, size_t N>
[[nodiscard]] constexpr size_t expm_ell(const Matrix<N, N, T>& A, int m) {
    using real_t = scalar_type_t<T>;
    const Matrix<N, N, real_t> Abs = expm_abs(A);
    const real_t               nrm = one_norm(Abs);
    if (!(nrm > real_t{0}) || !damp::isfinite(nrm)) {
        return 0;
    }

    size_t e = 0;
    real_t scaled = nrm;
    while (scaled > real_t{1}) {
        scaled *= real_t{0.5};
        ++e;
    }
    const Matrix<N, N, real_t> B = expm_div_pow2(Abs, e);
    const int                  p = (2 * m) + 1;
    const real_t               np = one_norm(expm_pow(B, p));
    if (!(np > real_t{0}) || !damp::isfinite(np)) {
        return 0;
    }

    const real_t two_m = static_cast<real_t>(2 * m);
    const real_t acc = (expm_log2(np) + (static_cast<real_t>(e) * static_cast<real_t>(p)) - expm_log2(nrm) - expm_log2_coeff<real_t>(m) - expm_log2_u<real_t>()) / two_m;
    if (!damp::isfinite(acc) || acc <= real_t{0}) {
        return 0;
    }
    const real_t up = damp::ceil(acc);
    if (up > static_cast<real_t>(1024)) {
        return 1024;
    }
    return static_cast<size_t>(up);
}

template<typename T>
[[nodiscard]] constexpr T expm_root(T nrm, int k) {
    if (!damp::isfinite(nrm) || !(nrm > T{0})) {
        return T{0};
    }
    return damp::pow(nrm, T{1} / static_cast<T>(k));
}

// ||A^k||_1^{1/k}, with A scaled to unit 1-norm before the power.
template<typename T, size_t N>
struct ExpmRoots {
    scalar_type_t<T> d4{};
    scalar_type_t<T> d6{};
    scalar_type_t<T> d8{};
    scalar_type_t<T> d10{};
};

template<typename T, size_t N>
[[nodiscard]] constexpr ExpmRoots<T, N> expm_power_roots(const Matrix<N, N, T>& A) {
    using real_t = scalar_type_t<T>;
    ExpmRoots<T, N> roots{real_t{0}, real_t{0}, real_t{0}, real_t{0}};
    const real_t    nrm = one_norm(A);
    if (!(nrm > real_t{0}) || !damp::isfinite(nrm)) {
        return roots;
    }

    size_t e = 0;
    real_t scaled = nrm;
    while (scaled > real_t{1}) {
        scaled *= real_t{0.5};
        ++e;
    }
    const Matrix<N, N, T> B = expm_div_pow2(A, e);
    const Matrix<N, N, T> B2 = B * B;
    const Matrix<N, N, T> B4 = B2 * B2;
    const Matrix<N, N, T> B6 = B4 * B2;
    const Matrix<N, N, T> B8 = B6 * B2;
    const Matrix<N, N, T> B10 = B8 * B2;
    const real_t          two_e = expm_pow2<real_t>(e);
    const real_t          cap = damp::numeric_limits<real_t>::max();
    const auto            lift = [two_e, cap](const Matrix<N, N, T>& M, int k) {
        const real_t root = expm_root(one_norm(M), k) * two_e;
        return damp::isfinite(root) ? root : cap;
    };
    roots.d4 = lift(B4, 4);
    roots.d6 = lift(B6, 6);
    roots.d8 = lift(B8, 8);
    roots.d10 = lift(B10, 10);
    return roots;
}

template<typename T, size_t N>
[[nodiscard]] constexpr ExpmScale expm_scale(const Matrix<N, N, T>& A) {
    using real_t = scalar_type_t<T>;
    // Θ_m from Al-Mohy & Higham 2009 (the same table SciPy uses for float and double).
    constexpr real_t theta3 = static_cast<real_t>(1.495585217958292e-2);
    constexpr real_t theta5 = static_cast<real_t>(2.539398330063230e-1);
    constexpr real_t theta7 = static_cast<real_t>(9.504178996162932e-1);
    constexpr real_t theta9 = static_cast<real_t>(2.097847961257068e+0);
    constexpr real_t theta13 = static_cast<real_t>(4.25);

    const ExpmRoots<T, N> roots = expm_power_roots(A);
    const real_t          eta1 = damp::max(roots.d4, roots.d6);
    if (eta1 < theta3 && expm_ell(A, 3) == 0) {
        return {3, 0};
    }
    if (eta1 < theta5 && expm_ell(A, 5) == 0) {
        return {5, 0};
    }
    const real_t eta3 = damp::max(roots.d6, roots.d8);
    if (eta3 < theta7 && expm_ell(A, 7) == 0) {
        return {7, 0};
    }
    if (eta3 < theta9 && expm_ell(A, 9) == 0) {
        return {9, 0};
    }

    const real_t eta5 = damp::min(eta3, damp::max(roots.d8, roots.d10));
    size_t       s = 0;
    if (eta5 > theta13) { // finite: expm_power_roots caps each root at max()
        const real_t lg = expm_log2(eta5 / theta13);
        if (lg > real_t{0}) {
            s = static_cast<size_t>(damp::ceil(lg));
        }
    }
    s += expm_ell(expm_div_pow2(A, s), 13);
    if (s > 2048U) {
        s = 2048U;
    }
    return {13, s};
}

template<typename T>
[[nodiscard]] constexpr damp::array<T, 10> expm_pade_coeff(int m) {
    damp::array<T, 10> b{};
    if (m == 3) {
        b[0] = static_cast<T>(120.0);
        b[1] = static_cast<T>(60.0);
        b[2] = static_cast<T>(12.0);
        b[3] = T{1};
    } else if (m == 5) {
        b[0] = static_cast<T>(30240.0);
        b[1] = static_cast<T>(15120.0);
        b[2] = static_cast<T>(3360.0);
        b[3] = static_cast<T>(420.0);
        b[4] = static_cast<T>(30.0);
        b[5] = T{1};
    } else if (m == 7) {
        b[0] = static_cast<T>(17297280.0);
        b[1] = static_cast<T>(8648640.0);
        b[2] = static_cast<T>(1995840.0);
        b[3] = static_cast<T>(277200.0);
        b[4] = static_cast<T>(25200.0);
        b[5] = static_cast<T>(1512.0);
        b[6] = static_cast<T>(56.0);
        b[7] = T{1};
    } else {
        b[0] = static_cast<T>(17643225600.0);
        b[1] = static_cast<T>(8821612800.0);
        b[2] = static_cast<T>(2075673600.0);
        b[3] = static_cast<T>(302702400.0);
        b[4] = static_cast<T>(30270240.0);
        b[5] = static_cast<T>(2162160.0);
        b[6] = static_cast<T>(110880.0);
        b[7] = static_cast<T>(3960.0);
        b[8] = static_cast<T>(90.0);
        b[9] = T{1};
    }
    return b;
}

template<typename T, size_t N>
[[nodiscard]] constexpr damp::optional<Matrix<N, N, T>> expm_pade_low(const Matrix<N, N, T>& B, int m) {
    const damp::array<T, 10> coeff = expm_pade_coeff<T>(m);
    const Matrix<N, N, T>    ident = Matrix<N, N, T>::identity();
    const Matrix<N, N, T>    B2 = B * B;
    Matrix<N, N, T>          power = ident;
    Matrix<N, N, T>          odd = ident * coeff[1];
    Matrix<N, N, T>          even = ident * coeff[0];
    for (int k = 1; k <= m / 2; ++k) {
        power = power * B2;
        even = even + (power * coeff[static_cast<size_t>(2 * k)]);
        if ((2 * k) + 1 <= m) {
            odd = odd + (power * coeff[static_cast<size_t>((2 * k) + 1)]);
        }
    }
    const Matrix<N, N, T> U = B * odd;
    const Matrix<N, N, T> V = even;
    return solve(V - U, V + U);
}

template<typename T, size_t N>
[[nodiscard]] constexpr damp::optional<Matrix<N, N, T>> expm_pade13(const Matrix<N, N, T>& B) {
    constexpr T b0 = static_cast<T>(64764752532480000.0);
    constexpr T b1 = static_cast<T>(32382376266240000.0);
    constexpr T b2 = static_cast<T>(7771770303897600.0);
    constexpr T b3 = static_cast<T>(1187353796428800.0);
    constexpr T b4 = static_cast<T>(129060195264000.0);
    constexpr T b5 = static_cast<T>(10559470521600.0);
    constexpr T b6 = static_cast<T>(670442572800.0);
    constexpr T b7 = static_cast<T>(33522128640.0);
    constexpr T b8 = static_cast<T>(1323241920.0);
    constexpr T b9 = static_cast<T>(40840800.0);
    constexpr T b10 = static_cast<T>(960960.0);
    constexpr T b11 = static_cast<T>(16380.0);
    constexpr T b12 = static_cast<T>(182.0);
    constexpr T b13 = T{1};

    const Matrix<N, N, T> ident = Matrix<N, N, T>::identity();
    const Matrix<N, N, T> B2 = B * B;
    const Matrix<N, N, T> B4 = B2 * B2;
    const Matrix<N, N, T> B6 = B4 * B2;
    const Matrix<N, N, T> U2 = B6 * ((B6 * b13) + (B4 * b11) + (B2 * b9));
    const Matrix<N, N, T> U = B * (U2 + (B6 * b7) + (B4 * b5) + (B2 * b3) + (ident * b1));
    const Matrix<N, N, T> V2 = B6 * ((B6 * b12) + (B4 * b10) + (B2 * b8));
    const Matrix<N, N, T> V = V2 + (B6 * b6) + (B4 * b4) + (B2 * b2) + (ident * b0);
    return solve(V - U, V + U);
}

template<typename T, size_t N>
[[nodiscard]] constexpr damp::optional<Matrix<N, N, T>> expm_pade(const Matrix<N, N, T>& B, int degree) {
    if (degree == 13) {
        return expm_pade13(B);
    }
    return expm_pade_low(B, degree);
}

template<typename T, size_t N>
[[nodiscard]] constexpr ExpmTriangle expm_triangle(const Matrix<N, N, T>& A) {
    bool upper = true;
    bool lower = true;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            if (!(damp::abs(A(i, j)) > scalar_type_t<T>{0})) {
                continue;
            }
            if (i > j) {
                upper = false;
            }
            if (i < j) {
                lower = false;
            }
        }
    }
    if (upper && lower) {
        return ExpmTriangle::diagonal;
    }
    if (upper) {
        return ExpmTriangle::upper;
    }
    if (lower) {
        return ExpmTriangle::lower;
    }
    return ExpmTriangle::none;
}

template<typename T, size_t N>
[[nodiscard]] constexpr Matrix<N, N, T> expm_diagonal(const Matrix<N, N, T>& A) {
    Matrix<N, N, T> out = Matrix<N, N, T>::zeros();
    for (size_t i = 0; i < N; ++i) {
        out(i, i) = expm_entry_exp(A(i, i));
    }
    return out;
}

// exp(a) * sinh(x) / x. Taylor when x is small so equal eigenvalues do not cancel.
template<typename T>
[[nodiscard]] constexpr T expm_sinch(T a, T x) {
    using real_t = scalar_type_t<T>;
    if (damp::abs(x) < static_cast<real_t>(0.0135)) {
        const T x2 = x * x;
        const T series = T{1} + ((x2 / static_cast<T>(6)) * (T{1} + ((x2 / static_cast<T>(20)) * (T{1} + (x2 / static_cast<T>(42))))));
        return expm_entry_exp(a) * series;
    }
    return (expm_entry_exp(a + x) - expm_entry_exp(a - x)) / (x * static_cast<T>(2));
}

template<typename T, size_t N>
constexpr void expm_write_band(Matrix<N, N, T>& X, const Matrix<N, N, T>& A, size_t shift, ExpmTriangle tri, bool with_off) {
    Matrix<N, N, T> scaled = expm_div_pow2(A, shift);
    for (size_t k = 0; k < N; ++k) {
        X(k, k) = expm_entry_exp(scaled(k, k));
    }
    if (!with_off) {
        return;
    }
    for (size_t k = 0; k + 1 < N; ++k) {
        const size_t row = (tri == ExpmTriangle::upper) ? k : k + 1;
        const size_t col = (tri == ExpmTriangle::upper) ? k + 1 : k;
        const T      lam1 = scaled(k, k);
        const T      lam2 = scaled(k + 1, k + 1);
        const T      half = static_cast<T>(0.5);
        X(row, col) = scaled(row, col) * expm_sinch((lam1 + lam2) * half, (lam1 - lam2) * half);
    }
}

template<typename T, size_t N>
[[nodiscard]] constexpr Matrix<N, N, T> expm_square(Matrix<N, N, T> R, size_t s) {
    for (size_t i = 0; i < s; ++i) {
        R = R * R;
    }
    return R;
}

template<typename T, size_t N>
[[nodiscard]] constexpr Matrix<N, N, T> expm_square_triangular(Matrix<N, N, T> X, const Matrix<N, N, T>& A, size_t s, ExpmTriangle tri) {
    expm_write_band(X, A, s, tri, false);
    for (size_t i = s; i-- > 0;) {
        X = X * X;
        expm_write_band(X, A, i, tri, true);
    }
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            const bool clear = (tri == ExpmTriangle::upper) ? (i > j) : (i < j);
            if (clear) {
                X(i, j) = T{0};
            }
        }
    }
    return X;
}

} // namespace detail

/**
 * @brief Matrix exponential by scaling and squaring (Al-Mohy & Higham 2009)
 *
 * Chooses a Padé degree m ∈ {3, 5, 7, 9, 13} and a power-of-two scale s from
 * ‖Aᵏ‖₁^{1/k}, then forms exp(A) = (r_m(A/2ˢ))^{2ˢ}. A triangular argument
 * rewrites its diagonal from the scalar exponential while squaring, so a huge
 * off-diagonal does not erase a small diagonal.
 *
 * Series definition:
 * @f[
 *   \exp(A) = I + A + \frac{A^2}{2!} + \frac{A^3}{3!} + \cdots
 * @f]
 *
 * For linear ODEs @f$ \dot x = A x @f$, the solution is
 * @f$ x(t) = \exp(A t)\, x(0) @f$. Used by ZOH discretization and exact LTI steps.
 *
 * The Padé step solves @f$ (V-U) R = V+U @f$ with mat::solve (no explicit
 * inverse). A non-finite ‖A‖₁ (a NaN or ±Inf entry, or a column sum past the
 * floating-point range) or a failed Padé solve returns damp::nullopt.
 *
 * @note Compare with MATLAB®'s expm(A).
 * @see Al-Mohy and Higham, "A New Scaling and Squaring Algorithm for the
 *      Matrix Exponential" (SIAM J. Matrix Anal. Appl., 2009)
 * @see Higham, "Functions of Matrices" (2008), §10.3
 *
 * @tparam T Element type
 * @tparam N Matrix dimension
 * @param A Square matrix
 * @return exp(A), or damp::nullopt if the input is non-finite or the Padé solve fails
 */
template<typename T, size_t N>
[[nodiscard]] constexpr damp::optional<Matrix<N, N, T>> expm(const Matrix<N, N, T>& A) {
    // one_norm's max skips a NaN column sum, so the entries are checked themselves.
    if (!mat::all_of(A, [](const T& v) { return damp::isfinite(damp::abs(v)); }) || !damp::isfinite(one_norm(A))) {
        return damp::nullopt;
    }

    const detail::ExpmTriangle tri = detail::expm_triangle(A);
    if (tri == detail::ExpmTriangle::diagonal) {
        return detail::expm_diagonal(A);
    }

    const detail::ExpmScale               scale = detail::expm_scale(A);
    const Matrix<N, N, T>                 scaled = detail::expm_div_pow2(A, scale.s);
    const damp::optional<Matrix<N, N, T>> R = detail::expm_pade(scaled, scale.degree);
    if (!R) {
        return damp::nullopt;
    }
    if (tri == detail::ExpmTriangle::none) {
        return detail::expm_square(R.value(), scale.s);
    }
    return detail::expm_square_triangular(R.value(), A, scale.s, tri);
}

namespace detail {

/**
 * @brief Solve R_ii·X + X·R_jj = C for the p×q block X (p, q ∈ {1, 2}) of a quasi-triangular square root.
 *
 * The Kronecker form (I_q ⊗ R_ii + R_jjᵀ ⊗ I_p)·vec(X) = vec(C), at most 4×4, by pivot_eliminate().
 * @return false on an exactly singular system (two eigenvalues of R sum to zero).
 */
template<typename T, size_t N>
[[nodiscard]] constexpr bool sqrtm_block_solve(Matrix<N, N, T>& R, const Matrix<N, N, T>& C, size_t i0, size_t p, size_t j0, size_t q) {
    Matrix<4, 4, T> K{};
    ColVec<4, T>    rhs{};
    const size_t    n = p * q;
    for (size_t c = 0; c < q; ++c) {
        for (size_t r = 0; r < p; ++r) {
            const size_t row = r + (c * p);
            rhs(row) = C(i0 + r, j0 + c);
            for (size_t s = 0; s < p; ++s) {
                K(row, s + (c * p)) += R(i0 + r, i0 + s);
            }
            for (size_t s = 0; s < q; ++s) {
                K(row, r + (s * p)) += R(j0 + s, j0 + c);
            }
        }
    }
    if (pivot_eliminate(K, rhs, n) == T{0}) {
        return false;
    }
    for (size_t k = n; k-- > 0;) {
        T sum = rhs(k);
        for (size_t m = k + 1; m < n; ++m) {
            sum -= K(k, m) * rhs(m);
        }
        rhs(k) = sum / K(k, k);
    }
    for (size_t c = 0; c < q; ++c) {
        for (size_t r = 0; r < p; ++r) {
            R(i0 + r, j0 + c) = rhs(r + (c * p));
        }
    }
    return true;
}

/**
 * @brief Principal square root of the diagonal block of @p Tq at @p i0 (size 1, or 2 for a complex pair) into @p R.
 *
 * A standardized 2×2 block [[a, b], [c, a]] with bc < 0 has eigenvalues a ± iμ, μ = √(−bc); its root is
 * αI + (T_ii − aI)/(2α) with α = Re √(a + iμ) > 0, taken from β = Im √(a + iμ) when a < 0 so that
 * neither is formed by cancellation.
 * @return false for a real eigenvalue that is not positive (no real principal root, or A singular).
 */
template<typename T, size_t N>
[[nodiscard]] constexpr bool sqrtm_diagonal_block(const Matrix<N, N, T>& Tq, Matrix<N, N, T>& R, size_t i0, size_t size) {
    if (size == 1) {
        if (!(Tq(i0, i0) > T{0})) {
            return false;
        }
        R(i0, i0) = damp::sqrt(Tq(i0, i0));
        return true;
    }
    const T a = Tq(i0, i0);
    const T b = Tq(i0, i0 + 1);
    const T c = Tq(i0 + 1, i0);
    const T mu = damp::sqrt(damp::abs(b)) * damp::sqrt(damp::abs(c));
    const T m = damp::max(damp::abs(a), mu);
    const T mod = m * damp::sqrt(((a / m) * (a / m)) + ((mu / m) * (mu / m))); // |a + iμ|, scaled
    T       alpha{};
    if (a >= T{0}) {
        alpha = damp::sqrt(static_cast<T>(0.5) * (mod + a));
    } else {
        const T beta = damp::sqrt(static_cast<T>(0.5) * (mod - a));
        alpha = mu / (T{2} * beta);
    }
    const T inv2a = T{1} / (T{2} * alpha);
    R(i0, i0) = alpha;
    R(i0 + 1, i0 + 1) = alpha;
    R(i0, i0 + 1) = b * inv2a;
    R(i0 + 1, i0) = c * inv2a;
    return true;
}

/**
 * @brief R = √Tq for a quasi-triangular @p Tq with standardized 2×2 blocks; R has the same block structure.
 *
 * Each 1×1 or 2×2 diagonal block takes its closed-form root (sqrtm_diagonal_block()), and each
 * off-diagonal block solves a Sylvester equation of at most 4 unknowns (sqrtm_block_solve()).
 * Returns damp::nullopt for a real eigenvalue that is not positive or an exactly singular
 * Sylvester system.
 *
 * @see "Computing real square roots of a real matrix" (Higham, 1987), Linear Algebra Appl. 88/89
 */
template<typename T, size_t N>
[[nodiscard]] constexpr damp::optional<Matrix<N, N, T>> sqrtm_quasi_triangular(const Matrix<N, N, T>& Tq) {
    // Diagonal blocks: a nonzero subdiagonal marks a standardized 2×2 block (a complex pair).
    damp::array<size_t, N> start{};
    damp::array<size_t, N> size{};
    size_t                 nb = 0;
    for (size_t i = 0; i < N; i += size[nb - 1]) {
        start[nb] = i;
        size[nb] = (i + 1 < N && Tq(i + 1, i) != T{0}) ? 2 : 1;
        ++nb;
    }

    Matrix<N, N, T> R{};
    Matrix<N, N, T> C{};
    for (size_t jb = 0; jb < nb; ++jb) {
        if (!sqrtm_diagonal_block(Tq, R, start[jb], size[jb])) {
            return damp::nullopt;
        }
        for (size_t ib = jb; ib-- > 0;) {
            const size_t i0 = start[ib];
            const size_t j0 = start[jb];
            for (size_t r = 0; r < size[ib]; ++r) {
                for (size_t c = 0; c < size[jb]; ++c) {
                    T sum = Tq(i0 + r, j0 + c);
                    for (size_t k = i0 + size[ib]; k < j0; ++k) {
                        sum -= R(i0 + r, k) * R(k, j0 + c);
                    }
                    C(i0 + r, j0 + c) = sum;
                }
            }
            if (!sqrtm_block_solve(R, C, i0, size[ib], j0, size[jb])) {
                return damp::nullopt;
            }
        }
    }
    return R;
}

/**
 * @brief Principal square root of a real A by the real Schur method; the float and double path of sqrt().
 *
 * Returns damp::nullopt for a non-finite entry, a real eigenvalue that is not positive, a Schur form
 * that did not converge, or a result whose residual fails the bound documented on sqrt().
 */
template<typename T, size_t N>
[[nodiscard]] constexpr damp::optional<Matrix<N, N, T>> sqrtm_real_schur(const Matrix<N, N, T>& A) {
    if (!mat::all_of(A, [](T v) { return damp::isfinite(v); })) {
        return damp::nullopt;
    }
    if constexpr (N == 1) {
        if (!(A(0, 0) > T{0})) {
            return damp::nullopt;
        }
        return Matrix<1, 1, T>{{damp::sqrt(A(0, 0))}};
    } else {
        // Balance first (a permutation and a diagonal similarity by powers of two, exact):
        // √A = D·√(D⁻¹AD)·D⁻¹, so the accuracy does not depend on the units of each state.
        const auto rs = real_schur(A);
        if (!rs.converged) {
            return damp::nullopt;
        }
        const auto&                           e = rs.e;
        const auto&                           Z = rs.Z;
        const Matrix<N, N, T>                 Ab = balanced(A, e);
        const damp::optional<Matrix<N, N, T>> R = sqrtm_quasi_triangular(rs.T_form);
        if (!R) {
            return damp::nullopt;
        }

        const Matrix<N, N, T> X = Z * R.value() * Z.t(); // √(D⁻¹AD)
        const T               xn = X.norm();
        const Matrix<N, N, T> residual = (X * X) - Ab;
        if (!damp::isfinite(xn) || !(residual.norm() <= static_cast<T>(N) * T{64} * eps<T>() * xn * xn)) {
            return damp::nullopt;
        }
        damp::array<int, N> neg{};
        for (size_t i = 0; i < N; ++i) {
            neg[i] = -e[i];
        }
        return balanced(X, neg);
    }
}

/**
 * @brief Principal square root by the Denman–Beavers iteration; the damp::complex path of sqrt().
 *
 * There is no complex Schur form here, so a complex A iterates Y₀ = A, Z₀ = I,
 * Yₖ₊₁ = (Yₖ + Zₖ⁻¹)/2, Zₖ₊₁ = (Zₖ + Yₖ⁻¹)/2 until Y → √A, at most 50 times. Returns
 * damp::nullopt on a singular iterate or when the iteration does not converge.
 *
 * @see Higham, "Functions of Matrices" (2008), §6.3
 */
template<typename T, size_t N>
[[nodiscard]] constexpr damp::optional<Matrix<N, N, T>> sqrtm_denman_beavers(const Matrix<N, N, T>& A) {
    using real_t = scalar_type_t<T>;
    Matrix<N, N, T> Y = A;
    Matrix<N, N, T> Z = Matrix<N, N, T>::identity();
    real_t          prev_diff = damp::numeric_limits<real_t>::max();

    for (int iter = 0; iter < 50; ++iter) {
        // Inverse is the deliverable of each Denman–Beavers half-step.
        auto Y_inv = Y.inverse();
        auto Z_inv = Z.inverse();

        if (!Y_inv || !Z_inv) {
            return damp::nullopt;
        }

        Matrix<N, N, T> Y_next = (Y + Z_inv.value()) * static_cast<T>(0.5);
        Matrix<N, N, T> Z_next = (Z + Y_inv.value()) * static_cast<T>(0.5);

        real_t diff = real_t{0};
        real_t size = real_t{0};
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                diff += damp::abs(Y_next(i, j) - Y(i, j));
                size += damp::abs(Y_next(i, j));
            }
        }

        Y = Y_next;
        Z = Z_next;

        //! Relative stop at roundoff level; quadratic convergence then stalls,
        //! so a step that no longer shrinks below √eps·‖Y‖ is converged too.
        if (diff <= static_cast<real_t>(N) * real_t{64} * eps<T>() * size
            || (diff <= damp::sqrt(eps<T>()) * size && diff >= prev_diff)) {
            return Y;
        }
        prev_diff = diff;
    }

    return damp::nullopt;
}

} // namespace detail

/**
 * @brief Principal matrix square root by the real Schur method
 *
 * Computes the principal square root @f$ X = \sqrt{A} @f$ with @f$ X\cdot X = A @f$ and every
 * eigenvalue of X in the open right half-plane. A real Schur form @f$ A = Z T Z^\top @f$
 * of the balanced matrix D⁻¹AD (balancing, Hessenberg reduction and Francis double-shift QR,
 * detail::real_schur as in compute_eigenvalues, with the permutation folded into Z; √A = D·√(D⁻¹AD)·D⁻¹, so graded units cost no accuracy) reduces the
 * problem to the quasi-triangular T: each 1×1 or 2×2 diagonal block has a closed-form root, and
 * each off-diagonal block of R solves a Sylvester equation of at most 4 unknowns,
 * @f[
 *   R_{ii} R_{ij} + R_{ij} R_{jj} = T_{ij} - \sum_{k=i+1}^{j-1} R_{ik} R_{kj},
 * @f]
 * after which @f$ X = Z R Z^\top @f$. The method is backward stable: the residual
 * @f$ \|X^2 - A\| @f$ is a small multiple of @f$ \varepsilon \|X\|^2 @f$, and the forward error
 * follows the condition number of the square root, not κ(A). Denman–Beavers, used here before
 * for real matrices, inverts iterates with condition κ(A), which failed in float once κ(A)
 * passed about 10⁴.
 *
 * Returns damp::nullopt for a non-finite entry, a real eigenvalue that is not positive (no real
 * principal root; A singular), a Schur form that did not converge, or a result whose residual
 * ‖X² − A‖_F, taken on the balanced matrix, exceeds N·64·ε·‖X‖²_F.
 *
 * A damp::complex matrix has no Schur form here and takes the Denman–Beavers iteration
 * (detail::sqrtm_denman_beavers), which returns damp::nullopt on a singular iterate or when it
 * does not converge within 50 steps; its iterates carry the condition number κ(A).
 *
 * @note Compare with MATLAB®'s sqrtm(A).
 * @see "Computing real square roots of a real matrix" (Higham, 1987), Linear Algebra Appl. 88/89
 * @see Higham, "Functions of Matrices" (2008), §6.2 (Schur method) and Alg. 6.7, §6.3 (Denman–Beavers)
 *
 * @tparam T float or double (real Schur method), or damp::complex thereof (Denman–Beavers iteration)
 * @tparam N Matrix dimension
 * @param A Square matrix (no eigenvalues on the closed negative real axis)
 * @return Principal √A, or damp::nullopt on failure
 */
template<typename T, size_t N>
[[nodiscard]] constexpr damp::optional<Matrix<N, N, T>> sqrt(const Matrix<N, N, T>& A) {
    if constexpr (damp::is_floating_point_v<T>) {
        return detail::sqrtm_real_schur(A);
    } else {
        return detail::sqrtm_denman_beavers(A);
    }
}

/**
 * @brief Principal matrix logarithm via inverse scaling and squaring
 *
 * Computes @f$ X = \log(A) @f$ such that @f$ \exp(X) = A @f$ (principal branch).
 * Balances A first (an exact diagonal similarity by powers of two,
 * @f$ \log A = D \log(D^{-1} A D) D^{-1} @f$, as in sqrt()), so graded units cost no accuracy.
 * A real A that is not already within 1/2 of I (∞-norm) then takes the real Schur form
 * @f$ D^{-1} A D = Z T Z^\top @f$ once (detail::real_schur, as in sqrt(); its D comes from a
 * balance that first isolates the eigenvalues a permutation exposes), and the roots (sqrt()'s quasi-triangular recurrence)
 * and the series run on T: @f$ \log A = D Z \log(T) Z^\top D^{-1} @f$. A Schur form per root
 * repeats its roundoff at every step, amplified by @f$ 2^s @f$, and a rotation of a
 * near-identity matrix mixes ε·I into its small entries. A damp::complex A takes sqrt() of the
 * balanced matrix at each step.
 * The roots continue until @f$ \|A_{\mathrm{scaled}} - I\|_\infty \le 1/2 @f$, at most 20 (at
 * the cap the series still converges, since @f$ \rho(Z) \lesssim 10^{-3} @f$ after 20 roots),
 * and the result evaluates the Gregory (atanh) series
 * @f[
 *   \log(A) = 2 \sum_{k\ \mathrm{odd}} \frac{Z^k}{k},\quad Z = (A + I)^{-1}(A - I)
 * @f]
 * and multiplies by @f$ 2^s @f$. With @f$ \|Z\| \lesssim 1/3 @f$ the odd-power
 * series reaches double precision in about 15 terms; the plain
 * @f$ \log(I + X) @f$ Taylor series at @f$ \|X\| = 1/2 @f$ stalls near 10⁻⁸.
 *
 * @p A must be invertible with no eigenvalues on the closed negative real axis.
 * Returns damp::nullopt for a non-finite entry, a Schur form that did not converge, or a
 * square-root step that fails (an eigenvalue on the closed negative real axis, or a complex
 * square root that does not converge).
 *
 * @note Compare with MATLAB®'s logm(A).
 * @see Higham, "Functions of Matrices" (2008), §11, §11.5 (inverse scaling and squaring, Schur)
 * @see sqrt() for the square-root step
 *
 * @tparam T Element type (real or damp::complex)
 * @tparam N Matrix dimension
 * @param A Square matrix
 * @return Principal log(A), or damp::nullopt on failure
 */
template<typename T, size_t N>
[[nodiscard]] constexpr damp::optional<Matrix<N, N, T>> log(const Matrix<N, N, T>& A) {
    // Balance first (exact, powers of two): log A = D·log(D⁻¹AD)·D⁻¹, so the number of square
    // roots, and with it the accuracy, does not depend on the units of each state. The exponents
    // come from the magnitudes, so a damp::complex A is balanced too. A non-finite entry is
    // refused here: its NaN norm would pass the stop test below and skip the roots.
    if (!mat::all_of(A, [](const T& v) { return detail::is_finite_entry(v); })) {
        return damp::nullopt;
    }
    Matrix<N, N, scalar_type_t<T>> mag{};
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            mag(i, j) = damp::abs(A(i, j));
        }
    }
    auto   e = detail::balance_exponents(mag);
    Matrix I = Matrix<N, N, T>::identity();
    Matrix A_scaled = detail::balanced(A, e);
    size_t s = 0;

    // Real A: when roots are needed, take the real Schur form of the balanced matrix once and
    // run the roots and the series on its quasi-triangular factor, log(A) = D·Z·log(T)·Zᵀ·D⁻¹.
    // A near-identity A skips the rotation, which would mix eps·I into its small entries.
    constexpr bool                   schur = damp::is_floating_point_v<T> && N > 1;
    [[maybe_unused]] Matrix<N, N, T> Z_schur = I;
    if constexpr (schur) {
        if (infinity_norm(A_scaled - I) > T{0.5}) {
            const auto rs = detail::real_schur(A); // its own balance: a permutation, then the scaling
            if (!rs.converged) {
                return damp::nullopt;
            }
            A_scaled = rs.T_form;
            Z_schur = rs.Z;
            e = rs.e;
        }
    }

    // A_scaled = A^(1/2^s) until ||A_scaled - I||_∞ ≤ 1/2; for a real A the loop runs only if the Schur step did
    while (infinity_norm(A_scaled - I) > scalar_type_t<T>{0.5} && s < 20) {
        damp::optional<Matrix<N, N, T>> S;
        if constexpr (schur) {
            S = detail::sqrtm_quasi_triangular(A_scaled);
        } else {
            S = sqrt(A_scaled);
        }
        if (!S) {
            return damp::nullopt;
        }
        A_scaled = S.value();
        ++s;
    }

    // Gregory series, Z = (A + I)⁻¹(A − I) (A and its rational functions commute)
    const auto Zopt = mat::solve(Matrix(A_scaled + I), Matrix(A_scaled - I));
    if (!Zopt) {
        return damp::nullopt;
    }
    const Matrix Z = Zopt.value();
    const Matrix Z2 = Z * Z;
    Matrix       Z_power = Z;
    Matrix       result = Z;
    for (size_t k = 3; k <= 41; k += 2) {
        Z_power = Z_power * Z2;
        const scalar_type_t<T> inv_k = scalar_type_t<T>{1} / static_cast<scalar_type_t<T>>(k); // real, for complex T too
        result = result + Z_power * T{inv_k};
    }
    result = result * T{2};

    // log(A) = 2^s * log(A^(1/2^s))
    const T scale{static_cast<scalar_type_t<T>>(size_t{1} << s)};
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            result(i, j) *= scale;
        }
    }
    if constexpr (schur) {
        result = Z_schur * result * Z_schur.t();
    }

    damp::array<int, N> neg{};
    for (size_t i = 0; i < N; ++i) {
        neg[i] = -e[i];
    }
    return detail::balanced(result, neg);
}

/**
 * @brief Integer matrix power via binary exponentiation
 *
 * Computes @f$ A^p @f$. For @f$ p < 0 @f$, returns @f$ (A^{-1})^{|p|} @f$: inverting A first
 * costs a condition number κ(A), where inverting @f$ A^{|p|} @f$ would cost up to
 * @f$ \kappa(A)^{|p|} @f$. Returns damp::nullopt if A is singular — no identity fallback.
 *
 * @note Compare with MATLAB®'s mpower(A, p) for integer @p p.
 *
 * @tparam T Element type
 * @tparam N Matrix dimension
 * @tparam I Integer type of the exponent (a floating exponent takes the real overload)
 * @param A Square matrix
 * @param p Integer exponent
 * @return Aᵖ, or damp::nullopt if p < 0 and A is singular
 */
template<typename T, size_t N, typename I>
    requires damp::is_integral_v<I>
[[nodiscard]] constexpr damp::optional<Matrix<N, N, T>> pow(const Matrix<N, N, T>& A, I p) {
    if (p == I{0}) {
        return Matrix<N, N, T>::identity();
    }

    bool negate = false;
    if constexpr (damp::is_signed_v<I>) {
        negate = p < I{0};
    }
    // Magnitude in unsigned arithmetic: −min() overflows a signed type.
    uintmax_t e = negate ? (uintmax_t{0} - static_cast<uintmax_t>(p)) : static_cast<uintmax_t>(p);

    Matrix<N, N, T> base = A;
    if (negate) {
        const auto inv = A.inverse();
        if (!inv) {
            return damp::nullopt;
        }
        base = *inv;
    }

    Matrix<N, N, T> result = Matrix<N, N, T>::identity();
    while (e > 0u) {
        if ((e & 1u) != 0u) {
            result = result * base;
        }
        base = base * base;
        e >>= 1u;
    }
    return result;
}

/**
 * @brief Matrix power Aᵖ with a real exponent, via exp(p log(A))
 *
 * An exponent within default_tol of the nearest integer (from either side) reduces to the
 * integer @ref pow path. Other exponents require a successful principal logarithm:
 * @f[
 *   A^p = \exp\bigl(p\,\log(A)\bigr)
 * @f]
 *
 * @note Compare with MATLAB®'s mpower(A, p).
 * @see log(), expm()
 *
 * @tparam T Element type
 * @tparam N Matrix dimension
 * @param A Square matrix (real or damp::complex)
 * @param p Real exponent, of A's real scalar type (not deduced: a double literal works on a float matrix)
 * @return Aᵖ, or damp::nullopt if @ref log or @ref expm fails (or integer inverse fails)
 */
template<typename T, size_t N>
[[nodiscard]] constexpr damp::optional<Matrix<N, N, T>> pow(const Matrix<N, N, T>& A, damp::type_identity_t<scalar_type_t<T>> p) {
    using R = scalar_type_t<T>;
    // Integer case (constexpr-safe, no std::modf). Casting an out-of-range
    // value to int is UB, so only exponents inside int's range take this path.
    constexpr R int_limit = static_cast<R>(damp::numeric_limits<int>::max() / 2);
    if (damp::abs(p) < int_limit) {
        const R nearest = damp::floor(p + static_cast<R>(0.5));
        if (damp::abs(p - nearest) < default_tol<R>()) {
            return pow(A, static_cast<int>(nearest));
        }
    }
    const auto L = log(A);
    if (!L) {
        return damp::nullopt;
    }
    return expm(Matrix<N, N, T>(L.value() * T{p}));
}

namespace detail {

/**
 * @brief The s with ‖A‖∞/2ˢ ≤ 1/2, for the double-angle recovery of sin/cos; sinh/cosh take the series only at s = 0.
 *
 * Bounded by the exponent range: a finite norm reaches 1/2 within max_exponent + 1 halvings. A
 * non-finite norm (a NaN or ±Inf entry, or a row sum past the range) returns 0, and the series
 * then carries the non-finite values into the result.
 */
template<typename T, size_t N>
[[nodiscard]] constexpr size_t trig_scale_exponent(const Matrix<N, N, T>& A) {
    using real_t = scalar_type_t<T>;
    real_t nrm = infinity_norm(A);
    if (!damp::isfinite(nrm)) {
        return 0;
    }
    constexpr size_t cap = static_cast<size_t>(damp::numeric_limits<real_t>::max_exponent) + 1;
    size_t           s = 0;
    while (nrm > real_t{0.5} && s < cap) {
        nrm *= real_t{0.5};
        ++s;
    }
    return s;
}

/** @brief The odd and even parts of a sin/cos or sinh/cosh pair. */
template<typename T, size_t N>
struct TrigPair {
    Matrix<N, N, T> odd{};  ///< sin or sinh
    Matrix<N, N, T> even{}; ///< cos or cosh
};

/**
 * @brief sin/cos (σ = −1) or sinh/cosh (σ = +1) of A by scaling, Taylor series and double-angle recovery.
 *
 * The series run to X²¹/21! and X²⁰/20! on X = A/2ˢ with ‖X‖∞ ≤ 1/2, then s steps of
 * odd(2X) = 2·odd(X)·even(X), even(2X) = 2·even(X)² − I, identities shared by both pairs. A small
 * A takes s = 0 and its series, so sinh(A) ≈ A keeps its relative accuracy; the difference
 * (exp(A) − exp(−A))/2 cancels to eps/‖A‖. The odd part is skipped when @p with_odd is false.
 */
template<typename T, size_t N>
[[nodiscard]] constexpr TrigPair<T, N> trig_pair(const Matrix<N, N, T>& A, T sigma, bool with_odd) {
    const Matrix<N, N, T> I = Matrix<N, N, T>::identity();
    const size_t          s = trig_scale_exponent(A);
    const Matrix<N, N, T> X = expm_div_pow2(A, s); // halve A itself: a scalar 2⁻ˢ is subnormal past s = 126 (float) and FTZ zeroes it
    const Matrix<N, N, T> X2 = X * X;

    TrigPair<T, N>  out{X, I};
    Matrix<N, N, T> odd_power = X;
    Matrix<N, N, T> even_power = I;
    T               odd_factorial = T{1};
    T               even_factorial = T{1};
    T               sign = sigma;
    for (size_t n = 2; n <= 20; n += 2) {
        even_factorial *= static_cast<T>(n - 1) * static_cast<T>(n);
        even_power = even_power * X2;
        out.even = out.even + (even_power * (sign / even_factorial));
        if (with_odd) {
            odd_factorial *= static_cast<T>(n) * static_cast<T>(n + 1);
            odd_power = odd_power * X2;
            out.odd = out.odd + (odd_power * (sign / odd_factorial));
        }
        sign *= sigma;
    }

    for (size_t i = 0; i < s; ++i) {
        if (with_odd) {
            out.odd = out.odd * out.even * T{2};
        }
        out.even = (out.even * out.even * T{2}) - I;
    }
    return out;
}

/**
 * @brief sinh/cosh of A: the series of trig_pair() when ‖A‖∞ ≤ 1/2, (exp(A) ∓ exp(−A))/2 otherwise.
 *
 * The exponential difference cancels only for a small A; above ‖A‖∞ = 1/2 it does not, and expm's
 * scaling (θ₁₃ ≈ 5.4) squares fewer times than the double-angle recovery, which loses a factor 2
 * per step. A failed expm (a NaN entry the norm skipped, a column sum past the range, or a failed
 * Padé solve) falls back to the series, which carries the non-finite values into the result.
 */
template<typename T, size_t N>
[[nodiscard]] constexpr TrigPair<T, N> hyperbolic_pair(const Matrix<N, N, T>& A, bool with_odd) {
    if (trig_scale_exponent(A) > 0) {
        const damp::optional<Matrix<N, N, T>> ep = expm(A);
        const damp::optional<Matrix<N, N, T>> em = expm(Matrix<N, N, T>(A * T{-1}));
        if (ep && em) {
            return {(ep.value() - em.value()) * static_cast<T>(0.5), (ep.value() + em.value()) * static_cast<T>(0.5)};
        }
    }
    return trig_pair(A, T{1}, with_odd);
}

} // namespace detail

/**
 * @brief Compute sin(A) and cos(A) together via scaling and double-angle reconstruction
 *
 * More efficient than calling sin() and cos() separately — computes both
 * with one scaling pass and shared Taylor series evaluation.
 *
 * Scales A down so ||A/2ˢ|| ≤ 0.5 (where the Taylor series converges
 * accurately), then recovers via repeated double-angle formulas:
 *
 *     sin(2A) = 2·sin(A)·cos(A)
 *     cos(2A) = 2·cos²(A) − I
 *
 * Since sin(A) and cos(A) are polynomials in A, they commute with each other,
 * so the scalar double-angle formulas apply directly. A non-finite ‖A‖∞ gives a
 * non-finite result.
 *
 * @see Higham, "Functions of Matrices" (2008), §12.3
 *
 * @param A Square matrix
 * @return {sin(A), cos(A)}
 */
template<typename T, size_t N>
[[nodiscard]] constexpr damp::pair<Matrix<N, N, T>, Matrix<N, N, T>> sincos(const Matrix<N, N, T>& A) {
    const auto p = detail::trig_pair(A, T{-1}, true);
    return {p.odd, p.even};
}

/**
 * @brief Matrix sine via scaling and double-angle reconstruction
 *
 * @note Compare with MATLAB®'s funm(A, @@sin).
 * @see sincos() to compute both sin(A) and cos(A) in one call
 * @see Higham, "Functions of Matrices" (2008), §12.3
 *
 * @param A Square matrix
 * @return Matrix sine sin(A)
 */
template<typename T, size_t N>
[[nodiscard]] constexpr Matrix<N, N, T> sin(const Matrix<N, N, T>& A) {
    return sincos(A).first;
}

/**
 * @brief Matrix cosine via scaling and double-angle reconstruction
 *
 * cos(2X) = 2·cos²(X) − I needs no sine, so the sine series and its doubling are skipped.
 *
 * @note Compare with MATLAB®'s funm(A, @@cos).
 * @see sincos() to compute both sin(A) and cos(A) in one call
 * @see Higham, "Functions of Matrices" (2008), §12.3
 *
 * @param A Square matrix
 * @return Matrix cosine cos(A)
 */
template<typename T, size_t N>
[[nodiscard]] constexpr Matrix<N, N, T> cos(const Matrix<N, N, T>& A) {
    return detail::trig_pair(A, T{-1}, false).even;
}

/**
 * @brief Matrix hyperbolic sine sinh(A) = (exp(A) − exp(−A))/2
 *
 * For ‖A‖∞ ≤ 1/2, the Taylor series, so a small A keeps its relative accuracy where the
 * exponential difference would cancel; for a larger A, (expm(A) − expm(−A))/2, which does
 * not cancel there and keeps expm's accuracy.
 *
 * @note Compare with MATLAB®'s funm(A, @@sinh).
 * @see expm(), cosh()
 *
 * @param A Square matrix
 * @return sinh(A)
 */
template<typename T, size_t N>
[[nodiscard]] constexpr Matrix<N, N, T> sinh(const Matrix<N, N, T>& A) {
    return detail::hyperbolic_pair(A, true).odd;
}

/**
 * @brief Matrix hyperbolic cosine cosh(A) = (exp(A) + exp(−A))/2
 *
 * For ‖A‖∞ ≤ 1/2, the even Taylor series; for a larger A, (expm(A) + expm(−A))/2.
 *
 * @note Compare with MATLAB®'s funm(A, @@cosh).
 * @see expm(), sinh()
 *
 * @param A Square matrix
 * @return cosh(A)
 */
template<typename T, size_t N>
[[nodiscard]] constexpr Matrix<N, N, T> cosh(const Matrix<N, N, T>& A) {
    return detail::hyperbolic_pair(A, false).even;
}
} // namespace mat
} // namespace damp