// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file functions.hpp
 * @brief Matrix functions (expm, sqrtm, funm, ...)
 */

#include <limits>

#include "core.hpp"
#include "damp/backend.hpp"
#include "solve.hpp" // expm()/sqrtm() use mat::solve; include directly so this
// does not depend on aggregation order (formatter sorts includes)

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
 * @brief Frobenius norm ‖A‖F = √(Σᵢⱼ |aᵢⱼ|²)
 *
 * Always returns a real scalar (@c scalar_type_t&lt;T&gt;), including for complex matrices.
 *
 * @f[
 *   \|A\|_F = \sqrt{\sum_{ij} |a_{ij}|^2}
 * @f]
 *
 * @note Compare with MATLAB®'s norm(A, 'fro').
 *
 * @tparam T Element type (real or damp::complex)
 * @tparam N Matrix dimension
 * @param A Square matrix
 * @return Frobenius norm of the matrix
 */
template<typename T, size_t N>
[[nodiscard]] constexpr scalar_type_t<T> frobenius_norm(const Matrix<N, N, T>& A) {
    using real_t = scalar_type_t<T>;
    real_t sum_squares = real_t{0};
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            const real_t abs_val = damp::abs(A(i, j));
            sum_squares += abs_val * abs_val;
        }
    }
    return damp::sqrt(sum_squares);
}

/**
 * @brief Spectral norm ‖A‖₂ = σₘₐₓ(A)
 *
 * @f[
 *   \|A\|_2 = \sigma_{\max}(A) = \sqrt{\lambda_{\max}(A^{\mathrm{H}} A)}
 * @f]
 * where @f$ A^{\mathrm{H}} @f$ is the conjugate transpose.
 *
 * Computed via power iteration on @f$ A^{\mathrm{H}} A @f$ (geometric rate
 * @f$ \sigma_1/\sigma_2 @f$). Intended for the small fixed-size matrices typical
 * in control design (@f$ n \lesssim 20 @f$). Always returns a real scalar.
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
    using real_t = scalar_type_t<T>;
    constexpr real_t tol = default_tol<T>();
    constexpr size_t max_iter = 100;

    // AᴴA (conjugate transpose) — required for complex matrices; equals AᵀA when real.
    Matrix<N, N, T> AtA = A.conjugate_transpose() * A;

    // Initial vector: all ones, then normalize
    ColVec<N, T> v;
    for (size_t i = 0; i < N; ++i) {
        v[i] = T{1};
    }

    {
        const real_t n = v.norm();
        if (n > tol) {
            v = v * (T{1} / static_cast<T>(n));
        }
    }

    real_t lambda = real_t{0};
    for (size_t iter = 0; iter < max_iter; ++iter) {
        ColVec<N, T> w = AtA * v;
        const real_t new_lambda = w.norm();
        if (new_lambda < tol) {
            return real_t{0};
        }

        v = w * (T{1} / static_cast<T>(new_lambda));

        if (damp::abs(new_lambda - lambda) < tol * damp::abs(new_lambda)) {
            return damp::sqrt(new_lambda);
        }
        lambda = new_lambda;
    }
    return damp::sqrt(lambda);
}

/**
 * @brief Matrix determinant det(A)
 *
 * Cofactor expansion for @f$ N \le 4 @f$ (exact, fast); LU for larger @f$ N @f$
 * (@f$ \det(A) = \mathrm{sign}(P)\prod_i U_{ii} @f$). Returns 0 if LU reports
 * singularity.
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
        // General case: det(A) = det(P) * det(L) * det(U) = sign * product(U_ii)
        auto lu = lu_decomposition(A);
        if (!lu) {
            return T{0}; // singular
        }

        const auto& [L, U, piv] = lu.value();

        // Compute sign from permutation parity
        auto   p = piv;
        size_t swaps = 0;
        for (size_t i = 0; i < N; ++i) {
            while (p[i] != i) {
                damp::swap(p[i], p[p[i]]);
                ++swaps;
            }
        }

        T d = (swaps % 2 == 0) ? T{1} : T{-1};
        for (size_t i = 0; i < N; ++i) {
            d *= U(i, i);
        }
        return d;
    }
}

/**
 * @brief Matrix rank via Gaussian elimination with partial pivoting
 *
 * Works for rectangular and square matrices, real or complex scalars (the
 * pivot magnitude uses damp::abs, which is defined for both). Single source for
 * the structural rank tests in stability.hpp and riccati.hpp.
 *
 * @tparam R   Rows
 * @tparam C   Columns
 * @tparam T   Element type (float/double or damp::complex thereof)
 * @param M    Input matrix
 * @param tol  Magnitude below which a pivot is treated as zero
 * @return Rank of the matrix (number of linearly independent rows/columns)
 */
template<size_t R, size_t C, typename T>
[[nodiscard]] constexpr size_t rank(const Matrix<R, C, T>& M, scalar_type_t<T> tol = default_tol<T>()) {
    Matrix<R, C, T> work = M;
    size_t          r = 0;
    for (size_t col = 0; col < C && r < R; ++col) {
        // Find the largest-magnitude pivot in this column at or below row r.
        size_t pivot = r;
        auto   max_val = damp::abs(work(r, col));
        for (size_t i = r + 1; i < R; ++i) {
            const auto val = damp::abs(work(i, col));
            if (val > max_val) {
                max_val = val;
                pivot = i;
            }
        }
        if (max_val < tol) {
            continue; // column is dependent on the ones already pivoted
        }
        if (pivot != r) {
            for (size_t j = 0; j < C; ++j) {
                damp::swap(work(r, j), work(pivot, j));
            }
        }
        for (size_t i = r + 1; i < R; ++i) {
            const T factor = work(i, col) / work(r, col);
            for (size_t j = col; j < C; ++j) {
                work(i, j) -= factor * work(r, j);
            }
        }
        ++r;
    }
    return r;
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
    int    degree;
    size_t s;
};

template<typename T>
[[nodiscard]] constexpr T expm_log2(T x) {
    return damp::log(x) * damp::numbers::log2e_v<T>;
}

// Unit roundoff u = 2^{-digits}: 2^{-53} for double, 2^{-24} for float.
template<typename T>
[[nodiscard]] constexpr T expm_log2_u() {
    return -static_cast<T>(std::numeric_limits<T>::digits);
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
            return std::numeric_limits<T>::max();
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
    scalar_type_t<T> d4;
    scalar_type_t<T> d6;
    scalar_type_t<T> d8;
    scalar_type_t<T> d10;
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
    const real_t          cap = std::numeric_limits<real_t>::max();
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
    if (eta5 > theta13 && damp::isfinite(eta5)) {
        const real_t lg = expm_log2(eta5 / theta13);
        if (lg > real_t{0}) {
            s = static_cast<size_t>(damp::ceil(lg));
        }
    } else if (!damp::isfinite(eta5)) {
        const real_t nrm = one_norm(A);
        if (nrm > theta13 && damp::isfinite(nrm)) {
            s = static_cast<size_t>(damp::ceil(expm_log2(nrm / theta13)));
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
 * inverse). On a failed solve, or a non-finite input, returns the identity —
 * a total function for the common discretize / integrator call sites.
 *
 * @note Compare with MATLAB®'s expm(A).
 * @see Al-Mohy and Higham, "A New Scaling and Squaring Algorithm for the
 *      Matrix Exponential" (SIAM J. Matrix Anal. Appl., 2009)
 * @see Higham, "Functions of Matrices" (2008), §10.3
 *
 * @tparam T Element type
 * @tparam N Matrix dimension
 * @param A Square matrix
 * @return exp(A), or I if the input is non-finite or the Padé solve fails
 */

template<typename T, size_t N>
[[nodiscard]] constexpr Matrix<N, N, T> expm(const Matrix<N, N, T>& A) {
    if (!damp::isfinite(one_norm(A))) {
        return Matrix<N, N, T>::identity();
    }

    const detail::ExpmTriangle tri = detail::expm_triangle(A);
    if (tri == detail::ExpmTriangle::diagonal) {
        return detail::expm_diagonal(A);
    }

    const detail::ExpmScale               scale = detail::expm_scale(A);
    const Matrix<N, N, T>                 scaled = detail::expm_div_pow2(A, scale.s);
    const damp::optional<Matrix<N, N, T>> R = detail::expm_pade(scaled, scale.degree);
    if (!R) {
        return Matrix<N, N, T>::identity();
    }
    if (tri == detail::ExpmTriangle::none) {
        return detail::expm_square(R.value(), scale.s);
    }
    return detail::expm_square_triangular(R.value(), A, scale.s, tri);
}

/**
 * @brief Matrix square root via Denman–Beavers iteration
 *
 * Computes the principal square root @f$ S = \sqrt{A} @f$ with @f$ S\cdot S = A @f$.
 *
 * Denman–Beavers iteration (inverses are the step deliverables, not
 * intermediate solves of @f$ Ax=b @f$):
 * @f[
 *   Y_{k+1} = \tfrac{1}{2}\bigl(Y_k + Z_k^{-1}\bigr),\quad
 *   Z_{k+1} = \tfrac{1}{2}\bigl(Z_k + Y_k^{-1}\bigr)
 * @f]
 * with @f$ Y_0 = A @f$, @f$ Z_0 = I @f$. Converges to
 * @f$ Y \to \sqrt{A} @f$, @f$ Z \to (\sqrt{A})^{-1} @f$.
 *
 * Returns damp::nullopt if a singular iterate appears or the iteration does not
 * converge (e.g. @f$ A @f$ has a real negative eigenvalue).
 *
 * @note Compare with MATLAB®'s sqrtm(A).
 * @see Higham, "Functions of Matrices" (2008), §6.3
 *
 * @tparam T Element type
 * @tparam N Matrix dimension
 * @param A Square matrix (no real negative eigenvalues)
 * @return Principal √A, or damp::nullopt on failure
 */
template<typename T, size_t N>
[[nodiscard]] constexpr damp::optional<Matrix<N, N, T>> sqrt(const Matrix<N, N, T>& A) {
    using real_t = scalar_type_t<T>;
    Matrix<N, N, T> Y = A;
    Matrix<N, N, T> Z = Matrix<N, N, T>::identity();

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
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                diff += damp::abs(Y_next(i, j) - Y(i, j));
            }
        }

        Y = Y_next;
        Z = Z_next;

        if (diff < default_tol<T>()) {
            return Y;
        }
    }

    return damp::nullopt;
}

/// @brief MATLAB®-style alias for @ref sqrt
template<typename T, size_t N>
[[nodiscard]] constexpr damp::optional<Matrix<N, N, T>> sqrtm(const Matrix<N, N, T>& A) {
    return sqrt(A);
}

/**
 * @brief Principal matrix logarithm via inverse scaling and squaring
 *
 * Computes @f$ X = \log(A) @f$ such that @f$ \exp(X) = A @f$ (principal branch).
 * Scales by repeated square roots until @f$ \|A_{\mathrm{scaled}} - I\|_\infty @f$
 * is small, evaluates the series
 * @f[
 *   \log(I + X) = X - \frac{X^2}{2} + \frac{X^3}{3} - \cdots
 * @f]
 * and multiplies by @f$ 2^s @f$.
 *
 * @p A must be invertible with no eigenvalues on the closed negative real axis.
 * Returns damp::nullopt if a square-root step fails (singular iterate / bad spectrum).
 *
 * @note Compare with MATLAB®'s logm(A).
 * @see Higham, "Functions of Matrices" (2008), §11
 * @see sqrt() for the Denman–Beavers square-root step
 *
 * @tparam T Element type
 * @tparam N Matrix dimension
 * @param A Square matrix
 * @return Principal log(A), or damp::nullopt on failure
 */
template<typename T, size_t N>
[[nodiscard]] constexpr damp::optional<Matrix<N, N, T>> log(const Matrix<N, N, T>& A) {
    Matrix I = Matrix<N, N, T>::identity();
    Matrix A_scaled = A;
    size_t s = 0;

    // A_scaled = A^(1/2^s) until ||A_scaled - I||_∞ < 1/2
    while (infinity_norm(A_scaled - I) > scalar_type_t<T>{0.5} && s < 20) {
        auto S = sqrt(A_scaled);
        if (!S) {
            return damp::nullopt;
        }
        A_scaled = S.value();
        ++s;
    }

    // log(I + X) series, X = A_scaled - I
    Matrix X = A_scaled - I;
    Matrix result = X;
    Matrix X_power = X;

    for (size_t n = 2; n <= 20; ++n) {
        X_power = X_power * X;
        T sign = (n % 2 == 0) ? T{-1} : T{1};
        result = result + X_power * (sign / static_cast<T>(n));
    }

    // log(A) = 2^s * log(A^(1/2^s))
    T scale = static_cast<T>(size_t{1} << s);
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            result(i, j) *= scale;
        }
    }

    return result;
}

/// @brief MATLAB®-style alias for @ref log
template<typename T, size_t N>
[[nodiscard]] constexpr damp::optional<Matrix<N, N, T>> logm(const Matrix<N, N, T>& A) {
    return log(A);
}

/**
 * @brief Integer matrix power via binary exponentiation
 *
 * Computes @f$ A^p @f$. For @f$ p < 0 @f$, returns @f$ (A^{|p|})^{-1} @f$ via
 * Matrix::inverse (the inverse is the deliverable). Returns damp::nullopt
 * if that inverse is singular — no identity fallback.
 *
 * @note Compare with MATLAB®'s mpower(A, p) for integer @p p.
 *
 * @tparam T Element type
 * @tparam N Matrix dimension
 * @param A Square matrix
 * @param p Integer exponent
 * @return Aᵖ, or damp::nullopt if a required inverse is singular
 */
template<typename T, size_t N>
[[nodiscard]] constexpr damp::optional<Matrix<N, N, T>> pow(const Matrix<N, N, T>& A, int p) {
    if (p == 0) {
        return Matrix<N, N, T>::identity();
    }

    bool negate = p < 0;
    if (negate) {
        p = -p;
    }

    Matrix<N, N, T> result = Matrix<N, N, T>::identity();
    Matrix<N, N, T> base = A;

    while (p > 0) {
        if (p & 1) {
            result = result * base;
        }
        base = base * base;
        p >>= 1;
    }

    if (negate) {
        return result.inverse();
    }

    return result;
}

/**
 * @brief Real matrix power Aᵖ via exp(p log(A))
 *
 * Near-integer exponents reduce to the integer @ref pow path. Non-integer
 * exponents require a successful principal logarithm:
 * @f[
 *   A^p = \exp\bigl(p\,\log(A)\bigr)
 * @f]
 *
 * @note Compare with MATLAB®'s mpower(A, p).
 * @see log(), expm()
 *
 * @tparam T Element type
 * @tparam N Matrix dimension
 * @param A Square matrix
 * @param p Real exponent
 * @return Aᵖ, or damp::nullopt if @ref log fails (or integer inverse fails)
 */
template<typename T, size_t N>
[[nodiscard]] constexpr damp::optional<Matrix<N, N, T>> pow(const Matrix<N, N, T>& A, T p) {
    // Integer case (constexpr-safe, no std::modf)
    int p_int = static_cast<int>(p);
    T   p_frac = p - static_cast<T>(p_int);
    if (damp::abs(p_frac) < default_tol<T>()) {
        return pow(A, p_int);
    }

    auto L = log(A);
    if (!L) {
        return damp::nullopt;
    }
    return expm(L.value() * p);
}

/**
 * @brief Compute sin(A) and cos(A) together via scaling and double-angle reconstruction
 *
 * More efficient than calling sin() and cos() separately — computes both
 * with one scaling pass and shared Taylor series evaluation.
 *
 * Scales A down so ||A/2ˢ|| < 0.5 (where the Taylor series converges
 * accurately), then recovers via repeated double-angle formulas:
 *
 *     sin(2A) = 2·sin(A)·cos(A)
 *     cos(2A) = 2·cos²(A) − I
 *
 * Since sin(A) and cos(A) are polynomials in A, they commute with each other,
 * so the scalar double-angle formulas apply directly.
 *
 * @see Higham, "Functions of Matrices" (2008), §12.3
 *
 * @param A Square matrix
 * @return {sin(A), cos(A)}
 */
template<typename T, size_t N>
[[nodiscard]] constexpr damp::pair<Matrix<N, N, T>, Matrix<N, N, T>> sincos(const Matrix<N, N, T>& A) {
    using real_t = scalar_type_t<T>;
    Matrix<N, N, T> I = Matrix<N, N, T>::identity();

    // Scale down until ||A/2^s|| < 0.5
    size_t s = 0;
    real_t scaled_norm = infinity_norm(A);
    while (scaled_norm > real_t{0.5}) {
        scaled_norm *= real_t{0.5};
        ++s;
    }

    T scale = T{1};
    for (size_t i = 0; i < s; ++i) {
        scale *= static_cast<T>(0.5);
    }

    Matrix<N, N, T> As = A * scale;
    Matrix<N, N, T> As2 = As * As;

    // Taylor series for sin(As) = As - As³/3! + As⁵/5! - ...
    Matrix<N, N, T> sinA = As;
    {
        Matrix<N, N, T> A_power = As;
        T               factorial = T{1};
        T               sign = T{-1};
        for (size_t n = 3; n <= 21; n += 2) {
            factorial *= static_cast<T>(n - 1) * static_cast<T>(n);
            A_power = A_power * As2;
            sinA = sinA + A_power * (sign / factorial);
            sign = -sign;
        }
    }

    // Taylor series for cos(As) = I - As²/2! + As⁴/4! - ...
    Matrix<N, N, T> cosA = I;
    {
        Matrix<N, N, T> A_power = I;
        T               factorial = T{1};
        T               sign = T{-1};
        for (size_t n = 2; n <= 20; n += 2) {
            factorial *= static_cast<T>(n - 1) * static_cast<T>(n);
            A_power = A_power * As2;
            cosA = cosA + A_power * (sign / factorial);
            sign = -sign;
        }
    }

    // Double-angle reconstruction: s iterations
    for (size_t i = 0; i < s; ++i) {
        Matrix<N, N, T> new_sin = sinA * cosA * T{2};
        Matrix<N, N, T> new_cos = cosA * cosA * T{2} - I;
        sinA = new_sin;
        cosA = new_cos;
    }

    return {sinA, cosA};
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
 * @note Compare with MATLAB®'s funm(A, @@cos).
 * @see sincos() to compute both sin(A) and cos(A) in one call
 * @see Higham, "Functions of Matrices" (2008), §12.3
 *
 * @param A Square matrix
 * @return Matrix cosine cos(A)
 */
template<typename T, size_t N>
[[nodiscard]] constexpr Matrix<N, N, T> cos(const Matrix<N, N, T>& A) {
    // Unlike sin (whose double-angle step sin(2x)=2·sin·cos needs cos), cos can
    // be reconstructed from cos alone via cos(2x)=2·cos²−1, so this avoids the
    // sin Taylor series and sin doubling that sincos(A).second would discard.
    using real_t = scalar_type_t<T>;
    Matrix<N, N, T> I = Matrix<N, N, T>::identity();

    // Scale down until ||A/2^s|| < 0.5
    size_t s = 0;
    real_t scaled_norm = infinity_norm(A);
    while (scaled_norm > real_t{0.5}) {
        scaled_norm *= real_t{0.5};
        ++s;
    }

    T scale = T{1};
    for (size_t i = 0; i < s; ++i) {
        scale *= static_cast<T>(0.5);
    }

    Matrix<N, N, T> As = A * scale;
    Matrix<N, N, T> As2 = As * As;

    // Taylor series for cos(As) = I - As²/2! + As⁴/4! - ...
    Matrix<N, N, T> cosA = I;
    {
        Matrix<N, N, T> A_power = I;
        T               factorial = T{1};
        T               sign = T{-1};
        for (size_t n = 2; n <= 20; n += 2) {
            factorial *= static_cast<T>(n - 1) * static_cast<T>(n);
            A_power = A_power * As2;
            cosA = cosA + A_power * (sign / factorial);
            sign = -sign;
        }
    }

    // Double-angle reconstruction: cos(2x) = 2·cos²(x) − 1
    for (size_t i = 0; i < s; ++i) {
        cosA = (cosA * cosA * T{2}) - I;
    }

    return cosA;
}

/**
 * @brief Matrix hyperbolic sine sinh(A) = (exp(A) − exp(−A))/2
 *
 * @note Compare with MATLAB®'s funm(A, @@sinh).
 * @see expm(), cosh()
 *
 * @param A Square matrix
 * @return sinh(A)
 */
template<typename T, size_t N>
[[nodiscard]] constexpr Matrix<N, N, T> sinh(const Matrix<N, N, T>& A) {
    Matrix<N, N, T> exp_A = expm(A);
    Matrix<N, N, T> exp_neg_A = expm(A * T{-1});
    return (exp_A - exp_neg_A) * static_cast<T>(0.5);
}

/**
 * @brief Matrix hyperbolic cosine cosh(A) = (exp(A) + exp(−A))/2
 *
 * @note Compare with MATLAB®'s funm(A, @@cosh).
 * @see expm(), sinh()
 *
 * @param A Square matrix
 * @return cosh(A)
 */
template<typename T, size_t N>
[[nodiscard]] constexpr Matrix<N, N, T> cosh(const Matrix<N, N, T>& A) {
    Matrix<N, N, T> exp_A = expm(A);
    Matrix<N, N, T> exp_neg_A = expm(A * T{-1});
    return (exp_A + exp_neg_A) * static_cast<T>(0.5);
}
} // namespace mat
} // namespace damp