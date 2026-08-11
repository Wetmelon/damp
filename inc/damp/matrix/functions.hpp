// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file functions.hpp
 * @brief Matrix functions (expm, sqrtm, funm, ...)
 */

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

/**
 * @brief Matrix exponential via scaling and squaring with Padé approximant of degree 13
 *
 * Computes @f$ \exp(A) = \bigl(\exp(A/2^s)\bigr)^{2^s} @f$ where @f$ s @f$ is
 * chosen so @f$ \|A/2^s\|_\infty @f$ is below the Padé-13 threshold (Higham 2005).
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
 * inverse). On solve failure (pathologically singular @f$ V-U @f$), returns the
 * identity — a total function for the common discretize / integrator call sites.
 * Well-scaled control matrices do not hit this path.
 *
 * @note Compare with MATLAB®'s expm(A).
 * @see Higham, "The Scaling and Squaring Method for the Matrix Exponential"
 *      (SIAM J. Matrix Anal. Appl., 2005)
 * @see Higham, "Functions of Matrices" (2008), §10
 *
 * @tparam T Element type
 * @tparam N Matrix dimension
 * @param A Square matrix
 * @return exp(A), or I if the Padé linear solve fails
 */

template<typename T, size_t N>
[[nodiscard]] constexpr Matrix<N, N, T> expm(const Matrix<N, N, T>& A) {
    using real_t = scalar_type_t<T>;
    Matrix<N, N, T> I = Matrix<N, N, T>::identity();

    // Compute infinity norm (always real)
    const real_t norm = mat::infinity_norm(A);

    // Tiny/nilpotent matrix shortcut (Taylor series)
    if (norm <= default_tol<T>()) {
        Matrix A2 = A * A;
        Matrix A3 = A2 * A;
        Matrix A4 = A3 * A;
        Matrix A5 = A4 * A;
        Matrix A6 = A5 * A;
        return I + A
             + A2 * (T{1} / T{2})
             + A3 * (T{1} / T{6})
             + A4 * (T{1} / T{24})
             + A5 * (T{1} / T{120})
             + A6 * (T{1} / T{720});
    }

    // Scaling for Padé-13
    constexpr real_t theta13 = real_t(2.097847961257068); // from Higham 2005
    size_t           s = 0;
    real_t           scaled_norm = norm;
    while (scaled_norm > theta13) {
        scaled_norm *= real_t(0.5);
        ++s;
    }

    T scale = T(1);
    for (size_t i = 0; i < s; ++i) {
        scale *= T(0.5);
    }

    Matrix<N, N, T> A_scaled = A;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            A_scaled(i, j) *= scale;
        }
    }

    // Precompute powers
    Matrix<N, N, T> A2 = A_scaled * A_scaled;
    Matrix<N, N, T> A4 = A2 * A2;
    Matrix<N, N, T> A6 = A4 * A2;
    Matrix<N, N, T> A8 = A6 * A2;
    Matrix<N, N, T> A10 = A8 * A2;
    Matrix<N, N, T> A12 = A10 * A2;

    // Padé-13 coefficients
    constexpr T b0 = T(64764752532480000.0);
    constexpr T b1 = T(32382376266240000.0);
    constexpr T b2 = T(7771770303897600.0);
    constexpr T b3 = T(1187353796428800.0);
    constexpr T b4 = T(129060195264000.0);
    constexpr T b5 = T(10559470521600.0);
    constexpr T b6 = T(670442572800.0);
    constexpr T b7 = T(33522128640.0);
    constexpr T b8 = T(1323241920.0);
    constexpr T b9 = T(40840800.0);
    constexpr T b10 = T(960960.0);
    constexpr T b11 = T(16380.0);
    constexpr T b12 = T(182.0);
    constexpr T b13 = T(1.0);

    // Compute U and V
    Matrix<N, N, T> U = A_scaled * (b1 * I + b3 * A2 + b5 * A4 + b7 * A6 + b9 * A8 + b11 * A10 + b13 * A12);
    Matrix<N, N, T> V = b0 * I + b2 * A2 + b4 * A4 + b6 * A6 + b8 * A8 + b10 * A10 + b12 * A12;

    // Solve (V-U) * R = V+U
    auto R_opt = solve(V - U, V + U);
    if (!R_opt) {
        // Documented fallback: Padé denominator singular / ill-conditioned.
        return I;
    }
    Matrix<N, N, T> R = R_opt.value();

    // Iterative refinement: solve for residual, apply two refinements
    for (int iter = 0; iter < 2; ++iter) {
        Matrix<N, N, T> residual = (V + U) - (V - U) * R;
        auto            delta_opt = solve(V - U, residual);
        if (!delta_opt) {
            break;
        }
        R = R + delta_opt.value();
    }

    // 7️⃣ Squaring phase
    for (size_t i = 0; i < s; ++i) {
        R = R * R;
    }

    return R;
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