// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file model_inverse.hpp
 * @brief Feedforward from the model's inverse: design::state_feedforward → StateFeedforward (u* from a
 *        state reference) and design::output_inverse → OutputInverse (x* and u* from an output reference)
 *
 * Both run on a discrete model x⁺ = A·x + B·u. StateFeedforward takes a state reference and returns
 * the input that moves the model along it, u*ₖ = B⁺(x*ₖ₊₁ − A·x*ₖ). OutputInverse takes the
 * controlled outputs' reference y* over a fixed preview and returns the state and input of the
 * model's trajectory through it: the reference point (x*, u*) that the LQ laws and MPC take.
 *
 * @code
 * inline constexpr auto kInverse = design::output_inverse<14>(model, damp::array<size_t, 2>{0, 1});
 * static_assert(kInverse.success);
 *
 * Law<kInverse, float> reference;      // OutputInverse
 * const auto p = reference.step(preview); // preview[j] = y* j samples ahead, preview[0] = now
 * const auto u = law.control({}, p.x, p.u);
 * @endcode
 */

#include <cstddef>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/controllers/law_detail.hpp"
#include "damp/design/disturbance.hpp"
#include "damp/design/riccati.hpp"
#include "damp/design/status.hpp"
#include "damp/design/steady_state_target.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/block.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/const_matrix.hpp"
#include "damp/matrix/decomposition.hpp"
#include "damp/matrix/eigen.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/solve.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {

namespace design {

// ---- State reference → input ------------------------------------------------------------

/**
 * @struct StateFeedforwardResult
 * @brief The plant-inversion feedforward of a discrete model: u*ₖ = B⁺·(x*ₖ₊₁ − x*ₖ) − B⁺(A − I)·x*ₖ
 *
 * The two maps give B⁺(x*ₖ₊₁ − A·x*ₖ) without forming A·x*ₖ: a state reference in metres at 0.1 m
 * steps by micrometres per sample, and x*ₖ₊₁ − x*ₖ keeps those digits where x*ₖ₊₁ − A·x*ₖ would
 * cancel them.
 */
template<size_t NX, size_t NU, typename T = double>
struct [[nodiscard]] StateFeedforwardResult {
    Matrix<NU, NX, T> B_pinv{};                            ///< B⁺ = (BᵀB)⁻¹Bᵀ: u* per unit state step x*ₖ₊₁ − x*ₖ
    Matrix<NU, NX, T> drift{};                             ///< B⁺(A − I): u* per unit x*ₖ that the model's own motion (A − I)·x*ₖ takes away
    bool              success{false};                      ///< status == DesignStatus::ok
    DesignStatus      status{DesignStatus::not_converged}; ///< ok, or why the feedforward is refused

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /** @brief The same feedforward in the scalar @p U. */
    template<typename U>
    [[nodiscard]] constexpr StateFeedforwardResult<NX, NU, U> as() const {
        return {B_pinv.template as<U>(), drift.template as<U>(), success, status};
    }
};

/**
 * @brief The plant-inversion feedforward of a discrete model: the input that moves x⁺ = A·x + B·u
 *        along a state reference, u*ₖ the least-squares solution of B·u = x*ₖ₊₁ − A·x*ₖ
 *
 * B⁺ comes from the Householder QR of B, B = QR, as R⁻¹Qᵀ by back substitution: no normal
 * equations and no inverse. When x*ₖ₊₁ − A·x*ₖ lies in the range of B (a reference the model
 * can follow) u* is exact; otherwise it is the input whose step comes nearest in the 2-norm of
 * the state, in the model's units.
 *
 * Refused with DesignStatus::invalid_sample_time for a continuous model (Ts = 0), with
 * DesignStatus::bad_weights when an entry of A or B is not finite, and with
 * DesignStatus::singular when B does not have full column rank (two inputs that act alike).
 *
 * @see "Matrix Computations" (Golub & Van Loan, 2013), §5.3 (the full-rank least-squares problem)
 * @see "Computer-Controlled Systems" (Åström & Wittenmark, 1997): feedforward from the inverse of the process model
 *
 * @param sys Discrete model (Ts > 0)
 * @return StateFeedforwardResult: deploy with `StateFeedforward ff{res}`
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr StateFeedforwardResult<NX, NU, T> state_feedforward(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    using Result = StateFeedforwardResult<NX, NU, T>;
    static_assert(NU <= NX, "state_feedforward: more inputs than states cannot have a full-column-rank B.");
    if (!damp::finite_positive(sys.Ts)) {
        return ::damp::detail::refused<Result>(DesignStatus::invalid_sample_time);
    }
    const auto finite = [](T v) { return damp::isfinite(v); };
    if (!mat::all_of(sys.A, finite) || !mat::all_of(sys.B, finite)) {
        return ::damp::detail::refused<Result>(DesignStatus::bad_weights);
    }
    const auto qr = mat::qr_decompose(sys.B);
    if (!qr.full_rank()) {
        return ::damp::detail::refused<Result>(DesignStatus::singular); // B loses rank: an input does nothing the others cannot
    }
    const auto B_pinv = mat::solve(qr.R.upper_triangle(), Matrix<NU, NX, T>(qr.Q.transpose()));
    if (!B_pinv) {
        return ::damp::detail::refused<Result>(DesignStatus::singular);
    }
    Result out{};
    out.B_pinv = B_pinv.value();
    out.drift = out.B_pinv * (sys.A - Matrix<NX, NX, T>::identity());
    return ::damp::detail::set_status(out, DesignStatus::ok);
}

// ---- Output reference → state and input -------------------------------------------------

/**
 * @struct OutputInverseResult
 * @brief The stable inverse of a discrete model over a preview of NP output samples
 *
 * With d the preview steps stacked, dₗ = y*ₖ₊ₗ − y*ₖ for l = 1 … NP − 1, and s the state of
 * the stable zero dynamics in the model's coordinates (zero at rest),
 * @f[
 *   x^\star_k = X_h y^\star_k + X_p d + s_k,\qquad
 *   u^\star_k = U_h y^\star_k + U_p d + U_s s_k,\qquad
 *   s_{k+1} = S s_k + S_p d .
 * @f]
 * At rest d = 0 and s = 0, so x* = X_h·y* is the steady state of design::steady_state_target.
 *
 * @tparam NX States of the model
 * @tparam NU Inputs, and as many controlled outputs
 * @tparam NP Preview samples y*ₖ … y*ₖ₊NP₋₁
 */
template<size_t NX, size_t NU, size_t NP, typename T = double>
struct [[nodiscard]] OutputInverseResult {
    static_assert(NP >= 2, "output_inverse: the preview holds the reference now and at least one sample ahead.");
    static constexpr size_t ND = NU * (NP - 1); ///< Preview steps y*ₖ₊ₗ − y*ₖ, l = 1 … NP − 1, stacked

    Matrix<NX, NU, T>       x_hold{};                            ///< X_h: x* per unit y*ₖ, the steady state that holds the outputs there
    Matrix<NU, NU, T>       u_hold{};                            ///< U_h: u* per unit y*ₖ, the input that holds it
    Matrix<NX, ND, T>       x_preview{};                         ///< X_p: x* per unit preview step
    Matrix<NU, ND, T>       u_preview{};                         ///< U_p: u* per unit preview step
    Matrix<NX, NX, T>       zero_dynamics{};                     ///< S: the stable zero dynamics in the model's coordinates
    Matrix<NX, ND, T>       zero_preview{};                      ///< S_p: the zero dynamics' input per unit preview step
    Matrix<NU, NX, T>       u_zero{};                            ///< U_s: u* per unit zero-dynamics state (x* takes it as is)
    damp::array<size_t, NU> relative_degree{};                   ///< rᵢ: samples from an input to controlled output i
    size_t                  unstable_zeros{0};                   ///< Zeros outside the unit circle, solved over the preview
    T                       truncation{0};                       ///< Weight of the first sample past the preview relative to the first one in it, 1/|z|^(NP − max rᵢ) for the zero z nearest the unit circle from outside (0: none)
    bool                    success{false};                      ///< status == DesignStatus::ok
    DesignStatus            status{DesignStatus::not_converged}; ///< ok, or why the inverse is refused

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    /** @brief The same inverse in the scalar @p U. */
    template<typename U>
    [[nodiscard]] constexpr OutputInverseResult<NX, NU, NP, U> as() const {
        return {
            x_hold.template as<U>(),
            u_hold.template as<U>(),
            x_preview.template as<U>(),
            u_preview.template as<U>(),
            zero_dynamics.template as<U>(),
            zero_preview.template as<U>(),
            u_zero.template as<U>(),
            relative_degree,
            unstable_zeros,
            static_cast<U>(truncation),
            success,
            status,
        };
    }
};

namespace detail {

/** @brief |M| entry by entry: the size of each term a product of absolute values bounds. */
template<size_t R, size_t C, typename T>
[[nodiscard]] constexpr Matrix<R, C, T> entrywise_abs(const Matrix<R, C, T>& M) {
    Matrix<R, C, T> out{};
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            out(i, j) = damp::abs(M(i, j));
        }
    }
    return out;
}

/** @brief A·M summed over the inner indices [@p j0, @p j1) only, for an M whose other rows are zero. */
template<size_t R, size_t K, size_t C, typename T>
[[nodiscard]] constexpr Matrix<R, C, T> partial_product(const Matrix<R, K, T>& A, const Matrix<K, C, T>& M, size_t j0, size_t j1) {
    Matrix<R, C, T> out{};
    for (size_t i = 0; i < R; ++i) {
        for (size_t c = 0; c < C; ++c) {
            T sum = T{0};
            for (size_t j = j0; j < j1; ++j) {
                sum += A(i, j) * M(j, c);
            }
            out(i, c) = sum;
        }
    }
    return out;
}

/** @brief The relative degree of each controlled output and the decoupling matrix: output i at sample k + rᵢ is C*ᵢ·xₖ + D*ᵢ·uₖ. */
template<size_t NX, size_t NU, typename T>
struct Decoupling {
    damp::array<size_t, NU> r{};                            ///< rᵢ per controlled output
    Matrix<NU, NU, T>       D{};                            ///< D*: row i is Cᵢ·A^(rᵢ−1)·B (Dᵢ when rᵢ = 0)
    Matrix<NU, NX, T>       C{};                            ///< C*: row i is Cᵢ·A^rᵢ
    DesignStatus            status{DesignStatus::singular}; ///< ok, or singular
};

/**
 * @brief The first nonzero Markov parameter of each controlled output, and the decoupling matrix
 *
 * A Markov row m = Cᵢ·A^(j−1)·B counts as zero when every entry is within 16(j + 1)·n·ε of the
 * same product taken in absolute values, |Cᵢ|·|A|^(j−1)·|B|. That bound does not change when
 * the states are rescaled, so the relative degree does not depend on their units. An output no
 * input reaches within n samples is refused, and so is a decoupling matrix that loses rank once
 * its rows are scaled to unit size.
 */
template<size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr Decoupling<NX, NU, T> decoupling(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NU, NX, T>& Cz,
    const Matrix<NU, NU, T>& Dz
) {
    Decoupling<NX, NU, T>   out{};
    const Matrix<NX, NX, T> A_abs = entrywise_abs(A);
    const Matrix<NX, NU, T> B_abs = entrywise_abs(B);
    for (size_t i = 0; i < NU; ++i) {
        bool direct = false;
        for (size_t c = 0; c < NU; ++c) {
            direct = direct || (Dz(i, c) != T{0});
        }
        if (direct) {
            out.r[i] = 0;
            out.D.template block<1, NU>(i, 0) = Dz.template block<1, NU>(i, 0);
            out.C.template block<1, NX>(i, 0) = Cz.template block<1, NX>(i, 0);
            continue;
        }
        Matrix<1, NX, T> row = Cz.template block<1, NX>(i, 0);
        Matrix<1, NX, T> row_abs = entrywise_abs(row);
        bool             found = false;
        for (size_t j = 1; j <= NX && !found; ++j) {
            const Matrix<1, NU, T> m = row * B;
            const Matrix<1, NU, T> bound = row_abs * B_abs;
            const T                tol = T{16} * static_cast<T>((j + 1) * NX) * eps<T>();
            for (size_t c = 0; c < NU; ++c) {
                found = found || (damp::abs(m(0, c)) > tol * bound(0, c));
            }
            if (found) {
                out.r[i] = j;
                out.D.template block<1, NU>(i, 0) = m;
                out.C.template block<1, NX>(i, 0) = row * A;
            } else {
                row = row * A;
                row_abs = row_abs * A_abs;
            }
        }
        if (!found) {
            return out; // no input reaches this output: singular
        }
    }
    Matrix<NU, NU, T> unit_rows = out.D;
    for (size_t i = 0; i < NU; ++i) {
        T peak = T{0};
        for (size_t c = 0; c < NU; ++c) {
            peak = damp::max(peak, damp::abs(out.D(i, c)));
        }
        for (size_t c = 0; c < NU; ++c) {
            unit_rows(i, c) = out.D(i, c) / peak;
        }
    }
    if (mat::qr_decompose(unit_rows).full_rank()) {
        out.status = DesignStatus::ok;
    }
    return out;
}

/** @brief A real Schur form of the balanced matrix, eigenvalues inside the unit circle first: A = V·T·V⁻¹, V = D·Z. */
template<size_t NX, typename T>
struct DiskSchur {
    Matrix<NX, NX, T>    T_form{};                            ///< Quasi-triangular T: |λ| < 1 in the leading `stable` rows
    Matrix<NX, NX, T>    Z{};                                 ///< Orthogonal Schur vectors of the balanced matrix
    damp::array<int, NX> e{};                                 ///< Balancing exponents: D = diag(2^e)
    size_t               stable{0};                           ///< Eigenvalues with |λ| < 1, counting a pair as two
    T                    nearest_outside{0};                  ///< Smallest |λ| among the others (0: none)
    DesignStatus         status{DesignStatus::not_converged}; ///< ok, singular (an eigenvalue on the unit circle) or not_converged

    /** @brief V = D·Z: Schur coordinates to the model's, column j the model-coordinate image of ζⱼ. */
    [[nodiscard]] constexpr Matrix<NX, NX, T> basis() const {
        Matrix<NX, NX, T> out{};
        for (size_t i = 0; i < NX; ++i) {
            for (size_t j = 0; j < NX; ++j) {
                out(i, j) = mat::detail::scale_by_pow2(Z(i, j), e[i]);
            }
        }
        return out;
    }

    /** @brief V⁻¹ = Zᵀ·D⁻¹, exact: Z is orthogonal and D a diagonal of powers of two. */
    [[nodiscard]] constexpr Matrix<NX, NX, T> inverse_basis() const {
        Matrix<NX, NX, T> out{};
        for (size_t i = 0; i < NX; ++i) {
            for (size_t j = 0; j < NX; ++j) {
                out(j, i) = mat::detail::scale_by_pow2(Z(i, j), -e[i]);
            }
        }
        return out;
    }

    /** @brief V⁻¹·M = Zᵀ·D⁻¹·M: from the model's coordinates to Schur coordinates. */
    template<size_t NC>
    [[nodiscard]] constexpr Matrix<NX, NC, T> to_schur(const Matrix<NX, NC, T>& M) const {
        Matrix<NX, NC, T> scaled = M;
        for (size_t i = 0; i < NX; ++i) {
            for (size_t c = 0; c < NC; ++c) {
                scaled(i, c) = mat::detail::scale_by_pow2(M(i, c), -e[i]);
            }
        }
        return Z.transpose() * scaled;
    }
};

/**
 * @brief The real Schur form of @p A with the eigenvalues inside the unit circle leading
 *
 * Balanced first (powers of two, exact), so the split does not depend on the units of the
 * states. An eigenvalue within √ε of the unit circle has neither a causal nor an anticausal
 * bounded solution and is refused as DesignStatus::singular.
 *
 * @see "Matrix Computations" (Golub & Van Loan, 2013), §7.6 (ordered Schur form)
 */
template<size_t NX, typename T>
[[nodiscard]] constexpr DiskSchur<NX, T> disk_ordered_schur(const Matrix<NX, NX, T>& A) {
    DiskSchur<NX, T> out{};
    out.e = mat::detail::balance_exponents(A);
    out.T_form = mat::detail::balanced(A, out.e);
    out.Z = Matrix<NX, NX, T>::identity();
    damp::array<T, NX> wr{};
    damp::array<T, NX> wi{};
    if constexpr (NX > 1) {
        mat::detail::hessenberg_reduce(out.T_form, out.Z);
        if (!mat::detail::francis_qr(out.T_form, out.Z, wr, wi)) {
            return out; // not_converged
        }
    } else {
        wr[0] = out.T_form(0, 0);
    }
    const T band = damp::sqrt(eps<T>());
    size_t  inside = 0;
    for (size_t i = 0; i < NX; ++i) {
        const T modulus = damp::sqrt((wr[i] * wr[i]) + (wi[i] * wi[i]));
        if (damp::abs(modulus - T{1}) <= band) {
            out.status = DesignStatus::singular; // a zero on the unit circle: no bounded inverse
            return out;
        }
        if (modulus < T{1}) {
            ++inside;
        } else if (out.nearest_outside == T{0} || modulus < out.nearest_outside) {
            out.nearest_outside = modulus;
        }
    }
    out.stable = ::damp::detail::reorder_schur(out.T_form, out.Z, [](const damp::complex<T>& lambda) { return lambda.abs() < T{1}; });
    if (out.stable == inside) {
        out.status = DesignStatus::ok;
    }
    return out;
}

/**
 * @brief Solve (T₂₂ − σI)·X₂ = X₂ in place on the rows from @p first on of a real Schur form, by
 *        back substitution over its 1×1 and 2×2 diagonal blocks
 *
 * The rows above @p first are not read or written. A 2×2 block is solved by Cramer's rule; the
 * form is standardized, so its determinant |λ − σ|² is that of a complex pair. False when a
 * diagonal block is singular.
 */
template<size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr bool trailing_solve(const Matrix<NX, NX, T>& Tf, size_t first, T shift, Matrix<NX, NU, T>& X) {
    size_t end = NX; // rows [end, NX) are solved
    while (end > first) {
        const bool   pair = (end - first >= 2) && ::damp::detail::schur_subdiagonal(Tf, end - 2);
        const size_t top = pair ? end - 2 : end - 1;
        for (size_t c = 0; c < NU; ++c) {
            for (size_t i = top; i < end; ++i) {
                for (size_t j = end; j < NX; ++j) {
                    X(i, c) -= Tf(i, j) * X(j, c);
                }
            }
        }
        if (pair) {
            const T a = Tf(top, top) - shift;
            const T b = Tf(top, top + 1);
            const T c2 = Tf(top + 1, top);
            const T d = Tf(top + 1, top + 1) - shift;
            const T det = (a * d) - (b * c2);
            if (!(damp::abs(det) > T{0})) {
                return false;
            }
            for (size_t c = 0; c < NU; ++c) {
                const T r0 = X(top, c);
                const T r1 = X(top + 1, c);
                X(top, c) = ((d * r0) - (b * r1)) / det;
                X(top + 1, c) = ((a * r1) - (c2 * r0)) / det;
            }
        } else {
            const T a = Tf(top, top) - shift;
            if (!(damp::abs(a) > T{0})) {
                return false;
            }
            for (size_t c = 0; c < NU; ++c) {
                X(top, c) = X(top, c) / a;
            }
        }
        end = top;
    }
    return true;
}

/**
 * @brief The anticausal part of the inverse in Schur coordinates: ζ₂ per unit y*ₖ₊ₗ, one NX × NU
 *        block per preview sample (rows below `stable` only)
 *
 * The unstable block T₂₂ runs backwards: ζ₂,ₖ = −Σⱼ T₂₂^−(j+1)·B₂·ŷₖ₊ⱼ with ŷᵢ,ₖ₊ⱼ = y*ᵢ,ₖ₊ⱼ₊rᵢ.
 * The terms the preview holds are exact; past it the reference is held at its last sample,
 * whose geometric tail is −T₂₂^−j₀·(T₂₂ − I)⁻¹·B₂ with j₀ = NP − 1 − rᵢ. Each power is a back
 * substitution on T₂₂ alone (trailing_solve), so no inverse is formed and the stable rows cost
 * nothing. nullopt if a diagonal block is singular.
 */
template<size_t NP, size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr damp::optional<damp::array<Matrix<NX, NU, T>, NP>> anticausal_preview(
    const DiskSchur<NX, T>&        sch,
    const Matrix<NX, NU, T>&       B_schur,
    const damp::array<size_t, NU>& r
) {
    Matrix<NX, NU, T> power{}; // T₂₂^−(j+1)·B₂
    for (size_t row = sch.stable; row < NX; ++row) {
        for (size_t c = 0; c < NU; ++c) {
            power(row, c) = B_schur(row, c);
        }
    }
    Matrix<NX, NU, T> tail = power; // T₂₂^−j·(T₂₂ − I)⁻¹·B₂
    if (!trailing_solve(sch.T_form, sch.stable, T{0}, power) || !trailing_solve(sch.T_form, sch.stable, T{1}, tail)) {
        return damp::nullopt;
    }
    damp::array<Matrix<NX, NU, T>, NP> coeff{};
    for (size_t j = 0; j < NP; ++j) {
        for (size_t i = 0; i < NU; ++i) {
            const size_t j0 = NP - 1 - r[i];
            for (size_t row = sch.stable; row < NX; ++row) {
                if (j < j0) {
                    coeff[j + r[i]](row, i) -= power(row, i);
                } else if (j == j0) {
                    coeff[NP - 1](row, i) -= tail(row, i);
                }
            }
        }
        if (!trailing_solve(sch.T_form, sch.stable, T{0}, power) || !trailing_solve(sch.T_form, sch.stable, T{0}, tail)) {
            return damp::nullopt;
        }
    }
    return coeff;
}

/**
 * @brief The steady-state certificate: [X_h; U_h] solves the target system
 *        [A − I, B; C_z, D_z]·[X_h; U_h] = [0; I] with a normwise backward error below √ε
 *
 * Per unit output, the residual r of column c is checked against ‖M‖·‖z_c‖ + 1 (∞-norms, M the
 * target matrix, z_c the column). An entrywise test of the state rows would compare roundoff with
 * roundoff on a plant with an integrator, whose steady state leaves whole terms exactly zero. The
 * output rows C_z·X_h + D_z·U_h = I are also checked entry by entry, against √ε·(sᵢ + 1) with sᵢ
 * the largest entry of row i of |C_z||X_h| + |D_z||U_h|: a large ‖A‖ in the normwise bound would
 * otherwise pass an output that misses its reference at rest by percents.
 */
template<size_t NX, size_t NU, typename T>
[[nodiscard]] constexpr bool holds_steady_state(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NU, NX, T>& Cz,
    const Matrix<NU, NU, T>& Dz,
    const Matrix<NX, NU, T>& x_hold,
    const Matrix<NU, NU, T>& u_hold
) {
    Matrix<NX + NU, NX + NU, T> M{};
    M.template block<NX, NX>(0, 0) = A - Matrix<NX, NX, T>::identity();
    M.template block<NX, NU>(0, NX) = B;
    M.template block<NU, NX>(NX, 0) = Cz;
    M.template block<NU, NU>(NX, NX) = Dz;
    Matrix<NX + NU, NU, T> z{};
    z.template block<NX, NU>(0, 0) = x_hold;
    z.template block<NU, NU>(NX, 0) = u_hold;
    Matrix<NX + NU, NU, T> residual = M * z;
    for (size_t i = 0; i < NU; ++i) {
        residual(NX + i, i) -= T{1};
    }
    T M_norm = T{0};
    for (size_t i = 0; i < NX + NU; ++i) {
        T row = T{0};
        for (size_t j = 0; j < NX + NU; ++j) {
            row += damp::abs(M(i, j));
        }
        M_norm = damp::max(M_norm, row);
    }
    const T tol = damp::sqrt(eps<T>());
    for (size_t c = 0; c < NU; ++c) {
        T z_norm = T{0};
        T r_norm = T{0};
        for (size_t i = 0; i < NX + NU; ++i) {
            z_norm = damp::max(z_norm, damp::abs(z(i, c)));
            r_norm = damp::max(r_norm, damp::abs(residual(i, c)));
        }
        if (!(r_norm <= tol * ((M_norm * z_norm) + T{1}))) {
            return false;
        }
    }
    const Matrix<NU, NU, T> out_scale = (entrywise_abs(Cz) * entrywise_abs(x_hold)) + (entrywise_abs(Dz) * entrywise_abs(u_hold));
    for (size_t i = 0; i < NU; ++i) {
        T row_scale = T{0};
        for (size_t c = 0; c < NU; ++c) {
            row_scale = damp::max(row_scale, out_scale(i, c));
        }
        for (size_t c = 0; c < NU; ++c) {
            if (!(damp::abs(residual(NX + i, c)) <= tol * (row_scale + T{1}))) {
                return false;
            }
        }
    }
    return true;
}

/**
 * @brief The runtime maps of the inverse from its parts, in the model's coordinates and on the
 *        preview's steps, certified on the steady state
 *
 * The causal part runs ζ₁ ← T₁₁·ζ₁ + T₁₂·ζ₂ + B₁·ŷ; in the model's coordinates its state is
 * s = V·[ζ₁; 0] with S = V·blkdiag(T₁₁, 0)·V⁻¹. Taking s about its value at rest, Σ·y*ₖ with
 * Σ = (I − S)⁻¹·Σₗ Sₗ, leaves a state driven by the steps y*ₖ₊ₗ − y*ₖ alone.
 */
template<size_t NP, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr OutputInverseResult<NX, NU, NP, T> assemble_inverse(
    const StateSpace<NX, NU, NY, T, NW, NV>&  sys,
    const Matrix<NU, NX, T>&                  Cz,
    const Matrix<NU, NU, T>&                  Dz,
    const Decoupling<NX, NU, T>&              dec,
    const Matrix<NU, NX, T>&                  K_dec,
    const DiskSchur<NX, T>&                   sch,
    const Matrix<NX, NU, T>&                  B_schur,
    const damp::array<Matrix<NX, NU, T>, NP>& anti
) {
    using Result = OutputInverseResult<NX, NU, NP, T>;
    Result                  out{};
    const size_t            s = sch.stable;
    const Matrix<NX, NX, T> V = sch.basis();
    Matrix<NX, NX, T>       T_lead{}; // blkdiag(T₁₁, 0)
    for (size_t i = 0; i < s; ++i) {
        for (size_t j = 0; j < s; ++j) {
            T_lead(i, j) = sch.T_form(i, j);
        }
    }
    out.zero_dynamics = partial_product(V, partial_product(T_lead, sch.inverse_basis(), 0, s), 0, s);
    Matrix<NX, NU, T> Gs_sum{};
    Matrix<NX, NU, T> Gx_sum{};
    Matrix<NU, NU, T> Gu_sum{};
    for (size_t l = 0; l < NP; ++l) {
        //! ζ₂ per unit y*ₖ₊ₗ (anti[l]) has only unstable rows, ζ₁'s input only stable ones.
        Matrix<NU, NU, T> lead{}; // ŷ per unit y*ₖ₊ₗ: the outputs whose relative degree is l
        bool              leads = false;
        for (size_t i = 0; i < NU; ++i) {
            if (dec.r[i] == l) {
                lead(i, i) = T{1};
                leads = true;
            }
        }
        Matrix<NX, NU, T> causal = partial_product(sch.T_form, anti[l], s, NX); // T₁₂·ζ₂ (+ T₂₂·ζ₂ below s, dropped)
        for (size_t i = 0; i < NX; ++i) {
            for (size_t c = 0; c < NU; ++c) {
                causal(i, c) = (i < s) ? causal(i, c) + (B_schur(i, c) * lead(c, c)) : T{0};
            }
        }
        const Matrix<NX, NU, T> Gs = partial_product(V, causal, 0, s);
        const Matrix<NX, NU, T> Gx = partial_product(V, anti[l], s, NX);
        Matrix<NU, NU, T>       Gu = -(K_dec * Gx);
        if (leads) {
            const auto lead_u = mat::lu_solve(dec.D, lead);
            if (!lead_u) {
                return ::damp::detail::refused<Result>(DesignStatus::singular);
            }
            Gu += lead_u.value();
        }
        Gs_sum += Gs;
        Gx_sum += Gx;
        Gu_sum += Gu;
        if (l > 0) {
            out.x_preview.template block<NX, NU>(0, NU * (l - 1)) = Gx;
            out.u_preview.template block<NU, NU>(0, NU * (l - 1)) = Gu;
            out.zero_preview.template block<NX, NU>(0, NU * (l - 1)) = Gs;
        }
    }
    const auto rest = mat::lu_solve(Matrix<NX, NX, T>(Matrix<NX, NX, T>::identity() - out.zero_dynamics), Gs_sum); // Σ
    if (!rest) {
        return ::damp::detail::refused<Result>(DesignStatus::singular);
    }
    out.x_hold = Gx_sum + rest.value();
    out.u_hold = Gu_sum - (K_dec * rest.value());
    out.u_zero = -K_dec;
    const Matrix<NX, NU, T> first = out.zero_preview.template block<NX, NU>(0, 0);
    out.zero_preview.template block<NX, NU>(0, 0) = first - rest.value(); // s about rest: −Σ·(y*ₖ₊₁ − y*ₖ)
    size_t r_max = 0;
    for (size_t i = 0; i < NU; ++i) {
        out.relative_degree[i] = dec.r[i];
        r_max = damp::max(r_max, dec.r[i]);
    }
    out.unstable_zeros = NX - sch.stable;
    if (out.unstable_zeros > 0) {
        out.truncation = T{1};
        for (size_t j = r_max; j < NP; ++j) {
            out.truncation /= sch.nearest_outside;
        }
    }
    if (!holds_steady_state(sys.A, sys.B, Cz, Dz, out.x_hold, out.u_hold)) {
        return ::damp::detail::refused<Result>(DesignStatus::not_converged);
    }
    return ::damp::detail::set_status(out, DesignStatus::ok);
}

} // namespace detail

/**
 * @brief The stable inverse of a discrete model: the state and input of its trajectory through
 *        a reference y* on the controlled outputs, from y*ₖ … y*ₖ₊NP₋₁
 *
 * For each controlled output i the relative degree rᵢ is the first sample at which an input
 * reaches it, and the decoupling matrix D* (rows Cᵢ·A^(rᵢ−1)·B) must be nonsingular. Then
 * z_{i,k+rᵢ} = C*ᵢ·xₖ + D*ᵢ·uₖ gives uₖ = D*⁻¹(ŷₖ − C*·xₖ), and every trajectory through y* obeys
 * @f[
 *   x_{k+1} = A_z x_k + B D^{\star-1} \hat y_k,\qquad A_z = A - B D^{\star-1} C^\star ,
 * @f]
 * whose eigenvalues are the model's zeros and Σrᵢ at 0. The ordered Schur form splits them at
 * the unit circle. The zeros inside run forward as the runtime's state, so a minimum-phase model
 * needs NP = max rᵢ + 1. The zeros outside (sampling zeros of a fast-sampled plant of relative
 * degree three or more, a non-minimum-phase plant) have one bounded solution, which runs
 * backwards from the future: the preview supplies NP − rᵢ of its terms exactly, and the
 * reference held at its last sample supplies the rest. `truncation` reports how much the held
 * part weighs; a longer preview shrinks it geometrically. At rest the inverse is exact for any NP.
 *
 * The runtime works on the preview's steps y*ₖ₊ₗ − y*ₖ, so in float it spends its digits on
 * the motion, not on the absolute position.
 *
 * Refused with DesignStatus::invalid_sample_time for a continuous model; DesignStatus::bad_weights
 * for a controlled output out of range or listed twice; DesignStatus::singular when an output no
 * input reaches, a singular decoupling matrix (the model has no vector relative degree for these
 * outputs, or is not invertible), or a zero within √ε of the unit circle (the zero at −1 of a
 * sampled double integrator, a zero at 1); DesignStatus::invalid_parameter when NP ≤ max rᵢ;
 * DesignStatus::not_converged when the Schur form does not converge or the inverse fails its
 * steady-state certificate.
 *
 * @see "Decoupling in the design and synthesis of multivariable control systems" (Falb &
 *      Wolovich, 1967), IEEE TAC 12(6): the relative degrees and the decoupling matrix
 * @see "Inversion of multivariable linear systems" (Silverman, 1969), IEEE TAC 14(3)
 * @see "Nonlinear inversion-based output tracking" (Devasia, Chen & Paden, 1996), IEEE TAC 41(7):
 *      the stable inverse, causal on the stable zeros and anticausal on the unstable ones
 * @see "Preview-based stable-inversion for output tracking of linear systems" (Zou & Devasia,
 *      1999), ASME J. Dyn. Sys., Meas., Control 121(4): the finite preview
 * @see "Zeros of sampled systems" (Åström, Hagander & Sternby, 1984), Automatica 20(1): the
 *      sampling zeros outside the unit circle
 *
 * @tparam NP Preview samples y*ₖ … y*ₖ₊NP₋₁ (at least max rᵢ + 1)
 * @param sys     Discrete model (Ts > 0)
 * @param outputs Rows of C (and D) that are the controlled outputs, as many as inputs, each once
 * @return OutputInverseResult: deploy with `OutputInverse inverse{res}` or `Law<kRes, float>`
 */
template<size_t NP, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV, size_t NS>
[[nodiscard]] constexpr OutputInverseResult<NX, NU, NP, T> output_inverse(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const damp::array<size_t, NS>&           outputs
) {
    static_assert(((NS == 0) ? NY : NS) == NU, "output_inverse: list as many controlled outputs as the model has inputs.");
    using Result = OutputInverseResult<NX, NU, NP, T>;
    if (!damp::finite_positive(sys.Ts)) {
        return ::damp::detail::refused<Result>(DesignStatus::invalid_sample_time);
    }
    const auto rows = detail::selected_outputs<NY>(outputs);
    if (!rows) {
        return ::damp::detail::refused<Result>(DesignStatus::bad_weights);
    }
    Matrix<NU, NX, T> Cz{};
    Matrix<NU, NU, T> Dz{};
    for (size_t i = 0; i < NU; ++i) {
        Cz.template block<1, NX>(i, 0) = sys.C.template block<1, NX>((*rows)[i], 0);
        Dz.template block<1, NU>(i, 0) = sys.D.template block<1, NU>((*rows)[i], 0);
    }

    //! The decoupling structure, then the inverse dynamics A_z and its input B·D*⁻¹.
    const auto dec = detail::decoupling(sys.A, sys.B, Cz, Dz);
    if (dec.status != DesignStatus::ok) {
        return ::damp::detail::refused<Result>(dec.status);
    }
    size_t r_max = 0;
    for (size_t i = 0; i < NU; ++i) {
        r_max = damp::max(r_max, dec.r[i]);
    }
    if (r_max + 1 > NP) {
        return ::damp::detail::refused<Result>(DesignStatus::invalid_parameter); // the preview does not reach the inputs' effect
    }
    const auto K_dec = mat::lu_solve(dec.D, dec.C);                        // D*⁻¹·C*
    const auto Bz_t = mat::lu_solve(dec.D.transpose(), sys.B.transpose()); // (B·D*⁻¹)ᵀ
    if (!K_dec || !Bz_t) {
        return ::damp::detail::refused<Result>(DesignStatus::singular);
    }
    //! A − B·D*⁻¹C* cancels whole entries to roundoff; balancing would scale that noise up to the
    //! size of the real entries and spread the exponents over decades, so it is zeroed first, by
    //! the entrywise bound of the product (as decoupling() bounds C*·B).
    Matrix<NX, NX, T>       Az = sys.A - (sys.B * K_dec.value());
    const Matrix<NX, NX, T> Az_bound = detail::entrywise_abs(sys.A) + (detail::entrywise_abs(sys.B) * detail::entrywise_abs(K_dec.value()));
    const T                 cancel_tol = T{16} * static_cast<T>(NX) * eps<T>();
    for (size_t i = 0; i < NX; ++i) {
        for (size_t j = 0; j < NX; ++j) {
            if (damp::abs(Az(i, j)) <= cancel_tol * Az_bound(i, j)) {
                Az(i, j) = T{0};
            }
        }
    }
    const auto sch = detail::disk_ordered_schur(Az);
    if (sch.status != DesignStatus::ok) {
        return ::damp::detail::refused<Result>(sch.status);
    }
    const Matrix<NX, NU, T> B_schur = sch.to_schur(Matrix<NX, NU, T>(Bz_t->transpose()));
    const auto              anti = detail::anticausal_preview<NP>(sch, B_schur, dec.r);
    if (!anti) {
        return ::damp::detail::refused<Result>(DesignStatus::singular);
    }
    return detail::assemble_inverse<NP>(sys, Cz, Dz, dec, K_dec.value(), sch, B_schur, anti.value());
}

/** @brief design::output_inverse with every output of @p sys controlled (NY = NU). */
template<size_t NP, size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr OutputInverseResult<NX, NU, NP, T> output_inverse(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    return output_inverse<NP>(sys, damp::array<size_t, 0>{});
}

} // namespace design

namespace detail {

/** @brief Where a StateFeedforward keeps its maps: in the object, any value at run time. */
template<size_t NX, size_t NU, typename T>
struct HeldFeedforward {
    static constexpr bool constant = false; ///< The matrices are members
    using B_pinv = Matrix<NU, NX, T>;       ///< B⁺
    using drift = Matrix<NU, NX, T>;        ///< B⁺(A − I)
};

/** @brief A StateFeedforward's maps as constants of the type. */
template<const auto& Res, typename T>
struct ConstantFeedforward {
    static constexpr bool constant = true;               ///< The matrices are constants of the type
    static constexpr auto design = Res.template as<T>(); ///< The design, in the runtime's scalar
    static constexpr auto kB_pinv = design.B_pinv;       ///< B⁺
    static constexpr auto kDrift = design.drift;         ///< B⁺(A − I)
    using B_pinv = ConstMatrix<kB_pinv>;                 ///< B⁺
    using drift = ConstMatrix<kDrift>;                   ///< B⁺(A − I)
};

/** @brief Where an OutputInverse keeps its maps: in the object, any value at run time. */
template<size_t NX, size_t NU, size_t NP, typename T>
struct HeldOutputInverse {
    static constexpr bool   constant = false;   ///< The matrices are members
    static constexpr size_t ND = NU * (NP - 1); ///< Preview steps, stacked
    using x_hold = Matrix<NX, NU, T>;           ///< x* per unit y*ₖ
    using u_hold = Matrix<NU, NU, T>;           ///< u* per unit y*ₖ
    using x_preview = Matrix<NX, ND, T>;        ///< x* per unit preview step
    using u_preview = Matrix<NU, ND, T>;        ///< u* per unit preview step
    using zero_dynamics = Matrix<NX, NX, T>;    ///< The stable zero dynamics
    using zero_preview = Matrix<NX, ND, T>;     ///< Their input per unit preview step
    using u_zero = Matrix<NU, NX, T>;           ///< u* per unit zero-dynamics state
};

/** @brief An OutputInverse's maps as constants of the type. */
template<const auto& Res, typename T>
struct ConstantOutputInverse {
    static constexpr bool constant = true;                       ///< The matrices are constants of the type
    static constexpr auto design = Res.template as<T>();         ///< The design, in the runtime's scalar
    static constexpr auto kx_hold = design.x_hold;               ///< x* per unit y*ₖ
    static constexpr auto ku_hold = design.u_hold;               ///< u* per unit y*ₖ
    static constexpr auto kx_preview = design.x_preview;         ///< x* per unit preview step
    static constexpr auto ku_preview = design.u_preview;         ///< u* per unit preview step
    static constexpr auto kzero_dynamics = design.zero_dynamics; ///< The stable zero dynamics
    static constexpr auto kzero_preview = design.zero_preview;   ///< Their input per unit preview step
    static constexpr auto ku_zero = design.u_zero;               ///< u* per unit zero-dynamics state
    using x_hold = ConstMatrix<kx_hold>;                         ///< x* per unit y*ₖ
    using u_hold = ConstMatrix<ku_hold>;                         ///< u* per unit y*ₖ
    using x_preview = ConstMatrix<kx_preview>;                   ///< x* per unit preview step
    using u_preview = ConstMatrix<ku_preview>;                   ///< u* per unit preview step
    using zero_dynamics = ConstMatrix<kzero_dynamics>;           ///< The stable zero dynamics
    using zero_preview = ConstMatrix<kzero_preview>;             ///< Their input per unit preview step
    using u_zero = ConstMatrix<ku_zero>;                         ///< u* per unit zero-dynamics state
};

} // namespace detail

/**
 * @brief The plant-inversion feedforward at run time: u* = B⁺·(x*ₖ₊₁ − x*ₖ) − B⁺(A − I)·x*ₖ, two
 *        matrix–vector products
 *
 * Runs design::state_feedforward. It holds no state.
 */
template<size_t NX, size_t NU, typename T = float, typename Store = detail::HeldFeedforward<NX, NU, T>>
struct StateFeedforward {
    [[no_unique_address]] typename Store::B_pinv B_pinv{}; ///< B⁺: u* per unit state step
    [[no_unique_address]] typename Store::drift  drift{};  ///< B⁺(A − I): u* per unit state

    constexpr StateFeedforward() = default;

    /** @brief From design::state_feedforward (any scalar). A failed design does not compile when constant-evaluated. */
    template<typename U>
        requires(!Store::constant)
    constexpr StateFeedforward(const design::StateFeedforwardResult<NX, NU, U>& result) // NOLINT
        : B_pinv(result.B_pinv.template as<T>()), drift(result.drift.template as<T>()) {
        if (!result.success && std::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(result.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
    }

    /**
     * @brief The input that moves the model from @p x_ref to @p x_ref_next in one sample (least squares)
     * @param x_ref      State reference now, x*ₖ
     * @param x_ref_next State reference one sample ahead, x*ₖ₊₁
     */
    [[nodiscard]] constexpr ColVec<NU, T> control(const ColVec<NX, T>& x_ref, const ColVec<NX, T>& x_ref_next) const {
        return ColVec<NU, T>((B_pinv * ColVec<NX, T>(x_ref_next - x_ref)) - (drift * x_ref));
    }
};

/** @brief Deduce the feedforward from its design, keeping the result's scalar: StateFeedforward ff{res.as<float>()}. */
template<size_t NX, size_t NU, typename T>
StateFeedforward(const design::StateFeedforwardResult<NX, NU, T>&) -> StateFeedforward<NX, NU, T>;

/**
 * @brief The stable inverse at run time: the reference point (x*, u*) from the output reference
 *        over the preview, seven matrix–vector products and the zero dynamics' state
 *
 * Runs design::output_inverse. Each tick takes y*ₖ … y*ₖ₊NP₋₁ and returns x*ₖ and u*ₖ, the
 * state and the input now of the model's trajectory through them. A law takes them as its
 * reference point: `law.control({}, p.x, p.u)` for a law with a target (OFLQG), whose r is then
 * zero. The preview advances one sample per tick. The zero dynamics start at rest, and reset()
 * returns them there.
 *
 * @tparam NX States of the model; NU Inputs and controlled outputs; NP Preview samples
 */
template<size_t NX, size_t NU, size_t NP, typename T = float, typename Store = detail::HeldOutputInverse<NX, NU, NP, T>>
struct OutputInverse {
    static constexpr size_t ND = NU * (NP - 1);     ///< Preview steps, stacked
    using Preview = damp::array<ColVec<NU, T>, NP>; ///< y*ₖ … y*ₖ₊NP₋₁: preview[0] is now

    [[no_unique_address]] typename Store::x_hold        x_hold{};              ///< x* per unit y*ₖ
    [[no_unique_address]] typename Store::u_hold        u_hold{};              ///< u* per unit y*ₖ
    [[no_unique_address]] typename Store::x_preview     x_preview{};           ///< x* per unit preview step
    [[no_unique_address]] typename Store::u_preview     u_preview{};           ///< u* per unit preview step
    [[no_unique_address]] typename Store::zero_dynamics zero_dynamics{};       ///< The stable zero dynamics
    [[no_unique_address]] typename Store::zero_preview  zero_preview{};        ///< Their input per unit preview step
    [[no_unique_address]] typename Store::u_zero        u_zero{};              ///< u* per unit zero-dynamics state
    ColVec<NX, T>                                       zero_dynamics_state{}; ///< The stable zero dynamics' state (0 at rest)

    constexpr OutputInverse() = default;

    /** @brief From design::output_inverse (any scalar). A failed design does not compile when constant-evaluated. */
    template<typename U>
        requires(!Store::constant)
    constexpr OutputInverse(const design::OutputInverseResult<NX, NU, NP, U>& result) // NOLINT
        : x_hold(result.x_hold.template as<T>()),
          u_hold(result.u_hold.template as<T>()),
          x_preview(result.x_preview.template as<T>()),
          u_preview(result.u_preview.template as<T>()),
          zero_dynamics(result.zero_dynamics.template as<T>()),
          zero_preview(result.zero_preview.template as<T>()),
          u_zero(result.u_zero.template as<T>()) {
        if (!result.success && std::is_constant_evaluated()) {
            ::damp::detail::refuse_failed_design(result.status); // constant-evaluated (constinit / constexpr): a failed design does not compile
        }
    }

    /**
     * @brief One tick: the reference point now from the preview, then the zero dynamics advance
     * @param y_ref The controlled outputs' reference, y_ref[j] = y*ₖ₊ⱼ (y_ref[0] now)
     * @return x*ₖ and u*ₖ
     */
    [[nodiscard]] constexpr TargetPoint<NX, NU, T> step(const Preview& y_ref) {
        ColVec<ND, T> steps{};
        for (size_t l = 1; l < NP; ++l) {
            for (size_t i = 0; i < NU; ++i) {
                steps((NU * (l - 1)) + i) = y_ref[l](i) - y_ref[0](i);
            }
        }
        TargetPoint<NX, NU, T> p{};
        p.x = ColVec<NX, T>((x_hold * y_ref[0]) + (x_preview * steps) + zero_dynamics_state);
        p.u = ColVec<NU, T>((u_hold * y_ref[0]) + (u_preview * steps) + (u_zero * zero_dynamics_state));
        zero_dynamics_state = ColVec<NX, T>((zero_dynamics * zero_dynamics_state) + (zero_preview * steps));
        return p;
    }

    /** @brief Return the zero dynamics to rest. */
    constexpr void reset() { zero_dynamics_state = ColVec<NX, T>{}; }
};

/** @brief Deduce the inverse from its design, keeping the result's scalar: OutputInverse inverse{res.as<float>()}. */
template<size_t NX, size_t NU, size_t NP, typename T>
OutputInverse(const design::OutputInverseResult<NX, NU, NP, T>&) -> OutputInverse<NX, NU, NP, T>;

namespace detail {

/** @brief The runtime a model-inverse design runs as (none for any other result). */
template<typename R, typename T>
struct inverse_runtime {};

/** @brief design::state_feedforward runs as StateFeedforward. */
template<size_t NX, size_t NU, typename U, typename T>
struct inverse_runtime<design::StateFeedforwardResult<NX, NU, U>, T> {
    using type = StateFeedforward<NX, NU, T>; ///< StateFeedforward
    template<const auto& Res>
    using on_constants = StateFeedforward<NX, NU, T, ConstantFeedforward<Res, T>>; ///< StateFeedforward on the design's constants
};

/** @brief design::output_inverse runs as OutputInverse. */
template<size_t NX, size_t NU, size_t NP, typename U, typename T>
struct inverse_runtime<design::OutputInverseResult<NX, NU, NP, U>, T> {
    using type = OutputInverse<NX, NU, NP, T>; ///< OutputInverse
    template<const auto& Res>
    using on_constants = OutputInverse<NX, NU, NP, T, ConstantOutputInverse<Res, T>>; ///< OutputInverse on the design's constants
};

/** @brief A model-inverse design. */
template<const auto& Res, typename T>
concept inverse_design = requires { typename inverse_runtime<std::remove_cvref_t<decltype(Res)>, T>::type; };

} // namespace detail

/** @brief Law<kRes, T> of a `constexpr` model inverse: the same runtime on the design's matrices as constants of the type. */
template<const auto& Res, typename T>
    requires(detail::inverse_design<Res, T> && detail::constant_design<Res>)
struct Law<Res, T> : detail::inverse_runtime<std::remove_cvref_t<decltype(Res)>, T>::template on_constants<Res> {
    using Runtime = detail::inverse_runtime<std::remove_cvref_t<decltype(Res)>, T>::type; ///< The runtime with members: the type of a member a redesign assigns

    constexpr Law() {
        if (!Res.success) {
            ::damp::detail::refuse_failed_design(Res.status);
        }
    }
};

/** @brief Law<kRes, T> of a model inverse initialized at start-up: StateFeedforward or OutputInverse, default-constructed from the design. */
template<const auto& Res, typename T>
    requires detail::inverse_design<Res, T>
struct Law<Res, T> : detail::inverse_runtime<std::remove_cvref_t<decltype(Res)>, T>::type {
    using Runtime = detail::inverse_runtime<std::remove_cvref_t<decltype(Res)>, T>::type; ///< The runtime the design runs as

    /** @brief The design's runtime. */
    constexpr Law() : Runtime{Res.template as<T>()} {}

    using Runtime::operator=;
};

} // namespace damp
