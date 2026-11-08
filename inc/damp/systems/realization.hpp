// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file realization.hpp
 * @brief State-coordinate changes of a StateSpace: similarity, balancing, canonical and staircase forms
 *
 * Every function here keeps the input-output behaviour and changes only the state basis,
 * z = T·x (T is MATLAB®'s transform): a general similarity, diagonal balancing, the real
 * modal form, the SISO controllable / observable canonical and companion forms, and the
 * controllability / observability staircase forms of any model. TransferFunction realizes
 * the same canonical forms directly from coefficients (to_controllable_canonical_form,
 * to_observable_canonical_form).
 *
 * Companion-type forms (canonical, companion) carry the characteristic polynomial's
 * coefficients, so their sensitivity to rounding grows quickly with the order: use them for
 * small models and analysis; deploy a physical, modal or balanced realization.
 */


#include "damp/backend.hpp"
#include "damp/design/status.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "state_space.hpp"

namespace damp {

/**
 * @brief A model in balanced state coordinates, and the diagonal scaling that took it there
 *
 * The scaling is x = diag(2^exponent)·x̃: powers of two, so it is exact in floating point.
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
struct [[nodiscard]] StateScaling {
    StateSpace<NX, NU, NY, T, NW, NV> sys{};      ///< The model in the scaled state x̃
    damp::array<int, NX>              exponent{}; ///< x_i = 2^exponent_i · x̃_i (0: state unscaled)

    /** @brief The scale factors 2^exponent: x = diag(scale)·x̃, so MATLAB®'s T (x̃ = T x) is diag(1 / scale). */
    [[nodiscard]] constexpr ColVec<NX, T> scale() const {
        ColVec<NX, T> s{};
        for (size_t i = 0; i < NX; ++i) {
            s[i] = mat::detail::scale_by_pow2(T{1}, exponent[i]);
        }
        return s;
    }

    template<typename U>
    [[nodiscard]] constexpr StateScaling<NX, NU, NY, U, NW, NV> as() const {
        return StateScaling<NX, NU, NY, U, NW, NV>{sys.template as<U>(), exponent};
    }
};

/**
 * @brief Balance a model by a diagonal state scaling, so A, B and C have comparable magnitudes
 *
 * Applies the Parlett–Reinsch / LAPACK xGEBAL iteration (mat::detail::balance_exponents) to
 * @f$ \left[\begin{smallmatrix} A & B \\ C & 0 \end{smallmatrix}\right] @f$ with the inputs
 * and outputs held fixed: only the states are scaled, x = D·x̃ with D = diag(2^e), giving
 * @f$ (D^{-1} A D,\; D^{-1} B,\; C D,\; D) @f$ and @f$ D^{-1} G @f$ for the noise input.
 * Each state's row (its A row and B row) and column (its A column and C column) then have
 * 2-norms within a factor of two of each other. A state in counts beside one in metres stops
 * costing digits in eigenvalue, Riccati and frequency-response computations. The scaling is
 * a similarity by powers of two, so it is exact and the transfer function does not change.
 *
 * A model with a non-finite entry is returned unscaled.
 *
 * @note Compare with MATLAB®'s ssbal(sys), and with prescale(sys), which also scales the
 *       states but chooses the factors over a frequency band; matlab::prescale is this function.
 * @see "Balancing a Matrix for Calculation of Eigenvalues and Eigenvectors" (Parlett & Reinsch, 1969)
 * @see LAPACK xGEBAL (JOB = 'S')
 *
 * @param sys  Model, continuous or discrete
 * @return The balanced model and its scaling exponents
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr StateScaling<NX, NU, NY, T, NW, NV> balance_states(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    //! [|A| ‖B row‖; ‖C column‖ 0], the I/O folded into one fixed index (NX): its row and
    //! column norms are those of each state's full row and column.
    Matrix<NX + 1, NX + 1, T> mag{};
    for (size_t i = 0; i < NX; ++i) {
        for (size_t j = 0; j < NX; ++j) {
            mag(i, j) = (i == j) ? T{0} : damp::abs(sys.A(i, j));
        }
        mag(i, NX) = damp::detail::scaled_norm2<T>(NU, [&](size_t j) { return damp::abs(sys.B(i, j)); });
        mag(NX, i) = damp::detail::scaled_norm2<T>(NY, [&](size_t k) { return damp::abs(sys.C(k, i)); });
    }
    damp::array<int, NX + 1>            e{};
    StateScaling<NX, NU, NY, T, NW, NV> out{.sys = sys};
    if (!mat::detail::gebal_log2_scale(mag, e)) {
        return out; // a non-finite entry: leave the model unscaled
    }
    //! A common shift leaves D⁻¹AD unchanged; shifting the I/O index to 0 makes it a pure state scaling.
    for (size_t i = 0; i < NX; ++i) {
        out.exponent[i] = e[i] - e[NX];
    }
    for (size_t i = 0; i < NX; ++i) {
        const int ei = out.exponent[i];
        for (size_t j = 0; j < NX; ++j) {
            out.sys.A(i, j) = mat::detail::scale_by_pow2(sys.A(i, j), out.exponent[j] - ei);
        }
        for (size_t j = 0; j < NU; ++j) {
            out.sys.B(i, j) = mat::detail::scale_by_pow2(sys.B(i, j), -ei);
        }
        for (size_t j = 0; j < NW; ++j) {
            out.sys.G(i, j) = mat::detail::scale_by_pow2(sys.G(i, j), -ei);
        }
        for (size_t k = 0; k < NY; ++k) {
            out.sys.C(k, i) = mat::detail::scale_by_pow2(sys.C(k, i), ei);
        }
    }
    return out;
}

/**
 * @brief A model in new state coordinates z = transform·x, and the transform that took it there
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
struct [[nodiscard]] StateTransformResult {
    StateSpace<NX, NU, NY, T, NW, NV> sys{};                                       ///< The model in the new coordinates
    Matrix<NX, NX, T>                 transform{};                                 ///< z = transform·x (MATLAB®'s T)
    bool                              success{false};                              ///< status == DesignStatus::ok
    design::DesignStatus              status{design::DesignStatus::not_converged}; ///< ok, or why there is no such form

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    template<typename U>
    [[nodiscard]] constexpr StateTransformResult<NX, NU, NY, U, NW, NV> as() const {
        return StateTransformResult<NX, NU, NY, U, NW, NV>{sys.template as<U>(), transform.template as<U>(), success, status};
    }
};

namespace detail {

/**
 * @brief Solve the small Sylvester equation S_ii·X − X·S_jj = R for a p × q block (p, q ∈ {1, 2}),
 * through its Kronecker form padded to 4 × 4. nullopt when the blocks share an eigenvalue.
 */
template<size_t N, typename T>
[[nodiscard]] constexpr damp::optional<Matrix<2, 2, T>> modal_sylvester(
    const Matrix<N, N, T>& S,
    size_t                 i0,
    size_t                 p,
    size_t                 j0,
    size_t                 q,
    const Matrix<2, 2, T>& R
) {
    //! vec(X) column-major: index r + p·c; K = I_q ⊗ S_ii − S_jjᵀ ⊗ I_p.
    Matrix<4, 4, T> K = Matrix<4, 4, T>::identity();
    Matrix<4, 1, T> rhs{};
    for (size_t c = 0; c < q; ++c) {
        for (size_t r = 0; r < p; ++r) {
            const size_t row = r + (p * c);
            rhs(row, 0) = R(r, c);
            for (size_t c2 = 0; c2 < q; ++c2) {
                for (size_t r2 = 0; r2 < p; ++r2) {
                    T k = T{0};
                    if (c == c2) {
                        k += S(i0 + r, i0 + r2);
                    }
                    if (r == r2) {
                        k -= S(j0 + c2, j0 + c);
                    }
                    K(row, r2 + (p * c2)) = k;
                }
            }
        }
    }
    const auto x = mat::lu_solve(K, rhs);
    if (!x) {
        return damp::nullopt;
    }
    Matrix<2, 2, T> X{};
    for (size_t c = 0; c < q; ++c) {
        for (size_t r = 0; r < p; ++r) {
            X(r, c) = (*x)(r + (p * c), 0);
        }
    }
    return X;
}

/** @brief Diagonal 1×1 / 2×2 blocks of a real Schur form: a nonzero subdiagonal opens a pair. */
template<size_t N>
struct SchurBlocks {
    damp::array<size_t, N> start{};  ///< First row of each block
    damp::array<size_t, N> size{};   ///< 1 (real mode) or 2 (complex pair)
    size_t                 count{0}; ///< Number of blocks
};

template<size_t N, typename T>
[[nodiscard]] constexpr SchurBlocks<N> schur_blocks(const Matrix<N, N, T>& S) {
    SchurBlocks<N> b{};
    for (size_t i = 0; i < N;) {
        const size_t s = (i + 1 < N && S(i + 1, i) != T{0}) ? 2 : 1;
        b.start[b.count] = i;
        b.size[b.count] = s;
        ++b.count;
        i += s;
    }
    return b;
}

/**
 * @brief Block-diagonalize a real Schur form S: V unit block upper triangular with S·V = V·blkdiag(S).
 * false when two blocks share an eigenvalue and the coupling between them is not zero (a
 * defective, Jordan-type, mode that has no modal form).
 */
template<size_t N, typename T>
[[nodiscard]] constexpr bool block_diagonalize(const Matrix<N, N, T>& S, const SchurBlocks<N>& b, Matrix<N, N, T>& V) {
    V = Matrix<N, N, T>::identity();
    const T tol = static_cast<T>(16 * N) * eps<T>() * S.norm();
    for (size_t bj = 1; bj < b.count; ++bj) {
        const size_t j0 = b.start[bj];
        const size_t q = b.size[bj];
        for (size_t bi = bj; bi-- > 0;) {
            const size_t i0 = b.start[bi];
            const size_t p = b.size[bi];
            //! R = −Σ_{l = i+1 … j} S_il·V_lj (V_jj = I)
            Matrix<2, 2, T> R{};
            T               rmax = T{0};
            for (size_t r = 0; r < p; ++r) {
                for (size_t c = 0; c < q; ++c) {
                    T acc = T{0};
                    for (size_t m = i0 + p; m < j0 + q; ++m) {
                        acc += S(i0 + r, m) * V(m, j0 + c);
                    }
                    R(r, c) = -acc;
                    rmax = damp::max(rmax, damp::abs(acc));
                }
            }
            auto X = modal_sylvester(S, i0, p, j0, q, R);
            if (!X) {
                if (rmax > tol) {
                    return false; // shared eigenvalue, coupled: defective
                }
                X = Matrix<2, 2, T>{}; // shared eigenvalue, uncoupled: already split
            }
            for (size_t r = 0; r < p; ++r) {
                for (size_t c = 0; c < q; ++c) {
                    V(i0 + r, j0 + c) = (*X)(r, c);
                }
            }
        }
    }
    return true;
}

} // namespace detail

/**
 * @brief Real modal canonical form: A block diagonal, one block per mode
 *
 * A real eigenvalue λ is a 1 × 1 block; a complex pair σ ± jω is the 2 × 2 block
 * @f$ \left[\begin{smallmatrix} \sigma & \omega \\ -\omega & \sigma \end{smallmatrix}\right] @f$,
 * ω > 0. The modes come from the real Schur form of the balanced A (Francis QR, the same
 * path as mat::compute_eigenvalues), which a unit block upper-triangular similarity then
 * splits into its diagonal blocks by solving one small Sylvester equation per block pair
 * (Bavely–Stewart), and a diagonal scaling brings each 2 × 2 block to the form above. With
 * x = X·z, the result is @f$ (X^{-1} A X,\; X^{-1} B,\; C X,\; D) @f$ and
 * @f$ X^{-1} G @f$, the inverses applied by LU solves; @c transform is @f$ X^{-1} @f$
 * (z = transform·x), MATLAB®'s T.
 *
 * The block structure is exact by construction, and the transform is certified by the
 * residual of A·X = X·A_modal. Close eigenvalues make X ill-conditioned, as in any modal
 * decomposition.
 *
 * Statuses: DesignStatus::not_converged when the QR iteration does not deflate or the
 * certificate fails; DesignStatus::singular for a defective A (a repeated eigenvalue with
 * fewer eigenvectors than its multiplicity, such as a double integrator), which has no
 * modal form, or a transform that cannot be solved against.
 *
 * @note Compare with MATLAB®'s [csys, T] = canon(sys, 'modal'); matlab::canon is the alias.
 * @see "An Algorithm for Computing Reducing Subspaces by Block Diagonalization" (Bavely & Stewart, 1979), SIAM J. Numer. Anal. 16(2)
 * @see "Matrix Computations" (Golub & Van Loan, 2013), §7.6 (invariant subspace computations)
 *
 * @param sys  Model, continuous or discrete
 * @return The modal model, the transform z = transform·x, and status
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr StateTransformResult<NX, NU, NY, T, NW, NV> modal_form(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    using Status = design::DesignStatus;
    StateTransformResult<NX, NU, NY, T, NW, NV> out{};

    //! Real Schur form of the balanced A: D⁻¹AD = Z·S·Zᵀ, D = diag(2^e).
    const auto        e = mat::detail::balance_exponents(sys.A);
    Matrix<NX, NX, T> S = mat::detail::balanced(sys.A, e);
    Matrix<NX, NX, T> Z = Matrix<NX, NX, T>::identity();
    if constexpr (NX > 1) {
        mat::detail::hessenberg_reduce(S, Z);
        damp::array<T, NX> wr{};
        damp::array<T, NX> wi{};
        if (!mat::detail::francis_qr(S, Z, wr, wi)) {
            return detail::set_status(out, Status::not_converged);
        }
    }
    const auto        blocks = detail::schur_blocks(S);
    Matrix<NX, NX, T> V{};
    if (!detail::block_diagonalize(S, blocks, V)) {
        return detail::set_status(out, Status::singular); // defective: no modal form
    }

    //! X = D·Z·V·W, W scaling each pair [m b; c m] to [σ ω; −ω σ] (column k+1 times ω/b).
    Matrix<NX, NX, T> DZ{};
    for (size_t i = 0; i < NX; ++i) {
        for (size_t j = 0; j < NX; ++j) {
            DZ(i, j) = mat::detail::scale_by_pow2(Z(i, j), e[i]);
        }
    }
    Matrix<NX, NX, T> X = DZ * V;
    Matrix<NX, NX, T> Am{};
    for (size_t k = 0; k < blocks.count; ++k) {
        const size_t i = blocks.start[k];
        if (blocks.size[k] == 1) {
            Am(i, i) = S(i, i);
            continue;
        }
        const T b = S(i, i + 1);
        const T w = damp::sqrt(damp::abs(b)) * damp::sqrt(damp::abs(S(i + 1, i)));
        const T sigma = static_cast<T>(0.5) * (S(i, i) + S(i + 1, i + 1));
        for (size_t r = 0; r < NX; ++r) {
            X(r, i + 1) = X(r, i + 1) * (w / b);
        }
        Am(i, i) = sigma;
        Am(i, i + 1) = w;
        Am(i + 1, i) = -w;
        Am(i + 1, i + 1) = sigma;
    }

    //! Certify A·X = X·A_modal against the data's scale.
    const Matrix<NX, NX, T> residual = (sys.A * X) - (X * Am);
    const T                 scale = sys.A.norm() * X.norm();
    if (!(residual.norm() <= static_cast<T>(256 * NX) * eps<T>() * scale)) {
        return detail::set_status(out, Status::not_converged);
    }
    const auto Xinv = mat::lu_solve(X, Matrix<NX, NX, T>::identity()); // the transform is a deliverable
    const auto Bm = mat::lu_solve(X, sys.B);
    const auto Gm = mat::lu_solve(X, sys.G);
    if (!Xinv || !Bm || !Gm) {
        return detail::set_status(out, Status::singular);
    }
    out.transform = *Xinv;
    out.sys = StateSpace<NX, NU, NY, T, NW, NV>{
        .A = Am,
        .B = *Bm,
        .C = sys.C * X,
        .D = sys.D,
        .G = *Gm,
        .H = sys.H,
        .Ts = sys.Ts,
    };
    return detail::set_status(out, Status::ok);
}

/**
 * @brief Change the state basis to z = T·x: (T·A·T⁻¹, T·B, C·T⁻¹, D), and T·G for the noise input
 *
 * T⁻¹ is never formed: one LU factorization of Tᵀ gives A·T⁻¹ and C·T⁻¹ as solves of
 * @f$ T^\top Y^\top = (\cdot)^\top @f$. The transfer function does not change. A singular T
 * (one LU cannot pivot through) is refused with DesignStatus::singular.
 *
 * @note Compare with MATLAB®'s ss2ss(sys, T); matlab::ss2ss is the alias.
 *
 * @param sys        Model, continuous or discrete
 * @param transform  The transform T, z = T·x
 * @return The model in z, the transform, and status
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr StateTransformResult<NX, NU, NY, T, NW, NV>
similarity_transform(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const Matrix<NX, NX, T>& transform) {
    StateTransformResult<NX, NU, NY, T, NW, NV> out{};
    const Matrix<NX, NX, T>                     Tt = transform.transpose();
    const auto                                  ATinvT = mat::lu_solve(Tt, Matrix<NX, NX, T>((transform * sys.A).transpose())); // (T·A·T⁻¹)ᵀ
    const auto                                  CTinvT = mat::lu_solve(Tt, Matrix<NX, NY, T>(sys.C.transpose()));               // (C·T⁻¹)ᵀ
    if (!ATinvT || !CTinvT) {
        return detail::set_status(out, design::DesignStatus::singular);
    }
    out.transform = transform;
    out.sys = StateSpace<NX, NU, NY, T, NW, NV>{
        .A = ATinvT->transpose(),
        .B = transform * sys.B,
        .C = CTinvT->transpose(),
        .D = sys.D,
        .G = transform * sys.G,
        .H = sys.H,
        .Ts = sys.Ts,
    };
    return detail::set_status(out, design::DesignStatus::ok);
}

namespace detail {

/** @brief [b, A·b, …, Aⁿ⁻¹·b] for one column b, and Aⁿ·b in @p next. */
template<size_t NX, typename T>
[[nodiscard]] constexpr Matrix<NX, NX, T> krylov(const Matrix<NX, NX, T>& A, const ColVec<NX, T>& b, ColVec<NX, T>& next) {
    Matrix<NX, NX, T> K{};
    ColVec<NX, T>     v = b;
    for (size_t j = 0; j < NX; ++j) {
        for (size_t i = 0; i < NX; ++i) {
            K(i, j) = v[i];
        }
        v = ColVec<NX, T>(A * v);
    }
    next = v;
    return K;
}

/** @brief The controllable canonical A of sⁿ + Σ αⱼ·sʲ: ones on the superdiagonal, −α in the last row. */
template<size_t NX, typename T>
[[nodiscard]] constexpr Matrix<NX, NX, T> phase_variable_A(const Matrix<NX, 1, T>& alpha) {
    Matrix<NX, NX, T> A{};
    for (size_t i = 0; i + 1 < NX; ++i) {
        A(i, i + 1) = T{1};
    }
    for (size_t j = 0; j < NX; ++j) {
        A(NX - 1, j) = -alpha(j, 0);
    }
    return A;
}

/** @brief eₖ as an NX × 1 column. */
template<size_t NX, typename T>
[[nodiscard]] constexpr Matrix<NX, 1, T> unit_column(size_t k) {
    Matrix<NX, 1, T> e{};
    e(k, 0) = T{1};
    return e;
}

/**
 * @brief Whether two SISO models share their first NX inverse moments C·A⁻ᵏ·B, to √ε of their size
 *
 * The forward Markov parameters of a stiff model are dominated by its fastest mode; the inverse
 * moments by its slowest, which sets the low-frequency and DC gain. Skipped (true) when @p a's A
 * is singular to working precision, ‖A‖·‖A⁻¹‖ ≥ 1/(64·NX·ε): its inverse moments are then
 * roundoff of a mode at the origin.
 */
template<size_t NX, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr bool same_inverse_moments(const StateSpace<NX, 1, 1, T, NW, NV>& a, const StateSpace<NX, 1, 1, T, NW, NV>& b) {
    const auto a_inv = mat::lu_solve(a.A, Matrix<NX, NX, T>::identity()); // for its norm: a condition estimate
    if (!a_inv) {
        return true;
    }
    const T condition = a.A.norm() * a_inv->norm();
    if (!damp::isfinite(condition) || !(condition * static_cast<T>(NX) * T{64} * eps<T>() < T{1})) {
        return true;
    }
    Matrix<NX, 1, T> va = a.B;
    Matrix<NX, 1, T> vb = b.B;
    T                diff{0};
    T                size{0};
    for (size_t k = 0; k < NX; ++k) {
        const auto na = mat::lu_solve(a.A, va);
        const auto nb = mat::lu_solve(b.A, vb);
        if (!na || !nb) {
            return false; // b's A is singular where a's is not
        }
        va = *na;
        vb = *nb;
        const T ha = (a.C * va)(0, 0);
        const T hb = (b.C * vb)(0, 0);
        diff = damp::max(diff, damp::abs(ha - hb));
        size = damp::max(size, damp::abs(ha));
    }
    return damp::isfinite(diff) && diff <= damp::sqrt(eps<T>()) * size;
}

/**
 * @brief Whether two SISO models share their first 2·NX Markov parameters C·Aᵏ·B and their first
 *        NX inverse moments C·A⁻ᵏ·B, each sequence to √ε of its size
 */
template<size_t NX, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr bool same_markov(const StateSpace<NX, 1, 1, T, NW, NV>& a, const StateSpace<NX, 1, 1, T, NW, NV>& b) {
    ColVec<NX, T> va = a.B.col(0);
    ColVec<NX, T> vb = b.B.col(0);
    T             diff{0};
    T             size{0};
    for (size_t k = 0; k < 2 * NX; ++k) {
        const T ha = (a.C * va)(0, 0);
        const T hb = (b.C * vb)(0, 0);
        diff = damp::max(diff, damp::abs(ha - hb));
        size = damp::max(size, damp::abs(ha));
        va = ColVec<NX, T>(a.A * va);
        vb = ColVec<NX, T>(b.A * vb);
    }
    return damp::isfinite(diff) && diff <= damp::sqrt(eps<T>()) * size && same_inverse_moments(a, b);
}

} // namespace detail

/**
 * @brief Controllable canonical form of a SISO model: ones on A's superdiagonal, −α in its last row, B = eₙ
 *
 * The same matrices as TransferFunction::to_controllable_canonical_form of the model's
 * transfer function, reached from the state-space model by z = T·x with
 * @f$ T = W_{cc} W_c^{-1} @f$, where @f$ W_c = [b, Ab, \ldots, A^{n-1}b] @f$ is the
 * controllability matrix and @f$ W_{cc} @f$ that of the canonical pair. The characteristic
 * coefficients come from the same matrix (Cayley–Hamilton: @f$ W_c \alpha = -A^n b @f$), so A
 * and B are exact by construction; C and G are transformed by solves. The result is certified
 * by its first 2n Markov parameters C·Aᵏ·B, which the fastest mode dominates, and, when A is
 * invertible to working precision, its first n inverse moments C·A⁻ᵏ·B, which the slowest mode
 * dominates, each to √ε of its size: a companion form loses digits with the order and the
 * spread of the modes, and one that has lost too many is refused.
 *
 * Statuses: DesignStatus::singular when (A, b) is not controllable (W_c singular);
 * DesignStatus::not_converged when the Markov certificate fails.
 *
 * @note Compare with MATLAB®'s ss(tf(sys)) for a SISO model.
 * @see "Linear System Theory and Design" (Chen, 1999), §4.3.1
 *
 * @param sys  SISO model, continuous or discrete
 * @return The canonical model, the transform z = transform·x, and status
 */
template<size_t NX, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr StateTransformResult<NX, 1, 1, T, NW, NV> controllable_canonical_form(const StateSpace<NX, 1, 1, T, NW, NV>& sys) {
    StateTransformResult<NX, 1, 1, T, NW, NV> out{};
    ColVec<NX, T>                             Anb{};
    const Matrix<NX, NX, T>                   Wc = detail::krylov(sys.A, ColVec<NX, T>(sys.B.col(0)), Anb);
    const auto                                alpha = mat::lu_solve(Wc, Matrix<NX, 1, T>(Anb * T{-1}));
    if (!alpha) {
        return detail::set_status(out, design::DesignStatus::singular); // (A, b) is not controllable
    }
    const Matrix<NX, NX, T> Ac = detail::phase_variable_A(*alpha);
    const Matrix<NX, 1, T>  en = detail::unit_column<NX, T>(NX - 1);
    ColVec<NX, T>           unused{};
    const Matrix<NX, NX, T> Wcc = detail::krylov(Ac, ColVec<NX, T>(en.col(0)), unused);
    //! T = Wcc·Wc⁻¹ and C·T⁻¹ = C·Wc·Wcc⁻¹, both as transposed solves.
    const auto Tt = mat::lu_solve(Matrix<NX, NX, T>(Wc.transpose()), Matrix<NX, NX, T>(Wcc.transpose()));
    const auto CcT = mat::lu_solve(Matrix<NX, NX, T>(Wcc.transpose()), Matrix<NX, 1, T>((sys.C * Wc).transpose()));
    if (!Tt || !CcT) {
        return detail::set_status(out, design::DesignStatus::singular);
    }
    out.transform = Tt->transpose();
    out.sys = StateSpace<NX, 1, 1, T, NW, NV>{
        .A = Ac,
        .B = en,
        .C = CcT->transpose(),
        .D = sys.D,
        .G = out.transform * sys.G,
        .H = sys.H,
        .Ts = sys.Ts,
    };
    const bool certified = detail::same_markov(sys, out.sys);
    return detail::set_status(out, certified ? design::DesignStatus::ok : design::DesignStatus::not_converged);
}

/**
 * @brief Observable canonical form of a SISO model: ones on A's subdiagonal, −α in its last column, C = eₙᵀ
 *
 * The dual of controllable_canonical_form and the same matrices as
 * TransferFunction::to_observable_canonical_form: z = T·x with @f$ T = W_{oo}^{-1} W_o @f$,
 * @f$ W_o = [c; cA; \ldots; cA^{n-1}] @f$ the observability matrix, and the coefficients from
 * @f$ W_o^\top \alpha = -(cA^n)^\top @f$. A and C are exact by construction, B and G are
 * transformed, and the result is certified by its Markov parameters as the controllable form.
 *
 * Statuses: DesignStatus::singular when (A, c) is not observable; DesignStatus::not_converged
 * when the Markov certificate fails.
 *
 * @see "Linear System Theory and Design" (Chen, 1999), §5.3
 *
 * @param sys  SISO model, continuous or discrete
 * @return The canonical model, the transform z = transform·x, and status
 */
template<size_t NX, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr StateTransformResult<NX, 1, 1, T, NW, NV> observable_canonical_form(const StateSpace<NX, 1, 1, T, NW, NV>& sys) {
    StateTransformResult<NX, 1, 1, T, NW, NV> out{};
    ColVec<NX, T>                             cAn{};
    const Matrix<NX, NX, T>                   WoT = detail::krylov(Matrix<NX, NX, T>(sys.A.transpose()), ColVec<NX, T>(sys.C.transpose().col(0)), cAn);
    const auto                                alpha = mat::lu_solve(WoT, Matrix<NX, 1, T>(cAn * T{-1}));
    if (!alpha) {
        return detail::set_status(out, design::DesignStatus::singular); // (A, c) is not observable
    }
    const Matrix<NX, NX, T> Ao = detail::phase_variable_A(*alpha).transpose();
    const Matrix<NX, 1, T>  en = detail::unit_column<NX, T>(NX - 1);
    ColVec<NX, T>           unused{};
    const Matrix<NX, NX, T> WooT = detail::krylov(Matrix<NX, NX, T>(Ao.transpose()), ColVec<NX, T>(en.col(0)), unused);
    const auto              Tm = mat::lu_solve(Matrix<NX, NX, T>(WooT.transpose()), Matrix<NX, NX, T>(WoT.transpose())); // T = Woo⁻¹·Wo
    if (!Tm) {
        return detail::set_status(out, design::DesignStatus::singular);
    }
    out.transform = *Tm;
    out.sys = StateSpace<NX, 1, 1, T, NW, NV>{
        .A = Ao,
        .B = out.transform * sys.B,
        .C = en.transpose(),
        .D = sys.D,
        .G = out.transform * sys.G,
        .H = sys.H,
        .Ts = sys.Ts,
    };
    const bool certified = detail::same_markov(sys, out.sys);
    return detail::set_status(out, certified ? design::DesignStatus::ok : design::DesignStatus::not_converged);
}

/**
 * @brief Companion form of a SISO model, MATLAB®'s canon(sys, 'companion'): B = e₁, −α in A's last column
 *
 * z = W_c⁻¹·x, with W_c the controllability matrix: A becomes @f$ W_c^{-1} A W_c @f$, whose
 * first n − 1 columns are e₂ … eₙ and whose last is −α (the characteristic polynomial
 * sⁿ + Σ αⱼ·sʲ), B becomes e₁, and C becomes C·W_c, the Markov parameters. A and B are exact
 * by construction; the result is certified by its Markov parameters and inverse moments, as
 * controllable_canonical_form.
 *
 * Statuses: DesignStatus::singular when (A, b) is not controllable; DesignStatus::not_converged
 * when the Markov certificate fails.
 *
 * @note Compare with MATLAB®'s [csys, T] = canon(sys, 'companion'); matlab::canon is the alias.
 *
 * @param sys  SISO model, continuous or discrete
 * @return The companion model, the transform z = transform·x, and status
 */
template<size_t NX, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr StateTransformResult<NX, 1, 1, T, NW, NV> companion_form(const StateSpace<NX, 1, 1, T, NW, NV>& sys) {
    StateTransformResult<NX, 1, 1, T, NW, NV> out{};
    ColVec<NX, T>                             Anb{};
    const Matrix<NX, NX, T>                   Wc = detail::krylov(sys.A, ColVec<NX, T>(sys.B.col(0)), Anb);
    const auto                                alpha = mat::lu_solve(Wc, Matrix<NX, 1, T>(Anb * T{-1}));
    const auto                                Winv = mat::lu_solve(Wc, Matrix<NX, NX, T>::identity()); // the transform is a deliverable
    const auto                                Gz = mat::lu_solve(Wc, sys.G);
    if (!alpha || !Winv || !Gz) {
        return detail::set_status(out, design::DesignStatus::singular); // (A, b) is not controllable
    }
    const Matrix<NX, NX, T> Az = detail::phase_variable_A(*alpha).transpose();
    out.transform = *Winv;
    out.sys = StateSpace<NX, 1, 1, T, NW, NV>{
        .A = Az,
        .B = detail::unit_column<NX, T>(0),
        .C = sys.C * Wc,
        .D = sys.D,
        .G = *Gz,
        .H = sys.H,
        .Ts = sys.Ts,
    };
    const bool certified = detail::same_markov(sys, out.sys);
    return detail::set_status(out, certified ? design::DesignStatus::ok : design::DesignStatus::not_converged);
}

/**
 * @brief A model in staircase coordinates, the orthogonal transform that took it there, and its rank
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
struct [[nodiscard]] StaircaseFormResult {
    StateSpace<NX, NU, NY, T, NW, NV> sys{};                                       ///< The model in staircase coordinates
    Matrix<NX, NX, T>                 transform{};                                 ///< z = transform·x, orthogonal
    size_t                            rank{0};                                     ///< Controllable (observable) states: the last rank states of z
    bool                              success{false};                              ///< status == DesignStatus::ok
    design::DesignStatus              status{design::DesignStatus::not_converged}; ///< ok, or why the staircase failed

    /** @brief True when the design succeeded (success). */
    [[nodiscard]] constexpr explicit operator bool() const { return success; }

    template<typename U>
    [[nodiscard]] constexpr StaircaseFormResult<NX, NU, NY, U, NW, NV> as() const {
        return StaircaseFormResult<NX, NU, NY, U, NW, NV>{sys.template as<U>(), transform.template as<U>(), rank, success, status};
    }
};

namespace detail {

/** @brief M ← H·M on rows [s, R), H = I − β·v·vᵀ. */
template<size_t R, size_t C, typename T>
constexpr void reflect_rows(Matrix<R, C, T>& M, const ColVec<R, T>& v, size_t s, T beta) {
    for (size_t j = 0; j < C; ++j) {
        T w{0};
        for (size_t i = s; i < R; ++i) {
            w += v[i] * M(i, j);
        }
        for (size_t i = s; i < R; ++i) {
            M(i, j) -= beta * v[i] * w;
        }
    }
}

/** @brief M ← M·H on columns [s, C), H = I − β·v·vᵀ. */
template<size_t R, size_t C, typename T>
constexpr void reflect_cols(Matrix<R, C, T>& M, const ColVec<C, T>& v, size_t s, T beta) {
    for (size_t i = 0; i < R; ++i) {
        T w{0};
        for (size_t j = s; j < C; ++j) {
            w += M(i, j) * v[j];
        }
        for (size_t j = s; j < C; ++j) {
            M(i, j) -= beta * w * v[j];
        }
    }
}

/**
 * @brief Controllability staircase in place, controllable states first; returns their count
 *
 * Each step takes the current "input" (B, then the coupling block of A just split off),
 * reduces it over the unreached rows [p, n) by column-pivoted Householder reflections, applied
 * to the whole model as similarities (A ← H·A·H, B ← H·B, G ← H·G, C ← C·H, Q ← Q·H), and
 * stops at the first pivot column whose norm is at most the tolerance of the matrix it reads
 * (@p tol_b for B, @p tol_a for a block of A): that is the step's rank.
 * The remainder below the new block is then at most @p tol and is set to zero, so the
 * staircase structure is exact. The reach stops when a step has rank 0.
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
constexpr size_t controllable_staircase(StateSpace<NX, NU, NY, T, NW, NV>& m, Matrix<NX, NX, T>& Q, T tol_b, T tol_a) {
    constexpr size_t NS = (NU > NX) ? NU : NX;
    size_t           p0 = 0;   // first unreached row
    size_t           prev = 0; // first column of the source block, when it lies in A
    size_t           width = NU;
    bool             from_b = true;
    auto             src = [&](size_t i, size_t c) -> T& { return from_b ? m.B(i, c) : m.A(i, prev + c); };
    while (p0 < NX && width > 0) {
        const T               tol = from_b ? tol_b : tol_a;
        damp::array<bool, NS> used{};
        size_t                r = 0;
        while (p0 + r < NX && r < width) {
            const size_t s = p0 + r;
            size_t       best = width;
            T            best_norm{0};
            for (size_t c = 0; c < width; ++c) {
                const T nc = used[c] ? T{0} : damp::detail::scaled_norm2<T>(NX - s, [&](size_t k) { return damp::abs(src(s + k, c)); });
                if (nc > best_norm) {
                    best_norm = nc;
                    best = c;
                }
            }
            if (best == width || !(best_norm > tol)) {
                break;
            }
            used[best] = true;
            const T       x0 = src(s, best);
            const T       alpha = (x0 < T{0}) ? best_norm : -best_norm;
            ColVec<NX, T> v{};
            v[s] = x0 - alpha;
            for (size_t k = s + 1; k < NX; ++k) {
                v[k] = src(k, best);
            }
            const T beta = T{1} / (best_norm * (best_norm + damp::abs(x0))); // 2 / vᵀv
            reflect_rows(m.A, v, s, beta);
            reflect_cols(m.A, v, s, beta);
            reflect_rows(m.B, v, s, beta);
            reflect_rows(m.G, v, s, beta);
            reflect_cols(m.C, v, s, beta);
            reflect_cols(Q, v, s, beta);
            for (size_t k = s + 1; k < NX; ++k) {
                src(k, best) = T{0};
            }
            ++r;
        }
        for (size_t c = 0; c < width; ++c) {
            for (size_t k = p0 + r; k < NX; ++k) {
                src(k, c) = T{0}; // below the rank tolerance: the staircase is exact
            }
        }
        if (r == 0) {
            break;
        }
        prev = p0;
        p0 += r;
        width = r;
        from_b = false;
    }
    return p0;
}

} // namespace detail

/**
 * @brief Controllability staircase form: the controllable part split off by an orthogonal transform
 *
 * z = T·x with T orthogonal, in MATLAB®'s layout
 * @f$ \bar A = \left[\begin{smallmatrix} A_{nc} & 0 \\ A_{21} & A_c \end{smallmatrix}\right] @f$,
 * @f$ \bar B = \left[\begin{smallmatrix} 0 \\ B_c \end{smallmatrix}\right] @f$,
 * @f$ \bar C = [C_{nc}\; C_c] @f$, with the last @c rank states controllable and the
 * uncontrollable modes in @f$ A_{nc} @f$; A_c is itself in staircase (block Hessenberg)
 * form. Any number of inputs and outputs. Built by column-pivoted Householder reflections, so
 * it is backward stable: the input step's rank is decided against n·ε·‖B‖_F and each coupling
 * step's against n·ε·‖A‖_F, each relative to the matrix it reads, so the rank does not depend
 * on the units of the input; the entries below that are set to zero. The transform is
 * certified by A and B reconstructing to those tolerances.
 *
 * Statuses: DesignStatus::bad_weights for a model with a non-finite entry;
 * DesignStatus::not_converged when the certificate fails.
 *
 * @note Compare with MATLAB®'s [Abar, Bbar, Cbar, T, k] = ctrbf(A, B, C), rank = sum(k);
 *       matlab::ctrbf is the alias.
 * @see "The Generalized Eigenstructure Problem in Linear System Theory" (Van Dooren, 1981), IEEE TAC 26(1)
 *
 * @param sys  Model, continuous or discrete
 * @return The staircase model, the orthogonal transform z = transform·x, the controllable rank, and status
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr StaircaseFormResult<NX, NU, NY, T, NW, NV> controllability_staircase(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    StaircaseFormResult<NX, NU, NY, T, NW, NV> out{};
    const T                                    norm_a = sys.A.norm();
    const T                                    norm_b = sys.B.norm();
    if (!damp::isfinite(norm_a) || !damp::isfinite(norm_b)) {
        return detail::set_status(out, design::DesignStatus::bad_weights);
    }
    const T                           tol_a = static_cast<T>(NX) * eps<T>() * norm_a;
    const T                           tol_b = static_cast<T>(NX) * eps<T>() * norm_b;
    StateSpace<NX, NU, NY, T, NW, NV> m = sys;
    Matrix<NX, NX, T>                 Q = Matrix<NX, NX, T>::identity();
    out.rank = detail::controllable_staircase(m, Q, tol_b, tol_a);

    //! MATLAB®'s layout: reverse the state order, uncontrollable states first, z = P·Qᵀ·x.
    out.sys = StateSpace<NX, NU, NY, T, NW, NV>{.D = sys.D, .H = sys.H, .Ts = sys.Ts};
    for (size_t i = 0; i < NX; ++i) {
        const size_t ri = NX - 1 - i;
        for (size_t j = 0; j < NX; ++j) {
            out.sys.A(i, j) = m.A(ri, NX - 1 - j);
            out.transform(i, j) = Q(j, ri);
        }
        for (size_t j = 0; j < NU; ++j) {
            out.sys.B(i, j) = m.B(ri, j);
        }
        for (size_t j = 0; j < NW; ++j) {
            out.sys.G(i, j) = m.G(ri, j);
        }
        for (size_t k = 0; k < NY; ++k) {
            out.sys.C(k, i) = m.C(k, ri);
        }
    }
    //! Certify: Tᵀ·Ā·T and Tᵀ·B̄ give A and B back, to the tolerance the rank decisions used.
    const Matrix<NX, NX, T> Tt = out.transform.transpose();
    const T                 resA = Matrix<NX, NX, T>((Tt * out.sys.A * out.transform) - sys.A).norm();
    const T                 resB = Matrix<NX, NU, T>((Tt * out.sys.B) - sys.B).norm();
    const T                 bound_a = static_cast<T>(64 * NX) * (tol_a + (eps<T>() * norm_a));
    const T                 bound_b = static_cast<T>(64 * NX) * (tol_b + (eps<T>() * norm_b));
    const bool              certified = (resA <= bound_a) && (resB <= bound_b);
    return detail::set_status(out, certified ? design::DesignStatus::ok : design::DesignStatus::not_converged);
}

/**
 * @brief Observability staircase form: the observable part split off by an orthogonal transform
 *
 * The dual of controllability_staircase, in MATLAB®'s layout
 * @f$ \bar A = \left[\begin{smallmatrix} A_{no} & A_{12} \\ 0 & A_o \end{smallmatrix}\right] @f$,
 * @f$ \bar B = \left[\begin{smallmatrix} B_{no} \\ B_o \end{smallmatrix}\right] @f$,
 * @f$ \bar C = [0\; C_o] @f$, with the last @c rank states observable: the controllability
 * staircase of (Aᵀ, Cᵀ), transposed back.
 *
 * Statuses: as controllability_staircase.
 *
 * @note Compare with MATLAB®'s [Abar, Bbar, Cbar, T, k] = obsvf(A, B, C), rank = sum(k);
 *       matlab::obsvf is the alias.
 *
 * @param sys  Model, continuous or discrete
 * @return The staircase model, the orthogonal transform z = transform·x, the observable rank, and status
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr StaircaseFormResult<NX, NU, NY, T, NW, NV> observability_staircase(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    const StateSpace<NX, NY, NU, T, NW, NV> dual{
        .A = sys.A.transpose(),
        .B = sys.C.transpose(),
        .C = sys.B.transpose(),
        .D = sys.D.transpose(),
        .G = sys.G,
        .Ts = sys.Ts,
    };
    const auto                                 d = controllability_staircase(dual);
    StaircaseFormResult<NX, NU, NY, T, NW, NV> out{};
    out.sys = StateSpace<NX, NU, NY, T, NW, NV>{
        .A = d.sys.A.transpose(),
        .B = d.sys.C.transpose(),
        .C = d.sys.B.transpose(),
        .D = sys.D,
        .G = d.sys.G,
        .H = sys.H,
        .Ts = sys.Ts,
    };
    out.transform = d.transform;
    out.rank = d.rank;
    return detail::set_status(out, d.status);
}

} // namespace damp
