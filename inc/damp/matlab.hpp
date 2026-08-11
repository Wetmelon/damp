// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file matlab.hpp
 * @brief MATLAB®-style free-function aliases for Damp APIs (host-oriented)
 *
 * Thin spellings over Damp design/analysis/systems. Prefer the descriptive
 * Damp names in new code; these exist for Control System Toolbox muscle memory.
 *
 * Optional returns match the underlying APIs after the systems fail-closed
 * pass: @c c2d / @c place / @c acker / @c pidtune / @c gram / @c lyap / @c dlyap
 * / @c norm return @c damp::optional when the design or conversion can fail.
 * Frequency aliases (@c bode, @c nyquist, @c sigma, @c rlocus, …) skip contour
 * poles honestly (shorter grids) via @c analysis::.
 */

#include <cstddef>
#include <limits>
#include <type_traits>
#include <vector>

#include "damp/analysis/analysis.hpp"
#include "damp/backend.hpp"
#include "damp/controllers/lqg.hpp"
#include "damp/controllers/lqgi.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/design/linearization.hpp"
#include "damp/design/minreal.hpp"
#include "damp/design/pole_placement.hpp"
#include "damp/design/qp.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/filters/iir_design.hpp"
#include "damp/math/complex.hpp"
#include "damp/matrix/eigen.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/matrix/svd.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"
#include "damp/systems/zpk.hpp"
#include "design/stability.hpp"

namespace damp {
// MATLAB®-style matrix functions
namespace matlab {

/**
 * @brief MATLAB®-style transfer function constructor
 *
 * Coefficients are in ascending powers of s or z.
 */
template<size_t Nnum, size_t Nden, typename T = double>
[[nodiscard]] constexpr TransferFunction<Nnum, Nden, T>
tf(const damp::array<T, Nnum>& num, const damp::array<T, Nden>& den) noexcept {
    return TransferFunction<Nnum, Nden, T>{.num = num, .den = den};
}

/**
 * @brief MATLAB®-style transfer function constructor from braced coefficient lists
 *
 * Enables calls like: tf({1.0, 2.0}, {3.0, 4.0, 5.0})
 * Coefficients are in ascending powers of s or z.
 */
template<typename TNum, size_t Nnum, typename TDen, size_t Nden>
[[nodiscard]] constexpr TransferFunction<Nnum, Nden, std::common_type_t<TNum, TDen>>
tf(const TNum (&num)[Nnum], const TDen (&den)[Nden]) noexcept {
    using T = std::common_type_t<TNum, TDen>;
    TransferFunction<Nnum, Nden, T> result{};

    for (size_t i = 0; i < Nnum; ++i) {
        result.num[i] = static_cast<T>(num[i]);
    }
    for (size_t i = 0; i < Nden; ++i) {
        result.den[i] = static_cast<T>(den[i]);
    }

    return result;
}

/**
 * @brief MATLAB®-style state-space model constructor.
 *
 * Builds a StateSpace from (A, B, C, D) with an optional sampling period
 * (Ts = 0 for continuous, > 0 for discrete). Process/measurement noise
 * channels are omitted (NW = NV = 0); construct StateSpace directly for those.
 *
 * @note Compare with MATLAB®'s sys = ss(A, B, C, D) / ss(A, B, C, D, Ts).
 */
template<size_t NX, size_t NU, size_t NY, typename T = double>
[[nodiscard]] constexpr StateSpace<NX, NU, NY, T> ss(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NY, NX, T>& C,
    const Matrix<NY, NU, T>& D = Matrix<NY, NU, T>::zeros(),
    T                        Ts = T{0}
) noexcept {
    return StateSpace<NX, NU, NY, T>{
        .A = A,
        .B = B,
        .C = C,
        .D = D,
        .Ts = Ts,
    };
}

/**
 * @brief MATLAB®-style zero-pole-gain model constructor.
 *
 * @param z Zeros, @param p Poles, @param k Scalar gain.
 * @note Compare with MATLAB®'s sys = zpk(z, p, k).
 */
template<size_t Nz, size_t Np, typename T = double>
[[nodiscard]] constexpr ZPK<Nz, Np, T> zpk(
    const damp::array<damp::complex<T>, Nz>& z,
    const damp::array<damp::complex<T>, Np>& p,
    T                                        k
) noexcept {
    return ZPK<Nz, Np, T>{z, p, k};
}

/**
 * @brief MATLAB®-style TF → ZPK conversion.
 * @note Compare with MATLAB®'s tf2zpk / zpk(tf).
 * @see damp::to_zpk
 */
using damp::tf2zpk;
using damp::to_zpk;

/**
 * @brief MATLAB®-style ZPK → TF conversion.
 * @note Compare with MATLAB®'s zpk2tf / tf(zpk).
 */
using damp::zpk2tf;

/**
 * @brief MATLAB®-style SISO SS → ZPK / SS → TF.
 * @note Compare with MATLAB®'s ss2zpk, ss2tf.
 */
using damp::ss2tf;
using damp::ss2zpk;

/**
 * @brief MATLAB®-style parallel-form continuous PID constructor.
 *
 * Returns a @ref ContinuousPID with the given continuous-time gains; @p Tf is
 * the derivative filter time constant (0 = unfiltered). Deploy fixed-rate with
 * `design::pid(...).discretize(Ts)` → PIDController; this alias matches
 * MATLAB®'s continuous `pid` object before `c2d`.
 *
 * @note Compare with MATLAB®'s C = pid(Kp, Ki, Kd, Tf). MATLAB®'s filter
 *       coefficient N maps to Tf = 1/N.
 */
template<typename T = double>
[[nodiscard]] constexpr ContinuousPID<T> pid(T Kp, T Ki = T{0}, T Kd = T{0}, T Tf = T{0}) noexcept {
    ContinuousPID<T> c{};
    c.Kp = Kp;
    c.Ki = Ki;
    c.Kd = Kd;
    c.Tf = Tf;
    return c;
}

/**
 * @brief Standard-form continuous PID constructor (1-DOF).
 *
 * Maps industrial standard form
 * @f[
 *   C(s) = K_p\Bigl(1 + \frac{1}{T_i s} + \frac{T_d s}{1+(T_d/N)s}\Bigr)
 * @f]
 * to the shipped parallel @ref design::PIDResult with
 * @f$K_i = K_p/T_i@f$, @f$K_d = K_p T_d@f$, @f$T_f = T_d/N@f$ when @p N > 0
 * (else unfiltered). Default @p N = 10 matches MATLAB®'s pidstd when @p Td > 0.
 * @p Ti <= 0 yields @f$K_i = 0@f$ (P/PD). Pass @p N = 0 for unfiltered D.
 *
 * Deploy with `.discretize(Ts)` → PIDController, same as other PID designs.
 *
 * @note Compare with MATLAB®'s C = pidstd(Kp, Ti, Td, N) (default N = 10).
 * @see design::pid for the parallel-form constructor
 * @see pidstd2 for the 2-DOF standard-form variant
 */
template<typename T = double>
[[nodiscard]] constexpr design::PIDResult<T> pidstd(T Kp, T Ti, T Td = T{0}, T N = T{10}) noexcept {
    const T Ki = (Ti > T{0}) ? (Kp / Ti) : T{0};
    const T Kd = Kp * Td;
    const T Tf = (N > T{0} && Td > T{0}) ? (Td / N) : T{0};
    return design::pid(Kp, Ki, Kd, -std::numeric_limits<T>::max(), std::numeric_limits<T>::max(), -std::numeric_limits<T>::max(), std::numeric_limits<T>::max(), T{0}, T{1}, T{1}, Tf);
}

/**
 * @brief Standard-form continuous 2-DOF PID constructor.
 *
 * Same Ti/Td → Ki/Kd map as @ref pidstd, with proportional/derivative setpoint
 * weights @p b and @p c (see @ref design::PIDResult).
 *
 * @note Compare with MATLAB®'s C = pidstd2(Kp, Ti, Td, N, b, c).
 * @see pidstd, make1DOF, make2DOF
 */
template<typename T = double>
[[nodiscard]] constexpr design::PIDResult<T>
pidstd2(T Kp, T Ti, T Td = T{0}, T N = T{10}, T b = T{1}, T c = T{0}) noexcept {
    auto r = pidstd(Kp, Ti, Td, N);
    r.b = b;
    r.c = c;
    return r;
}

/**
 * @brief Force 1-DOF setpoint weights on a PID design result (b=c=1).
 *
 * @note Compare with MATLAB®'s make1DOF (Control System Toolbox).
 */
template<typename T = double>
[[nodiscard]] constexpr design::PIDResult<T> make1DOF(design::PIDResult<T> r) noexcept {
    r.b = T{1};
    r.c = T{1};
    return r;
}

/**
 * @brief Apply 2-DOF setpoint weights on a PID design result.
 *
 * Defaults match the common PI-D industrial form (@f$b=1@f$, @f$c=0@f$: derivative
 * on measurement only).
 *
 * @note Compare with MATLAB®'s make2DOF (Control System Toolbox).
 */
template<typename T = double>
[[nodiscard]] constexpr design::PIDResult<T>
make2DOF(design::PIDResult<T> r, T b = T{1}, T c = T{0}) noexcept {
    r.b = b;
    r.c = c;
    return r;
}

/**
 * @brief MATLAB® short alias for design::minreal
 *
 * @note Compare with MATLAB®'s minreal(sys, tol).
 * @see design::minreal(), design::minimal_realization()
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr design::MinimalRealizationResult<NX, NU, NY, T, NW, NV>
minreal(const StateSpace<NX, NU, NY, T, NW, NV>& sys, T tol = default_tol<T>()) {
    return design::minreal(sys, tol);
}

/**
 * @brief MATLAB® short alias for controllability_matrix
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr Matrix<NX, NX * NU, T> ctrb(const Matrix<NX, NX, T>& A, const Matrix<NX, NU, T>& B) noexcept {
    return stability::controllability_matrix(A, B);
}

/**
 * @brief MATLAB® short alias for controllability_matrix taking StateSpace
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr Matrix<NX, NX * NU, T> ctrb(const StateSpace<NX, NU, NY, T, NW, NV>& sys) noexcept {
    return stability::controllability_matrix(sys.A, sys.B);
}

/**
 * @brief MATLAB® short alias for observability_matrix
 */
template<size_t NX, size_t NY, typename T = double>
[[nodiscard]] constexpr Matrix<NX * NY, NX, T> obsv(const Matrix<NX, NX, T>& A, const Matrix<NY, NX, T>& C) noexcept {
    return stability::observability_matrix(A, C);
}

/**
 * @brief MATLAB® short alias for observability_matrix
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr Matrix<NX * NY, NX, T> obsv(const StateSpace<NX, NU, NY, T, NW, NV>& sys) noexcept {
    return stability::observability_matrix(sys.A, sys.C);
}

/**
 * @brief MATLAB® interface function c2d to discretize a continuous-time state-space system
 *
 * @param sys           Continuous-time state-space model (Ts should be 0)
 * @param sampling_time Desired sampling period for discrete system
 * @param method        Discretization method (ZOH or Tustin)
 *
 * @return Discrete StateSpace, or @c nullopt if the method fails (e.g. singular Tustin)
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr damp::optional<StateSpace<NX, NU, NY, T, NW, NV>> c2d(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    T                                        sampling_time,
    DiscretizationMethod                     method = DiscretizationMethod::ZOH
) {
    return discretize(sys, sampling_time, method);
}

/**
 * @brief MATLAB®-style c2d for SISO transfer functions
 *
 * Realizes the TF then discretizes. Returns @c nullopt if realization or c2d fails.
 */
template<size_t Nnum, size_t Nden, typename T = double>
[[nodiscard]] constexpr damp::optional<StateSpace<Nden - 1, 1, 1, T>> c2d(
    const TransferFunction<Nnum, Nden, T>& tf_sys,
    T                                      sampling_time,
    DiscretizationMethod                   method = DiscretizationMethod::ZOH
) {
    const auto ss = tf_sys.to_state_space();
    if (!ss) {
        return damp::nullopt;
    }
    return discretize(*ss, sampling_time, method);
}

/**
 * @brief MATLAB®-style nonlinear linearization about an operating point
 *
 * Similar in spirit to MATLAB®'s linmod workflow: linearize nonlinear
 * dynamics and output maps, then return continuous-time state-space.
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, typename Dynamics, typename Output>
[[nodiscard]] constexpr StateSpace<NX, NU, NY, T> linmod(
    const Dynamics&      dynamics,
    const Output&        output,
    const ColVec<NX, T>& x_op,
    const ColVec<NU, T>& u_op,
    T                    epsilon = T{0}
) {
    if (epsilon > T{0}) {
        return design::linearize<NX, NU, NY, T>(dynamics, output, x_op, u_op, epsilon).to_state_space();
    }
    return design::linearize<NX, NU, NY, T>(dynamics, output, x_op, u_op).to_state_space();
}

/**
 * @brief Block diagonal matrix construction
 * @ingroup linear_algebra
 * @tparam Args Variadic list of matrix types
 * @return Block diagonal matrix composed of the input matrices
 */
template<typename... Args>
[[nodiscard]] constexpr auto blkdiag(Args... args) noexcept
    requires(sizeof...(Args) > 0 && (is_matrix_type<Args>::value && ...))
{
    constexpr size_t total_rows = (args.rows() + ...);
    constexpr size_t total_cols = (args.cols() + ...);
    using T = std::common_type_t<typename std::decay_t<Args>::value_type...>;
    Matrix<total_rows, total_cols, T> result = Matrix<total_rows, total_cols, T>::zeros();

    size_t row_offset = 0;
    size_t col_offset = 0;
    ((result.template block<args.rows(), args.cols()>(row_offset, col_offset) = args, row_offset += args.rows(), col_offset += args.cols()), ...);

    return result;
}

/**
 * @brief Returns a square diagonal matrix from the given array
 *
 * @tparam k Diagonal offset: 0 main diagonal, &gt;0 superdiagonal, &lt;0 subdiagonal
 * @param diag Diagonal entries
 * @return Matrix&lt;N, N, T&gt;
 */
template<size_t N, typename T, int k = 0>
[[nodiscard]] constexpr auto diag(const damp::array<T, N>& diag) noexcept {
    Matrix<N, N, T> result = Matrix<N, N, T>::zeros();
    for (size_t i = 0; i < N; ++i) {
        size_t row = i;
        size_t col = i + k;
        if (col < N) {
            result(row, col) = diag[i];
        }
    }
    return result;
}

/**
 * @brief Returns the diagonal elements of a square matrix as a column vector
 *
 * @param A Input square matrix
 *
 * @return ColVec<N, T>
 */
template<size_t N, typename T, int k = 0>
[[nodiscard]] constexpr ColVec<N, T> diag(const Matrix<N, N, T>& A) noexcept {
    ColVec<N, T> result;
    for (size_t i = 0; i < N; ++i) {
        size_t row = i;
        size_t col = i + k;
        if (col < N) {
            result(i, 0) = A(row, col);
        } else {
            result(i, 0) = T{0};
        }
    }
    return result;
}

/**
 * @brief Create an identity matrix of size n x n
 *
 * @tparam N Size of the identity matrix
 * @return constexpr Matrix<N, N, T>
 */
template<size_t N, typename T = double>
[[nodiscard]] constexpr auto eye() noexcept {
    return Matrix<N, N, T>::identity();
}

template<size_t N, typename T = double>
[[nodiscard]]
constexpr auto ones() noexcept {
    Matrix<N, N, T> result;
    for (auto& row : result.data_) {
        row.fill(T{1});
    }
    return result;
}

/**
 * @brief MATLAB® short alias for the singular value decomposition.
 *
 * Returns a mat::SVDResult{singular_U, singular_values, singular_V} rather than
 * the MATLAB® [U, S, V] tuple; the singular values are a descending array, not a
 * diagonal matrix.
 *
 * @note Compare with MATLAB®'s [U, S, V] = svd(A).
 */
template<size_t M, size_t N, typename T>
[[nodiscard]] constexpr auto svd(const Matrix<M, N, T>& A) {
    return mat::svd(A);
}

/**
 * @brief MATLAB® short alias for the Moore–Penrose pseudoinverse.
 * @note Compare with MATLAB®'s pinv(A).
 */
template<size_t M, size_t N, typename T>
[[nodiscard]] constexpr Matrix<N, M, T> pinv(const Matrix<M, N, T>& A) {
    return mat::pseudo_inverse(A);
}

/**
 * @brief MATLAB® short alias for an orthonormal null-space basis.
 *
 * Returns a mat::NullSpace whose trailing `dim` columns are the kernel basis
 * (MATLAB®'s null returns those columns directly as an N×dim matrix).
 *
 * @note Compare with MATLAB®'s null(A).
 */
template<size_t M, size_t N, typename T>
[[nodiscard]] constexpr mat::NullSpace<N, T> null(const Matrix<M, N, T>& A) {
    return mat::null_space(A);
}

/**
 * @brief MATLAB® short alias for the eigenvalues of a square matrix.
 *
 * Returns the eigenvalues as a complex column vector. Closed-form for N ≤ 4,
 * Francis double-shift QR for larger systems. For eigenvectors and the
 * convergence flag, call mat::compute_eigenvalues() directly.
 *
 * @note Compare with MATLAB®'s e = eig(A).
 */
template<size_t N, typename T>
[[nodiscard]] constexpr ColVec<N, damp::complex<T>> eig(const Matrix<N, N, T>& A) {
    return mat::compute_eigenvalues(A).values;
}

/**
 * @brief Form state estimator from system and estimator gain
 *
 * est = estim(sys,L) produces a state/output estimator est given the plant state-space model sys
 *       and the estimator gain L. All inputs w of sys are assumed stochastic (process and/or measurement noise),
 *       and all outputs y are measured. The estimator est is returned in state-space form.
 *
 * @param sys State-space system (discrete or continuous)
 * @param L   Estimator gain matrix
 *
 * @return StateSpace<NX, NY, NX, T, NW, NV>
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
constexpr auto estim(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const Matrix<NX, NY, T>& L) noexcept {
    Matrix<NX, NX, T> A_est = sys.A - L * sys.C;
    Matrix<NX, NY, T> B_est = L;
    Matrix<NX, NX, T> C_est = Matrix<NX, NX, T>::identity();
    Matrix<NX, NY, T> D_est = Matrix<NX, NY, T>::zeros();
    Matrix<NX, NW, T> G_est = sys.G - L * sys.H;
    Matrix<NX, NV, T> H_est = L;

    return StateSpace<NX, NY, NX, T, NW, NV>{A_est, B_est, C_est, D_est, G_est, H_est};
}

/**
 * @brief Form dynamic regulator from system, state-feedback gain, and estimator gain
 *
 * rsys = reg(sys,K,L) forms a dynamic regulator or compensator rsys given a state-space model sys of the plant,
 *         a state-feedback gain matrix K, and an estimator gain matrix L.
 *         The gains K and L are typically designed using pole placement or LQG techniques.
 *         The function reg handles both continuous- and discrete-time cases.
 *
 *    This syntax assumes that all inputs of sys are controls, and all outputs are measured.
 *    The regulator rsys is obtained by connecting the state-feedback law u = –Kx and the state estimator with gain matrix L (see estim).
 *
 * @param sys State-space system (discrete or continuous)
 * @param K   State-feedback gain matrix
 * @param L   Estimator gain matrix
 *
 * @return StateSpace<2*NX, NW+NV, NY, T>
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
constexpr auto reg(const StateSpace<NX, NU, NY, T, NW, NV>& sys, const Matrix<NU, NX, T>& K, const Matrix<NX, NY, T>& L) noexcept {
    constexpr size_t N_reg = 2 * NX;
    constexpr size_t M_reg = NW + NV;
    constexpr size_t P_reg = NY;

    Matrix<N_reg, N_reg, T> A_reg = Matrix<N_reg, N_reg, T>::zeros();
    A_reg.template block<NX, NX>(0, 0) = sys.A;
    A_reg.template block<NX, NX>(0, NX) = -sys.B * K;
    A_reg.template block<NX, NX>(NX, 0) = L * sys.C;
    A_reg.template block<NX, NX>(NX, NX) = sys.A - L * sys.C - L * sys.D * K;

    Matrix<N_reg, M_reg, T> B_reg = Matrix<N_reg, M_reg, T>::zeros();
    B_reg.template block<NX, NW>(0, 0) = sys.G;
    B_reg.template block<NX, NV>(NX, NW) = L * sys.H;

    Matrix<P_reg, N_reg, T> C_reg = Matrix<P_reg, N_reg, T>::zeros();
    C_reg.template block<NY, NX>(0, 0) = sys.C;
    C_reg.template block<NY, NX>(0, NX) = -sys.D * K;

    Matrix<P_reg, M_reg, T> D_reg = Matrix<P_reg, M_reg, T>::zeros();
    D_reg.template block<NY, NV>(0, NW) = sys.H;

    return StateSpace<N_reg, M_reg, P_reg, T>{A_reg, B_reg, C_reg, D_reg};
}

/**
 * @brief Pole placement for state-feedback control
 *
 * @param A State matrix
 * @param B Input matrix
 * @param p Desired poles (damp::array / std::array of damp::complex or std::complex)
 *
 * @return damp::optional<Matrix<NU, NX, T>> State-feedback gain K, or nullopt if not implementable
 */
template<size_t NX, size_t NU, typename T = double>
constexpr damp::optional<Matrix<NU, NX, T>> acker(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const auto&              p
) {
    static_assert(NU == 1, "acker is single-input only; use place for multi-input systems");
    const auto poles = std::apply(
        [](const auto&... elems) { return damp::array<damp::complex<T>, sizeof...(elems)>{elems...}; }, p
    );
    return design::ackermann(A, B, poles);
}

/**
 * @brief Robust multi-input pole placement (MATLAB®'s place)
 *
 * Thin alias for design::place — Kautsky–Nichols–Van Dooren robust eigenvalue
 * assignment, spending the multi-input freedom to minimize eigenvector
 * conditioning. Prefer this for multi-input systems; acker remains the
 * single-input Ackermann path.
 *
 * @param A State matrix
 * @param B Input matrix
 * @param p Desired poles (damp::array / std::array of damp::complex or std::complex)
 * @return State-feedback gain K (NU×NX), or nullopt if not assignable.
 *
 * @note Compare with MATLAB®'s K = place(A, B, p).
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NU, NX, T>> place(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const auto&              p
) {
    const auto poles = std::apply(
        [](const auto&... elems) { return damp::array<damp::complex<T>, sizeof...(elems)>{elems...}; }, p
    );
    return design::place(A, B, poles);
}

/**
 * @brief Continuous-time LQR design (MATLAB®'s lqr)
 *
 * Thin alias for design::continuous_lqr — the optimal gain K for u = −Kx
 * minimizing ∫(xᵀQx + uᵀRu + 2xᵀNu) dt via the continuous ARE (care()).
 *
 * @note Compare with MATLAB®'s K = lqr(A, B, Q, R, N). For the Riccati solution
 *       S and the closed-loop poles, call design::continuous_lqr directly.
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr Matrix<NU, NX, T> lqr(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    const Matrix<NX, NU, T>& N = Matrix<NX, NU, T>{}
) {
    return design::continuous_lqr(A, B, Q, R, N).K; // Only return K
}

/**
 * @brief Discrete-time Linear-Quadratic Regulator design
 * @note Alias for design::dlqr (same as design::discrete_lqr). Compare with MATLAB®'s dlqr.
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr auto dlqr(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    const Matrix<NX, NU, T>& N = Matrix<NX, NU, T>{}
) {
    return design::dlqr(A, B, Q, R, N);
}

/**
 * @brief Design discrete LQR from continuous-time system via discretization
 * @note Alias for design::lqrd (same as design::discrete_lqr_from_continuous). Compare with MATLAB®'s lqrd.
 */
template<size_t NX, size_t NU, typename T = double>
[[nodiscard]] constexpr auto lqrd(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NX, NX, T>& Q,
    const Matrix<NU, NU, T>& R,
    T                        Ts,
    const Matrix<NX, NU, T>& N = Matrix<NX, NU, T>{}
) {
    return design::lqrd(A, B, Q, R, Ts, N);
}

/**
 * @brief Design discrete LQR from continuous state-space system via discretization
 * @note Alias for design::lqrd. Compare with MATLAB®'s lqrd(sys, Q, R, Ts).
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr auto lqrd(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NX, NX, T>&                 Q,
    const Matrix<NU, NU, T>&                 R,
    T                                        Ts,
    const Matrix<NX, NU, T>&                 N = Matrix<NX, NU, T>{}
) {
    return design::lqrd(sys, Q, R, Ts, N);
}

// kalmd (discrete Kalman from continuous plant + continuous noise) is deferred:
// design::kalman requires a discrete plant, and there is no Van Loan process-
// noise discretizer yet. Use discretize(sys, Ts) then design::kalman when Q/R
// are already discrete intensities. See roadmap #29.

/**
 * @brief Output-weighted continuous LQR (state cost Q = Cᵀ Q_y C).
 *
 * Builds the equivalent state-feedback LQR problem from output weights:
 * @f[
 *   Q_x = C^\top Q_y C,\quad
 *   N = C^\top Q_y D,\quad
 *   R_u = R + D^\top Q_y D
 * @f]
 * then calls @ref design::continuous_lqr. When @p D is zero this reduces to
 * pure output weighting @f$Q_x = C^\top Q_y C@f$.
 *
 * @note Compare with MATLAB®'s [K,S,e] = lqry(A, B, C, D, Qy, R).
 * @see design::continuous_lqr, dlqr, lqrd
 */
template<size_t NX, size_t NU, size_t NY, typename T = double>
[[nodiscard]] constexpr design::LQRResult<NX, NU, T> lqry(
    const Matrix<NX, NX, T>& A,
    const Matrix<NX, NU, T>& B,
    const Matrix<NY, NX, T>& C,
    const Matrix<NY, NU, T>& D,
    const Matrix<NY, NY, T>& Qy,
    const Matrix<NU, NU, T>& R
) {
    const Matrix<NX, NX, T> Qx = C.transpose() * Qy * C;
    const Matrix<NX, NU, T> N = C.transpose() * Qy * D;
    const Matrix<NU, NU, T> Ru = R + (D.transpose() * Qy * D);
    return design::continuous_lqr(A, B, Qx, Ru, N);
}

/**
 * @brief Output-weighted LQR for a state-space plant.
 *
 * Continuous plants (@p sys.Ts == 0) use @ref design::continuous_lqr; discrete
 * plants use @ref design::discrete_lqr. Same @f$Q_x/N/R_u@f$ map as the
 * matrix overload.
 *
 * @note Compare with MATLAB®'s [K,S,e] = lqry(sys, Qy, R).
 * @see lqry(A,B,C,D,Qy,R), design::discrete_lqr, design::continuous_lqr
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr design::LQRResult<NX, NU, T> lqry(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NY, NY, T>&                 Qy,
    const Matrix<NU, NU, T>&                 R
) {
    const Matrix<NX, NX, T> Qx = sys.C.transpose() * Qy * sys.C;
    const Matrix<NX, NU, T> N = sys.C.transpose() * Qy * sys.D;
    const Matrix<NU, NU, T> Ru = R + (sys.D.transpose() * Qy * sys.D);
    if (sys.is_discrete()) {
        return design::discrete_lqr(sys.A, sys.B, Qx, Ru, N);
    }
    return design::continuous_lqr(sys.A, sys.B, Qx, Ru, N);
}

/**
 * @brief Linear-Quadratic Integral design for tracking
 * @note Alias for design::discrete_lqi. Compare with MATLAB®'s lqi(sys, Q, R).
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr auto lqi(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NX + NY, NX + NY, T>&       Q,
    const Matrix<NU, NU, T>&                 R
) {
    return design::discrete_lqi(sys, Q, R);
}

/**
 * @brief Linear-Quadratic-Gaussian regulator design
 * @note Alias for design::discrete_lqg. Compare with MATLAB®'s lqg(sys, ...).
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr auto lqg(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NX, NX, T>&                 Q_lqr,
    const Matrix<NU, NU, T>&                 R_lqr,
    const Matrix<NW, NW, T>&                 Q_kf,
    const Matrix<NV, NV, T>&                 R_kf,
    const Matrix<NX, NU, T>&                 N = Matrix<NX, NU, T>{}
) {
    return design::discrete_lqg(sys, Q_lqr, R_lqr, Q_kf, R_kf, N);
}

/**
 * @brief Combine separate Kalman filter and LQR designs into an LQG controller
 * @note Alias for design::lqg_from_parts. Compare with MATLAB®'s lqgreg(kest, k).
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = NX, size_t NV = NY>
[[nodiscard]] constexpr auto lqgreg(
    const design::KalmanResult<NX, NU, NY, T, NW, NV>& kest,
    const design::LQRResult<NX, NU, T>&                lqr_result
) {
    return design::lqg_from_parts(kest, lqr_result);
}

/**
 * @brief Linear-Quadratic-Gaussian design with integral action for tracking
 * @note Alias for design::discrete_lqgi. Compare with MATLAB®'s lqgtrack(...).
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr auto lqgtrack(
    const StateSpace<NX, NU, NY, T, NW, NV>& sys,
    const Matrix<NX + NY, NX + NY, T>&       Q_aug,
    const Matrix<NU, NU, T>&                 R,
    const Matrix<NW, NW, T>&                 Q_kf,
    const Matrix<NV, NV, T>&                 R_kf
) {
    return design::discrete_lqgi(sys, Q_aug, R, Q_kf, R_kf);
}

/**
 * @brief PID controller tuning using frequency domain method
 *
 * Tunes a PID controller for a given plant to achieve a specified crossover frequency wc.
 * Uses the method similar to MATLAB®'s pidtune, aiming for 60 degrees phase margin.
 *
 * @param sys Plant state-space system (SISO, continuous-time)
 * @param wc Desired crossover frequency (rad/s)
 * @return PIDResult with tuned gains
 */
template<size_t NX, typename T = double>
[[nodiscard]] constexpr damp::optional<damp::design::PIDResult<T>>
pidtune(const StateSpace<NX, 1, 1, T>& sys, T wc) noexcept {
    using Cplx = damp::complex<T>;
    constexpr T pi = damp::numbers::pi_v<T>;
    Cplx        jwc{0, wc};
    auto        G_frf_opt = eval_frf(sys, jwc);
    if (!G_frf_opt) {
        return damp::nullopt;
    }
    Cplx G = (*G_frf_opt)(0, 0);
    T    mag_G = damp::abs(G);
    if (!(mag_G > T{0})) {
        return damp::nullopt;
    }
    T arg_G = damp::arg(G);
    // Desired phase margin: 60 degrees = pi/3 radians
    T desired_phase = -pi + pi / 3 - arg_G;
    T mag_C = T{1} / mag_G;
    const auto [s, c] = damp::sincos(desired_phase);
    T real_C = mag_C * c;
    T imag_C = mag_C * s;
    T Kp = real_C;
    T Ki = -imag_C * wc; // imag_C = -Ki/wc for PI
    // For PID, set Td = Ti/4
    T Ti = Kp / Ki;
    T Td = Ti / T{4};
    T Kd = Kp * Td;
    T Kbc = Ki; // Back-calculation gain

    return damp::design::PIDResult<T>{
        Kp, Ki, Kd, T{0},
        -std::numeric_limits<T>::max(), std::numeric_limits<T>::max(),
        -std::numeric_limits<T>::max(), std::numeric_limits<T>::max(),
        Kbc
    };
}

// ===========================================================================
// Frequency-domain analysis — MATLAB® spellings over damp::analysis
// ===========================================================================

// These already carry their MATLAB® names in analysis::; surface them under
// matlab:: too so a MATLAB®-style call site finds them in one namespace.
using analysis::arange;
using analysis::bode;
using analysis::damp;
using analysis::dcgain;
using analysis::geomspace;
using analysis::impulse;
using analysis::initial;
using analysis::linspace;
using analysis::logspace;
using analysis::lsim;
using analysis::lsiminfo;
using analysis::nichols;
using analysis::nyquist;
using analysis::pzmap;
using analysis::rlocus;
using analysis::sigma;
using analysis::step;
using analysis::stepinfo;

/**
 * @brief MATLAB® short alias for the open-loop poles of a system.
 * @note Compare with MATLAB®'s p = pole(sys). (analysis names it `poles`.)
 */
template<size_t NX, typename T = double>
[[nodiscard]] constexpr ColVec<NX, damp::complex<T>> pole(const Matrix<NX, NX, T>& A) {
    return analysis::poles(A);
}

/**
 * @brief Gain/phase margins and their crossover frequencies.
 *
 * @note Compare with MATLAB®'s [Gm, Pm, Wcg, Wcp] = margin(...). Gm is a linear
 *       ratio; +inf marks a missing crossover.
 *
 * @c Wcg is the phase-crossover frequency (GM site); @c Wcp is the
 * gain-crossover frequency (PM site).
 */
template<typename T = double>
struct MarginResult {
    T Gm{};  ///< Gain margin (linear ratio); +inf if phase never crosses -180°
    T Pm{};  ///< Phase margin (degrees); +inf if gain never crosses 0 dB
    T Wcg{}; ///< Phase-crossover frequency (rad/s, where phase = -180°) — GM site
    T Wcp{}; ///< Gain-crossover frequency (rad/s, where |G| = 1) — PM site
};

/**
 * @brief Gain and phase margins of a SISO loop over a frequency grid.
 *
 * Thin composition of analysis::bode and BodeResult::gain_margin /
 * phase_margin. Unlike MATLAB®'s margin(sys), the grid is explicit — pass
 * e.g. matlab::logspace(...) — since damp does no automatic frequency gridding.
 *
 * @param sys   SISO state-space loop (continuous or discrete)
 * @param omega Frequency grid (rad/s)
 * @return MarginResult{Gm, Pm, Wcg, Wcp}
 *
 * @note Compare with MATLAB®'s [Gm, Pm, Wcg, Wcp] = margin(sys, w).
 */
template<size_t NX, size_t NW, size_t NV, typename T>
[[nodiscard]] MarginResult<T> margin(const StateSpace<NX, 1, 1, T, NW, NV>& sys, const std::vector<T>& omega) {
    const auto      bode_data = analysis::bode(sys, omega);
    constexpr T     inf = std::numeric_limits<T>::max();
    MarginResult<T> m{inf, inf, T{0}, T{0}};
    if (const auto gm = bode_data.gain_margin()) {
        m.Gm = damp::pow(T{10}, gm->first / T{20}); // dB -> linear ratio
        m.Wcg = gm->second;
    }
    if (const auto pm = bode_data.phase_margin()) {
        m.Pm = pm->first;
        m.Wcp = pm->second;
    }
    return m;
}

/**
 * @brief All classical margins including delay margin (superset of @ref margin).
 *
 * Delay margin @f$D_m = P_m\mathrm{(rad)} / \omega_{gc}@f$ where @f$\omega_{gc}@f$
 * is the gain-crossover frequency stored in @ref MarginResult::Wcp. Missing
 * crossovers keep the corresponding margin at +inf (same sentinel as @ref margin).
 *
 * @note Compare with MATLAB®'s S = allmargin(sys) (subset: first GM/PM/DM only;
 *       multi-crossing arrays are not collected).
 * @see margin
 */
template<typename T = double>
struct AllMarginResult {
    T Gm{};  ///< Gain margin (linear ratio); +inf if no phase crossover
    T Pm{};  ///< Phase margin (degrees); +inf if no gain crossover
    T Dm{};  ///< Delay margin (seconds); +inf if no gain crossover or ω = 0
    T Wcg{}; ///< Phase-crossover frequency (rad/s) — GM site
    T Wcp{}; ///< Gain-crossover frequency (rad/s) — PM / DM site
};

/**
 * @brief Gain, phase, and delay margins of a SISO loop over a frequency grid.
 *
 * Thin extension of margin: delay margin is the classical conversion of
 * phase margin at the gain crossover. Grid is explicit (pass @ref logspace).
 *
 * @param sys   SISO state-space loop (continuous or discrete)
 * @param omega Frequency grid (rad/s)
 * @return AllMarginResult{Gm, Pm, Dm, Wcg, Wcp}
 *
 * @note Compare with MATLAB®'s S = allmargin(sys).
 */
template<size_t NX, size_t NW, size_t NV, typename T>
[[nodiscard]] AllMarginResult<T> allmargin(const StateSpace<NX, 1, 1, T, NW, NV>& sys, const std::vector<T>& omega) {
    const MarginResult<T> m = margin(sys, omega);
    constexpr T           inf = std::numeric_limits<T>::max();
    AllMarginResult<T>    a{m.Gm, m.Pm, inf, m.Wcg, m.Wcp};
    if (m.Pm < inf && m.Wcp > T{0}) {
        a.Dm = (m.Pm * damp::numbers::pi_v<T> / T{180}) / m.Wcp;
    }
    return a;
}

/**
 * @brief Continuous-time stability predicate on a state matrix.
 *
 * True when all eigenvalues of @p A satisfy @f$\mathrm{Re}(\lambda) < 0@f$.
 *
 * @note Compare with MATLAB®'s isstable(sys) for continuous models.
 * @see analysis::is_stable_continuous, stability::is_stable_discrete
 */
template<size_t NX, typename T = double>
[[nodiscard]] constexpr bool isstable(const Matrix<NX, NX, T>& A) {
    return analysis::is_stable_continuous(A);
}

/**
 * @brief Stability predicate for a state-space model.
 *
 * Continuous (@p sys.Ts == 0): left half-plane. Discrete (@p sys.Ts > 0): unit disk.
 *
 * @note Compare with MATLAB®'s tf = isstable(sys).
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr bool isstable(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    if (sys.is_discrete()) {
        return stability::is_stable_discrete(sys.A);
    }
    return analysis::is_stable_continuous(sys.A);
}

/**
 * @brief Stability predicate for a continuous transfer function.
 *
 * Realizes @p tf as companion state-space and tests the poles.
 *
 * @note Compare with MATLAB®'s isstable(sys).
 */
template<size_t Nnum, size_t Nden, typename T = double>
[[nodiscard]] constexpr bool isstable(const TransferFunction<Nnum, Nden, T>& tf) {
    const auto ss = tf.to_state_space();
    return ss ? isstable(*ss) : false;
}

/**
 * @brief First-order Padé approximation of pure delay e^{−sT}.
 *
 * Thin alias for @ref design::pade_delay_1st.
 *
 * @note Compare with MATLAB®'s pade(T, 1) (numerator/denominator of the
 *       rational approximation; Damp returns a TransferFunction).
 * @see pade2, design::pade_delay_1st
 */
template<typename T = double>
[[nodiscard]] constexpr TransferFunction<2, 2, T> pade(T T_delay) {
    return design::pade_delay_1st(T_delay);
}

/**
 * @brief First-order Padé delay, Tustin-discretized at @p Ts.
 * @note Compare with MATLAB®'s c2d(pade(T,1), Ts, 'tustin') path.
 */
template<typename T = double>
[[nodiscard]] constexpr TransferFunction<2, 2, T> pade(T T_delay, T Ts) {
    return design::pade_delay_1st(T_delay, Ts);
}

/**
 * @brief Second-order Padé approximation of pure delay e^{−sT}.
 *
 * Thin alias for @ref design::pade_delay_2nd.
 *
 * @note Compare with MATLAB®'s pade(T, 2).
 * @see pade, design::pade_delay_2nd
 */
template<typename T = double>
[[nodiscard]] constexpr TransferFunction<3, 3, T> pade2(T T_delay) {
    return design::pade_delay_2nd(T_delay);
}

/**
 * @brief Second-order Padé delay, Tustin-discretized at @p Ts.
 */
template<typename T = double>
[[nodiscard]] constexpr TransferFunction<3, 3, T> pade2(T T_delay, T Ts) {
    return design::pade_delay_2nd(T_delay, Ts);
}

/**
 * @brief -3 dB bandwidth of a SISO system over a frequency grid.
 *
 * @param sys   SISO state-space system
 * @param omega Frequency grid (rad/s), ascending from ~DC
 * @return Bandwidth (rad/s), or nullopt if the response never drops 3 dB.
 *
 * @note Compare with MATLAB®'s fb = bandwidth(sys). Grid is explicit here.
 */
template<size_t NX, size_t NW, size_t NV, typename T>
[[nodiscard]] damp::optional<T> bandwidth(const StateSpace<NX, 1, 1, T, NW, NV>& sys, const std::vector<T>& omega) {
    return analysis::bode(sys, omega).bandwidth();
}

/**
 * @brief MATLAB® alias for the continuous Lyapunov solve AX+XAᵀ+Q=0.
 * @note Compare with MATLAB®'s @c X=lyap(A,Q). Returns nullopt if no unique solution.
 */
template<size_t NX, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NX, T>>
lyap(const Matrix<NX, NX, T>& A, const Matrix<NX, NX, T>& Q) {
    return damp::lyap(A, Q);
}

/**
 * @brief MATLAB® alias for the discrete Lyapunov solve AXAᵀ−X+Q=0.
 * @note Compare with MATLAB®'s @c X=dlyap(A,Q). Returns nullopt if no unique solution.
 */
template<size_t NX, typename T = double>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NX, T>>
dlyap(const Matrix<NX, NX, T>& A, const Matrix<NX, NX, T>& Q) {
    return damp::dlyap(A, Q);
}

/**
 * @brief MATLAB® alias for the dense inequality-constrained QP solve.
 * @note Compare with MATLAB®'s @c x=quadprog(H,f,A,b). Returns only the
 *       minimizer; for multipliers, objective, and status call
 *       design::solve_qp directly.
 */
template<size_t NV, size_t NI, typename T = double>
[[nodiscard]] constexpr ColVec<NV, T> quadprog(
    const Matrix<NV, NV, T>& H,
    const ColVec<NV, T>&     f,
    const Matrix<NI, NV, T>& A,
    const ColVec<NI, T>&     b
) {
    return design::solve_qp(H, f, A, b).x; // Only return x
}

/**
 * @brief MATLAB® alias for the equality-constrained QP solve.
 * @note Compare with MATLAB®'s @c x=quadprog(H,f,A,b,Aeq,beq). Returns only the
 *       minimizer; call design::solve_qp directly for multipliers and status.
 */
template<size_t NV, size_t NI, size_t NE, typename T = double>
[[nodiscard]] constexpr ColVec<NV, T> quadprog(
    const Matrix<NV, NV, T>& H,
    const ColVec<NV, T>&     f,
    const Matrix<NI, NV, T>& A,
    const ColVec<NI, T>&     b,
    const Matrix<NE, NV, T>& Aeq,
    const ColVec<NE, T>&     beq
) {
    return design::solve_qp(H, f, A, b, Aeq, beq).x; // Only return x
}

/**
 * @brief MATLAB® alias for the controllability/observability Gramian of a system.
 *
 * @param sys  State-space system (continuous or discrete, dispatched on Ts).
 * @param kind @c 'c' for the controllability Gramian, @c 'o' for observability.
 * @return Gramian (NX × NX), or nullopt if the Lyapunov solve fails.
 * @note Compare with MATLAB®'s @c gram(sys,'c') / @c gram(sys,'o').
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr damp::optional<Matrix<NX, NX, T>>
gram(const StateSpace<NX, NU, NY, T, NW, NV>& sys, char kind = 'c') {
    if (kind == 'o' || kind == 'O') {
        return stability::observability_gramian(sys.A, sys.C, sys.is_discrete());
    }
    return stability::controllability_gramian(sys.A, sys.B, sys.is_discrete());
}

/**
 * @brief MATLAB® alias for the H2 system norm @c norm(sys,2).
 * @note Returns nullopt if infinite (non-strictly-proper continuous system) or the
 *       Gramian solve fails. @see analysis::norm_h2.
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] constexpr damp::optional<T> norm(const StateSpace<NX, NU, NY, T, NW, NV>& sys) {
    return analysis::norm_h2(sys);
}

/**
 * @brief MATLAB® alias for the H∞ system norm @c norm(sys,Inf) / @c hinfnorm(sys).
 * @note Frequency-sweep estimate; see analysis::norm_hinf for the accuracy ceiling.
 */
template<size_t NX, size_t NU, size_t NY, typename T = double, size_t NW = 0, size_t NV = 0>
[[nodiscard]] T hinfnorm(const StateSpace<NX, NU, NY, T, NW, NV>& sys, size_t n_points = 1024) {
    return analysis::norm_hinf(sys, n_points);
}

} // namespace matlab

} // namespace damp