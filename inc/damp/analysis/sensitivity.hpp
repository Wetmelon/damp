// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file sensitivity.hpp
 * @brief Gang of four and gang of seven: the closed-loop maps of a plant and its controller, as state-space models (host-only)
 */

#include <cstddef>

#include "damp/backend.hpp"
#include "damp/design/stability.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {
namespace analysis {

/**
 * @brief Gang of four of a loop u = C(r − y): S, T, PS and CS, each a state-space model on the closed loop's states
 *
 * Output-side forms for a MIMO loop. Every member shares the closed loop's A, so @c stable holds for all four.
 */
template<size_t NX, size_t NU, size_t NY, typename T = double>
struct [[nodiscard]] GangOfFourResult {
    StateSpace<NX, NY, NY, T> sensitivity{};               ///< S = (I + PC)⁻¹: r → e, and an output disturbance → y
    StateSpace<NX, NY, NY, T> complementary_sensitivity{}; ///< T = PC(I + PC)⁻¹: r → y, and measurement noise → −y
    StateSpace<NX, NU, NY, T> load_sensitivity{};          ///< PS = P(I + CP)⁻¹: a load disturbance at the plant input → y
    StateSpace<NX, NY, NU, T> noise_sensitivity{};         ///< CS = C(I + PC)⁻¹: measurement noise → −u
    bool                      stable{false};               ///< Closed loop stable
    bool                      success{false};              ///< Sample times match and the loop is well posed
};

/**
 * @brief Gang of seven of a two-degree-of-freedom loop u = K[r; y]: the gang of four of its feedback path plus the reference responses
 *
 * With the controller split as u = C_r r − C_y y, the four are those of C_y. For a one-degree-of-freedom
 * controller (C_r = C_y) the reference responses equal T, CS and S.
 */
template<size_t NX, size_t NU, size_t NY, typename T = double>
struct [[nodiscard]] GangOfSevenResult {
    StateSpace<NX, NY, NY, T> sensitivity{};               ///< S = (I + PC_y)⁻¹: an output disturbance → y
    StateSpace<NX, NY, NY, T> complementary_sensitivity{}; ///< T = PC_y(I + PC_y)⁻¹: measurement noise → −y
    StateSpace<NX, NU, NY, T> load_sensitivity{};          ///< PS = P(I + C_y P)⁻¹: a load disturbance at the plant input → y
    StateSpace<NX, NY, NU, T> noise_sensitivity{};         ///< CS = C_y(I + PC_y)⁻¹: measurement noise → −u
    StateSpace<NX, NY, NY, T> reference_to_output{};       ///< r → y: (I + PC_y)⁻¹ P C_r
    StateSpace<NX, NY, NU, T> reference_to_control{};      ///< r → u: (I + C_y P)⁻¹ C_r
    StateSpace<NX, NY, NY, T> reference_to_error{};        ///< r → e = r − y: I minus reference_to_output
    bool                      stable{false};               ///< Closed loop stable
    bool                      success{false};              ///< Sample times match and the loop is well posed
};

namespace detail {

/**
 * @brief Generalized plant of the loop: inputs [r; d; n; u], outputs [y; u; r; y + n]
 *
 * d is a load disturbance added to the plant input u, n measurement noise added to y. The last
 * two output blocks are what the controller reads; lft() closes u = K[r; y + n].
 */
template<size_t NX, size_t NU, size_t NY, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr StateSpace<NX, (2 * NY) + (2 * NU), (3 * NY) + NU, T> loop_plant(const StateSpace<NX, NU, NY, T, NW, NV>& P) {
    constexpr size_t d0 = NY;            // columns: r, d, n, u
    constexpr size_t n0 = NY + NU;       //
    constexpr size_t u0 = (2 * NY) + NU; //
    constexpr size_t zu = NY;            // rows: y, u, r, y + n
    constexpr size_t mr = NY + NU;       //
    constexpr size_t my = (2 * NY) + NU; //

    StateSpace<NX, (2 * NY) + (2 * NU), (3 * NY) + NU, T> G{.A = P.A, .Ts = P.Ts};
    if constexpr (NX > 0) {
        G.B.template block<NX, NU>(0, d0) = P.B;
        G.B.template block<NX, NU>(0, u0) = P.B;
        G.C.template block<NY, NX>(0, 0) = P.C;
        G.C.template block<NY, NX>(my, 0) = P.C;
    }
    G.D.template block<NY, NU>(0, d0) = P.D;
    G.D.template block<NY, NU>(0, u0) = P.D;
    G.D.template block<NU, NU>(zu, u0) = Matrix<NU, NU, T>::identity();
    G.D.template block<NY, NY>(mr, 0) = Matrix<NY, NY, T>::identity();
    G.D.template block<NY, NU>(my, d0) = P.D;
    G.D.template block<NY, NY>(my, n0) = Matrix<NY, NY, T>::identity();
    G.D.template block<NY, NU>(my, u0) = P.D;
    return G;
}

/** @brief Rows [r0, r0 + NO) and columns [c0, c0 + NI) of a closed loop, negated when @p negate, plus @p shift on D's diagonal. */
template<size_t NO, size_t NI, size_t NX, size_t NUC, size_t NYC, typename T, size_t NW, size_t NV>
[[nodiscard]] constexpr StateSpace<NX, NI, NO, T>
loop_channel(const StateSpace<NX, NUC, NYC, T, NW, NV>& cl, size_t r0, size_t c0, bool negate = false, T shift = T{0}) {
    StateSpace<NX, NI, NO, T> out{.A = cl.A, .Ts = cl.Ts};
    if constexpr (NX > 0) {
        out.B = cl.B.template block<NX, NI>(0, c0);
        out.C = cl.C.template block<NO, NX>(r0, 0);
    }
    out.D = cl.D.template block<NO, NI>(r0, c0);
    if (negate) {
        out.C = -out.C;
        out.D = -out.D;
    }
    if constexpr (NO == NI) {
        for (size_t i = 0; i < NO; ++i) {
            out.D(i, i) += shift;
        }
    }
    return out;
}

} // namespace detail

/**
 * @brief Gang of seven of a plant and a two-degree-of-freedom controller u = K[r; y]
 *
 * The controller is the map from the reference and the measurement to the control, the way a
 * law's @c control(r, y) reads: its first NY inputs are r, its last NY inputs y. Split as
 * @f$ u = C_r r - C_y y @f$, the loop closes through @f$ C_y @f$ alone, so S, T, PS and CS are
 * those of @f$ C_y @f$, and the reference reaches y, u and e = r − y through @f$ C_r @f$:
 * @f[
 *   y = (I + P C_y)^{-1} P C_r\, r, \qquad u = (I + C_y P)^{-1} C_r\, r .
 * @f]
 * A setpoint weight on P or D changes the reference responses and leaves the four alone.
 *
 * The loop is built as the generalized plant with inputs [r; d; n; u] and outputs
 * [y; u; r; y + n], closed by lft(), so every member is a block of one closed loop on the
 * NXP + NXK states [x_P; x_K]: T is −(n → y), S = I − T, PS is d → y and CS is −(n → u). The
 * plant's and controller's noise channels (G, H) are not carried.
 *
 * Returns @c success = false when the sample times differ or the loop is not well posed (lft()
 * refuses I − D₂₂D_K). An unstable closed loop is returned with @c stable = false; norm_hinf()
 * refuses its members.
 *
 * @note Compare with MATLAB®'s loopsens(P, C) (Robust Control Toolbox), which returns the
 *       one-degree-of-freedom set with input- and output-side forms.
 * @see "Feedback Systems: An Introduction for Scientists and Engineers" (Åström & Murray, 2008),
 *      §12.1 (the gang of four and the gang of six)
 * @see "Multivariable Feedback Control" (Skogestad & Postlethwaite, 2005), §3.2 (input and
 *      output sensitivities)
 * @see gang_of_four(), norm_hinf() for the peaks
 *
 * @param plant      Plant P, NU inputs and NY outputs (continuous or discrete)
 * @param controller Controller K, inputs [r; y] (2·NY), NU outputs, on the plant's sample time
 * @return The seven closed-loop maps, closed-loop stability and success
 */
template<size_t NXP, size_t NU, size_t NY, size_t NWP, size_t NVP, size_t NXK, size_t NUK, size_t NYK, size_t NWK, size_t NVK, typename T>
    requires(NUK == 2 * NY) && (NYK == NU)
[[nodiscard]] constexpr GangOfSevenResult<NXP + NXK, NU, NY, T> gang_of_seven(
    const StateSpace<NXP, NU, NY, T, NWP, NVP>&   plant,
    const StateSpace<NXK, NUK, NYK, T, NWK, NVK>& controller
) {
    GangOfSevenResult<NXP + NXK, NU, NY, T> out{};
    const auto                              cl = lft(detail::loop_plant(plant), controller);
    if (!cl) {
        return out;
    }
    constexpr size_t y0 = 0; // closed-loop rows: y, u
    constexpr size_t u0 = NY;
    constexpr size_t r0 = 0; // closed-loop columns: r, d, n
    constexpr size_t d0 = NY;
    constexpr size_t n0 = NY + NU;

    out.sensitivity = detail::loop_channel<NY, NY>(*cl, y0, n0, false, T{1});
    out.complementary_sensitivity = detail::loop_channel<NY, NY>(*cl, y0, n0, true);
    out.load_sensitivity = detail::loop_channel<NY, NU>(*cl, y0, d0);
    out.noise_sensitivity = detail::loop_channel<NU, NY>(*cl, u0, n0, true);
    out.reference_to_output = detail::loop_channel<NY, NY>(*cl, y0, r0);
    out.reference_to_control = detail::loop_channel<NU, NY>(*cl, u0, r0);
    out.reference_to_error = detail::loop_channel<NY, NY>(*cl, y0, r0, true, T{1});
    out.stable = stability::is_stable(*cl);
    out.success = true;
    return out;
}

/**
 * @brief Gang of four of a plant P and a controller C on the error, u = C(r − y)
 *
 * @f[
 *   S = (I + PC)^{-1},\quad T = PC(I + PC)^{-1},\quad PS = P(I + CP)^{-1},\quad CS = C(I + PC)^{-1}.
 * @f]
 * S and T alone do not show a plant pole the controller cancels: it stays in PS, the response
 * to a load disturbance, and CS is the gain from measurement noise to the actuator. The peaks
 * are norm_hinf() of each member, e.g. @f$ M_s = \lVert S \rVert_\infty @f$.
 *
 * Runs gang_of_seven() with K = [C, −C] on [r; y], so the members are blocks of one closed loop
 * on the NXP + NXC states [x_P; x_C].
 *
 * @code
 * const auto g = analysis::gang_of_four(plant, *design::pid(2.0, 3.0, 0.0).to_ss());
 * const auto Ms = analysis::norm_hinf(g.sensitivity);      // peak sensitivity and its frequency
 * const auto Mn = analysis::norm_hinf(g.noise_sensitivity); // noise → actuator gain
 * @endcode
 *
 * @note Compare with MATLAB®'s loopsens(P, C) (Robust Control Toolbox): So, To, PSi and CSo.
 * @see "Feedback Systems: An Introduction for Scientists and Engineers" (Åström & Murray, 2008),
 *      §12.1 (the gang of four)
 * @see gang_of_seven() for a controller with separate setpoint and feedback paths
 *
 * @param plant      Plant P, NU inputs and NY outputs (continuous or discrete)
 * @param controller Controller C from the error r − y (NY inputs) to u (NU outputs)
 * @return S, T, PS, CS, closed-loop stability and success
 */
template<size_t NXP, size_t NU, size_t NY, size_t NWP, size_t NVP, size_t NXC, size_t NUC, size_t NYC, size_t NWC, size_t NVC, typename T>
    requires(NUC == NY) && (NYC == NU)
[[nodiscard]] constexpr GangOfFourResult<NXP + NXC, NU, NY, T> gang_of_four(
    const StateSpace<NXP, NU, NY, T, NWP, NVP>&   plant,
    const StateSpace<NXC, NUC, NYC, T, NWC, NVC>& controller
) {
    StateSpace<NXC, 2 * NY, NU, T> two_dof{.A = controller.A, .Ts = controller.Ts};
    if constexpr (NXC > 0) {
        two_dof.B.template block<NXC, NY>(0, 0) = controller.B;
        two_dof.B.template block<NXC, NY>(0, NY) = -controller.B;
        two_dof.C = controller.C;
    }
    two_dof.D.template block<NU, NY>(0, 0) = controller.D;
    two_dof.D.template block<NU, NY>(0, NY) = -controller.D;

    const auto g = gang_of_seven(plant, two_dof);
    return GangOfFourResult<NXP + NXC, NU, NY, T>{
        .sensitivity = g.sensitivity,
        .complementary_sensitivity = g.complementary_sensitivity,
        .load_sensitivity = g.load_sensitivity,
        .noise_sensitivity = g.noise_sensitivity,
        .stable = g.stable,
        .success = g.success,
    };
}

} // namespace analysis
} // namespace damp
