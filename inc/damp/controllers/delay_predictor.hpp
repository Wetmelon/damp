// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

/**
 * @file delay_predictor.hpp
 * @brief Input delay by prediction: the delay-free law on the state the command will meet
 */

#include <cstddef>

#include "damp/backend.hpp"
#include "damp/matrix/colvec.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

namespace damp {

/**
 * @brief Predictor form of an input delay of D whole samples.
 *
 * A command sent now reaches the plant D periods later, so the law that decides it should
 * see the state it will meet then, not the state now. With the D commands already in
 * flight that state is known exactly:
 * @f[
 *   x_{k+D} = A^D x_k + \sum_{j=0}^{D-1} A^{D-1-j} B\, u_{k-D+j}.
 * @f]
 * The optimal law on the delay-augmented model, with no cost on its register states, is the
 * delay-free law applied to this prediction, and the augmented Kalman filter is the
 * delay-free filter driven by the input the plant actually received, since the register
 * states carry no uncertainty. So the design stays on the plant (plus one register state for
 * a fractional remainder, see split_input_delay), the D commands live in this buffer, and the
 * A matrix keeps the plant's size: O(D·NX²) per tick and D·NU numbers of memory, where the
 * register form (discretize with InputDelay<D>) costs (NX + D·NU)² in both. The two agree
 * to roundoff (tests/test_delay_predictor.cpp). For MPC the register form spends its first D
 * prediction steps on outputs no decision can change, so a predictor-form horizon of NP
 * matches a register-form horizon of NP + D; with a preview, the plan sees the reference D
 * steps further on (OffsetFreeMPC shifts it).
 *
 * Buffer at tick k, after the previous tick's push: u_{k−D−1} … u_{k−1}.
 * - applied_last(): u_{k−D−1}, what the plant received over the last period (the estimator's predict)
 * - applied_now():  u_{k−D}, what it receives over this one (the feedthrough in the update)
 * - predict(x̂_k):  x_{k+D} from u_{k−D} … u_{k−1}
 * - push(u_k):      at the end of the tick, after the estimator's predict
 *
 * Give it the model the estimator runs on: for an offset-free design that is the
 * disturbance-augmented model, so the estimated disturbance acts during the delay too.
 * OFLQG and OffsetFreeMPC do all of this internally: give the design call the delay in
 * periods (discrete_oflqg<θ / Ts>, design::mpc<NP, NC, θ / Ts>) and the runtime deduces D
 * from the result.
 *
 * @see "Linear systems with delayed controls: A reduction" (Artstein, 1982), IEEE TAC 27(4)
 * @see "Digital Control of Dynamic Systems" (Franklin, Powell, Workman, 1998), §4.3.3
 * @note Compare with MATLAB®'s absorbDelay, the register form this is equivalent to.
 */
template<size_t NX, size_t NU, size_t D, typename T = float>
struct InputDelayPredictor {
    static_assert(D >= 1, "a zero delay needs no predictor");

    Matrix<NX, NX, T>                 A{};      ///< State map of the model the estimate lives on
    Matrix<NX, NU, T>                 B{};      ///< Input map
    damp::array<ColVec<NU, T>, D + 1> buffer{}; ///< the last D+1 commands, oldest at @ref head
    size_t                            head{0};  ///< index of the oldest entry

    constexpr InputDelayPredictor() = default;

    /// From the estimator's model (any scalar; converted with .as<T>()).
    template<size_t NY, typename U, size_t NW, size_t NV>
    constexpr explicit InputDelayPredictor(const StateSpace<NX, NU, NY, U, NW, NV>& model)
        : A(model.A.template as<T>()), B(model.B.template as<T>()) {}

    /// u_{k−D−1}: the command the plant received over the last period.
    [[nodiscard]] constexpr const ColVec<NU, T>& applied_last() const { return buffer[head]; }

    /// u_{k−D}: the command the plant receives over this period.
    [[nodiscard]] constexpr const ColVec<NU, T>& applied_now() const { return buffer[(head + 1) % (D + 1)]; }

    /// x_{k+D}: the state the command decided now will meet.
    [[nodiscard]] constexpr ColVec<NX, T> predict(ColVec<NX, T> x) const {
        for (size_t j = 0; j < D; ++j) {
            x = (A * x) + (B * buffer[(head + 1 + j) % (D + 1)]);
        }
        return x;
    }

    /// The same prediction in deviation coordinates: @p x relative to an equilibrium held by
    /// @p u_offset, each buffered command taken relative to it.
    [[nodiscard]] constexpr ColVec<NX, T> predict(ColVec<NX, T> x, const ColVec<NU, T>& u_offset) const {
        for (size_t j = 0; j < D; ++j) {
            x = (A * x) + (B * ColVec<NU, T>(buffer[(head + 1 + j) % (D + 1)] - u_offset));
        }
        return x;
    }

    /// End of the tick: the command just sent replaces the oldest entry.
    constexpr void push(const ColVec<NU, T>& u) {
        buffer[head] = u;
        head = (head + 1) % (D + 1);
    }

    /// Forget the commands in flight (zero).
    constexpr void reset() {
        for (auto& u : buffer) {
            u = ColVec<NU, T>{};
        }
        head = 0;
    }
};

} // namespace damp
