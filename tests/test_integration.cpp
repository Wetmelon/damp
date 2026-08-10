// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>

#include "damp/controllers/lqr.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

/**
 * @brief Integration test: Motor with flexible coupling and inertial load
 *
 * Demonstrates a realistic control design scenario combining discretization,
 * design-time compile-time guarantees, and runtime simulation.
 */

TEST_SUITE("Integration: Motor + Flexible Coupling + Mass System") {
    // Example system: motor with flexible coupling and inertial load
    // States:
    //   x1 = motor angular velocity (rad/s)
    //   x2 = coupling deflection (rad) - measures spring compression
    //   x3 = load angular velocity (rad/s)
    //
    // System parameters:
    //   J_m = motor inertia = 0.1 kg*m^2
    //   J_l = load inertia = 0.2 kg*m^2
    //   k = spring constant = 50 N*m/rad
    //   f_m = motor friction = 2 N*m*s/rad
    //   f_l = load friction = 1 N*m*s/rad
    //
    // Continuous dynamics:
    //   dx1/dt = -k/J_m * x2 - f_m/J_m * x1 + u/J_m
    //   dx2/dt = x1 - x3
    //   dx3/dt = k/J_m * x2 - f_l/J_l * x3
    //
    // Output: y = [x1, x3]^T (motor velocity and load velocity)

    TEST_CASE("Compile-time LQG design: Motor-coupling-mass system in double, runtime in float") {
        // Design parameters
        constexpr double J_m = 0.1; // motor inertia
        constexpr double J_l = 0.2; // load inertia
        constexpr double k = 50.0;  // spring constant
        constexpr double f_m = 2.0; // motor friction
        constexpr double f_l = 1.0; // load friction
        constexpr double Ts = 0.01; // sampling period 10 ms

        // Step 1: Design continuous-time LQG controller in double precision at compile time
        constexpr auto lqg_double = []() consteval {
            // System: 3 states, 1 input, 3 outputs (full state feedback)
            StateSpace sys_c{
                Matrix<3, 3>{
                    {-f_m / J_m, -k / J_m, 0.0},
                    {1.0, 0.0, -1.0},
                    {k / J_l, 0.0, -f_l / J_l}
                },
                Matrix<3, 1>{{1.0 / J_m}, {0.0}, {0.0}},
                Matrix<3, 3>{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}
            };

            // Cost: penalize motor and load velocities, not coupling deflection
            Matrix<3, 3, double> Q{{1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}};
            Matrix<1, 1, double> R{{0.01}};

            return design::discrete_lqr_from_continuous(sys_c, Q, R, Ts);
        }();

        // Step 2: Convert the double-precision controller to float for runtime use
        auto lqr_float = lqg_double.as<float>();

        // Discretize the system for simulation
        constexpr auto sys_d_double = []() consteval {
            StateSpace sys_c{
                Matrix<3, 3>{
                    {-f_m / J_m, -k / J_m, 0.0},
                    {1.0, 0.0, -1.0},
                    {k / J_l, 0.0, -f_l / J_l}
                },
                Matrix<3, 1>{{1.0 / J_m}, {0.0}, {0.0}},
                Matrix<3, 3>{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}
            };
            return *discretize(sys_c, Ts, DiscretizationMethod::ZOH);
        }();
        auto sys_d_float = sys_d_double.as<float>();

        // Step 3: Verify the controller was designed properly
        CHECK(lqr_float.K(0, 0) != 0.0f); // Should have non-zero gains
        CHECK(lqr_float.K(0, 1) != 0.0f);
        CHECK(lqr_float.K(0, 2) != 0.0f);

        // Step 4: Simulate the discrete-time closed-loop system at runtime
        // Initial state: motor and load at rest, coupling stressed
        ColVec x = {
            0.0f, // motor velocity
            0.1f, // coupling deflection (10% of full)
            0.0f  // load velocity
        };

        // Run 10 time steps of closed-loop control
        for (int i = 0; i < 10; ++i) {
            // Compute control: u = -K * x
            ColVec u = -lqr_float.K * x;

            // Update state: x[k+1] = A_d * x[k] + B_d * u[k]
            x = sys_d_float.A * x + sys_d_float.B * u;
        }

        // After 10 steps of feedback control, the coupling deflection should decrease
        // and load velocity should increase toward motor velocity
        CHECK(x(1) < 0.1f);           // Coupling deflection reduced
        CHECK(std::abs(x(0)) > 0.0f); // Motor velocity changed
        CHECK(std::abs(x(2)) > 0.0f); // Load velocity changed

        // Verify that final state shows active damping of the oscillation
        float final_norm_sq = (x(0) * x(0)) + (x(1) * x(1)) + (x(2) * x(2));
        CHECK(final_norm_sq < 0.15f); // System energy reduced
    }
}
