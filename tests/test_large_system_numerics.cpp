// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstddef>

#include "damp/backend.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/design/riccati.hpp"
#include "damp/design/stability.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;
using namespace damp::mat;

/**
 * @brief Tests specifically for 4x4 systems to ensure numerical accuracy
 *        for cart-pole, quad-rotor attitude, and other 4-state systems
 */

TEST_SUITE("4x4 Eigenvalue Computation") {
    TEST_CASE("4x4 diagonal matrix - exact eigenvalues") {
        Matrix<4, 4> D = {
            {1.0, 0.0, 0.0, 0.0},
            {0.0, 2.0, 0.0, 0.0},
            {0.0, 0.0, 3.0, 0.0},
            {0.0, 0.0, 0.0, 4.0}
        };

        auto result = compute_eigenvalues(D);
        CHECK(result.converged);

        // Extract eigenvalues and sort
        std::array<double, 4> evs = {};
        for (size_t i = 0; i < 4; ++i) {
            evs[i] = result.values[i].real();
            CHECK(std::abs(result.values[i].imag()) < 1e-10);
        }
        std::ranges::sort(evs);

        CHECK(evs[0] == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(evs[1] == doctest::Approx(2.0).epsilon(1e-6));
        CHECK(evs[2] == doctest::Approx(3.0).epsilon(1e-6));
        CHECK(evs[3] == doctest::Approx(4.0).epsilon(1e-6));
    }

    TEST_CASE("4x4 identity matrix") {
        auto I = Matrix<4, 4, double>::identity();

        auto result = compute_eigenvalues(I);
        CHECK(result.converged);

        for (size_t i = 0; i < 4; ++i) {
            CHECK(result.values[i].real() == doctest::Approx(1.0).epsilon(1e-10));
            CHECK(std::abs(result.values[i].imag()) < 1e-10);
        }
    }

    TEST_CASE("4x4 symmetric positive definite matrix") {
        // SPD matrix with known eigenvalues
        Matrix<4, 4> A = {
            {4.0, 1.0, 0.5, 0.2},
            {1.0, 3.0, 0.5, 0.3},
            {0.5, 0.5, 2.0, 0.1},
            {0.2, 0.3, 0.1, 1.0}
        };

        auto result = compute_eigenvalues(A);
        CHECK(result.converged);

        // All eigenvalues should be positive (SPD)
        for (size_t i = 0; i < 4; ++i) {
            CHECK(result.values[i].real() > 0.0);
            CHECK(std::abs(result.values[i].imag()) < 1e-6);
        }

        // Sum of eigenvalues equals trace
        double trace = A(0, 0) + A(1, 1) + A(2, 2) + A(3, 3);
        double sum_ev = 0.0;
        for (size_t i = 0; i < 4; ++i) {
            sum_ev += result.values[i].real();
        }
        CHECK(sum_ev == doctest::Approx(trace).epsilon(1e-6));
    }

    TEST_CASE("4x4 matrix with complex eigenvalues") {
        // Rotation-like matrix with complex eigenvalues
        Matrix<4, 4> A = {
            {0.0, -1.0, 0.0, 0.0},
            {1.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, -2.0},
            {0.0, 0.0, 2.0, 0.0}
        };

        auto result = compute_eigenvalues(A);
        CHECK(result.converged);

        // Should have complex conjugate pairs: ±i and ±2i
        // Count pairs with imaginary parts
        int complex_count = 0;
        for (size_t i = 0; i < 4; ++i) {
            if (std::abs(result.values[i].imag()) > 0.1) {
                complex_count++;
            }
        }
        CHECK(complex_count == 4);
    }

    TEST_CASE("4x4 cart-pole system matrix (typical control problem)") {
        // Linearized cart-pole system at upright equilibrium
        // State: [x, x_dot, theta, theta_dot]
        constexpr double g = 9.81;
        constexpr double m = 0.1; // pole mass
        constexpr double M = 1.0; // cart mass
        constexpr double L = 0.5; // pole length

        // Continuous A matrix for cart-pole
        Matrix<4, 4> A = {
            {0.0, 1.0, 0.0, 0.0},
            {0.0, 0.0, -m * g / M, 0.0},
            {0.0, 0.0, 0.0, 1.0},
            {0.0, 0.0, (M + m) * g / (M * L), 0.0}
        };

        auto result = compute_eigenvalues(A);
        CHECK(result.converged);

        // Cart-pole at upright has one positive real eigenvalue (unstable)
        // and one negative real eigenvalue, plus two zeros
        bool has_positive = false;
        bool has_negative = false;
        for (size_t i = 0; i < 4; ++i) {
            if (result.values[i].real() > 0.1) {
                has_positive = true;
            }
            if (result.values[i].real() < -0.1) {
                has_negative = true;
            }
        }
        CHECK(has_positive);
        CHECK(has_negative);
    }

    TEST_CASE("4x4 consteval eigenvalue computation") {
        constexpr Matrix<4, 4> D = {
            {1.0, 0.0, 0.0, 0.0},
            {0.0, 2.0, 0.0, 0.0},
            {0.0, 0.0, 3.0, 0.0},
            {0.0, 0.0, 0.0, 4.0}
        };

        constexpr auto result = compute_eigenvalues(D);
        static_assert(result.converged, "4x4 eigenvalue computation should converge");

        // Verify at compile time that eigenvalues sum to trace (1+2+3+4=10)
        constexpr double sum_ev = result.values[0].real() + result.values[1].real() + result.values[2].real() + result.values[3].real();
        static_assert(sum_ev > 9.9 && sum_ev < 10.1, "Sum of eigenvalues should equal trace");
    }
}

TEST_SUITE("4x4 DARE and LQR") {
    TEST_CASE("4-state discrete LQR (stable test system)") {
        // A well-conditioned stable 4-state system
        Matrix<4, 4> A_d = {
            {0.9, 0.05, 0.0, 0.0},
            {0.0, 0.9, 0.05, 0.0},
            {0.0, 0.0, 0.9, 0.05},
            {0.0, 0.0, 0.0, 0.9}
        };
        Matrix<4, 1> B_d = {
            {0.1},
            {0.05},
            {0.02},
            {0.01}
        };

        Matrix<4, 4> Q = Matrix<4, 4>::identity();
        Matrix<1, 1> R = {{1.0}};

        auto result = design::discrete_lqr(A_d, B_d, Q, R);
        CHECK(result.success);

        if (result.success) {
            // Verify closed-loop stability
            Matrix<4, 4> A_cl = A_d - B_d * result.K;
            CHECK(stability::is_stable_discrete(A_cl));
        }
    }

    TEST_CASE("4-state unstable system stabilized by LQR") {
        // An open-loop unstable 4-state system (one eigenvalue > 1)
        Matrix<4, 4> A_d = {
            {1.05, 0.05, 0.0, 0.0}, // Marginally unstable mode
            {0.0, 0.9, 0.05, 0.0},
            {0.0, 0.0, 0.9, 0.05},
            {0.0, 0.0, 0.0, 0.9}
        };
        Matrix<4, 1> B_d = {
            {0.2}, // Strong control authority on first state
            {0.05},
            {0.02},
            {0.01}
        };

        Matrix<4, 4> Q = Matrix<4, 4>::identity();
        Matrix<1, 1> R = {{0.1}}; // Lower control cost to enable stabilization

        auto result = design::discrete_lqr(A_d, B_d, Q, R);
        CHECK(result.success);

        if (result.success) {
            // Verify closed-loop stability
            Matrix<4, 4> A_cl = A_d - B_d * result.K;
            CHECK(stability::is_stable_discrete(A_cl));
        }
    }

    TEST_CASE("4-state DARE solution quality") {
        // Well-conditioned 4x4 system
        Matrix<4, 4> A = {
            {0.9, 0.1, 0.0, 0.0},
            {0.0, 0.9, 0.1, 0.0},
            {0.0, 0.0, 0.9, 0.1},
            {0.0, 0.0, 0.0, 0.9}
        };
        Matrix<4, 1> B = {
            {0.1},
            {0.0},
            {0.0},
            {0.0}
        };
        Matrix<4, 4> Q = Matrix<4, 4>::identity();
        Matrix<1, 1> R = {{1.0}};

        auto P_opt = dare(A, B, Q, R);
        REQUIRE(P_opt.has_value());

        auto P = P_opt.value();

        // Verify P is symmetric
        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = i + 1; j < 4; ++j) {
                CHECK(P(i, j) == doctest::Approx(P(j, i)).epsilon(1e-8));
            }
        }

        // Verify P is positive definite (all eigenvalues > 0)
        auto P_eigen = compute_eigenvalues(P);
        CHECK(P_eigen.converged);
        for (size_t i = 0; i < 4; ++i) {
            CHECK(P_eigen.values[i].real() > 0.0);
        }

        // Verify DARE equation: P = A'PA - A'PB(R + B'PB)^{-1}B'PA + Q
        auto S = R + B.transpose() * P * B;
        auto S_inv_opt = S.inverse();
        REQUIRE(S_inv_opt.has_value());
        auto K = S_inv_opt.value() * B.transpose() * P * A;
        auto P_check = Q + A.transpose() * P * A - A.transpose() * P * B * K;

        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                CHECK(P(i, j) == doctest::Approx(P_check(i, j)).epsilon(1e-6));
            }
        }
    }

    TEST_CASE("4-state compile-time LQR design") {
        constexpr Matrix<4, 4> A = {
            {0.9, 0.1, 0.0, 0.0},
            {0.0, 0.9, 0.1, 0.0},
            {0.0, 0.0, 0.9, 0.1},
            {0.0, 0.0, 0.0, 0.9}
        };
        constexpr Matrix<4, 1> B = {
            {0.1},
            {0.0},
            {0.0},
            {0.0}
        };
        constexpr Matrix<4, 4> Q = Matrix<4, 4>::identity();
        constexpr Matrix<1, 1> R = {{1.0}};

        constexpr auto result = design::discrete_lqr(A, B, Q, R);

        static_assert(result.success, "Compile-time 4x4 LQR should succeed");
        static_assert(result.K(0, 0) != 0.0, "K should have non-zero elements");

        // Runtime verification
        CHECK(result.success);
        CHECK(result.K(0, 0) != 0.0);
    }

    TEST_CASE("4-state 2-input LQR") {
        // 4 states, 2 inputs
        Matrix<4, 4> A = {
            {0.95, 0.05, 0.0, 0.0},
            {0.0, 0.95, 0.0, 0.0},
            {0.0, 0.0, 0.95, 0.05},
            {0.0, 0.0, 0.0, 0.95}
        };
        Matrix<4, 2> B = {
            {0.1, 0.0},
            {0.05, 0.0},
            {0.0, 0.1},
            {0.0, 0.05}
        };
        Matrix<4, 4> Q = Matrix<4, 4>::identity();
        Matrix<2, 2> R = Matrix<2, 2>::identity();

        auto result = design::discrete_lqr(A, B, Q, R);
        CHECK(result.success);

        // Verify closed-loop stability
        Matrix<4, 4> A_cl = A - B * result.K;
        CHECK(stability::is_stable_discrete(A_cl));
    }
}

TEST_SUITE("4x4 Kalman Filter Design") {
    TEST_CASE("4-state Kalman filter") {
        // 4 states, 1 input, 2 outputs
        StateSpace<4, 1, 2, double, 4, 2> sys;
        sys.A = {
            {0.95, 0.05, 0.0, 0.0},
            {0.0, 0.95, 0.0, 0.0},
            {0.0, 0.0, 0.95, 0.05},
            {0.0, 0.0, 0.0, 0.95}
        };
        sys.B = {{0.1}, {0.0}, {0.0}, {0.0}};
        sys.C = {
            {1.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 1.0, 0.0}
        };
        sys.D = Matrix<2, 1>::zeros();
        sys.G = Matrix<4, 4>::identity();
        sys.H = Matrix<2, 2>::identity();
        sys.Ts = 0.01;

        Matrix<4, 4> Q = Matrix<4, 4>::identity() * 0.01;
        Matrix<2, 2> R = Matrix<2, 2>::identity() * 0.1;

        auto result = design::kalman(sys, Q, R);
        CHECK(result.success);

        // Verify L is non-zero
        bool L_nonzero = false;
        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < 2; ++j) {
                if (result.L(i, j) != 0.0) {
                    L_nonzero = true;
                    break;
                }
            }
        }
        CHECK(L_nonzero);
    }
}

TEST_SUITE("4x4 Matrix Inverse") {
    TEST_CASE("4x4 identity inverse") {
        auto I = Matrix<4, 4, double>::identity();
        auto I_inv_opt = I.inverse();

        REQUIRE(I_inv_opt.has_value());
        auto I_inv = I_inv_opt.value();

        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                double expected = (i == j) ? 1.0 : 0.0;
                CHECK(I_inv(i, j) == doctest::Approx(expected).epsilon(1e-12));
            }
        }
    }

    TEST_CASE("4x4 diagonal inverse") {
        Matrix<4, 4> D = {
            {2.0, 0.0, 0.0, 0.0},
            {0.0, 3.0, 0.0, 0.0},
            {0.0, 0.0, 4.0, 0.0},
            {0.0, 0.0, 0.0, 5.0}
        };

        auto D_inv_opt = D.inverse();
        REQUIRE(D_inv_opt.has_value());
        auto D_inv = D_inv_opt.value();

        CHECK(D_inv(0, 0) == doctest::Approx(0.5).epsilon(1e-12));
        CHECK(D_inv(1, 1) == doctest::Approx(1.0 / 3.0).epsilon(1e-12));
        CHECK(D_inv(2, 2) == doctest::Approx(0.25).epsilon(1e-12));
        CHECK(D_inv(3, 3) == doctest::Approx(0.2).epsilon(1e-12));
    }

    TEST_CASE("4x4 general inverse verification") {
        Matrix<4, 4> A = {
            {1.0, 2.0, 3.0, 4.0},
            {5.0, 6.0, 7.0, 8.0},
            {9.0, 10.0, 12.0, 12.0},
            {13.0, 14.0, 15.0, 17.0}
        };

        auto A_inv_opt = A.inverse();
        REQUIRE(A_inv_opt.has_value());
        const auto& A_inv = A_inv_opt.value();

        // Verify A * A^{-1} = I
        auto I_check = A * A_inv;
        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                double expected = (i == j) ? 1.0 : 0.0;
                CHECK(I_check(i, j) == doctest::Approx(expected).epsilon(1e-8));
            }
        }
    }

    TEST_CASE("4x4 singular matrix detection") {
        // Matrix with linearly dependent rows
        Matrix<4, 4> S = {
            {1.0, 2.0, 3.0, 4.0},
            {2.0, 4.0, 6.0, 8.0}, // 2x row 0
            {5.0, 6.0, 7.0, 8.0},
            {9.0, 10.0, 11.0, 12.0}
        };

        auto S_inv_opt = S.inverse();
        CHECK_FALSE(S_inv_opt.has_value());
    }
}

TEST_SUITE("4x4 Matrix Exponential") {
    TEST_CASE("4x4 zero matrix exponential is identity") {
        auto Z = Matrix<4, 4, double>::zeros();
        auto exp_Z = mat::expm(Z);

        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                double expected = (i == j) ? 1.0 : 0.0;
                CHECK(exp_Z(i, j) == doctest::Approx(expected).epsilon(1e-10));
            }
        }
    }

    TEST_CASE("4x4 diagonal matrix exponential") {
        Matrix<4, 4> D = {
            {1.0, 0.0, 0.0, 0.0},
            {0.0, 2.0, 0.0, 0.0},
            {0.0, 0.0, -1.0, 0.0},
            {0.0, 0.0, 0.0, 0.0}
        };

        auto exp_D = mat::expm(D);

        // Improved Padé approximation should achieve ~1e-6 accuracy
        CHECK(exp_D(0, 0) == doctest::Approx(std::exp(1.0)).epsilon(1e-6));
        CHECK(exp_D(1, 1) == doctest::Approx(std::exp(2.0)).epsilon(1e-6));
        CHECK(exp_D(2, 2) == doctest::Approx(std::exp(-1.0)).epsilon(1e-6));
        CHECK(exp_D(3, 3) == doctest::Approx(1.0).epsilon(1e-10));
    }

    TEST_CASE("4x4 nilpotent matrix exponential") {
        // Upper triangular nilpotent (only superdiagonals)
        Matrix<4, 4> N = {
            {0.0, 1.0, 0.0, 0.0},
            {0.0, 0.0, 1.0, 0.0},
            {0.0, 0.0, 0.0, 1.0},
            {0.0, 0.0, 0.0, 0.0}
        };

        auto exp_N = mat::expm(N);

        // exp(N) = I + N + N²/2! + N³/3!
        // N² has ones on 2nd superdiagonal
        // N³ has one in corner (0,3)
        CHECK(exp_N(0, 0) == doctest::Approx(1.0).epsilon(1e-10));
        CHECK(exp_N(0, 1) == doctest::Approx(1.0).epsilon(1e-10));
        CHECK(exp_N(0, 2) == doctest::Approx(0.5).epsilon(1e-10));
        CHECK(exp_N(0, 3) == doctest::Approx(1.0 / 6.0).epsilon(1e-10));
    }
}

TEST_CASE("Matrix exponential (exp) - compile-time and runtime verification") {

    // Test 1: exp(0) = I (identity matrix)
    constexpr Matrix<2, 2, double> A_zero = Matrix<2, 2, double>::zeros();
    constexpr auto                 exp_zero = expm(A_zero);
    CHECK(exp_zero(0, 0) == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(exp_zero(0, 1) == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(exp_zero(1, 0) == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(exp_zero(1, 1) == doctest::Approx(1.0).epsilon(1e-12));
    static_assert(exp_zero == Matrix<2, 2, double>::identity(), "exp(0) should be identity");

    // Test 2: exp(scalar * I) = scalar * I (for small scalar)
    constexpr double scalar = 0.1;
    constexpr Matrix A_scalar = Matrix<2, 2, double>::identity() * scalar;
    constexpr auto   exp_scalar = expm(A_scalar);
    static_assert(exp_scalar(0, 0) - damp::exp(scalar) < 1e-12, "exp(scalar*I) diagonal should match scalar exp");

    // Test 3: Runtime test with 4x4 matrix (golden reference from scipy.linalg.expm)
    // Example matrix: [[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0]] (cyclic permutation)
    Matrix<4, 4, double> A_4x4 = {
        {0.0, 1.0, 0.0, 0.0},
        {0.0, 0.0, 1.0, 0.0},
        {0.0, 0.0, 0.0, 1.0},
        {1.0, 0.0, 0.0, 0.0},
    };

    auto exp_A_4x4 = expm(A_4x4);

    // Golden reference (computed via scipy.linalg.expm in Python)
    Matrix<4, 4, double> expected_exp_A_4x4 = {
        {1.04169147, 1.00833609, 0.50138916, 0.1668651},
        {0.1668651, 1.04169147, 1.00833609, 0.50138916},
        {0.50138916, 0.1668651, 1.04169147, 1.00833609},
        {1.00833609, 0.50138916, 0.1668651, 1.04169147},
    };

    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            CHECK(exp_A_4x4(i, j) == doctest::Approx(expected_exp_A_4x4(i, j)).epsilon(1e-8));
        }
    }

    // Test 4: Verify exp(A) * exp(-A) ≈ I for small A
    Matrix<2, 2, double> A_small = {
        {0.1, 0.05},
        {-0.05, 0.1}
    };

    auto exp_A_small = expm(A_small);
    auto exp_neg_A_small = expm(A_small * -1.0);
    auto product = exp_A_small * exp_neg_A_small;
    auto identity = Matrix<2, 2, double>::identity();
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            CHECK(product(i, j) == doctest::Approx(identity(i, j)).epsilon(1e-8));
        }
    }
}

TEST_SUITE("4x4 Stability Analysis") {
    TEST_CASE("4x4 discrete stability check - stable system") {
        // All eigenvalues inside unit circle
        Matrix<4, 4> A = {
            {0.5, 0.1, 0.0, 0.0},
            {0.0, 0.6, 0.1, 0.0},
            {0.0, 0.0, 0.7, 0.1},
            {0.0, 0.0, 0.0, 0.8},
        };

        CHECK(stability::is_stable_discrete(A));
    }

    TEST_CASE("4x4 discrete stability check - unstable system") {
        // At least one eigenvalue outside unit circle
        Matrix<4, 4> A = {
            {1.1, 0.0, 0.0, 0.0}, // eigenvalue 1.1 > 1
            {0.0, 0.5, 0.0, 0.0},
            {0.0, 0.0, 0.5, 0.0},
            {0.0, 0.0, 0.0, 0.5},
        };

        CHECK_FALSE(stability::is_stable_discrete(A));
    }

    TEST_CASE("4x4 stability margin") {
        Matrix<4, 4> A = {
            {0.5, 0.0, 0.0, 0.0},
            {0.0, 0.6, 0.0, 0.0},
            {0.0, 0.0, 0.7, 0.0},
            {0.0, 0.0, 0.0, 0.8},
        };

        double margin = stability::stability_margin_discrete(A);
        CHECK(margin > 0.0);                                 // Stable
        CHECK(margin == doctest::Approx(0.2).epsilon(0.01)); // 1 - 0.8 = 0.2
    }
}

TEST_SUITE("Larger Systems (5x5+)") {
    TEST_CASE("5x5 DARE solution works via QR iteration") {
        // 5-state system - DARE should work because it uses QR iteration
        Matrix<5, 5> A = {
            {0.9, 0.05, 0.0, 0.0, 0.0},
            {0.0, 0.9, 0.05, 0.0, 0.0},
            {0.0, 0.0, 0.9, 0.05, 0.0},
            {0.0, 0.0, 0.0, 0.9, 0.05},
            {0.0, 0.0, 0.0, 0.0, 0.9},
        };
        Matrix<5, 1> B = {{0.1}, {0.05}, {0.02}, {0.01}, {0.005}};
        Matrix<5, 5> Q = Matrix<5, 5>::identity();
        Matrix<1, 1> R = {{1.0}};

        auto P_opt = dare(A, B, Q, R);
        CHECK(P_opt.has_value());

        if (P_opt.has_value()) {
            auto P = P_opt.value();
            // Verify P is symmetric
            for (size_t i = 0; i < 5; ++i) {
                for (size_t j = i + 1; j < 5; ++j) {
                    CHECK(P(i, j) == doctest::Approx(P(j, i)).epsilon(1e-8));
                }
            }
        }
    }

    TEST_CASE("6x6 matrix inverse and operations") {
        /* This test is particularly good, as the matrix is non-SPD / non-symmetric, which checks the Cholesky fallback path */

        // Verify 6x6 matrix operations work
        Matrix<6, 6> A = Matrix<6, 6>::identity();
        A(0, 1) = 0.1;
        A(1, 2) = 0.1;
        A(2, 3) = 0.1;
        A(3, 4) = 0.1;
        A(4, 5) = 0.1;

        for (size_t i = 0; i < 6; ++i) {
            for (size_t j = 0; j < 6; ++j) {
            }
        }

        auto A_inv_opt = A.inverse();
        REQUIRE(A_inv_opt.has_value());

        const auto& A_inv = A_inv_opt.value();
        for (size_t i = 0; i < 6; ++i) {
            for (size_t j = 0; j < 6; ++j) {
            }
        }

        auto I_check = A * A_inv;
        for (size_t i = 0; i < 6; ++i) {
            for (size_t j = 0; j < 6; ++j) {
                double expected = (i == j) ? 1.0 : 0.0;
                CHECK(I_check(i, j) == doctest::Approx(expected).epsilon(1e-8));
            }
        }
    }
}