// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <complex>
#include <cstddef>

#include "damp/backend.hpp"
#include "damp/controllers/lqg.hpp"
#include "damp/controllers/lqgi.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/design/stability.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/filters/iir_design.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/matlab.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"
#include "damp/systems/transfer_function.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

TEST_SUITE("matlab") {

    TEST_CASE("Blkdiag") {
        Matrix<2, 2> mat1 = {{1, 2}, {3, 4}};
        Matrix<3, 3> mat2 = {{5, 6, 7}, {8, 9, 10}, {11, 12, 13}};
        Matrix<1, 1> mat3 = {{14}};

        auto block_diag = matlab::blkdiag(mat1, mat2, mat3);

        CHECK(block_diag.rows() == 6);
        CHECK(block_diag.cols() == 6);

        // Check first block
        CHECK(block_diag(0, 0) == 1);
        CHECK(block_diag(0, 1) == 2);
        CHECK(block_diag(1, 0) == 3);
        CHECK(block_diag(1, 1) == 4);

        // Check second block
        CHECK(block_diag(2, 2) == 5);
        CHECK(block_diag(2, 3) == 6);
        CHECK(block_diag(2, 4) == 7);
        CHECK(block_diag(3, 2) == 8);
        CHECK(block_diag(3, 3) == 9);
        CHECK(block_diag(3, 4) == 10);
        CHECK(block_diag(4, 2) == 11);
        CHECK(block_diag(4, 3) == 12);
        CHECK(block_diag(4, 4) == 13);

        // Check third block
        CHECK(block_diag(5, 5) == 14);

        // Check off-diagonal zeros
        for (size_t r = 0; r < block_diag.rows(); ++r) {
            for (size_t c = 0; c < block_diag.cols(); ++c) {
                if ((r < 2 && c < 2) || (r >= 2 && r < 5 && c >= 2 && c < 5) || (r == 5 && c == 5)) {
                    continue; // Skip diagonal blocks
                }
                CHECK(block_diag(r, c) == 0);
            }
        }
    }

    TEST_CASE("Place - Pole placement for single input") {
        // Simple 2x2 system
        Matrix<2, 2> A = {{0, 1}, {-1, -2}};
        Matrix<2, 1> B = {{0}, {1}};

        auto poles = std::array<damp::complex<double>, 2>{
            damp::complex<double>{-3.0, 0.0},
            damp::complex<double>{-4.0, 0.0},
        }; // Desired closed-loop poles

        auto K_opt = matlab::acker(A, B, poles);
        REQUIRE(K_opt.has_value());
        const auto& K = *K_opt;

        CHECK(K.rows() == 1);
        CHECK(K.cols() == 2);

        // Check that the closed-loop poles are approximately the desired ones
        // Closed-loop A_cl = A - B*K
        Matrix<2, 2> A_cl = A - B * K;

        // Compute characteristic polynomial of A_cl
        // For 2x2, det(sI - A_cl) = s^2 - trace*s + det
        double trace = A_cl(0, 0) + A_cl(1, 1);
        double det = (A_cl(0, 0) * A_cl(1, 1)) - (A_cl(0, 1) * A_cl(1, 0));

        // Roots of s^2 + trace*s + det = 0 should be -3 and -4
        // Sum of roots = -trace = 3+4=7
        // Product = det = 12
        CHECK(doctest::Approx(trace).epsilon(1e-6) == -7.0);
        CHECK(doctest::Approx(det).epsilon(1e-6) == 12.0);
    }

    TEST_CASE("Place - Pole placement with std::complex poles") {
        // Simple 2x2 system
        Matrix<2, 2> A = {{0, 1}, {-1, -2}};
        Matrix<2, 1> B = {{0}, {1}};

        auto poles = std::array{std::complex<double>(-3, 0), std::complex<double>(-4, 0)}; // Desired closed-loop poles

        auto K_opt = matlab::acker(A, B, poles);
        REQUIRE(K_opt.has_value());
        const auto& K = *K_opt;

        CHECK(K.rows() == 1);
        CHECK(K.cols() == 2);

        // Check that the closed-loop poles are approximately the desired ones
        Matrix<2, 2> A_cl = A - B * K;
        double       trace = A_cl(0, 0) + A_cl(1, 1);
        double       det = (A_cl(0, 0) * A_cl(1, 1)) - (A_cl(0, 1) * A_cl(1, 0));
        CHECK(doctest::Approx(trace).epsilon(1e-6) == -7.0);
        CHECK(doctest::Approx(det).epsilon(1e-6) == 12.0);
    }

    TEST_CASE("Place - Pole placement with damp::complex poles") {
        // Simple 2x2 system
        Matrix<2, 2> A = {{0, 1}, {-1, -2}};
        Matrix<2, 1> B = {{0}, {1}};

        auto poles = std::array{damp::complex<double>(-3, 0), damp::complex<double>(-4, 0)}; // Desired closed-loop poles

        auto K_opt = matlab::acker(A, B, poles);
        REQUIRE(K_opt.has_value());
        const auto& K = *K_opt;

        CHECK(K.rows() == 1);
        CHECK(K.cols() == 2);

        // Check that the closed-loop poles are approximately the desired ones
        Matrix<2, 2> A_cl = A - B * K;
        double       trace = A_cl(0, 0) + A_cl(1, 1);
        double       det = (A_cl(0, 0) * A_cl(1, 1)) - (A_cl(0, 1) * A_cl(1, 0));
        CHECK(doctest::Approx(trace).epsilon(1e-6) == -7.0);
        CHECK(doctest::Approx(det).epsilon(1e-6) == 12.0);
    }

    TEST_CASE("PID Tune") {
        // Simple second-order plant: G(s) = 1/(s^2 + 0.5s + 1)
        Matrix<2, 2, double> A = {{0, 1}, {-1, -0.5}};
        Matrix<2, 1, double> B = {{0}, {1}};
        Matrix<1, 2, double> C = {{1, 0}};
        Matrix<1, 1, double> D = {{0}};

        StateSpace<2, 1, 1> sys{A, B, C, D};

        double wc = 1.0; // Desired crossover frequency

        auto pid_result = matlab::pidtune(sys, wc);
        REQUIRE(pid_result.has_value());

        CHECK(pid_result->Kp > 0);
        CHECK(pid_result->Ki > 0);
        CHECK(pid_result->Kd > 0);
        CHECK(pid_result->Kbc > 0); // Back-calculation gain
    }

    TEST_CASE("c2d returns optional and fails closed on bad TF") {
        TransferFunction<1, 2, double> good{{1.0}, {1.0, 1.0}}; // 1/(s+1)
        auto                           d = matlab::c2d(good, 0.01);
        REQUIRE(d.has_value());
        CHECK(d->Ts == doctest::Approx(0.01));
        CHECK(d->is_discrete());

        TransferFunction<1, 2, double> bad{{1.0}, {1.0, 0.0}}; // leading den = 0
        CHECK_FALSE(matlab::c2d(bad, 0.01).has_value());
        CHECK_FALSE(matlab::isstable(bad));
        CHECK(matlab::isstable(good));
    }

    TEST_CASE("MATLAB® structural analysis aliases") {
        Matrix<2, 2> A{{0.0, 1.0}, {-2.0, -3.0}};
        Matrix<2, 1> B{{0.0}, {1.0}};
        Matrix<1, 2> C{{1.0, 0.0}};

        auto Co_short = matlab::ctrb(A, B);
        auto Ob_short = matlab::obsv(A, C);
        auto Co_core = stability::controllability_matrix(A, B);
        auto Ob_core = stability::observability_matrix(A, C);

        CHECK(Co_short == Co_core);
        CHECK(Ob_short == Ob_core);
    }

    TEST_CASE("MATLAB® eig alias") {
        // A = [0 1; -2 -3] has eigenvalues -1 and -2.
        Matrix<2, 2> A{{0.0, 1.0}, {-2.0, -3.0}};
        auto         e = matlab::eig(A);
        // the 2x2 closed form returns (trace+sqrt)/2 first → -1, then -2
        CHECK(e[0].real() == doctest::Approx(-1.0));
        CHECK(e[1].real() == doctest::Approx(-2.0));
        CHECK(e[0].imag() == doctest::Approx(0.0));
    }

    TEST_CASE("MATLAB® tf helper") {
        auto tf_sys = matlab::tf(
            std::array<double, 2>{1.0, 1.0},
            std::array<double, 3>{1.0, 3.0, 2.0}
        );

        CHECK(tf_sys.num[0] == doctest::Approx(1.0));
        CHECK(tf_sys.num[1] == doctest::Approx(1.0));
        CHECK(tf_sys.den[0] == doctest::Approx(1.0));
        CHECK(tf_sys.den[1] == doctest::Approx(3.0));
        CHECK(tf_sys.den[2] == doctest::Approx(2.0));
    }

    TEST_CASE("MATLAB® tf helper accepts braced lists") {
        auto tf_sys = matlab::tf({1.0, 1.0}, {1.0, 3.0, 2.0});

        CHECK(tf_sys.num[0] == doctest::Approx(1.0));
        CHECK(tf_sys.num[1] == doctest::Approx(1.0));
        CHECK(tf_sys.den[0] == doctest::Approx(1.0));
        CHECK(tf_sys.den[1] == doctest::Approx(3.0));
        CHECK(tf_sys.den[2] == doctest::Approx(2.0));
    }

    TEST_CASE("MATLAB® pidstd maps Ti/Td to parallel Ki/Kd") {
        // Standard form: C = Kp*(1 + 1/(Ti s) + Td s) → Ki = Kp/Ti, Kd = Kp*Td
        constexpr double Kp = 4.0;
        constexpr double Ti = 2.0;
        constexpr double Td = 0.25;
        constexpr double N = 10.0; // filter coefficient → Tf = Td/N

        const auto r = matlab::pidstd(Kp, Ti, Td, N);
        CHECK(r.Kp == doctest::Approx(Kp));
        CHECK(r.Ki == doctest::Approx(Kp / Ti));
        CHECK(r.Kd == doctest::Approx(Kp * Td));
        CHECK(r.Tf == doctest::Approx(Td / N));
        CHECK(r.b == doctest::Approx(1.0));
        CHECK(r.c == doctest::Approx(1.0));

        // Ti <= 0 → no integral; default N=10 → Tf = Td/10
        const auto pd = matlab::pidstd(Kp, 0.0, Td);
        CHECK(pd.Ki == doctest::Approx(0.0));
        CHECK(pd.Kd == doctest::Approx(Kp * Td));
        CHECK(pd.Tf == doctest::Approx(Td / 10.0));

        // Explicit N = 0 → unfiltered D
        const auto unf = matlab::pidstd(Kp, Ti, Td, 0.0);
        CHECK(unf.Tf == doctest::Approx(0.0));
    }

    TEST_CASE("MATLAB® allmargin delay margin from phase margin") {
        // Unity-gain plant G = 1/(s+1); PM and gain-crossover should be finite.
        StateSpace<1, 1, 1> sys{};
        sys.A(0, 0) = -1.0;
        sys.B(0, 0) = 1.0;
        sys.C(0, 0) = 1.0;
        sys.Ts = 0.0;
        const auto w = matlab::logspace(-2.0, 2.0, 200);
        const auto a = matlab::allmargin(sys, w);
        // Dm = Pm_rad / Wcp when both finite
        if (a.Pm < 1e300 && a.Wcp > 0.0) {
            const double dm_ref = (a.Pm * damp::numbers::pi_v<double> / 180.0) / a.Wcp;
            CHECK(a.Dm == doctest::Approx(dm_ref).epsilon(1e-9));
        }
    }

    TEST_CASE("MATLAB® pidstd2 and make1DOF/make2DOF setpoint weights") {
        auto r2 = matlab::pidstd2(2.0, 1.0, 0.1, 0.0, 0.5, 0.0);
        CHECK(r2.Kp == doctest::Approx(2.0));
        CHECK(r2.Ki == doctest::Approx(2.0));
        CHECK(r2.Kd == doctest::Approx(0.2));
        CHECK(r2.b == doctest::Approx(0.5));
        CHECK(r2.c == doctest::Approx(0.0));

        auto one = matlab::make1DOF(r2);
        CHECK(one.b == doctest::Approx(1.0));
        CHECK(one.c == doctest::Approx(1.0));
        // original gains preserved
        CHECK(one.Kp == doctest::Approx(r2.Kp));
        CHECK(one.Ki == doctest::Approx(r2.Ki));

        auto two = matlab::make2DOF(one, 0.8, 0.0);
        CHECK(two.b == doctest::Approx(0.8));
        CHECK(two.c == doctest::Approx(0.0));
    }

    TEST_CASE("MATLAB® lqry matches dlqr with Q = C' Qy C") {
        // Double integrator (continuous), then discretize for discrete_lqr path
        Matrix<2, 2> A{{0.0, 1.0}, {0.0, 0.0}};
        Matrix<2, 1> B{{0.0}, {1.0}};
        Matrix<1, 2> C{{1.0, 0.0}};
        Matrix<1, 1> D{{0.0}};
        Matrix<1, 1> Qy{{4.0}};
        Matrix<1, 1> R{{1.0}};

        // Continuous output-weighted LQR
        auto cont = matlab::lqry(A, B, C, D, Qy, R);
        REQUIRE(cont.success);
        const Matrix<2, 2> Qx = C.transpose() * Qy * C;
        auto               cont_ref = design::continuous_lqr(A, B, Qx, R);
        REQUIRE(cont_ref.success);
        CHECK(cont.K(0, 0) == doctest::Approx(cont_ref.K(0, 0)).epsilon(1e-9));
        CHECK(cont.K(0, 1) == doctest::Approx(cont_ref.K(0, 1)).epsilon(1e-9));

        // Discrete StateSpace path: lqry(sys, Qy, R) ≡ discrete_lqr(A,B,C'Qy C,R)
        constexpr double Ts = 0.01;
        StateSpace       sys_c{A, B, C, D};
        auto             sys_d = *discretize(sys_c, Ts, DiscretizationMethod::ZOH);
        auto             disc = matlab::lqry(sys_d, Qy, R);
        REQUIRE(disc.success);
        auto disc_ref = design::discrete_lqr(sys_d.A, sys_d.B, Qx, R);
        REQUIRE(disc_ref.success);
        CHECK(disc.K(0, 0) == doctest::Approx(disc_ref.K(0, 0)).epsilon(1e-9));
        CHECK(disc.K(0, 1) == doctest::Approx(disc_ref.K(0, 1)).epsilon(1e-9));
        CHECK(disc.is_stable());
    }

    TEST_CASE("MATLAB® pade aliases design::pade_delay_*") {
        constexpr double T_delay = 0.05;
        const auto       p1 = matlab::pade(T_delay);
        const auto       ref1 = design::pade_delay_1st(T_delay);
        CHECK(p1.num[0] == doctest::Approx(ref1.num[0]));
        CHECK(p1.num[1] == doctest::Approx(ref1.num[1]));
        CHECK(p1.den[0] == doctest::Approx(ref1.den[0]));
        CHECK(p1.den[1] == doctest::Approx(ref1.den[1]));

        const auto p2 = matlab::pade2(T_delay);
        const auto ref2 = design::pade_delay_2nd(T_delay);
        CHECK(p2.num[0] == doctest::Approx(ref2.num[0]));
        CHECK(p2.num[1] == doctest::Approx(ref2.num[1]));
        CHECK(p2.num[2] == doctest::Approx(ref2.num[2]));
        CHECK(p2.den[0] == doctest::Approx(ref2.den[0]));
        CHECK(p2.den[1] == doctest::Approx(ref2.den[1]));
        CHECK(p2.den[2] == doctest::Approx(ref2.den[2]));
    }

    TEST_CASE("MATLAB® isstable continuous and discrete") {
        // Stable continuous: poles at -1, -2
        Matrix<2, 2> A_stable{{0.0, 1.0}, {-2.0, -3.0}};
        CHECK(matlab::isstable(A_stable));

        // Unstable continuous: pole at +1
        Matrix<2, 2> A_unstable{{0.0, 1.0}, {1.0, 0.0}};
        CHECK_FALSE(matlab::isstable(A_unstable));

        StateSpace sys_c{
            .A = A_stable,
            .B = Matrix<2, 1>{{0.0}, {1.0}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
            .D = Matrix<1, 1>{{0.0}},
            .Ts = 0.0
        };
        CHECK(matlab::isstable(sys_c));

        // Discrete stable: |λ| < 1
        Matrix<2, 2> Ad{{0.5, 0.1}, {0.0, 0.8}};
        StateSpace   sys_d{
              .A = Ad,
              .B = Matrix<2, 1>{{0.0}, {1.0}},
              .C = Matrix<1, 2>{{1.0, 0.0}},
              .D = Matrix<1, 1>{{0.0}},
              .Ts = 0.01
        };
        CHECK(matlab::isstable(sys_d));

        // Discrete unstable
        sys_d.A = Matrix<2, 2>{{1.2, 0.0}, {0.0, 0.5}};
        CHECK_FALSE(matlab::isstable(sys_d));
    }

} // TEST_SUITE("matlab")

TEST_SUITE("MATLAB®-Style Control Design API") {
    // Test dlqr: Discrete LQR from A, B, Q, R matrices
    TEST_CASE("dlqr: discrete LQR design") {
        // Double integrator (discrete): x[k+1] = [1 0.1; 0 1]*x[k] + [0.005; 0.1]*u[k]
        constexpr double Ts = 0.1;
        Matrix<2, 2>     Ad{{1.0, Ts}, {0.0, 1.0}};
        Matrix<2, 1>     Bd{{Ts * Ts / 2}, {Ts}};
        Matrix<2, 2>     Q{{1.0, 0.0}, {0.0, 1.0}};
        Matrix<1, 1>     R{{0.1}};

        auto result = matlab::dlqr(Ad, Bd, Q, R);

        REQUIRE(result.success);
        CHECK(result.K(0, 0) != 0.0);
        CHECK(result.K(0, 1) != 0.0);

        // S should be positive definite (check diagonal)
        CHECK(result.S(0, 0) > 0.0);
        CHECK(result.S(1, 1) > 0.0);
    }

    // Test lqrd: discrete LQR from a continuous-time system
    TEST_CASE("lqrd: discrete LQR from continuous design") {
        Matrix<2, 2> A{{0.0, 1.0}, {0.0, 0.0}}; // continuous double integrator
        Matrix<2, 1> B{{0.0}, {1.0}};
        Matrix<2, 2> Q = Matrix<2, 2>::identity();
        Matrix<1, 1> R{{1.0}};

        auto result = matlab::lqrd(A, B, Q, R, 0.1);

        REQUIRE(result.success);
        CHECK(result.is_stable());
        CHECK(result.K(0, 0) != 0.0);
    }

    // Test lqi: LQI controller design
    TEST_CASE("lqi: LQI controller design") {
        StateSpace<2, 1, 1, double, 2, 1> sys{
            Matrix<2, 2>{{0.0, 1.0}, {0.0, 0.0}}, // double integrator
            Matrix<2, 1>{{0.0}, {1.0}},
            Matrix<1, 2>{{1.0, 0.0}} // position output
        };

        // Augmented cost: [state; integral]
        Matrix<3, 3> Q_aug{};
        Q_aug(0, 0) = 1.0;  // position
        Q_aug(1, 1) = 1.0;  // velocity
        Q_aug(2, 2) = 10.0; // integral (high weight for tracking)
        Matrix<1, 1> R{{0.1}};

        auto lqi_result = matlab::lqi(sys, Q_aug, R);

        REQUIRE(lqi_result.success);
        CHECK(lqi_result.K(0, 0) != 0.0);
        CHECK(lqi_result.K(0, 1) != 0.0);
        CHECK(lqi_result.K(0, 2) != 0.0);
    }

    // Test lqg: LQG regulator design
    TEST_CASE("lqg: LQG regulator design") {
        StateSpace sys{
            .A = Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}}, // discrete double integrator
            .B = Matrix<2, 1>{{0.005}, {0.1}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
            .D = Matrix<1, 1>::zeros(),
            .G = Matrix<2, 2>::identity(),
            .H = Matrix<1, 1>::identity(),
        };

        Matrix<2, 2> Q_lqr = Matrix<2, 2>::identity();
        Matrix<1, 1> R_lqr{{0.1}};
        Matrix<2, 2> Q_kf{{0.01, 0.0}, {0.0, 0.01}};
        Matrix<1, 1> R_kf{{0.1}};

        auto lqg_result = matlab::lqg(sys, Q_lqr, R_lqr, Q_kf, R_kf);

        REQUIRE(lqg_result.success);
        CHECK(lqg_result.lqr.K(0, 0) != 0.0);
        CHECK(lqg_result.lqr.K(0, 1) != 0.0);
    }

    // Test lqgtrack: LQG servo controller
    TEST_CASE("lqgtrack: LQG servo design") {
        StateSpace sys{
            Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}},
            Matrix<2, 1>{{0.005}, {0.1}},
            Matrix<1, 2>{{1.0, 0.0}},
            Matrix<1, 1>::zeros(),
            Matrix<2, 2>::identity(),
            Matrix<1, 1>::identity(),
        };

        Matrix<3, 3> Q_aug{};
        Q_aug(0, 0) = 1.0;
        Q_aug(1, 1) = 1.0;
        Q_aug(2, 2) = 10.0;
        Matrix<1, 1> R{{0.1}};
        Matrix<2, 2> Q_kf{{0.01, 0.0}, {0.0, 0.01}};
        Matrix<1, 1> R_kf{{0.1}};

        auto servo_result = matlab::lqgtrack(sys, Q_aug, R, Q_kf, R_kf);

        // Should have computed gains
        CHECK(servo_result.lqi.K(0, 0) != 0.0);
        CHECK(servo_result.lqi.K(0, 1) != 0.0);
        CHECK(servo_result.lqi.K(0, 2) != 0.0);
    }

    // Test lqgreg: Form LQG from Kalman and LQR results
    TEST_CASE("lqgreg: combine Kalman result and LQR result") {
        StateSpace sys{
            .A = Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}},
            .B = Matrix<2, 1>{{0.005}, {0.1}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
            .D = Matrix<1, 1>::zeros(),
            .G = Matrix<2, 2>::identity(),
            .H = Matrix<1, 1>::identity(),
        };

        // Create Kalman result
        Matrix<2, 2> Q_kf{{0.01, 0.0}, {0.0, 0.01}};
        Matrix<1, 1> R_kf{{0.1}};
        auto         kf_result = design::kalman(sys, Q_kf, R_kf);

        // Create LQR result using dlqr free function
        Matrix<2, 2> Q_lqr = Matrix<2, 2>::identity();
        Matrix<1, 1> R_lqr{{0.1}};
        auto         lqr_result = matlab::dlqr(sys.A, sys.B, Q_lqr, R_lqr);

        // Combine into LQGResult
        auto lqg_result = matlab::lqgreg(kf_result, lqr_result);

        // Should preserve the LQR gain
        CHECK(lqg_result.lqr.K(0, 0) == lqr_result.K(0, 0));
        CHECK(lqg_result.lqr.K(0, 1) == lqr_result.K(0, 1));
    }

    // Test compile-time design (consteval)
    TEST_CASE("design:: consteval functions compile") {
        // This test verifies that design:: functions can be evaluated at compile time
        constexpr Matrix<2, 2> Ad{{1.0, 0.1}, {0.0, 1.0}};
        constexpr Matrix<2, 1> Bd{{0.005}, {0.1}};
        constexpr Matrix<2, 2> Q{{1.0, 0.0}, {0.0, 1.0}};
        constexpr Matrix<1, 1> R{{0.1}};

        // Note: We can't use consteval directly in runtime CHECK, but we can
        // verify the design functions work at compile time by using constexpr
        constexpr auto result = design::discrete_lqr(Ad, Bd, Q, R);

        // At runtime, verify the compile-time result
        CHECK(result.K(0, 0) != 0.0);
        CHECK(result.K(0, 1) != 0.0);
        CHECK(result.S(0, 0) > 0.0);
    }
}
