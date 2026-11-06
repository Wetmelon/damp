// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <type_traits>
#include <vector>

#include "damp/backend.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/controllers/lqr.hpp"
#include "damp/controllers/pid.hpp"
#include "damp/design/riccati.hpp"
#include "damp/design/stability.hpp"
#include "damp/design/status.hpp"
#include "damp/estimation/kalman.hpp"
#include "damp/filters/pade.hpp"
#include "damp/math/complex.hpp"
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
        REQUIRE(K_opt.success);
        const auto& K = K_opt.K;

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
        REQUIRE(K_opt.success);
        const auto& K = K_opt.K;

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
        REQUIRE(K_opt.success);
        const auto& K = K_opt.K;

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
        REQUIRE(pid_result.success);

        CHECK(pid_result.Kp > 0);
        CHECK(pid_result.Ki > 0);
        CHECK(pid_result.Kd > 0);
        CHECK(pid_result.Kbc > 0); // Back-calculation gain
    }

    TEST_CASE("c2d returns a DiscretizationResult and fails closed on bad TF") {
        TransferFunction<1, 2, double> good{{1.0}, {1.0, 1.0}}; // 1/(s+1)
        auto                           d = matlab::c2d(good, 0.01);
        REQUIRE(d.has_value());
        CHECK(d.status == design::DesignStatus::ok);
        CHECK(d->Ts == doctest::Approx(0.01));
        CHECK(d->is_discrete());

        TransferFunction<1, 2, double> bad{{1.0}, {0.0, 1.0}}; // leading den = 0
        CHECK_FALSE(matlab::c2d(bad, 0.01).has_value());
        CHECK(matlab::c2d(bad, 0.01).status == design::DesignStatus::singular);
        const StateSpace<1, 1, 1> pole_at_2_over_Ts{.A = {{4.0}}, .B = {{1.0}}, .C = {{1.0}}};
        CHECK(matlab::c2d(pole_at_2_over_Ts, 0.5, DiscretizationMethod::Tustin).status == design::DesignStatus::singular);
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
        const auto   e = matlab::eig(A);
        REQUIRE(e.has_value());
        // the 2x2 closed form returns (trace+sqrt)/2 first → -1, then -2
        CHECK((*e)[0].real() == doctest::Approx(-1.0));
        CHECK((*e)[1].real() == doctest::Approx(-2.0));
        CHECK((*e)[0].imag() == doctest::Approx(0.0));
        CHECK_FALSE(matlab::eig(Matrix<2, 2>{{std::numeric_limits<double>::quiet_NaN(), 0.0}, {0.0, 1.0}}).has_value());
        CHECK_FALSE(matlab::pole(Matrix<2, 2>{{std::numeric_limits<double>::quiet_NaN(), 0.0}, {0.0, 1.0}}).has_value());
    }

    TEST_CASE("MATLAB® expm, sqrtm and logm aliases") {
        // A = diag(1, 4): expm = diag(e, e⁴), sqrtm = diag(1, 2), logm = diag(0, ln 4).
        const Matrix<2, 2> A{{1.0, 0.0}, {0.0, 4.0}};
        const auto         E = matlab::expm(A);
        const auto         S = matlab::sqrtm(A);
        const auto         L = matlab::logm(A);
        REQUIRE(E);
        REQUIRE(S);
        REQUIRE(L);
        CHECK(*E == *mat::expm(A));
        CHECK(*S == *mat::sqrt(A));
        CHECK(*L == *mat::log(A));
        CHECK((*S)(1, 1) == doctest::Approx(2.0));
        CHECK((*L)(1, 1) == doctest::Approx(std::log(4.0)));
        CHECK_FALSE(matlab::sqrtm(Matrix<2, 2>{{-1.0, 0.0}, {0.0, 1.0}}));
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

        // MATLAB® tf([1 1], [1 3 2]) = (s+1)/((s+1)(s+2)): zero −1, poles −1 and −2
        const auto zpk = matlab::tf2zpk(tf_sys);
        REQUIRE(zpk.success);
        CHECK(zpk.zpk.gain == doctest::Approx(1.0));
        CHECK(zpk.zpk.zeros[0].real() == doctest::Approx(-1.0));
        CHECK(damp::max(zpk.zpk.poles[0].real(), zpk.zpk.poles[1].real()) == doctest::Approx(-1.0));
        CHECK(damp::min(zpk.zpk.poles[0].real(), zpk.zpk.poles[1].real()) == doctest::Approx(-2.0));
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

        // Ti <= 0 → no integral; default N is MATLAB®'s ∞ → unfiltered D (Tf = 0)
        const auto pd = matlab::pidstd(Kp, 0.0, Td);
        CHECK(pd.Ki == doctest::Approx(0.0));
        CHECK(pd.Kd == doctest::Approx(Kp * Td));
        CHECK(pd.Tf == doctest::Approx(0.0));

        // Explicit N = 0 → unfiltered D
        const auto unf = matlab::pidstd(Kp, Ti, Td, 0.0);
        CHECK(unf.Tf == doctest::Approx(0.0));
    }

    TEST_CASE("MATLAB® logspace takes decade exponents") {
        // MATLAB® logspace(-2, 2, 5) = 0.01 0.1 1 10 100; the value-endpoint reading gave 1 … 2.
        const auto w = matlab::logspace(-2.0, 2.0, 5);
        REQUIRE(w.size() == 5);
        CHECK(w[0] == doctest::Approx(0.01));
        CHECK(w[1] == doctest::Approx(0.1));
        CHECK(w[2] == doctest::Approx(1.0));
        CHECK(w[4] == doctest::Approx(100.0));
    }

    TEST_CASE("MATLAB® allmargin delay margin from phase margin") {
        // L = 4/(s+1)³: |L| crosses 0 dB at ω = √(4^(2/3) − 1) ≈ 1.233 rad/s, PM ≈ 27.1°.
        const auto L = TransferFunction<1, 4>{.num = {4.0}, .den = {1.0, 3.0, 3.0, 1.0}}.to_state_space();
        REQUIRE(L.has_value());
        const auto a = matlab::allmargin(*L, matlab::logspace(-2.0, 2.0, 2000));
        CHECK(a.Wcp == doctest::Approx(1.2328).epsilon(1e-3));
        CHECK(a.Pm == doctest::Approx(27.1).epsilon(0.01));
        const double dm_ref = (a.Pm * damp::numbers::pi_v<double> / 180.0) / a.Wcp;
        CHECK(a.Dm == doctest::Approx(dm_ref).epsilon(1e-9));
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

        // MATLAB® defaults b = c = 1: make2DOF(C) and pidstd2 without weights keep the 1-DOF response.
        CHECK(matlab::make2DOF(one).c == doctest::Approx(1.0));
        CHECK(matlab::pidstd2(2.0, 1.0, 0.1).c == doctest::Approx(1.0));
    }

    TEST_CASE("MATLAB® pid is a design result, so make2DOF and discretize take it") {
        const auto c = matlab::pid(1.0, 2.0, 0.1, 0.01);
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(c)>, design::PIDGains<double>>);
        CHECK(c.Tf == doctest::Approx(0.01));
        const auto two = matlab::make2DOF(c, 1.0, 0.0);
        CHECK(two.Kd == doctest::Approx(0.1));
        CHECK(two.c == doctest::Approx(0.0));
        const ContinuousPID<double> runtime{c};
        CHECK(runtime.Kp == doctest::Approx(1.0));
    }

    TEST_CASE("MATLAB® ones, eye and minreal's default tolerance") {
        const auto o = matlab::ones<3>();
        const auto I = matlab::eye<3>();
        CHECK(o(1, 2) == doctest::Approx(1.0));
        CHECK(I(1, 2) == doctest::Approx(0.0));
        CHECK(I(2, 2) == doctest::Approx(1.0));
        // A mode coupled at 1e-10 on both B and C (1e-10 after balancing too) sits between √ε and
        // 1e-12: the alias once kept it while design::minreal removed it.
        const StateSpace<2, 1, 1> weak{
            .A = Matrix<2, 2>{{-1.0, 0.0}, {0.0, -2.0}},
            .B = Matrix<2, 1>{{1.0}, {1e-10}},
            .C = Matrix<1, 2>{{1.0, 1e-10}},
        };
        CHECK(matlab::minreal(weak).order == design::minreal(weak).order);
        CHECK(matlab::minreal(weak).order == 1);
    }

    TEST_CASE("MATLAB® estim carries Ts and places the estimator poles at eig(A − LC)") {
        const StateSpace<2, 1, 1> plant{
            .A = Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}},
            .B = Matrix<2, 1>{{0.005}, {0.1}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
            .Ts = 0.1,
        };
        const Matrix<2, 1> L{{0.9}, {2.0}};
        const auto         est = matlab::estim(plant, L);
        CHECK(est.Ts == doctest::Approx(0.1));
        CHECK(est.is_discrete());
        const Matrix<2, 2> Acl = plant.A - (L * plant.C);
        CHECK(est.A(0, 0) == doctest::Approx(Acl(0, 0)));
        CHECK(est.A(1, 0) == doctest::Approx(Acl(1, 0)));
        CHECK(est.C(1, 1) == doctest::Approx(1.0));
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

    TEST_CASE("MATLAB® pade<N> aliases design::pade_delay_*") {
        constexpr double T_delay = 0.05;
        const auto       p1 = matlab::pade<1>(T_delay).tf;
        const auto       ref1 = design::pade_delay_1st(T_delay).tf;
        CHECK(p1.num[0] == doctest::Approx(ref1.num[0]));
        CHECK(p1.num[1] == doctest::Approx(ref1.num[1]));
        CHECK(p1.den[0] == doctest::Approx(ref1.den[0]));
        CHECK(p1.den[1] == doctest::Approx(ref1.den[1]));

        const auto p2 = matlab::pade<2>(T_delay).tf;
        const auto ref2 = design::pade_delay_2nd(T_delay).tf;
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
        CHECK(matlab::isstable(StateSpace<2, 1, 1>{.A = A_stable}));

        // Unstable continuous: pole at +1
        Matrix<2, 2> A_unstable{{0.0, 1.0}, {1.0, 0.0}};
        CHECK_FALSE(matlab::isstable(StateSpace<2, 1, 1>{.A = A_unstable}));

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

    TEST_CASE("tf2sos reproduces MATLAB's example") {
        // MATLAB: sos = tf2sos([1 3 3 1]/6, [3 0 1 0]/3)
        //   0.1667  0.1667  0       1  0  0
        //   1       2       1       1  0  0.3333
        const auto sos = matlab::tf2sos({1.0 / 6, 0.5, 0.5, 1.0 / 6}, {1.0, 0.0, 1.0 / 3, 0.0});
        REQUIRE(sos.success);
        CHECK(sos.sections[0].b0 == doctest::Approx(1.0 / 6));
        CHECK(sos.sections[0].b1 == doctest::Approx(1.0 / 6));
        CHECK(sos.sections[0].b2 == doctest::Approx(0.0).epsilon(1e-12).scale(1.0));
        CHECK(sos.sections[0].a1 == doctest::Approx(0.0).epsilon(1e-12).scale(1.0));
        CHECK(sos.sections[0].a2 == doctest::Approx(0.0).epsilon(1e-12).scale(1.0));
        CHECK(sos.sections[1].b0 == doctest::Approx(1.0));
        CHECK(sos.sections[1].b1 == doctest::Approx(2.0));
        CHECK(sos.sections[1].b2 == doctest::Approx(1.0));
        CHECK(sos.sections[1].a1 == doctest::Approx(0.0).epsilon(1e-12).scale(1.0));
        CHECK(sos.sections[1].a2 == doctest::Approx(1.0 / 3));
    }

    TEST_CASE("zp2sos matches second_order_sections of the TF") {
        const damp::array<damp::complex<double>, 2> z{damp::complex<double>{-1.0, 0.0}, damp::complex<double>{-1.0, 0.0}};
        const damp::array<damp::complex<double>, 2> p{damp::complex<double>{0.5, 0.3}, damp::complex<double>{0.5, -0.3}};
        const auto                                  sos = matlab::zp2sos(z, p, 0.2);
        REQUIRE(sos.success);
        CHECK(sos.sections[0].b0 == doctest::Approx(0.2));
        CHECK(sos.sections[0].b1 == doctest::Approx(0.4));
        CHECK(sos.sections[0].b2 == doctest::Approx(0.2));
        CHECK(sos.sections[0].a1 == doctest::Approx(-1.0));
        CHECK(sos.sections[0].a2 == doctest::Approx(0.34));
    }


    TEST_CASE("matlab::reg builds the y→u compensator (separation principle)") {
        const Matrix<2, 2> A{{0.0, 1.0}, {2.0, -1.0}};
        const Matrix<2, 1> B{{0.0}, {1.0}};
        const Matrix<1, 2> C{{1.0, 0.0}};
        const Matrix<1, 1> D{{0.5}};
        const auto         sys = matlab::ss(A, B, C, D);
        const Matrix<1, 2> K{{8.0, 4.0}};    // eig(A − BK) = {−2, −3}
        const Matrix<2, 1> L{{9.0}, {17.0}}; // eig(A − LC) = {−4, −6}
        const auto         rsys = matlab::reg(sys, K, L);
        // Close the loop: plant u = C_r ξ + D_r y, compensator input y.
        Matrix<4, 4> Acl{};
        Acl.template block<2, 2>(0, 0) = A;
        Acl.template block<2, 2>(0, 2) = B * rsys.C;
        Acl.template block<2, 2>(2, 0) = rsys.B * C;
        Acl.template block<2, 2>(2, 2) = rsys.A + rsys.B * D * rsys.C;
        auto                ev = mat::compute_eigenvalues(Acl).values;
        std::vector<double> re;
        for (std::size_t i = 0; i < 4; ++i) {
            CHECK(std::abs(ev[i].imag()) < 1e-9);
            re.push_back(ev[i].real());
        }
        std::sort(re.begin(), re.end());
        CHECK(re[0] == doctest::Approx(-6.0));
        CHECK(re[1] == doctest::Approx(-4.0));
        CHECK(re[2] == doctest::Approx(-3.0));
        CHECK(re[3] == doctest::Approx(-2.0));
    }

    TEST_CASE("matlab::lqi / lqg refuse a continuous plant") {
        const auto sys = matlab::ss(Matrix<1, 1>{{-1.0}}, Matrix<1, 1>{{1.0}}, Matrix<1, 1>{{1.0}});
        CHECK_FALSE(matlab::lqi(sys, Matrix<2, 2>::identity(), Matrix<1, 1>{{1.0}}).success);
    }

    TEST_CASE("matlab::pidstd default N is MATLAB's (no derivative filter)") {
        const auto r = matlab::pidstd(2.0, 1.0, 0.5);
        CHECK(r.Tf == 0.0);
    }

    TEST_CASE("matlab::lqr is fail-closed") {
        // Unstable, uncontrollable mode: no stabilizing gain exists.
        const Matrix<2, 2> A{{1.0, 0.0}, {0.0, -1.0}};
        const Matrix<2, 1> B{{0.0}, {1.0}};
        CHECK_FALSE(matlab::lqr(A, B, Matrix<2, 2>::identity(), Matrix<1, 1>{{1.0}}).success);
    }
} // TEST_SUITE("matlab")

TEST_SUITE("MATLAB®-Style Control Design API") {
    TEST_CASE("lqr returns the continuous_lqr result") {
        const Matrix<2, 2> A{{0.0, 1.0}, {0.0, 0.0}};
        const Matrix<2, 1> B{{0.0}, {1.0}};
        const Matrix<2, 2> Q = Matrix<2, 2>::identity();
        const Matrix<1, 1> R{{1.0}};

        const auto result = matlab::lqr(A, B, Q, R);
        const auto ref = design::continuous_lqr(A, B, Q, R);
        REQUIRE(result.success);
        CHECK(result.status == design::DesignStatus::ok);
        CHECK(result.status == ref.status);
        CHECK(result.K(0, 0) == doctest::Approx(ref.K(0, 0)));
        CHECK(result.K(0, 1) == doctest::Approx(ref.K(0, 1)));
        CHECK(result.S(0, 0) == doctest::Approx(ref.S(0, 0)));
        CHECK(result.S(1, 1) == doctest::Approx(ref.S(1, 1)));

        // Uncontrollable unstable mode: no stabilizing gain. The Result fails closed.
        const Matrix<2, 2> Au{{1.0, 0.0}, {0.0, -1.0}};
        const Matrix<2, 1> Bu{{0.0}, {1.0}};
        const auto         bad = matlab::lqr(Au, Bu, Matrix<2, 2>::identity(), Matrix<1, 1>{{1.0}});
        const auto         bad_ref = design::continuous_lqr(Au, Bu, Matrix<2, 2>::identity(), Matrix<1, 1>{{1.0}});
        CHECK_FALSE(bad.success);
        CHECK(bad.status == bad_ref.status);
        CHECK(bad.status != design::DesignStatus::ok);
    }

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

    TEST_CASE("dlqe: matrix form forwards to design::kalman with G and H = I") {
        const Matrix<2, 2> A{{0.9, 0.1}, {0.0, 0.8}};
        const Matrix<2, 2> G = Matrix<2, 2>::identity();
        const Matrix<1, 2> C{{1.0, 0.0}};
        const Matrix<2, 2> Q = Matrix<2, 2>::diagonal({0.1, 0.1});
        const Matrix<1, 1> R{{0.5}};

        const StateSpace<2, 1, 1, double, 2, 1> sys{.A = A, .C = C, .G = G, .H = Matrix<1, 1>{{1.0}}, .Ts = 0.01};
        const auto                              via_kalman = design::kalman(sys, {.process = Q, .measurement = R});
        const auto                              via_dlqe = matlab::dlqe(A, G, C, Q, R);
        REQUIRE(via_kalman.success);
        REQUIRE(via_dlqe.success);
        CHECK(via_dlqe.M(0, 0) == doctest::Approx(via_kalman.M(0, 0)));
        CHECK(via_dlqe.M(1, 0) == doctest::Approx(via_kalman.M(1, 0))); // dlqe's M, the innovation gain
        CHECK(via_dlqe.P(0, 0) == doctest::Approx(via_kalman.P(0, 0)));
        CHECK(via_dlqe.P(1, 1) == doctest::Approx(via_kalman.P(1, 1)));
    }

    TEST_CASE("d2d / upsample: forward to resample / upsample") {
        constexpr StateSpace<1, 1, 1> sys{.A = Matrix<1, 1>{{0.9}}, .B = Matrix<1, 1>{{0.1}}, .C = Matrix<1, 1>{{1.0}}, .Ts = 0.1};
        constexpr auto                slow = matlab::d2d(sys, 0.2);
        static_assert(slow.has_value() && slow->Ts == 0.2);
        CHECK(slow->A(0, 0) == doctest::Approx(0.81));
        CHECK(slow->B(0, 0) == doctest::Approx(0.19));
        static_assert(matlab::d2d(sys, -0.2).status == design::DesignStatus::invalid_sample_time);
        constexpr auto up = matlab::upsample<2>(sys);
        static_assert(up.Ts == 0.05 && up.A(0, 1) == 0.9 && up.A(1, 0) == 1.0);
    }

    TEST_CASE("kalmd: MATLAB's Rn/Ts, then forwards to design::kalman(sys, Ts, noise)") {
        constexpr StateSpace<2, 1, 1, double, 1, 1> sys{
            .A = Matrix<2, 2>{{0.0, 1.0}, {-2.0, -0.5}},
            .B = Matrix<2, 1>{{0.0}, {1.0}},
            .C = Matrix<1, 2>{{1.0, 0.0}},
            .G = Matrix<2, 1>{{0.0}, {1.0}},
            .H = Matrix<1, 1>{{1.0}},
        };
        constexpr auto kest = matlab::kalmd(sys, Matrix<1, 1>{{0.3}}, Matrix<1, 1>{{0.02}}, 0.1);
        static_assert(kest.success);
        CHECK(kest.R(0, 0) == doctest::Approx(0.2).epsilon(1e-14));                // Rd = Rn/Ts, as kalmd
        CHECK(kest.M(0, 0) == doctest::Approx(0.15876599438771183).epsilon(1e-9)); // SciPy, see test_kalman.cpp
        CHECK(kest.M(1, 0) == doctest::Approx(0.13651425878842288).epsilon(1e-9));
        CHECK((kest.L - (kest.sys.A * kest.M)).norm() == 0.0); // MATLAB®'s L, the predictor gain

        //! design::kalman takes the per-reading variance itself: kalmd's Rd = 0.2 gives the same filter.
        constexpr auto same = design::kalman(sys, 0.1, {.process = 0.3, .measurement = 0.2});
        CHECK(same.M(0, 0) == doctest::Approx(kest.M(0, 0)).epsilon(1e-12));
        //! ...while passing kalmd's Rn there is a 10× more trusted sensor, a different (larger) gain.
        constexpr auto unscaled = design::kalman(sys, 0.1, {.process = 0.3, .measurement = 0.02});
        CHECK(unscaled.M(0, 0) > kest.M(0, 0));

        //! A bad Ts is refused, not divided by.
        constexpr auto bad = matlab::kalmd(sys, Matrix<1, 1>{{0.3}}, Matrix<1, 1>{{0.02}}, 0.0);
        static_assert(bad.status == design::DesignStatus::invalid_sample_time);
    }

    // Test lqi: LQI controller design
    TEST_CASE("lqi: LQI controller design") {
        // lqi is a discrete design: ZOH double integrator at Ts = 0.1.
        StateSpace<2, 1, 1, double, 2, 1> sys{
            .A = Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}},
            .B = Matrix<2, 1>{{0.005}, {0.1}},
            .C = Matrix<1, 2>{{1.0, 0.0}}, // position output
            .Ts = 0.1,
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
            .Ts = 0.1,
        };

        Matrix<2, 2> Q_lqr = Matrix<2, 2>::identity();
        Matrix<1, 1> R_lqr{{0.1}};
        Matrix<2, 2> Q_kf{{0.01, 0.0}, {0.0, 0.01}};
        Matrix<1, 1> R_kf{{0.1}};

        auto lqg_result = matlab::lqg(sys, Q_lqr, R_lqr, Q_kf, R_kf);

        REQUIRE(lqg_result.success);
        CHECK(lqg_result.K(0, 0) != 0.0);
        CHECK(lqg_result.K(0, 1) != 0.0);
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
            0.1, // Ts: discrete double integrator
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
        CHECK(servo_result.K(0, 0) != 0.0);
        CHECK(servo_result.K(0, 1) != 0.0);
        CHECK(servo_result.K(0, 2) != 0.0);
    }

    TEST_CASE("lqi / lqgtrack integrate a subset of the outputs, same argument as the design calls") {
        // Position and velocity sensors; integrate the position only.
        const StateSpace<2, 1, 2, double, 2, 2> sys{
            .A = Matrix<2, 2>{{1.0, 0.1}, {0.0, 1.0}},
            .B = Matrix<2, 1>{{0.005}, {0.1}},
            .C = Matrix<2, 2>::identity(),
            .G = Matrix<2, 2>::identity(),
            .H = Matrix<2, 2>::identity(),
            .Ts = 0.1,
        };
        const Matrix<3, 3> Q_aug{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 10.0}};
        const Matrix<1, 1> R{{0.1}};
        const Matrix<2, 2> Q_kf{{0.01, 0.0}, {0.0, 0.01}};
        const Matrix<2, 2> R_kf{{0.1, 0.0}, {0.0, 0.1}};
        const Matrix<1, 2> pick_position{{1.0, 0.0}};

        const auto direct = design::lqi(sys, {.Q = Matrix<2, 2>::identity(), .R = R, .Q_integral = 10.0}, design::IntegralAction<1>{{0}});
        const auto by_index = matlab::lqi(sys, Q_aug, R, {0});
        const auto by_selector = matlab::lqi(sys, Q_aug, R, pick_position);
        const auto track_index = matlab::lqgtrack(sys, Q_aug, R, Q_kf, R_kf, {0});
        const auto track_selector = matlab::lqgtrack(sys, Q_aug, R, Q_kf, R_kf, pick_position);
        REQUIRE(direct.success);
        REQUIRE(by_index.success);
        REQUIRE(by_selector.success);
        REQUIRE(track_index.success);
        REQUIRE(track_selector.success);
        for (size_t j = 0; j < 3; ++j) {
            CHECK(by_index.K(0, j) == doctest::Approx(direct.K(0, j)));
            CHECK(by_selector.K(0, j) == doctest::Approx(direct.K(0, j)));
            CHECK(track_index.K(0, j) == doctest::Approx(direct.K(0, j)));
            CHECK(track_selector.K(0, j) == doctest::Approx(direct.K(0, j)));
        }
        CHECK(by_index.integral.Ci(0, 1) == 0.0);
        CHECK(track_index.integral.Ci(0, 0) == 1.0);
        CHECK_FALSE(matlab::lqi(sys, Q_aug, R, {2}).success);
        CHECK_FALSE(matlab::lqgtrack(sys, Q_aug, R, Q_kf, R_kf, {2}).success);

        // Continuous plants are refused, with the subset as without it.
        auto cont = sys;
        cont.Ts = 0.0;
        CHECK_FALSE(matlab::lqi(cont, Q_aug, R, {0}).success);
        CHECK_FALSE(matlab::lqgtrack(cont, Q_aug, R, Q_kf, R_kf, pick_position).success);
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
            .Ts = 0.1,
        };

        // Create Kalman result
        Matrix<2, 2> Q_kf{{0.01, 0.0}, {0.0, 0.01}};
        Matrix<1, 1> R_kf{{0.1}};
        auto         kf_result = design::kalman(sys, {.process = Q_kf, .measurement = R_kf});

        // The regulator on the same discrete plant
        Matrix<2, 2> Q_lqr = Matrix<2, 2>::identity();
        Matrix<1, 1> R_lqr{{0.1}};
        auto         lqr_result = design::lqr(sys, {.Q = Q_lqr, .R = R_lqr});

        auto lqg_result = matlab::lqgreg(kf_result, lqr_result);
        CHECK(lqg_result.success);

        // Should preserve the LQR gain
        CHECK(lqg_result.K(0, 0) == lqr_result.K(0, 0));
        CHECK(lqg_result.K(0, 1) == lqr_result.K(0, 1));

        // MATLAB®'s lqgreg(kest, dlqr(A, B, Q, R)): the bare-matrix gain takes the filter's Ts
        const auto from_dlqr = matlab::lqgreg(kf_result, matlab::dlqr(sys.A, sys.B, Q_lqr, R_lqr));
        CHECK(from_dlqr.success);
        CHECK(from_dlqr.Ts == kf_result.sys.Ts);
        CHECK(from_dlqr.K(0, 0) == doctest::Approx(lqr_result.K(0, 0)));
        CHECK(from_dlqr.K(0, 1) == doctest::Approx(lqr_result.K(0, 1)));
        const auto from_care = matlab::lqgreg(kf_result, design::continuous_lqr(sys.A, sys.B, Q_lqr, R_lqr));
        CHECK(from_care.status == design::DesignStatus::invalid_sample_time);
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
