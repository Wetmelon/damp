// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cstddef>
#include <limits>
#include <type_traits>

#include "damp/backend.hpp"
#include "damp/controllers/adrc.hpp"
#include "damp/design/status.hpp"
#include "damp/estimation/eso.hpp"
#include "damp/math/math.hpp"
#include "damp/systems/zpk.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

namespace {

constexpr auto                 kEso2 = design::eso<2>(50.0, 2.0, 1e-3);
static constinit ESO<2, float> kObs2{kEso2};

/// Coefficients of (z − λ)ⁿ, highest power first.
template<size_t N>
damp::array<double, N + 1> repeated_root_poly(double lambda) {
    damp::array<double, N + 1> c{};
    c[0] = 1.0;
    for (size_t k = 0; k < N; ++k) {
        for (size_t i = k + 1; i > 0; --i) {
            c[i] = c[i] - (lambda * c[i - 1]);
        }
    }
    return c;
}

/// ÿ = a, advanced exactly over one period of held a.
struct DoubleIntegrator {
    double y{0.0};
    double v{0.0};

    void step(double a, double Ts) {
        y += (v * Ts) + (0.5 * a * Ts * Ts);
        v += a * Ts;
    }
};

} // namespace

TEST_SUITE("Extended State Observer (ESO)") {
    TEST_CASE("1st-order binomial β") {
        constexpr auto r = design::eso<1>(20.0, 1.0);
        static_assert(r.success);
        CHECK(r.beta[0] == doctest::Approx(40.0).epsilon(1e-12));
        CHECK(r.beta[1] == doctest::Approx(400.0).epsilon(1e-12));
        CHECK(r.Phi(0, 0) == 0.0);
        CHECK(r.L[0] == 0.0);
    }

    TEST_CASE("2nd-order binomial β") {
        constexpr auto r = design::eso<2>(20.0, 1.0);
        static_assert(r.success);
        CHECK(r.beta[0] == doctest::Approx(60.0).epsilon(1e-12));
        CHECK(r.beta[1] == doctest::Approx(1200.0).epsilon(1e-12));
        CHECK(r.beta[2] == doctest::Approx(8000.0).epsilon(1e-12));
    }

    TEST_CASE("Φ, Γ and L match the closed form") {
        const double h = 1e-3, wo = 50.0, b0 = 2.0;
        const double l = damp::exp(-wo * h), m = 1.0 - l;

        const auto r1 = design::eso<1>(wo, b0, h);
        REQUIRE(r1.success);
        CHECK(r1.Phi(0, 1) == doctest::Approx(h));
        CHECK(r1.Gamma[0] == doctest::Approx(b0 * h));
        CHECK(r1.Gamma[1] == 0.0);
        CHECK(r1.L[0] == doctest::Approx(1.0 - (l * l)));
        CHECK(r1.L[1] == doctest::Approx(m * m / h));

        const auto r2 = design::eso<2>(wo, b0, h);
        REQUIRE(r2.success);
        CHECK(r2.Phi(0, 2) == doctest::Approx(0.5 * h * h));
        CHECK(r2.Phi(1, 2) == doctest::Approx(h));
        CHECK(r2.Gamma[0] == doctest::Approx(0.5 * b0 * h * h));
        CHECK(r2.Gamma[1] == doctest::Approx(b0 * h));
        CHECK(r2.Gamma[2] == 0.0);
        CHECK(r2.L[0] == doctest::Approx(1.0 - (l * l * l)));
        CHECK(r2.L[1] == doctest::Approx(1.5 * m * m * (1.0 + l) / h));
        CHECK(r2.L[2] == doctest::Approx(m * m * m / (h * h)));
    }

    TEST_CASE("design::adrc copies the ESO matrices") {
        constexpr double wo = 50.0, b0 = 2.0, Ts = 1e-3;
        const auto       a = design::adrc<2>(10.0, wo, b0, Ts);
        const auto       e = design::eso<2>(wo, b0, Ts);
        REQUIRE(a.success);
        REQUIRE(e.success);
        CHECK(a.wo == doctest::Approx(e.wo));
        CHECK(a.b0 == doctest::Approx(e.b0));
        CHECK(a.Ts == doctest::Approx(e.Ts));
        for (size_t i = 0; i < 3; ++i) {
            CHECK(a.beta[i] == doctest::Approx(e.beta[i]).epsilon(1e-15));
            CHECK(a.Gamma[i] == doctest::Approx(e.Gamma[i]).epsilon(1e-15));
            CHECK(a.L[i] == doctest::Approx(e.L[i]).epsilon(1e-15));
            for (size_t j = 0; j < 3; ++j) {
                CHECK(a.Phi(i, j) == doctest::Approx(e.Phi(i, j)).epsilon(1e-15));
            }
        }
        const auto from_adrc = a.as_eso();
        CHECK(from_adrc.L[0] == doctest::Approx(e.L[0]).epsilon(1e-15));
        CHECK(from_adrc.success);
    }

    TEST_CASE("ADRCController is ESO plus PD") {
        const auto                res = design::adrc<2>(10.0, 50.0, 2.0, 1e-3);
        ADRCController<2, double> law{res};
        ESO<2, double>            obs{res.as_eso()};
        double                    y = 0.0, v = 0.0;
        const double              r = 1.0, Ts = 1e-3, b = 2.0, d = 0.3;
        for (int k = 0; k < 400; ++k) {
            const double u_law = law.control(r, y);
            const auto&  z = obs.step(y);
            const double u_obs = ((res.Kp * (r - z[0])) - (res.Kd * z[1]) - z[2]) / res.b0;
            obs.back_calculate(u_obs);
            CHECK(u_law == doctest::Approx(u_obs).epsilon(1e-14));
            CHECK(law.observer()[2] == doctest::Approx(obs.disturbance()).epsilon(1e-14));
            v += ((b * u_law) + d) * Ts;
            y += v * Ts;
        }
    }

    TEST_CASE("an exact ZOH plant with u ≠ 0: every state converges, for both orders") {
        const double h = 1e-3, b0 = 20.0, d = -3.5;

        SUBCASE("NX = 2") {
            ESO<2, double>   obs{design::eso<2>(100.0, b0, h)};
            DoubleIntegrator p{.y = 0.2, .v = -0.1};
            double           u = 0.0;
            for (int k = 0; k < 6000; ++k) {
                (void)obs.step(p.y, u);
                u = 0.1 + (0.3 * damp::sin(0.01 * k));
                p.step((b0 * u) + d, h);
            }
            (void)obs.step(p.y, u);
            CHECK(obs.output() == doctest::Approx(p.y).epsilon(1e-9));
            CHECK(obs.state()[1] == doctest::Approx(p.v).epsilon(1e-9));
            CHECK(obs.disturbance() == doctest::Approx(d).epsilon(1e-9));
        }

        SUBCASE("NX = 1") {
            ESO<1, double> obs{design::eso<1>(100.0, b0, h)};
            double         y = 0.2, u = 0.0;
            for (int k = 0; k < 6000; ++k) {
                (void)obs.step(y, u);
                u = 0.1 + (0.3 * damp::sin(0.01 * k));
                y += ((b0 * u) + d) * h;
            }
            (void)obs.step(y, u);
            CHECK(obs.output() == doctest::Approx(y).epsilon(1e-9));
            CHECK(obs.disturbance() == doctest::Approx(d).epsilon(1e-9));
        }
    }

    TEST_CASE("a plant-gain error lands in f̂") {
        const double     h = 1e-3, b0 = 20.0, b = 26.0, d = -3.5, u = 0.5;
        ESO<2, double>   obs{design::eso<2>(100.0, b0, h)};
        DoubleIntegrator p{};
        for (int k = 0; k < 6000; ++k) {
            (void)obs.step(p.y, u);
            p.step((b * u) + d, h);
        }
        CHECK(obs.disturbance() == doctest::Approx(d + ((b - b0) * u)).epsilon(1e-9));
    }

    TEST_CASE("predict coasts through missed samples") {
        const double     h = 1e-3, a = 2.0;
        ESO<2, double>   obs{design::eso<2>(100.0, 1.0, h)};
        DoubleIntegrator p{};
        for (int k = 0; k < 3000; ++k) {
            (void)obs.step(p.y, 0.0);
            p.step(a, h);
        }
        for (int k = 0; k < 50; ++k) { // 50 ms with no reading: the model carries on
            obs.predict(ColVec<1>{0.0});
            p.step(a, h);
        }
        obs.predict(ColVec<1>{0.0});
        obs.update(ColVec<1>{p.y});
        CHECK(obs.output() == doctest::Approx(p.y).epsilon(1e-6));
        CHECK(obs.state()[1] == doctest::Approx(p.v).epsilon(1e-6));
        CHECK(obs.disturbance() == doctest::Approx(a).epsilon(1e-6));
    }

    TEST_CASE("update corrects without advancing time") {
        ESO<2, double> obs{design::eso<2>(100.0, 1.0, 1e-3)};
        obs.reset(ColVec<3>{1.0, 5.0, 2.0});
        obs.update(ColVec<1>{1.0}); // zero innovation: nothing to correct, and no time passes
        CHECK(obs.state()[0] == 1.0);
        CHECK(obs.state()[1] == 5.0);
        CHECK(obs.state()[2] == 2.0);
    }

    TEST_CASE("the low-latency ISR order matches the fused step") {
        const auto     res = design::eso<2>(100.0, 2.0, 1e-3);
        ESO<2, double> fused{res};
        ESO<2, double> isr{res};
        double         y = 0.1, u = 0.0;
        for (int k = 0; k < 200; ++k) {
            (void)fused.step(y, u);
            isr.update(ColVec<1>{y});
            CHECK(isr.disturbance() == doctest::Approx(fused.disturbance()).epsilon(1e-12));
            u = 0.2 + (0.1 * damp::sin(0.05 * k));
            isr.back_calculate(u);
            isr.predict();
            y += 1e-3 * u;
        }
    }

    TEST_CASE("first tick is L y when z and u are zero") {
        const auto     res = design::eso<2>(40.0, 1.5, 1e-3);
        ESO<2, double> obs{res};
        const double   y = 0.4;
        const auto&    z = obs.step(y, 0.0);
        CHECK(z[0] == doctest::Approx(res.L[0] * y).epsilon(1e-14));
        CHECK(z[1] == doctest::Approx(res.L[1] * y).epsilon(1e-14));
        CHECK(z[2] == doctest::Approx(res.L[2] * y).epsilon(1e-14));
    }

    TEST_CASE("back_calculate then step(y) matches step(y, u)") {
        const auto     res = design::eso<1>(80.0, 1.0, 5e-4);
        ESO<1, double> a{res};
        ESO<1, double> b{res};
        const double   y = 0.25, u = 0.1;
        (void)a.step(y, u);
        b.back_calculate(u);
        (void)b.step(y);
        CHECK(a.state()[0] == doctest::Approx(b.state()[0]).epsilon(1e-15));
        CHECK(a.state()[1] == doctest::Approx(b.state()[1]).epsilon(1e-15));
        CHECK(a.last_command() == doctest::Approx(u));
        CHECK(b.last_command() == doctest::Approx(u));
    }

    TEST_CASE("predict then update equals fused step") {
        const auto     res = design::eso<2>(40.0, 1.0, 1e-3);
        ESO<2, double> fused{res};
        ESO<2, double> manual{res};
        for (int k = 0; k < 40; ++k) {
            const ColVec<1> y{0.05 * k};
            const ColVec<1> u{0.2};
            const auto&     zf = fused.step(y, u);
            manual.predict(u);
            manual.update(y);
            CHECK(zf[0] == doctest::Approx(manual.state()[0]).epsilon(1e-15));
            CHECK(zf[1] == doctest::Approx(manual.state()[1]).epsilon(1e-15));
            CHECK(zf[2] == doctest::Approx(manual.state()[2]).epsilon(1e-15));
        }
    }

    TEST_CASE("float at a large |y|: the innovation form leaves no DC bias in the estimates") {
        ESO<2, float> obs{design::eso<2>(40.0, 20.0, 1e-4)};
        const float   y = 1000.0f;
        obs.reset(ColVec<3, float>{y, 0.0f, 0.0f});
        double    v_sum = 0.0, f_sum = 0.0;
        const int N = 20000;
        for (int k = 0; k < N; ++k) {
            const auto& z = obs.step(y, 0.0f);
            v_sum += static_cast<double>(z[1]);
            f_sum += static_cast<double>(z[2]);
        }
        // True values are 0. F ẑ + G_y y gave means of about −0.14 and −0.34 here.
        CHECK(damp::abs(v_sum / N) < 0.02);
        CHECK(damp::abs(f_sum / N) < 0.05);
    }

    TEST_CASE("any finite b0: a reverse-acting plant and a plain differentiator") {
        const double h = 1e-3;

        SUBCASE("negative b0") {
            const double     b0 = -5.0, d = 1.5, u = 0.4;
            ESO<2, double>   obs{design::eso<2>(100.0, b0, h)};
            DoubleIntegrator p{};
            for (int k = 0; k < 6000; ++k) {
                (void)obs.step(p.y, u);
                p.step((b0 * u) + d, h);
            }
            CHECK(obs.disturbance() == doctest::Approx(d).epsilon(1e-9));
        }

        SUBCASE("b0 = 0 with no input differentiates y") {
            const double   slope = 3.0;
            ESO<2, double> obs{design::eso<2>(100.0, 0.0, h)};
            for (int k = 0; k < 6000; ++k) {
                (void)obs.step(slope * h * k, 0.0);
            }
            CHECK(obs.state()[1] == doctest::Approx(slope).epsilon(1e-9));
            CHECK(obs.disturbance() == doctest::Approx(0.0).scale(1.0).epsilon(1e-6));
        }
    }

    TEST_CASE("invalid knobs fail closed") {
        CHECK(design::eso<1>(0.0, 1.0).status == design::DesignStatus::invalid_parameter);
        CHECK(design::eso<2>(-10.0, 1.0).status == design::DesignStatus::invalid_parameter);
        CHECK(design::eso<1>(10.0, std::numeric_limits<double>::infinity()).status == design::DesignStatus::invalid_parameter);
        CHECK(design::eso<1>(10.0, std::numeric_limits<double>::quiet_NaN()).status == design::DesignStatus::invalid_parameter);
        CHECK(design::eso<2>(10.0, 1.0, -1e-3).status == design::DesignStatus::invalid_sample_time);
        CHECK(design::eso<1>(10.0, 1.0, std::numeric_limits<double>::quiet_NaN()).status == design::DesignStatus::invalid_sample_time);
        CHECK_FALSE(design::eso<1>(std::numeric_limits<double>::infinity(), 1.0).success);

        ESO<1, double> bad{design::eso<1>(0.0, 1.0, 1e-3)};
        CHECK_FALSE(bad.valid());
        CHECK(bad.step(1.0, 1.0)[0] == 0.0);
        CHECK(bad.disturbance() == 0.0);
    }

    TEST_CASE("gains that T cannot represent are refused, not reported as success") {
        const auto bad = design::DesignStatus::invalid_parameter;
        CHECK(design::eso<2>(1e150, 1.0).status == bad);               // ω_o³ overflows
        CHECK(design::eso<2, float>(1e13f, 1.0f).status == bad);       // ω_o³ overflows float
        CHECK(design::eso<2>(1.0, 1.0, 1e-200).status == bad);         // λ = 1, h² underflows
        CHECK(design::eso<2>(1e-20, 1.0, 1e-3).status == bad);         // λ = 1: L = 0 never corrects
        CHECK(design::eso<2, float>(1.0f, 1.0f, 1e-9f).status == bad); // λ = 1 in float
        CHECK(design::eso<2>(1e200, 1.0, 1e-200).status == bad);       // L₂ overflows through h²
        CHECK(design::eso<1>(1e6, 1.0, 1.0).success);                  // deadbeat is fine
    }

    TEST_CASE("a continuous design (Ts = 0) leaves the runtime invalid") {
        auto obs = ESO<2, double>(design::eso<2>(40.0, 1.0));
        CHECK_FALSE(obs.valid());
        CHECK(obs.step(1.0, 0.5)[0] == 0.0);
        CHECK(obs.output() == 0.0);
    }

    TEST_CASE("ESOResult::as and converting ctor") {
        constexpr auto ok = design::eso<1>(20.0, 1.5, 1e-3);
        static_assert(ok.success);
        constexpr auto ok_f = ok.as<float>();
        static_assert(ok_f.success);
        CHECK(ok_f.wo == doctest::Approx(static_cast<float>(ok.wo)));

        constexpr auto bad = design::eso<1>(-20.0, 1.0);
        static_assert(!bad.success);
        static_assert(!bad.as<float>().success);

        ESO<1, double> d{ok};
        (void)d.step(0.2, 0.1);
        (void)d.step(0.4, 0.1);
        ESO<1, float> f{d};
        CHECK(f.valid());
        const auto zd = d.step(0.5, 0.1);
        const auto zf = f.step(0.5f, 0.1f);
        CHECK(static_cast<double>(zf[0]) == doctest::Approx(zd[0]).epsilon(1e-5));
    }

    TEST_CASE("CTAD keeps the result's scalar; constinit deploys float") {
        ESO obs{kEso2};
        static_assert(std::is_same_v<decltype(obs), ESO<2, double>>);
        CHECK(kObs2.valid());
        CHECK(obs.valid());
        CHECK(static_cast<double>(kObs2.step(0.1f, 0.0f)[0]) == doctest::Approx(obs.step(0.1, 0.0)[0]).epsilon(1e-5));
    }

    TEST_CASE("reset clears the observer state") {
        const auto     res = design::eso<1>(40.0, 1.0, 1e-3);
        ESO<1, double> obs{res};
        for (int k = 0; k < 50; ++k) {
            (void)obs.step(0.3, 0.1);
        }
        obs.reset();
        CHECK(obs.state()[0] == 0.0);
        CHECK(obs.state()[1] == 0.0);
        CHECK(obs.last_command() == 0.0);
        ESO<1, double> fresh{res};
        CHECK(obs.step(0.3, 0.1)[0] == doctest::Approx(fresh.step(0.3, 0.1)[0]));
    }

    TEST_CASE("reset(z0) starts on the measurement with no transient") {
        const auto     res = design::eso<2>(50.0, 1.0, 1e-3);
        ESO<2, double> cold{res};
        ESO<2, double> warm{res};
        warm.reset(ColVec<3>{1.0, 0.0, 0.0});
        double cold_peak = 0.0, warm_peak = 0.0;
        for (int k = 0; k < 2000; ++k) {
            cold_peak = damp::max(cold_peak, damp::abs(cold.step(1.0, 0.0)[2]));
            warm_peak = damp::max(warm_peak, damp::abs(warm.step(1.0, 0.0)[2]));
        }
        CHECK(cold_peak > 100.0);
        CHECK(warm_peak < 1e-9);
    }

    TEST_CASE("the discrete ESO's poles are all at e^{−ω_o T_s}") {
        for (const double woTs : {1e-3, 0.1, 1.0, 3.0}) {
            const double Ts = 1e-3;
            const double lambda = damp::exp(-woTs);
            CAPTURE(woTs);
            const auto r1 = design::eso<1>(woTs / Ts, 2.0, Ts);
            const auto r2 = design::eso<2>(woTs / Ts, 2.0, Ts);
            const auto c1 = zpk_detail::leverrier(r1.error_dynamics()).den;
            const auto c2 = zpk_detail::leverrier(r2.error_dynamics()).den;
            const auto e1 = repeated_root_poly<2>(lambda);
            const auto e2 = repeated_root_poly<3>(lambda);
            for (size_t i = 0; i < 3; ++i) {
                CHECK(c1[i] == doctest::Approx(e1[i]).epsilon(1e-12).scale(1.0));
            }
            for (size_t i = 0; i < 4; ++i) {
                CHECK(c2[i] == doctest::Approx(e2[i]).epsilon(1e-12).scale(1.0));
            }
        }
        constexpr auto rc = design::eso<2>(10.0, 50.0, 0.0);
        static_assert(rc.Phi(0, 0) == 0.0 && rc.L(0) == 0.0);
        constexpr auto rd = design::eso<2>(10.0, 50.0, 1e-3);
        static_assert(rd.L(0) > 0.0 && rd.L(0) < 1.0);
    }
}
