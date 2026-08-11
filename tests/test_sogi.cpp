// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <numbers>

#include "damp/backend.hpp"
#include "damp/filters/sogi.hpp"
#include "damp/math/math.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/discretization.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

TEST_CASE("SOGI design") {
    constexpr double omega_0 = 2 * std::numbers::pi * 50.0; // 50 Hz
    constexpr double alpha = 1.414;
    constexpr auto   sogi_sys = design::sogi<double>(omega_0, alpha);

    // Check matrix dimensions
    CHECK(sogi_sys.A.rows() == 2);
    CHECK(sogi_sys.A.cols() == 2);
    CHECK(sogi_sys.B.rows() == 2);
    CHECK(sogi_sys.B.cols() == 1);
    CHECK(sogi_sys.C.rows() == 2);
    CHECK(sogi_sys.C.cols() == 2);

    // Check A matrix structure
    CHECK(sogi_sys.A(0, 0) == doctest::Approx(-alpha * omega_0));
    CHECK(sogi_sys.A(0, 1) == doctest::Approx(-omega_0));
    CHECK(sogi_sys.A(1, 0) == doctest::Approx(omega_0));
    CHECK(sogi_sys.A(1, 1) == doctest::Approx(0.0));

    // Check B matrix
    CHECK(sogi_sys.B(0, 0) == doctest::Approx(alpha * omega_0));
    CHECK(sogi_sys.B(1, 0) == doctest::Approx(0.0));

    // Check C matrix (identity for direct state output)
    CHECK(sogi_sys.C(0, 0) == doctest::Approx(1.0));
    CHECK(sogi_sys.C(0, 1) == doctest::Approx(0.0));
    CHECK(sogi_sys.C(1, 0) == doctest::Approx(0.0));
    CHECK(sogi_sys.C(1, 1) == doctest::Approx(1.0));
}

TEST_CASE("SOGI discrete design is the ZOH of the continuous SOGI and band-passes at w0") {
    constexpr double w0 = 2 * std::numbers::pi * 50.0;
    constexpr double alpha = 1.414;
    constexpr double Ts = 1.0 / 10000.0;

    constexpr auto sys = design::sogi<double>(w0, alpha, Ts);

    // The discrete overload is exactly discretize(continuous, Ts, ZOH).
    constexpr auto zoh = *discretize(design::sogi<double>(w0, alpha), Ts, DiscretizationMethod::ZOH);
    CHECK(sys.A(0, 0) == doctest::Approx(zoh.A(0, 0)));
    CHECK(sys.A(1, 1) == doctest::Approx(zoh.A(1, 1)));
    CHECK(sys.B(0, 0) == doctest::Approx(zoh.B(0, 0)));

    // Step the SS with a 50 Hz tone; the band-pass state (y0) reaches unity gain.
    damp::ColVec<2, double> x = {};
    double                  max_bp = 0.0;
    for (int i = 0; i < 5000; ++i) {
        const double                  in = std::sin(w0 * Ts * static_cast<double>(i));
        const damp::ColVec<2, double> y = sys.C * x;
        x = (sys.A * x) + (sys.B * in);
        if (i > 2500) {
            max_bp = damp::max(max_bp, std::abs(y(0)));
        }
    }
    CHECK(max_bp == doctest::Approx(1.0).epsilon(0.05));
}

TEST_CASE("MSTOGI design - TOGI structure with DC-rejecting quadrature") {
    constexpr double omega_0 = 2 * std::numbers::pi * 50.0;
    constexpr double alpha = 1.414;
    constexpr auto   mstogi_sys = design::mstogi<double>(omega_0, alpha);

    // Check dimensions
    CHECK(mstogi_sys.A.rows() == 3);
    CHECK(mstogi_sys.A.cols() == 3);
    CHECK(mstogi_sys.B.rows() == 3);
    CHECK(mstogi_sys.B.cols() == 1);
    CHECK(mstogi_sys.C.rows() == 2);
    CHECK(mstogi_sys.C.cols() == 3);

    // SOGI core (rows 0,1): in-phase + quadrature integrators.
    CHECK(mstogi_sys.A(0, 0) == doctest::Approx(-alpha * omega_0));
    CHECK(mstogi_sys.A(0, 1) == doctest::Approx(-omega_0));
    CHECK(mstogi_sys.A(1, 0) == doctest::Approx(omega_0));
    CHECK(mstogi_sys.A(1, 1) == doctest::Approx(0.0));

    // TOGI stage (row 2): v̇‴ = alpha·ω₀·(v − v′) − ω₀·v‴ — first-order, self-damped,
    // fed from the same post-gain bus as the SOGI (NOT the old −ω₀²·x₂ double
    // integrator, which made neither a washout nor a unity quadrature pair).
    CHECK(mstogi_sys.A(2, 0) == doctest::Approx(-alpha * omega_0));
    CHECK(mstogi_sys.A(2, 1) == doctest::Approx(0.0));
    CHECK(mstogi_sys.A(2, 2) == doctest::Approx(-omega_0));
    CHECK(mstogi_sys.B(2, 0) == doctest::Approx(alpha * omega_0));

    // Output map: v_o = v′ (band-pass), q·v_o = v″ − v‴ (DC-rejecting quadrature).
    CHECK(mstogi_sys.C(0, 0) == doctest::Approx(1.0));
    CHECK(mstogi_sys.C(1, 1) == doctest::Approx(1.0));
    CHECK(mstogi_sys.C(1, 2) == doctest::Approx(-1.0));
}

TEST_CASE("MSTOGI discrete design is the ZOH of the continuous MSTOGI and band-passes at w0") {
    constexpr double w0 = 2 * std::numbers::pi * 50.0;
    constexpr double alpha = 1.414;
    constexpr double Ts = 1.0 / 10000.0;

    constexpr auto sys = design::mstogi<double>(w0, alpha, Ts);

    // The discrete overload is exactly discretize(continuous, Ts, ZOH).
    constexpr auto zoh = *discretize(design::mstogi<double>(w0, alpha), Ts, DiscretizationMethod::ZOH);
    CHECK(sys.A(0, 0) == doctest::Approx(zoh.A(0, 0)));
    CHECK(sys.A(2, 2) == doctest::Approx(zoh.A(2, 2)));
    CHECK(sys.B(0, 0) == doctest::Approx(zoh.B(0, 0)));

    // C/D output map carried through from the continuous design: v_o = v′ (band-pass),
    // q·v_o = v″ − v‴ (DC-rejecting quadrature).
    CHECK(sys.C(0, 0) == doctest::Approx(1.0));
    CHECK(sys.C(1, 1) == doctest::Approx(1.0));
    CHECK(sys.C(1, 2) == doctest::Approx(-1.0));

    // Step the SS with a 50 Hz tone; the band-pass output (y0) reaches unity gain.
    damp::ColVec<3, double> x = {};
    double                  max_bp = 0.0;
    for (int i = 0; i < 5000; ++i) {
        const double                  in = std::sin(w0 * Ts * static_cast<double>(i));
        const damp::ColVec<2, double> y = sys.C * x;
        x = (sys.A * x) + (sys.B * in);
        if (i > 2500) {
            max_bp = damp::max(max_bp, std::abs(y(0)));
        }
    }
    CHECK(max_bp == doctest::Approx(1.0).epsilon(0.05));
}

TEST_CASE("MSTOGI runtime - rejects DC offset on quadrature, tracks fundamental") {
    constexpr float f0 = 50.0f;
    constexpr float Ts = 1.0f / 10000.0f; // 10 kHz
    constexpr float alpha = 1.414f;
    MSTOGI<float>   mstogi;

    const float omega = 2.0f * std::numbers::pi_v<float> * f0;
    const float dc_offset = 0.5f; // bias that a plain SOGI-QSG would leak into qv′

    float bp_amp = 0.0f;
    float q_amp = 0.0f;
    float q_mean = 0.0f;
    int   n_mean = 0;

    // Run ~30 cycles; measure steady-state behaviour over the last 10.
    const int n = 6000;
    for (int i = 0; i < n; ++i) {
        const float t = static_cast<float>(i) * Ts;
        const float v = dc_offset + std::sin(omega * t);
        const auto [bp, q] = mstogi(v, f0, alpha, Ts);

        if (i > n - 2000) {
            bp_amp = damp::max(bp_amp, std::abs(bp));
            q_amp = damp::max(q_amp, std::abs(q));
            q_mean += q;
            ++n_mean;
        }
    }
    q_mean /= static_cast<float>(n_mean);

    // Band-pass tracks the fundamental (unity gain at ω₀), and the quadrature
    // channel has comparable amplitude...
    CHECK(bp_amp == doctest::Approx(1.0f).epsilon(0.05f));
    CHECK(q_amp == doctest::Approx(1.0f).epsilon(0.05f));
    // ...but its DC mean is driven to ~0 by the TOGI even though the input has a
    // 0.5 DC bias. That is the whole point of the MSTOGI vs a plain SOGI-QSG.
    CHECK(std::abs(q_mean) < 0.02f);
}

TEST_CASE("MSTOGI runtime - frequency retune (center frequency moves mid-run)") {
    // The motivating use case: the center frequency is changed every tick (e.g.
    // an inverter PLL feeding its estimated frequency into the band-pass). Lock
    // at 50 Hz, step the commanded frequency to 60 Hz, and confirm the band-pass
    // re-locks to the new input — the per-call (in, freq, alpha, Ts) interface
    // rebuilds the resonator from sincos(w*Ts) each sample, so no rediscretize.
    constexpr float Ts = 1.0f / 10000.0f;
    constexpr float alpha = 1.414f;
    MSTOGI<float>   mstogi;

    // Phase-continuous input that switches frequency at the retune instant, so
    // the only thing under test is the filter retune, not an input phase jump.
    float       phase = 0.0f;
    const float two_pi = 2.0f * std::numbers::pi_v<float>;

    // Settle at 50 Hz.
    for (int i = 0; i < 2000; ++i) {
        phase += two_pi * 50.0f * Ts;
        (void)mstogi(std::sin(phase), 50.0f, alpha, Ts);
    }

    // Step the commanded (and input) frequency to 60 Hz; measure after re-lock.
    float bp_amp = 0.0f;
    float q_amp = 0.0f;
    for (int i = 0; i < 4000; ++i) {
        phase += two_pi * 60.0f * Ts;
        const auto [bp, q] = mstogi(std::sin(phase), 60.0f, alpha, Ts);
        if (i > 2000) {
            bp_amp = damp::max(bp_amp, std::abs(bp));
            q_amp = damp::max(q_amp, std::abs(q));
        }
    }
    CHECK(bp_amp == doctest::Approx(1.0f).epsilon(0.05f)); // re-locked to 60 Hz
    CHECK(q_amp == doctest::Approx(1.0f).epsilon(0.05f));
}

TEST_CASE("SOGI runtime - resonator form") {
    constexpr float f0 = 50.0f;
    constexpr float Ts = 1.0f / 10000.0f;
    constexpr float alpha = 1.414f;
    SOGI<float>     sogi;

    const float omega_input = 2.0f * std::numbers::pi_v<float> * f0;

    float max_bp = 0.0f;
    float max_quad = 0.0f;

    for (int i = 0; i < 5000; ++i) {
        const float t = static_cast<float>(i) * Ts;
        const float input = std::sin(omega_input * t);

        const auto [bp, quad] = sogi(input, f0, alpha, Ts);

        if (i > 2500) {
            max_bp = damp::max(max_bp, std::abs(bp));
            max_quad = damp::max(max_quad, std::abs(quad));
        }
    }

    CHECK(max_bp == doctest::Approx(1.0f).epsilon(0.05f));
    CHECK(max_quad == doctest::Approx(1.0f).epsilon(0.05f));
}

TEST_CASE("SOGI runtime - frequency retuning") {
    constexpr float Ts = 1.0f / 10000.0f;
    constexpr float alpha = 1.414f;
    SOGI<float>     sogi;
    const float     omega_50 = 2.0f * std::numbers::pi_v<float> * 50.0f;
    for (int i = 0; i < 2000; ++i) {
        const float t = static_cast<float>(i) * Ts;
        (void)sogi(std::sin(omega_50 * t), 50.0f, alpha, Ts);
    }

    const float omega_60 = 2.0f * std::numbers::pi_v<float> * 60.0f;
    float       max_bp = 0.0f;
    for (int i = 0; i < 4000; ++i) {
        const float t = static_cast<float>(i) * Ts;
        const auto [bp, q] = sogi(std::sin(omega_60 * t), 60.0f, alpha, Ts);
        (void)q;
        if (i > 2000) {
            max_bp = damp::max(max_bp, std::abs(bp));
        }
    }

    CHECK(max_bp == doctest::Approx(1.0f).epsilon(0.05f));
}

TEST_CASE("SOGI runtime exposes derived notch output") {
    constexpr float f0 = 50.0f;
    constexpr float Ts = 1.0f / 10000.0f;
    constexpr float alpha = 1.414f;
    SOGI<float>     sogi;

    const float omega = 2.0f * std::numbers::pi_v<float> * f0;
    float       max_notch = 0.0f;
    for (int i = 0; i < 5000; ++i) {
        const float t = static_cast<float>(i) * Ts;
        const float input = std::sin(omega * t);
        const auto [bp, q] = sogi(input, f0, alpha, Ts);
        (void)q;

        if (i > 2500) {
            const float notch = input - bp;
            max_notch = damp::max(max_notch, std::abs(notch));
        }
    }

    CHECK(max_notch < 0.15f);
}

TEST_CASE("SOGI runtime matches resonator realization") {
    constexpr float Ts = 1.0f / 10000.0f;
    constexpr float f0 = 50.0f;
    constexpr float alpha = 1.414f;

    SOGI<float> sogi;
    const float w0 = 2.0f * std::numbers::pi_v<float> * f0;
    const float wT = w0 * Ts;
    const auto [sin_wt, cos_wt] = damp::sincos(wT);

    float x_r1 = 0.0f; // quadrature state
    float x_r2 = 0.0f; // band-pass state

    for (int i = 0; i < 32; ++i) {
        const float in = std::sin(2.0f * std::numbers::pi_v<float> * f0 * Ts * static_cast<float>(i));

        const float bp_ref = x_r2;
        const float q_ref = x_r1;
        const float u = alpha * (in - bp_ref);

        const auto [bp, q] = sogi(in, f0, alpha, Ts);

        CHECK(bp == doctest::Approx(bp_ref).epsilon(1e-6f));
        CHECK(q == doctest::Approx(q_ref).epsilon(1e-6f));

        const float x_r1_next = (cos_wt * x_r1) + (sin_wt * x_r2) + ((1.0f - cos_wt) * u);
        const float x_r2_next = (-sin_wt * x_r1) + (cos_wt * x_r2) + (sin_wt * u);

        x_r1 = x_r1_next;
        x_r2 = x_r2_next;
    }
}

TEST_CASE("MSTOGI runtime wrapper matches FE washout realization") {
    constexpr float Ts = 1.0f / 10000.0f;
    constexpr float f0 = 50.0f;
    constexpr float alpha = 1.414f;

    MSTOGI<float> mstogi;
    const float   w0 = 2.0f * std::numbers::pi_v<float> * f0;
    const float   wT = w0 * Ts;
    const auto [sin_wt, cos_wt] = damp::sincos(wT);

    float x_r1 = 0.0f;
    float x_r2 = 0.0f;
    float x_t = 0.0f;

    for (int i = 0; i < 32; ++i) {
        const float in = std::sin(2.0f * std::numbers::pi_v<float> * f0 * Ts * static_cast<float>(i));

        const float bp_ref = x_r2;
        const float u = alpha * (in - bp_ref);
        x_t += (u - x_t) * wT;
        const float q_ref = x_r1 - x_t;

        const auto [bp, q] = mstogi(in, f0, alpha, Ts);

        CHECK(bp == doctest::Approx(bp_ref).epsilon(1e-6f));
        CHECK(q == doctest::Approx(q_ref).epsilon(1e-6f));

        const float x_r1_next = (cos_wt * x_r1) + (sin_wt * x_r2) + ((1.0f - cos_wt) * u);
        const float x_r2_next = (-sin_wt * x_r1) + (cos_wt * x_r2) + (sin_wt * u);

        x_r1 = x_r1_next;
        x_r2 = x_r2_next;
    }
}

TEST_CASE("SOGI runtime - simple per-call frequency interface") {
    constexpr float f0 = 50.0f;
    constexpr float Ts = 1.0f / 10000.0f;
    constexpr float alpha = 1.414f;
    SOGI<float>     sogi;

    float       max_bp = 0.0f;
    float       max_q = 0.0f;
    const float omega = 2.0f * std::numbers::pi_v<float> * f0;

    for (int i = 0; i < 5000; ++i) {
        const float t = static_cast<float>(i) * Ts;
        const float in = std::sin(omega * t);
        const auto [bp, q] = sogi(in, f0, alpha, Ts);
        if (i > 2500) {
            max_bp = damp::max(max_bp, std::abs(bp));
            max_q = damp::max(max_q, std::abs(q));
        }
    }

    CHECK(max_bp == doctest::Approx(1.0f).epsilon(0.05f));
    CHECK(max_q == doctest::Approx(1.0f).epsilon(0.05f));
}

TEST_CASE("MSTOGI runtime - simple per-call frequency interface") {
    constexpr float f0 = 50.0f;
    constexpr float Ts = 1.0f / 10000.0f;
    constexpr float alpha = 1.414f;
    MSTOGI<float>   mstogi;

    const float omega = 2.0f * std::numbers::pi_v<float> * f0;
    float       bp_amp = 0.0f;
    float       q_amp = 0.0f;

    for (int i = 0; i < 5000; ++i) {
        const float t = static_cast<float>(i) * Ts;
        const float v = std::sin(omega * t);
        const auto [bp, q] = mstogi(v, f0, alpha, Ts);
        if (i > 2500) {
            bp_amp = damp::max(bp_amp, std::abs(bp));
            q_amp = damp::max(q_amp, std::abs(q));
        }
    }

    CHECK(bp_amp == doctest::Approx(1.0f).epsilon(0.05f));
    CHECK(q_amp == doctest::Approx(1.0f).epsilon(0.05f));
}

TEST_SUITE("sogi_fll") {
    TEST_CASE("locks onto a tone offset from the initial guess") {
        const double    fs = 2000.0, dt = 1.0 / fs;
        const double    f_true = 57.0;
        SogiFll<double> fll(50.0, dt); // start 7 Hz low

        for (int k = 0; k < 6000; ++k) { // 3 s
            fll.update(std::sin(2.0 * std::numbers::pi * f_true * (k * dt)));
        }
        CHECK(fll.frequency_hz() == doctest::Approx(f_true).epsilon(0.01)); // locked within 1%
    }

    TEST_CASE("recovers the tone amplitude after lock") {
        const double fs = 2000.0, dt = 1.0 / fs;
        const double f_true = 45.0, A = 2.5;

        SogiFll<double> fll(48.0, dt);
        for (int k = 0; k < 8000; ++k) {
            fll.update(A * std::sin(2.0 * std::numbers::pi * f_true * (k * dt)));
        }
        CHECK(fll.frequency_hz() == doctest::Approx(f_true).epsilon(0.02));
        CHECK(fll.amplitude() == doctest::Approx(A).epsilon(0.05));
    }

    TEST_CASE("tracks a frequency step") {
        const double    fs = 2000.0, dt = 1.0 / fs;
        SogiFll<double> fll(50.0, dt);

        // Phase-continuous tone that steps 50 -> 62 Hz halfway.
        double phase = 0.0;
        double f_final = 0.0;
        for (int k = 0; k < 12000; ++k) {
            const double f = (k < 6000) ? 50.0 : 62.0;
            phase += 2.0 * std::numbers::pi * f * dt;
            fll.update(std::sin(phase));
            f_final = f;
        }
        CHECK(f_final == 62.0);
        CHECK(fll.frequency_hz() == doctest::Approx(62.0).epsilon(0.02)); // re-locked
    }

    TEST_CASE("locks through broadband noise (it is a narrow band-pass)") {
        const double    fs = 2000.0, dt = 1.0 / fs;
        const double    f_true = 53.0;
        SogiFll<double> fll(50.0, dt);
        auto            noise = [](double t) {
            return 0.3 * (std::sin(617.0 * t) + std::sin(1303.0 * t) + std::sin(2999.0 * t)) / 3.0;
        };
        for (int k = 0; k < 8000; ++k) {
            const double t = k * dt;
            fll.update(std::sin(2.0 * std::numbers::pi * f_true * t) + noise(t));
        }
        CHECK(fll.frequency_hz() == doctest::Approx(f_true).epsilon(0.02));
    }

    TEST_CASE("clamps the tracked frequency to the configured band") {
        const double fs = 2000.0, dt = 1.0 / fs;
        // Drive at 90 Hz but clamp the search to [40, 70] Hz.
        SogiFll<double> fll(55.0, dt, damp::numbers::sqrt2_v<double>, 2.0, 40.0, 70.0);
        for (int k = 0; k < 6000; ++k) {
            fll.update(std::sin(2.0 * std::numbers::pi * 90.0 * (k * dt)));
        }
        CHECK(fll.frequency_hz() <= 70.0 + 1e-6);
        CHECK(fll.frequency_hz() >= 40.0 - 1e-6);
    }

    TEST_CASE("invalid config is inert; reset re-seeds the guess") {
        SogiFll<double> bad(50.0, -1.0); // Ts <= 0
        CHECK_FALSE(bad.valid());
        bad.update(1.0);
        CHECK(bad.frequency_hz() == doctest::Approx(50.0)); // unchanged

        SogiFll<double> fll(50.0, 1.0 / 2000.0);
        for (int k = 0; k < 3000; ++k) {
            fll.update(std::sin(2.0 * std::numbers::pi * 58.0 * (k / 2000.0)));
        }
        fll.reset(50.0);
        CHECK(fll.frequency_hz() == doctest::Approx(50.0));
        CHECK(fll.amplitude() == doctest::Approx(0.0));
    }

    TEST_CASE("float specialization locks") {
        SogiFll<float> fll(50.0f, 1.0f / 2000.0f);
        for (int k = 0; k < 8000; ++k) {
            fll.update(std::sin(2.0f * 3.14159265f * 56.0f * (static_cast<float>(k) / 2000.0f)));
        }
        CHECK(fll.frequency_hz() == doctest::Approx(56.0f).epsilon(0.03));
    }

    TEST_CASE("SOGI-FLL is constexpr-evaluable") {
        constexpr double f = []() consteval {
            SogiFll<double> fll(50.0, 1.0 / 2000.0);
            for (int k = 0; k < 6000; ++k) {
                fll.update(damp::sin(2.0 * 3.14159265358979323846 * 55.0 * (k / 2000.0)));
            }
            return fll.frequency_hz();
        }();
        static_assert(f > 54.0 && f < 56.0, "SOGI-FLL must lock at compile time");
        CHECK(f == doctest::Approx(55.0).epsilon(0.02));
    }
}
