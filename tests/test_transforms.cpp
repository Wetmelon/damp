// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <cstddef>
#include <numbers>

#include "damp/backend.hpp"
#include "damp/math/complex.hpp"
#include "damp/math/math.hpp"
#include "damp/math/transforms.hpp"
#include "damp/matrix/colvec.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

/**
 * @brief Tests for reference-frame transforms (Clarke/Park, zero sequence, Fortescue)
 */

TEST_SUITE("Transforms") {
    TEST_CASE("Clarke transform") {
        const ColVec<3, float> abc = {1.0f, -0.5f, -0.5f};
        const auto [alpha, beta] = clarke_transform(abc);

        CHECK(alpha == doctest::Approx(1.0f));
        CHECK(beta == doctest::Approx(0.0f).epsilon(1e-6f));
    }

    TEST_CASE("Inverse Clarke transform (three-wire)") {
        const ColVec<3, float> abc_orig = {1.0f, -0.5f, -0.5f};
        const auto             ab = clarke_transform(abc_orig);
        const auto             abc = inverse_clarke_transform(ab);

        CHECK(abc[0] == doctest::Approx(abc_orig[0]).epsilon(1e-6f));
        CHECK(abc[1] == doctest::Approx(abc_orig[1]).epsilon(1e-6f));
        CHECK(abc[2] == doctest::Approx(abc_orig[2]).epsilon(1e-6f));
    }

    TEST_CASE("Park transform identity at theta=0") {
        const AlphaBeta<float> ab = {.alpha = 1.0f, .beta = 0.5f};
        const auto [d, q] = park_transform(ab, 0.0f);

        CHECK(d == doctest::Approx(ab.alpha));
        CHECK(q == doctest::Approx(ab.beta));
    }

    TEST_CASE("Park d-axis convention: space vector at theta maps to (|v|, 0)") {
        const float            theta = std::numbers::pi_v<float> / 3.0f;
        const float            mag = 2.0f;
        const AlphaBeta<float> ab = {
            .alpha = mag * std::cos(theta),
            .beta = mag * std::sin(theta),
        };
        const auto [d, q] = park_transform(ab, theta);
        CHECK(d == doctest::Approx(mag).epsilon(1e-6f));
        CHECK(q == doctest::Approx(0.0f).epsilon(1e-6f));
    }

    TEST_CASE("Park round-trip") {
        const DirectQuadrature<float> dq_orig = {.d = 0.8f, .q = 0.3f};
        const float                   theta = std::numbers::pi_v<float> / 6.0f;
        const auto                    ab = inverse_park_transform(dq_orig, theta);
        const auto [d, q] = park_transform(ab, theta);

        CHECK(d == doctest::Approx(dq_orig.d).epsilon(1e-6f));
        CHECK(q == doctest::Approx(dq_orig.q).epsilon(1e-6f));
    }

    TEST_CASE("Fused Clarke-Park round-trip") {
        const DirectQuadrature<float> dq = {.d = 1.0f, .q = 0.5f};
        const float                   theta = std::numbers::pi_v<float> / 4.0f;
        const auto                    abc = inverse_park_clarke_transform(dq, theta);
        const auto [d2, q2] = clarke_park_transform(abc, theta);

        CHECK(d2 == doctest::Approx(dq.d).epsilon(1e-6f));
        CHECK(q2 == doctest::Approx(dq.q).epsilon(1e-6f));
    }

    TEST_CASE("Zero sequence: common-mode DC lands in zero, not αβ") {
        const float            d = 0.7f;
        const ColVec<3, float> abc = {d, d, d};
        const auto             ab = clarke_transform<float>(abc);

        CHECK(ab.alpha == doctest::Approx(0.0f).epsilon(1e-6f));
        CHECK(ab.beta == doctest::Approx(0.0f).epsilon(1e-6f));
        CHECK(zero_sequence<float>(abc) == doctest::Approx(d).epsilon(1e-6f));
    }

    TEST_CASE("Clarke round-trip with scalar zero") {
        const ColVec<3, float> abc = {1.3f, -0.4f, 0.2f};
        const auto             ab = clarke_transform<float>(abc);
        const auto             z = zero_sequence<float>(abc);
        const auto             rt = inverse_clarke_transform(ab, z);

        for (size_t i = 0; i < 3; ++i) {
            CHECK(rt[i] == doctest::Approx(abc[i]).epsilon(1e-6f));
        }
    }

    TEST_CASE("Fused Park-Clarke round-trip with scalar zero passthrough") {
        const ColVec<3, float> abc = {1.3f, -0.4f, 0.2f};
        const float            theta = 0.9f;
        const auto             dq = clarke_park_transform<float>(abc, theta);
        const auto             z = zero_sequence<float>(abc);
        const auto             rt = inverse_park_clarke_transform(dq, theta, z);

        for (size_t i = 0; i < 3; ++i) {
            CHECK(rt[i] == doctest::Approx(abc[i]).epsilon(1e-5f));
        }
        CHECK(z == doctest::Approx(zero_sequence<float>(rt)).epsilon(1e-5f));
    }

    TEST_CASE("Power-invariant Clarke scaling and round-trip") {
        using damp::Convention;
        const ColVec<3, float> abc = {1.0f, -0.4f, 0.3f};
        const auto             amp = clarke_transform<float>(abc);
        const auto             pwr = clarke_transform<float, Convention::PowerInvariant>(abc);
        const float            ratio = std::sqrt(1.5f);

        CHECK(pwr.alpha == doctest::Approx(amp.alpha * ratio).epsilon(1e-6f));
        CHECK(pwr.beta == doctest::Approx(amp.beta * ratio).epsilon(1e-6f));

        const auto z = zero_sequence<float, Convention::PowerInvariant>(abc);
        const auto rt = inverse_clarke_transform(pwr, z);
        for (size_t i = 0; i < 3; ++i) {
            CHECK(rt[i] == doctest::Approx(abc[i]).epsilon(1e-6f));
        }
    }

    TEST_CASE("Symmetrical components: balanced set is pure positive sequence") {
        using Cplx = damp::complex<float>;
        const Cplx a = {1.0f, 0.0f};
        const Cplx b = {-0.5f, -std::numbers::sqrt3_v<float> / 2.0f};
        const Cplx c = {-0.5f, std::numbers::sqrt3_v<float> / 2.0f};
        const auto seq = symmetrical_components<float>({a, b, c});

        CHECK(seq.zero.abs() == doctest::Approx(0.0f).epsilon(1e-6f));
        CHECK(seq.negative.abs() == doctest::Approx(0.0f).epsilon(1e-6f));
        CHECK(seq.positive.real() == doctest::Approx(1.0f).epsilon(1e-6f));
        CHECK(seq.positive.imag() == doctest::Approx(0.0f).epsilon(1e-6f));
    }

    TEST_CASE("Instantaneous power: balanced load, αβ and dq agree") {
        const float            tp = 2.0f * std::numbers::pi_v<float> / 3.0f;
        const float            th = 0.37f;
        const ColVec<3, float> v = {std::cos(th), std::cos(th - tp), std::cos(th + tp)};
        const ColVec<3, float> i = v;

        const auto p_ab = instantaneous_power<float>(clarke_transform<float>(v), clarke_transform<float>(i));
        CHECK(p_ab.p == doctest::Approx(1.5f).epsilon(1e-5f));
        CHECK(p_ab.q == doctest::Approx(0.0f).epsilon(1e-5f));

        const auto p_dq =
            instantaneous_power<float>(clarke_park_transform<float>(v, th), clarke_park_transform<float>(i, th));
        CHECK(p_dq.p == doctest::Approx(p_ab.p).epsilon(1e-5f));
        CHECK(p_dq.q == doctest::Approx(p_ab.q).epsilon(1e-5f));
    }
}
