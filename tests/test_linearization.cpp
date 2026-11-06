// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <limits>

#include "damp/design/linearization.hpp"
#include "damp/matlab.hpp"
#include "damp/matrix/matrix.hpp"
#include "damp/systems/state_space.hpp"

#define DOCTEST_CONFIG_INCLUDE_TYPE_TRAITS
#include "doctest.h"

using namespace damp;

TEST_SUITE("Linearization") {
    TEST_CASE("linearize computes A B C D at operating point") {
        auto dynamics = [](const ColVec<2>& x, const ColVec<1>& u) -> ColVec<2> {
            return ColVec<2>{
                x(1, 0),
                std::sin(x(0, 0)) + (2.0 * u(0, 0))
            };
        };

        auto output = [](const ColVec<2>& x, const ColVec<1>& u) -> ColVec<1> {
            return ColVec<1>{x(0, 0) + (3.0 * u(0, 0))};
        };

        const ColVec<2> x_op{0.0, 0.0};
        const ColVec<1> u_op{0.0};

        const auto lin = design::linearize<2, 1, 1>(dynamics, output, x_op, u_op, 1e-6);

        CHECK(lin.A(0, 0) == doctest::Approx(0.0).epsilon(1e-6));
        CHECK(lin.A(0, 1) == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(lin.A(1, 0) == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(lin.A(1, 1) == doctest::Approx(0.0).epsilon(1e-6));

        CHECK(lin.B(0, 0) == doctest::Approx(0.0).epsilon(1e-6));
        CHECK(lin.B(1, 0) == doctest::Approx(2.0).epsilon(1e-6));

        CHECK(lin.C(0, 0) == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(lin.C(0, 1) == doctest::Approx(0.0).epsilon(1e-6));
        CHECK(lin.D(0, 0) == doctest::Approx(3.0).epsilon(1e-6));
    }

    TEST_CASE("linearize captures nonzero-point Jacobian") {
        auto dynamics = [](const ColVec<2>& x, const ColVec<1>& u) -> ColVec<2> {
            return ColVec<2>{
                x(1, 0),
                std::sin(x(0, 0)) + u(0, 0)
            };
        };

        auto output = [](const ColVec<2>& x, const ColVec<1>& /*u*/) -> ColVec<1> {
            return ColVec<1>{x(0, 0)};
        };

        const ColVec<2> x_op{0.3, 0.0};
        const ColVec<1> u_op{0.0};

        const auto lin = design::linearize<2, 1, 1>(dynamics, output, x_op, u_op, 1e-6);

        CHECK(lin.A(1, 0) == doctest::Approx(std::cos(0.3)).epsilon(1e-6));
        CHECK(lin.C(0, 0) == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(lin.D(0, 0) == doctest::Approx(0.0).epsilon(1e-6));
    }

    TEST_CASE("identity-output overload and to_state_space") {
        auto dynamics = [](const ColVec<2>& x, const ColVec<1>& u) -> ColVec<2> {
            return ColVec<2>{x(1, 0), x(0, 0) + u(0, 0)};
        };

        const ColVec<2> x_op{0.0, 0.0};
        const ColVec<1> u_op{0.0};

        const auto lin = design::linearize<2, 1>(dynamics, x_op, u_op, 1e-6);
        const auto sys = lin.to_state_space();

        CHECK(lin.C(0, 0) == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(lin.C(1, 1) == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(lin.D(0, 0) == doctest::Approx(0.0).epsilon(1e-6));
        CHECK(lin.D(1, 0) == doctest::Approx(0.0).epsilon(1e-6));

        CHECK(sys.A(0, 1) == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(sys.B(1, 0) == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(sys.Ts == doctest::Approx(0.0));
    }

    TEST_CASE("matlab linmod wrapper returns state-space") {
        auto dynamics = [](const ColVec<2>& x, const ColVec<1>& u) -> ColVec<2> {
            return ColVec<2>{x(1, 0), x(0, 0) + (2.0 * u(0, 0))};
        };

        auto output = [](const ColVec<2>& x, const ColVec<1>& u) -> ColVec<1> {
            return ColVec<1>{x(0, 0) + u(0, 0)};
        };

        const ColVec<2> x_op{0.0, 0.0};
        const ColVec<1> u_op{0.0};

        const auto sys = matlab::linmod<2, 1, 1>(dynamics, output, x_op, u_op, 1e-6);

        CHECK(sys.A(0, 1) == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(sys.A(1, 0) == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(sys.B(1, 0) == doctest::Approx(2.0).epsilon(1e-6));
        CHECK(sys.C(0, 0) == doctest::Approx(1.0).epsilon(1e-6));
        CHECK(sys.D(0, 0) == doctest::Approx(1.0).epsilon(1e-6));
    }

    TEST_CASE("linearize is accurate at a large operating point") {
        const auto f = [](const ColVec<1>& x, const ColVec<1>& u) {
            return ColVec<1>{(x(0) * x(0)) + (u(0) * u(0))};
        };
        const auto lin = design::linearize<1, 1>(f, ColVec<1>{1e10}, ColVec<1>{-3e10});
        CHECK(lin.A(0, 0) == doctest::Approx(2e10).epsilon(1e-6));
        CHECK(lin.B(0, 0) == doctest::Approx(-6e10).epsilon(1e-6));
    }

    TEST_CASE("linearize reports its status; a non-finite model or step is refused") {
        const auto dynamics = [](const ColVec<1>& x, const ColVec<1>& u) { return ColVec<1>{-x(0) + u(0)}; };
        const auto output = [](const ColVec<1>& x, const ColVec<1>&) { return ColVec<1>{x(0)}; };
        const auto ok = design::linearize<1, 1, 1>(dynamics, output, ColVec<1>{1.0}, ColVec<1>{0.0});
        CHECK(ok.success);
        CHECK(ok.status == design::DesignStatus::ok);
        const auto nan_model = [](const ColVec<1>& x, const ColVec<1>&) { return ColVec<1>{x(0) * std::numeric_limits<double>::quiet_NaN()}; };
        CHECK(design::linearize<1, 1, 1>(nan_model, output, ColVec<1>{1.0}, ColVec<1>{0.0}).status == design::DesignStatus::bad_weights);
        CHECK(design::linearize<1, 1, 1>(dynamics, output, ColVec<1>{1.0}, ColVec<1>{0.0}, std::numeric_limits<double>::quiet_NaN()).status == design::DesignStatus::invalid_parameter);
    }
}
