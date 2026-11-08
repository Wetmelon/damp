// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file lqi_sil.cpp
 * @brief Position LQI — host plant, same tick as the sketch; clamp-aware vs naive
 */

#include "damp/backend.hpp"
#include "damp/controllers/lqi.hpp"
#include "damp/matrix/colvec.hpp"
#include "fmt/base.h"
#include "fmt/core.h"
#include "lqi_controller.hpp"

using namespace damp;
using namespace damp::examples_lqi;

static constinit Law<lqi_result, float> aware;
static constinit Law<lqi_result, float> naive;

static constexpr auto kModel = *sysd;

int main() {
    fmt::print("===== LQI + clamp (host SIL, same tick as the sketch) =====\n\n");
    fmt::print("Ts={:.0f} Hz, |u|≤{:.2f}, r=1\n\n", 1.0 / Ts, u_max);

    ColVec<2> x_a{};
    ColVec<2> x_n{};
    double    peak_a = 0.0;
    double    peak_n = 0.0;
    for (int k = 0; k < 800; ++k) {
        const float ua = control_period(aware, x_a.as<float>(), 1.0f);
        naive.update(x_n.as<float>());

        const float u_req = naive.control(ColVec<1, float>{1.0f})[0];
        const float un = damp::clamp(u_req, static_cast<float>(u_min), static_cast<float>(u_max));
        naive.back_calculate(ColVec<1, float>{u_req}); // unsaturated request
        naive.predict();

        x_a = kModel.A * x_a + kModel.B * ColVec<1>{static_cast<double>(ua)};
        x_n = kModel.A * x_n + kModel.B * ColVec<1>{static_cast<double>(un)};
        peak_a = damp::max(peak_a, x_a[0]);
        peak_n = damp::max(peak_n, x_n[0]);
    }
    fmt::print("aware  overshoot={:.3f}  final={:.3f}\n", peak_a - 1.0, x_a[0]);
    fmt::print("naive  overshoot={:.3f}  final={:.3f}\n", peak_n - 1.0, x_n[0]);
    fmt::print("aware uses back_calculate(u_sat); naive reports the unsaturated request.\n");
    return 0;
}
