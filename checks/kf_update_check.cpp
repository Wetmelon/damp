// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "damp/estimation/kalman.hpp"
#include "damp/matrix/colvec.hpp"

using KF = damp::KalmanFilter<2, 1, 1, float, 2, 1>;

bool kf_update(KF& kf, const damp::ColVec<1, float>& y) {
    return kf.update(y);
}
