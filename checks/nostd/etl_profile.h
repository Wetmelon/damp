// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

// What README.md tells a toolchain with no C++ standard headers to put in etl_profile.h.
#pragma once

#define ETL_NO_STL
#define ETL_NO_CPP_NAN_SUPPORT // quiet_NaN() from the NAN macro, never libm's nanf()
