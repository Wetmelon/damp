// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

// Host profile for tests: stdlib containers + std:: math backend (the library
// defaults). damp_profile.hpp is MACRO-ONLY configuration — see damp/config.hpp
// for the recognized macros. The host defaults need none, so this file is empty
// apart from existing: its presence selects an explicit, warning-free config.

// #define DAMP_BACKEND_ETL
// #define DAMP_MATH_BACKEND_DAMP
