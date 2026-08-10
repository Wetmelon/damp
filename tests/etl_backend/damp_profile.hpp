// Copyright Paul Guénette and Damp contributors.
// SPDX-License-Identifier: BSL-1.0
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE or copy at http://www.boost.org/LICENSE_1_0.txt)

#pragma once

// ETL profile for this standalone suite: select the Embedded Template Library
// container backend instead of the stdlib default. damp_profile.hpp is MACRO-ONLY
// configuration — see damp/config.hpp. Found ahead of the host tests/damp_profile.hpp
// because this directory is on the include path first (-I.).

#define DAMP_BACKEND_ETL
