-- Compiler selection (tup.config / tup.config.default):
--   CONFIG_DAMP_SELF_BUILD=y     — enable this repo's Tupfiles (unset when a submodule)
--   CONFIG_COMPILER_PATH=...    — directory containing the compiler (optional)
--   CONFIG_CXX=clang++           — full driver override (optional; wins over the path)
--
-- These build and run host programs. Cross builds are PlatformIO's job: see targets/.
--
-- Official API: getconfig("FOO") reads CONFIG_FOO. Default driver is g++ on PATH.

local function with_slash(path)
    if path == "" then
        return ""
    end
    local last = path:sub(-1)
    if last ~= "/" and last ~= "\\" then
        return path .. "/"
    end
    return path
end

local cxx_override = tup.getconfig("CXX")
local compiler_path = with_slash(tup.getconfig("COMPILER_PATH"))

if cxx_override ~= "" then
    -- Full override for clang tooling or a non-g++ driver:
    --   CONFIG_CXX=clang++
    --   CONFIG_CXX=C:/Tools/llvm/bin/clang++
    CXX = cxx_override
else
    -- Default: g++, optionally under COMPILER_PATH.
    CXX = compiler_path .. "g++"
end

local plat = tup.getconfig("TUP_PLATFORM")
local on_windows = plat == "win32" or plat == "win64"

-- clang on Windows targets MSVC by default, and Damp does not support MSVC: use MinGW's libstdc++.
-- Emulated TLS too: Clang 19's native TLS crashes doctest's thread_local on that target.
if on_windows and string.find(CXX, "clang", 1, true) then
    CXX = CXX .. " --target=x86_64-w64-mingw32 -femulated-tls"
end

-- Warning set: Jason Turner / cppbestpractices "GCC / Clang" list
-- (https://github.com/cpp-best-practices/cppbestpractices/blob/master/02-Use_the_Tools_Available.md#compilers)
-- plus -Wfloat-conversion (subset of -Wconversion; kept explicit) and driver-specific silences.
WARNINGS_COMMON =
    "-Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wold-style-cast -Wcast-align " ..
    "-Wunused -Woverloaded-virtual -Wpedantic -Wconversion -Wsign-conversion " ..
    "-Wnull-dereference -Wdouble-promotion -Wfloat-conversion -Wformat=2 -Wimplicit-fallthrough"
-- GCC-only extras from the same list (unknown on clang → hard error without -Wno-unknown-warning-option).
WARNINGS_GCC =
    "-Wmisleading-indentation -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wuseless-cast -Wno-psabi"
-- clang: -ffast-math makes isfinite a no-op; silence the "nan/inf disabled" note on our bit-pattern isfinite.
WARNINGS_CLANG = "-Wno-nan-infinity-disabled"

if string.find(CXX, "clang", 1, true) then
    WARNINGS = WARNINGS_COMMON .. " " .. WARNINGS_CLANG
else
    WARNINGS = WARNINGS_COMMON .. " " .. WARNINGS_GCC
end

-- Tune for the host (-mtune=native), never target it (-march=native): on an AVX-512 host
-- (znver4), g++ 14 (xPack MinGW) under -O3 -march=native stores structs passed by value
-- with 64-byte-aligned vmovdqa64 onto a stack Windows x64 only aligns to 16 bytes. The
-- process dies (0xC0000005) whenever ASLR leaves the slot misaligned: the adaptive ODE
-- tests hit it. Clang is clean either way; -mtune keeps the host scheduling without that instruction selection.
CXXFLAGS_COMMON = "-O3 -std=c++20 -mtune=native -ffunction-sections -fdata-sections "
-- clang: raise fold-expression nesting (default 2048). Matrix copy historically
-- used an index_sequence fold over Rows*Cols; large fixed-N LU copies exceed that.
CXXFLAGS_CLANG = "-fbracket-depth=16384 "
-- GCC: compile-time designs of about 10 states or more exceed the default 2^25 constexpr ops.
CXXFLAGS_GCC = "-fconstexpr-ops-limit=67108864 "
if string.find(CXX, "clang", 1, true) then
    CXXFLAGS = CXXFLAGS_COMMON .. CXXFLAGS_CLANG
else
    CXXFLAGS = CXXFLAGS_COMMON .. CXXFLAGS_GCC
end

-- clang's flags whatever CXX is: clang-tests/Tupfile.lua builds the test suites with clang too
-- (make clang-tests), with the flags a CONFIG_CXX=clang++ build gives them.
CLANG_CXXFLAGS = CXXFLAGS_COMMON .. CXXFLAGS_CLANG
CLANG_WARNINGS = WARNINGS_COMMON .. " " .. WARNINGS_CLANG
-- The main suite's on clang (tests/Tupfile.lua, clang-tests/Tupfile.lua). Its consteval designs
-- exceed clang's default budget (.clangd matches). Product sketches must not: examples build at
-- the default, and tools/clang_budget_check.sh checks them. -fno-finite-math-only because these
-- tests spell NaN and Inf as numeric_limits constants, which clang folds to poison under
-- finite-math-only (a NaN the design then never sees). tests/damp_backend keeps the flag and
-- checks data that arrives at run time, which is what firmware sees.
CLANG_TEST_CXXFLAGS = "-fconstexpr-steps=100000000 -fno-finite-math-only"

LDFLAGS = "-O3 -mtune=native -ffunction-sections -fdata-sections -Wl,--gc-sections"
-- example_drive_governor --live needs Winsock (and joy optional); unused on other exes
if on_windows then
    LDFLAGS = LDFLAGS .. " -lws2_32"
end
