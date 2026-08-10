-- Compiler selection (tup.config / tup.config.default):
--   CONFIG_DAMP_SELF_BUILD=y     — enable this repo's Tupfiles (unset when a submodule)
--   CONFIG_COMPILER_PATH=...    — directory containing the compiler (optional)
--   CONFIG_COMPILER_PREFIX=...  — e.g. arm-none-eabi- (optional; empty = host)
--   CONFIG_CXX=clang++           — full driver override (optional; wins over path/prefix)
--   CONFIG_BACKEND=ETL          — ETL container + freestanding math (optional)
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
local compiler_prefix = tup.getconfig("COMPILER_PREFIX")

if cxx_override ~= "" then
    -- Full override for clang tooling or a non-g++ driver:
    --   CONFIG_CXX=clang++
    --   CONFIG_CXX=C:/Tools/llvm/bin/clang++
    CXX = cxx_override
else
    -- Default: g++ (optionally under COMPILER_PATH, optionally prefixed for cross).
    CXX = compiler_path .. compiler_prefix .. "g++"
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

CXXFLAGS = "-O3 -std=c++20 -march=native -ffunction-sections -fdata-sections "
-- clang: raise fold-expression nesting (default 2048). Matrix copy historically
-- used an index_sequence fold over Rows*Cols; large fixed-N LU copies exceed that.
if string.find(CXX, "clang", 1, true) then
    CXXFLAGS = CXXFLAGS .. "-fbracket-depth=16384 "
end
LDFLAGS = "-O3 -march=native -ffunction-sections -fdata-sections -Wl,--gc-sections"
-- example_drive_governor --live needs Winsock (and joy optional); unused on other exes
local plat = tup.getconfig("TUP_PLATFORM")
if plat == "win32" or plat == "win64" then
    LDFLAGS = LDFLAGS .. " -lws2_32"
end

-- Backend selection (mirrors damp/config.hpp). Default = C++ stdlib.
-- Variant with CONFIG_BACKEND=ETL: ETL containers + freestanding series math.
BACKEND = tup.getconfig("BACKEND")
if BACKEND == "ETL" then
    CXXFLAGS = CXXFLAGS .. "-DDAMP_BACKEND_ETL -DDAMP_MATH_BACKEND_FREESTANDING -I../libs/etl/include "
end
