-- Parent tup projects: leave CONFIG_DAMP_SELF_BUILD unset so these no-op.
-- Standalone Damp: tup.config sets CONFIG_DAMP_SELF_BUILD=y (see docs/submodule_tup.md).
-- Official Lua API: getconfig("FOO") reads CONFIG_FOO from tup.config.
if tup.getconfig("DAMP_SELF_BUILD") ~= "y" then
    return
end

INCLUDES = '-I.'
INCLUDES += '-I../inc'
INCLUDES += '-I..'   -- sil/ host plant ODEs (not embeddable)
-- Third-party: -isystem so their headers are not held to our -W* set.
INCLUDES += '-isystem../libs'
INCLUDES += '-isystem../libs/fmt/include'
INCLUDES += '-isystem../libs/plotlypp/include'
INCLUDES += '-isystem../libs/json/single_include'

-- The test runner builds with -ffast-math by default. Constexpr static_asserts
-- across the library act as a tripwire for the property that compile-time
-- evaluation stays IEEE-strict even when the optimizer is set to fast-math. If
-- a future toolchain leaks fast-math into constant evaluation, the suite will
-- fail on the next build. Drop the flag here to spot-check strict-IEEE runtime
-- behavior, or lower -O for precision investigations.
--
-- Drop -march=native for tests: g++ 14 (xPack MinGW) miscompiles multi-call
-- adaptive ODE paths under -O3 -march=native (intermittent 0xC0000005). Clang
-- and plain -O3 without -march=native are fine. Keep -mtune=native for scheduling.
TEST_CXXFLAGS = (CXXFLAGS:gsub('%-march=native', '-mtune=native'))..' -ffast-math'
TEST_LDFLAGS = (LDFLAGS:gsub('%-march=native', '-mtune=native'))

-- Compile all .cpp files in the tests directory
objs = tup.foreach_rule('*.cpp', '^j^'..CXX..' -c %f '..TEST_CXXFLAGS..' '..INCLUDES..' '..WARNINGS..' -o %o', 'build/objs/%B.o')

-- Compile all fmt source files (third-party: do not apply the Damp warning set)
fmt_sources = {'../libs/fmt/src/format.cc', '../libs/fmt/src/os.cc'}
objs += tup.foreach_rule(fmt_sources, '^j^'..CXX..' -c %f '..TEST_CXXFLAGS..' '..INCLUDES..' -Wno-error -o %o', 'build/objs/fmt/%B.o')

-- Link with g++
tup.rule(objs, CXX..' '..TEST_CXXFLAGS..' '..TEST_LDFLAGS..' -static %f -o %o', 'build/test_runner.exe')

-- Run test executable
-- tup.frule{inputs = 'build/test_runner.exe', command = './%f', outputs = {}}
