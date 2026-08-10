-- Parent tup projects: leave CONFIG_DAMP_SELF_BUILD unset so these no-op.
-- Standalone Damp: tup.config sets CONFIG_DAMP_SELF_BUILD=y (see docs/submodule_tup.md).
-- Official Lua API: getconfig("FOO") reads CONFIG_FOO from tup.config.
if tup.getconfig("DAMP_SELF_BUILD") ~= "y" then
    return
end

-- Standalone executable for the DAMP_MATH_BACKEND_DAMP profile. Kept out of the
-- main tests/ glob (and its single test_runner.exe) because the Damp
-- MathBackend<float> specialization is ODR-incompatible with the default std
-- backend every other test object links against. tup.foreach_rule globs are
-- per-directory, so the parent tests/Tupfile.lua never picks this file up.
INCLUDES = '-I. -I.. -I../../inc'

TEST_CXXFLAGS = CXXFLAGS..' -ffast-math'

obj = tup.foreach_rule('*.cpp', '^j^'..CXX..' '..TEST_CXXFLAGS..' '..INCLUDES..' '..WARNINGS..' -c %f -o %o', 'build/objs/%B.o')
tup.rule(obj, CXX..' '..TEST_CXXFLAGS..' '..LDFLAGS..' -static %f -o %o', 'build/test_damp_backend.exe')
