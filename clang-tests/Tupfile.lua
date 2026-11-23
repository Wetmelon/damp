-- Parent tup projects: leave CONFIG_DAMP_SELF_BUILD unset so these no-op.
-- Standalone Damp: tup.config sets CONFIG_DAMP_SELF_BUILD=y (see docs/submodule_tup.md).
-- Official Lua API: getconfig("FOO") reads CONFIG_FOO from tup.config.
if tup.getconfig("DAMP_SELF_BUILD") ~= "y" then
    return
end

-- The three test runners again, built with clang whatever CXX is, with the flags a
-- CONFIG_CXX=clang++ build gives them (Tuprules.lua): clang deduces, folds and evaluates at
-- compile time differently enough to fail where g++ passes. Not in the default build: make builds
-- tests, examples and checks, and `make clang-tests` builds this directory and runs the runners.
-- Without clang++ 19+ on PATH every output is empty and make reports the suites skipped
-- (tools/clang_cxx.sh). Not marked for compile_commands.json, which holds the g++ commands.
-- A clang build has already built and run them.
if string.find(CXX, 'clang', 1, true) then
    return
end

CLANG = 'sh ../tools/clang_cxx.sh '..CXX..' %o'

-- The main suite (tests/Tupfile.lua).
local main_includes = '-I../tests -I../inc -I.. -isystem../libs -isystem../libs/fmt/include ' ..
    '-isystem../libs/plotlypp/include -isystem../libs/json/single_include'
local backend_flags = CLANG_CXXFLAGS..' -ffast-math'
local main_flags = backend_flags..' '..CLANG_TEST_CXXFLAGS
objs = tup.foreach_rule('../tests/*.cpp', CLANG..' -c %f '..main_flags..' '..main_includes..' '..CLANG_WARNINGS..' -o %o', 'build/objs/%B.o')
objs += tup.foreach_rule({'../libs/fmt/src/format.cc', '../libs/fmt/src/os.cc'}, CLANG..' -c %f '..main_flags..' '..main_includes..' -o %o', 'build/objs/fmt/%B.o')
tup.rule(objs, CLANG..' '..backend_flags..' '..LDFLAGS..' -static %f -o %o', 'build/test_runner.exe')

-- The two backend runners (tests/damp_backend, tests/etl_backend), each with its own profile first.
local backends = {
    {dir = 'damp_backend', exe = 'test_damp_backend', includes = '-I../tests/damp_backend -I../tests -I../inc'},
    {dir = 'etl_backend', exe = 'test_etl_backend', includes = '-I../tests/etl_backend -I../tests -I../inc -isystem../libs/etl/include'},
}
for _, b in ipairs(backends) do
    local bobjs = tup.foreach_rule('../tests/'..b.dir..'/*.cpp', CLANG..' -c %f '..backend_flags..' '..b.includes..' '..CLANG_WARNINGS..' -o %o', 'build/'..b.dir..'/%B.o')
    tup.rule(bobjs, CLANG..' '..backend_flags..' '..LDFLAGS..' -static %f -o %o', 'build/'..b.exe..'.exe')
end
