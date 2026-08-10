-- Host contract stamps + Cortex-M7 asm dumps (incremental).
-- Parent projects: leave CONFIG_DAMP_SELF_BUILD unset so these no-op.
-- Standalone Damp: tup.config sets CONFIG_DAMP_SELF_BUILD=y (see docs/submodule_tup.md).
if tup.getconfig("DAMP_SELF_BUILD") ~= "y" then
    return
end

-- ---------------------------------------------------------------------------
-- Contract stamps (embedded umbrella / freestanding ETL)
-- ---------------------------------------------------------------------------

-- One-level globs (tup does not expand **). List only directories that exist in
-- *this* tree — tup.glob hard-errors on a missing dir (pcall cannot catch it),
-- and public v0.1 is deliberately sparse (no motor/power packs; kinematics is
-- pose.hpp only). When adding a pack directory, append it here.
local damp_dirs = {
    '../inc/damp',
    '../inc/damp/analysis',
    '../inc/damp/controllers',
    '../inc/damp/design',
    '../inc/damp/estimation',
    '../inc/damp/filters',
    '../inc/damp/kinematics',
    '../inc/damp/math',
    '../inc/damp/matrix',
    '../inc/damp/simulation',
    '../inc/damp/systems',
    '../inc/damp/toolbox',
    '../inc/damp/trajectory',
}

local headers = { '../tools/contract_check.sh' }
for _, d in ipairs(damp_dirs) do
    for _, f in ipairs(tup.glob(d .. '/*.hpp')) do
        headers[#headers + 1] = f
    end
end

-- freestanding: ETL + freestanding math compile + no hosted-std leaks in damp/
-- Headers are *inputs* (rebuild when they change). On Windows, g++/grep often
-- open those headers with write-sharing; tup then reports FILE OVERWRITTEN even
-- though content is unchanged. ^inc/.* ignores those false-positive writes
-- (same idea as examples' ^plots/.*). Probe temps live under build/.
tup.rule(
    headers,
    'sh ../tools/contract_check.sh freestanding ' .. CXX .. ' %o',
    {
        'build/freestanding.ok',
        '^build/probe.*',
        '^inc/.*',
        '^\\.\\./inc/.*',
    }
)

-- embedded: default profile must not pull <vector> via control.hpp
tup.rule(
    headers,
    'sh ../tools/contract_check.sh embedded ' .. CXX .. ' %o',
    {
        'build/embedded.ok',
        '^build/probe.*',
        '^inc/.*',
        '^\\.\\./inc/.*',
    }
)

-- ---------------------------------------------------------------------------
-- Asm dumps (arm-none-eabi, Cortex-M7) — cycle inspection via cycle_guesser.py
-- ---------------------------------------------------------------------------

local asm_includes = '-I. -I../inc -isystem../libs -isystem../libs/etl/include'
    .. ' -isystem../libs/fmt/include -isystem../libs/plotlypp/include -isystem../libs/json/single_include'

local asm_cxxflags =
    '-std=c++20 -O3 -ffast-math -mcpu=cortex-m7 -mfloat-abi=hard -mfpu=fpv5-d16'
    .. ' -fno-exceptions -fno-rtti -fno-threadsafe-statics'

local asm_objs = tup.foreach_rule(
    '*.cpp',
    'arm-none-eabi-g++ -c %f ' .. asm_cxxflags .. ' ' .. WARNINGS .. ' ' .. asm_includes
        .. ' -o %o -ffunction-sections -fdata-sections',
    'build/obj/%B.o'
)

tup.foreach_rule(asm_objs, 'arm-none-eabi-objdump -drC %f > %o', 'build/%B.asm')
