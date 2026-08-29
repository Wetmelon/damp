-- Parent tup projects: leave CONFIG_DAMP_SELF_BUILD unset so these no-op.
-- Standalone Damp: tup.config sets CONFIG_DAMP_SELF_BUILD=y (see docs/submodule_tup.md).
-- Official Lua API: getconfig("FOO") reads CONFIG_FOO from tup.config.
if tup.getconfig("DAMP_SELF_BUILD") ~= "y" then
    return
end

-- The example programs are hosted demos (fmt / plotlypp). They are not part of
-- the freestanding contract, so skip them under the ETL backend variant.
if BACKEND == 'ETL' then
    return
end

INCLUDES = '-I.'
INCLUDES += '-I../inc'
INCLUDES += '-I..'   -- repo root (examples may include shared helpers)
INCLUDES += '-isystem../libs'
INCLUDES += '-isystem../libs/fmt/include'
INCLUDES += '-isystem../libs/plotlypp/include'
INCLUDES += '-isystem../libs/json/single_include'

-- Compile all fmt source files into the <fmt> bin so the per-example link can
-- reference them with %<fmt> (bin paths carry the variant prefix correctly;
-- a bare $(var) group does not).
fmt_sources = {'../libs/fmt/src/format.cc', '../libs/fmt/src/os.cc'}
-- Third-party fmt: compile without the Damp warning set
tup.foreach_rule(fmt_sources, '^j^'..CXX..' -c %f '..CXXFLAGS..' '..INCLUDES..' -o %o', {'build/fmt/%B.o', '<fmt>'})

-- Compile every example .cpp: root demos + nested product folders.
-- Basenames must stay unique across the tree (cart_pole_sketch.cpp, …).
-- tup.glob hard-errors on a missing directory (cannot pcall-skip), so this list
-- is only product folders present in *this* tree. Append when adding a demo.
sources = tup.glob('*.cpp')
local product_dirs = {
    'control/cart_pole',
    'control/pendulum',
    'control/reaction_wheel',
    'control/pid',
    'control/lpf',
    'control/workflow_end_to_end',
    'control/adrc',
    'estimation/eskf',
    'estimation/ins_navigator',
    'estimation/ins_eskf',
    'estimation/ins_mechanization',
    'estimation/imu_pose',
    'estimation/encoder_velocity',
    'estimation/fo_plant_inertia',
    'motor/foc',
}
for _, dir in ipairs(product_dirs) do
    for _, f in ipairs(tup.glob(dir .. '/*.cpp')) do
        sources[#sources + 1] = f
    end
end
ex_objs = tup.foreach_rule(sources, '^j^'..CXX..' -c %f '..CXXFLAGS..' '..INCLUDES..' '..WARNINGS..' -o %o', 'build/obj/%B.o')
ex_objs.extra_inputs = {'<fmt>'}

-- Link with g++ (objects before LDFLAGS so -lws2_32 etc. resolve on MinGW)
tup.foreach_rule(ex_objs, CXX..' '..CXXFLAGS..' %f %<fmt> '..LDFLAGS..' -o %o', 'build/%B.exe')

-- ---------------------------------------------------------------------------
-- Run each example. Plot HTML under plots/<domain>/ is a *real* tup output of
-- the run rule — up-to-date → no re-run; delete a plot (or the stamp) → that
-- example re-runs. When an example gains/renames/drops a plot, update
-- plot_outputs below.
--
-- plots/js/plotly.min.js is a shared plotlypp side product (written once, reused):
-- exclude with ^ so it is not an undeclared output of every plotter.
--
-- CWD = examples/. Runner: tools/run_example.sh (portable sh; no Python).
-- ---------------------------------------------------------------------------

-- basename (no .exe / path) → files the default main() always writes under plots/
-- Public v0.1 product set only (match product_dirs above).
local plot_outputs = {
    encoder_velocity_sil = {
        'plots/estimation/encoder_velocity_absolute.html',
        'plots/estimation/encoder_velocity_hall.html',
    },
    imu_pose_sil = {'plots/estimation/imu_pose_3d.html'},
    ins_eskf_sil = {'plots/estimation/ins_eskf_3d.html'},
    ins_mechanization_sil = {'plots/estimation/ins_mechanization_3d.html'},
    cart_pole_sil = {
        'plots/control/cart_pole_lqr.html',
    },
    adrc_sil = {
        'plots/control/adrc_vs_pid.html',
    },
    pendulum_sil = {
        'plots/control/pendulum_sim.html',
        'plots/control/pendulum_sim_high_q.html',
        'plots/control/pendulum_phase.html',
    },
    foc_sil = {
        'plots/motor/foc_current_loop.html',
        'plots/motor/foc_id_iq.html',
    },
}

local function cpp_basename(path)
    local name = path:match('([^/\\]+)$') or path
    return (name:gsub('%.cpp$', ''))
end

for _, src in ipairs(sources) do
    local base = cpp_basename(src)
    local exe = 'build/' .. base .. '.exe'
    local stamp = 'build/run/' .. base .. '.ok'
    local outs = {
        stamp,
        '^plots/js/.*',       -- shared plotly.min.js at plots/js/
        '^plots/.*/js/.*',    -- legacy per-domain copies (write_html now centralizes)
    }
    local plots = plot_outputs[base]
    if plots then
        for _, p in ipairs(plots) do
            outs[#outs + 1] = p
        end
    end
    -- %f is not available outside foreach; name inputs/outputs explicitly.
    -- run script is an input so edits re-trigger.
    tup.rule(
        {exe, '../tools/run_example.sh'},
        'sh ../tools/run_example.sh ' .. exe .. ' ' .. stamp,
        outs
    )
end
