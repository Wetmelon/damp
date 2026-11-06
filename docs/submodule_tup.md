# Using Damp as a Tup submodule

Tup discovers every `Tupfile` / `Tupfile.lua` under the parent's project root (the directory with `.tup/` / `Tupfile.ini`). If Damp is cloned as `libs/damp` (or similar) into a Tup-using app, the parent's `tup` run would also parse Damp's `tests/`, `examples/`, and `checks/` Tupfiles and try to build the whole library suite, which is wrong for a header-only consumer.

All Damp Tupfiles are therefore no-ops unless `tup.config` sets:

```text
CONFIG_DAMP_SELF_BUILD=y
```

Lua reads it as `tup.getconfig("DAMP_SELF_BUILD")`, without the `CONFIG_` prefix. The flag is set in this repo's `tup.config.default` (and your local `tup.config` after the copy). Parent projects must not set it; Damp's Tupfiles then parse as empty and contribute no rules.

Consumer checklist:

1. The submodule (or tree) only needs headers on the include path, e.g. `-I libs/damp/inc`.
2. Do not set `CONFIG_DAMP_SELF_BUILD` in the parent `tup.config`.
3. Optional: git sparse-checkout to only `inc/` for zero Tupfiles in the tree (not required once the guard is in place).

Developing Damp itself: `tup.config` contains `CONFIG_DAMP_SELF_BUILD=y` (from `tup.config.default`); `make` / `tup` then builds tests and examples.

Tup does not support nested independent databases or an official switch to exclude a subtree: one project root means one database, and every Tupfile under it is eligible. The config guard is the portable, in-tree solution for header-only submodules.
