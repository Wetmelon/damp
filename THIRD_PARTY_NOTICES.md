# Third-party notices

Damp itself is licensed under the Boost Software License 1.0 — see
[LICENSE](LICENSE). That license applies to Damp sources under `inc/damp/`, and to
first-party tests, examples, and tooling written for this project.

This file lists third-party software that may appear in this repository or
in builds that opt into optional dependencies. None of these are part of the
BSL-licensed Damp library sources under `inc/damp/`. Their licenses remain in
force for those components only.

## Summary

| Component | Path | Role | License | In default MCU binary? |
| --------- | ---- | ---- | ------- | ---------------------- |
| Damp | `inc/damp/` | Product library | BSL-1.0 | Yes (if you use it) |
| {fmt} | `libs/fmt/` | Host formatting (tests/examples) | MIT-style (see `libs/fmt/LICENSE`) | No |
| nlohmann/json | `libs/json/` | Host JSON for plotting stack | MIT (`libs/json/LICENSE.MIT`) | No |
| plotlypp | `libs/plotlypp/` | Host Plotly HTML figures | MIT (`libs/plotlypp/LICENSE`) | No |
| ETL | `libs/etl/` | Optional freestanding backend | MIT (`libs/etl/LICENSE`) | Only if built with `DAMP_BACKEND_ETL` |
| tl::expected | `libs/expected.hpp` | Available header (not wired into core) | CC0-1.0 | No (not included by `inc/damp`) |
| doctest | `tests/doctest.h` | Unit-test framework | MIT (header comment) | No (tests only) |

The C++ standard library (or another toolchain runtime) used by a build is
covered by the compiler / standard-library license, not by Damp’s BSL.

## Optional ETL backend

When `DAMP_BACKEND_ETL` is defined, `inc/damp/backend.hpp` aliases containers and
utilities to the [Embedded Template Library](https://www.etlcpp.com). That build
does compile ETL headers into the program; you must comply with ETL’s MIT
license for that configuration. The default backend uses the hosted
standard library only.

## Host plotting stack

Including `damp/simulation/plot_plotly.hpp` (or building examples that write
Plotly HTML) requires plotlypp and nlohmann/json on the include path.
Those are host-only; they are not pulled by `damp/control.hpp`.

## Where to find full license texts

| Component | Full text |
| --------- | --------- |
| Damp | [LICENSE](LICENSE) |
| {fmt} | [libs/fmt/LICENSE](libs/fmt/LICENSE) |
| nlohmann/json | [libs/json/LICENSE.MIT](libs/json/LICENSE.MIT) |
| plotlypp | [libs/plotlypp/LICENSE](libs/plotlypp/LICENSE) |
| ETL | [libs/etl/LICENSE](libs/etl/LICENSE) |
| tl::expected | CC0 dedication in the header of [libs/expected.hpp](libs/expected.hpp) |
| doctest | Copyright / MIT notice at the top of [tests/doctest.h](tests/doctest.h) |

## Submodules

Third-party trees under `libs/` are git submodules (see `.gitmodules`). They are
not dual-licensed as Damp; each keeps its upstream license.

---

If you redistribute Damp source, keep Damp’s BSL notices and this notice file
(or equivalent attribution) for any third-party trees you also distribute.  
If you redistribute only machine-executable object code generated from Damp
(and the standard library), BSL does not require reproducing the BSL text in
that binary alone — see [LICENSE](LICENSE).
