# Known issues

Compiler bugs, and compiler behavior, that change the results of code built with Damp. Each entry gives the compilers and flags, what goes wrong, what Damp does about it, and what your own code needs. Checked 2026-11-08. Scope limits of the library itself are in [known_limitations.md](known_limitations.md).

| Issue | Compilers | Damp | Your code |
| --- | --- | --- | --- |
| `pow(0, y)` folded to 0 | GCC 15.1 to 16.2 and trunk, `-ffast-math` | Not affected | Affected |
| AVX stack-alignment fault | MinGW GCC 14 to 16, AVX enabled | Exposed; this repo's build avoids it | Affected with `-march=native` |
| Vectorizer wrong code | GCC 15.1 and 15.2, `-O3 -ffast-math` | Worked around | Possible |
| Fused multiply-add and fast-math folds | Any FMA target | Written for it | Exact comparisons break |

## `pow(0, y)` folded to 0

When the base is a constant zero after inlining, GCC 15.1 and later replace `pow(0.0, y)` and `powf(0.0f, y)` by the constant 0, without calling `pow`, whatever `y` is. `pow(x, ±0)` is 1 for every x (C Annex F, F.10.4.4), so `pow(0, 0)` comes out 0.

- Compilers: GCC 15.1, 15.2, 15.3, 16.1, 16.2 and trunk 17.0.0 20261108. GCC 14 is correct. Seen on x86-64, x86_64-w64-mingw32, and arm-none-eabi 15.2 and 16.1 for Cortex-M4 with hard float.
- Trigger: `-funsafe-math-optimizations`, `-ffinite-math-only` and `-fno-math-errno` all on. `-ffast-math` turns on all three.
- Where it bites: a power basis evaluated at an end point. Examples are tᵏ at t = 0, (1 − t)ᵏ at t = 1 (Bernstein, Bézier), a polynomial's derivative at a knot, and Taylor terms at the expansion point. Every k = 0 term reads 0.
- Damp: not affected. `damp::pow(x, y)` returns 1 for y = 0 and handles a zero base itself, before it calls the math backend. `damp::pow(x, n)` with an integer `n` is a product.
- Your code: call `damp::pow`, or form integer powers as products. Otherwise add any one of `-fno-finite-math-only`, `-fsigned-zeros`, `-fmath-errno` or `-ftrapping-math`.
- Upstream: not reported yet.

```cpp
#include <cmath>
#include <cstdio>

int main(int argc, char**) {
    const double y = argc - 1; // 0, unknown to the compiler
    std::printf("pow(0, %g) = %g\n", y, std::pow(0.0, y)); // 1; GCC 15+ with -ffast-math: 0
}
```

## AVX stack-alignment fault

Windows x64 aligns the stack to 16 bytes. With AVX enabled, MinGW GCC copies some 32- and 64-byte aggregates passed by value into stack slots with aligned vector moves (`vmovdqa`, `vmovdqa64`). Examples are a struct holding a `double[8]` or a `char[32]`. The program faults (`0xC0000005`) on about half the runs, depending on where the stack lands. This is [GCC PR 54412](https://gcc.gnu.org/bugzilla/show_bug.cgi?id=54412).

- Compilers: x86_64-w64-mingw32 GCC 14.2, 15.2 and 16.2 were checked. Trunk is fixed by [8aaa5f6e62](https://github.com/gcc-mirror/gcc/commit/8aaa5f6e62), so GCC 17.1 is the first release with the fix unless it is backported. Only the Windows x64 target is affected: elsewhere GCC realigns the stack.
- Trigger: a 64-byte aggregate with AVX-512, for example `-march=native` on Zen 4 or an AVX-512 Intel host. A 32-byte one with `-mavx` on GCC 15 and 16. GCC 14 copies the 32-byte case as two 16-byte halves.
- Damp: its headers pass `Matrix` and option structs by value, so they are exposed like any code. This repo's build tunes for the host (`-mtune=native`) and never targets it (`-march=native`). See `Tuprules.lua`.
- Your code: with MinGW GCC 16 or older, leave out `-march=native`. If you need AVX, add `-Wa,-muse-unaligned-vector-move`: the assembler rewrites every aligned vector move as unaligned. `-mprefer-vector-width=128` also clears the 64-byte case.

```cpp
struct A { double d[8]; }; // 64 bytes, alignment 8

__attribute__((noinline)) double f(A a) { return a.d[0]; }

int main() {
    volatile double s = f(A{1, 2, 3, 4, 5, 6, 7, 8}); // g++ -O3 -march=native on an AVX-512 host
    (void)s;
}
```

## Vectorizer wrong code

GCC 15.1 and 15.2 drop a whole product term from a hand-written complex multiply-accumulate `h += M · prev`, when each `h` feeds the next iteration. This happens at `-O3 -ffast-math` with the baseline x86-64 instruction set (SSE2). The error is O(1), not rounding.

- Compilers: GCC 15.1 and 15.2, on x86-64 and x86_64-w64-mingw32. Fixed in 15.3; 14.x and 16.x are correct. Other targets were not checked.
- Clean with GCC 15.2: `-O2`, `-fno-tree-loop-vectorize`, or `-march=x86-64-v3`.
- Damp: `design::place_jordan` hit it while building eigenvector chains. They are `Matrix` products now, which compile correctly.
- Your code: use GCC 15.3 or later, or add `-fno-tree-loop-vectorize` on 15.1 and 15.2.

## Fused multiply-add and fast-math folds

On FMA hardware, GCC rounds a·b + c once instead of twice. This hardware includes Cortex-M4F, M7 and M33 with an FPU, and x86-64-v3, the default target of Ubuntu 26.04's GCC 15. In C++, GCC fuses by default (`-ffp-contract=fast`), whatever `-std` is, not only under `-ffast-math`. GCC 15 also folds more under `-ffast-math`. For example, `1 − exp(−x) > 0` becomes `x > 0`, so a guard against the difference rounding to 0 goes away.

- Two inlined copies of one formula can round differently, so they need not compare equal. Which copies fuse can change with `-mtune`.
- A difference that cancels exactly in real arithmetic, such as x·y − x·y, can keep the rounding error of one product.
- Damp: no decision waits for a cancellation to round to exactly 0, and the tests compare to a few ulps of the operands.
- Your code: compare to a tolerance scaled to the operands. Use `-ffp-contract=off` where you need bit-for-bit results across builds.
