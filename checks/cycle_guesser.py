import re
import sys
from collections import defaultdict

# Approximate Cortex-M7 (Thumb-2 + FPv5-D16 hardware FPU) cycle costs, keyed by the
# *base* mnemonic (suffixes like .f32/.w/.n and condition codes are stripped first).
# These are rough throughput numbers: the M7 is a dual-issue superscalar pipeline, so a
# real count needs the hardware. This is a guess, not a model.
CYCLE_LOOKUP = {
    # Integer ALU / move / shift / bitfield — 1 cycle
    "add": 1, "adc": 1, "sub": 1, "sbc": 1, "rsb": 1, "adr": 1,
    "mov": 1, "movw": 1, "movt": 1, "mvn": 1, "neg": 1,
    "cmp": 1, "cmn": 1, "tst": 1, "teq": 1,
    "and": 1, "orr": 1, "orn": 1, "eor": 1, "bic": 1,
    "lsl": 1, "lsr": 1, "asr": 1, "ror": 1, "rrx": 1,
    "clz": 1, "rev": 1, "rev16": 1, "revsh": 1,
    "uxtb": 1, "uxth": 1, "sxtb": 1, "sxth": 1,
    "bfi": 1, "bfc": 1, "ubfx": 1, "sbfx": 1,
    "nop": 1, "it": 1, "ite": 1, "itt": 1, "itte": 1, "itee": 1, "ittt": 1,

    # Integer multiply / divide
    "mul": 1, "mla": 2, "mls": 2,
    "umull": 1, "smull": 1, "umlal": 1, "smlal": 1,
    "sdiv": 12, "udiv": 12,  # 2-12 cycles, data-dependent; assume worst-ish

    # Branches (taken branch refills the pipeline; M7 has prediction but assume a refill)
    "b": 2, "bl": 3, "bx": 2, "blx": 3, "cbz": 2, "cbnz": 2,

    # Single load / store
    "ldr": 2, "ldrb": 2, "ldrh": 2, "ldrsb": 2, "ldrsh": 2, "ldrd": 3,
    "str": 2, "strb": 2, "strh": 2, "strd": 3,
    # Load/store multiple: handled specially (1 + register count); table value is the base.
    "ldm": 1, "stm": 1, "push": 1, "pop": 1,

    # --- VFP (FPv5 single precision) ---
    "vldr": 2, "vstr": 2, "vldm": 1, "vstm": 1, "vpush": 1, "vpop": 1,
    "vmov": 1, "vmrs": 1, "vmsr": 1,
    "vadd": 1, "vsub": 1, "vmul": 1, "vnmul": 1,
    "vmla": 3, "vmls": 3, "vnmla": 3, "vnmls": 3,
    "vfma": 3, "vfms": 3, "vfnma": 3, "vfnms": 3,  # fused multiply-accumulate, ~3 cyc latency
    "vneg": 1, "vabs": 1, "vcmp": 1, "vcmpe": 1, "vminnm": 1, "vmaxnm": 1,
    "vcvt": 1, "vcvtr": 1, "vrintr": 1, "vrintz": 1, "vrintx": 1, "vsel": 1,
    "vdiv": 14, "vsqrt": 14,  # FPv5 single: ~14 cycles
}

# Addressing-mode spellings that objdump emits for the load/store-multiple family,
# folded onto the base mnemonic (which carries the cost + register-count handling).
LSM_ALIASES = {
    "ldmia": "ldm", "ldmdb": "ldm", "ldmib": "ldm", "ldmda": "ldm",
    "stmia": "stm", "stmdb": "stm", "stmib": "stm", "stmda": "stm",
    "vldmia": "vldm", "vstmia": "vstm", "vldmdb": "vldm", "vstmdb": "vstm",
}

# ARM condition-code suffixes, longest first so we don't strip a 1-letter prefix of a 2-letter code.
CONDITIONS = ["eq", "ne", "cs", "cc", "mi", "pl", "vs", "vc",
              "hi", "ls", "ge", "lt", "gt", "le", "al", "hs", "lo"]

# An instruction line is "  <hex addr>:\t<opcode bytes>\t<mnemonic>...".
INSN_RE = re.compile(r"^\s*[0-9a-fA-F]+:\t")
# A function symbol label: "00000000 <damp::sincos(float)>:".
LABEL_RE = re.compile(r"^[0-9a-fA-F]+ <(.+)>:\s*$")
# A new function section starts here; the next label line names it.
SECTION_RE = re.compile(r"^Disassembly of section ")
# Branch-and-link relocation carrying the *real* callee (the bl operand in a .o is bogus).
# "\t\t\t 3a: R_ARM_THM_CALL\tdamp::detail::sin_poly(float)"
RELOC_CALL_RE = re.compile(r"R_ARM\S*CALL\s+(.+?)\s*$")


def base_mnemonic(mnemonic):
    """Reduce 'vfma.f32' / 'ldrb.w' / 'addeq' / 'movs' / 'stmia.w' to a lookup key."""
    m = mnemonic.lower().split(".", 1)[0]  # drop .f32 / .w / .n / .i32 ...
    m = LSM_ALIASES.get(m, m)
    if m in CYCLE_LOOKUP:
        return m
    # Strip a trailing condition code, then a flag-setting 's', only if what remains is a
    # known mnemonic (so 'mls' isn't mangled into 'm' via the 'ls' condition, etc.).
    for suffix in CONDITIONS + ["s"]:
        if m.endswith(suffix) and m[: -len(suffix)] in CYCLE_LOOKUP:
            return m[: -len(suffix)]
    return None


def count_reglist(operands):
    """Number of registers in a {...} list, expanding ranges like r4-r11 / d8-d11."""
    m = re.search(r"\{([^}]*)\}", operands)
    if not m:
        return 1
    n = 0
    for part in (p.strip() for p in m.group(1).split(",")):
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            ai = re.sub(r"\D", "", a)
            bi = re.sub(r"\D", "", b)
            n += (int(bi) - int(ai) + 1) if (ai and bi) else 1
        else:
            n += 1
    return max(n, 1)


def insn_cost(mnemonic, operands):
    """(cycles, base) for one instruction; base is None for unknown mnemonics."""
    base = base_mnemonic(mnemonic)
    if base is None:
        return 1, None
    cycles = CYCLE_LOOKUP[base]
    if base in ("ldm", "stm", "push", "pop", "vldm", "vstm", "vpush", "vpop"):
        cycles = 1 + count_reglist(operands)
    return cycles, base


def simplify(name):
    """Drop template args, parameter lists, and the damp:: prefix for a readable label."""
    out, depth = [], 0
    for ch in name:
        if ch in "<(":
            depth += 1
        elif ch in ">)":
            depth = max(0, depth - 1)
        elif depth == 0:
            out.append(ch)
    return " ".join("".join(out).replace("damp::", "").split())


def parse_functions(text):
    """Split objdump -drC into per-function bodies with self cost and call edges.

    Returns {name: {'cyc','instr','vdiv','callees':[name,...]}}. Falls back to a single
    synthetic '<input>' function when there are no section/label markers (piped raw blob).
    """
    funcs = {}

    def get(name):
        return funcs.setdefault(name, {"cyc": 0, "instr": 0, "vdiv": 0, "callees": []})

    cur = None
    expect_label = False
    for line in text.splitlines():
        if SECTION_RE.match(line):
            expect_label = True
            continue
        m = LABEL_RE.match(line)
        if m and expect_label:
            cur = m.group(1)
            get(cur)
            expect_label = False
            continue
        rm = RELOC_CALL_RE.search(line)
        if rm:
            get(cur if cur is not None else "<input>")["callees"].append(rm.group(1))
            continue
        if INSN_RE.match(line):
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            mnemonic = fields[2].strip()
            if not mnemonic or mnemonic.startswith("."):
                continue  # assembler directive / constant-pool data
            operands = fields[3] if len(fields) > 3 else ""
            cycles, base = insn_cost(mnemonic, operands)
            f = get(cur if cur is not None else "<input>")
            f["cyc"] += cycles
            f["instr"] += 1
            if base in ("vdiv", "vsqrt"):
                f["vdiv"] += 1
    return funcs


def inclusive(funcs):
    """Inclusive (self + callees) cost per function. External callees (no body in this
    TU) are skipped — their cost runs in code we can't see. Recursion is broken at the
    back edge (self counted once) and flagged."""
    memo = {}

    def rec(name, stack):
        f = funcs.get(name)
        if f is None:
            return None  # external: body not in this translation unit
        if name in stack:
            return f["cyc"], f["vdiv"], True  # recursion back edge: self only
        if name in memo:
            return memo[name]
        cyc, vdiv, recursive = f["cyc"], f["vdiv"], False
        for callee in f["callees"]:
            r = rec(callee, stack | {name})
            if r is None:
                continue
            cyc += r[0]
            vdiv += r[1]
            recursive = recursive or r[2]
        result = (cyc, vdiv, recursive)
        if not recursive:
            memo[name] = result
        return result

    return {n: rec(n, frozenset()) for n in funcs}


def report(text):
    funcs = parse_functions(text)
    if not funcs:
        print("No instructions found.")
        return

    incl = inclusive(funcs)
    called = {c for f in funcs.values() for c in f["callees"]}
    externals = sorted(c for c in called if c not in funcs)

    rows = sorted(funcs, key=lambda n: -incl[n][0])
    name_w = min(54, max(len(simplify(n)) for n in rows))

    print(f"Per-function cycle estimate ({len(funcs)} functions).")
    print("  *  = entry point (never called within this TU)   r = recurses\n")
    print(f"  {'function':<{name_w}}  {'self':>6} {'incl':>7}  {'div(s/i)':>9}  calls")
    print(f"  {'-' * name_w}  {'-'*6} {'-'*7}  {'-'*9}  -----")
    for n in rows:
        self_c, incl_c = funcs[n]["cyc"], incl[n][0]
        self_v, incl_v = funcs[n]["vdiv"], incl[n][1]
        root = "*" if n not in called else " "
        rflag = "r" if incl[n][2] else " "
        ncalls = len(funcs[n]["callees"])
        print(f"{root} {simplify(n):<{name_w}}  {self_c:>6} {incl_c:>7}  "
              f"{self_v:>3}/{incl_v:<3}{rflag:>2}  {ncalls}")

    if externals:
        print("\nExternal calls (cost NOT included — body not in this TU):")
        for name in externals:
            print(f"  {simplify(name)}")

    print("-" * 50)
    total = sum(f["instr"] for f in funcs.values())
    print(f"Total instructions across all functions: {total}")
    print("Entry points (*) carry the per-call cost; multiply by each loop's rate.")


def _selftest():
    """Exercise the call-graph composition: a multi-call edge (A->B twice), an external
    callee (skipped), and a recursive back edge (broken, flagged). Run: --selftest."""
    sample = "\n".join([
        "Disassembly of section .text.A:",
        "",
        "00000000 <A>:",
        "   0:\tbf00      \tnop",
        "   4:\tf7ff fffe \tbl\t0 <B>",
        "\t\t\t4: R_ARM_THM_CALL\tB",
        "   8:\tf7ff fffe \tbl\t0 <B>",
        "\t\t\t8: R_ARM_THM_CALL\tB",
        "   c:\tf7ff fffe \tbl\t0 <ext>",
        "\t\t\tc: R_ARM_THM_CALL\text",
        "",
        "Disassembly of section .text.B:",
        "",
        "00000000 <B>:",
        "   0:\tee80 0a00 \tvdiv.f32\ts0, s0, s0",
        "   4:\tf7ff fffe \tbl\t0 <B>",
        "\t\t\t4: R_ARM_THM_CALL\tB",
    ])
    funcs = parse_functions(sample)
    assert funcs["A"]["cyc"] == 10, funcs["A"]                       # nop(1) + 3*bl(3)
    assert funcs["A"]["callees"] == ["B", "B", "ext"], funcs["A"]
    assert funcs["B"]["cyc"] == 17 and funcs["B"]["vdiv"] == 1, funcs["B"]  # vdiv(14)+bl(3)

    incl = inclusive(funcs)
    # B expands to self(17) + recursion back-edge self(17) = 34; A = self(10) + 2*34; ext skipped.
    assert incl["A"] == (78, 4, True), incl["A"]
    assert incl["B"] == (34, 2, True), incl["B"]

    called = {c for f in funcs.values() for c in f["callees"]}
    assert "ext" in called and "ext" not in funcs, "external not detected"
    assert "A" not in called, "A should be an entry point"
    print("selftest OK")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--selftest":
        _selftest()
        sys.exit(0)

    if len(sys.argv) > 1:
        with open(sys.argv[1], "r", encoding="utf-8", errors="replace") as f:
            text = f.read()
    elif not sys.stdin.isatty():
        text = sys.stdin.read()
    else:
        print("Usage: python cycle_guesser.py <objdump -drC output>   (or pipe it in)")
        print("       objdump must be run with -r so bl targets resolve.")
        sys.exit(1)

    print(f"Analyzing {sys.argv[1] if len(sys.argv) > 1 else 'stdin'}...\n")
    report(text)
