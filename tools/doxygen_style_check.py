#!/usr/bin/env python3
"""House-style greps for Doxygen comments under inc/damp.

Fails (exit 1) if any of the following appear in inc/damp/**/*.hpp:

  - //! < member-doc dialect (use ///<)
  - Markdown **bold** used as *editorial emphasis* (not math symbols)
  - @note Equivalent to MATLAB… / @note Mirrors MATLAB… (use Compare with)

Also requires every header to contain exactly one @file tag.

**bold** policy
  Allowed — math symbols / vectors / short notation tokens, e.g. **v**, **x**,
  **R0**, **i_q**, a single Greek letter. Prefer @f$\\mathbf{v}@f$ for real
  display math; markdown bold is fine for a short symbol in prose.
  Forbidden — emphasis on English prose: **not**, **Host-only.**, **out of
  scope**, multi-word phrases, section labels. State the fact without markup.

Usage (from repo root):
  python3 tools/doxygen_style_check.py
  # Windows without python3 on PATH:  py -3 tools/doxygen_style_check.py
  # Or:  make doxygen-style-check
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WET = ROOT / "inc" / "damp"

BANG_MEMBER = re.compile(r"//!<")
NOTE_EQUIV = re.compile(r"@note\s+Equivalent to MATLAB")
NOTE_MIRROR = re.compile(r"@note\s+Mirrors MATLAB")
# **bold** that is not the /** comment opener or **/ closer
BOLD = re.compile(r"(?<!/)\*\*(?!\*)([^*\n]+)\*\*")

# Math-symbol shapes (no spaces, no sentence punctuation):
#   v, x, A, R0, x1, i_q, v_d, ω (single Greek)
_MATH_SYMBOL = re.compile(
    r"^(?:"
    r"[A-Za-z]"  # single Latin letter (vector / matrix / scalar)
    r"|[A-Za-z][0-9]+"  # R0, x1
    # Multi-letter only with underscore parts so English words (not, only) fail:
    # i_q, x_d, omega_m
    r"|[A-Za-z][A-Za-z0-9]*(?:_[A-Za-z0-9]+)+"
    r"|[\u0370-\u03FF]"  # single Greek letter
    r")$"
)


def is_math_bold(inner: str) -> bool:
    """True if **inner** is a math symbol, not editorial emphasis."""
    s = inner.strip()
    if not s:
        return False
    if any(c.isspace() for c in s):
        return False
    if any(c in s for c in ".,:;!?"):
        return False
    return _MATH_SYMBOL.fullmatch(s) is not None


def check_file(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8")
    rel = path.relative_to(ROOT).as_posix()
    issues: list[str] = []

    n_file = text.count("@file")
    if n_file == 0:
        issues.append(f"{rel}: missing @file")
    elif n_file > 1:
        issues.append(f"{rel}: multiple @file tags ({n_file})")

    for i, line in enumerate(text.splitlines(), 1):
        if BANG_MEMBER.search(line):
            issues.append(f"{rel}:{i}: use ///< not //! < for member docs")
        if NOTE_EQUIV.search(line):
            issues.append(f"{rel}:{i}: use '@note Compare with MATLAB' not 'Equivalent to'")
        if NOTE_MIRROR.search(line):
            issues.append(f"{rel}:{i}: use '@note Compare with MATLAB' not 'Mirrors'")
        # strip /** and **/ noise for bold search
        stripped = line.replace("/**", "").replace("**/", "")
        for m in BOLD.finditer(stripped):
            inner = m.group(1)
            if not is_math_bold(inner):
                issues.append(
                    f"{rel}:{i}: editorial **bold** (use for math symbols only, "
                    f"not emphasis; saw **{inner}**)"
                )

    return issues


def main() -> int:
    if not WET.is_dir():
        print(f"doxygen_style_check: {WET} not found", file=sys.stderr)
        return 2

    issues: list[str] = []
    for path in sorted(WET.rglob("*.hpp")):
        issues.extend(check_file(path))

    if issues:
        print("doxygen-style-check FAILED:", file=sys.stderr)
        for msg in issues:
            print(f"  {msg}", file=sys.stderr)
        print(f"({len(issues)} issue(s))", file=sys.stderr)
        return 1

    print("doxygen-style-check OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
