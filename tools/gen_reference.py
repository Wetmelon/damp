#!/usr/bin/env python3
"""Scrape inc/damp headers for public types/functions + @brief into a categorized
markdown reference. Gate on column-0 doc comments so members/converters are skipped."""
import re, pathlib, collections

ROOT = pathlib.Path("inc/damp")

# Pretty names + ordering (core rings for public tree).
# Grouping: substrate → design-time → runtime → helpers → host.
CATS = [
    # --- Substrate ---
    ("",            "Core, configuration & backend vocabulary"),
    ("math",        "Scalar math, complex & frames"),
    ("matrix",      "Linear algebra"),
    ("systems",     "LTI systems (SS / TF / ZPK / discretize)"),
    # --- Design-time (not PWM-rate) ---
    ("design",      "Design-time synthesis (not PWM-rate)"),
    ("plants",      "Plant builders (optional grade B/C)"),
    # --- Runtime (generic) ---
    ("controllers", "Runtime controllers"),
    ("estimation",  "Observers & estimators"),
    ("filters",     "Filters & signal conditioning"),
    # --- Shared helpers (public core; full packs deferred) ---
    ("trajectory",  "Trajectory value types"),
    ("kinematics",  "Kinematics / pose"),
    ("motor",       "Motor control pack (if present)"),
    ("power",       "Power electronics pack (if present)"),
    # --- Embedded helpers ---
    ("toolbox",     "Embedded helpers (controls-adjacent utilities)"),
    # --- Host ---
    ("analysis",    "Frequency-domain analysis (host)"),
    ("simulation",  "Simulation / SIL harness (host)"),
    ("matlab",      "MATLAB®-style aliases (host)"),
]
CAT_ORDER = {k: i for i, (k, _) in enumerate(CATS)}
CAT_NAME = dict(CATS)

# /** ... */ block that starts at column 0, plus the code that follows it.
BLOCK = re.compile(r"^/\*\*(.*?)\*/\n(.*?)(?=\n/\*\*|\n///|\n\}|\Z)", re.S | re.M)
# Column-0 /// one-liner run, plus the code that follows it (make_*_servo, enums, ...).
LINE_DOC = re.compile(r"^((?:///[^\n]*\n)+)(.*?)(?=\n/\*\*|\n///|\n\}|\Z)", re.S | re.M)
# Brief terminates at the next @tag, a blank comment line, the comment end, or the
# end of the captured innards — BLOCK strips the trailing */, so a doc comment that
# is only an @brief line has no other terminator ($ catches it).
BRIEF = re.compile(r"@brief\s+(.*?)(?:\n\s*\*\s*@|\n\s*\*\s*\n|\*/|$)", re.S)

# Internal, compile-time-selected math backends. Each re-declares the same
# damp:: scalar surface (sin/cos/sqrt/...), so listing their functions is pure
# duplication — they get one file-level section instead. math.hpp (the public
# dispatcher) and complex.hpp stay in the scalar-math tables.
BACKEND_FILES = {"trig.hpp", "math_backend.hpp", "damp_backend.hpp",
                 "series_backend.hpp", "std_fallback.hpp", "constexpr_math.hpp"}


def doxygen_math_to_md(t):
    """Doxygen formula mode → markdown math $`...`$ (KaTeX-friendly tables).

    @f$...@f$ (inline) and @f[...@f] (display) become $`...`$ so REFERENCE.md
    renders in GFM viewers that accept that form. Collapse internal newlines so
    table cells stay single-line. Escape | inside the math body once for tables.
    """
    def squash(m):
        body = re.sub(r"\s+", " ", m.group(1).strip())
        # Table-safe: escape bare | only (keep LaTeX \| norms as single backslash).
        body = re.sub(r"(?<!\\)\|", r"\\|", body)
        return f"$`{body}`$"

    # Display before inline so nested edge cases don't double-match.
    t = re.sub(r"@f\[\s*(.*?)\s*@f\]", squash, t, flags=re.S)
    t = re.sub(r"@f\$\s*(.*?)\s*@f\$", squash, t, flags=re.S)
    return t


def escape_table_pipes(t):
    """Escape | outside $`...`$ segments (math bodies are already escaped)."""
    parts = re.split(r"(\$`.*?`\$)", t, flags=re.S)
    out = []
    for p in parts:
        if p.startswith("$`") and p.endswith("`$"):
            out.append(p)
        else:
            out.append(p.replace("|", "\\|"))
    return "".join(out)


def clean_brief(t):
    t = t.strip().rstrip(".")
    t = doxygen_math_to_md(t)
    t = re.sub(r"@(?:ref|c|p|a)\s+", "", t)  # drop Doxygen inline tags
    return escape_table_pipes(t)


def brief_text(body):
    m = BRIEF.search(body)
    if not m:
        return None
    return clean_brief(re.sub(r"\s*\*\s*", " ", m.group(1)))


def line_doc_text(comment):
    """First sentence-ish of a column-0 /// run."""
    lines = [re.sub(r"^///<?\s*", "", l.strip()) for l in comment.strip().split("\n")]
    return clean_brief(" ".join(l for l in lines if l))


def decl(code):
    """Return (kind, name, line_offset) for the declaration following a doc block,
    or None. line_offset is the 0-based line index (within code) of the declaration."""
    lines = [l for l in code.split("\n")]
    # skip blank/attribute/template prefix lines to reach the named declaration
    buf = []
    tmpl_depth = 0
    for i, l in enumerate(lines):
        s = l.strip()
        if not s:
            continue
        if s.startswith("#") or s.startswith("namespace") or s.startswith("using"):
            return None  # file-level / namespace brief, not a public entity
        # Consume template parameter headers (possibly line-wrapped) and requires-
        # clauses whole: their `class X` / `requires(...)` tokens would otherwise be
        # picked up as the entity name (DsogiPll's Resonator param, "requires" fns).
        if tmpl_depth > 0 or s.startswith("template"):
            tmpl_depth = max(0, tmpl_depth + s.count("<") - s.count(">"))
            continue
        if s.startswith("requires"):
            continue
        buf.append(s)
        joined = " ".join(buf)
        m = re.search(r"\bconcept\s+(\w+)", joined)
        if m:
            return ("concept", m.group(1), i)
        m = re.search(r"\b(?<!enum )(struct|class)\s+(\w+)", joined)
        if m and "{" in joined or (m and ";" not in joined and len(buf) <= 3):
            return ("type", m.group(2), i)
        m = re.search(r"\benum\s+(?:class\s+)?(\w+)", joined)
        if m:
            return ("enum", m.group(1), i)
        # function: identifier immediately before '(' once we've cleared prefixes
        if "(" in joined:
            fm = re.search(r"(\w+)\s*\(", joined)
            if fm and fm.group(1) not in ("if", "for", "while", "return", "as"):
                # require it looks like a free function decl (has a return-type token)
                if re.search(r"\b(struct|class|enum)\b", joined):
                    continue
                return ("fn", fm.group(1), i)
        if len(buf) > 4:
            break
    return None


def category(path):
    # Taxonomy: folder-first with D25 special cases (plants vs design core; AC math).
    # TODO(D20/D25): prefer declared domain namespace over folder when migration is done.
    if path.name == "matlab.hpp":
        return "matlab"  # all thin wrappers — keep them out of the core tables
    rel = path.relative_to(ROOT)
    top = rel.parts[0] if len(rel.parts) > 1 else ""
    # Grade B/C plant builders (D25) — folder plants/ or legacy filenames
    if top == "plants" or path.name in ("models.hpp", "power_converters.hpp", "plants.hpp"):
        return "plants"
    return top if top in CAT_NAME else "toolbox"


def entity_category(path, kind, name):
    """Folder category, with co-located design Results / free functions retargeted."""
    cat = category(path)
    # Dual headers under controllers/: runtime types stay here; design free functions
    # and *Result structs are the non-PWM-rate shelf (D25 design:: policy).
    if cat == "controllers":
        if kind == "fn":
            return "design"
        if kind in ("type", "enum") and name.endswith("Result"):
            return "design"
    return cat


entries = collections.defaultdict(list)  # cat -> list[(kind,name,brief,header)]
seen = collections.defaultdict(set)
overloads = collections.Counter()  # (cat, name) -> extra same-name doc blocks seen


def detail_spans(text):
    """Character ranges of `namespace detail { ... } // namespace detail` blocks —
    implementation helpers, excluded from the public reference."""
    spans = []
    for m in re.finditer(r"^namespace detail \{", text, re.M):
        close = re.compile(r"^\} // namespace detail", re.M).search(text, m.end())
        spans.append((m.start(), close.end() if close else len(text)))
    return spans


def collect(path, text, header, matches, get_brief, hidden):
    for m in matches:
        if any(a <= m.start(2) < b for a, b in hidden):
            continue  # inside namespace detail
        b = get_brief(m.group(1))
        if not b:
            continue
        d = decl(m.group(2))
        if not d:
            continue
        kind, name, off = d
        cat = entity_category(path, kind, name)
        if name in seen[cat]:
            overloads[(cat, name)] += 1  # keep the first (primary) overload's brief
            continue
        seen[cat].add(name)
        # line of the declaration itself: first line of the code block + its offset
        line = text.count("\n", 0, m.start(2)) + 1 + off
        entries[cat].append((kind, name, b, f"inc/{header}#L{line}"))


for hpp in sorted(ROOT.rglob("*.hpp")):
    if hpp.name in BACKEND_FILES:
        continue  # internal math backends — covered by their own section below
    text = hpp.read_text(encoding="utf-8", errors="replace")
    header = str(hpp.relative_to("inc")).replace("\\", "/")
    hidden = detail_spans(text)
    collect(hpp, text, header, BLOCK.finditer(text), brief_text, hidden)
    collect(hpp, text, header, LINE_DOC.finditer(text), line_doc_text, hidden)

# Same-name doc blocks are overload sets — say so instead of pretending the first
# overload's brief is the whole story (e.g. mat::solve's triangular/LU/Cholesky family).
for cat, rows in entries.items():
    entries[cat] = [(k, n, b + (f" (+{overloads[(cat, n)]} more overload{'s' if overloads[(cat, n)] > 1 else ''})"
                               if overloads[(cat, n)] else ""), link)
                    for k, n, b, link in rows]

out = [
    "# API Reference\n",
    "Auto-generated from `@brief` doc comments in `inc/damp/`. "
    "Regenerate with `python tools/gen_reference.py`. "
    "Flat A→Z view: [REFERENCE_INDEX.md](REFERENCE_INDEX.md).\n",
    "Compile-time (or init-time) control design in the firmware tree "
    "(variant gains as `constexpr`); same runtime objects in SIL. "
    "Scope: [docs/known_limitations.md](docs/known_limitations.md).\n",
    "`damp::design::` is the design-time shelf (not for PWM-rate): heavy solvers and "
    "Result structs. Domain packs (motor / power / full motion kits) are not assumed "
    "present in the public core tree.\n",
]

# Group labels before the first category of each ring (blockquote, not ## — keeps TOC clean).
SECTION_BREAKS = {
    "": "Substrate — types every path uses",
    "design": "Design-time — not for the PWM interrupt",
    "controllers": "Runtime — tick objects and signal processing",
    "trajectory": "Pose / trajectory helpers (core)",
    "toolbox": "Embedded helpers — controls-adjacent utilities (not ETL)",
    "analysis": "Host — analysis, sim, MATLAB® aliases",
}


def table(title, rows):
    if not rows:
        return
    out.append(f"\n**{title}**\n")
    out.append("| Name | Description |")
    out.append("| ---- | ----------- |")
    for _, name, b, link in sorted(rows, key=lambda r: r[1].lower()):
        out.append(f"| [`{name}`]({link}) | {b} |")


for key, _ in CATS:
    cat = entries.get(key)
    if not cat:
        continue
    if key in SECTION_BREAKS:
        out.append(f"\n> **{SECTION_BREAKS[key]}**")
    out.append(f"\n## {CAT_NAME[key]}")
    table("Blocks (structs, classes, enums, concepts)", [r for r in cat if r[0] != "fn"])
    table("Functions", [r for r in cat if r[0] == "fn"])

# --- Math backends --------------------------------------------------------
out.append("\n## Math backends")
out.append("\nInternal, compile-time-selected implementations of the `damp::` scalar-math "
           "surface (`sin`/`cos`/`sqrt`/`exp`/…), chosen via `damp/config.hpp`. Every "
           "backend exposes the same functions, so they're listed once here as files "
           "rather than repeated in the tables above. The public dispatcher is "
           "[`damp/math/math.hpp`](inc/damp/math/math.hpp).\n")
out.append("| File | Role |")
out.append("| ---- | ---- |")
for name in ("math_backend.hpp", "std_fallback.hpp", "damp_backend.hpp",
             "series_backend.hpp", "constexpr_math.hpp", "trig.hpp"):
    f = ROOT / "math" / name
    role = brief_text(f.read_text(encoding="utf-8", errors="replace")[:1500])
    if not role:
        role = f"`{name}` math backend"
    out.append(f"| [`damp/math/{name}`](inc/damp/math/{name}) | {role} |")

# --- Examples -------------------------------------------------------------
ex_rows = []
for cpp in sorted(pathlib.Path("examples").rglob("*.cpp")):
    head = cpp.read_text(encoding="utf-8", errors="replace")[:1500]
    m = BRIEF.search(head)
    desc = brief_text(head) if m else None
    if not desc:  # no @brief: humanize the filename
        stem = cpp.stem.removeprefix("example_").replace("_", " ")
        desc = stem[:1].upper() + stem[1:]
    rel = str(cpp).replace("\\", "/")
    ex_rows.append((cpp.name, rel, desc))

out.append("\n## Examples")
out.append(f"\nRunnable programs in `examples/` ({len(ex_rows)} total). "
           "Build with `make` (or `tup --quiet examples`); outputs go to `examples/build/`.\n")
out.append("| Example | Description |")
out.append("| ------- | ----------- |")
for name, rel, desc in sorted(ex_rows, key=lambda r: r[0].lower()):
    out.append(f"| [`{name}`]({rel}) | {desc} |")

# --- Table of contents (GitHub-slug anchors) ------------------------------
def slug(title):
    s = re.sub(r"[^\w\s-]", "", title.lower())  # drop punctuation, keep spaces/hyphens
    return s.replace(" ", "-")                   # per-space so '&' leaves a double hyphen


toc = ["- [API Reference](#api-reference)"]
for line in out:
    s = line.lstrip()
    if s.startswith("## "):
        t = s[3:].strip()
        toc.append(f"  - [{t.replace('&', chr(92) + '&')}](#{slug(t)})")
# Splice TOC after the intro blurb (before ring labels / first ##).
insert_at = next(
    (i for i, line in enumerate(out)
     if line.lstrip().startswith("## ") or line.lstrip().startswith("> **")),
    len(out),
)
out[insert_at:insert_at] = ["\n" + "\n".join(toc) + "\n"]

pathlib.Path("REFERENCE.md").write_text("\n".join(out) + "\n", encoding="utf-8")
total = sum(len(v) for v in entries.values())
print(f"Wrote REFERENCE.md: {total} entries across {len(entries)} categories")
for key, _ in CATS:
    if entries.get(key):
        print(f"  {CAT_NAME[key]:40s} {len(entries[key])}")

# --- Flat alphabetical index (REFERENCE_INDEX.md) -------------------------
# Same entries, every public symbol A→Z with its domain, for quick lookup. A name
# appearing twice with different links flags a possible canonical-implementation
# violation (D20) — surfaced here on purpose.
KIND_LABEL = {"fn": "function", "type": "block", "enum": "enum", "concept": "concept"}
flat = [(name, KIND_LABEL.get(kind, kind), CAT_NAME[cat], b, link)
        for cat, rows in entries.items() for kind, name, b, link in rows]

idx = ["# API Reference — Alphabetical Index\n",
       "Auto-generated from `@brief` doc comments in `inc/damp/`. "
       "Regenerate with `python tools/gen_reference.py`. "
       "Grouped-by-domain view: [REFERENCE.md](REFERENCE.md).\n",
       "| Name | Kind | Domain | Description |",
       "| ---- | ---- | ------ | ----------- |"]
for name, kind, domain, b, link in sorted(flat, key=lambda r: (r[0].lower(), r[2])):
    idx.append(f"| [`{name}`]({link}) | {kind} | {domain} | {b} |")

pathlib.Path("REFERENCE_INDEX.md").write_text("\n".join(idx) + "\n", encoding="utf-8")
print(f"Wrote REFERENCE_INDEX.md: {len(flat)} entries (flat A-Z)")
