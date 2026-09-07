#!/usr/bin/env python3
"""Vendor FASPR 1:1 into imp.bff: flat src/Faspr*.h + src/Faspr*.cpp.

Mechanical transforms only (numerics untouched):
  * wrap declarations/definitions in namespace IMP::bff::faspr
  * prefix include guards (FASPR_*) and sibling includes ("FasprX.h")
  * drop MSVC #pragma warning lines
  * exit(0) -> faspr_fail()  (throws; defined in Utility)
  * sprintf -> snprintf (macOS deprecation, identical behaviour)

FASPR.cpp (the CLI main) is not ported; its pipeline becomes src/Faspr.cpp.
"""
import re
import sys
from pathlib import Path

SRC = Path("/Users/tpeulen/dev/imp.bff/junk/FASPR/src")
DST_H = Path("/Users/tpeulen/dev/imp.bff/src")
DST_C = Path("/Users/tpeulen/dev/imp.bff/src")
DST_H.mkdir(parents=True, exist_ok=True)

UNITS = ["Utility", "AAName", "Structure", "RotamerBuilder",
         "SelfEnergy", "PairEnergy", "Search"]

BANNER = """/* Vendored 1:1 from FASPR (github.com/tommyhuangthu/FASPR, version
 * 20200309) into IMP.bff by tools/port_faspr.py -- behaviour-identical port
 * wrapped in namespace IMP::bff::faspr. Only mechanical changes: namespace,
 * include guards/prefixes, exit() -> throwing faspr_fail(), MSVC pragmas
 * dropped, sprintf -> snprintf. See IMP/bff/Faspr.h for the license
 * preserved below and the citation. */
"""

VENDORED_NOTE = "/* FASPR original header follows verbatim. */\n"

NS_OPEN = "\nnamespace IMP {\nnamespace bff {\nnamespace faspr {\n\n"
NS_CLOSE = "\n} // namespace faspr\n} // namespace bff\n} // namespace IMP\n"


def last_include_idx(lines):
    last = -1
    for i, ln in enumerate(lines[:40]):
        if ln.startswith("#include"):
            last = i
    return last


def transform(text, is_header):
    lines = text.splitlines(keepends=True)
    out = []
    for ln in lines:
        s = ln
        if s.strip().startswith("#pragma warning"):
            continue
        # sibling include prefix -- in .cpp only. Headers must keep plain
        # "X.h" so file-relative lookup finds them without src/ on the
        # include path (the amalgamated bff_all.cpp build has no -Isrc).
        if not is_header:
            for u in UNITS:
                s = s.replace(f'#include "{u}.h"', f'#include "Faspr{u}.h"')
        # guard prefix (headers)
        if is_header:
            s = re.sub(r"#ifndef (\w+_H)$", r"#ifndef FASPR_\1", s, flags=re.M)
            s = re.sub(r"#define (\w+_H)$", r"#define FASPR_\1", s, flags=re.M)
        out.append(s)
    lines = out
    # exit(0) -> faspr_fail()
    lines = [("        faspr_fail();\n" if ln.strip() == "exit(0);" else ln)
             for ln in lines]
    # inline `cerr<<...;exit(0);` compound statements
    joined = "".join(lines)
    joined = joined.replace(";exit(0);", ";faspr_fail();")
    lines = joined.splitlines(keepends=True)
    # unqualified swap is hidden by IMP::bff-scope declarations when the
    # amalgamated bff_all.cpp build pulls IMP headers in first
    lines = [ln.replace("swap(i1,i2)", "std::swap(i1,i2)") for ln in lines]
    # sprintf -> snprintf
    lines = [ln.replace("sprintf(inpath,", "snprintf(inpath, sizeof(inpath),")
             for ln in lines]
    # namespace wrap after the last #include of the prologue
    k = last_include_idx(lines)
    lines.insert(k + 1, NS_OPEN)
    text = "".join(lines)
    if is_header:
        # close before the final #endif
        idx = text.rstrip().rfind("#endif")
        text = text[:idx] + NS_CLOSE + text[idx:]
    else:
        text = text.rstrip() + "\n" + NS_CLOSE
    return text


for u in UNITS:
    orig = (SRC / f"{u}.h").read_text()
    (DST_H / f"{u}.h").write_text(BANNER + VENDORED_NOTE + transform(orig, True))
    orig = (SRC / f"{u}.cpp").read_text()
    (DST_C / f"Faspr{u}.cpp").write_text(BANNER + VENDORED_NOTE + transform(orig, False))
    print("ported", u)

# faspr_fail declaration into Utility.h (before namespace close)
uh = DST_H / "Utility.h"
t = uh.read_text()
t = t.replace("#include <algorithm>", "#include <algorithm>\n#include <stdexcept>")
t = t.replace(NS_CLOSE,
              "\n//! FASPR's exit(0) error paths become exceptions in-process.\n"
              "[[noreturn]] inline void faspr_fail() {\n"
              "  throw std::runtime_error(\n"
              "      \"IMP.bff FASPR port: fatal error (see stderr)\");\n}\n"
              + NS_CLOSE)
uh.write_text(t)
print("faspr_fail injected into Utility.h")
