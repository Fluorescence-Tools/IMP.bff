#!/usr/bin/env python3
"""Merge related headers or sources into one, flat, the same way every time.

PRD-138: the layout stays flat (IMP's module tooling globs ``include/*.h`` and
``src/*.cpp``), so consolidation means fewer files, not directories. Doing
that by hand a dozen times is a dozen chances to drop an include or leave a
dangling ``%include``; this does it as one reviewable operation::

    utility/merge_sources.py --into include/FCS.h \\
        include/FcsMdf.h include/FcsSaturation.h include/FcsSaturationCurve.h
    utility/merge_sources.py --into src/FCS.cpp \\
        src/FcsMdf.cpp src/FcsSaturation.cpp src/FcsSaturationCurve.cpp

What it does, in order:

1. Appends each source's body to the target under a ``// ---- from X ----``
   marker, so the former file boundaries stay findable by eye.
   For headers the merged files' own include guards are stripped and the body
   goes *inside* the target's guard; the target's guard is kept.
2. Drops ``#include`` lines that name a file being merged away (they would be
   self-includes now) and ones the target already has.
3. Deletes the merged files.
4. Rewrites every ``#include <IMP/bff/X.h>`` / ``"X.h"`` and every SWIG
   ``%include "IMP/bff/X.h"`` under ``include/``, ``src/`` and ``pyext/`` to
   name the target, then removes duplicate include lines that creates.
5. Drops merged ``.cpp`` names from ``src/Files.cmake``.

Order matters and is not left to the caller: the merged bodies are emitted in
**dependency order**, computed from the includes among the files being merged
(a header that includes another goes after it), so a typedef reaches the code
that uses it. The first merge attempted alphabetically put ``FasprUtility.h``
-- the one with the typedefs -- last, and the build stopped at ``IV2``. When
the includes form a cycle (``PathMap.h`` / ``PathMapTile.h`` include each other
and get by on forward declarations) the tool refuses and asks for ``--order``,
which takes the list verbatim; the target itself may then appear in that list
to place an existing header inside the fresh guard.

What it deliberately does not do: reason about symbol collisions. The module
compiles as a unity build (``bff_all.cpp`` includes every source into one
translation unit), so two sources that clashed would already have failed to
build. Merging them into one file cannot introduce a clash that did not exist.

It also does not touch the per-topic ``pyext/IMP_bff.*.i`` files beyond the
``%include`` rewrite: which topic file absorbs which is a judgment per merge,
made in the commit that does it. ``--dry-run`` prints the plan and changes
nothing.
"""
import argparse
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCAN_DIRS = ("include", "src", "pyext")

_GUARD_IFNDEF = re.compile(r"^\s*#ifndef\s+(\w+_H)\s*$", re.M)
_INCLUDE = re.compile(r'^\s*#include\s+[<"]([^>"]+)[>"]\s*$', re.M)
_SWIG_INCLUDE = re.compile(r'^(\s*%include\s+")([^"]+)(")\s*$', re.M)


def _read(p):
    with open(p, encoding="utf-8", errors="replace") as fh:
        return fh.read()


def _write(p, s):
    with open(p, "w", encoding="utf-8") as fh:
        fh.write(s)


def _module_name(path):
    """`include/FcsMdf.h` -> `FcsMdf.h`; `include/internal/CifReader.h` -> `internal/CifReader.h`."""
    rel = os.path.relpath(path, ROOT)
    for d in ("include", "src"):
        if rel.startswith(d + os.sep):
            return rel[len(d) + 1:]
    return rel


def _strip_guard(body):
    """Remove the outer `#ifndef X_H / #define X_H ... #endif` of a header body."""
    m = _GUARD_IFNDEF.search(body)
    if not m:
        return body
    name = m.group(1)
    body = body[:m.start()] + body[m.end():]
    body = re.sub(r"^\s*#define\s+" + re.escape(name) + r"\s*$\n?", "", body, count=1, flags=re.M)
    # the closing #endif is the last one in the file
    last = body.rfind("#endif")
    if last >= 0:
        # take the whole line: `#endif // X_H` and `#endif  /* X_H */` both
        # leave a dangling comment otherwise, and the first merges did
        eol = body.find("\n", last)
        body = body[:last] + (body[eol:] if eol >= 0 else "")
    return body.rstrip() + "\n"


def _dedupe_includes(text, drop_names):
    """Drop includes naming `drop_names`; keep the first of any repeated include."""
    seen = set()
    out = []
    for line in text.split("\n"):
        m = _INCLUDE.match(line)
        if m:
            inc = m.group(1)
            base = inc.split("/")[-1]
            if base in drop_names or inc in drop_names:
                continue
            if inc in seen:
                continue
            seen.add(inc)
        out.append(line)
    return "\n".join(out)


def _topo_order(paths):
    """`paths` in dependency order: a file goes after every set member it includes.

    Edges come from `#include <IMP/bff/X.h>` and `#include "X.h"` lines that
    name another member of the set. Ties keep the given order. A cycle is an
    error: the tool cannot know which forward declaration the authors relied on.
    """
    by_name = {}
    for p in paths:
        by_name[os.path.basename(p)] = p
        by_name[_module_name(p)] = p
    deps = {p: set() for p in paths}
    for p in paths:
        for inc in _INCLUDE.findall(_read(p)):
            q = by_name.get(inc) or by_name.get(inc.split("/")[-1])
            if q and q != p:
                deps[p].add(q)
    order, done = [], set()
    remaining = list(paths)
    while remaining:
        ready = [p for p in remaining if deps[p] <= done]
        if not ready:
            cyc = ", ".join(os.path.relpath(p, ROOT) for p in remaining)
            sys.exit("include cycle among: %s -- pass --order to give the order by hand" % cyc)
        p = ready[0]
        order.append(p)
        done.add(p)
        remaining.remove(p)
    return order


def merge(target, sources, dry_run, explicit_order=False):
    target = os.path.abspath(target)
    sources = [os.path.abspath(s) for s in sources]
    is_header = target.endswith(".h")
    for s in sources:
        if s.endswith(".h") != is_header:
            sys.exit("cannot merge headers and sources into one target: %s" % s)
        if not os.path.exists(s):
            sys.exit("no such file: %s" % s)
    tname = os.path.basename(target)
    target_in_list = target in sources
    if not explicit_order:
        sources = _topo_order(sources)
    merged_names = {os.path.basename(s) for s in sources if s != target}
    merged_mods = {_module_name(s) for s in sources if s != target}
    tmod = _module_name(target)

    # With the target in the ordered list it is rebuilt from scratch: its old
    # body is one of the sections, its guard is stripped like the others, and
    # the file gets a fresh guard around everything.
    body = "" if target_in_list else (_read(target) if os.path.exists(target) else "")
    if is_header and body:
        end = body.rfind("#endif")
        head, tail = body[:end], body[end:]
    elif is_header:
        # A new header gets one guard of its own, in the module's spelling,
        # and the merged bodies go inside it with their guards stripped.
        guard = "IMPBFF_%s_H" % re.sub(r"\W", "_", os.path.splitext(tname)[0]).upper()
        head = "#ifndef %s\n#define %s\n" % (guard, guard)
        tail = "#endif  // %s\n" % guard
    else:
        head, tail = body, ""

    parts = [head.rstrip() + "\n"]
    for s in sources:
        b = _read(s)
        if is_header:
            b = _strip_guard(b)
        if s == target:
            # The target's own leading doc block describes the whole file now,
            # so it goes to the top, ahead of every section.
            m = re.match(r"\s*(/\*\*.*?\*/\n)", b, re.S)
            if m:
                parts.insert(1, "\n" + m.group(1))
                b = b[m.end():]
        else:
            # A `\file X.h` tag in a merged section would name a file that no
            # longer exists; say what it was instead, so Doxygen stays quiet
            # and a reader still learns where the section came from.
            b = re.sub(r"(\*\s+)\\file ((?:IMP/bff/)?[\w/]+\.(?:h|cpp))",
                       r"\1(formerly \2, now a section of this file)", b)
        parts.append("\n// %s from %s %s\n%s" % ("-" * 8, os.path.basename(s), "-" * 8, b.rstrip() + "\n"))
    new = "".join(parts) + ("\n" + tail if tail else "")
    new = _dedupe_includes(new, merged_names | merged_mods | {tname})

    # repo-wide include rewrite
    rewrites = []
    for d in SCAN_DIRS:
        for dp, _, fn in os.walk(os.path.join(ROOT, d)):
            for f in fn:
                if not f.endswith((".h", ".cpp", ".i", ".i-in")):
                    continue
                p = os.path.join(dp, f)
                if p == target or p in sources:
                    continue
                txt = _read(p)
                orig = txt
                for old in merged_mods:
                    txt = re.sub(r'(#include\s+<IMP/bff/)' + re.escape(old) + r'(>)', r'\g<1>' + tmod + r'\g<2>', txt)
                    txt = re.sub(r'(#include\s+")' + re.escape(os.path.basename(old)) + r'(")', r'\g<1>' + tname + r'\g<2>', txt)
                    txt = re.sub(r'(%include\s+"IMP/bff/)' + re.escape(old) + r'(")', r'\g<1>' + tmod + r'\g<2>', txt)
                if txt != orig:
                    txt = _dedupe_includes(txt, set())
                    txt = _dedupe_swig_includes(txt, os.path.relpath(p, ROOT))
                    rewrites.append((p, txt))

    # Files.cmake
    files_cmake = os.path.join(ROOT, "src", "Files.cmake")
    fc_new = None
    if not is_header and os.path.exists(files_cmake):
        fc = _read(files_cmake)
        fc_new = fc
        for n in merged_names:
            fc_new = re.sub(r'(?<![\w/])' + re.escape(n) + r';?', '', fc_new)
        fc_new = fc_new.replace(';;', ';').replace('";', '"').replace(';"', '"')

    print("merge %s <- %s  [%s order]" % (os.path.relpath(target, ROOT),
          ", ".join(os.path.relpath(s, ROOT) for s in sources),
          "given" if explicit_order else "dependency"))
    print("  rewrite includes in %d file(s)" % len(rewrites))
    for p, _ in rewrites:
        print("    ", os.path.relpath(p, ROOT))
    if fc_new is not None and fc_new != fc:
        print("  src/Files.cmake: drop %s" % ", ".join(sorted(merged_names)))
    if dry_run:
        print("  (dry run; nothing written)")
        return
    _write(target, new)
    for s in sources:
        if s != target:
            os.remove(s)
    for p, txt in rewrites:
        _write(p, txt)
    if fc_new is not None and fc_new != fc:
        _write(files_cmake, fc_new)
    print("  done")


def _dedupe_swig_includes(text, where=""):
    """Keep the first `%include` of a path; drop the rest -- and say so.

    SWIG wraps in file order and a header can only be wrapped where every
    type it names is already known. When a merge folds `%include`s that sat
    at different positions into one, keeping the *first* is right only if
    nothing between the old positions needed the later ones first. The tool
    cannot know; it reports the collapsed positions so the merge commit
    decides where the one survivor belongs.
    """
    seen = {}
    out = []
    for n, line in enumerate(text.split("\n"), 1):
        m = _SWIG_INCLUDE.match(line)
        if m:
            path = m.group(2)
            if path in seen:
                seen[path].append(n)
                continue
            seen[path] = [n]
        out.append(line)
    for path, lines in seen.items():
        if len(lines) > 1 and max(lines) - min(lines) > 1:
            print("  NOTE %s: %%include \"%s\" collapsed from lines %s to line %d -- "
                  "check nothing in between needed the later ones first"
                  % (where, path, lines, lines[0]))
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--into", required=True, help="target file (created if absent)")
    ap.add_argument("sources", nargs="+", help="files merged into the target and deleted")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--order", action="store_true",
                    help="emit sources in the order given (needed when their includes form a cycle); "
                         "the target may appear in the list")
    a = ap.parse_args()
    merge(a.into, a.sources, a.dry_run, explicit_order=a.order)


if __name__ == "__main__":
    main()
