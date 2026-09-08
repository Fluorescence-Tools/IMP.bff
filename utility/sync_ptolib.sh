#!/usr/bin/env sh
# Refresh include/internal/ptolib.h from a sibling ptolib checkout: a verbatim
# copy (what is committed -- a clone has no checkout beside it), or a symlink
# for a working tree while editing both sides.
#
#   utility/sync_ptolib.sh [path-to-ptolib]         # a copy (default ../ptolib)
#   utility/sync_ptolib.sh --link [path-to-ptolib]  # a symlink, never committed
#
# ptolib (https://github.com/tpeulen/ptolib) is the source of truth for the PTO
# container and the DataStore; test/test_vendored_headers.py compares the file
# against the checkout either way. One-way: fix upstream.
set -eu
mode=copy
if [ "${1:-}" = "--link" ]; then mode=link; shift; fi
if [ "${1:-}" = "--copy" ]; then mode=copy; shift; fi
here=$(cd "$(dirname "$0")/.." && pwd)
src=${1:-$here/../ptolib}
[ -f "$src/include/ptolib/ptolib.h" ] || { echo "no ptolib checkout at $src" >&2; exit 1; }
rm -f "$here/include/internal/ptolib.h"
if [ "$mode" = link ]; then
  rel=$(python3 -c "import os,sys; print(os.path.relpath(sys.argv[1], sys.argv[2]))" "$src/include/ptolib/ptolib.h" "$here/include/internal")
  ln -s "$rel" "$here/include/internal/ptolib.h"
else
  cp "$src/include/ptolib/ptolib.h" "$here/include/internal/ptolib.h"
fi
tag=$(git -C "$src" describe --tags --always --dirty 2>/dev/null || echo unknown)
echo "include/internal/ptolib.h <- ptolib $tag ($mode)"
