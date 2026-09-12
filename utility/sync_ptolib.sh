#!/usr/bin/env sh
# Refresh the source package; include/internal/ptolib.h is a maintained forwarder.
# Usage: utility/sync_ptolib.sh [--copy|--link] [path-to-ptolib]
set -eu
mode=--copy
if [ "${1:-}" = "--link" ]; then mode=--link; shift; fi
if [ "${1:-}" = "--copy" ]; then shift; fi
here=$(cd "$(dirname "$0")/.." && pwd)
src=${1:-$here/../ptolib}
[ -f "$src/scripts/vendor.sh" ] || { echo "no ptolib checkout at $src" >&2; exit 1; }
sh "$src/scripts/vendor.sh" "$mode" "$here/thirdparty/ptolib"
