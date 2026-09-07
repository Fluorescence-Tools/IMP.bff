python-ihm's C mmCIF/BinaryCIF parser (`ihm_format.c/h`, with the `cmp`
MessagePack library it uses), vendored verbatim from IMP's copy of
python-ihm (`modules/core/dependency/python-ihm/src`, MIT, LICENSE beside
this file) for the standalone build. The IMP build compiles against IMP's
copy; `include/internal/Cif.h` includes `"ihm_format.h"` unqualified, and
each build puts the copy it wants on the include path. Refresh from the
same place; do not edit.
