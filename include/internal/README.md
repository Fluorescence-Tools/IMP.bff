Place the private header files in this directory. They will be
available to your code with

     #include <IMP/bff/internal/myheader.h>

All headers should include `IMP/bff/bff_config.h` as their
first include and surround all code with `IMPBFF_BEGIN_INTERNAL_NAMESPACE`
and `IMPBFF_END_INTERNAL_NAMESPACE` to put it in the
IMP::bff::internal namespace and manage compiler warnings.
