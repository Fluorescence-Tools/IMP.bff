Place the public header files in this directory. They will be
available to your code (and other modules) with

     #include <IMP/bff/myheader.h>

All headers should include `IMP/bff/bff_config.h` as their
first include and surround all code with `IMPBFF_BEGIN_NAMESPACE`
and `IMPBFF_END_NAMESPACE` to put it in the IMP::bff namespace
and manage compiler warnings.

Headers should also be exposed to SWIG in the `pyext/swig.i-in` file.
