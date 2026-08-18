Place the private header files in this directory. They will be
available to your code with

     #include <IMP/bff/internal/myheader.h>

All headers should include `IMP/bff/bff_config.h` as their
first include and surround all code with `IMPBFF_BEGIN_INTERNAL_NAMESPACE`
and `IMPBFF_END_INTERNAL_NAMESPACE` to put it in the
IMP::bff::internal namespace and manage compiler warnings.

Vendored copies live here too and are the exception to the rule above: they
keep their upstream namespaces and include guards so they can be refreshed by
a plain `cp`. `pcg_*.h` (PCG random numbers), `json.h` (nlohmann/json), and
`MlpCore.h` — the differentiable MLP kernels shared with tttrlib
(`../tttrlib/modules/math/include/MlpCore.h` is the source; never edit the copy;
`test/test_vendored_mlpcore.py` fails when the two diverge). `MlpCore.h`
honours `TTTRLIB_MLPCORE_NAMESPACE`, so a bff header may define it as
`IMP::bff::internal` before including if the `tttrlib` namespace is unwanted.
