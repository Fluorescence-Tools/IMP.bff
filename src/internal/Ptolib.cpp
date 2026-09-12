/** \file Ptolib.cpp
 * Vendored PTO implementation, compiled once for IMP.bff.
 * Include the single header directly so its implementation remains reachable
 * when declarations have already been included in IMP's unity build.
 */
#define PTOLIB_JSON_INCLUDE <IMP/bff/internal/json.h>
#define PTOLIB_IMPLEMENTATION
#include <IMP/bff/internal/ptolib.h>
#undef PTOLIB_IMPLEMENTATION
