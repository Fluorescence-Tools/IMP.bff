/*
 * The dye as a species.
 *
 * `Dye` and `Spectrum` were frozen dataclasses holding numpy arrays, a CIF
 * reader built on `ihm.format`'s Python parser with a hand-rolled mtime cache,
 * and a Förster radius computed with `np.trapezoid`. None of that is arithmetic
 * that needed C++; what needed it is the same thing `LifetimeSpectrum` needed,
 * which is that the *object* be a C++ value the way `atom`'s and `core`'s are.
 *
 * Everything that reads or computes is C++ -- the name-keyed front door the
 * rotamer code uses, the `_cgdye_metadata` read, and the two flrCIF item maps,
 * which move into DyeLibrary.cpp as plain data exposed over SWIG.
 */

IMP_SWIG_VALUE(IMP::bff, Spectrum, Spectrums);
IMP_SWIG_VALUE(IMP::bff, Dye, Dyes);

%include "IMP/bff/DyeLibrary.h"

%template(DyeMap) std::map<std::string, IMP::bff::Dye>;

// The two flrCIF item maps are data about a dictionary, now returned from C++.
// They declare themselves in DyeLibrary.h; over SWIG they surface as Python
// dicts via the MapStringString template. The absent-item keys (a bff-native
// category with no flrCIF item) simply have no key, where the Python once
// mapped them to `None`.
//   DYE_FLRCIF_ITEMS                -> dye_flrcif_items()
//   FORSTER_RADIUS_FLRCIF_ITEMS     -> forster_radius_flrcif_items()
