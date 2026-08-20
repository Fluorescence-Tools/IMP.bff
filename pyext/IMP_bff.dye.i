/*
 * The dye as a species.
 *
 * `Dye` and `Spectrum` were frozen dataclasses holding numpy arrays, a CIF
 * reader built on `ihm.format`'s Python parser with a hand-rolled mtime cache,
 * and a Förster radius computed with `np.trapezoid`. None of that is arithmetic
 * that needed C++; what needed it is the same thing `LifetimeSpectrum` needed,
 * which is that the *object* be a C++ value the way `atom`'s and `core`'s are.
 *
 * Everything that reads or computes is C++, including the name-keyed front door
 * the rotamer code uses and the `_cgdye_metadata` read that folds a dye
 * template's dipole atoms into the species. What is left here is the two
 * flrCIF item maps, which are data about a dictionary rather than code.
 */

IMP_SWIG_VALUE(IMP::bff, Spectrum, Spectrums);
IMP_SWIG_VALUE(IMP::bff, Dye, Dyes);

%include "IMP/bff/DyeLibrary.h"

%template(DyeMap) std::map<std::string, IMP::bff::Dye>;

%pythoncode %{
#: flrCIF item for each :class:`Dye` field, or ``None`` where the dictionary has
#: none.
#:
#: **flrCIF's word for a dye is "probe", and for the fluorescent moiety
#: "chromophore".** A *probe* is the labelling reagent
#: (``reactive_probe_name``, e.g. a maleimide) and its *chromophore* is what
#: fluoresces (``chromophore_name``). ``Dye`` is the chromophore plus the
#: photophysics.
#:
#: **Note what is missing**: no dictionary in the stack has an item for quantum
#: yield, extinction coefficient or a spectrum. Those three are bff-native and
#: marked ``None`` deliberately, not by omission -- see the header.
DYE_FLRCIF_ITEMS = {
    "name": "_flr_probe_list.chromophore_name",
    "reactive_probe_name": "_flr_probe_list.reactive_probe_name",
    "probe_origin": "_flr_probe_list.probe_origin",
    "probe_link_type": "_flr_probe_list.probe_link_type",
    "chromophore_center_atom": "_flr_probe_descriptor.chromophore_center_atom",
    "lifetime": "_flr_reference_measurement_lifetime.lifetime",
    # no flrCIF item exists for these
    "spectrum": None,
    "extinction_coefficient": None,
    "quantum_yield": None,
    "dipole_atoms": None,
    "radius": None,
    "hydrodynamic_radius": None,
    "positive_atoms": None,
    "negative_atoms": None,
}

#: flrCIF items for the derivation's inputs and output. ``index_of_refraction``
#: and ``kappa_squared`` are **not** in the upstream IHM-FLR dictionary -- they
#: are added by ``mmfdb_flr_ext.dic``, which is what makes a stored R0
#: reproducible rather than a bare number.
FORSTER_RADIUS_FLRCIF_ITEMS = {
    "forster_radius": "_flr_fret_forster_radius.forster_radius",
    "k2": "_flr_fret_forster_radius.kappa_squared",
    "refractive_index": "_flr_fret_forster_radius.index_of_refraction",
    "donor": "_flr_fret_forster_radius.donor_probe_id",
    "acceptor": "_flr_fret_forster_radius.acceptor_probe_id",
    # no dictionary in the stack has an item for this
    "spectral_overlap": None,
}

%}
