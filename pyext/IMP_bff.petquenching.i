/*
 * Photoinduced electron transfer: which moieties quench, and how hard.
 *
 * The three values -- `Quencher` (identity), `PETParameters` (the pair's rate)
 * and `ResidueQuenching` (what a diffusing dye sees) -- were frozen dataclasses
 * split across `label.py` and `quenching.py`, with the tables in one and the
 * types that read them in the other. They are one header now, because a rate
 * that takes two partners to define cannot live in only one of them.
 *
 * `quench_radius` and `attenuation_length` are **NaN** where the Python had
 * `None`: NaN is what a C++ double carries, and for the radius it means
 * "inherit the model-wide critical distance" exactly as `None` did.
 *
 * Every function here returns a **plain Python dict** rather than the SWIG map
 * proxy. A table is looked up with `.get()`, iterated, copied and handed back
 * out by half the quenching model, and a proxy that supports `[]` but not
 * `.get()` is a boundary that shows through.
 */

IMP_SWIG_VALUE(IMP::bff, Quencher, Quenchers);
IMP_SWIG_VALUE(IMP::bff, PETParameters, PETParametersList);
IMP_SWIG_VALUE(IMP::bff, ResidueQuenching, ResidueQuenchings);
IMP_SWIG_VALUE(IMP::bff, PETReference, PETReferences);
IMP_SWIG_VALUE(IMP::bff, ResidueSites, ResidueSitesList);

// Two output views from one call: `atomic_quenching_parameters` returns kQ and
// rC together because they are read together and looked up once.
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** out_kQ, int* n_out_kQ)};
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** out_rC, int* n_out_rC)};

// The atoms, straight from numpy.
%apply(double* IN_ARRAY2, int DIM1, int DIM2) {
    (double* coords, int n_atoms, int n_dim)
};

%attribute_py(IMP::bff::Quencher, bool, is_typed, get_is_typed);
%attribute_py(IMP::bff::PETParameters, bool, is_transferred, get_is_transferred);

// The C++ names are taken by the dict-returning surface below.
%rename(_quencher_atoms) IMP::bff::quencher_atoms;
%rename(_pet_quenching_reference) IMP::bff::pet_quenching_reference;
%rename(_reference_quenchers) IMP::bff::reference_quenchers;
%rename(_reference_pet_parameters) IMP::bff::reference_pet_parameters;
%rename(_normalize_amino_acid_quenching) IMP::bff::normalize_amino_acid_quenching;
%rename(_amino_acid_quenching_defaults) IMP::bff::amino_acid_quenching_defaults;
%rename(_slow_factors_for_residues) IMP::bff::slow_factors_for_residues;
%rename(_quenching_rates_for_residues) IMP::bff::quenching_rates_for_residues;
%rename(_quench_radii_for_residues) IMP::bff::quench_radii_for_residues;
%rename(_residue_sites) IMP::bff::residue_sites;

// A build helper over raw pointers; nothing Python-side appends to a result.
%ignore IMP::bff::ResidueSites::add;

%include "IMP/bff/PETQuenching.h"

%extend IMP::bff::ResidueSites {
    %pythoncode %{
        @property
        def slow_centers(self):
            """``(n, 3)`` -- CB, or CA, or the residue's first atom."""
            return _IMP_bff.ResidueSites_get_slow_centers(self).reshape(-1, 3)

        @property
        def quench_centers(self):
            """``(n, 3)`` -- the centroid of the redox-active atoms."""
            return _IMP_bff.ResidueSites_get_quench_centers(self).reshape(-1, 3)

        @property
        def residue_names(self):
            return list(self.get_residue_names())

        def __len__(self):
            return int(self.size())
    %}
}

%template(QuencherMap) std::map<std::string, IMP::bff::Quencher>;
%template(PETParametersMap) std::map<std::string, IMP::bff::PETParameters>;
%template(ResidueQuenchingMap) std::map<std::string, IMP::bff::ResidueQuenching>;
%template(PETReferenceMap) std::map<std::string, IMP::bff::PETReference>;
%template(MapStringVectorString) std::map<std::string, std::vector<std::string> >;
%template(VectorInt2) std::vector<int>;

%pythoncode %{
def _pet_table(table):
    """A per-residue table as `{name: ResidueQuenching}`.

    `None` is the empty table -- the two model defaults spell "use the
    defaults" that way, and a null does not convert to a `std::map`.
    """
    return {} if not table else dict(table)


def quencher_atoms():
    """The redox-active atoms of each residue type, as ``{name: [atoms]}``.

    Not CB: electron transfer happens at the indole ring, the phenol, the
    thioether or the thiol, and stamping a rate on CB puts it up to 4 A from
    the chemistry.
    """
    return {k: list(v) for k, v in _IMP_bff._quencher_atoms().items()}


def pet_quenching_reference():
    """The published PET chemistry, as ``{residue: PETReference}``."""
    return dict(_IMP_bff._pet_quenching_reference())


def reference_quenchers():
    """The redox-active moieties, as ``{residue: Quencher}``."""
    return dict(_IMP_bff._reference_quenchers())


def reference_pet_parameters(dye=REFERENCE_DYE, rate_scale=1.0,
                             attenuation_length=None):
    """The published PET chemistry for one dye, as ``{residue: PETParameters}``.

    ``attenuation_length`` of ``None`` keeps the hard contact-sphere model.
    """
    return dict(_IMP_bff._reference_pet_parameters(
        str(dye), float(rate_scale),
        float("nan") if attenuation_length is None else float(attenuation_length)))


def normalize_amino_acid_quenching(table=None):
    """The full per-residue interaction table, with defaults filled in.

    Accepts a partial ``{residue: ResidueQuenching}``. Unknown residue names are
    kept, so a non-standard residue can be given a rate.
    """
    return dict(_IMP_bff._normalize_amino_acid_quenching(_pet_table(table)))


def amino_acid_quenching_defaults(kQ_scale=1.0, slow_factor=1.0,
                                  dye_radius=DEFAULT_DYE_RADIUS):
    """A full interaction table built from :func:`pet_quenching_reference`."""
    return dict(_IMP_bff._amino_acid_quenching_defaults(
        float(kQ_scale), float(slow_factor), float(dye_radius)))


def slow_factors_for_residues(residue_names, table):
    """The stickiness factor of each residue, in ``residue_names`` order."""
    return np.asarray(_IMP_bff._slow_factors_for_residues(
        [str(r) for r in residue_names], _pet_table(table)), dtype=np.float64)


def quenching_rates_for_residues(residue_names, table):
    """The quenching rate (1/ns) of each residue."""
    return np.asarray(_IMP_bff._quenching_rates_for_residues(
        [str(r) for r in residue_names], _pet_table(table)), dtype=np.float64)


def quench_radii_for_residues(residue_names, table, critical_distance=0.0):
    """The contact radius of each residue, inheriting ``critical_distance``.

    A NaN ``quench_radius`` in the table means "use the model-wide critical
    distance", which is how a project sets one radius for everything and
    overrides it per residue type where it matters.
    """
    return np.asarray(_IMP_bff._quench_radii_for_residues(
        [str(r) for r in residue_names], _pet_table(table),
        float(critical_distance or 0.0)), dtype=np.float64)


def residue_sites(chains, res_ids, res_names, atom_names, coords, table=None):
    """Group atoms by residue and locate its slow and quench centres.

    Residues are keyed by ``(chain, res_id, res_name)``. **Keying on ``res_id``
    alone is wrong** and was a real defect in QuEst: residue numbers restart per
    chain, so in a homodimer every number occurs twice and two residues' atoms
    were folded into one centre.
    """
    return _IMP_bff._residue_sites(
        [str(c) for c in chains], [int(r) for r in res_ids],
        [str(r) for r in res_names], [str(a) for a in atom_names],
        np.ascontiguousarray(np.asarray(coords, dtype=np.float64)),
        _pet_table(table))
%}
