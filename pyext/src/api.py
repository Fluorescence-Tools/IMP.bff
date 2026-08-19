"""The public surface of ``IMP.bff``: domain-scoped, with a flat view over it.

The code is organised by **domain** -- ``dye``, ``label``, ``representation``,
``photophysics``, ``sampling``, ``scoring``, ``analysis``, ``observables``,
``io``, ``restraints`` -- and that organisation is what a reader should learn.
A domain is a *name*, not necessarily a directory: most are one flat module,
and the handful that are still packages earned it (see ``pyext/src/README.md``). ``IMP.bff.observables.LifetimeSpectrum``
says where a name lives and what it is about; ``IMP.bff.LifetimeSpectrum`` is a
convenience over it.

So the domains are the authored thing here (:data:`BY_DOMAIN`) and the flat map
(:data:`EXPORTS`) is **derived** from them. That inversion is the point of
PRD-113 stage 7. Before it, the flat map was authored and the grouping existed
only as comments, which meant nothing could check that a name was filed
correctly -- and the test that stood in for it matched export *names* against a
regex of naming families, a list that had to grow with every feature and could
only ever assert that a name looked plausible.

What is checked now instead, in ``test/test_public_api_names.py``:

* every export resolves and carries a docstring;
* every export's module actually lies inside the domain it is filed under, so a
  misfiled name is an error rather than a comment that drifted;
* no name appears in two domains -- the flat view would silently keep one.

Names resolve **lazily** through the ``__getattr__`` that ``pyext/swig.i-in``
installs on the package, so ``import IMP.bff`` imports almost nothing and
nothing here needs ``click``. Retired names are removed, never aliased.

``IMP.bff.cgdye`` is **not** a domain and has no names here. It is explicit
all-atom dye modelling under a force field -- molecular mechanics, which is
IMP's own territory rather than a spectroscopy library's -- and as of
2026-08-18 it sits outside the layout, reachable only by module path. Its 38
flat names were removed. What was *not* molecular mechanics came out of it
first: the rotamer library into ``representation``, attachment into ``label``,
its readers into ``io``. See that package's docstring for the table.
"""

from __future__ import annotations

import importlib

#: domain -> {public name -> the module that defines it}. **The authored map.**
BY_DOMAIN = {
    # -- dye -- the species: spectra, photophysics parameters, topology, the library
    "dye": {
        # IMP.bff.dye
        "read_dye_library": "IMP.bff.dye",
        # IMP.bff.dye
        "available_dyes": "IMP.bff.dye",
        "find_dye": "IMP.bff.dye",
        # IMP.bff.dye
        "Dye": "IMP.bff.dye",
        "Spectrum": "IMP.bff.dye",
        # IMP.bff.dye
        "forster_radius": "IMP.bff.dye",
        "forster_radius_from_spectra": "IMP.bff.dye",
        "spectral_overlap": "IMP.bff.dye",
    },
    # -- io -- the formats: fps.json, its legacy ancestors, structures, templates
    "io": {
        # IMP.bff.io.structure
        "read_dcd": "IMP.bff.io.structure",
        # IMP.bff.io.fps
        "fps_positions_for_docking": "IMP.bff.io.fps",
        "read_fps_json": "IMP.bff.io.fps",
        "write_fps_json": "IMP.bff.io.fps",
        # IMP.bff.io.fps
        "fps_schema_validate": "IMP.bff.io.fps",
        # IMP.bff.io.cif
        "read_nmr_restraints": "IMP.bff.io.cif",
        "write_nmr_restraints": "IMP.bff.io.cif",
        # IMP.bff.io.structure
        "read_rotamer_library_rmf": "IMP.bff.io.structure",
        "write_rotamer_library_rmf": "IMP.bff.io.structure",
    },
    # -- label -- a dye attached at a site, its linker and frame, and the quenchers it sees
    "label": {
        # IMP.bff.label
        "SITE_KEEP_ATOM_NAMES": "IMP.bff.label",
        "align_hierarchies": "IMP.bff.label",
        "attach_dyes": "IMP.bff.label",
        "place_dye_from_coords": "IMP.bff.label",
        "place_dye_from_rotamer_cb": "IMP.bff.label",
        "resolve_dye_site": "IMP.bff.label",
        "strip_sidechain_at_site": "IMP.bff.label",
        # IMP.bff.label
        "backbone_frame": "IMP.bff.label",
        "backbone_frame_from_coords": "IMP.bff.label",
        "backbone_transformation": "IMP.bff.label",
        "backbone_transformation_from_coords": "IMP.bff.label",
        # IMP.bff.label
        "PETParameters": "IMP.bff.label",
        "Quencher": "IMP.bff.label",
        "reference_pet_parameters": "IMP.bff.label",
        "reference_quenchers": "IMP.bff.label",
        # IMP.bff.label
        "FLUOROPHORE_TYPES": "IMP.bff.label",
        "Label": "IMP.bff.label",
        # IMP.bff.label
        "default_strip_mask": "IMP.bff.label",
        "parse_strip_mask": "IMP.bff.label",
        "select_atoms": "IMP.bff.label",
        "site_strip_mask": "IMP.bff.label",
        "strip_hierarchy": "IMP.bff.label",
        "strip_obstacles": "IMP.bff.label",
        "strip_pdb_lines": "IMP.bff.label",
    },
    # -- observables -- stage 4b -- the projection onto an experiment: (amplitude, rate) pairs
    "observables": {
        # IMP.bff.observables
        "lifetime_spectrum_from_rates": "IMP.bff.observables",
        "lifetime_spectrum_from_states": "IMP.bff.observables",
        # IMP.bff.observables
        "LifetimeSpectrum": "IMP.bff.observables",
        "fret_efficiency_from_lifetimes": "IMP.bff.observables",
    },
    # -- photophysics -- interaction terms, rate fields, and the orientation factor
    "photophysics": {
        # IMP.bff.photophysics
        "kappa2_from_dipoles": "IMP.bff.photophysics",
        "kappa2_isotropic": "IMP.bff.photophysics",
        # IMP.bff.photophysics
        "kappa2_distribution_all": "IMP.bff.photophysics",
        "kappa2_distribution_dynamic": "IMP.bff.photophysics",
        "kappa2_isotropic_distribution": "IMP.bff.photophysics",
        "kappa2_order_parameters": "IMP.bff.photophysics",
        "kappa2_to_distance_ratio": "IMP.bff.photophysics",
        # IMP.bff.photophysics
        "FRETTerm": "IMP.bff.photophysics",
        "InteractionTerm": "IMP.bff.photophysics",
        "PETTerm": "IMP.bff.photophysics",
        "RadiativeTerm": "IMP.bff.photophysics",
        "total_rate": "IMP.bff.photophysics",
    },
    # -- quenching -- the PET model, assembled for one labelling site
    "quenching": {
        # IMP.bff.quenching.asa
        "solvent_accessible_surface": "IMP.bff.quenching.asa",
        # IMP.bff.quenching.dynamic
        "DynamicAccessibleVolume": "IMP.bff.quenching.dynamic",
        # IMP.bff.quenching.fret_trace
        "fret_rate_pair_trace": "IMP.bff.quenching.fret_trace",
        "fret_rate_trace": "IMP.bff.quenching.fret_trace",
        # IMP.bff.quenching.grids
        "av_contact_mask": "IMP.bff.quenching.grids",
        "grid_center_index": "IMP.bff.quenching.grids",
        "quenching_rate_grid": "IMP.bff.quenching.grids",
        "slow_factor_grid": "IMP.bff.quenching.grids",
        # IMP.bff.quenching.model
        "atomic_quenching_parameters": "IMP.bff.quenching.model",
        "diffusion_coefficient_map": "IMP.bff.quenching.model",
        "fret_rate_map": "IMP.bff.quenching.model",
        "quenching_rate_map": "IMP.bff.quenching.model",
        "radial_diffusion_map": "IMP.bff.quenching.model",
        "slow_diffusion_near_atoms": "IMP.bff.quenching.model",
        # IMP.bff.quenching.model
        "DyeDiffusionSimulation": "IMP.bff.quenching.model",
        "QuenchedDonorDecay": "IMP.bff.quenching.model",
        # IMP.bff.quenching.pet
        "PET_QUENCHING_REFERENCE": "IMP.bff.quenching.pet",
        "QUENCHER_ATOMS": "IMP.bff.quenching.pet",
        "amino_acid_quenching_defaults": "IMP.bff.quenching.pet",
        "normalize_amino_acid_quenching": "IMP.bff.quenching.pet",
        "quencher_atom_indices": "IMP.bff.quenching.pet",
        "quencher_centers": "IMP.bff.quenching.pet",
        # IMP.bff.quenching.model
        "ResidueSites": "IMP.bff.quenching.model",
        "quench_radii_for_residues": "IMP.bff.quenching.model",
        "quenching_rates_for_residues": "IMP.bff.quenching.model",
        "residue_sites": "IMP.bff.quenching.model",
        "slow_factors_for_residues": "IMP.bff.quenching.model",
    },
    # -- representation -- stage 1 -- where the dye can be: accessible volume, rotamer library, distributions
    "representation": {
        # IMP.bff.representation.av
        "compute_av": "IMP.bff.representation.av",
        "compute_av_from_arrays": "IMP.bff.representation.av",
        "compute_avs_for_structure": "IMP.bff.representation.av",
        # IMP.bff.representation.compare
        "compare_av_and_rotamer_pairs": "IMP.bff.representation.compare",
        "compare_av_and_rotamer_positions": "IMP.bff.representation.compare",
        # IMP.bff.representation.distance
        "av_pair_statistics": "IMP.bff.representation.distance",
        "distance_from_fret_efficiency": "IMP.bff.representation.distance",
        "fret_efficiency": "IMP.bff.representation.distance",
        "fret_pair_distribution": "IMP.bff.representation.distance",
        "fret_pair_efficiencies": "IMP.bff.representation.distance",
        "fret_pair_geometry": "IMP.bff.representation.distance",
        # IMP.bff.representation.distance
        "DyeDistributionNormal": "IMP.bff.representation.distance",
        "LabelDistribution": "IMP.bff.representation.distance",
        "LabelDistributionAV": "IMP.bff.representation.distance",
        # IMP.bff.representation.rotamer
        "RotamerEnsemble": "IMP.bff.representation.rotamer",
        "SIMULATION_TYPE_R1": "IMP.bff.representation.rotamer",
        "rotamer_ensembles_from_fps": "IMP.bff.representation.rotamer",
        # IMP.bff.representation.rotamer
        "RotamerDistance": "IMP.bff.representation.rotamer",
        "RotamerPosition": "IMP.bff.representation.rotamer",
        "distances_from_ensembles": "IMP.bff.representation.rotamer",
        "read_rotamer_fps": "IMP.bff.representation.rotamer",
        "rotamer_ensemble_payload": "IMP.bff.representation.rotamer",
        "rotamer_fret_from_fps": "IMP.bff.representation.rotamer",
        "rotamer_position_payload": "IMP.bff.representation.rotamer",
        "write_rotamer_fps": "IMP.bff.representation.rotamer",
        # IMP.bff.representation.rotamer
        "RotamerFRET": "IMP.bff.representation.rotamer",
        # IMP.bff.representation.rotamer
        "load_protein_frames": "IMP.bff.representation.rotamer",
        "load_rotamer_library": "IMP.bff.representation.rotamer",
        "resolve_rotamer_library_path": "IMP.bff.representation.rotamer",
        "rotamer_library_metadata": "IMP.bff.representation.rotamer",
        "rotamer_library_registry": "IMP.bff.representation.rotamer",
        # IMP.bff.representation.distance
        "States": "IMP.bff.representation.distance",
        # IMP.bff.representation.distance
        "AccessibleVolume": "IMP.bff.representation.distance",
    },
    # -- sampling -- stage 3 -- how configurations are drawn: walks, densities, library screening
    "sampling": {
        # IMP.bff.sampling
        "DyeDiffusionTrajectory": "IMP.bff.sampling",
        "simulate_dye_diffusion": "IMP.bff.sampling",
        # IMP.bff.sampling
        "simulate_photon_trace": "IMP.bff.sampling",
        "simulate_quenched_decay": "IMP.bff.sampling",
        # IMP.bff.sampling
        "apply_rotamer_coordinates": "IMP.bff.sampling",
        "load_rotamer_library_dcd": "IMP.bff.sampling",
        "sample_rotamer_index": "IMP.bff.sampling",
        # IMP.bff.sampling
        "GridDiffusionGradient": "IMP.bff.sampling",
        "GridDiffusionResult": "IMP.bff.sampling",
        "GridDiffusionSolver": "IMP.bff.sampling",
        "diffusion_stability_limit": "IMP.bff.sampling",
        "equilibrium_occupancy": "IMP.bff.sampling",
    },
    # -- scoring -- stage 2 -- is this configuration allowed, and how heavily weighted
    "scoring": {
        # IMP.bff.scoring
        "compute_rotamer_score": "IMP.bff.scoring",
    },
}

#: The flat view: ``IMP.bff.<Name>``. Derived from :data:`BY_DOMAIN`, never
#: edited directly -- a name added to two domains must be an error, not a
#: silently surviving duplicate.
EXPORTS = {}
for _domain, _members in BY_DOMAIN.items():
    for _name, _module in _members.items():
        if _name in EXPORTS:
            raise RuntimeError(
                f"{_name!r} is exported by two domains: "
                f"{EXPORTS[_name]} and {_module}")
        EXPORTS[_name] = _module
del _domain, _members, _name, _module

#: a few exports are attributes with a different name in their module
_SOURCE_NAME = {
    "kappa2_distribution_dynamic": "kappasq_dwt",
    "kappa2_distribution_all": "kappasq_all",
    "kappa2_order_parameters": "s2delta",
    "kappa2_isotropic_distribution": "p_isotropic_orientation_factor",
    "fps_schema_validate": "validate_fps",
    "compare_av_and_rotamer_positions": "compare_positions",
    "compare_av_and_rotamer_pairs": "compare_pairs",
    "make_langevin_simulator": "make_simulator",
    # `IMP.bff.compute_av` has always been the *structure* front door -- a PDB
    # plus an fps position definition -- while `IMP.bff.representation.av`
    # exported the *array* one under the same name. The two lived in different
    # modules of `representation/av/` and were only distinguishable by the route
    # taken to them; merging those modules made the collision visible and forced
    # the second a name of its own. The flat name keeps its meaning.
    "compute_av": "compute_av_from_structure",
    "compute_av_from_arrays": "compute_av",
}


def domain_of(name: str) -> str:
    """Which domain a public name belongs to (raises ``KeyError`` if none)."""
    for domain, members in BY_DOMAIN.items():
        if name in members:
            return domain
    raise KeyError(name)


def resolve(name: str):
    """Import and return the object behind a public name (raises AttributeError)."""
    target = EXPORTS.get(name)
    if target is None:
        raise AttributeError(f"module 'IMP.bff' has no attribute {name!r}")
    module = importlib.import_module(target)
    return getattr(module, _SOURCE_NAME.get(name, name))


def module_getattr(name: str):
    """PEP 562 hook installed on ``IMP.bff``: resolve lazily, then cache."""
    import IMP.bff
    value = resolve(name)
    setattr(IMP.bff, name, value)
    return value


def module_dir():
    """``dir(IMP.bff)`` including the lazy exports."""
    import IMP.bff
    return sorted(set(vars(IMP.bff)) | set(EXPORTS))


def public_names():
    """The sorted list of flat public names."""
    return sorted(EXPORTS)
