"""The public surface of ``IMP.bff``: domain-scoped, with a flat view over it.

The code is organised by **domain** -- ``dye``, ``label``, ``representation``,
``photophysics``, ``dynamics``, ``observables``, ``io``, ``restraints`` -- and
that organisation is what a reader should learn. ``IMP.bff.observables.LifetimeSpectrum``
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
        # IMP.bff.dye.cif
        "read_dye_library": "IMP.bff.dye.cif",
        # IMP.bff.dye.library
        "available_dyes": "IMP.bff.dye.library",
        "find_dye": "IMP.bff.dye.library",
        # IMP.bff.dye.species
        "Dye": "IMP.bff.dye.species",
        "Spectrum": "IMP.bff.dye.species",
        # IMP.bff.dye.spectra
        "forster_radius": "IMP.bff.dye.spectra",
        "forster_radius_from_spectra": "IMP.bff.dye.spectra",
        "spectral_overlap": "IMP.bff.dye.spectra",
    },
    # -- io -- the formats: fps.json, its legacy ancestors, structures, templates
    "io": {
        # IMP.bff.io.dcd
        "read_dcd": "IMP.bff.io.dcd",
        # IMP.bff.io.fps
        "fps_positions_for_docking": "IMP.bff.io.fps",
        "read_fps_json": "IMP.bff.io.fps",
        "write_fps_json": "IMP.bff.io.fps",
        # IMP.bff.io.fps_schema
        "fps_schema_validate": "IMP.bff.io.fps_schema",
        # IMP.bff.io.nmr_cif
        "read_nmr_restraints": "IMP.bff.io.nmr_cif",
        "write_nmr_restraints": "IMP.bff.io.nmr_cif",
        # IMP.bff.io.rotamer_rmf
        "read_rotamer_library_rmf": "IMP.bff.io.rotamer_rmf",
        "write_rotamer_library_rmf": "IMP.bff.io.rotamer_rmf",
    },
    # -- label -- a dye attached at a site, its linker and frame, and the quenchers it sees
    "label": {
        # IMP.bff.label.attachment
        "SITE_KEEP_ATOM_NAMES": "IMP.bff.label.attachment",
        "align_hierarchies": "IMP.bff.label.attachment",
        "attach_dyes": "IMP.bff.label.attachment",
        "place_dye_from_coords": "IMP.bff.label.attachment",
        "place_dye_from_rotamer_cb": "IMP.bff.label.attachment",
        "resolve_dye_site": "IMP.bff.label.attachment",
        "strip_sidechain_at_site": "IMP.bff.label.attachment",
        # IMP.bff.label.backbone_frame
        "backbone_frame": "IMP.bff.label.backbone_frame",
        "backbone_frame_from_coords": "IMP.bff.label.backbone_frame",
        "backbone_transformation": "IMP.bff.label.backbone_frame",
        "backbone_transformation_from_coords": "IMP.bff.label.backbone_frame",
        # IMP.bff.label.quencher
        "PETParameters": "IMP.bff.label.quencher",
        "Quencher": "IMP.bff.label.quencher",
        "reference_pet_parameters": "IMP.bff.label.quencher",
        "reference_quenchers": "IMP.bff.label.quencher",
        # IMP.bff.label.site
        "FLUOROPHORE_TYPES": "IMP.bff.label.site",
        "Label": "IMP.bff.label.site",
        # IMP.bff.label.strip
        "default_strip_mask": "IMP.bff.label.strip",
        "parse_strip_mask": "IMP.bff.label.strip",
        "select_atoms": "IMP.bff.label.strip",
        "site_strip_mask": "IMP.bff.label.strip",
        "strip_hierarchy": "IMP.bff.label.strip",
        "strip_obstacles": "IMP.bff.label.strip",
        "strip_pdb_lines": "IMP.bff.label.strip",
    },
    # -- observables -- stage 4b -- the projection onto an experiment: (amplitude, rate) pairs
    "observables": {
        # IMP.bff.observables.reduce
        "lifetime_spectrum_from_rates": "IMP.bff.observables.reduce",
        "lifetime_spectrum_from_states": "IMP.bff.observables.reduce",
        # IMP.bff.observables.spectrum
        "LifetimeSpectrum": "IMP.bff.observables.spectrum",
        "fret_efficiency_from_lifetimes": "IMP.bff.observables.spectrum",
    },
    # -- photophysics -- interaction terms, rate fields, and the orientation factor
    "photophysics": {
        # IMP.bff.photophysics.kappa2
        "kappa2_from_dipoles": "IMP.bff.photophysics.kappa2",
        "kappa2_isotropic": "IMP.bff.photophysics.kappa2",
        # IMP.bff.photophysics.orientation
        "kappa2_distribution_all": "IMP.bff.photophysics.orientation",
        "kappa2_distribution_dynamic": "IMP.bff.photophysics.orientation",
        "kappa2_isotropic_distribution": "IMP.bff.photophysics.orientation",
        "kappa2_order_parameters": "IMP.bff.photophysics.orientation",
        "kappa2_to_distance_ratio": "IMP.bff.photophysics.orientation",
        # IMP.bff.photophysics.terms
        "FRETTerm": "IMP.bff.photophysics.terms",
        "InteractionTerm": "IMP.bff.photophysics.terms",
        "PETTerm": "IMP.bff.photophysics.terms",
        "RadiativeTerm": "IMP.bff.photophysics.terms",
        "total_rate": "IMP.bff.photophysics.terms",
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
        # IMP.bff.quenching.maps
        "atomic_quenching_parameters": "IMP.bff.quenching.maps",
        "diffusion_coefficient_map": "IMP.bff.quenching.maps",
        "fret_rate_map": "IMP.bff.quenching.maps",
        "quenching_rate_map": "IMP.bff.quenching.maps",
        "radial_diffusion_map": "IMP.bff.quenching.maps",
        "slow_diffusion_near_atoms": "IMP.bff.quenching.maps",
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
        # IMP.bff.quenching.sites
        "ResidueSites": "IMP.bff.quenching.sites",
        "quench_radii_for_residues": "IMP.bff.quenching.sites",
        "quenching_rates_for_residues": "IMP.bff.quenching.sites",
        "residue_sites": "IMP.bff.quenching.sites",
        "slow_factors_for_residues": "IMP.bff.quenching.sites",
    },
    # -- representation -- stage 1 -- where the dye can be: accessible volume, rotamer library, distributions
    "representation": {
        # IMP.bff.representation.av.structure
        "compute_av": "IMP.bff.representation.av.structure",
        "compute_avs_for_structure": "IMP.bff.representation.av.structure",
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
        # IMP.bff.representation.label_distribution
        "DyeDistributionNormal": "IMP.bff.representation.label_distribution",
        "LabelDistribution": "IMP.bff.representation.label_distribution",
        "LabelDistributionAV": "IMP.bff.representation.label_distribution",
        # IMP.bff.representation.rotamer.ensemble
        "RotamerEnsemble": "IMP.bff.representation.rotamer.ensemble",
        "SIMULATION_TYPE_R1": "IMP.bff.representation.rotamer.ensemble",
        "rotamer_ensembles_from_fps": "IMP.bff.representation.rotamer.ensemble",
        # IMP.bff.representation.rotamer.fps
        "RotamerDistance": "IMP.bff.representation.rotamer.fps",
        "RotamerPosition": "IMP.bff.representation.rotamer.fps",
        "distances_from_ensembles": "IMP.bff.representation.rotamer.fps",
        "read_rotamer_fps": "IMP.bff.representation.rotamer.fps",
        "rotamer_ensemble_payload": "IMP.bff.representation.rotamer.fps",
        "rotamer_fret_from_fps": "IMP.bff.representation.rotamer.fps",
        "rotamer_position_payload": "IMP.bff.representation.rotamer.fps",
        "write_rotamer_fps": "IMP.bff.representation.rotamer.fps",
        # IMP.bff.representation.rotamer.fret
        "RotamerFRET": "IMP.bff.representation.rotamer.fret",
        # IMP.bff.representation.rotamer.io
        "load_protein_frames": "IMP.bff.representation.rotamer.io",
        "load_rotamer_library": "IMP.bff.representation.rotamer.io",
        "resolve_rotamer_library_path": "IMP.bff.representation.rotamer.io",
        "rotamer_library_metadata": "IMP.bff.representation.rotamer.io",
        "rotamer_library_registry": "IMP.bff.representation.rotamer.io",
        # IMP.bff.representation.states
        "States": "IMP.bff.representation.states",
        # IMP.bff.representation.types
        "AccessibleVolume": "IMP.bff.representation.types",
    },
    # -- sampling -- stage 3 -- how configurations are drawn: walks, densities, library screening
    "sampling": {
        # IMP.bff.sampling.brownian
        "DyeDiffusionTrajectory": "IMP.bff.sampling.brownian",
        "simulate_dye_diffusion": "IMP.bff.sampling.brownian",
        # IMP.bff.sampling.excited_state
        "simulate_photon_trace": "IMP.bff.sampling.excited_state",
        "simulate_quenched_decay": "IMP.bff.sampling.excited_state",
        # IMP.bff.sampling.rotamer_library
        "apply_rotamer_coordinates": "IMP.bff.sampling.rotamer_library",
        "load_rotamer_library_dcd": "IMP.bff.sampling.rotamer_library",
        "sample_rotamer_index": "IMP.bff.sampling.rotamer_library",
        # IMP.bff.sampling.smoluchowski
        "GridDiffusionGradient": "IMP.bff.sampling.smoluchowski",
        "GridDiffusionResult": "IMP.bff.sampling.smoluchowski",
        "GridDiffusionSolver": "IMP.bff.sampling.smoluchowski",
        "diffusion_stability_limit": "IMP.bff.sampling.smoluchowski",
        "equilibrium_occupancy": "IMP.bff.sampling.smoluchowski",
    },
    # -- scoring -- stage 2 -- is this configuration allowed, and how heavily weighted
    "scoring": {
        # IMP.bff.scoring.rotamer
        "compute_rotamer_score": "IMP.bff.scoring.rotamer",
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
    "fps_schema_validate": "validate",
    "compare_av_and_rotamer_positions": "compare_positions",
    "compare_av_and_rotamer_pairs": "compare_pairs",
    "make_langevin_simulator": "make_simulator",
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
