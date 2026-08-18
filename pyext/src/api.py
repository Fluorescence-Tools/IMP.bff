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
installs on the package: ``import IMP.bff`` does not import ``cgdye``, and
nothing here needs ``click``. Retired names are removed, never aliased.
"""

from __future__ import annotations

import importlib

#: domain -> {public name -> the module that defines it}. **The authored map.**
BY_DOMAIN = {
    # -- cgdye -- explicit coarse-grained dyes: topology, force field, rotamers, sampling
    "cgdye": {
        # IMP.bff.cgdye.analysis.density
        "analyze_dye_density": "IMP.bff.cgdye.analysis.density",
        # IMP.bff.cgdye.analysis.fret
        "fret_efficiency_regimes": "IMP.bff.cgdye.analysis.fret",
        "fret_efficiency_exact_kinetic": "IMP.bff.cgdye.analysis.fret",
        "fret_efficiency_exact_kinetic_pair": "IMP.bff.cgdye.analysis.fret",
        # IMP.bff.cgdye.io.cif
        "read_dye_forcefield_cif": "IMP.bff.cgdye.io.cif",
        "write_dye_forcefield_cif": "IMP.bff.cgdye.io.cif",
        # IMP.bff.cgdye.io.dcd
        "read_dcd": "IMP.bff.cgdye.io.dcd",
        # IMP.bff.cgdye.io.nmr_cif
        "read_nmr_restraints": "IMP.bff.cgdye.io.nmr_cif",
        "write_nmr_restraints": "IMP.bff.cgdye.io.nmr_cif",
        # IMP.bff.cgdye.io.rotamer_rmf
        "read_rotamer_library_rmf": "IMP.bff.cgdye.io.rotamer_rmf",
        "write_rotamer_library_rmf": "IMP.bff.cgdye.io.rotamer_rmf",
        # IMP.bff.cgdye.io.template_cif
        "read_component_template_cif": "IMP.bff.cgdye.io.template_cif",
        "write_component_template_cif": "IMP.bff.cgdye.io.template_cif",
        "write_dye_template_cif": "IMP.bff.cgdye.io.template_cif",
        # IMP.bff.cgdye.labeling.attachment
        "attach_dyes": "IMP.bff.cgdye.labeling.attachment",
        "resolve_dye_site": "IMP.bff.cgdye.labeling.attachment",
        "place_dye_from_coords": "IMP.bff.cgdye.labeling.attachment",
        "place_dye_from_rotamer_cb": "IMP.bff.cgdye.labeling.attachment",
        "align_hierarchies": "IMP.bff.cgdye.labeling.attachment",
        "strip_sidechain_at_site": "IMP.bff.cgdye.labeling.attachment",
        "SITE_KEEP_ATOM_NAMES": "IMP.bff.cgdye.labeling.attachment",
        # IMP.bff.cgdye.labeling.backbone_frame
        "backbone_frame": "IMP.bff.cgdye.labeling.backbone_frame",
        "backbone_frame_from_coords": "IMP.bff.cgdye.labeling.backbone_frame",
        "backbone_transformation": "IMP.bff.cgdye.labeling.backbone_frame",
        "backbone_transformation_from_coords": "IMP.bff.cgdye.labeling.backbone_frame",
        # IMP.bff.cgdye.rotamer.compare_av
        "compare_av_and_rotamer_positions": "IMP.bff.cgdye.rotamer.compare_av",
        "compare_av_and_rotamer_pairs": "IMP.bff.cgdye.rotamer.compare_av",
        # IMP.bff.cgdye.rotamer.ensemble
        "RotamerEnsemble": "IMP.bff.cgdye.rotamer.ensemble",
        "rotamer_ensembles_from_fps": "IMP.bff.cgdye.rotamer.ensemble",
        "SIMULATION_TYPE_R1": "IMP.bff.cgdye.rotamer.ensemble",
        # IMP.bff.cgdye.rotamer.fps
        "RotamerPosition": "IMP.bff.cgdye.rotamer.fps",
        "RotamerDistance": "IMP.bff.cgdye.rotamer.fps",
        "read_rotamer_fps": "IMP.bff.cgdye.rotamer.fps",
        "write_rotamer_fps": "IMP.bff.cgdye.rotamer.fps",
        "rotamer_position_payload": "IMP.bff.cgdye.rotamer.fps",
        "rotamer_ensemble_payload": "IMP.bff.cgdye.rotamer.fps",
        "distances_from_ensembles": "IMP.bff.cgdye.rotamer.fps",
        "rotamer_fret_from_fps": "IMP.bff.cgdye.rotamer.fps",
        # IMP.bff.cgdye.rotamer.fret
        "RotamerFRET": "IMP.bff.cgdye.rotamer.fret",
        # IMP.bff.cgdye.rotamer.io
        "load_rotamer_library": "IMP.bff.cgdye.rotamer.io",
        "resolve_rotamer_library_path": "IMP.bff.cgdye.rotamer.io",
        "rotamer_library_registry": "IMP.bff.cgdye.rotamer.io",
        "rotamer_library_metadata": "IMP.bff.cgdye.rotamer.io",
        "load_protein_frames": "IMP.bff.cgdye.rotamer.io",
        # IMP.bff.cgdye.rotamer.scoring
        "compute_rotamer_score": "IMP.bff.cgdye.rotamer.scoring",
        # IMP.bff.cgdye.sampling.boltzmann
        "rotamer_cluster_weights": "IMP.bff.cgdye.sampling.boltzmann",
        "boltzmann_weights": "IMP.bff.cgdye.sampling.boltzmann",
        # IMP.bff.cgdye.sampling.clustering
        "cluster_frames_leader": "IMP.bff.cgdye.sampling.clustering",
        "assign_frames_to_clusters": "IMP.bff.cgdye.sampling.clustering",
        # IMP.bff.cgdye.sampling.kinetic
        "rotamer_transition_matrix": "IMP.bff.cgdye.sampling.kinetic",
        "rotamer_correlation_times": "IMP.bff.cgdye.sampling.kinetic",
        "rotamer_rotational_correlation_time": "IMP.bff.cgdye.sampling.kinetic",
        "reconstruct_rotamer_trajectory": "IMP.bff.cgdye.sampling.kinetic",
        # IMP.bff.cgdye.sampling.langevin
        "LangevinDyeSampler": "IMP.bff.cgdye.sampling.langevin",
        "LangevinTrajectory": "IMP.bff.cgdye.sampling.langevin",
        "make_langevin_simulator": "IMP.bff.cgdye.sampling.langevin",
        # IMP.bff.cgdye.sampling.library_gen
        "LinkerSampler": "IMP.bff.cgdye.sampling.library_gen",
        "generate_linker_rotamers": "IMP.bff.cgdye.sampling.library_gen",
        # IMP.bff.cgdye.sampling.mean_field
        "rotamer_mean_field_weights": "IMP.bff.cgdye.sampling.mean_field",
        "rotamer_mean_field_weights_multi_dye": "IMP.bff.cgdye.sampling.mean_field",
        # IMP.bff.cgdye.sampling.rotamer
        "load_rotamer_library_dcd": "IMP.bff.cgdye.sampling.rotamer",
        "apply_rotamer_coordinates": "IMP.bff.cgdye.sampling.rotamer",
        "sample_rotamer_index": "IMP.bff.cgdye.sampling.rotamer",
        # IMP.bff.cgdye.sampling.rrt
        "run_rigid_body_rrt": "IMP.bff.cgdye.sampling.rrt",
        "run_torsion_rrt": "IMP.bff.cgdye.sampling.rrt",
        # IMP.bff.cgdye.sampling.scoring
        "DyeInternalEnergyEvaluator": "IMP.bff.cgdye.sampling.scoring",
        "dye_internal_system": "IMP.bff.cgdye.sampling.scoring",
        # IMP.bff.cgdye.sim.dye_restraints
        "build_dye_restraints": "IMP.bff.cgdye.sim.dye_restraints",
        # IMP.bff.cgdye.system
        "DyeForceFieldSystem": "IMP.bff.cgdye.system",
        # IMP.bff.cgdye.topology.builder
        "parse_dye_mol2": "IMP.bff.cgdye.topology.builder",
        # IMP.bff.cgdye.topology.combined
        "build_dye_protein_system": "IMP.bff.cgdye.topology.combined",
        "dye_forcefield_system": "IMP.bff.cgdye.topology.combined",
        # IMP.bff.cgdye.topology.dye
        "build_dye_topology": "IMP.bff.cgdye.topology.dye",
        "CHARMM36_LJ": "IMP.bff.cgdye.topology.dye",
        "lj_cross": "IMP.bff.cgdye.topology.dye",
        "lj_energy": "IMP.bff.cgdye.topology.dye",
        "torsion_cosine": "IMP.bff.cgdye.topology.dye",
    },
    # -- dynamics -- integrators: Brownian, Smoluchowski, excited-state kMC
    "dynamics": {
        # IMP.bff.dynamics.brownian
        "DyeDiffusionTrajectory": "IMP.bff.dynamics.brownian",
        "simulate_dye_diffusion": "IMP.bff.dynamics.brownian",
        # IMP.bff.dynamics.excited_state
        "simulate_photon_trace": "IMP.bff.dynamics.excited_state",
        "simulate_quenched_decay": "IMP.bff.dynamics.excited_state",
        # IMP.bff.dynamics.smoluchowski
        "GridDiffusionSolver": "IMP.bff.dynamics.smoluchowski",
        "GridDiffusionResult": "IMP.bff.dynamics.smoluchowski",
        "diffusion_stability_limit": "IMP.bff.dynamics.smoluchowski",
        "equilibrium_occupancy": "IMP.bff.dynamics.smoluchowski",
    },
    # -- dye -- the species: spectra, photophysics parameters, topology, the library
    "dye": {
        # IMP.bff.dye.cif
        "read_dye_library": "IMP.bff.dye.cif",
        # IMP.bff.dye.library
        "find_dye": "IMP.bff.dye.library",
        "available_dyes": "IMP.bff.dye.library",
        # IMP.bff.dye.species
        "Dye": "IMP.bff.dye.species",
        "Spectrum": "IMP.bff.dye.species",
        # IMP.bff.dye.spectra
        "forster_radius_from_spectra": "IMP.bff.dye.spectra",
        "forster_radius": "IMP.bff.dye.spectra",
        "spectral_overlap": "IMP.bff.dye.spectra",
    },
    # -- fret -- the FRET engine, distances, and the strip tool
    "fret": {
        # IMP.bff.fret.av
        "compute_av": "IMP.bff.fret.av",
        "compute_avs_for_structure": "IMP.bff.fret.av",
        # IMP.bff.fret.strip
        "parse_strip_mask": "IMP.bff.fret.strip",
        "default_strip_mask": "IMP.bff.fret.strip",
        "site_strip_mask": "IMP.bff.fret.strip",
        "strip_hierarchy": "IMP.bff.fret.strip",
        "strip_pdb_lines": "IMP.bff.fret.strip",
        "strip_obstacles": "IMP.bff.fret.strip",
        "select_atoms": "IMP.bff.fret.strip",
    },
    # -- io -- the formats: fps.json, its legacy ancestors, and structures
    "io": {
        # IMP.bff.io.fps
        "read_fps_json": "IMP.bff.io.fps",
        "write_fps_json": "IMP.bff.io.fps",
        "fps_positions_for_docking": "IMP.bff.io.fps",
        # IMP.bff.io.fps_schema
        "fps_schema_validate": "IMP.bff.io.fps_schema",
    },
    # -- label -- a dye at a site, and the quenchers it sees
    "label": {
        # IMP.bff.label.quencher
        "Quencher": "IMP.bff.label.quencher",
        "reference_quenchers": "IMP.bff.label.quencher",
        "PETParameters": "IMP.bff.label.quencher",
        "reference_pet_parameters": "IMP.bff.label.quencher",
        # IMP.bff.label.site
        "Label": "IMP.bff.label.site",
        "FLUOROPHORE_TYPES": "IMP.bff.label.site",
    },
    # -- observables -- the output contract: (amplitude, rate) pairs and rate constants
    "observables": {
        # IMP.bff.observables.reduce
        "lifetime_spectrum_from_rates": "IMP.bff.observables.reduce",
        "lifetime_spectrum_from_states": "IMP.bff.observables.reduce",
        # IMP.bff.observables.spectrum
        "LifetimeSpectrum": "IMP.bff.observables.spectrum",
        "fret_efficiency_from_lifetimes": "IMP.bff.observables.spectrum",
    },
    # -- photophysics -- interaction terms and the orientation factor
    "photophysics": {
        # IMP.bff.photophysics.kappa2
        "kappa2_from_dipoles": "IMP.bff.photophysics.kappa2",
        "kappa2_isotropic": "IMP.bff.photophysics.kappa2",
        # IMP.bff.photophysics.orientation
        "kappa2_distribution_dynamic": "IMP.bff.photophysics.orientation",
        "kappa2_distribution_all": "IMP.bff.photophysics.orientation",
        "kappa2_order_parameters": "IMP.bff.photophysics.orientation",
        "kappa2_isotropic_distribution": "IMP.bff.photophysics.orientation",
        "kappa2_to_distance_ratio": "IMP.bff.photophysics.orientation",
        # IMP.bff.photophysics.terms
        "InteractionTerm": "IMP.bff.photophysics.terms",
        "RadiativeTerm": "IMP.bff.photophysics.terms",
        "PETTerm": "IMP.bff.photophysics.terms",
        "FRETTerm": "IMP.bff.photophysics.terms",
        "total_rate": "IMP.bff.photophysics.terms",
    },
    # -- quenching -- PET quenching: the field model, the particle model, the maps
    "quenching": {
        # IMP.bff.quenching.asa
        "solvent_accessible_surface": "IMP.bff.quenching.asa",
        # IMP.bff.quenching.dynamic
        "DynamicAccessibleVolume": "IMP.bff.quenching.dynamic",
        # IMP.bff.quenching.fret_trace
        "fret_rate_trace": "IMP.bff.quenching.fret_trace",
        "fret_rate_pair_trace": "IMP.bff.quenching.fret_trace",
        # IMP.bff.quenching.grids
        "grid_center_index": "IMP.bff.quenching.grids",
        "slow_factor_grid": "IMP.bff.quenching.grids",
        "quenching_rate_grid": "IMP.bff.quenching.grids",
        "av_contact_mask": "IMP.bff.quenching.grids",
        # IMP.bff.quenching.maps
        "atomic_quenching_parameters": "IMP.bff.quenching.maps",
        "diffusion_coefficient_map": "IMP.bff.quenching.maps",
        "radial_diffusion_map": "IMP.bff.quenching.maps",
        "slow_diffusion_near_atoms": "IMP.bff.quenching.maps",
        "quenching_rate_map": "IMP.bff.quenching.maps",
        "fret_rate_map": "IMP.bff.quenching.maps",
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
        "residue_sites": "IMP.bff.quenching.sites",
        "slow_factors_for_residues": "IMP.bff.quenching.sites",
        "quenching_rates_for_residues": "IMP.bff.quenching.sites",
        "quench_radii_for_residues": "IMP.bff.quenching.sites",
    },
    # -- representation -- where the dye can be: accessible volumes, rotamers, distributions
    "representation": {
        # IMP.bff.representation.distance
        "av_pair_statistics": "IMP.bff.representation.distance",
        "fret_pair_geometry": "IMP.bff.representation.distance",
        "fret_pair_efficiencies": "IMP.bff.representation.distance",
        "fret_pair_distribution": "IMP.bff.representation.distance",
        "fret_efficiency": "IMP.bff.representation.distance",
        "distance_from_fret_efficiency": "IMP.bff.representation.distance",
        # IMP.bff.representation.distribution
        "LabelDistribution": "IMP.bff.representation.distribution",
        "LabelDistributionAV": "IMP.bff.representation.distribution",
        "DyeDistributionNormal": "IMP.bff.representation.distribution",
        # IMP.bff.representation.states
        "States": "IMP.bff.representation.states",
        # IMP.bff.representation.types
        "AccessibleVolume": "IMP.bff.representation.types",
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
