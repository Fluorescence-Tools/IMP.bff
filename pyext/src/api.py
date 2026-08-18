"""The flat, user-facing surface of ``IMP.bff`` (PRD-107).

Every public name of the explicit-dye (cgdye) and FRET (fret) Python layers
is reachable directly as ``IMP.bff.<Name>``; the sub-packages
(``IMP.bff.cgdye.*``, ``IMP.bff.fret.*``) are where the code is *organised*,
not what a user has to remember. Names are resolved lazily -- ``import
IMP.bff`` does not import cgdye, and nothing here needs ``click`` -- through
the module-level ``__getattr__`` that ``pyext/swig.i-in`` installs.

Naming rule (asserted by ``test/test_public_api_names.py``): every export has a
docstring and says what it is -- ``Rotamer*``, ``Dye*``, ``Linker*``,
``Langevin*``, ``AccessibleVolume``, ``fret_*``, ``fps_*``, ``forster_*``,
``kappa2_*``, ``strip_*``, ``rotamer_*``, ``dye_*``, ``linker_*``, ... No
aliases of retired names are kept.
"""

from __future__ import annotations

import importlib

#: name -> module that defines it. Grouped by family; keep sorted within a group.
EXPORTS = {
    # --- labelling: put an explicit dye on a structure ----------------------
    "attach_dyes": "IMP.bff.cgdye.labeling.attachment",
    "resolve_dye_site": "IMP.bff.cgdye.labeling.attachment",
    "place_dye_from_coords": "IMP.bff.cgdye.labeling.attachment",
    "place_dye_from_rotamer_cb": "IMP.bff.cgdye.labeling.attachment",
    "align_hierarchies": "IMP.bff.cgdye.labeling.attachment",
    "strip_sidechain_at_site": "IMP.bff.cgdye.labeling.attachment",
    "SITE_KEEP_ATOM_NAMES": "IMP.bff.cgdye.labeling.attachment",
    "backbone_frame": "IMP.bff.cgdye.labeling.backbone_frame",
    "backbone_frame_from_coords": "IMP.bff.cgdye.labeling.backbone_frame",
    "backbone_transformation": "IMP.bff.cgdye.labeling.backbone_frame",
    "backbone_transformation_from_coords": "IMP.bff.cgdye.labeling.backbone_frame",
    # --- dye topology / force field ------------------------------------------
    "DyeForceFieldSystem": "IMP.bff.cgdye.system",
    "parse_dye_mol2": "IMP.bff.cgdye.topology.builder",
    "build_dye_topology": "IMP.bff.cgdye.topology.dye",
    "build_dye_protein_system": "IMP.bff.cgdye.topology.combined",
    "build_dye_restraints": "IMP.bff.cgdye.sim.dye_restraints",
    "CHARMM36_LJ": "IMP.bff.cgdye.topology.dye",
    "lj_cross": "IMP.bff.cgdye.topology.dye",
    "lj_energy": "IMP.bff.cgdye.topology.dye",
    "DyeInternalEnergyEvaluator": "IMP.bff.cgdye.sampling.scoring",
    "dye_internal_system": "IMP.bff.cgdye.sampling.scoring",
    "read_dye_forcefield_cif": "IMP.bff.cgdye.io.cif",
    "write_dye_forcefield_cif": "IMP.bff.cgdye.io.cif",
    "read_component_template_cif": "IMP.bff.cgdye.io.template_cif",
    "write_component_template_cif": "IMP.bff.cgdye.io.template_cif",
    "write_dye_template_cif": "IMP.bff.cgdye.io.template_cif",
    # --- rotamer libraries and rotamer FRET ----------------------------------
    "RotamerFRET": "IMP.bff.cgdye.rotamer.fret",
    "RotamerEnsemble": "IMP.bff.cgdye.rotamer.ensemble",
    "rotamer_ensembles_from_fps": "IMP.bff.cgdye.rotamer.ensemble",
    "compare_av_and_rotamer_positions": "IMP.bff.cgdye.rotamer.compare_av",
    "compare_av_and_rotamer_pairs": "IMP.bff.cgdye.rotamer.compare_av",
    "SIMULATION_TYPE_R1": "IMP.bff.cgdye.rotamer.ensemble",
    "RotamerPosition": "IMP.bff.cgdye.rotamer.fps",
    "RotamerDistance": "IMP.bff.cgdye.rotamer.fps",
    "read_rotamer_fps": "IMP.bff.cgdye.rotamer.fps",
    "write_rotamer_fps": "IMP.bff.cgdye.rotamer.fps",
    "rotamer_position_payload": "IMP.bff.cgdye.rotamer.fps",
    "rotamer_ensemble_payload": "IMP.bff.cgdye.rotamer.fps",
    "distances_from_ensembles": "IMP.bff.cgdye.rotamer.fps",
    "rotamer_fret_from_fps": "IMP.bff.cgdye.rotamer.fps",
    "load_rotamer_library": "IMP.bff.cgdye.rotamer.io",
    "resolve_rotamer_library_path": "IMP.bff.cgdye.rotamer.io",
    "rotamer_library_registry": "IMP.bff.cgdye.rotamer.io",
    "rotamer_library_metadata": "IMP.bff.cgdye.rotamer.io",
    "load_protein_frames": "IMP.bff.cgdye.rotamer.io",
    "load_rotamer_library_dcd": "IMP.bff.cgdye.sampling.rotamer",
    "read_rotamer_library_rmf": "IMP.bff.cgdye.io.rotamer_rmf",
    "write_rotamer_library_rmf": "IMP.bff.cgdye.io.rotamer_rmf",
    "apply_rotamer_coordinates": "IMP.bff.cgdye.sampling.rotamer",
    "sample_rotamer_index": "IMP.bff.cgdye.sampling.rotamer",
    "compute_rotamer_score": "IMP.bff.cgdye.rotamer.scoring",
    "rotamer_mean_field_weights": "IMP.bff.cgdye.sampling.mean_field",
    "rotamer_mean_field_weights_multi_dye": "IMP.bff.cgdye.sampling.mean_field",
    "rotamer_cluster_weights": "IMP.bff.cgdye.sampling.boltzmann",
    "rotamer_transition_matrix": "IMP.bff.cgdye.sampling.kinetic",
    "rotamer_correlation_times": "IMP.bff.cgdye.sampling.kinetic",
    "rotamer_rotational_correlation_time": "IMP.bff.cgdye.sampling.kinetic",
    "reconstruct_rotamer_trajectory": "IMP.bff.cgdye.sampling.kinetic",
    # --- linker sampling ------------------------------------------------------
    "LinkerSampler": "IMP.bff.cgdye.sampling.library_gen",
    "generate_linker_rotamers": "IMP.bff.cgdye.sampling.library_gen",
    "boltzmann_weights": "IMP.bff.cgdye.sampling.boltzmann",
    "cluster_frames_leader": "IMP.bff.cgdye.sampling.clustering",
    "assign_frames_to_clusters": "IMP.bff.cgdye.sampling.clustering",
    "run_rigid_body_rrt": "IMP.bff.cgdye.sampling.rrt",
    "LangevinDyeSampler": "IMP.bff.cgdye.sampling.langevin",
    "LangevinTrajectory": "IMP.bff.cgdye.sampling.langevin",
    "make_langevin_simulator": "IMP.bff.cgdye.sampling.langevin",
    "dye_forcefield_system": "IMP.bff.cgdye.topology.combined",
    "torsion_cosine": "IMP.bff.cgdye.topology.dye",
    "run_torsion_rrt": "IMP.bff.cgdye.sampling.rrt",
    # --- FRET efficiencies from ensembles and kinetics ------------------------
    "fret_efficiency_regimes": "IMP.bff.cgdye.analysis.fret",
    "fret_efficiency_exact_kinetic": "IMP.bff.cgdye.analysis.fret",
    "fret_efficiency_exact_kinetic_pair": "IMP.bff.cgdye.analysis.fret",
    "analyze_dye_density": "IMP.bff.cgdye.analysis.density",
    # --- label-pair physics (fret) --------------------------------------------
    "forster_radius_from_spectra": "IMP.bff.dye.spectra",
    "forster_radius": "IMP.bff.dye.spectra",
    "spectral_overlap": "IMP.bff.dye.spectra",
    "Dye": "IMP.bff.dye.species",
    "Spectrum": "IMP.bff.dye.species",
    "find_dye": "IMP.bff.dye.library",
    # --- system: what is attached where (PRD-113 stage 2) --------------------
    "Label": "IMP.bff.label.site",
    "FLUOROPHORE_TYPES": "IMP.bff.label.site",
    "Quencher": "IMP.bff.label.quencher",
    "reference_quenchers": "IMP.bff.label.quencher",
    "PETParameters": "IMP.bff.label.quencher",
    "reference_pet_parameters": "IMP.bff.label.quencher",
    "available_dyes": "IMP.bff.dye.library",
    "read_dye_library": "IMP.bff.dye.cif",
    "kappa2_from_dipoles": "IMP.bff.fret.kappa2",
    "kappa2_isotropic": "IMP.bff.fret.kappa2",
    # --- accessible volumes, fps.json, distances (fret) -----------------------
    "AccessibleVolume": "IMP.bff.representation.types",
    "States": "IMP.bff.representation.states",
    "compute_av": "IMP.bff.fret.av",
    "compute_avs_for_structure": "IMP.bff.fret.av",
    "read_fps_json": "IMP.bff.fret.io",
    "write_fps_json": "IMP.bff.fret.io",
    "fps_schema_validate": "IMP.bff.fret.fps_schema",
    "fps_positions_for_docking": "IMP.bff.fret.io",
    "av_pair_statistics": "IMP.bff.fret.distance",
    "fret_pair_geometry": "IMP.bff.fret.distance",
    "fret_pair_efficiencies": "IMP.bff.fret.distance",
    "fret_pair_distribution": "IMP.bff.fret.distance",
    "fret_efficiency": "IMP.bff.fret.distance",
    "distance_from_fret_efficiency": "IMP.bff.fret.distance",
    # --- PET quenching of a dye diffusing in its AV (quenching) ---------------
    "PET_QUENCHING_REFERENCE": "IMP.bff.quenching.pet",
    "QUENCHER_ATOMS": "IMP.bff.quenching.pet",
    "amino_acid_quenching_defaults": "IMP.bff.quenching.pet",
    "normalize_amino_acid_quenching": "IMP.bff.quenching.pet",
    "quencher_atom_indices": "IMP.bff.quenching.pet",
    "quencher_centers": "IMP.bff.quenching.pet",
    "sphere_points": "IMP.bff.quenching.asa",
    "solvent_accessible_surface": "IMP.bff.quenching.asa",
    "grid_center_index": "IMP.bff.quenching.grids",
    "slow_factor_grid": "IMP.bff.quenching.grids",
    "quenching_rate_grid": "IMP.bff.quenching.grids",
    "av_contact_mask": "IMP.bff.quenching.grids",
    "DyeDiffusionTrajectory": "IMP.bff.quenching.diffusion",
    "simulate_dye_diffusion": "IMP.bff.quenching.diffusion",
    "simulate_photon_trace": "IMP.bff.quenching.photon",
    "simulate_quenched_decay": "IMP.bff.quenching.photon",
    "fret_rate_trace": "IMP.bff.quenching.fret_trace",
    "fret_rate_pair_trace": "IMP.bff.quenching.fret_trace",
    # --- the field formulation: rate/mobility maps and a grid solver ----------
    "atomic_quenching_parameters": "IMP.bff.quenching.maps",
    "diffusion_coefficient_map": "IMP.bff.quenching.maps",
    "radial_diffusion_map": "IMP.bff.quenching.maps",
    "slow_diffusion_near_atoms": "IMP.bff.quenching.maps",
    "quenching_rate_map": "IMP.bff.quenching.maps",
    "fret_rate_map": "IMP.bff.quenching.maps",
    "GridDiffusionSolver": "IMP.bff.quenching.solver",
    "GridDiffusionResult": "IMP.bff.quenching.solver",
    "diffusion_stability_limit": "IMP.bff.quenching.solver",
    "equilibrium_occupancy": "IMP.bff.quenching.solver",
    "DynamicAccessibleVolume": "IMP.bff.quenching.dynamic",
    # --- the particle model, end to end --------------------------------------
    "DyeDiffusionSimulation": "IMP.bff.quenching.model",
    "QuenchedDonorDecay": "IMP.bff.quenching.model",
    "ResidueSites": "IMP.bff.quenching.sites",
    "residue_sites": "IMP.bff.quenching.sites",
    "slow_factors_for_residues": "IMP.bff.quenching.sites",
    "quenching_rates_for_residues": "IMP.bff.quenching.sites",
    "quench_radii_for_residues": "IMP.bff.quenching.sites",
    # --- strip engine (fret) ----------------------------------------------------
    "parse_strip_mask": "IMP.bff.fret.strip",
    "default_strip_mask": "IMP.bff.fret.strip",
    "site_strip_mask": "IMP.bff.fret.strip",
    "strip_hierarchy": "IMP.bff.fret.strip",
    "strip_pdb_lines": "IMP.bff.fret.strip",
    "strip_obstacles": "IMP.bff.fret.strip",
    "select_atoms": "IMP.bff.fret.strip",
    # --- other readers ----------------------------------------------------------
    "read_dcd": "IMP.bff.cgdye.io.dcd",
    "read_nmr_restraints": "IMP.bff.cgdye.io.nmr_cif",
    "write_nmr_restraints": "IMP.bff.cgdye.io.nmr_cif",
}

#: a few exports are attributes with a different name in their module
_SOURCE_NAME = {
    "fps_schema_validate": "validate",
    "compare_av_and_rotamer_positions": "compare_positions",
    "compare_av_and_rotamer_pairs": "compare_pairs",
    "make_langevin_simulator": "make_simulator",
}


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
