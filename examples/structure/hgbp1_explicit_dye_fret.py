## \example structure/hgbp1_explicit_dye_fret.py
# Explicit dyes on hGBP1: rotamer ensembles, FRET between two sites, an fps.json
# with R1 positions, and a short Langevin run -- all through the flat IMP.bff API.
#
# hGBP1 (PDB 1DG3, chain A) labelled at 481 (donor, Alexa488 C1R library) and
# 496 (acceptor, Alexa594 C1R library). Every name used here is reachable as
# IMP.bff.<Name>: the namespace is flat, and there are no sub-packages. What
# the numbers mean, and how they compare with the AV model, is recorded in
# okf/validation/av_vs_rotamer.md (PRD-108).

import json

import IMP
import IMP.atom
import IMP.bff
from IMP.bff import get_structure_dir

pdb = str(get_structure_dir("1DG3.pdb"))

# --- 1. rotamer ensembles: a screened library at each site (fps type R1) ---
donor = IMP.bff.RotamerEnsemble.from_site(pdb, "A", 481, "AlexaFluor 488 C1R cutoff30", position_name="A481")
acceptor = IMP.bff.RotamerEnsemble.from_site(pdb, "A", 496, "AlexaFluor 594 C1R cutoff30", position_name="A496")
print(f"donor: {donor.n_rotamers} rotamers, Z={donor.partition:.3f}, mean position {donor.mean_position.round(1)}")
print(f"acceptor: {acceptor.n_rotamers} rotamers, Z={acceptor.partition:.3f}")

# --- 2. FRET between them: R0 from the spectra at the pair's <kappa2>, all regimes ---
eff = donor.pair_distribution_from_probes(acceptor, "AlexaFluor 488", "AlexaFluor 594")
print(f"R0 = {eff.forster_radius / 10:.2f} nm at <kappa2> = {eff.kappa2_avg:.3f}; "
      f"E_static {eff.static_efficiency:.3f}, E_dynamic1 {eff.dynamic1:.3f}, E_dynamic2 {eff.dynamic2:.3f}")
# the same pair as an AV-style statistic (Rmp, <R_DA>, <R_DA>_E, sigma) -- an ensemble is an AccessibleVolume
print("av_pair_statistics:", [round(x, 1) for x in IMP.bff.av_pair_statistics(donor, acceptor, forster_radius=52.0)])

# --- 3. an fps.json carrying the two R1 positions and the predicted <R_DA>_E ---
positions = IMP.bff.rotamer_positions_payload({"A481": donor, "A496": acceptor})
distances = IMP.bff.distances_from_ensembles({"A481": donor, "A496": acceptor}, [("A481", "A496")], 52.0)
IMP.bff.write_rotamer_fps("hgbp1_rotamer.fps.json", positions, distances)
print("wrote hgbp1_rotamer.fps.json:",
      {k: round(v["distance"], 1) for k, v in json.loads(distances).items()})
# (ProbeNetworkRestraint scores AV positions only -- IMP.bff.fps_positions_for_docking() strips R1 ones)

# --- 4. an explicit dye in motion: Langevin dynamics of Alexa488 at 481 ---
model = IMP.Model()
protein = IMP.atom.read_pdb(pdb, model, IMP.atom.NonWaterPDBSelector())
dye = IMP.atom.read_mol2(str(get_structure_dir("alexa488_r48.mol2")), model)
IMP.bff.attach_probes(protein, [IMP.bff.ProbeAttachment(dye, "A", 481)], strip_site_sidechain=True)
sampler = IMP.bff.AttachedProbeDynamics(protein, dye, str(get_structure_dir("alexa488_r48.mol2")), "A", 481,
                                     integrator="md", temperature=300.0, seed=1)
sampler.minimize(200)
traj = sampler.run(2000, write_every=100)
t_kin = sum(sampler.kinetic_temperature(k) for k in traj.kinetic_energy) / traj.n_frames
print(f"Langevin md: {traj.n_frames} frames, <T_kin> = {t_kin:.0f} K, <E_pot> = {traj.potential_energy.mean():.1f} kcal/mol")
