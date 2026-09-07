"""
FRET-restrained docking of HIV-1 reverse transcriptase and its DNA
=================================================================
Two rigid bodies, twenty measured distances, and the four things you do with
them: convert the data, score the structure you have, dock, and read the
result honestly.

The data is FPS's own docking test case -- HIV-1 reverse transcriptase (p66/p51,
PDB 1R0A) with its DNA primer/template -- shipped here as the original C# FPS
``LabelingPositions.txt`` and ``Distances.txt`` beside the converted
``hiv_rt.fps.json``. Eleven labelling positions: eight AV1 sites on the
protein, three AV3 sites on the DNA.

Everything below is also a command line; see ``imp_bff_fps --help``.
"""
# sphinx_gallery_thumbnail_path = 'img/example_structure_fret_docking.png'
import json
import math
import os
import tempfile

import IMP
import IMP.bff

IMP.set_log_level(IMP.SILENT)

HERE = IMP.bff.get_example_path("structure/HIV_RT")
PROTEIN = os.path.join(HERE, "protein_1R0A.pdb")
DNA = os.path.join(HERE, "dna.pdb")
FPS_JSON = os.path.join(HERE, "hiv_rt.fps.json")

# %%
# The labelling data
# ------------------
# A legacy FPS position names its site by **atom serial number**, which means
# nothing without the structure it was numbered in. Reading the file resolves
# each serial to a chain, residue and atom name against the PDBs given -- so
# the structures are an argument, not an afterthought.
doc = IMP.bff.read_old_lps_txt(
    os.path.join(HERE, "LabelingPositions.txt"), [PROTEIN, DNA])
positions = json.loads(doc.positions)
distances = json.loads(IMP.bff.read_old_distances_txt(
    os.path.join(HERE, "Distances.txt")))

print("molecules: %s" % ", ".join(doc.molecules))
print("positions: %d, distances: %d" % (len(positions), len(distances)))
for name in ("p66_Q6C", "p_1bp"):
    p = positions[name]
    print("  %-10s chain %s residue %s atom %-4s %s"
          % (name, p["chain_identifier"], p["residue_seq_number"],
             p["atom_name"], p["simulation_type"]))

# %%
# The shipped ``hiv_rt.fps.json`` is this, plus one thing the legacy format
# cannot say: which rigid body each position belongs to. The protein is body 0
# and the DNA body 1, matching the order the PDBs are passed in.

# %%
# Scoring the structure you have
# ------------------------------
# The first thing to run. It says whether the labelling file resolves against
# the structures at all, how many volumes were built, and which measured
# distance the model misses.
#
# Note the score set. ``all`` is all twenty distances and scores ``inf``,
# because the two involving ``p66_K287C`` have no model value at all: that
# site is buried at the protein-DNA interface and its accessible volume comes
# out **empty**. That is a real property of this complex, not a failure to
# hide -- so it is a named score set rather than a silently shorter table.
result = IMP.bff.score_structures([PROTEIN, DNA], FPS_JSON,
                                  score_set="resolved")
print("\nscore set 'resolved': chi2 = %.3f over %d distances, %d volumes"
      % (result.score, result.n_distances, result.n_avs))

with_all = IMP.bff.score_structures([PROTEIN, DNA], FPS_JSON, score_set="all")
missing = [p.name for p in with_all.pairs
           if not math.isfinite(p.distance_model)]
print("score set 'all':       chi2 = %s   (no model value: %s)"
      % (with_all.score, ", ".join(missing)))

# %%
# The worst-fitting pairs are where the model and the measurement disagree
# most, and are what to look at before believing a structure.
worst = sorted((p for p in result.pairs if math.isfinite(p.distance_model)),
               key=lambda p: -p.get_chi2())
print("\n%-24s %8s %8s %8s" % ("pair", "exp", "model", "chi2"))
for p in worst[:5]:
    print("%-24s %8.2f %8.2f %8.2f"
          % (p.name, p.distance_exp, p.distance_model, p.get_chi2()))

# %%
# Docking
# -------
# Each PDB becomes one rigid body. The mobile bodies are shuffled, then driven
# by conjugate gradients under the measured distances and an excluded-volume
# term. ``fixed_body = 0`` holds the protein still, so the DNA is what moves
# and the result is reported in the protein's frame.
#
# The optimiser is IMP's, not FPS's damped rigid-body dynamics, so **scores are
# comparable with FPS and coordinates are not**.
out_dir = tempfile.mkdtemp(prefix="fret_docking_")
IMP.random_number_generator.seed(1)

params = IMP.bff.DockingParameters()
params.n_frames = 200               # minimiser iterations
params.score_set = "resolved"
params.fixed_body = 0               # the protein is the frame
params.shuffle_max_translation = 8.0

docked = IMP.bff.dock_minimize([PROTEIN, DNA], FPS_JSON, out_dir, params)
print("\ndocked: chi2 = %.3f  (started from a random displacement)" % docked.score)
print("wrote %s" % ", ".join(sorted(os.listdir(out_dir))))

# %%
# One run is one local minimum. Repeating from different random starts is the
# only way to see whether a pose is the answer or an accident -- the spread of
# the scores is the result, not the best of them.
scores = []
for seed in (1, 2, 3):
    IMP.random_number_generator.seed(seed)
    trial = tempfile.mkdtemp(prefix="fret_docking_%d_" % seed)
    scores.append(IMP.bff.dock_minimize([PROTEIN, DNA], FPS_JSON, trial,
                                        params).score)
print("three independent starts: %s"
      % ", ".join("%.2f" % s for s in scores))

# %%
# Refinement
# ----------
# The local settle, with no shuffle: start from a pose that is already roughly
# right -- the deposited complex here -- and let it relax against the data.
refined = IMP.bff.refine_docking([PROTEIN, DNA], FPS_JSON,
                                 tempfile.mkdtemp(prefix="fret_refine_"),
                                 "resolved", 200)
print("\ndeposited complex: chi2 = %.3f -> refined %.3f"
      % (result.score, refined.score))

# %%
# Screening
# ---------
# The complementary question: rather than moving anything, score every
# candidate structure as it stands and rank the library, best first. Each path
# is one **candidate**, not one body -- screening asks "which of these
# structures fits the data", so a candidate has to carry every labelling site
# on its own.
#
# That is worth showing rather than asserting. Screening the protein alone
# against a file whose sites include the DNA gives a NaN, because four of the
# positions are on a chain this candidate does not have:
ranked = IMP.bff.screen_structures([PROTEIN, DNA], FPS_JSON, "resolved")
for entry in ranked:
    print("%10s  %s"
          % ("%.3f" % entry.score if math.isfinite(entry.score) else "nan",
             os.path.basename(entry.path)))

# %%
# A structure that cannot be scored keeps its row rather than being dropped --
# a table with a row missing looks like a smaller experiment. So **read the
# NaNs**: a table that is entirely NaN, as here, is not a library of bad models
# but a labelling file that does not match the library at all.
