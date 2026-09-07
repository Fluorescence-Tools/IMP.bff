"""
Choosing a labelling site: the Labelizer score, natively
========================================================

Which residue of a protein should carry the dye? The Labelizer score
(Gebhardt *et al.*, *Nat. Commun.* **16**, 3305, 2025) answers it by combining
four per-residue quantities into one number: how conserved the position is, how
exposed it is, what secondary structure it sits in, and how much the residue
already resembles a cysteine.

This walks the whole model on mouse BID (PDB 1DDB, model 39, chain A), the
worked example the reference ships, and then asks the second question the
reference asks: which *pair* of those sites would make the most informative
FRET measurement.

Two things worth knowing before reading the numbers:

* **A score is a likelihood ratio, not a probability.** Each table holds
  ``P(labelable | observation) / P(labelable)``, so a value above 1 means the
  observation makes a site more likely to be labelable than the base rate, and
  the combined score is unbounded above.
* **The published arithmetic is the default, defects and all**, because that is
  what the paper's numbers were computed with. ``LL_MODEL_CORRECTED`` selects
  the arithmetic the reference's own documentation describes.

The companion example, ``plot_labelizer_fret_pair.py``, takes the second
question further: given *two* conformations, which pair reports on the change.
The model itself is described in the "Choosing a labelling site" page of the
API documentation.
"""

import os
import tempfile

import matplotlib.pyplot as plt
import numpy as np

import IMP.bff as bff

# --- 1. the structure, and the conservation grades it needs -----------------
#
# Conservation is imported, never computed here: an alignment and a rate
# estimate are a different program. The reference's own example file is used,
# which is also how its published output is reproducible -- see the caveat in
# section 5.
data = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "..", "..", "test", "input", "labelizer")
pdb = os.path.join(data, "1DDB-39.pdb")
conservation = os.path.join(data, "1DDB-39_cs.pdb")

structure = bff.ll_read_structure(pdb)
print("%d residues, %d atoms" % (len(structure.residues), len(structure.vdw)))

# --- 2. the structural features, each on its own ----------------------------
#
# These are the quantities the model scores. `ll_dssp` is a native
# Kabsch-Sander assignment -- the reference shells out to the DSSP binary and
# will not start on macOS at all -- and it reproduces that binary exactly on
# every residue of this structure.
secondary = bff.ll_dssp(structure)
depth = np.asarray(bff.ll_residue_depth(structure, 1.4, 590))
rsa = np.asarray(bff.ll_relative_solvent_accessibility(
    structure, bff.LL_MAXASA_WILKE, 1.4, 590))
hse = np.asarray(bff.ll_half_sphere_exposure(structure, 13.0)).reshape(-1, 2)

print("secondary structure: " + "".join(
    "%s=%d " % (s, secondary.count(s)) for s in sorted(set(secondary))))
print("residue depth  %.2f - %.2f A" % (depth.min(), depth.max()))

# --- 3. the score ------------------------------------------------------------
model = bff.ll_model_paper()
print("model: " + ", ".join("%s(w=%d)" % (p.tag, p.weight) for p in model))

scores = bff.ll_score_structure(pdb, model, bff.LlOptions(), conservation)

# The rows are tidy -- one per (position, score_type) -- which is the shape the
# container stores and the shape a dataframe wants.
combined = {bff.ll_residue_key(r.asym_id, r.seq_id): r.value
            for r in scores if r.score_type == "combined" and r.status == "scored"}
best = sorted(combined.items(), key=lambda kv: -kv[1])[:8]
print("\nbest sites")
for key, value in best:
    print("   %-6s %.4f" % (key, value))

# --- 4. which pair to measure ------------------------------------------------
#
# A pair is informative when both sites are labelable *and* the dye-dye
# distance sits where FRET responds to it, which is near R0. The screen places
# the dye with the reference's analytic "alpha cone"; `n_refine` then rebuilds
# the best few with a real accessible volume, which is IMP.bff's own -- the
# reference calls LabelLib, which this package does not use.
options = bff.LlFretOptions()
options.forster_radius = 52.0
options.n_refine = 5

pairs = list(bff.ll_pair_scores(pdb, combined, options))
print("\n%d pairs above the label-score threshold" % len(pairs))
print("best pairs")
for p in pairs[:8]:
    print("   %s%d-%s%d  %.4f  d=%.1f A"
          % (p.asym_id_1, p.seq_id_1, p.asym_id_2, p.seq_id_2,
             p.value, p.distance))

# The refinement is worth seeing on its own, because it is the reason the two
# tiers exist: the cheap cone systematically over-reaches, so the pairs it
# ranked highest are not the ones that survive a real dye cloud.
refined = [p for p in pairs if p.probe_model == bff.PROBE_MODEL_ACCESSIBLE_VOLUME]
print("\nthe %d pairs the cone ranked best, after rebuilding the dye clouds"
      % len(refined))
for p in refined:
    print("   %s%d-%s%d  %.4f  d=%.1f A"
          % (p.asym_id_1, p.seq_id_1, p.asym_id_2, p.seq_id_2,
             p.value, p.distance))

# --- 5. one file out ---------------------------------------------------------
#
# Not six CSVs and a zip: one container holding the structure verbatim, the
# scores with their columns named by MMFDB dictionary items, and the complete
# settings. A position that was not scored carries a status and no number.
# Into a scratch directory: an example that writes into the source tree
# leaves derived files behind that look like inputs the next time.
work = tempfile.mkdtemp()
out = os.path.join(work, "1DDB-39.mmfdb.pto")
settings = bff.ll_settings_json(model, bff.LlOptions(), options, conservation)
bff.ll_write_pto(out, pdb, scores, pairs[:200], settings)
print("\nwrote %s (%d bytes)" % (os.path.basename(out), os.path.getsize(out)))
print("recovered structure sha256 %s"
      % bff.ll_extract_pto_structure(
          out, os.path.join(work, "1DDB-39-recovered.pdb")))

# A caveat the reference's own example carries: it passes its conservation PDB
# as both input and output, so re-running it feeds the previous run's scores
# back in as grades. The lookup has a two-cycle, which is why this structure's
# conservation term takes only two distinct values. The score machinery is
# exact; the input is not. See okf/validation/labelizer_ab.md.

# --- 6. the picture ----------------------------------------------------------
seq = [r.seq_id for r in structure.residues]
values = [combined.get(bff.ll_residue_key(r.chain, r.seq_id), np.nan)
          for r in structure.residues]

fig, (ax, bx) = plt.subplots(2, 1, figsize=(9, 5.5), sharex=True,
                             gridspec_kw={"height_ratios": [2, 1]})
ax.plot(seq, values, lw=1.0, color="#3b6ea5")
ax.axhline(1.0, color="0.6", lw=0.8, ls="--")
ax.axhline(0.5, color="#b5651d", lw=0.8, ls=":")
for key, value in best[:5]:
    i = seq.index(int(key[1:]))
    ax.annotate(key, (seq[i], value), textcoords="offset points",
                xytext=(0, 5), ha="center", fontsize=8)
ax.set_ylabel("combined label score")
ax.set_title("Labelizer score along mouse BID (1DDB model 39, chain A)")
ax.text(0.99, 0.05, "dashed: base rate   dotted: pairing threshold",
        transform=ax.transAxes, ha="right", fontsize=8, color="0.4")

bx.plot(seq, depth, lw=1.0, color="#5a5a5a", label="residue depth (A)")
bx.plot(seq, rsa * 4.0, lw=1.0, color="#7aa457", label="RSA (x4)")
bx.set_xlabel("residue")
bx.set_ylabel("exposure")
bx.legend(fontsize=8, loc="upper right")
fig.tight_layout()
plt.show()
