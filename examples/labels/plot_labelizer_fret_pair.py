"""
Which pair to measure: a FRET assay for a conformational change
===============================================================

`plot_labelizer_score.py` asks where a dye can go. This asks the second
question, which is the one that decides an experiment: **given two
conformations of a protein, which pair of sites would actually report on the
change?**

The case is maltose-binding protein, which closes around its ligand by a hinge
bend — 1OMP open, 1ANF closed. A good pair has to satisfy three things at once,
and the score is the product of them:

1. both sites are labelable at all (the per-residue score);
2. the dye–dye distance sits where FRET responds, near :math:`R_0`;
3. that distance **changes** between the two states.

Along the way this shows the parts the first example does not: taking
:math:`R_0` from named dyes rather than asserting a number, the dye container,
and reading a result back out of one.
"""

import os
import tempfile

import matplotlib.pyplot as plt
import numpy as np

import IMP.bff as bff

data = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "..", "..", "test", "input", "labelizer")
apo_pdb = os.path.join(data, "1OMP.pdb")     # open, no maltose
holo_pdb = os.path.join(data, "1anf.pdb")    # closed, maltose bound

# --- 1. R0 from the dyes, not from a number ---------------------------------
#
# A Forster radius is derived from the donor's quantum yield and the overlap of
# its emission with the acceptor's absorption. Naming the dyes makes the number
# traceable; asserting `52.0` does not. Names resolve loosely, so the spelling
# an experimenter uses works: `Alexa488`, not `AlexaFluor488`.
donor, acceptor = "Alexa488", "Alexa647"
r0 = bff.forster_radius(bff.get_probe(donor), bff.get_probe(acceptor))
print("R0(%s -> %s) = %.1f A" % (donor, acceptor, r0))

# The same dyes, and the same R0, out of a container whose every column name is
# an MMFDB dictionary item. This is what `imp_bff_labelizer --dyes` reads.
# Into a scratch directory: an example that writes into the source tree
# leaves derived files behind that look like inputs the next time.
work = tempfile.mkdtemp()
container = os.path.join(work, "dyes.mmfdb.pto")
print("%d dyes written to %s"
      % (bff.probe_library_to_pto(container), os.path.basename(container)))
print("R0 via the container = %.1f A"
      % bff.probe_pto_forster_radius(container, donor, acceptor))

# --- 2. score both conformations --------------------------------------------
#
# No ConSurf grades ship for MalE, and under the published model a missing term
# makes the whole combined score unavailable -- correctly, but it leaves
# nothing to pair. These three terms need only the coordinates.
model = bff.LabelizerParameterList()
for tag, table in (("se", "N_SE11_MEAN_SURFACE_DIST"),
                   ("cr", "C_CR1_Name"),
                   ("ss", "C_SS1_SS")):
    model.append(bff.LabelizerParameter(tag, table, 1))

options = bff.LabelizerOptions()
apo_scores = bff.labelizer_score_structure(apo_pdb, model, options, "")
holo_scores = bff.labelizer_score_structure(holo_pdb, model, options, "")
apo = bff.labelizer_combined_by_key(apo_scores)
holo = bff.labelizer_combined_by_key(holo_scores)
print("\nscored %d positions in each conformation" % len(apo))

# A property of the whole molecule rather than of a site: MBP is acidic, and a
# charged dye is not indifferent to that.
net, positive, negative = bff.labelizer_global_charge(bff.labelizer_read_structure(holo_pdb))
print("net formal charge %+.0f (%+.0f / %.0f)" % (net, positive, negative))

# --- 3. rank the pairs -------------------------------------------------------
fret = bff.LabelizerFRETOptions()
fret.forster_radius = r0
fret.n_refine = 0          # the cone screen alone; see the caveat below

pairs = list(bff.labelizer_pair_scores_two_states(
    apo_pdb, holo_pdb, apo, holo, fret))
print("\n%d pairs; the ten worth measuring:" % len(pairs))
print("   %-14s %8s %8s %8s   %s" % ("pair", "d(apo)", "d(holo)", "delta", "score"))
for p in pairs[:10]:
    print("   %-14s %7.1f A %7.1f A %+7.1f A   %.4f"
          % ("%s%d-%s%d" % (p.asym_id_1, p.seq_id_1, p.asym_id_2, p.seq_id_2),
             p.distance, p.distance_2, p.distance_2 - p.distance, p.value))

# Every one of them contracts: the score has found the pairs that span the two
# domains, which is what a hinge closure looks like from the outside.
shifts = np.array([p.distance_2 - p.distance for p in pairs[:20]])
print("\ntop 20: %d of 20 contract, by %.1f to %.1f A"
      % (int((shifts < 0).sum()), -shifts.max(), -shifts.min()))

# A caveat worth carrying: this ranking used the cheap analytic dye position.
# It over-reaches by 2-5 A, and because the score peaks sharply at R = R0 that
# reorders the top of the list. Set `n_refine` to rebuild the best few with
# real accessible volumes before committing to an experiment.

# --- 4. the whole result as one file, and back again -------------------------
out = os.path.join(work, "MalE_apo_holo.mmfdb.pto")
settings = bff.labelizer_settings_json(model, options, fret, "")
bff.labelizer_write_pto(out, holo_pdb, holo_scores, pairs[:500], settings)
print("\nwrote %s (%.2f MB)"
      % (os.path.basename(out), os.path.getsize(out) / 1048576.0))

# Reading it back needs nothing but the file: the scores, the pairs and the
# settings that produced them all travel together.
back_scores = bff.labelizer_read_pto_scores(out)
back_pairs = bff.labelizer_read_pto_pairs(out)
print("read back: %d score rows, %d pair rows" % (len(back_scores),
                                                  len(back_pairs)))
print("recovered structure sha256 %s..."
      % bff.labelizer_extract_pto_structure(
          out, os.path.join(work, "MalE-recovered.pdb"))[:16])

# --- 5. the picture ----------------------------------------------------------
#
# The Cbeta difference map is the change itself, before any dye is involved:
# how much further apart every pair of positions ends up. The pair score is
# essentially this map, weighted by where FRET can see it.
first, flat = bff.labelizer_cbeta_difference_map(bff.labelizer_read_structure(apo_pdb),
                                          bff.labelizer_read_structure(holo_pdb))
matrix = np.asarray(flat)
n = int(round(np.sqrt(matrix.size)))
matrix = matrix.reshape(n, n)

fig, (ax, bx) = plt.subplots(1, 2, figsize=(11, 4.6))
limit = np.abs(matrix).max()
image = ax.imshow(matrix, cmap="RdBu", vmin=-limit, vmax=limit,
                  origin="lower", extent=[first, first + n, first, first + n])
ax.set_title("C$\\beta$ distance change, holo $-$ apo (\\AA)")
ax.set_xlabel("residue")
ax.set_ylabel("residue")
fig.colorbar(image, ax=ax, fraction=0.046)

for p in pairs[:15]:
    ax.plot(p.seq_id_2, p.seq_id_1, "k.", ms=4)
    ax.plot(p.seq_id_1, p.seq_id_2, "k.", ms=4)
ax.text(0.02, 0.97, "dots: the 15 best pairs", transform=ax.transAxes,
        va="top", fontsize=8)

bx.scatter([p.distance for p in pairs], [p.distance_2 for p in pairs],
           s=3, alpha=0.25, color="#3b6ea5", label="all pairs")
bx.scatter([p.distance for p in pairs[:20]], [p.distance_2 for p in pairs[:20]],
           s=22, color="#b5651d", label="top 20")
span = [min(p.distance for p in pairs), max(p.distance for p in pairs)]
bx.plot(span, span, "k--", lw=0.8, label="no change")
bx.axvline(r0, color="0.7", lw=0.8)
bx.axhline(r0, color="0.7", lw=0.8)
bx.set_xlabel("dye-dye distance, apo (\\AA)")
bx.set_ylabel("dye-dye distance, holo (\\AA)")
bx.set_title("the pairs that move, near $R_0$")
bx.legend(fontsize=8, loc="upper left")
fig.tight_layout()
plt.show()
