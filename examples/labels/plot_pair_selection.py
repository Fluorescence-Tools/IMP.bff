"""
Which FRET pair to measure next
===============================
An ensemble of candidate structures and a list of labelling pairs you *could*
measure: which one do you measure first?

``IMP.bff.select_informative_pairs`` answers it the way Olga does (Dimura
*et al.*, *Nat. Commun.* **11**, 5394, 2020). A measurement is worth making if
it separates the ensemble -- if the conformers it cannot tell apart are also
the ones that are structurally alike. The selector adds the pair that leaves
the smallest *expected* RMSD between the true structure and the one the
measurements would single out, then repeats on what is left. What comes back is
an ordered list and a decay curve: how much precision each further measurement
buys, and where buying more stops paying.

Two inputs are needed, and both are computed here rather than assumed: the FRET
efficiency of every candidate pair in every frame, and the pairwise RMSD
between the frames.
"""
import numpy as np
import pylab as plt

import RMF
import IMP
import IMP.atom
import IMP.core
import IMP.rmf
import IMP.bff

# %%
# The ensemble. A hundred docked T4 lysozyme conformers, and the labelling
# positions and pairs of an fps.json file.
m = IMP.Model()
rmf_fn = IMP.bff.get_example_path("structure/T4L/t4l_docking.rmf3")
fps_json_path = IMP.bff.get_example_path("structure/T4L/fret.fps.json")

f = RMF.open_rmf_file_read_only(rmf_fn)
hier = IMP.rmf.create_hierarchies(f, m)[0]
IMP.rmf.load_frame(f, RMF.FrameID(0))

# The beads the RMSD is measured over, taken before the restraint decorates the
# hierarchy with accessible volumes -- those move with the dye, not the fold.
beads = [IMP.core.XYZ(p) for p in IMP.atom.get_leaves(hier)]
frames = list(f.get_root_frames())
print(f"{len(frames)} frames, {len(beads)} beads")

# %%
# The candidate pairs. A score set names the subset of the file's distances a
# given experiment covers; its pairs are the candidates to choose among.
restraint = IMP.bff.ProbeNetworkRestraint(
    hier, fps_json_path, score_set="chi2_C1_33p")
pair_names = restraint.get_pair_names()
print(f"{len(pair_names)} candidate pairs")

# %%
# One pass over the ensemble collects both inputs. ``get_pair_efficiencies``
# re-evaluates the accessible volumes and reports the mean FRET efficiency of
# every candidate, in the order ``get_pair_names`` gives -- so the columns of
# ``effs`` are the candidates and its rows are the frames.
effs = np.empty((len(frames), len(pair_names)))
coords = np.empty((len(frames), len(beads), 3))
for i, frame in enumerate(frames):
    IMP.rmf.load_frame(f, frame)
    effs[i] = restraint.get_pair_efficiencies()
    coords[i] = [b.get_coordinates() for b in beads]

# %%
# How far apart the ensemble places its conformers, after superposition. This
# is the yardstick the selection is scored against: a pair that separates
# structures which are already alike has told you nothing.
rmsds = IMP.bff.pairwise_rmsd(coords, True)
print(f"RMSD spread {rmsds[np.triu_indices_from(rmsds, 1)].mean():.2f} A mean, "
      f"{rmsds.max():.2f} A max")

# %%
# The selection. ``err`` is the expected absolute error of a FRET efficiency
# measurement -- the selector's only notion of what an experiment can resolve,
# and what stops it from believing a hair's-breadth separation.
selected, decay = IMP.bff.select_informative_pairs(
    effs, rmsds, err=0.06, max_pairs=10)

start = IMP.bff.expected_rmsd(
    rmsds.ravel(), np.zeros(rmsds.size), 1, 0.99, len(frames))
print(f"{'#':>2}  {'pair':<12}  {'expected RMSD':>13}  {'gain':>6}")
print(f"{'-':>2}  {'(no data)':<12}  {start:>11.2f} A  {'':>6}")
previous = start
for rank, (index, value) in enumerate(zip(selected, decay), start=1):
    print(f"{rank:>2}  {pair_names[index]:<12}  {value:>11.2f} A  "
          f"{previous - value:>5.2f} A")
    previous = value

# %%
# The decay curve. The first pairs are worth several angstrom each; the curve
# flattens once the remaining candidates report on distances the earlier ones
# have already pinned down.
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(range(len(decay) + 1), np.concatenate([[start], decay]), "o-")
ax.set_xticks(range(len(decay) + 1))
ax.set_xticklabels(["none"] + [pair_names[i] for i in selected],
                   rotation=45, ha="right")
ax.set_xlabel("measurements made, in the order chosen")
ax.set_ylabel("expected mean RMSD / $\\AA$")
ax.set_title("Precision bought by each further FRET pair")
fig.tight_layout()
plt.show()
