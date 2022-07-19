"""
Scoring of structures
=====================
This examples illustrates how to use  ``IMP.bff.AVNetworkRestraint`` for
scoring structures in a trajectory.
"""
import json
import pathlib
import tqdm
import pylab as plt

import RMF
import IMP
import IMP.rmf
import IMP.atom
import IMP.bff

m = IMP.Model()
root_path = pathlib.Path(IMP.bff.get_example_path("structure")) / "T4L"

# %%
# Load a RMF file and create a hierarchy
rmf_fn = str(root_path / "t4l_docking.rmf3")
frame_index = 0
f = RMF.open_rmf_file_read_only(rmf_fn)
hier = IMP.rmf.create_hierarchies(f, m)[frame_index]
IMP.rmf.load_frame(f, RMF.FrameID(frame_index))

# %%
# Load fps.json file and create restraint for scoring
fps_json_path = root_path / "screening.fps.json"
with open(fps_json_path) as fp:
    fps_json = json.load(fp)
score_set_c1 = "chi2_C1_33p"
fret_restraint = IMP.bff.AVNetworkRestraint(hier, str(fps_json_path), score_set=score_set_c1)
v = fret_restraint.unprotected_evaluate(None)

# %%
# Score each frame in the RMF file and plot score
scores = list()
for frame in tqdm.tqdm(f.get_root_frames()):
    IMP.rmf.load_frame(f, frame)
    v = fret_restraint.unprotected_evaluate(None)
    scores.append(v)

plt.plot(scores, "o-")
plt.show()

