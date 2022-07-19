"""
Accessible Volume Restraint
===========================
This example illustrates how to use the accessible volume restraint. Note,
there is another version of the AV restraint wrapped for PMI that can operate
of mean dye positions. The ``IMP.bff.AVNetworkRestraint`` recomputes the AVs
used for scoring in each evaluation. This can be costly.

The use of the PMI wrapped ``IMP.bff.AVNetworkRestraint`` is illustrated in the
structural modeling use cases.
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

rmf_fn = str(root_path / "t4l_docking.rmf3")
frame_index = 0
f = RMF.open_rmf_file_read_only(rmf_fn)
hier = IMP.rmf.create_hierarchies(f, m)[frame_index]
IMP.rmf.load_frame(f, RMF.FrameID(frame_index))

# %%
# Scoring
# -------
# The AV network restraint

fps_json_path = root_path / "screening.fps.json"
fps_json = ""
with open(fps_json_path) as fp:
    fps_json = json.load(fp)

score_set_c1 = "chi2_C1_33p"
fret_restraint = IMP.bff.AVNetworkRestraint(hier, str(fps_json_path), score_set=score_set_c1)
v = fret_restraint.unprotected_evaluate(None)
print("Score: {:.1f}".format(v))

