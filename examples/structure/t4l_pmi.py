## \example bff/structure/T4L.py
"""
T4 lysozyme: Using PMI
======================
T4 lysozme

This example illustrate modeling of  by rigid bodies (RBs) in PMI.
Note, this example serves illustration and not scientific purposes.
"""
# sphinx_gallery_thumbnail_path = 'img/example_structure_t4l.png'

import sys
import pathlib

import IMP
import IMP.core
import RMF
import IMP.atom
import IMP.bff
import IMP.bff

import IMP.rmf
import IMP.pmi
import IMP.pmi.tools
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.stereochemistry

IMP.setup_from_argv(sys.argv, "Simulation of an atomic T4L system")
if IMP.get_is_quick_test():
    print("This example is too slow to test in debug mode - run without")
    print("internal tests enabled, or without the --run-quick-test flag")
    sys.exit(0)


output_objects = list()
root_dir = pathlib.Path(IMP.bff.get_example_path('structure')) / "T4L/"
output_dir = root_dir / 'modeling/'

# %%
# System setup
# -------------
# The TopologyReader reads the text file, and the BuildSystem macro constructs it.
# The RB decomposition of T4L is defined in the topology tile.
mdl = IMP.Model()
topol_fn = str(root_dir / 'topology.dat')
reader = IMP.pmi.topology.TopologyReader(
    topol_fn,
    pdb_dir=root_dir,
    fasta_dir=root_dir
)
bs = IMP.pmi.macros.BuildSystem(mdl, resolutions=[1])

# %%
# note you can call this multiple times to create a multi-state system
bs.add_state(reader)
hier, dof = bs.execute_macro()
state = hier.get_children()[0]
mol = state.get_children()[0]
rb = mol.get_children()[2]

# %%
# General restraints
# --------------------
# Connectivity keeps things connected along the backbone (ignores if inside same rigid body)
crs = []
moldict = bs.get_molecules()[0]
mols = []
for molname in moldict:
    for mol in moldict[molname]:
        IMP.pmi.tools.display_bonds(mol)
        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol)
        cr.add_to_model()
        output_objects.append(cr)
        crs.append(cr)
        mols.append(mol)

# %%
# Restrain in sphere
barrier_restraint = IMP.pmi.restraints.basic.ExternalBarrier(hier, 30)
barrier_restraint.add_to_model()
output_objects.append(barrier_restraint)

# %%
# The library builds the restraints; PMI's bookkeeping is the caller's. A
# `RestraintBase` subclass is what `output_objects` and the replica-exchange
# macro report through, and `IMP.pmi.restraints.RestraintBase` is a Python
# class in a module `IMP.bff` does not depend on -- so the adapter is here,
# in the script that already imports PMI, and the restraints it wraps come
# from `IMP.bff.probe_network_restraint_set` (C++).
class ProbeNetworkRestraintWrapper(IMP.pmi.restraints.RestraintBase):
    """The fps.json AV network restraint, as a PMI restraint."""

    def __init__(self, hier, fps_json_fn, score_set="", weight=1.0,
                 mean_position_restraint=False, sigma_DA=6.0,
                 label="ProbeNetworkRestraint", occupy_volume=True):
        super().__init__(hier.get_model(), label=label, weight=weight)
        self.rs = IMP.bff.probe_network_restraint_set(
            hier, str(fps_json_fn), name=self.name, score_set=score_set,
            mean_position_restraint=mean_position_restraint,
            sigma_DA=sigma_DA, occupy_volume=occupy_volume, weight=weight)
        self.restraint_sets = [self.rs]


# %%
# FRET restraint
# --------------
# The experimental information is contained in an ``fps.json`` file.
# Fps.json files can contain multiple experimental distance dataset.
# Different combinations of distances can be used for scoring. The set
# of distances that is used for scoring needs to be specified. There are
# two modes of operation of the PMI AV network restraint. The restraint
# can operate on the mean position of AVs. If the restraint operates on
# mean positions, the AV is treated as rigid body member and the distances
# between the mean dye positions computed for the initial RB arrangement
# are used for scoring (with the use of a transfer function). Elsewise,
# the AVs are recalculated on each restraint evaluation.
fps_json_fn = str(root_dir / "fret.fps.json")
score_set = "chi2_C3"  # molecule in close c2, go from c2 -> open c1
fret_restraint = ProbeNetworkRestraintWrapper(
    hier, fps_json_fn,
    mean_position_restraint=True,
    occupy_volume=True,
    score_set=score_set
)
fret_restraint.add_to_model()
output_objects.append(fret_restraint)

# %%
# Excluded volume - automatically more efficient due to rigid bodies
# Note: add after AV restraint if the dyes should occupy some volume
ev_weight = 1.0
evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=mols)
evr.add_to_model()
evr.set_weight(ev_weight)
output_objects.append(evr)


# %%
# Sampling
# --------
# First shuffle the system
for state in bs.system.get_states():
    IMP.pmi.tools.shuffle_configuration(state.get_hierarchy(), max_translation=10)
# Quickly move all flexible beads into place
dof.optimize_flexible_beads(10)

# %%
# Monte carlo sampling. For better results increase the number of frames
num_frames = 1000
rex = IMP.pmi.macros.ReplicaExchange(
    mdl,
    simulated_annealing=False,
    root_hier=hier,  # pass the root hierarchy
    monte_carlo_sample_objects=dof.get_movers(),  # pass MC movers
    output_objects=output_objects,
    monte_carlo_steps=100,
    global_output_directory='./out/',
    rmf_output_objects=None,
    monte_carlo_temperature=1.0,
    number_of_best_scoring_models=0,  # set >0 to store best PDB files (but this is slow to do online)
    number_of_frames=num_frames       # increase number of frames to get better results!
)
rex.execute_macro()
