"""
Conformer modeling
==================

This example illustrate modeling of T4 lysozme by rigid bodies (RBs) in PMI.
Note, this example serves illustration and not scientific purposes.
"""
import pathlib

import IMP
import IMP.core
import RMF
import IMP.atom
import IMP.bff
import IMP.bff.tools
import IMP.bff.restraints

import IMP.rmf
import IMP.pmi
import IMP.pmi.tools
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.stereochemistry


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
        cr.set_label(molname)
        output_objects.append(cr)
        crs.append(cr)
        mols.append(mol)

# %%
# Excluded volume - automatically more efficient due to rigid bodies
ev_weight = 1.0
evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=mols)
evr.add_to_model()
evr.set_weight(ev_weight)
output_objects.append(evr)

# %%
# Restrain in sphere
barrier_restraint = IMP.pmi.restraints.basic.ExternalBarrier(hier, 30)
barrier_restraint.add_to_model()
output_objects.append(barrier_restraint)

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
fps_json_fn = str(root_dir / "screening.fps.json")
score_set = "chi2_C1_33p"  # molecule in close c2, go from c2 -> open c1
fret_restraint = IMP.bff.restraints.AVNetworkRestraintWrapper(
    hier, fps_json_fn,
    mean_position_restraint=True,
    score_set=score_set
)
fret_restraint.add_to_model()
output_objects.append(fret_restraint)

# %%
# Adds the AV to the atom Hierarchy to display the dye mean position:
used_avs = fret_restraint.av_network_restraint.get_used_avs()
IMP.bff.tools.display_mean_av_positions(used_avs)


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
num_frames = 100
rex = IMP.pmi.macros.ReplicaExchange0(
    mdl,
    simulated_annealing=False,
    root_hier=hier,  # pass the root hierarchy
    monte_carlo_sample_objects=dof.get_movers(),  # pass MC movers
    output_objects=output_objects,
    monte_carlo_steps=10,
    global_output_directory='./out/',
    rmf_output_objects=None,
    monte_carlo_temperature=1.0,
    number_of_best_scoring_models=0,  # set >0 to store best PDB files (but this is slow to do online)
    number_of_frames=num_frames       # increase number of frames to get better results!
)
rex.execute_macro()
