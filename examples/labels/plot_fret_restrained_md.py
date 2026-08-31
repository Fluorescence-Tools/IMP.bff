"""
FRET-restrained molecular dynamics: driving T4L from C2 to C1
=============================================================
An accessible-volume network can *score* a structure. To make it *move* one, the
restraint has to have a gradient and a shape an integrator can live with, which
is what :class:`IMP.bff.AVFlatBottomRestraint` is: zero inside the experimental
error bars, harmonic outside them, and **linear** past that so the force is
capped and a badly-placed start cannot blow up the first step. It is the same
well an AMBER ``&rst`` record describes, so a restraint written here and one
written for an MD engine are the same restraint.

``md_flat_bottom_restraints`` builds the whole thing from an fps.json: it
computes both volumes, converts each measured distance and its asymmetric errors
into *mean-position* bounds with ``rmp_flat_bottom_bounds`` -- a conversion that
is not a constant offset, because it depends on the shapes of both volumes --
creates one probe particle per labelling position tethered to its site, and puts
a well between every measured pair.

The T4L file carries two independent sets of 33 measured distances, ``C1`` and
``C2``, for the same 33 position pairs: two conformational states of the same
protein. Here the structure is equilibrated under the C2 distances and then the
restraints are switched to C1, so the transition is driven by the data alone.
"""
import numpy as np
import pylab as plt

import IMP
import IMP.algebra
import IMP.atom
import IMP.core
import IMP.bff

IMP.set_log_level(IMP.SILENT)

# %%
# The structure. All heavy atoms, because the accessible volumes are computed
# against them; only the alpha carbons will be moved.
m = IMP.Model()
hier = IMP.atom.read_pdb(IMP.bff.get_example_path("structure/T4L/3GUN.pdb"), m,
                         IMP.atom.NonWaterNonHydrogenPDBSelector())
fps_json = IMP.bff.get_example_path("structure/T4L/fret.fps.json")

cas = IMP.atom.Selection(hier, atom_type=IMP.atom.AT_CA).get_selected_particles()
for p in cas:
    if not IMP.core.XYZR.get_is_setup(p):
        IMP.core.XYZR.setup_particle(p, 2.5)
    if not IMP.atom.Mass.get_is_setup(p):
        IMP.atom.Mass.setup_particle(p, 110.0)
    IMP.core.XYZ(p).set_coordinates_are_optimized(True)
print(f"{len(cas)} alpha carbons")

# %%
# A C-alpha elastic network: stiff along the backbone, soft for every other
# contact within 10 A. It holds the fold together without deciding what the
# fold is, so the FRET restraints are what pick the conformation.
def elastic_network(particles, cutoff=10.0, k_backbone=20.0, k_contact=0.3):
    xyz = [IMP.core.XYZ(p).get_coordinates() for p in particles]
    out = []
    for i in range(len(particles)):
        for j in range(i + 1, len(particles)):
            d = IMP.algebra.get_distance(xyz[i], xyz[j])
            if d < cutoff:
                k = k_backbone if j == i + 1 else k_contact
                out.append(IMP.core.DistanceRestraint(
                    m, IMP.core.Harmonic(d, k), particles[i], particles[j]))
    return out

network = elastic_network(cas)
print(f"{len(network)} elastic-network restraints")

# %%
# The FRET restraints. One call per score set: the volumes are built on the
# structure as it stands, so calling it again later re-derives the bounds from
# the geometry the trajectory has reached.
#
# ``tether_atom="CA"`` puts each probe on the alpha carbon of its labelling
# residue rather than on the C-beta the fps.json names, because only the alpha
# carbons move here -- a probe tethered to an atom nothing optimises cannot
# follow the structure.
def fret_restraints(score_set):
    # f_max = 15 kcal/mol/A at one error bar, tether 30: strong enough to bend
    # the network above, which is the balance any restrained refinement has to
    # strike. Too soft and the data is decoration; too stiff and the force
    # field is.
    return IMP.bff.md_flat_bottom_restraints(
        hier, fps_json, score_set, 15.0, 30.0, 100.0, "CA")

state_c2 = fret_restraints("chi2_C2_33p")
# The same 33 pairs under the C1 distances, on the *starting* structure, purely
# as the baseline the run is measured against.
baseline_c1 = fret_restraints("chi2_C1_33p")
names = list(state_c2.get_pair_names())
print(f"{len(names)} restrained pairs, "
      f"{len(state_c2.get_probes())} probes")

# %%
# Run. The reporter records every restrained probe-probe distance so the
# trajectory can be read as the experiment reads it.
def simulate(system, n_steps, temperature=300.0, report_every=250,
             rebuild_every=2000):
    """Run MD under the elastic network plus `system`, reporting distances.

    ``AVRebuildOptimizerState`` resamples the volumes every ``rebuild_every``
    steps and re-derives every well from the geometry the trajectory has
    reached. Without it the wells keep the bounds derived from the *starting*
    structure for the whole run, and on this system those bounds move by
    several angstrom as the fold opens.
    """
    sf = IMP.core.RestraintsScoringFunction(
        network + [system.get_restraints()])
    md = IMP.atom.MolecularDynamics(m)
    md.set_scoring_function(sf)
    md.set_maximum_time_step(2.0)
    md.add_optimizer_state(
        IMP.atom.VelocityScalingOptimizerState(m, cas, temperature))
    md.add_optimizer_state(
        IMP.bff.AVRebuildOptimizerState(m, system, rebuild_every))
    wells = system.get_wells()
    trace = []
    for _ in range(0, n_steps, report_every):
        md.optimize(report_every)
        trace.append([w.get_distance() for w in wells])
    # MD explores; a short minimisation settles into the well it found.
    cg = IMP.core.ConjugateGradients(m)
    cg.set_scoring_function(sf)
    cg.optimize(200)
    trace.append([w.get_distance() for w in wells])
    return np.array(trace)


N = 12000
start_c2 = [w.get_distance() for w in state_c2.get_wells()]
trace_c2 = simulate(state_c2, N)
print(f"under C2: {trace_c2.shape[0]} reports of {trace_c2.shape[1]} distances")

# %%
# Switch to C1. The volumes are rebuilt on the structure the C2 run reached, so
# the new bounds are derived from the geometry actually in hand.
state_c1 = fret_restraints("chi2_C1_33p")
# The network is re-anchored on the C2 conformer too. It stands for *local*
# structure -- what a fold will not give up -- and not for one conformation;
# left on the crystal coordinates it becomes a memory of where the run started,
# and the C1 phase stalls at 29 of 33 because the network is pulling back.
network = elastic_network(cas)
trace_c1 = simulate(state_c1, N)
print(f"under C1: {trace_c1.shape[0]} reports")

# %%
# Did the data move the structure? A distance is *satisfied* when it sits inside
# its flat bottom. Count how many of the 33 are satisfied under each set, before
# and after the switch.
def satisfied(wells, distances):
    return sum(1 for w, d in zip(wells, distances)
               if w.get_bounds()[1] <= d <= w.get_bounds()[2])

wells_c2, wells_c1 = state_c2.get_wells(), state_c1.get_wells()
base = baseline_c1.get_wells()
print(f"C1 distances satisfied by the crystal structure: "
      f"{satisfied(base, [w.get_distance() for w in base])} of {len(base)}")
print(f"C2 distances satisfied, start -> end of the C2 run: "
      f"{satisfied(wells_c2, start_c2)} -> {satisfied(wells_c2, trace_c2[-1])}"
      f" of {len(wells_c2)}")
# Note the counts do not reach 33 of 33. With the volumes held fixed they do --
# but only because the wells are then a target derived from a structure the run
# has already left. Rebuilding moves the target onto the geometry in hand, which
# is a harder and a truer one.
print(f"C1 distances satisfied, start -> end of the C1 run: "
      f"{satisfied(wells_c1, trace_c1[0])} -> {satisfied(wells_c1, trace_c1[-1])}"
      f" of {len(wells_c1)}")

# %%
# The trajectory, read as the experiment reads it. Each line is one restrained
# probe-probe distance; the dashed line marks the switch from the C2 targets to
# the C1 targets.
trace = np.vstack([trace_c2, trace_c1])
fig, ax = plt.subplots(figsize=(8, 4.5))
for k in range(trace.shape[1]):
    ax.plot(trace[:, k], lw=0.8, alpha=0.7)
ax.axvline(len(trace_c2) - 0.5, color="k", ls="--", lw=1.2)
ax.text(len(trace_c2) - 0.5, ax.get_ylim()[1], " restraints switched C2 -> C1",
        va="top", fontsize=9)
ax.set_xlabel("report (250 MD steps each)")
ax.set_ylabel(r"probe-probe distance / $\AA$")
ax.set_title("T4 lysozyme driven from the C2 distances to the C1 distances")
fig.tight_layout()
plt.show()


# %%
# Out to OpenMM
# ------------
# IMP has no bridge to an MD engine -- ``IMP.modeller`` is its only external
# package interface and Modeller is not one -- so the way to run these
# restraints in OpenMM is to export them. ``write_openmm_restraints`` writes a
# document with everything an OpenMM script has to add to a ``System`` built
# from the same structure, **in OpenMM's units**: nanometres and kJ/mol.
import json
import tempfile
import os

work = tempfile.mkdtemp()
out = os.path.join(work, "t4l_c1.openmm.json")
IMP.bff.write_openmm_restraints(state_c1, out, hier)
doc = json.load(open(out))

# ...and the OpenMM around it, as a script that runs. The restraint table is
# embedded, so the generated file is the whole input.
script = os.path.join(work, "run_t4l_c1.py")
IMP.bff.write_openmm_script(
    state_c1, script, hier,
    IMP.bff.get_example_path("structure/T4L/3GUN.pdb"))
print(f"generated {os.path.basename(script)}, "
      f"{len(open(script).read().splitlines())} lines")
print(f"{len(doc['probes'])} probes, {len(doc['restraints'])} restraints, "
      f"units {doc['units']['length']}/{doc['units']['energy']}")
print(json.dumps(doc["restraints"][0], indent=2))

# %%
# What consumes it. Attachment atoms are named by chain, residue and atom
# rather than indexed, because an index depends on how the reader built its
# topology -- whether it kept hydrogens, waters, altlocs -- and a spec that
# names atoms by position in someone else's file silently restrains the wrong
# ones.
#
# .. code-block:: python
#
#     import openmm, openmm.unit as unit
#
#     doc = json.load(open("t4l_c1.openmm.json"))
#     index_of = {(c.id, r.id, a.name): a.index
#                 for c in topology.chains() for r in c.residues()
#                 for a in r.atoms()}
#
#     probe_index = {}
#     tether = openmm.HarmonicBondForce()
#     for probe in doc["probes"]:
#         probe_index[probe["name"]] = system.addParticle(probe["mass"])
#         at = probe["attachment"]
#         tether.addBond(index_of[(at["chain"], at["residue"], at["atom"])],
#                        probe_index[probe["name"]],
#                        probe["tether"]["length"], probe["tether"]["k"])
#     system.addForce(tether)
#
#     wells = openmm.CustomBondForce(doc["energy_expression"])
#     for name in doc["per_bond_parameters"]:
#         wells.addPerBondParameter(name)
#     for r in doc["restraints"]:
#         wells.addBond(probe_index[r["probe_1"]], probe_index[r["probe_2"]],
#                       [r[p] for p in doc["per_bond_parameters"]])
#     system.addForce(wells)
#
# The probe positions in the document are the starting coordinates those new
# particles need, appended to whatever the topology's own positions are.

# %%
# The exported well *is* the well. ``openmm_flat_bottom_energy`` is one string
# used by both sides, and it can be checked here without OpenMM installed --
# OpenMM's ``select`` and ``step`` are two lines of Python, and ``^`` is
# Python's ``**``.
def openmm_eval(expression, **values):
    ns = dict(values)
    ns["select"] = lambda c, a, b: a if c else b
    ns["step"] = lambda x: 1.0 if x >= 0 else 0.0
    return eval(expression.replace("^", "**"), {"__builtins__": {}}, ns)


A_TO_NM, KCAL_TO_KJ = 0.1, 4.184
example = state_c1.get_wells()[0]
r1, r2, r3, r4 = example.get_bounds()
k2, k3 = example.get_force_constants()

grid = np.linspace(0.6 * r1, 1.4 * r4, 300)
here = np.array([openmm_eval(
    "select(step(r1-r), k2*(r1-r2)^2 + 2*k2*(r1-r2)*(r-r1),"
    "select(step(r2-r), k2*(r-r2)^2,"
    "select(step(r-r4), k3*(r4-r3)^2 + 2*k3*(r4-r3)*(r-r4),"
    "select(step(r-r3), k3*(r-r3)^2, 0))))",
    r=d, r1=r1, r2=r2, r3=r3, r4=r4, k2=k2, k3=k3) for d in grid])
there = np.array([openmm_eval(
    doc["energy_expression"], r=d * A_TO_NM, r1=r1 * A_TO_NM, r2=r2 * A_TO_NM,
    r3=r3 * A_TO_NM, r4=r4 * A_TO_NM, k2=k2 * KCAL_TO_KJ * 100.0,
    k3=k3 * KCAL_TO_KJ * 100.0) for d in grid])
print("worst |exported - local| over the well: "
      f"{np.max(np.abs(there - here * KCAL_TO_KJ)):.2e} kJ/mol")

fig, ax = plt.subplots(figsize=(6.5, 3.6))
ax.plot(grid, there, lw=2, label="as exported to OpenMM")
ax.plot(grid, here * KCAL_TO_KJ, lw=1, ls="--", color="k",
        label="as evaluated by IMP.bff")
ax.axvspan(r2, r3, color="0.9", label="flat bottom (within the error bars)")
ax.set_xlabel(r"probe-probe distance / $\AA$")
ax.set_ylabel("energy / kJ mol$^{-1}$")
ax.set_title(f"{example.get_name()}: harmonic walls, capped tails")
ax.legend(fontsize=8)
fig.tight_layout()
plt.show()
