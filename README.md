\brief Bayesian Fluorescence Framework

# imp_bff {#imp_bff}

Bayesian Fluorescence Framework (BFF) processes and analyzes 
fluorescence data. BFF contains functions to computed label 
distributions and forward models. Among other data the program 
imp_bff samples can sample biomolecular conformations restrained 
by experimental data by  simulating fluorophore distributions 
around attachment sites and  by comparing simulated observables 
to experimental data.

Labels are modelled two ways: as **accessible volumes** (`IMP.bff.AV`,
`ProbeNetworkRestraint`, the `fret` docking layer) or as **explicit dyes** — an
atomistic dye + linker placed on a residue from a rotamer library (the
FRETpredict libraries ship as module data) or sampled over its linker degrees
of freedom, screened against the protein and turned into R0/κ²/FRET
efficiencies (`IMP.bff.FRETRotamer`, `IMP.bff.attach_dyes`, ...). Every public
name is reachable flat as `IMP.bff.<Name>`; the manual page
`doc/manual/structure/structure_cgprobe.ipynb` walks the explicit route, and
`okf/cgprobe.md` records how the code is organised.

## Install

Two packages, one import name. Both give you `import IMP.bff`; they own the
same files, so install one or the other, not both.

**`bff` — the core, no IMP.** Label distributions and accessible volumes,
explicit dyes from a rotamer library (`FRETRotamer`), dye diffusion on a
grid (`ProbeDiffusionSimulation`, `DynamicAccessibleVolume`), linker
sampling (`RRTTree`), probe force fields (`create_forcefield_system`),
side-chain packing (`pack_protein_sidechains`), the whole fitting stack (`Minimizer`,
`TCSPCDecay`, `ChiSquared`, FCS, κ²), the labelling-site score, and
structure and trajectory IO.

```bash
pip install bff             # wheels for Linux and macOS
conda install -c conda-forge bff
```

**`imp.bff` — the IMP module.** The same library built inside the Integrative
Modeling Platform, with the connection layer: restraints IMP's optimizers
score, decorators on IMP particles, `IMP.atom` hierarchies — and the two
dye roads that run on them, attaching a dye onto a structure
(`attach_probes`, `place_probe`) and driving it with Langevin dynamics
(`make_langevin_simulator`, `AttachedProbeDynamics`), which use IMP's own
integrators.

```bash
conda install -c conda-forge imp.bff
```

`IMP.bff.get_build()` says which one you have: `"core"` or `"imp"`.

The large data (rotamer libraries, coarse-grained probe inputs — 62 MB) is
downloaded on first use into a per-user cache. `imp_bff_fetch_data` takes it
all at once, which is what you want before going offline or in a container
image.

## Coming from LabelLib

`IMP.bff.labellib` is [LabelLib](https://github.com/Fluorescence-Tools/LabelLib)'s
interface, name for name — the same functions, the same argument order, the
same array conventions — so a script moves over by changing its import:

```python
# import LabelLib as ll
from IMP.bff import labellib as ll

av = ll.dyeDensityAV1(atoms_xyzr, source_xyz, 20.0, 2.0, 3.5, 0.9)
d = ll.meanDistance(av1, av2)
```

`dyeDensityAV1`, `dyeDensityAV3`, `minLinkerLength`, `addWeights`,
`meanDistance`, `meanEfficiency`, `sampleDistanceDistInv` and `Grid3D` are all
there, with their `_arr` twins. The numbers are bff's own lattice search, so
they agree with LabelLib to within the discretisation rather than bit for bit;
`test/test_labellib_interface.py` measures that against LabelLib itself
wherever it is installed. bff's own spelling of the same computation is
`IMP.bff.get_av` (arrays) and `IMP.bff.get_av_from_pdb` (a structure file),
which return an `AccessibleVolume` carrying the cloud, the density and the
grid.

## A first volume

```python
import IMP.bff

av = IMP.bff.get_av_from_pdb("structure.pdb", chain="A", resseq=132,
                             atom_name="CB", linker_length=20.0,
                             linker_width=2.0, r1=3.5)
print(av.get_mean_position(), av.get_n_points())
```

Two volumes give a distance distribution and a FRET observable:

```python
av2 = IMP.bff.get_av_from_pdb("structure.pdb", "A", 55, "CB", 20.0, 2.0, 3.5)
print(IMP.bff.average_distance(av.get_points(), av2.get_points()))
print(IMP.bff.mean_fret_distance(av.get_points(), av2.get_points(), 52.0))
```

Explicit dyes instead of volumes, with no IMP either:

```python
fret = IMP.bff.FRETRotamer()              # dye + linker from a rotamer library
sim = IMP.bff.ProbeDiffusionSimulation()  # a dye diffusing on its linker
```

Attaching a dye onto a structure and running Langevin dynamics on it
(`attach_probes`, `make_langevin_simulator`) is the one road that needs
`imp.bff`: those use IMP's integrators and hierarchies rather than
reimplementing them.

`ProbeSimulation` is the common interface: `get_positions`/`set_positions`,
`minimize`, `step` and `run`. `MolecularProbeSimulation` implements it using
IMP, while `ProbeDiffusionSimulation` implements a grid walk. A recorded run
returns `ProbeSimulationTrajectory`. The grid-specific `simulate(...)` method
retains its duration/diffusion arguments; `run(...)` takes a step count.

For C++ callers, include the matching `ProbeSimulation.h`,
`MolecularProbeSimulation.h` or `ProbeDiffusionSimulation.h`. The standalone
volume is `ProbeAccessibleVolume` in `ProbeAccessibleVolume.h`; IMP particles
use `ProbeAccessibleVolumeDecorator`. Probe templates, topology, force-field
CIF, potential tables and data locations live in the corresponding
`ProbeComponentTemplate`, `ProbeTopology`, `ProbeForceFieldCIF`,
`ProbePotentialTables` and `ProbeDataPaths` headers. Rotamer scoring helpers
are in `RotamerScoring.h`.


## Inter-label distance score usage:

First, import the module:

```python
import IMP.bff
```

Then, select the "score set", i.e., a set of distances that are used for 
score calculation from a FPS.JSON file:

```python
fps_json_fn = str(root_dir / "screening.fps.json")
score_set = "inter"
```

Finally, create the restraint and add it to the model.

```python
fret_rs = IMP.bff.probe_network_restraint_set(
    hier, fps_json_fn,
    mean_position_restraint=True,
    score_set=score_set
)
IMP.pmi.tools.add_restraint_to_model(hier.get_model(), fret_rs)
```

`probe_network_restraint_set` returns an `IMP.RestraintSet`: the network
restraint itself, or -- with `mean_position_restraint=True` -- one
mean-position restraint per distance over volumes made rigid-body members.
PMI's own bookkeeping stays with PMI: a `RestraintBase` subclass around this
set, for `output_objects` and the replica-exchange macro, is fifteen lines and
is written out in `examples/structure/t4l_pmi.py`.


The command tree is one tree. `imp_bff --help` lists every command:
`flexfit` and `rmsd` fit against distance restraints, `select-pairs` ranks
labelling pairs before an experiment is done, `openmm` writes a restrained
OpenMM run and `av-export` writes one volume for a viewer, `decays` runs the automated decay
analysis, `dye` is explicit-dye labelling and sampling, `rotamer` is
rotamer-library FRET, and `av-vs-rotamer` regenerates the comparison note. Two of those groups used to live *inside* the package and
were reachable only as `python -m IMP.bff.cgprobe.cli` and
`python -m IMP.bff.cli`; a click command is a decorated function, so a library
module carrying one cannot be imported without click, which is why command
trees belong in `bin/`.

## FRET-restrained docking

`imp_bff dock` samples rigid bodies against the distances in an fps.json with
PMI's replica exchange, writing an RMF trajectory, the best-scoring PDBs and a
per-distance CSV:

```bash
imp_bff dock -p bodyA.pdb -p bodyB.pdb -j labels.fps.json -o out/ \
             -c chi2_set -n 500 --mc-steps 10
```

One PDB per rigid body, `--fixed-body` anchors the frame. The Monte-Carlo
sampler is a command rather than a library call because it drives IMP.pmi
(Python-only) and writes a directory of results; scoring, minimising and
screening stay in the library -- `IMP.bff.score`, `dock_minimize`, `refine`,
`screen` -- because an application drives those with a cancellation callback
and reads the `DockingResult` back.

## Which FRET pair to measure next

`imp_bff select-pairs` ranks the labelling pairs of an fps.json by how much
each would tell you, given an ensemble of candidate structures -- the greedy
selection Olga performs (Dimura *et al.*, *Nat. Commun.* **11**, 5394, 2020):

```bash
imp_bff select-pairs -j labels.fps.json -r ensemble.rmf3 -n 10 -o ranking.tsv
```

Each step adds the pair leaving the smallest expected RMSD between the true
structure and the one the measurements would single out, so the output is an
order and a decay curve -- how much precision each further measurement buys,
and where buying more stops paying. The two inputs it needs are computed from
the ensemble rather than assumed: `ProbeNetworkRestraint::get_pair_efficiencies`
gives one row of FRET efficiencies per frame, and `pairwise_rmsd` the matrix
they are weighed against. `examples/labels/plot_pair_selection.py` plots the
curve. The kernels underneath -- `select_probe_pairs`, `expected_rmsd`,
`expected_rmsd_after_adding`, `chi2_right_tail` -- are public in their own
right, for a caller who already has the two matrices.

## FRET-restrained molecular dynamics

An accessible-volume network can score a structure; to *move* one it has to have
a gradient and a shape an integrator can live with.
`IMP::bff::AVFlatBottomRestraint` is that shape -- zero inside the experimental
error bars, harmonic outside them, and linear past that so the force is capped
and a badly-placed start cannot blow up the first step. It is the same well an
AMBER `&rst` record describes.

`md_flat_bottom_restraints` builds the whole system from an fps.json in one
call:

```python
system = IMP.bff.md_flat_bottom_restraints(
    hierarchy, "labels.fps.json", "chi2_C1_33p",
    f_max=15.0, tether_k=30.0, probe_mass=100.0, tether_atom="CA")

sf = IMP.core.RestraintsScoringFunction(force_field + [system.get_restraints()])
md = IMP.atom.MolecularDynamics(model)
md.set_scoring_function(sf)
md.optimize(10000)

for well in system.get_wells():          # one per measured pair
    print(well.get_name(), well.get_distance(), list(well.get_bounds()))
```

It computes both volumes, converts each measured distance and its **asymmetric**
errors into mean-position bounds, creates one probe particle per labelling
position tethered to its site, and puts a well between every measured pair.
`tether_atom` puts the probes on a different atom of the same residue -- `"CA"`
for a run that moves only alpha carbons, since a probe tethered to an atom
nothing optimises cannot follow the structure.
`examples/labels/plot_fret_restrained_md.py` drives T4 lysozyme from one
measured state to the other and back through the data alone.

The volumes are built once, on the starting structure, and every well follows
from those shapes -- an approximation that decays as the structure moves.
`AVRebuildOptimizerState` refreshes it:

```python
md.add_optimizer_state(IMP.bff.AVRebuildOptimizerState(model, system, 2000))
```

Every 2000 steps it resamples the volumes and re-derives each well's bounds and
each probe's tether from the geometry the trajectory has reached. On the T4L
example those bounds move by several angstrom over a run. Expect the *satisfied*
count to go down when you switch it on: a fixed well is a target derived from a
structure the run has already left.

### Running them in OpenMM

IMP has no bridge to an MD engine -- `IMP.modeller` is its only external-package
interface and Modeller is not one -- so `write_openmm_restraints` exports
instead:

```python
IMP.bff.write_openmm_restraints(system, "restraints.openmm.json", hierarchy)
```

The document carries the probes, their tethers and one flat-bottom bond per
measured pair, **in OpenMM's units** (nm, kJ/mol), and `energy_expression` is
the same well as a `CustomBondForce` expression --
`IMP::bff::openmm_flat_bottom_energy()`, one string used by both sides, so the
restraint that runs there is the restraint that ran here. The example verifies
that numerically without OpenMM installed (its `select` and `step` are two lines
of Python) and shows the dozen lines that consume the document.

`write_openmm_script` goes one further and writes the OpenMM *around* it: a
self-contained Python script with the restraint table embedded, which loads the
PDB, builds the system, adds the probes, their tethers and the wells, minimises
and runs. `imp_bff openmm` does it from the command line:

```bash
imp_bff openmm -j labels.fps.json -p structure.pdb -s chi2_set \
               --tether-atom CA -o fret_restraints.py
```

Attachment atoms are named by **chain, residue and atom**, not indexed: an index
depends on how the reader built its topology -- hydrogens, waters, altlocs --
and a spec that names atoms by position in someone else's file silently
restrains the wrong ones.

`imp_bff av-export -p structure.pdb -c A -r 132 -o site.pqr` writes one
accessible volume in whichever format the extension names.

The conversion behind step two is public in its own right:
`rmp_from_model_distance` slides two volumes apart until the modelled observable
matches a target and reports the mean-position separation where that happens.
An experiment reports \f$\langle R_{DA}\rangle\f$; a two-point restraint needs
\f$R_{mp}\f$; the two differ by several angstrom in a way that depends on the
shape of both clouds, so it is not a constant offset.
`rmp_flat_bottom_bounds` does the same for a measurement's value and both of its
error bars at once.

## Measuring and exporting a volume

| call | answers |
|---|---|
| `minimum_distance(a, b)` | closest approach of two clouds -- a *bound*, where every other distance type is an average |
| `cloud_overlap(points, refs, r)`, `av_overlap(av, hier, sel, r)` | how much of a volume's weight lies within `r` of a selection |
| `cloud_model_distance(a, b, type, R0)` | any `ProbePairMeasures` convention over two clouds -- except `pRDA`, which is not one |
| `cloud_distance_distribution(a, b, axis)` | `p(R_DA)` as a histogram, which is what `pRDA` names |
| `write_av(av, path)` | the volume as `.xyz`, `.pqr`, `.dx`, `.mrc`/`.map`/`.ccp4` -- the format follows the extension |

`av_overlap` takes a selection in either dialect this module speaks
(`"resi 130-134"` or `"resid 130 to 134"`). Note that a volume already excludes
the structure's van der Waals envelope inflated by the dye radius, so a contact
radius under about 5 A returns zero for everything.

`chain_weighting` in an fps.json position now reaches the volume: grid points
are weighted by the linker's chain statistics rather than uniformly. Read
`okf/validation/chain_weighting.md` before turning it on -- the shipped table
does not cover a dye-length linker, and the code says so.

# imp_bff_traj2bcif: convert a trajectory to BinaryCIF {#imp_bff_traj2bcif}

Converts a DCD or XTC trajectory to a BinaryCIF `_atom_site` coordinate
category, which is the format the shipped rotamer libraries use. Lossless
float32 by default: the FRETpredict pins are sensitive to dipole *directions*
between atoms about 1.7 A apart, so a quantisation grid that looks harmless as
a displacement is not one as an angle. `--grid` opts into quantisation for
corpora where that does not hold. Reading a DCD needs only `IMP.bff`; reading
an XTC needs `mdtraj`, which `imp_bff` already imports.

# imp_bff_traj2drot: build a rotamer library {#imp_bff_traj2drot}

Converts a trajectory (`.bcif`, `.dcd`, or `.xtc` through `mdtraj`) plus a
template PDB into a `.drot` rotamer library -- the internal-coordinate store
the shipped libraries use (PRD-118), written as `<stem>.drot.pto` in the PTO
container this stack shares with tttrlib's photon streams and chimol's
structures. `--bundle` files several of them into one container per family,
which is how the libraries ship: `dyes.drot.pto` (95 dye+linker ensembles),
`spinlabels.drot.pto` (10 spin labels) and `sidechains.drot.pto` (the
Dunbrack-2010 backbone-dependent table), each library addressed as
`dyes.drot.pto::A48_C1R_cutoff10`. A `.drot` carries its own template, atom
names, residue names and elements, so nothing has to sit beside it, and each
conformer is stored as its own base coordinates plus a bond length, bond angle
and dihedral per Z-matrix row: lossless to ~1e-6 A at ~0.65x the BinaryCIF of
the same ensemble, which is what keeps the FRETpredict parity pins inside
their tolerance. `--grid` opts into int16 grids (a third the size, ~3e-3 A)
for corpora where dipole directions do not matter. `--cluster A` takes the
raw-MD path: leader clustering, leaders become rotamers and cluster
populations become weights. `--all DIR` re-encodes a whole library directory.

# imp_bff_fps: FRET-restrained docking and screening {#imp_bff_fps}

The run modes of FPS, the FRET Positioning and Screening toolkit (Kalinin
*et al.*, *Nat. Methods* **9**, 1218, 2012), on IMP. FPS is a Windows C#
application; the physics it needs -- accessible volumes, the asymmetric
chi-square, rigid bodies, excluded volume -- is native here, and this program
is the door onto it.

```bash
imp_bff_fps score  -p bodyA.pdb -p bodyB.pdb -j labels.fps.json -c set
imp_bff_fps dock   -p bodyA.pdb -p bodyB.pdb -j labels.fps.json -c set -o out/
imp_bff_fps refine -p bodyA.pdb -p bodyB.pdb -j labels.fps.json -c set -o out/
imp_bff_fps screen -s ./library -j labels.fps.json -o ranked.csv
imp_bff_fps convert --positions LabelingPositions.txt -p bodyA.pdb -o labels.fps.json
```

`-p` is repeated once per rigid body and the order given is the body order;
`--fixed-body 0` (the default) holds the first still, so the first structure is
the frame everything is reported in. Every mode also reads the **legacy C# FPS
pair** directly -- `--positions LabelingPositions.txt --distances
Distances.txt` -- with no conversion step; `convert` exists only for when the
fps.json should be kept. The legacy formats are read-only by design.

Worked examples run on the shipped data in
`examples/structure/HIV_RT/`: HIV-1 reverse transcriptase with its DNA
primer/template (1R0A), FPS's own docking test case, with the original C# FPS
`LabelingPositions.txt` and `Distances.txt` beside the converted
`hiv_rt.fps.json` -- 11 labelling positions, 20 measured distances, two rigid
bodies.

```bash
E=$IMP/examples/bff/structure/HIV_RT
imp_bff_fps score -p $E/protein_1R0A.pdb -p $E/dna.pdb -j $E/hiv_rt.fps.json -c resolved
imp_bff_fps dock  -p $E/protein_1R0A.pdb -p $E/dna.pdb -j $E/hiv_rt.fps.json -c resolved \
                  -o dock_out --seed 1
```

Two things that file is honest about, and worth knowing before trusting any
docking result. Its `all` score set carries all 20 distances and scores `inf`,
because the two involving `p66_K287C` have no model value: that site is buried
at the protein-DNA interface and its accessible volume is **empty**. The
`resolved` set is the other 18. And FPS builds each volume on its own subunit
in isolation, where `IMP.bff` builds them in the assembled complex -- which is
why a site can be open in FPS and buried here.

The optimiser underneath differs from FPS's damped rigid-body dynamics, so
**scores are comparable with FPS and coordinates are not**. A single run is one
local minimum: repeat with different `--seed`s and compare before believing a
pose. See [`okf/prds/prd-121.md`](okf/prds/prd-121.md) for what has been
measured against FPS itself.

# imp_bff_fps_av: compute one accessible volume {#imp_bff_fps_av}

FPS's standalone AV dialog, which FPS itself reaches as `FpsGui -av`. One site,
one dye, one volume, written where PyMOL or VMD can open it (`.xyz`, `.pqr`,
`.dx`, `.mrc`).

```bash
imp_bff_fps_av -p 3GUN.pdb -c A -r 55 -d alexa488-long -o d55.xyz
imp_bff_fps_av -p 3GUN.pdb -r 132 -l 20 -w 4.5 --radii 3.5 -o d132.xyz
```

The dye presets are FPS's own, from its `Fps/data/linker.txt`. Each row carries
a single AV1 radius **and** three AV3 radii; they are not interchangeable, so
both are kept and `--av1` selects the single-radius column. The grid spacing
defaults to FPS's rule, `max(min(0.2L, 0.2W, 0.4R_i), 0.4)`, so a preset
reproduces FPS's grid.

Leave `--clearance` unset unless you know why you are setting it. The path
search inflates every obstacle by half the linker width, and a clearance below
that walls the source in and returns an **empty** volume; an unset clearance
derives one that does not. FPS's `LinkerInitialSphere x width` is a different
quantity — its seed is unconditional — so an FPS number does not transfer.
An empty volume is a legitimate answer for a buried site, and this program says
so and exits non-zero rather than writing a file that looks computed.

# imp_bff_fps_distance: the distances between two dye clouds {#imp_bff_fps_distance}

FPS's distance calculator. Two accessible volumes in, and the numbers a FRET
measurement is compared against out:

```bash
imp_bff_fps_distance d55.xyz a86.xyz            # R0 = 52 A
imp_bff_fps_distance d55.xyz a86.xyz -r 60 --json
```

```
Rmp       =    59.43 A     distance between the mean positions
<RDA>     =    61.72 A     mean donor-acceptor distance
sigma_DA  =     7.45 A     width of that distribution
<E>       =    0.292       mean transfer efficiency (R0 = 52 A)
<RDA>_E   =    60.28 A     what a measured <E> reads as
```

The three distances are **not** interchangeable and an experiment measures one
of them; `<RDA>_E` is always the shortest, because efficiency weights close
pairs as `1/r^6`. Reporting one where the data mean another is a
several-Angstrom error.

It reads both `.xyz` dialects. One trap it handles rather than hides: an FPS
`.xyz` is **duplicate-expanded** — its AV3 export writes a voxel once per dye
radius that fits, so the line count is the sum of densities and not the volume
(5473 lines over 3187 voxels in the shipped example). That is a weighting, so
the duplicates are kept and counted, and the unique voxel count is reported
separately.

# imp_bff_labelizer: score the labelling sites of a structure {#imp_bff_labelizer}

Scores every residue of a structure for how good a fluorescent-labelling site
it is, and optionally every pair of good sites for how informative a FRET
measurement between them would be. A 1:1 native port of the Labelizer model
(Gebhardt *et al.*, *Nat. Commun.* **16**, 3305, 2025), whose own
implementation cannot run on a current machine: its secondary structure shells
out to DSSP and refuses to start on macOS, and its default solvent-exposure
term shells out to MSMS binaries that are 32-bit ppc/i386 Mach-O. Both are
native here. Agreement with the reference's published output is measured, not
asserted -- on its own 1DDB example the secondary-structure and
cysteine-resemblance terms reproduce it exactly on all 195 residues and the
native residue depth matches MSMS to 0.14 A rms with no bias
([`okf/validation/labelizer_ab.md`](okf/validation/labelizer_ab.md)).

```bash
imp_bff_labelizer structure.pdb --conservation grades.txt --pairs
imp_bff_labelizer structure.pdb --donor Alexa488 --acceptor Alexa647 --pairs
imp_bff_labelizer structure.mmfdb.pto --show
```

The model, what a score means and the traps worth knowing are in
[Choosing a labelling site](@ref labelizer).

The published arithmetic is the **default, defects included**, because that is
what the paper's numbers were computed with; `--corrected` selects the
arithmetic the reference documents, and the container records which was used.
Conservation is imported, never computed -- pass a ConSurf `.grades` table or a
PDB carrying the grade in its B-factor column.

Output is one `.mmfdb.pto` rather than the reference's six CSVs, four
score-carrying PDBs, heat-map JSON and zip: the structure verbatim, the scores
as tables whose columns are named by MMFDB dictionary items, and the complete
settings. A position that was not scored carries a status and no number, which
is the thing the reference's output cannot express -- it writes `-1` for
"excluded" and `0` for "no contribution" into the same column as real scores.

# imp_bff_fps_export: FPS's result files, and the errors that fill them {#imp_bff_fps_export}

The door onto FPS's `SaveForm` -- the dialog every result file is written from
-- and onto `ErrorEstimation`, the mode that produces a table with more than
one row in it.

```
imp_bff_fps_export errors  ...   perturb the docked model's own distances,
                                 re-fit, and report the spread (FPS's
                                 parametric bootstrap)
imp_bff_fps_export table   ...   write the five export files from a results table
imp_bff_fps_export screen  ...   rank a library and write FPS's two Filter-mode tables
```

`errors` and `table` each write five files: one PyMOL script per row
(`<prefix><N>.pml`), every row as its own object (`<prefix>Overlay.pml`), one
object with a state per row (`<prefix>OverlayStates.pml`), the model distances
(`<prefix>Rtable_<T>.txt`) and the results grid (`<prefix>chi2table.txt`).

FPS's sixth output, `SimulationResults.bin`, is deliberately not reproduced: it
is a .NET `BinaryFormatter` dump, unreadable outside a legacy .NET target and
meaningless away from the project beside it. `results.json` replaces it, and
carries the molecule manifest that FPS's reload-by-position needed and never
had.

Run `imp_bff_fps_export --help` for the flags, including the three that
reproduce known FPS defects on demand.

# imp_bff_probe_pdb2cif: convert a probe PDB to mmCIF {#imp_bff_probe_pdb2cif}

Writes the `_atom_site` records for a dye structure, deriving the element from
the atom name where the PDB does not carry one. The conversion itself is
`IMP.bff.io.structure.convert_pdb_to_cif`; this is its command-line driver.

# imp_bff_potentials2pto: build the potential-table container {#imp_bff_potentials2pto}

Writes `data/potentials.pto`, the one container carrying every parameter table
the coarse-grained potentials read: the Miyazawa-Jernigan contact matrix, the
UNRES side-chain centroid table, the four-channel hydrogen-bond lookup and the
Ramachandran maps.

The tables arrive as four loose `.npy` files in a ChiSurf checkout and the
conversion is not a copy -- NaNs, a short-range divergence and two channels
that are coordinate grids rather than maps are dealt with on the way in, and
what was done is written into the container's manifest. Miyazawa-Jernigan and
UNRES go in as text in IMP's PMF format, so `IMP.core.StatisticalPairScore`
reads them directly.

```bash
imp_bff_potentials2pto --database ../chisurf/chisurf/core/structure/potential/database
```

Rebuilding the container is a maintainer's job: the shipped one is in the
repository and a user never needs to run this.

# Info

_Author(s)_: Thomas-Otavio Peulen

_Maintainer_: `tpeulen`

_License_: [LGPL](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html)
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

_Publications_:
- None
