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
efficiencies (`IMP.bff.RotamerFRET`, `IMP.bff.attach_dyes`, ...). Every public
name is reachable flat as `IMP.bff.<Name>`; the manual page
`doc/manual/structure/structure_cgdye.ipynb` walks the explicit route, and
`okf/cgdye.md` records how the code is organised.


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
`flexfit` and `rmsd` fit against distance restraints, `decays` runs the
automated decay analysis, `dye` is explicit-dye labelling and sampling,
`rotamer` is rotamer-library FRET, and `av-vs-rotamer` regenerates the
comparison note. Two of those groups used to live *inside* the package and
were reachable only as `python -m IMP.bff.cgdye.cli` and
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

# imp_bff_dye_pdb2cif: convert a dye PDB to mmCIF {#imp_bff_dye_pdb2cif}

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
