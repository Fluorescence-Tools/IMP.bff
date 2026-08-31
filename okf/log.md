# Update Log

## 2026-08-31 (3) — the OpenMM run, rebuilt volumes, and pRDA as a distribution

Follow-ups on the six, all from the user reading the previous entry.

**A runnable OpenMM script, not just a document.** `write_openmm_script` writes
a self-contained Python file -- the restraint table embedded, so the generated
file is the whole input -- that loads the PDB, builds the system, adds the
probes, their tethers and the flat-bottom wells, minimises and runs, then
reports how many measured distances the run satisfied. `imp_bff openmm` does it
from the command line. OpenMM does not have to be installed to write it, and the
test parses the generated file with `ast` and checks its embedded table against
the system it came from.

**The volumes rebuild during a run.** `md_flat_bottom_restraints` derives every
well from volumes computed once on the starting structure, and that
approximation decays: as the fold moves, a volume's shape and its offset from
its attachment change, and so does the R_mp that reproduces a measured
<R_DA>. `AVRebuildOptimizerState` resamples every `period` steps and re-derives
each well and each tether. On the T4L example the bounds move **up to 47 A**
over 8000 steps, which is the size of the error being carried when they are
held fixed.

It makes the example's headline *worse* -- C2 goes 23 -> 31 of 33 rather than
23 -> 33 -- and that is the point: a fixed well is a target derived from a
structure the run has already left. The example says so.

**`pRDA` answers as a distribution.** Refusing it as a scalar was half the fix;
`cloud_distance_distribution` is the other half -- the weighted histogram of
pair distances over caller-supplied bin edges, which is what a decay is a
function of. Its mean agrees with the scalar mean-distance reduction to 0.002 A.

**Two commands**: `imp_bff openmm` and `imp_bff av-export` (one volume, format
from the extension).

One design fix found by a test that should have passed and did not: the rebuild
updated its **own copy** of the probe records, so `ProbeParticle::offset` went
stale while the tether it described moved. Two sources of truth for one number.
`offset` is gone; `MDRestraintSystem::get_tether_length(i)` reads the tether,
which is where it lives.

Also cleared, from the previous entry's audit: `pRDA` returning R_E, `Rmin`
reaching the vocabulary but not `av_distance`, `write_av`'s `.dx` branch writing
the obstacle raster rather than the volume, and the chain weighting reaching the
point cloud but not the exported grid.

Tests: `test/restraints/test_flat_bottom_md.py` 21 -> 29,
`test/label/test_av_kernels_and_io.py` 23 -> 30.

## 2026-08-31 (2) — the six things Olga had that this module did not

A survey of Olga's tree (17 evaluators, `AV/`, `av2restraints`, `screen-nox`,
`irmsd`) against this module found most of it covered or bettered here --
screening, the selection language, AV1/AV3, ACV (richer than their `freeSize`),
distance distributions, RMSD. Six gaps were real. All six are closed.

**1. FRET restraints an MD engine can integrate.** `AVMeanDistanceRestraint`
scores a chi-squared: it pulls at every separation and never stops, which is a
scoring function, not a potential. `AVFlatBottomRestraint` is AMBER's `&rst`
well -- zero inside the error bars, harmonic outside, **linear** past that so
the force is capped and a badly-placed start cannot blow up the first step.
`md_flat_bottom_restraints` builds the whole system from an fps.json: volumes,
bounds, one probe particle per position tethered to its site, one well per pair.
`tether_atom` moves the probe to another atom of the same residue, which a
C-alpha-only run needs.

The bug that would have been silent: the probes were created with mass and
coordinates and **not marked optimizable**, so every well pushed on something
that could not move and MD reported no error and no motion. The test now asserts
`get_coordinates_are_optimized()` per probe.

`examples/labels/plot_fret_restrained_md.py` drives T4 lysozyme through both
measured states. The crystal satisfies 18 of the 33 C1 distances and 23 of the
33 C2 ones; restrained MD under C2 reaches **33/33**, and after switching the
restraints to C1, **33/33** again. Two things had to be right for that: the
elastic network is re-anchored on the C2 conformer before the C1 phase (left on
the crystal it is a memory of where the run started, and C1 stalls at 29/33),
and f_max/tether have to out-pull the network -- too soft and the data is
decoration, too stiff and the force field is.

The export target is **OpenMM**, not AMBER (user, 2026-08-31), and IMP turns out
to have no MD bridge at all to hang it on: `IMP.modeller` is its only
external-package interface, Modeller is not an MD engine, and the only
`openmm` strings in the whole IMP tree are the two files this added.
`write_openmm_restraints` therefore writes a document -- probes, tethers, one
flat-bottom bond per pair -- in nm and kJ/mol, and
`openmm_flat_bottom_energy()` is the well as a `CustomBondForce` expression.

One string, used by both sides, and **checkable without OpenMM installed**:
its `select` and `step` are two lines of Python and `^` is `**`, so a test
evaluates the exported expression and compares it against the C++ restraint
over 5-120 A through both walls and both linear tails. They agree to
**2e-13 kJ/mol**, conversion included. That test is the thing that keeps the
two from drifting.

Attachment atoms are exported **by chain, residue and atom name**, never by
index: an index depends on how the reader built its topology -- hydrogens,
waters, altlocs -- and a spec that names atoms by position in someone else's
file silently restrains the wrong ones.

**2. `rmp_from_model_distance`.** The inverse of the model distance in the
separation of two clouds, which is the number that turns a measurement into a
two-point restraint. The pairs are sampled **once** and then translated, so the
function being bisected is smooth -- resampling per iteration would put
1/sqrt(n) noise on a root being found to 0.01 A. Bracket then bisect, which is
stricter than the reference's fixed-point step (it assumes dR_model/dR_mp = 1
and stops when it stops improving). A target no separation reproduces raises
rather than returning a restraint that pulls forever.
`rmp_flat_bottom_bounds` converts a value and both error bars, each separately,
because the relation is not linear and a symmetric error bar in the measured
quantity is not symmetric in R_mp.

**3. `chain_weighting` was declared and never read.** A schema field,
documented, with no AV attribute and one grep hit in the tree: a file asking for
it got the unweighted volume and no complaint. Implemented as a `LinkerWeighting`
folded into `density_soa_` once per raster, so every reader sees it. **And the
shipped table does not cover a dye linker** -- 0.02 % of its weight lies within
20.5 A, so it selects the volume's rim and moves the mean position 13.7 A. The
code measures that and warns.
[`chain_weighting.md`](validation/chain_weighting.md).

Three defects in the above came out of auditing the paths the first round of
tests did not reach, all of them silent:

* **`pRDA` returned a number.** `cloud_model_distance` with the distribution
  type fell through to the efficiency branch and returned R_E -- a value that
  looks like an answer and is a different quantity. It throws now, naming
  `av_distance_distribution`.
* **`Rmin` reached the vocabulary but not the distances.** An fps.json could
  declare it and `av_distance`/`av_distance_quadrature` would silently give the
  mean distance. Both know it now, exactly and without a quadrature, since a
  coarsened cloud gives a *wrong* minimum rather than a cheaper one. Inverting
  it is refused: a minimum from a sample is biased high.
* **`write_av`'s `.dx` branch wrote the obstacle raster.** `PathMap` derives
  from `em::SampledDensityMap` and its `get_value` is the map it was *built
  from*, not the accessible volume. Every OpenDX file this module wrote for an
  AV was the wrong map. It reads `PM_TILE_ACCESSIBLE_DENSITY` now, the same
  feature `write_map_feature` and `AVBuilder` use.

And the chain weighting was landing on the point cloud but not on the tile-based
grid, because a fold guarded by an "applied" flag cannot be sequenced reliably
against a multi-stage resample that copies tiles once. It is a multiply at the
point of read now (`PathMap::path_weight`), which has no state to desynchronise.
Details in [`chain_weighting.md`](validation/chain_weighting.md).

**4. `minimum_distance`.** Closest approach of two clouds -- a bound where every
other type here is an average, and `Rmin` joins the fps.json vocabulary. Exact,
not sampled: a minimum estimated from a sample is biased high. Zero-weight
points are skipped, which matters since radius-zero masking leaves them in.

**5. Export formats.** `write_av(av, path)`, format from the extension: `.xyz`,
`.pqr`, `.dx` (OpenDX), `.mrc`/`.map`/`.ccp4`. Two overloads, for the decorator
and for the value `compute_av` returns.

**6. `av_overlap`.** How much of a volume touches a selection, in either
dialect. Its documentation now carries the thing that cost a test: a volume
already excludes the van der Waals envelope inflated by the dye radius, so
**nothing is within about 5 A of an atom centre** and a smaller contact radius
returns zero for everything.

Structural: the `ProbePairMeasures` enum and the fps.json distance vocabulary
moved from `AV.h` to `AVDistance.h`, where the distance conventions belong, and
`avdistance.i` is now wrapped before `av.i` because AV.h's defaults name that
enum. `States`/`AccessibleVolume` gained `*_vector` accessors for C++ callers,
`%ignore`d so Python keeps only the numpy views.

Tests: `test/restraints/test_flat_bottom_md.py` (21) and
`test/label/test_av_kernels_and_io.py` (32). Suite 1366 -> **1471 passed**, 5 skipped, 3
xfailed, 180 subtests. Both new examples pass the docs sweep.

## 2026-08-31 — greedy Olga reaches a structure: `select-pairs`, an example, two kernels

The greedy pair selection was correct and unreachable. `select_informative_pairs`
and its three kernels had been in C++ since the `%pythoncode` sweep (PRD-117),
with 28 passing tests -- all of which handed it **synthetic matrices**. Nothing
in the tree turned an ensemble into the two things it takes, so a user with a
trajectory and a candidate pair list had no route in. No example, no bin.

Closed from both ends, and in C++ rather than twice in Python:

* **`ProbeNetworkRestraint::get_pair_names` / `get_pair_efficiencies`.** The
  score set's pairs and their mean FRET efficiency at the current coordinates.
  The order was the real content: a `std::map` iterates sorted, and *which
  column is which pair* was about to be rediscovered by every caller. It is
  stated on the class instead. `get_pair_efficiencies` re-evaluates first, so a
  freshly loaded frame answers for itself -- the failure it prevents is silent
  (the previous frame's volumes, read as this one's).
* **`pairwise_rmsd`** (`Clustering.h`), beside `rmsd_no_align` and on the same
  `(n_frames, n_atoms, 3)` typemap `cluster_frames_leader` already uses. It
  calls `compute_rmsd`, so the Kabsch superposition has one implementation, and
  copies each frame once rather than once per pair.

Then `imp_bff select-pairs` (RMF or PDBs in, ranking TSV out, `--stride` for a
long ensemble) and `examples/labels/plot_pair_selection.py`, which plots the
decay curve over the T4L docking ensemble and one 33-pair score set. Both call
the same two accessors; neither reimplements the marshalling.

`expected_rmsd` -- Olga's `rmsdMeanMean` -- had been exported with **no caller
and no test** since the port. It is now the "before any data" baseline both the
bin and the example print, which is what makes a gain column mean anything.

`test/restraints/test_pair_selection.py`: 14 tests, both entry paths (RMF and a
stack of PDBs), 4 s. Two things it found:

* **The decay is not monotone.** Both the bin and the example print a gain per
  measurement, and the obvious assertion -- that it never goes negative --
  holds on the 20-frame ensemble and fails on a four-member one, where the
  third selection gives back 0.03 A. The greedy scores each candidate against
  the chi-squared accumulated so far, so the best *remaining* candidate can
  leave a larger expectation than the step before it. The tests assert what is
  actually true (the first selection beats the uninformed baseline, the last
  beats the first) and say why the stronger claim is not asserted.
* **A dropped `IMP.Model` is a segfault, not an exception.** The first version
  of the restraint fixture returned the restraint and the hierarchy but not the
  model; the model is reference counted from Python, and `get_pair_efficiencies`
  faulted inside SWIG. Worth knowing before writing the next fixture.

`test_superposition_removes_the_rigid_motions_and_nothing_else` is the one that
earns its length: a frame rotated *and* translated must come back at RMSD 0
with `superpose=True` and far from it without -- the only check exercising the
reflection guard `compute_rmsd`'s Kabsch has.

**A/B against Olga itself** -- [`greedy_olga_ab.md`](validation/greedy_olga_ab.md),
`test/restraints/test_greedy_olga_ab.py`. Olga is Qt + pteros, but its selector
is not: two non-semantic edits (drop the pteros includes and everything after
`sys2xyz`; add the `template` disambiguator clang wants at `spline.hpp:165`)
compile `best_dist.h` standalone against Eigen, and the test builds it in a
temporary directory. On the T4L ensemble, 50 frames x 33 pairs: **the same ten
pairs, in the same order**, decays agreeing to 1.6e-4 A. Same on a
well-separated synthetic case.

Where they differ is worth knowing. Olga carries `MatrixXf` -- float32 -- and
does not evaluate the chi-squared right tail; it fits it with a 64-piece
quadratic spline, which at **one** degree of freedom is off by **0.11 in
probability** near the origin. `chiSqRTSpline`'s own accuracy assert is
commented out with an `|| ndof==1` escape for exactly this. And ndof is 1 for
the *first two* selections, because `bestPair` uses
`max(selPairs.size() - 1, 1)` -- so the one step where the decays visibly
disagree (3.2e-3 A at step 1) is the reference's approximation, not ours.

A third case was built to break the parity and does: 25 candidates cut from one
cloth, where the second selection's best two are 2.6e-4 A apart out of 1.70.
The two pick differently and swap back one step later -- same eight pairs,
different order. The test asserts what can be asserted, that a disagreement
happens only on a tie, and only checks the *first* divergence, since after it
the accumulated chi-squared differs and the scores are no longer comparable.

Still open, and the owner's call (PRD-120): whether the Labelizer's pair layer
should rank with this rather than with its own published score.


## 2026-08-26 (night, last +7) — a use-after-free that returned zeros

Chasing the last of the untested surface -- the two `bin/` programs nothing
covered, and the top-level commands -- turned up two real defects and one that
had been silently wrong for as long as the CIF writer has existed.

### A system CIF written anywhere else could not be read back

The writer records each component's MOL2 path *relative to the CIF*, and
computed the `..` count from the **separators** of the base directory rather
than its segments -- one short whenever the base has no trailing separator,
which is always:

```
base   /var/folders/cl/xxxx/T/tmp1        6 segments
mol2   /Users/me/dye.mol2
wrote  ../../../../../Users/me/dye.mol2   5 ups
read   /var/Users/me/dye.mol2             does not exist
```

Nothing noticed because every test writes its systems into the data tree,
where the common prefix is long and the miscount cancels. `imp_bff simulate`
was dead for any system built outside it. `relative_path` in
`TopologyBuild.cpp` counts segments now, and
`test/cgdye/test_system_paths.py` writes a system into a `tmp_path` -- which
is exactly the case that failed.

### A container read off a temporary was freed memory

Writing the test for that found something worse:

```python
IMP.bff.read_forcefield_cif(path).components.values()
```

returned `'\x00\x00\x00...'` where a path should be, and segfaulted under
pytest. An accessor returning `const std::map&` gives SWIG a **borrowed**
pointer into the owner, and the proxy does not keep the owner alive: the
system is a temporary, it dies as the expression unwinds, and the map proxy
outlives it. Silently, most of the time.

The accessors cannot return by value -- `Scoring.cpp` calls `get_bonds()` in
the *condition* of a loop over the bonds, so a copy per access is quadratic --
so the copy is made **at the language boundary**: `%owned_container_out` in
`IMP_bff.types.i` gives every `const container&` return an owning copy, which
costs one copy per Python access and hands the lifetime to Python. Eight
container types: the maps, and the vectors whose element is one of this
module's own values.

**Not** `std::vector<std::string>`, `<double>` or `<int>`. Those already come
back *converted* -- a Python list, a numpy array -- which is a copy, so they
were never borrowed; claiming them replaced the conversion with a raw proxy and
five tests stopped recognising their own results (`assert <RMF_HDF5.Strings> ==
('SD',)`). The first cut of this fix did claim them, and the suite said so.
`%naturalvar` is on beside it, for the same hazard in plain member variables.

`test/test_borrowed_containers.py` reads five kinds of container off a
temporary, and checks that the copy really is one.

### Two programs nothing had ever run

`imp_bff_dye_pdb2cif` passed `--dye-id`'s `None` default into a
`const std::string&`, so its documented default invocation raised a
`TypeError` before reading anything. `imp_bff_labelizer --show`, given a
structure instead of a container, failed with "does not begin with an EBML
header" -- true, and useless. Both fixed, both now in
`test/io/test_bin_programs.py`.

`imp_bff dock` runs (score 188.30 over 99 distances, 17 accessible volumes),
and so do `build-system`, `simulate` and the six `dye` commands.

## 2026-08-26 (night, last +6) — nine dead paths, found by running the docs

The unit suite was green at every step of the last three passes, and the CLI
and the example gallery were quietly rotting the whole time: **what breaks when
a C++ port changes a shape is the code nobody runs**, and nobody runs the
documentation.

Running every `bin/imp_bff` command and every example and notebook found nine.

### In the program

| command | what was wrong |
|---|---|
| `dye sample-dof-walk` | imported `LinkerSampler`, a name from before the linker sampler became C++; read `resolve_probe_site` as a dict; and then, once running, **accepted nothing** |
| `dye sample-langevin` | passed `timestep_fs=None` into a `double`, and asked `run()` to write an RMF, which it stopped doing when the sampler became C++ |
| `dye label-fusion` | `sys.exit` in a module that never imported `sys` |
| `flexfit` (FP branch) | `attach_probes` with four parallel lists, and another dict-style site |

`sample-dof-walk` is the interesting one. It rejected every proposal because
the probe is *bonded into* the site: its first atoms sit a bond length from
residues i-1 and i+1, and the walk counted those bonds as clashes. Excluding
the site and its two neighbours, and accepting a move that does not make the
clash count *worse* -- a walk that starts inside the protein can never leave,
otherwise -- gives 30/60 and 51/60 accepted with no clashes left. It also
built its geometry from a MOL2 while turning atoms read from a **PDB**, whose
order need not match; both come from the same file now.

### In the documentation

| page | what was wrong |
|---|---|
| `plot_convolution_routines`, `plot_pile_up` | `IMP.bff.decay_fconv*` -- decay convolution is **tttrlib's** by the placement rule; they call `tttrlib.fconv`/`fconv_per`/`fconv_simd` now |
| `hgbp1_label_and_sample` | `{"N", "CA", "C"} <= set(site)`, and `attached[0]["resnum"]` |
| `langevin_hgbp1_site481` | the `None` timestep and `out_rmf` |
| `plot_k2_uncertainty` | `k2_call` returned two values where its caller unpacked three -- it had *never* matched; `Kappa2Distribution` supplies the third |
| `t4l_pmi` | fetched AVs for a display helper that went with `representation/av.py` |
| `structure_cgdye` | `attach_probes` tuples, `SITE_KEEP_ATOM_NAMES` (a function now), a nested transition-count matrix (flat + `n` now), and a `RotamerLibrary` read as a dict |
| `structure_accessible_volumes` | stopped on a missing `ipyvolume`, which is a viewer and not a dependency |

### What holds them now

`test/expensive_test_docs_and_examples.py` runs **every** example and **every**
notebook
-- 42 of them -- with each notebook's own directory as the working directory,
which is what Jupyter does and what their relative paths assume. Opt-in, like
the other `expensive_test_` files. `test/label/test_dye_commands.py` runs six
commands for real, and asserts the walk *accepts* something: a sampler that
accepts nothing passes any "does it run" check.

### One decision left

`doc/manual/decays/decay_curves.ipynb`, `decay_forward_model.ipynb`,
`decay_objective_function.ipynb` and `doc/manual/programming/programming_imp_decorator.ipynb`
document `IMP.bff.DecayCurve`, `DecayConvolution`, `DecayLifetimeHandler`,
`DecayScale`, `DecayPattern` and `DecayLinearization` -- 50 references to an
API **deliberately deleted** by `93b7198` ("PRD-113 stage 0: delete the TCSPC
instrument layer"), and which exists nowhere in the stack now, tttrlib
included. They are not broken call sites; they are pages describing a layer
that left. Porting them to tttrlib's manual or dropping them is a call for the
maintainer, so they are skipped by name in the sweep, with the reason written
where the skip is.

## 2026-08-26 (night, last +5) — two landmines defused

**`.drot` could not read a payload that compressed well.** `DrotReader`'s
decompressor sized its output buffer as `compressed * 8` and grew it on
`NEEDS_MORE_OUTPUT`; brotli's one-shot decoder reports a too-small buffer as
`BROTLI_DECODER_RESULT_ERROR`, so the growth branch never fired and a good file
came back as "corrupt brotli stream". Every shipped `.drot` is under the ratio,
which is why nothing had noticed -- but a library of near-identical conformers
(a rigid dye on a short linker) is not. It streams now.

**Four `attach_probes` calls in `bin/imp_bff` were dead.** Three passed
`(hierarchy, chain, residue)` tuples and one passed four parallel lists; the
C++ takes `ProbeAttachment` values, so all four raised `TypeError` at the
binding. `dye sample-rotamer` then read `attached[0]["site"]` out of a list of
values. The command runs end to end now, and `resolve_probe_site` -- which
answers CA, N, C -- is where the site comes from.

Neither had a test, which is why neither was known to be dead. Both do now:
`test/io/test_drot.py::test_a_payload_that_compresses_hugely_still_reads` and
`test/label/test_attach_probes_cli.py` (3 tests, including the command itself
through `CliRunner`).

## 2026-08-26 (night, last +4) — the coarse-grained potentials, as IMP scores

The family that lived twice: `numba.njit` kernels in
`IMP.cgmol.{statpot,sterics,solvation}` (imp-tricks) and classes around them in
chisurf's `structure/potential/potentials.py`, which imported those kernels.
Coordinates are this module's business and numba is a prototyping tool in this
stack, so both are here and both are C++.

### The shape, corrected

The first cut of this was a `CoarseStructure` holding flat coordinate arrays, a
residue lookup table and a C-alpha distance matrix -- chisurf's design,
transliterated. That is a second scoring world beside IMP's, and it was thrown
away. **IMP has all three already**: particles in a `Model`, a `Hierarchy` over
them, and a `ClosePairContainer` for the pair loop with its cutoff.

| term | what it is here |
|---|---|
| Miyazawa-Jernigan | `IMP::core::StatisticalPairScore<atom::ResidueType>` + a table |
| UNRES centroid | the same, distance-binned |
| hydrogen bond | `HydrogenBondRestraint` (four-channel, no IMP form) |
| Go model | `GoRestraint` (truncated LJ, no IMP form) |
| generalized Born | `GeneralizedBornRestraint` (no IMP form) |
| Ramachandran | `RamachandranRestraint` (no IMP form) |
| LJ bead | `LennardJonesBeadPairScore` |
| clash | **IMP's** `SoftSpherePairScore`, k = 2/t^2, via `build_clash_restraint` |
| C-alpha internals | **IMP's** `Harmonic` distance/angle/dihedral, via `build_ca_internal_restraints` |
| accessible surface | this module's `solvent_accessible_surface_area`, per residue |
| radius of gyration | **IMP's** `atom::get_radius_of_gyration`; nothing added |

`add_residue_type_score_data` is `IMP::atom::add_dope_score_data` for residues:
it puts the type on one atom each and the scores read it through
`get_residue_type_key()`. The tables name residues and IMP's PMF reader maps
names to indices, so **nothing here carries a residue order** -- which is what
kept the two Python copies from agreeing.

The kernels stay public on flat arrays (`clash_energy`, `go_energy`,
`lennard_jones_bead_energy`, `generalized_born_energy`, `residue_asa`,
`ramachandran_energy`), because a caller holding numpy should not have to build
a model to score one, and because they are what the parity tests pin.

### Verified against the Python it replaces

Every kernel reproduces the unmodified numba kernel on a seeded 12-residue
structure, to the last bit or one ulp:

| | value |
|---|---|
| clash | 35.961679525439735 |
| clash (t 6, cov 1) | 6.521641123970415 |
| Go | -6.676356712583448 |
| residue ASA (64 points) | 288.3883881225005 |
| Miyazawa-Jernigan | 21 contacts, 3.95 |
| hydrogen bond | 3 bonds, -0.1324432499910165 |
| generalized Born | -0.024505222258838378 |
| LJ C-alpha | 91590.40855879632 |

The Miyazawa-Jernigan row is the interesting one: 3.95 comes back through
`ClosePairContainer` + `PairsRestraint` + `StatisticalPairScore`, so IMP's pair
loop and IMP's table reader agree with the hand-written double loop.

### Three defects, fixed and recorded

`okf/validation/hbond_ca_cutoff.md`: the H-bond kernel compared a plain
C-alpha distance against a **squared** cutoff, so its 8 A prefilter passed
every pair -- and unlike the same mismatch in the MJ and ASA kernels, there it
decides which bonds count. Applied here, with the difference pinned. Also:
chisurf's `go()` carried its accumulator across pairs, and `_asa_kernel` read
its representative atom from a column that holds `CG`.

### Not ported

`IMP/cgmol/martini/score.py` -- seven PMI `RestraintBase` subclasses over a
vermouth topology. Neither has a C++ spelling; their own pass.

### The tables: one container

`data/potentials.pto`, 680 KB, holding all four where four loose `.npy` files
(3.3 MB) used to sit in a chisurf directory. `MiyazawaJerniganPairScore()` and
`UNRESCentroidPairScore()` with no argument read it.

One file rather than four because a potential is not one table: the UNRES term
needs its grid *and* the residue order that indexes it, the Ramachandran map
needs to say which channel is glycine, and a reader that has to find four files
in agreement with each other will one day find three. The manifest carries what
each is, where it came from (with a sha256 of the source) and **what the
conversion did to it** -- because none of them were copies:

| table | as shipped | what the converter had to do |
|---|---|---|
| `mj` | PMF text | flat to 6.5 A, zero past it; two bins, because IMP's reader splines each pair |
| `unres` | PMF text | 6 NaNs in bin 0, and the flat repulsion below 3.5 A that the Python applied *in its loop*, written into the bins |
| `hbond` | 4x800 grid | 243 values above 1e5 below 1.3 A -- a divergence where no two atoms are -- held at the value there |
| `ramachandran` | 3x360x360 grid | two of the source's five channels are the phi and psi *coordinate* grids, not maps; and 114 556 cells of the proline map held a sentinel `1.0`, which read as a density would make the unmeasured region the most favourable place in the map |

The two contact potentials go in as **text in IMP's PMF format**, so
`IMP::core::StatisticalPairScore` reads them directly and the potential is the
file rather than a loop. The file names residues and IMP maps names to indices,
which is why nothing in the C++ carries a residue order.

### A fourth defect, found by the conversion

`unres.npy` is filled in its **upper triangle**: past the repulsion, entries
with `type_i < type_j` average -0.129 and their transposes -0.003. `centroid2`
indexes it by the residues' order *in the chain*, so about one contact in two
read the empty half and scored zero -- on T4 lysozyme, -207.198 where the
filled half gives **-371.031** over the same 4 983 pairs. A PMF file holds one
value per unordered pair, which is what a pair potential is; the container
stores the filled half and the manifest says so.

### The residue types are the potential's own, not IMP's

The first cut typed particles with `IMP::atom::ResidueType`'s index and wrote
tables as wide as IMP. That is wrong in a way that only shows up later: an IMP
`Key` table **grows at run time** -- read one PDB with an unusual residue in it
and the next index is one further out -- while a statistical table is a fixed
square. The full suite caught it: a test that passed alone failed after other
tests had widened the key table.

`IMP::bff::ResidueContactType` is a key family of its own, the shape
`IMP::atom::DopeType` has, and its names are the ones the table names.
`add_residue_type_score_data` types a residue only when the table covers it, so
no index can land outside the square. Two smaller things fell out: the key's
family number has to be one nobody else uses (`783462` is
`IMP::atom::ProteinLigandType`, and sharing it brought back two hundred names
that were not residues), and a `%template` base has to be spelled the way the
header spells it or SWIG quietly makes the class a bare `object` -- it built,
imported, and then would not go into a `PairsRestraint`.

### Two things the build taught, worth keeping

- **A payload that compresses better than eight times did not read.** Brotli's
  one-shot decoder reports a too-small output buffer as an *error*, not as
  `NEEDS_MORE_OUTPUT`, so a reader that guesses the size as `compressed * 8`
  fails outright -- and 2.4 MB of digits compresses far better than that. The
  container writes the decompressed length into the payload and allocates
  exactly. **`src/DrotReader.cpp` guessed the same way and is fixed too**: it
  uses brotli's *streaming* decoder now, which asks for more room when it runs
  out, so there is nothing to guess and nothing to re-decode. A library of
  near-identical conformers -- a rigid dye on a short linker, which is the
  common case -- would have been reported as a corrupt stream while being
  perfectly good; `test_a_payload_that_compresses_hugely_still_reads` writes
  one.
- **A new program in `bin/` silently disables the module.** IMP's configure
  wants a section title for every command-line tool in `README.md`; without one
  it prints "Module IMP.bff disabled" among a thousand lines and the next
  `ninja` fails with `library 'IMP.bff-lib' not found`, which says nothing about
  the cause. `imp_bff_potentials2pto` has its section now.

### Delivered

`include/Potentials.h`, `src/Potentials.cpp`, `pyext/IMP_bff.potentials.i`,
`test/potentials/test_potentials.py` (22 tests),
`examples/structure/coarse_grained_potentials.py`,
`doc/manual/structure/structure_potentials.ipynb`,
`okf/validation/hbond_ca_cutoff.md`.

## 2026-08-26 (night, last +3) — `AVNetwork` is `ProbeNetwork`; the template collision, fixed at the source

### The rename

| was | is |
|---|---|
| `AVNetworkRestraint` / `AVNetworkRestraints` | `ProbeNetworkRestraint` / `ProbeNetworkRestraints` |
| `SimpleAVNetworkRestraint` | `SimpleProbeNetworkRestraint` |
| `av_network_restraint_set` | `probe_network_restraint_set` |
| `AVNetworkRestraint.h`, `.cpp` | `ProbeNetworkRestraint.h`, `.cpp` |
| `test_AVNetworkRestraint.py`, `plot_AVNetworkRestraint.py` | `test_ProbeNetworkRestraint.py`, `plot_ProbeNetworkRestraint.py` |

Same reasoning as `ProbeSite` and `ProbeRestraints`: what the restraint scores
is a network of *probes* -- the volumes are how a dye's reach is modelled, not
what the restraint is about, and a spin label has no accessible volume in the
FRET sense at all. The four old names are on `test_public_api_names.py`'s
retired list. Prose in `okf/prds/`, `okf/handoff-*` and this log keeps the old
names where it is recording history; `AVNetworkRestraintWrapper` keeps its name
wherever a sentence says what it *was*, because that was its name.

`AVMeanDistanceRestraint`, `AVModel`, `AVBuilder` and the rest of the AV layer
are untouched: those really are about accessible volumes.

### The template collision

Yesterday's `rmf` dependency cost four `%template` names -- `rmf` brings `isd`,
`isd` brings `saxs`, and `saxs` declares
`%template(DistBase) std::vector<double>` first, so `types.i`'s
`%template(VectorDouble)` was silently skipped (warning 302 is filtered by
IMP's own `%warnfilter`). Three test files papered over it by importing
`IMP.saxs.DistBase`. The paper is off now, and the cause is gone:

**Nothing in this module's surface should ever have asked a caller to build a
wrapped vector.** Two functions did, through non-const `std::vector<double>&`
out-parameters, and both are fixed at the source:

- `wobbling_kappa2_distribution` and `wobbling_kappa2_distribution_delta`
  returned the values and *filled* `k2_scale` and `k2_hist`. They return one
  **`Kappa2Distribution`** value now -- `values`, `scale`, `hist`, three numpy
  attributes -- which is what `RotamerScoreResult` and every other result value
  in this module already does.
- `diffusion_propagate` filled a `std::vector<double>& fluorescence` beside its
  managed density view. It has **two managed views** now and hands Python a
  `(fluorescence, density)` tuple. The C++ callers in `GridDiffusionSolver`
  free both, which is the same contract they already had for the density.

`types.i` no longer *declares* `VectorString`, `VectorDouble`, `VectorFloat` or
`VectorInt`: a `%template` that SWIG always skips is a promise of a class that
never appears. The comment there says which module owns each type.
`VectorLong` stays -- nothing else wraps `std::vector<long>`.

`test_vector_input_typemap.py` still builds an `IMP.saxs.DistBase`, and should:
binding a wrapped vector is the property that test exists to check, and the
class is whatever SWIG has. `examples/spectroscopy/plot_k2_uncertainty.py`
imported `IMP.bff.VectorDouble` and had been **broken since the dependency
landed**; it uses the new value.

### Verification

`ninja IMP.bff-python` clean after a cmake reconfigure (renamed headers are new
headers as far as SWIG's dependency list is concerned); the full suite **1017
passed, 3 xfailed, 30 subtests** in 101 s.

## 2026-08-26 (night, last +2) — the last four `.i` files with Python in them

`avmeandistance.i`, `cif.i`, `interactionterms.i` and `structureio.i` were the
four left. They are `%include` files now, and there is no `%pythoncode` left in
`pyext/` outside `swig.i-in`'s `import numpy as np`.

### What each one was, and where it went

**`cif.i`** (4 defs). `read_dye_forcefield_cif(path)` was
`read_forcefield_cif(str(path))` -- a second name for the C++ reader, spelled
in Python so it could take a `pathlib.Path`; the callers say `read_forcefield_cif`
now. `write_dye_forcefield_cif` *shadowed* the C++ writer of the same name to
coerce its arguments, and had to call `_IMP_bff.` to avoid recursing into
itself. `forcefield_system_from_dict` and `as_forcefield_system` turned a dict
into a typed system through `json.dumps`: that is what a caller does with a
dict it built, so `bin/imp_bff` and the three tests that write systems as dict
literals coerce their own.

**`interactionterms.i`** (2 `%feature("shadow")` constructors). Three reasons,
all typemaps or features:

| the shadow did | now |
|---|---|
| keyword arguments with C++ defaults | `%feature("compactdefaultargs")` -- SWIG's default-argument *overloads* are what disable `kwargs`, and compact ones do not |
| `dict(...)`, `[str(x) for x in ...]`, `np.ascontiguousarray(...)` | std_map, std_vector and numpy.i typemaps, which already did all three |
| `kappa2=None` for "no orientation factor" | a `%typemap(in) double kappa2` that reads `None` as NaN |

**`avmeandistance.i`** (93 lines). `AVNetworkRestraintWrapper` subclassed
`IMP.pmi.restraints.RestraintBase`, so it was built lazily -- naming it at
import time would have made `import IMP.bff` require IMP.pmi. What it did
besides PMI's bookkeeping is `av_network_restraint_set` in
`AVMeanDistanceRestraint.h`: the network restraint, or one
`AVMeanDistanceRestraint` per distance over volumes made rigid-body members,
plus `add_avs_to_rigid_bodies` and `set_av_xyzr_mass`. **IMP.pmi stays
optional**: a `RestraintBase` subclass has no C++ spelling at all, so the
fifteen lines that make one are in the two PMI examples that want one. Nothing
had ever tested the wrapper; the set has five tests.

**`structureio.i`** (359 lines). Four RMF functions and an `IMP.rotamer` query,
all behind the same lazy door. They are C++ in `RmfIO.h` and
`ProbeAttachment.h`, which cost a dependency -- see below.
`place_probe_from_rotamer_cb` is not ported: it computed the
backbone-dependent C-beta, logged it, and then called
`place_probe_from_coords(label, ca, n, c)`, which is what a caller gets by
calling that function. The C-beta was never used.

With those gone there is **no `_LAZY` table and no module `__getattr__`**: every
public name is an ordinary attribute again, and a missing one is a miss.

### The dependency, and what it cost

`dependencies.py` gained `rmf` and `rotamer`. `rotamer` is free (`core:atom`).
`rmf` is not: it brings `isd`, and `isd` brings `saxs`, and SWIG wraps a type
once across a module and its imports -- so `saxs`'s
`%template(DistBase) std::vector<double>` is declared first and `types.i`'s
`%template(VectorDouble)` is **skipped**. `VectorDouble`, `VectorString`,
`VectorFloat` and `VectorInt` are no longer attributes of `IMP.bff`;
`IMP.saxs.DistBase` *is* `std::vector<double>`.

Only code that *names* those classes is affected -- a list or a numpy array
converts as it always did -- and in this tree that was three test files, which
now say where the class comes from. The two functions that made it visible
(`wobbling_kappa2_distribution`, `diffusion_propagate`) take
`std::vector<double>&` output collectors, which is the pattern that needs a
wrapped vector in the first place.

Two other things the import needed: `#include <RMF/TraverseHelper.h>` in the
wrapper (IMP.rmf's conversion table names it; an importer has to include it
itself), and `IMP_SHOWABLE_INLINE` on `LJSitePair`, `Mol2Atom`, `AlignedBlock`,
`ComponentTemplate::FeatureAtom` and `ComponentTemplate::Improper` -- RMF's
`Showable` is what SWIG reaches for when printing a `std::vector<T>`, and it
needs `operator<<` on `T`. Those five value types had none; they should have
had one anyway.

**This is reversible.** Dropping `rmf` from `required_modules` puts the four
RMF functions back in Python and restores the four template names; it is one
line either way.

### Verification

`ninja IMP.bff-python` clean; the full suite **1008 passed, 3 xfailed, 30
subtests** (one timing test, `test_a_short_walk_on_a_large_grid...`, failed
once under CPU contention with a second build and passes alone). New tests:
`test/io/test_rmf_io.py` (8) and `TestAVNetworkRestraintSet` in
`test/test_AVNetworkRestraint.py` (4).

## 2026-08-26 (night, last +1) — `ProbeSite`, `ProbeRestraints`; labelizer stays

Finishing the flrCIF alignment. The three names left behind last time were
`labelizer`, `LabelingSite` and `LabelingRestraints`; the first is a program
and keeps its name, the other two are API and did not.

| was | is |
|---|---|
| `LabelingSite` / `LabelingSites` | `ProbeSite` / `ProbeSites` |
| `LabelingRestraints.h`, `.cpp` | `ProbeRestraints.h`, `.cpp` |
| `DirectLabelingRestraint` | `DirectProbeRestraint` |
| `IMP_bff.labelingrestraints.i` | `IMP_bff.proberestraints.i` |
| `%template(LabelingSiteVector)` | `%template(ProbeSiteVector)` |

`ProbeSite` is a residue, an attachment atom and a measured distance -- the
same object whether what sits there fluoresces, carries an unpaired electron
or quenches. `DirectProbeRestraint` scores attachment-atom distances with no
volume at all, which is as true of a spin label as of a dye.

Two headers, two sources, one `.i` renamed; a cmake reconfigure before SWIG,
because a renamed header is a new header as far as the dependency list is
concerned.

Verification: `ninja IMP.bff-python` clean; **1004 passed, 3 xfailed, 30
subtests**, `medium_test_restraints.py` 18 passed.

## 2026-08-26 (night, last) — flrCIF says *probe*, so this package does

Asked which term aligns with the dictionary, `label` or `probe`. Checked
`mmcif_ihm_flr_ext.dic` rather than answering from memory:

- **probe appears 86 times, label 5.** The categories are `_flr_probe_list`,
  `_flr_probe_descriptor`, `_flr_poly_probe_position`,
  `_flr_poly_probe_conjugate`, `_flr_sample_probe_details`. There is no label
  category at all.
- **The five `label` hits mean something else.** `mm_atom_site_label` and
  `mm_poly_res_label` are mmCIF's `label_asym_id` / `label_seq_id` /
  `label_comp_id` convention -- the canonical-identifier naming. The other two
  are prose in an EPR example ("Spin label rotamer refinement using DEER").

The second point decided it: the value's fields **are** the `label_*` items,
so naming the type `Label` would have named it after its own field prefix, in
a dictionary where `label_` already means something.

| was | is | flrCIF |
|---|---|---|
| `Label` | `ProbePosition` | `_flr_poly_probe_position` |
| `label_from_source_info` | `probe_position_from_source_info` | |
| `label_flrcif_items` | `probe_position_flrcif_items` | |
| `LabelAttachment`, `attach_labels` | `ProbeAttachment`, `attach_probes` | |
| `place_label`, `resolve_label_site` | `place_probe`, `resolve_probe_site` | |
| `AttachedLabelDynamics` | `AttachedProbeDynamics` | |
| `label_forcefield_system` | `probe_forcefield_system` | |
| `Labeling.h` | `ProbeAttachment.h` | |

flrCIF also splits what this value merges -- the probe itself
(`_flr_probe_list`), where it attaches (`_flr_poly_probe_position`) and the
pairing that carries `fluorophore_type` (`_flr_sample_probe_details`). The
struct is the position, naming its probe; the split is recorded in the header
for whoever needs the other two categories.

"Labelling site" stays in prose -- it is what the field calls it, and
`labelizer`, `LabelingSite` and `LabelingRestraints` keep their names. What
changed is the identifiers a caller writes.

Verification: `ninja IMP.bff-python` clean; **1004 passed, 3 xfailed, 30
subtests**.

## 2026-08-26 (night, later) — a label is not a dye

The labelling layer described a *fluorophore* everywhere -- `resolve_dye_site`,
`place_dye`, `attach_dyes`, `DyeAttachment`, `LangevinDyeSampler`,
`dye_forcefield_system` -- while the package already ships
`spinlabels.drot.pto` beside `dyes.drot.pto` and the rotamer machinery
(`RotamerLibrary`, `RotamerEnsemble`) was general all along. The attachment
half was the only place that insisted.

**From a structure's point of view there is no difference.** A fluorophore at
a cysteine, a spin label for DEER and a quencher tryptophan are one object: a
molecule on a linker, anchored at a residue's backbone frame, with a rotamer
library behind it. What differs is what is *measured* -- spectra for FRET, a
dipole for DEER -- and that lives with the measurement, not with the site.

So the structural vocabulary is now labels:

| was | is |
|---|---|
| `resolve_dye_site` | `resolve_label_site` |
| `place_dye`, `place_dye_from_coords` | `place_label`, `place_label_from_coords` |
| `attach_dyes`, `DyeAttachment` | `attach_labels`, `LabelAttachment` |
| `LangevinDyeSampler` (`DyeDynamics.h`) | `AttachedLabelDynamics` (`LabelDynamics.h`) |
| `dye_forcefield_system` | `label_forcefield_system` |
| `dye_internal_system` | `internal_topology_system` |
| `DyeInternalEnergyEvaluator` | `IntramolecularEnergy` |

`Label` carries a `probe` **name** and, when the probe happens to be a
fluorophore this package knows, a `Dye` value with the photophysics.
`set_dye` sets both, so the two cannot disagree about what is attached.

**What was deliberately not added.** A first draft gave `Label` a
`probe_type` enumeration (`fluorophore` / `spin_label` / `quencher`) and
refused a FRET role to a non-fluorophore. That is over-modelling at the
structural layer, which cannot tell the three apart and should not try: the
field is gone and a probe with no spectra simply has no `dye`.

`Dye`, `find_dye`, `forster_radius_from_spectra` and the FRET distance
conventions keep their names. A dye *is* a dye where photophysics is meant.

### Two things the rename caught

1. **A name that had been retired came back.** `InternalEnergyEvaluator` is on
   `test_public_api_names`'s retired list, and renaming
   `DyeInternalEnergyEvaluator` to it un-retired something deliberately
   buried. It is `IntramolecularEnergy`, which says what it computes -- the
   intramolecular non-bonded energy of a conformation -- and was never
   anything else.
2. **A component name split in two.** The attached molecule's force-field
   component is named once and its anchor group is that name plus `_anchor`;
   renaming the first and not the second left every anchor unfound, so
   nothing was held fixed. Both come from one constant now, and the
   thermodynamic test caught it (`len(s.fixed) == 4` became `0 == 4`).

Verification: `ninja IMP.bff-python` clean; **1004 passed, 3 xfailed, 30
subtests** -- the new one asserting that a spin label is a label:
`Label("A", 41, "SG", "MTSSL")` has a site, a key and no `dye`.

## 2026-08-26 (night) — `label.i` is empty, and the backbone frame existed twice

986 lines gone. `Labeling.h` carries the frame, the site, the strip, the
attachment and the `Label` value; the fluorescent-protein detection, the pLDDT
parsing and the segmentation went to `SequenceAlignment.h`, beside the
Smith-Waterman they are built on.

**The find: the backbone frame existed twice.** Once as `backbone_rotation`
(`RotamerSite.h`), which every rotamer library is placed with, and once in the
Python labelling layer, which every explicit dye was placed with. The same
three cross products -- x along CA->N, z out of the peptide plane, y
completing the set -- in two languages, for the two halves of one package. One
`backbone_frame` now; the flat-matrix spelling is a view of it.

Three smaller ones:

- **`atom_name_of` was written out in three files.** The unity build caught it
  as a redefinition. It is `atom_name` in `HierarchyFrame.h`.
- **`strip_hierarchy` defaulted to `inplace=False`**, which cloned an entire
  structure on every call without saying so. Stripping is in place; a caller
  that wants the original clones it and can see the cost.
- **`attach_dyes` took three parallel vectors** in the first draft -- the exact
  smell being removed elsewhere. It takes `DyeAttachment` values, which also
  lets it report where each dye was placed and how many atoms were stripped.

`Label` now carries a `Dye` **value** rather than a dye name: a label that
carries only a name cannot derive an R0 or a correlation time without going
back to the library for what it already had.

### What could not be C++, and why

`get_anchor_cb_position` and `place_dye_from_rotamer_cb` need **`IMP.rotamer`**
-- a third optional IMP module, alongside `IMP.pmi` and `IMP.rmf`, that this
module does not depend on. They joined the lazy island in `structureio.i`,
which is why that file grew from 189 to 359 lines while the total fell.

### API changes that moved tests

A frame is a `Transformation3D`, not a `ReferenceFrame3D`; `resolve_dye_site`
answers CA, N, C **in order** rather than a dict; `select_atoms` returns
particles; `strip_obstacles` is `strip_keep_mask`, which says which rows
survive and leaves the slicing to numpy; `SITE_KEEP_ATOM_NAMES` is
`site_keep_atom_names()`. `Label` is constructed positionally, because an IMP
value needs a default constructor too and SWIG turns keyword arguments off for
anything overloaded.

The `_find_atom` test was **rebuilt, not deleted**: it pins *which* atom a
lookup returns -- first in hierarchy order, right chain and residue, matched
on the name as the structure spells it -- and it now checks `select_atoms`
against the same exhaustive scan.

### The same build trap, twice

A header added since the last cmake run is not in SWIG's dependency list, so
`wrap.cpp` is generated from the previous declaration and the wrapper compiles
against an API that no longer exists. It cost a cycle for `RRT.h` and another
for `Labeling.h`, with the lesson written down in between. **Reconfigure after
adding a header**, not after the error.

`label.i` 986 -> 0; the `%pythoncode` total is **509** (11,607 when PRD-117
started). What is left is almost entirely the optional-module islands:
`structureio.i` 359 (RMF, rotamer), `avmeandistance.i` 93 (the PMI wrapper),
`cif.i` 28, `swig.i-in` 24, `scoring.i` 5.

Verification: `ninja IMP.bff-python` clean; **1003 passed, 3 xfailed, 30
subtests**.

## 2026-08-26 (later) — `sampling.i` is empty

1,041 lines gone, in three headers and one addition to a fourth:

- **`DyeDynamics.h`** -- `make_langevin_simulator` (md gets a Langevin
  thermostat, bd gets Einstein's coefficient for each particle's own radius),
  `LangevinDyeSampler`, `LangevinTrajectory`.
- **`LinkerSampling.h`** -- the Metropolis sampler over a linker's torsions
  and bond angles, and the library it clusters into. **No IMP model at all**:
  the Python wrote every configuration into `XYZ` decorators and read the
  coordinates straight back, using the particles as a scratch buffer for
  numbers `LinkerGeometry::apply` had just returned.
- **`RRT.h`** -- the two planners. Kept as two functions rather than one
  generic tree because their metrics differ *in kind*: the rigid one adds a
  translation distance to a weighted rotation distance, which no single norm
  over six numbers expresses. Ångström and radian are not the same thing.
- **`markov_state_trajectory`** in `Clustering.h`, beside the transition
  matrix it draws from.

A caller's collision test crosses back through the `RRTCollision` director,
as a docking run's cancellation does through `DockingStop`.

### What the port changed, and what it cost

**The sampler pins moved, and had to.** `generate_linker_rotamers` is a
Metropolis walk; C++ draws from `boost::mt19937` where Python drew from its
own `random`. Same process, different stream, so every number downstream is a
different draw. The pins are re-recorded and their note says so -- drift
guard, not reference values. The invariants hold: 30 clusters, weights summing
to one, 29 transitions across 30 frames.

**An error the port caught in my own first draft**: I weighted each cluster by
its *representative's* Boltzmann factor. The original sums over members, and
that is right -- a broad shallow basin holds more of the ensemble than a
narrow deep one.

**Three tests are gone and their coverage is not.** The Python had unit tests
for `steer_transform`, `is_collision_sphere`, `rrt_grow_step` and
`tree.nearest`, all internal now. The tree-level assertions say the same
thing: no step longer than `step_size`, no node in a collision, a reachable
goal reached -- plus a case the old suite did not have, that an unreachable
goal reports `goal_node == -1` rather than failing, because a tree that did
not reach still explored. Suite 1006 -> 1003 for exactly that reason.

**A fourth copy of Boltzmann's constant** was in `test_physics_invariants`;
both tests take `kb_kcal()` now.

**A second van der Waals table** was in the dye sampler, keyed by symbol,
disagreeing with `AVBuilder`'s (by atomic number) about hydrogen -- 1.10
against 1.20, Rowland & Taylor against Bondi. One `vdw_radius(element)` now;
the thermodynamic pins did not move.

### Two build lessons, both paid for

1. **A shared build directory makes exit 0 meaningless.** `wrap.cpp`
   regenerated at 18:28 while `_IMP_bff.so` stayed at 18:13: this session's
   ninja planned its work, a peer's ninja regenerated the wrapper underneath
   it, and mine exited 0 believing it was done. Check the artefact
   timestamps, not the exit code.
2. **A new header needs a cmake reconfigure before SWIG sees it.** The swig
   rule's dependency list is fixed at configure time, so `RRT.h` changed and
   SWIG never re-ran -- the wrapper compiled against the previous
   declaration. `cmake -S ../imp -B .` then ninja.

Also: `%attribute_py(Class, T, name, name)` on a *public member* makes a
property of a property (`'property' object is not callable`). SWIG already
exposes a public member; `%attribute` with a real getter is for methods.

`sampling.i` 1,041 -> 0; the `%pythoncode` total is **1,325** (11,607 when
PRD-117 started). Left: `label.i` 986, `structureio.i` 189 (the RMF island),
`avmeandistance.i` 93, `cif.i` 28, `swig.i-in` 24, `scoring.i` 5.

Verification: `ninja IMP.bff-python` clean; **1003 passed, 3 xfailed, 30
subtests**; the 36 expensive all-dyes sampling tests green.

## 2026-08-26 (later) — `sim.i` is empty, and there is one repulsion

**The runner moved to `bin/imp_bff`.** All 1,319 remaining lines of it: it
reads a directory of force-field systems, propagates each one and writes
trajectories, and it was reached through a three-line `simulate` command that
forwarded twenty-five options into the library. It needed exactly eight names
from `IMP.bff` (one explicit import) and collided with nothing already in the
program.

What it *computes* stayed in C++ -- `build_dye_restraints`,
`build_steric_restraint`, `build_go_restraints`, `place_guest_by_score` --
because those are kernels rather than orchestration, and because the runner
had grown its own copy of the restraint builder that disagreed with the C++
one.

### One repulsion, for dynamics and for Monte Carlo

There were two implementations of the same physics: per-pair Lennard-Jones
lower bounds in `build_dye_restraints`, and soft spheres over a container in
the Monte-Carlo path. It is `build_steric_restraint` everywhere now --
differentiable, so a dynamics run can use it; one restraint over a container
rather than one per pair, so a long run can afford it; and dependent only on
where the spheres are, so a rigid move can be scored against it *alone*, which
is the whole reason the Monte-Carlo path wanted its own term. A rigid move
cannot change a bond, an angle or a torsion.

The Lennard-Jones parameters remain what `DyeInternalEnergyEvaluator`
evaluates. An energy of a conformation is a different question from keeping
two atoms apart.

**A hazard closed on the way**: a soft sphere is a sphere, so a site the
caller decorated without `XYZR` would have contributed nothing at all, and
silently -- the run would look fine and the molecule would pass through
itself. The builder gives such a site the radius the system says it has.

### The suite did not cover the term that changed

The Langevin test's `d.min() > 1.5` is a *dye-protein* assertion and comes
from a separate bipartite restraint, so it says nothing about the intra-dye
repulsion. Measured directly instead -- the steric term is the last of 697
restraints on the bundled dye system, scores 294 on an artificial line
layout, and rises when two sites are superposed -- and then pinned, because
"the suite is green" is not evidence when the suite does not look:
`test_the_repulsion_is_one_term_shared_by_dynamics_and_monte_carlo` (one
restraint, not thousands; absent when asked for), that it rises on overlap,
and that a site without a radius still repels.

`sim.i` 1,319 -> 0; the `%pythoncode` total is **2,366** (from 11,607 when
PRD-117 started). Left: `sampling.i` 1041, `label.i` 986, `structureio.i` 189
(the RMF island), `avmeandistance.i` 93, `cif.i` 28, `swig.i-in` 24,
`scoring.i` 5.

Verification: `ninja IMP.bff` clean; **1006 passed, 3 xfailed, 30 subtests**,
the three `imp_bff simulate` modes among them.

## 2026-08-26 — `sim.i`: three kernels out, and a restraint builder that existed twice

The cgdye MD runner, started. `bin/imp_bff simulate` was a three-line shim
forwarding twenty-five options into a 1,491-line body that lived in the
library -- the same shape docking had.

**The restraint builder existed twice.** `sim.i::_build_restraints` and the
C++ `build_dye_restraints` both walk a force-field system into IMP restraints,
and they disagreed in two ways:

- *The equilibrium fallback.* The Python restrained a bond about the geometry
  as it stands when the system carried no length; the C++ used the stored
  zero, which would pull the two sites together. A system built without a
  template carries exactly that. The C++ falls back the same way now, for
  bonds and for angles.
- *The steric term.* The C++ made one lower-bound harmonic **per pair** --
  thousands of restraints; the Python made one `PairsRestraint` over a
  container. The container form is what a molecular-dynamics run wants, and it
  is also what a *rigid-body* Monte-Carlo step should be scored against: a
  rigid move cannot change a bond, an angle or a torsion, so scoring those
  during one is work whose answer never changes.

So `build_dye_restraints` gained a `nonbonded` switch and `build_steric_restraint`
sits beside it. Both take the exclusions from
`DyeForceFieldSystem::get_exclusions`, so they cannot disagree about which
pairs are 1-2, 1-3 or 1-4.

Two more kernels moved:

- **`build_go_restraints`** -- native contacts holding a component in its own
  shape. The mobile component is restrained everywhere and the fixed one only
  where `fixed_flex_mode` released atoms, which is an `only_sites` filter
  rather than two near-identical loops. The `max(d, 1 A)` equilibrium is
  written down for what it is: two atoms a structure has placed on top of each
  other would otherwise get a harmonic at zero, and a minimiser walks straight
  into it.
- **`place_guest_by_score`** -- the rigid random placement search that gives a
  simulation somewhere plausible to start, rather than wherever two separate
  input files happened to put a dye and a protein (usually inside each other).
  Its uniform-sphere direction is now written down as the one construction
  that does not crowd the poles.

`_center` is gone in favour of `IMP.core.get_centroid`: a local mean nobody
outside that file could find, when IMP has had the function all along.

Unlike docking, this runner **is** covered -- `test_cgdye_integration` runs
`imp_bff simulate` end to end in three modes -- so the C++ builders are
exercised against real MD rather than only by construction.

`sim.i` 1,491 -> 1,323 lines; the `%pythoncode` total is **3,689**.

Verification: `ninja IMP.bff` clean; **1003 passed, 3 xfailed, 30 subtests**
after each of the three ports.

## 2026-08-25 (night, last) — `docking.i` is empty

1,427 lines of `%pythoncode` when the file was first opened this session;
**zero** now. It is declarations, `%feature("kwargs")` and attributes.

Where each kind of thing went, by what it *is*:

- **An API** -- `score_structures`, `dock_minimize`, `refine_docking`,
  `screen_structures`, the assembly, the pair table, the poses, the CSV -- is
  C++ in `Docking.h`. `refine` and `screen` were renamed on the way in:
  `IMP::bff::refine` says nothing in a namespace this size.
- **Programs** -- `imp_bff dock` (PMI replica exchange) and `imp_bff
  dock-errors` (forked trial driver) -- are in `bin/imp_bff`. Both write
  directories of results; one needs a Python-only dependency and the other
  forks workers.
- **The cancellation callback** is #IMP::bff::DockingStop, an `IMP::Object`
  with a SWIG director, rather than a bare Python callable. A caller
  subclasses it in either language and the C++ loop asks it between chunks --
  which is also the only reliable place to ask, since raising out of an
  `IMP::OptimizerState` mid-step is not.

Two things were **removed rather than translated**, both because they never
did anything:

- `IMP.pmi.tools.shuffle_configuration` was called inside a bare
  `except: pass`, so the initial shuffle silently did nothing whenever PMI was
  missing. It is a plain random displacement of the mobile bodies now, with no
  dependency.
- `ensure_fps_json` (the one with the undefined `_io`) is gone: `read_fps_json`
  already reads the legacy C# format, so `build_docking_assembly` converts a
  non-`.json` input once, in the one place that needs it.

### And one I introduced

Moving `estimate_errors` into `bin/imp_bff` carried its annotation
`params: Optional[DockingParameters] = None` with it. Annotations are
evaluated at `def` time and that name does not exist in a program file, so
**the whole `imp_bff` program stopped importing** -- 3 failures and 4 errors,
every test that shells out to it. Caught by the suite, fixed by dropping the
annotation. Worth remembering when moving code between a `.i` and `bin/`: the
library's names are not in scope there, and an annotation is not free.

`docking.i` 831 -> 0; the `%pythoncode` total is **3,857**. Left: `sim.i`
1491, `sampling.i` 1041, `label.i` 986, `structureio.i` 189 (the RMF island),
`avmeandistance.i` 93, `cif.i` 28, `swig.i-in` 24, `scoring.i` 5.

Verification: `ninja IMP.bff` clean; **1003 passed, 3 xfailed, 30 subtests**;
the twelve docking tests green against the C++ entry points.

## 2026-08-25 (night, later) — `dock` is a command; the rest of the engine is not a program

Asked whether the programs had been relocated: for docking they had not. The
samplers were still `%pythoncode` and this entry corrects that -- but only for
the one piece that is actually a program.

**`imp_bff dock`** (`bin/imp_bff`) is the PMI replica-exchange sampler. It
drives `IMP.pmi.macros.ReplicaExchange` -- IMP.pmi is Python-only and not one
of this module's `required_modules`, so it can never be C++ -- and it writes an
RMF trajectory, PMI stat files and best-scoring PDBs into a directory. That is
a program by any reading. Verified end to end against the bundled T4L fixture.

**`score`, `dock_minimize`, `refine` and `screen` stay in the library**, and
the reason is the consumer rather than taste: ChiSurf's FRET plugin calls
`dock_minimize(...)` and `estimate_errors(...)` **with a `stop_check`
callback from its GUI thread** and reads the returned `DockingResult`. Shelling
out to a command would lose both the cancellation and the structured result.
They are an API, so the rule that applies is "what can be C++ must be C++",
which is where their assembly, pair table, poses and CSV went earlier today.

`estimate_errors(method="mc")` now raises, naming `imp_bff dock`: its
Monte-Carlo branch called the function that moved. ChiSurf offers both
`minimize` and `mc`, so this is a real (if currently theoretical) capability
loss on that side.

### The consumer is already broken, and has been for a while

`chisurf/plugins/modelling/fret/core/imp_engine.py` is a forwarder to
`IMP.bff.fret.imp_engine` -- a sub-package the flat-namespace commit deleted.
It raises `ModuleNotFoundError` at import, independently of anything done
here. Both sides of this API were dead: the engine could not run (four
`NameError`s, one per entry point) and the caller could not import it. That is
why nothing noticed either.

Whoever picks the ChiSurf side up needs to know that `IMP.bff.dock` is gone
from the library and that the flat names (`IMP.bff.score`,
`IMP.bff.dock_minimize`, ...) are what the forwarder should point at.

`docking.i` 899 -> 831 lines.

Verification: `imp_bff dock --help` and a two-frame run on the T4L fixture;
**1003 passed, 3 xfailed, 30 subtests**.

## 2026-08-25 (late night) — the docking engine runs, and every entry point had to be walked

The `NameError` in `build_assembly` was the first of four. Each entry point
fails on its first unexercised line, so the only way through was to run them
one at a time:

1. `build_assembly` looked up the lazy `AVNetworkRestraintWrapper` as a bare
   global -- a module `__getattr__` answers attribute lookups only. Every
   entry point died here.
2. `ensure_fps_json` called `_io.read_fps_json(...)`; **`_io` is defined
   nowhere in the module**, so every legacy C#-format conversion raised
   `NameError`. Only a non-`.json` input reaches it, which is why the first
   fix did not surface it.
3. `dock_minimize` (and `dock`) indexed `asm.rigid_bodies` as a dict.
4. `screen` catches per-structure exceptions and records `NaN`, so a broken
   engine produced a full ranking table of NaNs and looked like a screening
   run with unscorable inputs.

All five now run against the bundled T4L fixture, and
`test/test_docking_values.py` pins them: `score` (22.69935116408601 -- the
same number `test_AVNetworkRestraint.py` pins for the quadrature score, so the
docking path and the restraint agree), `dock_minimize`, `refine`, `screen`,
and the pose round-trip. Twelve tests, four seconds.

### The assembly is C++ and needs no PMI

`build_docking_assembly` (`Docking.h`) does what the PMI wrapper did --
resample each volume, attach it to the rigid body of the atom it hangs off,
give it a radius and a mass, add the mean-distance restraints and the
excluded-volume term -- in plain IMP. The PMI wrapper was the only reason the
*scoring* path pulled a sampler's dependency, and the lazy name it was reached
through is what made the engine unrunnable. `dock` still drives PMI's replica
exchange, and registers the assembly's restraint set with one call.

`capture_poses` / `apply_poses` are C++ too, JSON in and out, so a docked
state can be stored and continued from without Python.

### What is still not tested, and should be said plainly

**The bundled fixture is one rigid body.** No rigid motion of a single body
changes the distance between two of its own dyes, so there is nothing for a
docking optimiser to optimise, and the score moving at all during
`dock_minimize` (16.02 -> 15.48) is the AV resampling between builds and not
the minimiser. Nothing in this tree exercises multi-body docking. That needs a
two-body fixture; the tests say so rather than implying coverage they do not
have.

`estimate_errors` and its forked-worker trial runner are still unexercised.

`docking.i` 1,119 -> 899 lines; the `%pythoncode` total is **4,756**.

Verification: `ninja IMP.bff` clean; **1003 passed, 3 xfailed, 30 subtests**.

## 2026-08-25 (night) — the docking engine could not run, and nothing said so

`IMP.bff.score`, `dock`, `dock_minimize`, `screen` and `estimate_errors` have
raised `NameError` on their first line since the PMI wrapper went lazy.
`build_assembly` looks up `AVNetworkRestraintWrapper` as a bare global; the
wrapper is provided by the module's `__getattr__` (it needs `IMP.pmi`, which
is not a dependency of this module), and **a module `__getattr__` answers
attribute lookups only** -- never a bare global inside a function defined in
that module. One `IMP.bff.` prefix fixes it.

Nothing noticed because the engine had no test of any kind, in this tree or
before the move: 1,400 lines driving IMP.pmi over a whole assembly, and no
caller here builds one. `test/test_docking_values.py` is its first coverage,
and the end-to-end case is worth more than it looks: `score()` on the bundled
T4L fixture returns **22.69935116408601**, which is exactly what
`test_AVNetworkRestraint.py` pins for the quadrature score of the same
fixture. Assembling a model and scoring it through the docking path gives what
the restraint gives on its own, and both numbers now move together or not at
all.

### What moved

`Docking.h` carries what a run is *told* and what it *reports* --
`DockingParameters`, `PairDistance`, `DockingResult` -- with the two readers
of a scored assembly beside them: `collect_pair_distances` (from the
network's own volumes) and `pair_distances_at_positions` (from the proxy
particles the minimisation path moves instead), plus `write_score_csv`.

**The mean-distance restraint defined in `docking.i` was a third copy.**
`AVMeanDistanceRestraint` has had the gradient in C++ since the AV batch --
its header records that the Python had two copies, "one with the gradient and
one without, and the one without was the one the wrapper used". This was the
third, and the only one still in Python. It is gone; the C++ class takes the
Förster radius and the distance type from the measurement rather than from a
re-parse of the JSON the measurement was built from.

**The fps.json distance vocabulary was private to a reader.**
`DyePairMeasure_name_to_type` sat in `internal/FPSReaderWriter.h`, so a
diagnostics table that wanted to report which convention it had scored
re-derived the names from the JSON -- and said `RDAMean` for an entry whose
absent `distance_type` the reader had defaulted to `RDAMeanE`. It is
`dye_pair_distance_type()` / `dye_pair_distance_type_name()` in `AV.h` now,
one table read both ways, used by the reader and by the table.

`docking.i` 1,427 -> 1,119 lines; the `%pythoncode` total is **4,976**. What
remains in it is the engine -- build the assembly, run the Monte-Carlo or the
minimisation, write the RMF and the PDBs -- which is IMP.pmi orchestration.
Now that it runs, it can be moved with a test watching.

Verification: `ninja IMP.bff` clean; **1000 passed, 3 xfailed, 30 subtests**
(999 before this file's own test).

## 2026-08-25 (evening) — `rotamer.i` is empty: the loaders, the fps layer and the FRETpredict driver

1,277 lines of `%pythoncode` when this batch started; **zero** now, a file
comment and one `%include`. Three new headers carry what left it.

**`RotamerFps.h` — fps.json and rotamer ensembles.** `RotamerPosition`,
`RotamerDistance` and `RotamerFpsSelection` (what `read_rotamer_fps` returned
as a six-tuple), the payload builders, `distances_from_ensembles`,
`write_rotamer_fps` and `rotamer_ensembles_from_fps`. Entries cross as JSON
text, which is what the rest of the fps layer does: an entry carries whatever
keys its writer put there, and a typed struct would either lose them or grow a
field per program. What *is* typed is the part this module reasons about --
which chain, which residue, which library, which dye.

The alias reading is the format's business and is done once now: a position's
chain is `chain_identifier`, `chain` or `segid`; a distance's donor is
`position1_name`, `donor_position`, `donor_position_name` or `dye1_position`.
`rotamer_ensembles_from_fps` had grown its *own* copy of that list -- the unity
build caught the duplicate as a redefinition -- and reads through the typed
parser now, one library load per library rather than one per position.

**`RotamerFret.h` — the FRETpredict driver.** Twenty options, five arrays and
the file writing, over kernels that were C++ already. The parameter names stay
FRETpredict's (`fixed_R0`, `ign_H`, `libname_1`, `r0lib`), deliberately: the
parity harness hands *one* keyword dictionary to this class and to FRETpredict
and compares the files they write, so a rename here would be a rename of the
experiment. `distance_distributions` was **not** ported: `trajectory_analysis`
allocated it whenever `calc_distr` was set and nothing ever wrote to it, so the
option and the `rmin`/`rmax`/`dr` axis it was shaped from are gone.

**The registry shims are gone too.** `rotamer_library_registry`,
`normalize_library_name`, `rotamer_library_metadata` and
`resolve_rotamer_library_path` turned the C++ answers back into dicts and
`Path`s. The C++ ones answer directly, and the registry itself is published as
JSON text (`rotamer_library_registry()`, read once and cached).

One more name of the kind the last two batches were about:
`rotamer_frame_weights` -> **`frame_weights_from_partitions`**. It takes two
partition functions per frame and weights the frames by their product; nothing
in it is particular to a rotamer.

`OwnedView` -- a malloc'ed out-view freed on scope exit, which is what a *C++*
caller of the numpy out-view protocol needs -- moved into
`internal/OutputView.h`, beside the protocol it belongs to.

### Two defects found, one fixed here

**The expensive FRETpredict parity harness has been feeding R0 = 5.5 to both
sides.** R0 became Ångström across this package on 2026-08-25; the pins file
was restated that day (5.5 -> 55.0) and
`expensive_test_fretpredict_parity.py` was missed, so it handed 5.5 Å to
IMP and 5.5 nm to FRETpredict out of one shared dictionary. It skips unless
FRETpredict is installed, which is why nothing said so. Fixed: `_imp_kwargs`
converts the one value whose unit differs, and says why.

**`RotamerFRET(rmf_path)` no longer works**, and that is a loss, not a fix:
`IMP.rmf` is not one of this module's `required_modules` and C++ cannot open
an RMF here. The trajectory goes through the lazy door instead --
`RotamerFRET.from_frames(protein_frames_from_rmf(path), ...)` -- which is the
same `ProteinFrame` values the C++ PDB reader returns. `from_frames` is a
named factory and not a second constructor because SWIG turns keyword
arguments off for anything overloaded, and this constructor is keyword
arguments in every caller.

### Standing

`rotamer.i` 1,277 -> 0 and `rotamer_ensemble.i` 332 -> 0. The `%pythoncode`
total is **5,284** (was 6,166 this morning, 7,695 when PRD-117 started). Left:
`sim.i` 1491, `docking.i` 1427, `sampling.i` 1041, `label.i` 986,
`structureio.i` 189 (the RMF island), `avmeandistance.i` 93, `cif.i` 28,
`swig.i-in` 24, `scoring.i` 5.

Verification: `ninja IMP.bff` clean; **991 passed, 3 xfailed, 30 subtests**
(988 before -- the three new are the API test's coverage of `RotamerFRET`,
`load_rotamer_library` and the lazy `protein_frames_from_rmf`); the
FRETpredict pins pass through the C++ driver unchanged; both
`examples/structure/` dye examples run end to end and the cgdye notebook
executes.

## 2026-08-25 (later) — `rotamer_ensemble.i` is empty, and three structs that were one library

**`rotamer_ensemble.i` holds no Python.** 332 lines of `%pythoncode` defining a
Python class that *subclassed a C++ value*; zero now. `RotamerEnsemble` is a
C++ `States` (`include/RotamerEnsemble.h`) carrying the atoms, the energies,
the partition function and the site, with `from_frame` / `from_site` placing a
library on a residue and `pair_geometry` / `pair_distribution` /
`pair_distribution_from_dyes` answering with the typed pair values instead of
rebuilding dicts out of them. It reproduces the FRETpredict parity pins
exactly -- E_static, E_dyn1, E_dyn2, ⟨κ²⟩, R0 and both partition functions,
`test/cgdye/rotamer/test_rotamer_ensemble.py`, all eight tests green.

**One rotamer library value, not three.** `RotamerLibraryData` (`CifIO.h`, the
numpy/text pair), `RotamerLibrary` (`DyeSampling.h`, a PDB plus a trajectory)
and `DrotLibrary` (`DrotReader.h`, a `.drot` container) held the same thing --
conformers, a weight each, the atom names the coordinates are ordered by -- and
differed only in which optional columns their own reader filled. A caller that
took one could not be handed another, so each reader grew its own consumers,
and `id` (always `1..n`) existed in exactly one of them. There is one
`RotamerLibrary` now (`include/RotamerLibrary.h`), with `resnames`, `elements`,
`metadata` and `path`, and three readers returning it.

Two names went with it, both of the kind the last batch was about:
`load_rotamer_library_dcd` reads a PDB plus *any* trajectory, so it is
`load_rotamer_library_trajectory`; and `DrotLibrary` named the container it
first came out of.

**The loaders are C++.** `load_rotamer_library(name, lib_dir)` (registry name,
locator or path -> the library with its metadata, resnames filled from the
container, from the `<stem>.pdb` beside it, or inferred) and
`load_protein_frames(path)` (a multi-MODEL PDB -> `ProteinFrame` values) were
the last Python between a file and an ensemble. `ProteinFrame` sits in
`HierarchyFrame.h`, beside the two functions that read one out of a hierarchy.
`rotamer_ensembles_from_fps` is C++ too, and reads each library once rather
than once per position.

**RMF stays out of the dependency graph.** `IMP.rmf` is not one of this
module's `required_modules`, so the RMF trajectory door is one more lazy
builder in `structureio.i` (`protein_frames_from_rmf`), six lines that open a
file and step it, returning the same `ProteinFrame` values the C++ PDB reader
returns.

**Shapes that were the caller's problem.** `FRETPairGeometry` and
`FRETPairEfficiencies` publish `R`, `kappa2`, `weight`, `E`, `rate_ratio` and
`k_fret` as `(n1, n2)` attributes now (`%attribute_np2v`, new); a library
answers `(n_rotamers, n_atoms, 3)` and a frame `(n, 3)`. Every consumer used to
fold those by hand out of `n1` and `n2` -- six sites in the ensemble code alone
-- and a matrix reshaped wrongly is a matrix transposed in silence.

`rotamer.i` 1277 -> 882 lines of Python; the `%pythoncode` total is **6,166**
(was 6,864). Left: `sim.i` 1491, `docking.i` 1427, `sampling.i` 1041,
`label.i` 986, `rotamer.i` 882, `structureio.i` 189 (the RMF island),
`avmeandistance.i` 93, `cif.i` 28, `swig.i-in` 24, `scoring.i` 5.

Verification: `ninja IMP.bff` clean; **988 passed, 3 xfailed, 30 subtests**;
both `examples/structure/` dye examples run end to end; the cgdye notebook
executes (it had been broken since the flat-namespace commit -- it imported
`IMP.bff.cgdye.utils`, called `kappa2_from_dipoles` and passed `R0=` to
`fret_efficiency_regimes`, all three gone; fixed, and its `.drot` paragraph
now describes the PTO container rather than the brotli+tar one).

**Found, not fixed** (not this session's model to re-pin):
`test/cgdye/rotamer/medium_test_av_vs_rotamer.py` fails on `av_n_points`
(2047 against a pinned 1137). The pins are from 2026-08-17 and the AV model
changed on 2026-08-21 (`9446fee`, "AV3 and the reference stencil"); the file is
named `medium_test_*`, which pytest does not collect, so nothing said so. Its
`bin/imp_bff` helper also still subscripted the typed pair values
(`g["weight"]`) -- that part is fixed here.

## 2026-08-25 (later still) — `topology.i` is empty, and eight names that described their first caller

**`topology.i` holds no Python.** 989 lines when PRD-117 started, 47 now, all
of it a file comment and a `%include`. The last of it: `build_forcefield_system`
turned dicts into the JSON string the C++ builder took, so the builder takes
`FFComponentSpec` values instead; `build_dye_protein_system` and
`dye_forcefield_system` are two- and one-component calls of it, in
`TopologyBuild.h` with the dye's `N`/`CA`/`C`/`O` anchor logic beside the
builder it belongs to; and `build_system_from_specs` was the `build-system`
command's body -- build, write, report -- which is a program's job and now sits
in `bin/imp_bff` with the `name=X,mol2=Y` parsing that is its own CLI syntax.
Five more helpers were dead: nothing called `sid`, `_alpha_suffix`,
`_serial_to_site_atom_names`, `_center_atom_serials_from_template` or
`_resolve_feature_ids`, and the site-id deduplication they duplicated is
`serial_to_site_atom_names` in C++.

**Then a naming pass, prompted by `write_path_map`.** The observation was that
a PINN will want to write a field to a map too, and the name says *path map*
when the operation is "write one voxel feature of a lattice to a density file".
Looking for the same defect -- a name that describes its first caller rather
than what it does -- found seven more:

| was | is | why |
|---|---|---|
| `write_path_map` | `write_map_feature` | picks a scalar field, writes a density file |
| `rotamer_transition_matrix` | `transition_matrix_from_counts` | counts to jump probabilities; Markov, not rotamer |
| `rotamer_correlation_times` | `relaxation_times` | eigenvalues of a transition matrix |
| `rotamer_rotational_correlation_time` | `slowest_relaxation_time` | it is *rotational* only when the states are orientations, which the arithmetic never checks |
| `rotamer_cluster_weights` | `cluster_weights` | aggregates frame weights by cluster label |
| `sample_rotamer_index` | `sample_weighted_index` | a weighted draw over an index |
| `apply_rotamer_coordinates` | `apply_coordinates` | puts a coordinate set on a hierarchy |
| `rotamer_pair_energy_matrix` | `pair_energy_matrix_kernel` | the kernel had the domain prefix while its general wrapper had the plain name -- backwards |

The docs moved with the names, which is the half that matters: a name a PINN
author does not recognise is as good as absent, and
`rotamer_rotational_correlation_time` additionally *claimed* something its
arithmetic does not do. 21 files, 988 tests pass.

7,245 -> 6,864 lines of Python in the `.i` files. `topology.i`, `observables.i`
and `av.i` are clear; `scoring.i` has four lines left.

## 2026-08-25 (later) — the dict view and the dict arithmetic

`parse_dye_mol2` copied `read_mol2_component`'s typed atoms into dicts of the
same fields under string keys, and `distance`/`angle_value` took those dicts
and did arithmetic IMP already does. All three are gone. A MOL2 atom is
`IMP.bff.Mol2Atom`, a distance is `IMP.algebra.get_distance`, and an angle is
`bond_angle_rad`.

**One angle primitive, two spellings, each with its unit in the name.** There
were three: `bond_angle_deg` in C++ returning 0.0 for a collapsed arm,
`angle_value` in `topology.i` returning the tetrahedral 1.9106 rad for the
same case, and `sim.i`'s `_angle`, which built coordinate dicts out of three
IMP particles so it could call the second one. `bond_angle_rad` is the formula
now and `bond_angle_deg` is `* 180/pi` over it; the fallback is the tetrahedral
angle, which is the choice this repository already recorded when it deleted the
*fourth* copy -- a harmonic whose minimum sits at a collapsed angle is not a
default anyone wants. It fires on nothing in the tree: 6,364 angles, no
zero-length arm.

The test that compared the dict view with the typed value went with the dict
view -- its subject no longer exists -- and what replaced it checks the reader
against real files, the 4,698-atom 1DG3 included.

7,295 -> 7,245 lines. 971 tests pass, built and run at `nice -n 19` with two
jobs.

## 2026-08-25 — R0 is Angstrom, the graph is C++, and a lesson about shared build trees

**R0 is Angstrom at the source.** It returned nanometres while every consumer
of it -- `fret_efficiency`, `av_distance`,
`AVPairDistanceMeasurement::forster_radius`, `LlFretOptions` -- worked in
Angstrom, so six call sites multiplied by ten and one nearly shipped without.
`forster_radius` converts once now; **the spectra stay nanometres**, which is
what a spectrometer, a datasheet and the overlap integral are in. Gone with it:
the x10 in `DyeContainer.cpp`, `InteractionTerms.cpp`, `bin/imp_bff_labelizer`,
`rotamer.i` and `rotamer_ensemble.i`, and the one in the container test that
existed to pin the conversion. `RotamerFRET(r0=...)` and `imp_bff --r0` are
Angstrom (54 by default), and the FRETpredict pin kwargs were restated
5.5 -> 55.0 and 5.68 -> 56.8 -- the same experiment in the package's unit, with
the reference E values untouched. The pins that cite FRETpredict's own R0 keep
it in *FRETpredict's* nanometres and convert at the comparison, because that
number is a citation.

**The molecular graph is C++.** `topology.i` had `build_graph` (a dict of
neighbour sets), `build_angles`, `build_dihedrals`, `find_cycles`,
`ring_atoms_from_graph`, `_component_without_edge` and two rotor helpers: 126
lines that numbered the caller's nodes, called `MolecularGraph`, and mapped the
answers back -- twice over, because a first-appearance numbering and a sorted
one give different orders out. The numbering belongs with the graph.
`MolecularGraph` gains `get_component_without_edge`, `get_bond_rotor`,
`get_angle_rotor` and `get_ring_atoms`; the new `LabelledGraph` does the same
for string-keyed nodes, which is what site ids are. `LinkerSampler`, `sim.i`
and four test files call the C++ directly. The sampler keeps **first-appearance
node order** explicitly rather than taking the C++ ascending order, because
that order reaches a seeded sampler whose weights are pinned -- and the pins
did not move.

7,381 -> 7,295 lines of Python in the `.i` files. 966 tests pass.

**And the process lesson, which cost more than the code.** Two sessions ran
`ninja` in the same build tree; ninja takes no lock, so the builds interleaved
writes to the same objects and `.ninja_log`. The visible symptom was a build
that **exited 0 having printed through [11/18]**, leaving a `_IMP_bff.so` that
did not contain the sources it had just built -- R0 still returning the old
value, a class missing from the module -- which reads exactly like a source
regression and was mis-diagnosed as one, twice, by two sessions. The machine
also went to load 45 on 8 cores. What replaces it: announce a build, one at a
time, `nice -n 19 ninja -j2`, and no `until ...; do sleep; done && ninja`
waiters -- an armed trigger is not a hold, and it cannot see that the tree went
quiet because somebody was rescuing the machine.

## 2026-08-24 (later still) — a third derivation of the dye topology, deleted rather than ported

`topology.i` had `build_dye_topology` and, under it, `_build_impropers` and
four `_build_*_impropers_from_template` wrappers: bonds, angles, dihedrals and
template impropers derived from a MOL2 graph. `build_forcefield_system` (C++)
derives the same four things, and `build_dye_protein_system` used to carry a
*second* copy of that derivation -- the one whose `impropers = []` cost 81
impropers and has its own page
([`validation/impropers_are_dropped.md`](validation/impropers_are_dropped.md)).
This was the third.

Nothing in the package called it. Its seven call sites were all in
`test/cgdye/test_dye_topology.py`, and what those tests assert about a
four-atom chain -- 3 bonds, 2 angles, 1 dihedral -- the three tests directly
above them assert against the graph itself; the only thing they added was
`distance > 0`. So it is deleted rather than ported, with the tests that
existed to keep it alive, and the note in its place says why. 100 lines of
Python and one more chance for the derivations to drift.

`dye_internal_system` went the other way, to `TopologyBuild.h`, and
`scoring.i` is down to four lines. 7,471 -> 7,381. One test is red on the
shared tree and it is the labelizer session's
(`test_dye_container.py::test_every_column_name_is_a_dictionary_item`, an
`mfdb_is_term` vocabulary check, file touched 21:03); they have been told.
958 pass otherwise.

## 2026-08-24 (late) — `dye_internal_system` to C++, and a wrong hypothesis measured before it was written down

`scoring.i` is down to four lines. `dye_internal_system` -- the topology-only
system `DyeInternalEnergyEvaluator` needs to know which site pairs are 1-2, 1-3
or 1-4 -- took `parse_dye_mol2`'s dicts, built site ids, and went out through
`forcefield_system_from_json`. It is `TopologyBuild.h` now and takes the typed
`Mol2Component` the C++ reader already returns, so the dict never exists: the
sites, the `MolecularGraph` angles and torsions and the exclusions are all
derived where the data is. Its three callers (two in `sampling.i`, one test)
took the typed value with it, and one of them lost a `re.match` that
re-derived the element from the atom name -- a rule the MOL2 reader had
already applied.

**The hypothesis that made this interesting was wrong, and measuring is what
said so.** The Python built its sites with `id` and `atom_name` only, so every
site reached the LJ table with an empty element: apparently a hydrogen
parameterised as carbon, in the dye's own internal energy. The port fills the
element from the MOL2 -- and the energy did not move, to four decimals. The
reason is that `site_element_map` never read `FFSite::element` at all: it
re-derived the element from the atom name in an inline loop, a **fourth** copy
of a rule that already exists as `element_from_atom_name` (C++), as
`structio::element_from_name`, and as that `re.match` in `sampling.i`. Nothing
was broken; something was written four times.

So the fix is the dedup, not the "bug": `site_element_map` now prefers the
site's own element and falls back to `element_from_atom_name`. Today those
agree on every structure in the tree -- 0 atoms differ across the shipped MOL2
set, because the reader derives the element with the same first-letter rule --
and they will part on a halogen, where the name rule calls `CL3` carbon and a
`_atom_site.type_symbol` does not. The dye's internal energy on alexa488_r48 is
18.8717 kcal/mol before and after.

7,495 -> 7,471 lines, 951 tests pass.

## 2026-08-24 (night) — the shadows are gone, and a typemap that lets a C++ signature say what it means

Three `%feature("shadow")` blocks re-implemented in Python what SWIG does
natively, and all three are gone.

`AVNetworkRestraint` (33 lines) and `LifetimeSpectrum` (26) hand-parsed
argument names, defaults and error messages because **SWIG will not generate
keyword arguments for an overloaded constructor** -- and both classes had a
second, default constructor that exists only for deserialisation.
`%ignore`ing that one (it is cereal's, and a restraint over no hierarchy is
not a restraint) leaves a single wrapped constructor, and
`%feature("kwargs")` then does the whole job, unknown-keyword `TypeError`
included. `LifetimeSpectrum` took the other route -- one constructor with
default arguments -- because its empty state is a real value the vector
template needs.

`AVMeanDistanceRestraint` (28) was not kwargs but a `hasattr` ladder turning a
Particle or an AV decorator into an index. That is `IMP::ParticleIndexAdaptor`,
which `IMP::core::AngleRestraint` and `DihedralRestraint` already take, so the
constructor takes it too and the ladder is deleted. **`av.i` and
`observables.i` now hold no Python at all.**

`rotamer_mean_field_weights_multi_dye` moved to `Scoring.h`: it was numpy
iteration over the C++ pair-energy matrix, which is the same fixed-point the
single-dye update already ran in C++, written a second time. The port needed
two typemaps, and they are the interesting part. A caller passes one
`(n_conf, n_atoms, 3)` array per dye, and SWIG's nested conversion walks the
outer sequence and then asks its traits for each element -- traits that only
know sequences, so a 3-D ndarray arrives as a sequence of 2-D ndarrays and
fails naming the whole nested type. `types.i` now converts a sequence of
arrays into `const std::vector<std::vector<double> >&` (the same
`PyArray_FROMANY` path the single-vector typemap has taken since 2026-08-19)
and returns one back as a **list of numpy arrays** rather than a tuple of
tuples of floats. The alternative was flattening at every call site, which is
what the Python version did and what a C++ signature should not require.

7,695 -> 7,495 lines. Cleared outright: `observables.i`, `av.i`. `scoring.i`
is down to `dye_internal_system`, which reads `parse_dye_mol2`'s dicts and can
only move when `topology.i` does. 951 tests pass.

## 2026-08-24 (evening) — dedup: four helpers with five copies, two bond rules, and the first program out of a `.i`

**Helpers.** `ends_with`, `file_exists`, `trimmed`/`trim` and `nan_value` had
between two and four definitions each across `src/`, and they were not all the
same function: `FPSIO`'s `ends_with` was case-sensitive and `TrajectoryIO`'s
was not, so `LIB.BCIF` was a trajectory and `X.PDB` was not a structure.
Extension dispatch should not depend on which file the caller sits in. One
copy each in [`include/internal/Text.h`](../include/internal/Text.h), pulled in
by `using` declarations so no call site moved; `trim`'s two callers took the
rename. The case-insensitive spelling won, which widens what `FPSIO` accepts.

**Two bond rules, measured before touching.** `perceive_bonds` (ZMatrix,
`r_i + r_j + 0.35`) and `infer_bonds` (StructureIO, `1.22 * (r_i + r_j)`) carry
*different covalent-radius tables* -- H 0.37 against 0.31, O 0.73 against 0.66
-- and different cutoff rules. On the 33 dye templates and on 1DG3's 4,698
atoms they return **identical bond sets**, and `perceive_bonds` reproduces the
A48_C1R MOL2 table exactly (87/87). So this is latent, not live: the cutoffs
differ by up to 0.25 A (C-H 1.49 against 1.31, S-S 2.39 against 2.56) and a
stretched bond would fall on different sides. Not merged, because merging is a
decision about which perception rule the package has, not a cleanup.

**`faspr` is flat in `src/`** -- `src/faspr/*.h` became `src/Faspr*.h` beside
the `.cpp` files they belong to. The vendored headers include each other by
their upstream names, which is what broke the shared build tree for the
rotamer session for ten minutes; `utility/port_faspr.py` writes the flat layout
now, so a re-vendor reproduces it.

**`observables.i` is clear** -- the first of the ten named files with no Python
left in it. Its 26-line `%feature("shadow")` re-derived argument names,
defaults and error messages by hand because SWIG cannot generate keyword
arguments for an *overloaded* constructor. `LifetimeSpectrum` has one
constructor with default arguments now, so `%feature("kwargs")` does it
natively -- including rejecting an unknown keyword with a proper `TypeError`.
The same trick is available for the shadows in `av.i` (31 lines) and
`avmeandistance.i` (111).

**`scoring.i` lost its program half.** `torsion_cosine` and
`build_dye_restraints` -- 73 lines that were *literally* IMP C++ API calls
written in Python (`DistanceRestraint`, `AngleRestraint`, `DihedralRestraint`,
`Harmonic`, `HarmonicLowerBound`) -- are `Scoring.h` now, and with them the
CHARMM-to-`Cosine` phase shift, which is physics that had no business in an
interface file. The builder takes parallel `site_ids`/`particles` rather than a
dict, so nothing has to teach SWIG about a map of particles. Two test
adjustments came with it and both are the platform's own contract: a torsion
type is the typed `FFTorsionType`, not an ad-hoc dict with a `phase_rad` key,
and restraints come back from an `IMP::Restraints` container as base
`IMP.Restraint` -- IMP's own `RestraintSet.get_restraints()` does exactly the
same, so counting by `isinstance` only ever worked because the objects had been
built in Python.

7,695 -> 7,633 lines of Python in the `.i` files, 886 tests passing. What is
left in `scoring.i` is the multi-dye mean-field loop and `dye_internal_system`,
which reads `parse_dye_mol2`'s dicts and moves when `topology.i` does.

**Then, once it was said the package is pre-release, the two things that had
been left for compatibility were done properly.** `score_model` is
`0.5 * chi2_score` -- what a Gaussian is worth, \(-\log L = z^2/2\), and what
`IMP::core::Harmonic` scores for \(k = 1/\sigma^2\). It read 0.25 because one
factor of a half was applied twice, which had this restraint entering a scoring
function at *half the weight of every IMP harmonic beside it*; the four pins
doubled. And the second covalent-radius table is gone: `infer_bonds` used
Cordero radii under a multiplicative rule while `perceive_bonds` used
Pauling-ish radii under an additive one, agreeing on every molecule in the tree
but not on the cutoffs (C-H 1.31 A against 1.49 A). `covalent_radius` and
`bond_tolerance` are published from `ZMatrix.h`, `infer_bonds` takes a
`tolerance` rather than a `scale`, and the bond sets are unchanged on all 33
templates and on 1DG3.

`structureio.i` gave up `load_structure_with_particles` and `read_angle_file`
(38 lines) for `LoadedStructure` and `FlexFitSelection` in `StructureIO.h` --
the flexfit block is JSON, which this module parses in C++ everywhere else, and
decorating a particle is not a reason to be in Python. Both are IMP value types
with by-value getters, because SWIG rejects a non-const reference to one and
says so in a compiler error naming the rule. 7,633 -> 7,595, 919 tests pass.

## 2026-08-24 — the rotamer path, 7x on loading and 2x on scoring, and what a site does to a library

Three fixes, each found by measuring rather than by guessing, and none of them
changes a number: the FRETpredict pins and the 883-test suite are untouched.

**`read_drot` read the whole container to look at four bytes.** It opened the
file, slurped every byte to check the EBML magic, and only then walked the
framing — harmless when the only envelope was a tar that had to be read whole
anyway, and *expensive* the moment libraries moved into a 19.7 MB family
container, because it read all of it to fetch one library. Four bytes now.
Pulling one library out of `dyes.drot.pto` went 20.5 -> 10.1 ms CPU, and a
one-conformer library 10.2 -> 1.7 ms.

**The reader's framing walk was syscall-bound.** Crossing 761 objects took
~9,000 reads of a dozen octets each. `PtoReader` now serves small reads from a
512-octet window — big enough for one `AttachedFile`'s whole header run, small
enough never to drag a payload in. Listing a catalog went 12.9 -> 1.4 ms.

**Scoring spent 46 % of its time formatting numpy arrays into strings.** The
profile said `arrayprint`, 240 calls, and the cause was
`metadata.get("library_name", metadata.get("name", str(library)))` in
`RotamerEnsemble.from_site`: `dict.get` evaluates its default whether or not
the key is present, and `library` here is the library *dict*, so every site of
every frame rendered 711x83x3 coordinates to text and threw the string away.
Twice. pp11's `trajectory_analysis` went 1200 -> 562 ms.

Together: the whole 95-library corpus loads in **0.51 s against 3.65 s**, one
library in 9.4 ms against 36.3.

**One negative result, and a caveat on it.** The decode walks column-major
grids with a stride of `n_rot` doubles, which looks like a cache miss per row;
blocking the conformer loop so each row reads a run of neighbours was measured
at block sizes 1, 4, 16, 64, 256 and 1024 (a block of 1 *is* the plain order)
and the whole spread was under 10 %. The blocked form is kept at 16 because
writing straight into the output drops a per-conformer buffer and its copy.

The caveat is the method, not the result: those six numbers were taken as six
*sequential* runs on a machine at load 13, so "under 10 %" cannot distinguish
no effect from an effect smaller than the drift between runs. The honest way
to A/B on shared hardware is to interleave the arms inside one run in short
alternating blocks and compare medians -- contention and thermal drift then
hit both arms equally, so the ratio survives even though the absolute numbers
do not. Re-take it that way before leaning on it. (The method is chimol's,
learned the same day from the same mistake.)

What is left in the scoring path is genuine all-pairs work — 409k donor x
acceptor pairs per frame for pp11 — and reducing it means changing the answer,
which belongs to the sampling question below rather than to optimisation.

## 2026-08-24 — a dye library is a fixed sample, and the site decides how much survives

Measured, because "the sampling is site dependent" turned out to be an
understatement. A library is one sample of the *free* dye; a site reweights it
and never resamples it, so the site's answer is worth what the overlap between
that sample and the site's allowed conformations is worth. Kish's effective
sample size measures exactly that and is now on every ensemble as
`RotamerEnsemble.effective_sample_size`.

**Hsp90 residue 637 with the default cutoff-30 library is an ensemble of 1.7
conformers**: seven exist, two survive the site, one holds 72 % of the weight.
Every orientational quantity that site reports is that one conformer's. The
same site at cutoff 10 carries 71.7.

It costs the answer: at that site cutoff30 gives `Es` 0.4677 where cutoff10
gives 0.3999 — **0.068 in E**, and `<kappa^2>` moves 0.97 -> 0.65. At pp11's
open sites the same comparison moves E by 0.008. That factor of ten between two
structures, same library and same code, *is* the site dependence. The
`<kappa^2>` of 0.97 at the tight site is not order in the dye; it is two
conformers being unable to represent a distribution of orientations.

Written up with the full table in
`okf/validation/site_dependent_rotamer_sampling.md`, including why the shipped
default is nonetheless right for parity (it is FRETpredict's default, and the
pins prove this package reproduces FRETpredict — including where FRETpredict is
under-sampled). The real fix is resampling rather than reweighting: a `.drot`
conformer is a dihedral vector, so conformers can be generated near the
survivors and rebuilt through `internal2cartesian`. That changes numbers by
design, so it is a decision and not a default.

## 2026-08-24 (PRD-120: the Labelizer ported, and what the A/B found)

* **The published label-site score is native, and it could not be run here
  before.** `../labelizer-backend`'s `ss` shells out to DSSP and raises
  `RuntimeError("Unknown platform")` on macOS; `se` and `me` shell out to MSMS
  binaries that are 32-bit ppc/i386 Mach-O; and it pins LabelLib, banned since
  2026-08-11. Meanwhile `IMP.bff` had every piece of physics underneath it —
  AVs, R₀, κ², PET quenchers, SASA, charge masks — and **no per-residue scoring
  layer at all**. `LabelizerFeatures.h`, `LabelizerScore.h`, `LabelizerFret.h`,
  `LabelizerIO.h`, `PtoProfile.h`, `internal/Sha256.h`, `bin/imp_bff_labelizer`.

* **The native DSSP reproduces the reference binary exactly — 195/195 on
  1DDB, zero score delta — and getting there found two real things.** First,
  DSSP marks the residues *between* an n-turn's hydrogen-bonded pair, `k+1 …
  k+n-1`, not the donor `k`; including `k` over-assigned `T` and cost nine
  residues. Second, the remaining four disagreements were all `I → H` in one
  stretch: **DSSP 3.0 reversed the π/α precedence** (Touw 2015), so a 5-turn
  now wins where a 4-turn overlaps it. Assigning `I` first took it to 195/195,
  which also dates the reference's CSVs to DSSP ≥ 3. Worth knowing if these
  numbers are ever compared against an older binary.

* **The shipped 1DDB example feeds its own output back in as its input.** Its
  `cs` reference has **two distinct values across 195 residues** where the
  table has ten bins (`cr` has 20, `se` has 10). The example passes
  `prot1_cs=[".../1DDB-39_cs.pdb"]` and `_save_pdb` writes the scores to that
  same path, so each run reads the last run's output as conservation grades —
  and the lookup has a two-cycle (1.63 → 2.3553…, 2.36 → 1.6322…, each
  rounding to the other's grade), which is the fixed point it is stuck on. Our
  lookup reproduces both values to sixteen digits, so the machinery is
  verified; the input that produced them is one generation back and is not in
  the distribution. The real ConSurf file for the entry correlates with the
  implied bins at **r = +0.07**, i.e. it is not the one either.

* **Residue depth against MSMS: no bias worth correcting.** Mean +0.015 Å,
  r = 0.976, slope 0.953, 92.3 % within one table bin, 74.9 % bin-exact. The
  0.163 Å scatter is close to the quantisation floor the binned reference
  imposes (0.077 Å); the genuine disagreement is ≈ 0.14 Å. The native surface
  is the **contact** surface only — MSMS includes the reentrant patches — which
  is where the rest should live, in crevices.

* **The cheap dye model over-reaches, and it matters more than it looks.** The
  analytic alpha cone screens 15 051 pairs in 0.01 s but places the dye 2–5 Å
  too far out; because the pair score peaks sharply at `R = R₀`, that is a
  ~30 % score error, and the top five pairs by cone score are **not** the top
  five after rebuilding the clouds (A3–A30: 1.8931 @ 52.4 Å → 1.3606 @ 47.0 Å).
  Do not leave `n_refine` at zero.

* **Absence replaced the sentinel, and the term layer earned its keep on the
  first run.** The reference writes `-1` for "excluded" and `0` for "no
  contribution" into the same column as real scores, so its CSVs cannot be read
  back; a position with no score now carries a `status` and no number.
  `mfdb_check_term` then refused to write `row_grain = label_site` because the
  dictionary had no such term — which is the layer working, not failing — so
  `label_site` went into `../mmfdb/.../mmfdb_flr_ext.dic` and the version to
  **1.8**. Left uncommitted there; that repo has other work in flight.

* **Attribution became dictionary terms, not header JSON.** Seven items on
  `_mmfdb_artifact` in the same 1.8 pass, at the `.drot` session's request:
  `author`, `citation`, `license`, `terms`, `terms_url`, `source`,
  `redistributed_via`, with `MfdbAttribution`/`mfdb_attribution_tags` over
  them. The reasoning is theirs and worth keeping: attribution is **per
  artifact** because one family container can hold GPL-3.0-only and
  academic-use libraries at once; **`license` must not be mandatory-SPDX**,
  since "free for academic use" has no identifier and coining a `LicenseRef`
  would state something the upstream did not — so at least one of
  `license`/`terms` is required and neither is preferred; and
  `redistributed_via` exists because Dunbrack-2010 reaches this stack through
  FASPR, which is MIT while the library is not, so a chain recording only the
  immediate source loses that the licences differ along it. Spelled `license`
  (US) to match the upstream mmCIF dictionaries.

* **Two kernels extended rather than copied.**
  `solvent_accessible_surface_area_per_atom` joins the existing SASA: the old
  one samples every atom at one radius and scales by `vdw²`, which is right for
  contact estimates and wrong for a *relative* accessibility, and the header now
  says which answers which. PRD-109 wrote that kernel and nothing ever called
  it; this is its first consumer.

* **Coordination.** `Pto.h`/`Pto.cpp` are `imp-bff-69`'s and were used as-is —
  no second container core, no tttrlib dependency. Kinds split by agreement:
  `drot.*` and `rot.bbdep.*` theirs, `label.*` here. They asked for the
  vocabulary layer to be reusable rather than labelizer-scoped, hence
  `PtoProfile.h`. Tags stay additive: a walker plus the framing reads a
  container with none of this code. `prototypes/fast_label_score/` — the
  *improved* model — is untouched; the original had to land first so there is
  something to measure an improvement against.

* **The container passes `../tttrlib/test/tools/pto_ebml_check.cpp`** — libebml
  itself, strict alignment included (5 payloads, 0 misaligned). The checker
  needs **libebml 2.0**, not the 1.4.7 Homebrew ships: `EbmlId::FromBuffer` is
  2.0-only and 1.4.7 puts `EMaxSizeLength`/`EDocType` behind an
  `EbmlSubHead.h` the tool does not include. Build it from the clone at
  `../tttrlib/junk/libebml` (recipe in `test/label/test_labelizer_pto.py`).
  I first read those two breaks as the tool being stale; that was backwards,
  and `Rotamer` corrected it — it targets the newer library and Homebrew is
  behind. The suite does not take that dependency, so the container is *also*
  walked by a from-scratch RFC 8794 parser in `test/label/test_labelizer_pto.py`
  that shares no lineage with either implementation.

* **A second A/B closed the biggest gap, and found two more defects.** 1DDB is
  all-helical, so the strand half of DSSP was never *compared* — only
  exercised. The paper's Supplementary Data 1 also carries maltose-binding
  protein in two conformations (1OMP apo, 1ANF holo, 370 residues, 72 `E` and
  5 `B` in the reference's own column), and against it `cr` is **370/370** and
  `ss` is **97.3 % on both**. Getting from 94.6 % found: **β-bulges were
  fragmenting sheets** — ladders were joined only between strictly adjacent
  bridges, so a bulge broke a strand into isolated bridges (`E → B`), now a
  union-find joins any two bridges of a type that advance together with one
  side stepping exactly one; and **`G` was a leftover stub** — assigning
  whatever remained after `H` took part of the span gave one- and two-residue
  3-10 helices at helix C-termini where the reference has `T`. The ten interior
  strand residues that remain are **not** a threshold artifact: −0.5 kcal/mol
  is optimal (94.3 % at −0.4, 94.9 % at −0.6), so what is left is DSSP's
  sheet-level bookkeeping. Fixture provenance: extracted from the published
  supplement, five decimals, so comparisons are at 1e-5 and that is the file's
  precision rather than a tolerance on the port.

* **The two-state layer got its first real test and passes a physical one.**
  It had no tests at all. On MBP it ranks the hinge closure first — the top
  twenty two-state pairs all move more than 5 Å and eighteen contract, by
  11–13 Å — and **A29, in the best pair, is a site the paper itself labelled**
  (Figure 4). The Cβ difference map peaks at 13.8 Å over the same motion, is
  symmetric, and is zero on the diagonal. A structure paired with itself scores
  exactly zero, which is the check that the two conformations are read
  independently.

* **`ll_dssp` no longer bridges a sequence gap.** Array-adjacent is not
  chain-adjacent: a dropped residue (an unnatural amino acid, a modified
  residue) or an unresolved loop leaves a hole in the author numbering, and the
  turn logic was indexing by array position — reading residues 99 and 103 as
  three apart when a missing 100 had made them four, and asserting a 3-turn
  that does not exist. Confident wrong output, not an error, and invisible on
  all three fixtures because none of them has a gap. Raised as a hypothesis by
  `fast-label-score-77` (whose bundle has 982 pAcF measurements over 75 sites);
  reproduced by deleting one mid-chain residue, then fixed and pinned. Scores
  were never at risk — they key on `(asym_id, seq_id)` — so only a caller
  zipping arrays by index could have been misled.

* Suite: **886 passed, 3 xfailed**. `okf/validation/labelizer_ab.md`,
  `okf/prds/prd-120.md`.


## 2026-08-24 — one container per rotamer family: dyes, spin labels, side chains

The store now holds all three rotamer families it should always have held.
`data/rotamer_library/` is three files: **`dyes.drot.pto`** (95 FRETpredict
dye+linker libraries, 19.7 MB), **`spinlabels.drot.pto`** (10 DEER-PREdict spin
labels, 0.24 MB) and **`sidechains.drot.pto`** (the Dunbrack-2010
backbone-dependent table, 3.3 MB). Before today only the dyes shipped, spin
labels merely converted, and side chains were not in the format at all.

**Bundling costs nothing and reading stays cheap.** `write_drot_bundle` copies
payload bytes across unread under `<library>/` names and adds a `drot.catalog`;
95 libraries bundle to 19,675,419 bytes against the 19,683,305 they occupied
separately. A library is addressed by a locator —
`dyes.drot.pto::A48_C1R_cutoff10` — and `PtoReader` now walks the framing by
seeking and reads only the objects asked for, so pulling one library out of the
family touches its own bytes. That is what makes one-file-per-family a
different proposition from the corpus tar rejected in the PRD a day earlier.

**Side chains are the case that tests the shared grammar.** A dye library is an
ensemble; a side-chain library is a distribution over (phi, psi). Same
envelope, new kinds — `rot.bbdep.header`, `rot.bbdep.records` — and no element
ID invented, which is the claim PRD-040 makes about a new domain. The records
keep FASPR's own 20-byte layout byte for byte, so `write_dunbrack_bin`
reproduces `dun2010bbdep.bin` with an identical SHA-256, `faspr_pack` accepts
the shipped container and packs byte-identically to packing from the binary,
and the decode agrees with the prototype's independent reader over 720
(residue, phi, psi) cases exactly. 13.76 MB → 3.31 MB.

The Dunbrack library is academic-use. It is redistributed here on the owner's
explicit instruction, after the terms were put to them; the attribution,
citation and terms ride *inside* the container in `bbdep.json`, where they
cannot be separated from the data.

## 2026-08-24 (later) — the property sugar goes native, and a 100x chi2 disagreement falls out of looking

`avmodel.i` had already decided how a C++ value shows itself to Python --
`%attribute` for scalars, `get_*()` views for arrays, "no property sugar, no
`%pythoncode`" -- and then `rotamer_ensemble.i` hand-wrote eight property pairs
over the same `States` getters anyway. They are gone: `points`, `orientations`,
`mu`, `attachment_point`, `mean_position`, `position_name`, `params` and
`n_points` are declared on `States` itself in `avmodel.i`, so an ensemble, an
`AccessibleVolume` and an `ACV` now answer alike and the two helpers that
branched on `isinstance(RotamerEnsemble)` are one line each. The shape is the
part worth having in one place: a cloud is `(n, 4)` and a dipole set is
`(n, 3)`, and this repository has paid twice for a caller reshaping by hand.
`%attribute_np2` (new, `types.i`) carries that reshape; `%attributestring` --
not `%attribute` -- carries `position_name`, because the getter returns the
string by value and `%attribute` would hand out the address of a temporary.
71 lines of Python left the `.i` files (7,676 -> 7,605) and 798 tests pass.

**The finding is worth more than the line count.** The next record in the
queue, `PairDistance` in `docking.i`, turned out to be the Python duplicate of
the C++ `AVPairDistanceMeasurement` -- and the two implement the asymmetric
chi2 in *opposite* directions. `score_model` (C++, what every restraint and the
docking engine minimise) scores a model above the experiment against
`error_neg`; `PairDistance.chi2` (Python, what the FPS tables report) scores it
against `error_pos`. On experiment 50 A with errors (1, 10) that is 25.0 versus
1.0 one way and 0.25 versus 100.0 the other -- and no test asserts either
number: the one test that calls `score_model` throws the result away.
Measurement and all three listings in
[`validation/two_chi2_conventions.md`](validation/two_chi2_conventions.md).

**Resolved the same day, on the maintainer's ruling: the residual is model
minus data -- a model distance that is too large is a positive deviation, and a
positive deviation is judged against `error_pos`.** Looking for the second
implementation turned up a *third*: `chi2_score` (`AVDistance.cpp`), which
`LabelingRestraints` and the `imp_bff` chi2 columns already use, and which had
it right all along. So the ruling cost `chi2_score` nothing but a comment, made
`AVPairDistanceMeasurement::score_model` delegate to it -- it had measured
`data - model` and then read the branch as if it were `model - data`, which is
how the bars came out swapped -- and turned `PairDistance`'s four properties
into calls to `chi2_score` and `fret_efficiency`. Three implementations became
one. The 0.25
`score_model` carries was left alone deliberately: it is a scale, and every
restraint weight in the package was set against it, whereas a branch reorders
models. Four of the package's own regression pins moved with the fix
(13.505 -> 11.350 and three like it, listed in the note) while the model
distances in the same tests did not -- which is what a branch fix should look
like -- and `test_distance_conventions.py` now asserts which bar divides in
each direction and that `score_model` is `0.25 * chi2_score` across a sweep.
863 tests pass.

## 2026-08-24 — `.drot` folds into `.pto`: one container for the stack

The rotamer libraries stopped having a private container. Their members —
header, template, Z-matrix rows, the four grids, the weights — are now the
attached objects of a **PTO** document (EBML, `DocType "pto"`), the envelope
tttrlib writes photon streams into and chimol writes structures into, and a
shipped library is `<stem>.drot.pto`. The suffix is the profile convention
PRD-040 settles: it names the primary payload kind, so a shell and a human see
which profile a container carries while a reader still asks each object what
it is. No element ID was invented — the eight members are `PtoKind` strings
(`drot.header`, `drot.template`, `drot.rows`, `drot.grid` ×4, `drot.weights`) —
which is the test of the claim that a new domain is a set of kinds, not a
dialect.

**Vendored, not depended on.** `include/Pto.h` + `src/Pto.cpp` (~570 lines)
were written from `../tttrlib/okf/specs/pto-binary-decoding.md`, not imported
from tttrlib: the convergence record's rule is one small core per repo with a
provenance note, and imp.bff acquires no build dependency for a container. It
is the third independent implementation of the same grammar, after tttrlib's
C++ and chimol's Python.

**The envelope is free.** Each object is brotli'd on its own instead of one
stream over a tar, which gives up the joint context across members and gets
tar's 512-byte headers back: over the 95 shipped libraries **19.60 MB against
19.64 MB**. Per file it is +0.5 %; the honest outlier is the one-rotamer
library, which pays 35 % of its 2 KB for the framing. Bought with it: payloads
on 8-byte boundaries for a mapping reader, a header that inflates without the
grids, and objects addressable by name and kind.

**Proved by readers that are not ours**, which is the only proof a shared
format admits — writer and reader agreeing on a misreading of RFC 8794 would
satisfy any test written here. `pto_ebml_check` walks all 95 with **libebml**
(`--aligned` enforced), and **chimol's `render/pto.py` reads imp.bff's
containers** and agrees on every object's kind, offset and payload — the first
time two implementations of this format have been pointed at each other's
files. Round-trip against the v9 files is bit-exact (0.0 Å), re-writes are
byte-identical, the FRETpredict pins are unmoved and the suite is 797 green.
v5–v9 still read; only v10 is written.

Also this session, before the fold: `resolve_rotamer_library_path` and
`load_rotamer_library` learned the `.drot` store (the loader normalises its
raw cluster populations, which the frame stores already did), `find_reference_
rotamer_files` prefers it, and the C++ `load_rotamer_library` now dispatches on
the extension — which incidentally repaired `.dcd`, dead in that path since it
was pointed at `read_bcif_trajectory` alone. The DEER-PREdict spin-label
libraries (10, in `junk/`) were run through the same path as a second family:
they encode, read back and load, worst error 2.1e-06 Å.

## 2026-08-24 — PRD-118: `.drot` ships, brotli is vendored both ways, and the pins hold

The rotamer-library store is finished end to end. `read_drot` was in the tree
but nothing could **write** one, the container linked against a brotli that was
not declared, and the shipped `.drot` set failed the FRETpredict parity gate --
10 red tests in `test/cgdye`. All three are closed.

**Vendoring, both directions.** `src/brotli/` is google/brotli master
(`8e10eeb`, MIT), decoder *and* encoder, compiled as the single hidden
translation unit `src/Brotli.cpp` -- `dependencies.py` declares nothing, no
`-lbrotli*` reaches the link line, and `nm -gU libimp_bff` shows zero `Brotli*`
symbols, so a process that also loads a real libbrotlidec keeps the two apart.
Three mechanical rules make it work in an IMP module and are written down in
`src/brotli/VENDORING.md`: implementation files carry `.inc` (a `.cpp` under
`src/` would be compiled twice -- once inside `bff_all.cpp`, once in per-file
mode -- because **`setup_all.py` globs `src/*.cpp` one level deep and ignores
`Files.cmake` in the default build**); `<brotli/x.h>` is rewritten to the
relative header so the vendored one always wins; and the three file-local
helpers two encoder units share are renamed around their `#include`. The
codec is byte-identical to the reference brotli in both directions.

**The format had to stop losing; the pins did not have to move.** v8 rebuilt
conformers on the *template's* bond lengths. Real MD flexes bonds by ~0.004 A
RMS, and measured against the shipped `.bcif` that was 0.016 A of the 0.020 A
total error -- 1.3e-3 in E against a 2e-5 tolerance. The fix is representational,
not numerical: store the conformer's own bond length per Z-matrix row and the
description is complete (`ZMatrix::frame_internals`/`decode_internals`, exact to
~1e-13). What is left is the grid, and **v9 stores float32 internals,
column-major, byte-plane shuffled**: 1.2e-6 A mean / 9.3e-6 A max, 0.66x the
`.bcif`. The 95 shipped libraries were re-encoded and verified on write (19.6 MB
against 30.0 MB), the parity pins pass untouched, and `test/cgdye` is green.
The int16 rung survives as `--grid` (0.30x, ~3e-3 A) carrying the warning
`imp_bff_traj2bcif` already carries: kappa^2 comes from dipole directions
between atoms ~1.7 A apart, and those errors do not average away.

**The builder is a program**: `bin/imp_bff_traj2drot` (`.bcif`/`.dcd`/`.xtc` +
template PDB -> `.drot`; `--cluster A` takes PRD-118's raw-MD path through
`cluster_frames_leader`, `--all DIR` re-encodes a corpus, and every write is
round-trip verified unless told otherwise). Two smaller things fell out of
using it: `_metadata_from_path` never matched a `<stem>_cutoff<N>` file, so an
explicit library path lost the dye's dipole and attachment selectors and
`RotamerEnsemble.from_site` raised `IndexError` -- a path is a library name
again, for `.bcif` as much as `.drot`. And `load_rotamer_library` in
`DyeSampling.cpp` called an `ends_with` that did not exist (a concurrent
session's in-flight edit), which is what "the module is disabled" looked like
that hour.

Documented in `README.md` (the tool section `setup_module.py` insists on),
`doc/manual/structure/structure_cgdye.ipynb`, the PRD-118 ledger, and
`examples/structure/drot_rotamer_library.py`, which runs trajectory -> library
-> two ensembles on hGBP1 -> FRET in one file. Tests: `test/io/test_drot.py`
(11: round trip, self-containedness, reproducible bytes, the shipped corpus
against its `.bcif`, the builder, the clustering path, the grid a molecule does
not fit). Suites: **796 passed, none failed** in the non-medium run,
`test/cgdye` 161, medium 43.

**Three things the corpus change made visible, all older than it.** The dye
CLI joined `str / str` in `find_dye_structure`, `find_dye_mol2` and
`build_lib` -- `get_structure_dir`/`get_template_dir` became C++ and return
`std::string`, so four commands of `imp_bff dye` had been dead since -- and
`sample-rotamer` still passed `rng=` to a `sample_rotamer_index` that takes a
seed. `test_AccessibleVolume.py` read its reference maps from
`./references/...`, so the AV feature test passed from `test/` and failed from
the repository root, where the suite is actually run; it had been recorded as
"a pre-existing data failure" for that reason. And
`expensive_test_dcd_reader.py` swept `data/rotamer_library` for `.dcd` files,
of which there have been none since 2026-08-19 -- it now finds the one kept at
`test/input/A56_C1R_cutoff30.dcd` and the MDAnalysis parity is checked again.

**And one test that could only ever have passed by luck.**
`expensive_test_langevin_equilibrium.py` (PRD-108, no rotamer library in
sight) asserted `KS < 0.3` on the md-vs-bd centroid-distance distributions and
failed at 0.3151 -- reproducibly, since it is seeded. That is not a near miss:
**two md runs differing only in their seed are 0.3151 apart by the same
statistic**, and md-vs-bd over three seed pairs gives 0.315 / 0.555 / 0.375.
At 0.4-1.2 ns the linker's centroid distribution has not converged, so a fixed
threshold below ~0.5 was testing the random number stream. The test now runs md
twice and asserts that bd sits no further from md than md sits from itself
(`ks <= 2 * floor`), with the measurement written into its docstring; whether
the two integrators agree once *converged* needs runs an order of magnitude
longer than a test should hold, and stays a PRD-108 question.

## 2026-08-22 — PRD-117 batch 18: `build_forcefield_system` to C++ (`TopologyBuild.h`)

The 410-line Python orchestration is one C++ function now:
`build_forcefield_system(components_json, ...force constants...)` reads the
MOL2 components, derives sites (alpha-suffix dedup of repeated MOL2 names --
`serial_to_site_atom_names` is public too), bonds/angles/dihedrals, the
template-driven impropers (ring/pi/flat/orient centers expanded against the
bond graph), the feature groups, the LJ table (via `build_lj_type_table`, one
source), and the type tables -- then builds the typed system through
`forcefield_system_from_json`, the one conversion path, which guarantees the
same semantics the Python produced.

The Python `build_forcefield_system` is a thin bridge (list-or-string specs ->
JSON -> `_IMP_bff`); `build_system_from_specs`, `build_dye_protein_system`,
`dye_forcefield_system` and the dict/`parse_dye_mol2` surface stay Python
(dict-shaped callers; the writer wrapper). `topology.i` 989 -> 628 lines.

**The bug the port caught**: `MolecularGraph` keys nodes by the raw pair
values (the serials), but my first improper expansion passed 0-based index
arrays beside a serial-keyed graph -- 58 impropers where the pins say 81.
Every count matches after passing serial spelling everywhere (138 sites, 146
bonds, 262 angles, 207 dihedrals, 81 impropers on the shipped two-component
system).

Verification: `ninja IMP.bff` clean; non-medium **624 passed, 3 xfailed** (only
the pre-existing `av_reference_0.mrc` data failure); cgdye **160 passed**;
expensive all-dyes **36 passed**; `imp_bff build-system --help` runs.

## 2026-08-23 — PRD-118 scope: spin labels (DEER/PRD + PRE) — DEER-PREdict cloned, `.drot` verified on its libraries

The label story now extends beyond fluorescent dyes: cloned
[KULL-Centre/DEERpredict](https://github.com/KULL-Centre/DEERpredict) (de-gitted,
176 MB of test MD stripped to 22 MB) into `junk/DEERpredict/` — Tesei, Martins,
Kunze, Wang, Crehuet & Lindorff-Larsen, *PLOS Comput Biol* 2021;17:e1008551,
GPLv3; **the same group and the same pdb+dcd+weights mdtraj triple as our
FRETpredict dye libraries**. Its `lib/` carries 10 spin-label rotamer
libraries (MTSSL 175K/298K variants, BASL, MA-proxyl) — and the `.drot`
converter handles them as-is: **MTSSL CaSd → `.drot` at 0.0023 Å RMSD
(rigid rotamers, cleaner than the dye MD), 17× smaller than its DCD**.
One gotcha found and recorded: their PDB templates are zero-geometry
placeholders, so the template comes from the first DCD frame
(`CONVERTING.md` gotchas). What this opens: EPR (DEER distance
distributions) and NMR (PRE rates) labels on the same rotamer machinery as
FRET dyes — one `.drot` store, one loader, two spectroscopies.

# Update Log

## 2026-08-26 (night, last +7) — a use-after-free that returned zeros

Chasing the last of the untested surface -- the two `bin/` programs nothing
covered, and the top-level commands -- turned up two real defects and one that
had been silently wrong for as long as the CIF writer has existed.

### A system CIF written anywhere else could not be read back

The writer records each component's MOL2 path *relative to the CIF*, and
computed the `..` count from the **separators** of the base directory rather
than its segments -- one short whenever the base has no trailing separator,
which is always:

```
base   /var/folders/cl/xxxx/T/tmp1        6 segments
mol2   /Users/me/dye.mol2
wrote  ../../../../../Users/me/dye.mol2   5 ups
read   /var/Users/me/dye.mol2             does not exist
```

Nothing noticed because every test writes its systems into the data tree,
where the common prefix is long and the miscount cancels. `imp_bff simulate`
was dead for any system built outside it. `relative_path` in
`TopologyBuild.cpp` counts segments now, and
`test/cgdye/test_system_paths.py` writes a system into a `tmp_path` -- which
is exactly the case that failed.

### A container read off a temporary was freed memory

Writing the test for that found something worse:

```python
IMP.bff.read_forcefield_cif(path).components.values()
```

returned `'\x00\x00\x00...'` where a path should be, and segfaulted under
pytest. An accessor returning `const std::map&` gives SWIG a **borrowed**
pointer into the owner, and the proxy does not keep the owner alive: the
system is a temporary, it dies as the expression unwinds, and the map proxy
outlives it. Silently, most of the time.

The accessors cannot return by value -- `Scoring.cpp` calls `get_bonds()` in
the *condition* of a loop over the bonds, so a copy per access is quadratic --
so the copy is made **at the language boundary**: `%owned_container_out` in
`IMP_bff.types.i` gives every `const container&` return an owning copy, which
costs one copy per Python access and hands the lifetime to Python. Eight
container types: the maps, and the vectors whose element is one of this
module's own values.

**Not** `std::vector<std::string>`, `<double>` or `<int>`. Those already come
back *converted* -- a Python list, a numpy array -- which is a copy, so they
were never borrowed; claiming them replaced the conversion with a raw proxy and
five tests stopped recognising their own results (`assert <RMF_HDF5.Strings> ==
('SD',)`). The first cut of this fix did claim them, and the suite said so.
`%naturalvar` is on beside it, for the same hazard in plain member variables.

`test/test_borrowed_containers.py` reads five kinds of container off a
temporary, and checks that the copy really is one.

### Two programs nothing had ever run

`imp_bff_dye_pdb2cif` passed `--dye-id`'s `None` default into a
`const std::string&`, so its documented default invocation raised a
`TypeError` before reading anything. `imp_bff_labelizer --show`, given a
structure instead of a container, failed with "does not begin with an EBML
header" -- true, and useless. Both fixed, both now in
`test/io/test_bin_programs.py`.

`imp_bff dock` runs (score 188.30 over 99 distances, 17 accessible volumes),
and so do `build-system`, `simulate` and the six `dye` commands.

## 2026-08-26 (night, last +6) — nine dead paths, found by running the docs

The unit suite was green at every step of the last three passes, and the CLI
and the example gallery were quietly rotting the whole time: **what breaks when
a C++ port changes a shape is the code nobody runs**, and nobody runs the
documentation.

Running every `bin/imp_bff` command and every example and notebook found nine.

### In the program

| command | what was wrong |
|---|---|
| `dye sample-dof-walk` | imported `LinkerSampler`, a name from before the linker sampler became C++; read `resolve_probe_site` as a dict; and then, once running, **accepted nothing** |
| `dye sample-langevin` | passed `timestep_fs=None` into a `double`, and asked `run()` to write an RMF, which it stopped doing when the sampler became C++ |
| `dye label-fusion` | `sys.exit` in a module that never imported `sys` |
| `flexfit` (FP branch) | `attach_probes` with four parallel lists, and another dict-style site |

`sample-dof-walk` is the interesting one. It rejected every proposal because
the probe is *bonded into* the site: its first atoms sit a bond length from
residues i-1 and i+1, and the walk counted those bonds as clashes. Excluding
the site and its two neighbours, and accepting a move that does not make the
clash count *worse* -- a walk that starts inside the protein can never leave,
otherwise -- gives 30/60 and 51/60 accepted with no clashes left. It also
built its geometry from a MOL2 while turning atoms read from a **PDB**, whose
order need not match; both come from the same file now.

### In the documentation

| page | what was wrong |
|---|---|
| `plot_convolution_routines`, `plot_pile_up` | `IMP.bff.decay_fconv*` -- decay convolution is **tttrlib's** by the placement rule; they call `tttrlib.fconv`/`fconv_per`/`fconv_simd` now |
| `hgbp1_label_and_sample` | `{"N", "CA", "C"} <= set(site)`, and `attached[0]["resnum"]` |
| `langevin_hgbp1_site481` | the `None` timestep and `out_rmf` |
| `plot_k2_uncertainty` | `k2_call` returned two values where its caller unpacked three -- it had *never* matched; `Kappa2Distribution` supplies the third |
| `t4l_pmi` | fetched AVs for a display helper that went with `representation/av.py` |
| `structure_cgdye` | `attach_probes` tuples, `SITE_KEEP_ATOM_NAMES` (a function now), a nested transition-count matrix (flat + `n` now), and a `RotamerLibrary` read as a dict |
| `structure_accessible_volumes` | stopped on a missing `ipyvolume`, which is a viewer and not a dependency |

### What holds them now

`test/expensive_test_docs_and_examples.py` runs **every** example and **every**
notebook
-- 42 of them -- with each notebook's own directory as the working directory,
which is what Jupyter does and what their relative paths assume. Opt-in, like
the other `expensive_test_` files. `test/label/test_dye_commands.py` runs six
commands for real, and asserts the walk *accepts* something: a sampler that
accepts nothing passes any "does it run" check.

### One decision left

`doc/manual/decays/decay_curves.ipynb`, `decay_forward_model.ipynb`,
`decay_objective_function.ipynb` and `doc/manual/programming/programming_imp_decorator.ipynb`
document `IMP.bff.DecayCurve`, `DecayConvolution`, `DecayLifetimeHandler`,
`DecayScale`, `DecayPattern` and `DecayLinearization` -- 50 references to an
API **deliberately deleted** by `93b7198` ("PRD-113 stage 0: delete the TCSPC
instrument layer"), and which exists nowhere in the stack now, tttrlib
included. They are not broken call sites; they are pages describing a layer
that left. Porting them to tttrlib's manual or dropping them is a call for the
maintainer, so they are skipped by name in the sweep, with the reason written
where the skip is.

## 2026-08-26 (night, last +5) — two landmines defused

**`.drot` could not read a payload that compressed well.** `DrotReader`'s
decompressor sized its output buffer as `compressed * 8` and grew it on
`NEEDS_MORE_OUTPUT`; brotli's one-shot decoder reports a too-small buffer as
`BROTLI_DECODER_RESULT_ERROR`, so the growth branch never fired and a good file
came back as "corrupt brotli stream". Every shipped `.drot` is under the ratio,
which is why nothing had noticed -- but a library of near-identical conformers
(a rigid dye on a short linker) is not. It streams now.

**Four `attach_probes` calls in `bin/imp_bff` were dead.** Three passed
`(hierarchy, chain, residue)` tuples and one passed four parallel lists; the
C++ takes `ProbeAttachment` values, so all four raised `TypeError` at the
binding. `dye sample-rotamer` then read `attached[0]["site"]` out of a list of
values. The command runs end to end now, and `resolve_probe_site` -- which
answers CA, N, C -- is where the site comes from.

Neither had a test, which is why neither was known to be dead. Both do now:
`test/io/test_drot.py::test_a_payload_that_compresses_hugely_still_reads` and
`test/label/test_attach_probes_cli.py` (3 tests, including the command itself
through `CliRunner`).

## 2026-08-23 — PRD-118: the lossy ladder measured — and the 1–3 MB/dye budget is already met lossless

Per-stem accounting first: v8's largest dye+linker stem is **0.41 MB**
(median 175 KB), so the stated 1–3 MB/dye budget is met 3–7× over with zero
loss. The lossy option was still measured so it is quantified, not guessed:
grids coarsened through the real pipeline (0.25° → 4°) with the transition-
dipole error as the acceptance metric, against the guard that a 1.58°
dipole error is what broke the FRETpredict pins in the bcif investigation.
Verdict: 0.25° is safe (−10 %, 0.27° median), 0.5° spends most of the
margin (−20 %, 0.52°/1.11° med/p95), 1° and beyond sit in or past the
pin-breaking regime. The rungs are now real options, not code changes:
`encode_drot.py --grid-deg/--grid-a`, the grid recorded per member in
`drot.json`, and the reader honours it (verified end-to-end on a 0.5°
encode: −21.7 %, errors land on the ladder rung). Also counted: the
**loss-free** diet of shipping cutoff20+30 only is 2.13 MB — the first cut
if package size ever matters. Ladder + guard in
`prototypes/drot_rotlib/CONVERTING.md`.

# Update Log

## 2026-08-26 (night, last +7) — a use-after-free that returned zeros

Chasing the last of the untested surface -- the two `bin/` programs nothing
covered, and the top-level commands -- turned up two real defects and one that
had been silently wrong for as long as the CIF writer has existed.

### A system CIF written anywhere else could not be read back

The writer records each component's MOL2 path *relative to the CIF*, and
computed the `..` count from the **separators** of the base directory rather
than its segments -- one short whenever the base has no trailing separator,
which is always:

```
base   /var/folders/cl/xxxx/T/tmp1        6 segments
mol2   /Users/me/dye.mol2
wrote  ../../../../../Users/me/dye.mol2   5 ups
read   /var/Users/me/dye.mol2             does not exist
```

Nothing noticed because every test writes its systems into the data tree,
where the common prefix is long and the miscount cancels. `imp_bff simulate`
was dead for any system built outside it. `relative_path` in
`TopologyBuild.cpp` counts segments now, and
`test/cgdye/test_system_paths.py` writes a system into a `tmp_path` -- which
is exactly the case that failed.

### A container read off a temporary was freed memory

Writing the test for that found something worse:

```python
IMP.bff.read_forcefield_cif(path).components.values()
```

returned `'\x00\x00\x00...'` where a path should be, and segfaulted under
pytest. An accessor returning `const std::map&` gives SWIG a **borrowed**
pointer into the owner, and the proxy does not keep the owner alive: the
system is a temporary, it dies as the expression unwinds, and the map proxy
outlives it. Silently, most of the time.

The accessors cannot return by value -- `Scoring.cpp` calls `get_bonds()` in
the *condition* of a loop over the bonds, so a copy per access is quadratic --
so the copy is made **at the language boundary**: `%owned_container_out` in
`IMP_bff.types.i` gives every `const container&` return an owning copy, which
costs one copy per Python access and hands the lifetime to Python. Eight
container types: the maps, and the vectors whose element is one of this
module's own values.

**Not** `std::vector<std::string>`, `<double>` or `<int>`. Those already come
back *converted* -- a Python list, a numpy array -- which is a copy, so they
were never borrowed; claiming them replaced the conversion with a raw proxy and
five tests stopped recognising their own results (`assert <RMF_HDF5.Strings> ==
('SD',)`). The first cut of this fix did claim them, and the suite said so.
`%naturalvar` is on beside it, for the same hazard in plain member variables.

`test/test_borrowed_containers.py` reads five kinds of container off a
temporary, and checks that the copy really is one.

### Two programs nothing had ever run

`imp_bff_dye_pdb2cif` passed `--dye-id`'s `None` default into a
`const std::string&`, so its documented default invocation raised a
`TypeError` before reading anything. `imp_bff_labelizer --show`, given a
structure instead of a container, failed with "does not begin with an EBML
header" -- true, and useless. Both fixed, both now in
`test/io/test_bin_programs.py`.

`imp_bff dock` runs (score 188.30 over 99 distances, 17 accessible volumes),
and so do `build-system`, `simulate` and the six `dye` commands.

## 2026-08-26 (night, last +6) — nine dead paths, found by running the docs

The unit suite was green at every step of the last three passes, and the CLI
and the example gallery were quietly rotting the whole time: **what breaks when
a C++ port changes a shape is the code nobody runs**, and nobody runs the
documentation.

Running every `bin/imp_bff` command and every example and notebook found nine.

### In the program

| command | what was wrong |
|---|---|
| `dye sample-dof-walk` | imported `LinkerSampler`, a name from before the linker sampler became C++; read `resolve_probe_site` as a dict; and then, once running, **accepted nothing** |
| `dye sample-langevin` | passed `timestep_fs=None` into a `double`, and asked `run()` to write an RMF, which it stopped doing when the sampler became C++ |
| `dye label-fusion` | `sys.exit` in a module that never imported `sys` |
| `flexfit` (FP branch) | `attach_probes` with four parallel lists, and another dict-style site |

`sample-dof-walk` is the interesting one. It rejected every proposal because
the probe is *bonded into* the site: its first atoms sit a bond length from
residues i-1 and i+1, and the walk counted those bonds as clashes. Excluding
the site and its two neighbours, and accepting a move that does not make the
clash count *worse* -- a walk that starts inside the protein can never leave,
otherwise -- gives 30/60 and 51/60 accepted with no clashes left. It also
built its geometry from a MOL2 while turning atoms read from a **PDB**, whose
order need not match; both come from the same file now.

### In the documentation

| page | what was wrong |
|---|---|
| `plot_convolution_routines`, `plot_pile_up` | `IMP.bff.decay_fconv*` -- decay convolution is **tttrlib's** by the placement rule; they call `tttrlib.fconv`/`fconv_per`/`fconv_simd` now |
| `hgbp1_label_and_sample` | `{"N", "CA", "C"} <= set(site)`, and `attached[0]["resnum"]` |
| `langevin_hgbp1_site481` | the `None` timestep and `out_rmf` |
| `plot_k2_uncertainty` | `k2_call` returned two values where its caller unpacked three -- it had *never* matched; `Kappa2Distribution` supplies the third |
| `t4l_pmi` | fetched AVs for a display helper that went with `representation/av.py` |
| `structure_cgdye` | `attach_probes` tuples, `SITE_KEEP_ATOM_NAMES` (a function now), a nested transition-count matrix (flat + `n` now), and a `RotamerLibrary` read as a dict |
| `structure_accessible_volumes` | stopped on a missing `ipyvolume`, which is a viewer and not a dependency |

### What holds them now

`test/expensive_test_docs_and_examples.py` runs **every** example and **every**
notebook
-- 42 of them -- with each notebook's own directory as the working directory,
which is what Jupyter does and what their relative paths assume. Opt-in, like
the other `expensive_test_` files. `test/label/test_dye_commands.py` runs six
commands for real, and asserts the walk *accepts* something: a sampler that
accepts nothing passes any "does it run" check.

### One decision left

`doc/manual/decays/decay_curves.ipynb`, `decay_forward_model.ipynb`,
`decay_objective_function.ipynb` and `doc/manual/programming/programming_imp_decorator.ipynb`
document `IMP.bff.DecayCurve`, `DecayConvolution`, `DecayLifetimeHandler`,
`DecayScale`, `DecayPattern` and `DecayLinearization` -- 50 references to an
API **deliberately deleted** by `93b7198` ("PRD-113 stage 0: delete the TCSPC
instrument layer"), and which exists nowhere in the stack now, tttrlib
included. They are not broken call sites; they are pages describing a layer
that left. Porting them to tttrlib's manual or dropping them is a call for the
maintainer, so they are skipped by name in the sweep, with the reason written
where the skip is.

## 2026-08-26 (night, last +5) — two landmines defused

**`.drot` could not read a payload that compressed well.** `DrotReader`'s
decompressor sized its output buffer as `compressed * 8` and grew it on
`NEEDS_MORE_OUTPUT`; brotli's one-shot decoder reports a too-small buffer as
`BROTLI_DECODER_RESULT_ERROR`, so the growth branch never fired and a good file
came back as "corrupt brotli stream". Every shipped `.drot` is under the ratio,
which is why nothing had noticed -- but a library of near-identical conformers
(a rigid dye on a short linker) is not. It streams now.

**Four `attach_probes` calls in `bin/imp_bff` were dead.** Three passed
`(hierarchy, chain, residue)` tuples and one passed four parallel lists; the
C++ takes `ProbeAttachment` values, so all four raised `TypeError` at the
binding. `dye sample-rotamer` then read `attached[0]["site"]` out of a list of
values. The command runs end to end now, and `resolve_probe_site` -- which
answers CA, N, C -- is where the site comes from.

Neither had a test, which is why neither was known to be dead. Both do now:
`test/io/test_drot.py::test_a_payload_that_compresses_hugely_still_reads` and
`test/label/test_attach_probes_cli.py` (3 tests, including the command itself
through `CliRunner`).

## 2026-08-22 — PRD-117 batch 17: rotamer site kernels + library registry to C++ (`RotamerSite.h`)

New header `include/RotamerSite.h` / `src/RotamerSite.cpp` carrying what was
Python in `representation/rotamer.py`:

- **The backbone frame**: `resolve_backbone_site` (CA/N/C by chain+residue,
  raising which atom is missing), `backbone_rotation` (x along CA->N, y in the
  N-CA-C plane, z = x cross y), `transform_library_to_site` -- nine
  multiplications a caller should not re-type, now out-views.
- **`selector_atom_indices`**: the FRETpredict selector index resolution, same
  `NAME and resname RES` grammar and ambiguity rules as scoring's -- the
  resname clause is honoured when resnames are given (atom names repeat
  between dye and linker residues), and a miss raises naming the selector.
- **The library registry**: `rotamer_library_metadata` (libraries.json read
  once in C++, entry returned as JSON text with name/library_name/cutoff
  added), `normalize_library_name` (cutoff-suffix grammar),
  `library_name_cutoff`, `library_filename`, `resolve_rotamer_library_path`
  (canonical FRETpredict set first so the *requested cutoff* loads; the
  only-cutoff-30-RMF mismatch raises).

`rotamer.i`'s Python functions with the same names became thin shape/bridge
wrappers (1420 -> 1252 lines); the fps payload dataclasses,
`load_rotamer_library` (bcif via C++ `load_rotamer_library_dcd`, RMF door
lazy), `load_protein_frames` (IMP.atom/IMP.rmf) and the `RotamerFRET` driver
stay Python -- RMF and IMP API glue. Gotcha repeated and fixed: a %pythoncode
def that shadows its own SWIG name must call `_IMP_bff.<name>`, or it
recurses into itself (all eight bridges).

Verification: `ninja IMP.bff` clean; non-medium **624 passed, 3 xfailed**
(only the pre-existing `av_reference_0.mrc` data failure); cgdye **160
passed**; expensive all-dyes **36 passed**.

## 2026-08-23 — PRD-118: `.drot` v8 — column-major grids + brotli, 5.31× vs `.bcif`

"Better encoding?" answered with the full codec matrix on the real payload
(95 libraries): zstd-22-long +11.5 % vs xz, PPMd worse, per-row centering
noise, xz-DELTA +9 % — and **brotli −6.6 %** on identical tar bytes with
**3.6× faster decode** (libbrotlidec already in the runtime env). Combined
with the other measured win — **column-major grids (−7.4 %: same-row values
cluster; contiguity is what the context model needs)** — v8 lands at
**30.05 MB `.bcif` → 5.66 MB `.drot` = 5.31× smaller (81 % saved)**, v5's
6.59 MB cut by another 14 %. Validation unchanged (93/95 exact, worst
0.018 Å; the two known template-defect libs), load + full reconstruction
of every conformer ~9 s, every file unpacks `brotli -dc | tar -t`
(GNU: `tar -I brotli -xf`). One **corpus archive** was measured too
(4.52 MB, 6.65×) and rejected as ship shape — couples all libraries,
kills per-file lazy loading. Full decision ledger (nine candidates,
eight rejected with numbers) in `prototypes/drot_rotlib/CONVERTING.md`;
reader still dispatches v2–v8.

# Update Log

## 2026-08-26 (night, last +7) — a use-after-free that returned zeros

Chasing the last of the untested surface -- the two `bin/` programs nothing
covered, and the top-level commands -- turned up two real defects and one that
had been silently wrong for as long as the CIF writer has existed.

### A system CIF written anywhere else could not be read back

The writer records each component's MOL2 path *relative to the CIF*, and
computed the `..` count from the **separators** of the base directory rather
than its segments -- one short whenever the base has no trailing separator,
which is always:

```
base   /var/folders/cl/xxxx/T/tmp1        6 segments
mol2   /Users/me/dye.mol2
wrote  ../../../../../Users/me/dye.mol2   5 ups
read   /var/Users/me/dye.mol2             does not exist
```

Nothing noticed because every test writes its systems into the data tree,
where the common prefix is long and the miscount cancels. `imp_bff simulate`
was dead for any system built outside it. `relative_path` in
`TopologyBuild.cpp` counts segments now, and
`test/cgdye/test_system_paths.py` writes a system into a `tmp_path` -- which
is exactly the case that failed.

### A container read off a temporary was freed memory

Writing the test for that found something worse:

```python
IMP.bff.read_forcefield_cif(path).components.values()
```

returned `'\x00\x00\x00...'` where a path should be, and segfaulted under
pytest. An accessor returning `const std::map&` gives SWIG a **borrowed**
pointer into the owner, and the proxy does not keep the owner alive: the
system is a temporary, it dies as the expression unwinds, and the map proxy
outlives it. Silently, most of the time.

The accessors cannot return by value -- `Scoring.cpp` calls `get_bonds()` in
the *condition* of a loop over the bonds, so a copy per access is quadratic --
so the copy is made **at the language boundary**: `%owned_container_out` in
`IMP_bff.types.i` gives every `const container&` return an owning copy, which
costs one copy per Python access and hands the lifetime to Python. Eight
container types: the maps, and the vectors whose element is one of this
module's own values.

**Not** `std::vector<std::string>`, `<double>` or `<int>`. Those already come
back *converted* -- a Python list, a numpy array -- which is a copy, so they
were never borrowed; claiming them replaced the conversion with a raw proxy and
five tests stopped recognising their own results (`assert <RMF_HDF5.Strings> ==
('SD',)`). The first cut of this fix did claim them, and the suite said so.
`%naturalvar` is on beside it, for the same hazard in plain member variables.

`test/test_borrowed_containers.py` reads five kinds of container off a
temporary, and checks that the copy really is one.

### Two programs nothing had ever run

`imp_bff_dye_pdb2cif` passed `--dye-id`'s `None` default into a
`const std::string&`, so its documented default invocation raised a
`TypeError` before reading anything. `imp_bff_labelizer --show`, given a
structure instead of a container, failed with "does not begin with an EBML
header" -- true, and useless. Both fixed, both now in
`test/io/test_bin_programs.py`.

`imp_bff dock` runs (score 188.30 over 99 distances, 17 accessible volumes),
and so do `build-system`, `simulate` and the six `dye` commands.

## 2026-08-26 (night, last +6) — nine dead paths, found by running the docs

The unit suite was green at every step of the last three passes, and the CLI
and the example gallery were quietly rotting the whole time: **what breaks when
a C++ port changes a shape is the code nobody runs**, and nobody runs the
documentation.

Running every `bin/imp_bff` command and every example and notebook found nine.

### In the program

| command | what was wrong |
|---|---|
| `dye sample-dof-walk` | imported `LinkerSampler`, a name from before the linker sampler became C++; read `resolve_probe_site` as a dict; and then, once running, **accepted nothing** |
| `dye sample-langevin` | passed `timestep_fs=None` into a `double`, and asked `run()` to write an RMF, which it stopped doing when the sampler became C++ |
| `dye label-fusion` | `sys.exit` in a module that never imported `sys` |
| `flexfit` (FP branch) | `attach_probes` with four parallel lists, and another dict-style site |

`sample-dof-walk` is the interesting one. It rejected every proposal because
the probe is *bonded into* the site: its first atoms sit a bond length from
residues i-1 and i+1, and the walk counted those bonds as clashes. Excluding
the site and its two neighbours, and accepting a move that does not make the
clash count *worse* -- a walk that starts inside the protein can never leave,
otherwise -- gives 30/60 and 51/60 accepted with no clashes left. It also
built its geometry from a MOL2 while turning atoms read from a **PDB**, whose
order need not match; both come from the same file now.

### In the documentation

| page | what was wrong |
|---|---|
| `plot_convolution_routines`, `plot_pile_up` | `IMP.bff.decay_fconv*` -- decay convolution is **tttrlib's** by the placement rule; they call `tttrlib.fconv`/`fconv_per`/`fconv_simd` now |
| `hgbp1_label_and_sample` | `{"N", "CA", "C"} <= set(site)`, and `attached[0]["resnum"]` |
| `langevin_hgbp1_site481` | the `None` timestep and `out_rmf` |
| `plot_k2_uncertainty` | `k2_call` returned two values where its caller unpacked three -- it had *never* matched; `Kappa2Distribution` supplies the third |
| `t4l_pmi` | fetched AVs for a display helper that went with `representation/av.py` |
| `structure_cgdye` | `attach_probes` tuples, `SITE_KEEP_ATOM_NAMES` (a function now), a nested transition-count matrix (flat + `n` now), and a `RotamerLibrary` read as a dict |
| `structure_accessible_volumes` | stopped on a missing `ipyvolume`, which is a viewer and not a dependency |

### What holds them now

`test/expensive_test_docs_and_examples.py` runs **every** example and **every**
notebook
-- 42 of them -- with each notebook's own directory as the working directory,
which is what Jupyter does and what their relative paths assume. Opt-in, like
the other `expensive_test_` files. `test/label/test_dye_commands.py` runs six
commands for real, and asserts the walk *accepts* something: a sampler that
accepts nothing passes any "does it run" check.

### One decision left

`doc/manual/decays/decay_curves.ipynb`, `decay_forward_model.ipynb`,
`decay_objective_function.ipynb` and `doc/manual/programming/programming_imp_decorator.ipynb`
document `IMP.bff.DecayCurve`, `DecayConvolution`, `DecayLifetimeHandler`,
`DecayScale`, `DecayPattern` and `DecayLinearization` -- 50 references to an
API **deliberately deleted** by `93b7198` ("PRD-113 stage 0: delete the TCSPC
instrument layer"), and which exists nowhere in the stack now, tttrlib
included. They are not broken call sites; they are pages describing a layer
that left. Porting them to tttrlib's manual or dropping them is a call for the
maintainer, so they are skipped by name in the sweep, with the reason written
where the skip is.

## 2026-08-26 (night, last +5) — two landmines defused

**`.drot` could not read a payload that compressed well.** `DrotReader`'s
decompressor sized its output buffer as `compressed * 8` and grew it on
`NEEDS_MORE_OUTPUT`; brotli's one-shot decoder reports a too-small buffer as
`BROTLI_DECODER_RESULT_ERROR`, so the growth branch never fired and a good file
came back as "corrupt brotli stream". Every shipped `.drot` is under the ratio,
which is why nothing had noticed -- but a library of near-identical conformers
(a rigid dye on a short linker) is not. It streams now.

**Four `attach_probes` calls in `bin/imp_bff` were dead.** Three passed
`(hierarchy, chain, residue)` tuples and one passed four parallel lists; the
C++ takes `ProbeAttachment` values, so all four raised `TypeError` at the
binding. `dye sample-rotamer` then read `attached[0]["site"]` out of a list of
values. The command runs end to end now, and `resolve_probe_site` -- which
answers CA, N, C -- is where the site comes from.

Neither had a test, which is why neither was known to be dead. Both do now:
`test/io/test_drot.py::test_a_payload_that_compresses_hugely_still_reads` and
`test/label/test_attach_probes_cli.py` (3 tests, including the command itself
through `CliRunner`).

## 2026-08-23 — PRD-118: `.drot` v6 — the data members become JSON

Final payload shape (user direction): the tar.xz stays, and every data
member is now JSON text — `base.json` / `theta.json` / `phi.json` (flat int
arrays, the same int16 grids spelled out) and `weights.json` (int array;
weights are cluster populations), alongside `drot.json`, `template.cif`
(CIF per the earlier explicit call) and `rows.json`. Unpack with `tar -xf`
and everything but the template parses in any JSON parser — no bespoke
binary decoding on the data path. **Measured cost: +15.6 % over v5's int16
members after LZMA — 7.68 vs 6.59 MB across the 95 libraries — still 3.91×
smaller than the shipped `.bcif` (74 % saved).** Validation unchanged:
93/95 reconstruct exactly (worst 0.018 Å; the two known template-defect
libs), all 95 files pass `xz -t` and `tar -tf`, JSON-member load of every
conformer ~7 s. Spec and full format-decision ledger (bit-pack rejected,
msgpack rejected, JSON accepted at its measured price) in
`prototypes/drot_rotlib/CONVERTING.md`.

# Update Log

## 2026-08-26 (night, last +7) — a use-after-free that returned zeros

Chasing the last of the untested surface -- the two `bin/` programs nothing
covered, and the top-level commands -- turned up two real defects and one that
had been silently wrong for as long as the CIF writer has existed.

### A system CIF written anywhere else could not be read back

The writer records each component's MOL2 path *relative to the CIF*, and
computed the `..` count from the **separators** of the base directory rather
than its segments -- one short whenever the base has no trailing separator,
which is always:

```
base   /var/folders/cl/xxxx/T/tmp1        6 segments
mol2   /Users/me/dye.mol2
wrote  ../../../../../Users/me/dye.mol2   5 ups
read   /var/Users/me/dye.mol2             does not exist
```

Nothing noticed because every test writes its systems into the data tree,
where the common prefix is long and the miscount cancels. `imp_bff simulate`
was dead for any system built outside it. `relative_path` in
`TopologyBuild.cpp` counts segments now, and
`test/cgdye/test_system_paths.py` writes a system into a `tmp_path` -- which
is exactly the case that failed.

### A container read off a temporary was freed memory

Writing the test for that found something worse:

```python
IMP.bff.read_forcefield_cif(path).components.values()
```

returned `'\x00\x00\x00...'` where a path should be, and segfaulted under
pytest. An accessor returning `const std::map&` gives SWIG a **borrowed**
pointer into the owner, and the proxy does not keep the owner alive: the
system is a temporary, it dies as the expression unwinds, and the map proxy
outlives it. Silently, most of the time.

The accessors cannot return by value -- `Scoring.cpp` calls `get_bonds()` in
the *condition* of a loop over the bonds, so a copy per access is quadratic --
so the copy is made **at the language boundary**: `%owned_container_out` in
`IMP_bff.types.i` gives every `const container&` return an owning copy, which
costs one copy per Python access and hands the lifetime to Python. Eight
container types: the maps, and the vectors whose element is one of this
module's own values.

**Not** `std::vector<std::string>`, `<double>` or `<int>`. Those already come
back *converted* -- a Python list, a numpy array -- which is a copy, so they
were never borrowed; claiming them replaced the conversion with a raw proxy and
five tests stopped recognising their own results (`assert <RMF_HDF5.Strings> ==
('SD',)`). The first cut of this fix did claim them, and the suite said so.
`%naturalvar` is on beside it, for the same hazard in plain member variables.

`test/test_borrowed_containers.py` reads five kinds of container off a
temporary, and checks that the copy really is one.

### Two programs nothing had ever run

`imp_bff_dye_pdb2cif` passed `--dye-id`'s `None` default into a
`const std::string&`, so its documented default invocation raised a
`TypeError` before reading anything. `imp_bff_labelizer --show`, given a
structure instead of a container, failed with "does not begin with an EBML
header" -- true, and useless. Both fixed, both now in
`test/io/test_bin_programs.py`.

`imp_bff dock` runs (score 188.30 over 99 distances, 17 accessible volumes),
and so do `build-system`, `simulate` and the six `dye` commands.

## 2026-08-26 (night, last +6) — nine dead paths, found by running the docs

The unit suite was green at every step of the last three passes, and the CLI
and the example gallery were quietly rotting the whole time: **what breaks when
a C++ port changes a shape is the code nobody runs**, and nobody runs the
documentation.

Running every `bin/imp_bff` command and every example and notebook found nine.

### In the program

| command | what was wrong |
|---|---|
| `dye sample-dof-walk` | imported `LinkerSampler`, a name from before the linker sampler became C++; read `resolve_probe_site` as a dict; and then, once running, **accepted nothing** |
| `dye sample-langevin` | passed `timestep_fs=None` into a `double`, and asked `run()` to write an RMF, which it stopped doing when the sampler became C++ |
| `dye label-fusion` | `sys.exit` in a module that never imported `sys` |
| `flexfit` (FP branch) | `attach_probes` with four parallel lists, and another dict-style site |

`sample-dof-walk` is the interesting one. It rejected every proposal because
the probe is *bonded into* the site: its first atoms sit a bond length from
residues i-1 and i+1, and the walk counted those bonds as clashes. Excluding
the site and its two neighbours, and accepting a move that does not make the
clash count *worse* -- a walk that starts inside the protein can never leave,
otherwise -- gives 30/60 and 51/60 accepted with no clashes left. It also
built its geometry from a MOL2 while turning atoms read from a **PDB**, whose
order need not match; both come from the same file now.

### In the documentation

| page | what was wrong |
|---|---|
| `plot_convolution_routines`, `plot_pile_up` | `IMP.bff.decay_fconv*` -- decay convolution is **tttrlib's** by the placement rule; they call `tttrlib.fconv`/`fconv_per`/`fconv_simd` now |
| `hgbp1_label_and_sample` | `{"N", "CA", "C"} <= set(site)`, and `attached[0]["resnum"]` |
| `langevin_hgbp1_site481` | the `None` timestep and `out_rmf` |
| `plot_k2_uncertainty` | `k2_call` returned two values where its caller unpacked three -- it had *never* matched; `Kappa2Distribution` supplies the third |
| `t4l_pmi` | fetched AVs for a display helper that went with `representation/av.py` |
| `structure_cgdye` | `attach_probes` tuples, `SITE_KEEP_ATOM_NAMES` (a function now), a nested transition-count matrix (flat + `n` now), and a `RotamerLibrary` read as a dict |
| `structure_accessible_volumes` | stopped on a missing `ipyvolume`, which is a viewer and not a dependency |

### What holds them now

`test/expensive_test_docs_and_examples.py` runs **every** example and **every**
notebook
-- 42 of them -- with each notebook's own directory as the working directory,
which is what Jupyter does and what their relative paths assume. Opt-in, like
the other `expensive_test_` files. `test/label/test_dye_commands.py` runs six
commands for real, and asserts the walk *accepts* something: a sampler that
accepts nothing passes any "does it run" check.

### One decision left

`doc/manual/decays/decay_curves.ipynb`, `decay_forward_model.ipynb`,
`decay_objective_function.ipynb` and `doc/manual/programming/programming_imp_decorator.ipynb`
document `IMP.bff.DecayCurve`, `DecayConvolution`, `DecayLifetimeHandler`,
`DecayScale`, `DecayPattern` and `DecayLinearization` -- 50 references to an
API **deliberately deleted** by `93b7198` ("PRD-113 stage 0: delete the TCSPC
instrument layer"), and which exists nowhere in the stack now, tttrlib
included. They are not broken call sites; they are pages describing a layer
that left. Porting them to tttrlib's manual or dropping them is a call for the
maintainer, so they are skipped by name in the sweep, with the reason written
where the skip is.

## 2026-08-26 (night, last +5) — two landmines defused

**`.drot` could not read a payload that compressed well.** `DrotReader`'s
decompressor sized its output buffer as `compressed * 8` and grew it on
`NEEDS_MORE_OUTPUT`; brotli's one-shot decoder reports a too-small buffer as
`BROTLI_DECODER_RESULT_ERROR`, so the growth branch never fired and a good file
came back as "corrupt brotli stream". Every shipped `.drot` is under the ratio,
which is why nothing had noticed -- but a library of near-identical conformers
(a rigid dye on a short linker) is not. It streams now.

**Four `attach_probes` calls in `bin/imp_bff` were dead.** Three passed
`(hierarchy, chain, residue)` tuples and one passed four parallel lists; the
C++ takes `ProbeAttachment` values, so all four raised `TypeError` at the
binding. `dye sample-rotamer` then read `attached[0]["site"]` out of a list of
values. The command runs end to end now, and `resolve_probe_site` -- which
answers CA, N, C -- is where the site comes from.

Neither had a test, which is why neither was known to be dead. Both do now:
`test/io/test_drot.py::test_a_payload_that_compresses_hugely_still_reads` and
`test/label/test_attach_probes_cli.py` (3 tests, including the command itself
through `CliRunner`).

## 2026-08-22 — PRD-117 batch 16: `scoring.i` to C++ — the whole stage-2 layer

667 lines of `%pythoncode` became 150 (the sanctioned remainder). Scoring.h
grew from parameter-table shims to the full stage-2 layer:

- `charmm36_lj`/`lj_cross` are public and return 2-value managed views
  (un-`%ignore`d; `std::array` -> the out-view pattern the repo already uses).
  New beside them: `atom_type`, `lj_score`, `scaled_parameters`,
  `cross_lj_params`, `lj_pairs_sum`, `pair_energy_matrix` (over the C++
  conformer-pair kernel), and `LJArrays` values.
- The **FRETpredict selector mini-language** is C++: `selector_matches`,
  `selector_resnames`, `selector_atom_names`, `hydrogen_mask`, `site_mask`,
  `protein_charge_mask`, `rotamer_charge_mask`.
- `BoundingBoxFilter` is a C++ class (`build`/`build_single`/`intersects`/
  `intersects_reference`/`filter_frames` -> `AABBFilterResult{coords, mask}`).
- `compute_rotamer_score` is **end-to-end** now: it takes the names/resnames/
  selectors the Python orchestration took, builds the masks itself, assembles
  the Lorentz-Berthelot parameters, calls the all-pairs kernel. The old
  pre-assembled-parameters entry point is gone (it had no caller).
- `rotamer_mean_field_weights` (single dye) iterates in C++ over the
  pair-energy matrix; `boltzmann_weights`/`rotamer_cluster_weights` already
  were C++ and now surface as such (they return tuples/lists -- callers wrap
  with `np.asarray`).
- Typed-system walkers: `site_element_map`, `compute_lj_pair_sites` (returns
  `LJSitePair` values), `build_lj_type_table` (map of `FFLJType`), and
  `DyeInternalEnergyEvaluator` as a C++ class taking the **typed system**
  (evaluate/evaluate_batch(flat, n, m)/evaluate_batch_filtered ->
  `EnergyMaskResult{energies, mask}`).

Still Python, deliberately: `rotamer_mean_field_weights_multi_dye` (numpy
matrix algebra over the C++ kernel), `dye_internal_system` (bridges
`parse_dye_mol2` dicts, moves with the topology batch), and
`torsion_cosine`/`build_dye_restraints` (create IMP.core objects).

Also: swig.i-in runs `forcefield.i` before `scoring.i` -- Scoring.h's
`DyeForceFieldSystem` parameters need DyeForceField.h parsed first or the
wrappers emit unqualified names and fail to compile; SWIG does not follow a
header's #include chain. rotamer_ensemble's from_site passes lists (a numpy
`or []` raised an ambiguous-truth error). topology.i's LJ-type table is
`build_lj_type_table` now (CHARMM36_LJ is gone with its dict). The
`_atom_type` lru_cache test became a plain answers test -- C++ needs no
memoisation.

Migrated tests: scoring/*, cgdye/test_dye_topology, test_rotamer_generation,
test_physics_invariants, dye/test_library_cache; sampling.i and
rotamer_ensemble.i call sites.

Verification: `ninja IMP.bff` clean; non-medium **624 passed, 3 xfailed**
(only the pre-existing `av_reference_0.mrc` data failure); cgdye **160
passed**; expensive all-dyes **36 passed**.

## 2026-08-23 — PRD-118: "why not msgpack?" — measured

Answered with the corpus rather than an argument: one msgpack document
(grids as `bin`, header as a map) through the same xz stream is **+0.13 %
larger** than v5's tar.xz (6.597 vs 6.588 MB over the 95 libraries) — tar's
repetitive headers compress to less than msgpack's framing saves. Two
presumptions corrected on the way: msgpack **is** already in the Python dep
tree (ships with python-ihm) and its C decoder **is** already compiled into
IMP (`ihm_format.c` includes `cmp.h`, and bff's own `read_bcif_trajectory`
drives it in binary mode) — so neither side lacks msgpack; the decision is
the size tie plus the requirement that the file unpack with zero tooling.
Recorded next to the bit-packing rejection in
`prototypes/drot_rotlib/CONVERTING.md`.

# Update Log

## 2026-08-26 (night, last +7) — a use-after-free that returned zeros

Chasing the last of the untested surface -- the two `bin/` programs nothing
covered, and the top-level commands -- turned up two real defects and one that
had been silently wrong for as long as the CIF writer has existed.

### A system CIF written anywhere else could not be read back

The writer records each component's MOL2 path *relative to the CIF*, and
computed the `..` count from the **separators** of the base directory rather
than its segments -- one short whenever the base has no trailing separator,
which is always:

```
base   /var/folders/cl/xxxx/T/tmp1        6 segments
mol2   /Users/me/dye.mol2
wrote  ../../../../../Users/me/dye.mol2   5 ups
read   /var/Users/me/dye.mol2             does not exist
```

Nothing noticed because every test writes its systems into the data tree,
where the common prefix is long and the miscount cancels. `imp_bff simulate`
was dead for any system built outside it. `relative_path` in
`TopologyBuild.cpp` counts segments now, and
`test/cgdye/test_system_paths.py` writes a system into a `tmp_path` -- which
is exactly the case that failed.

### A container read off a temporary was freed memory

Writing the test for that found something worse:

```python
IMP.bff.read_forcefield_cif(path).components.values()
```

returned `'\x00\x00\x00...'` where a path should be, and segfaulted under
pytest. An accessor returning `const std::map&` gives SWIG a **borrowed**
pointer into the owner, and the proxy does not keep the owner alive: the
system is a temporary, it dies as the expression unwinds, and the map proxy
outlives it. Silently, most of the time.

The accessors cannot return by value -- `Scoring.cpp` calls `get_bonds()` in
the *condition* of a loop over the bonds, so a copy per access is quadratic --
so the copy is made **at the language boundary**: `%owned_container_out` in
`IMP_bff.types.i` gives every `const container&` return an owning copy, which
costs one copy per Python access and hands the lifetime to Python. Eight
container types: the maps, and the vectors whose element is one of this
module's own values.

**Not** `std::vector<std::string>`, `<double>` or `<int>`. Those already come
back *converted* -- a Python list, a numpy array -- which is a copy, so they
were never borrowed; claiming them replaced the conversion with a raw proxy and
five tests stopped recognising their own results (`assert <RMF_HDF5.Strings> ==
('SD',)`). The first cut of this fix did claim them, and the suite said so.
`%naturalvar` is on beside it, for the same hazard in plain member variables.

`test/test_borrowed_containers.py` reads five kinds of container off a
temporary, and checks that the copy really is one.

### Two programs nothing had ever run

`imp_bff_dye_pdb2cif` passed `--dye-id`'s `None` default into a
`const std::string&`, so its documented default invocation raised a
`TypeError` before reading anything. `imp_bff_labelizer --show`, given a
structure instead of a container, failed with "does not begin with an EBML
header" -- true, and useless. Both fixed, both now in
`test/io/test_bin_programs.py`.

`imp_bff dock` runs (score 188.30 over 99 distances, 17 accessible volumes),
and so do `build-system`, `simulate` and the six `dye` commands.

## 2026-08-26 (night, last +6) — nine dead paths, found by running the docs

The unit suite was green at every step of the last three passes, and the CLI
and the example gallery were quietly rotting the whole time: **what breaks when
a C++ port changes a shape is the code nobody runs**, and nobody runs the
documentation.

Running every `bin/imp_bff` command and every example and notebook found nine.

### In the program

| command | what was wrong |
|---|---|
| `dye sample-dof-walk` | imported `LinkerSampler`, a name from before the linker sampler became C++; read `resolve_probe_site` as a dict; and then, once running, **accepted nothing** |
| `dye sample-langevin` | passed `timestep_fs=None` into a `double`, and asked `run()` to write an RMF, which it stopped doing when the sampler became C++ |
| `dye label-fusion` | `sys.exit` in a module that never imported `sys` |
| `flexfit` (FP branch) | `attach_probes` with four parallel lists, and another dict-style site |

`sample-dof-walk` is the interesting one. It rejected every proposal because
the probe is *bonded into* the site: its first atoms sit a bond length from
residues i-1 and i+1, and the walk counted those bonds as clashes. Excluding
the site and its two neighbours, and accepting a move that does not make the
clash count *worse* -- a walk that starts inside the protein can never leave,
otherwise -- gives 30/60 and 51/60 accepted with no clashes left. It also
built its geometry from a MOL2 while turning atoms read from a **PDB**, whose
order need not match; both come from the same file now.

### In the documentation

| page | what was wrong |
|---|---|
| `plot_convolution_routines`, `plot_pile_up` | `IMP.bff.decay_fconv*` -- decay convolution is **tttrlib's** by the placement rule; they call `tttrlib.fconv`/`fconv_per`/`fconv_simd` now |
| `hgbp1_label_and_sample` | `{"N", "CA", "C"} <= set(site)`, and `attached[0]["resnum"]` |
| `langevin_hgbp1_site481` | the `None` timestep and `out_rmf` |
| `plot_k2_uncertainty` | `k2_call` returned two values where its caller unpacked three -- it had *never* matched; `Kappa2Distribution` supplies the third |
| `t4l_pmi` | fetched AVs for a display helper that went with `representation/av.py` |
| `structure_cgdye` | `attach_probes` tuples, `SITE_KEEP_ATOM_NAMES` (a function now), a nested transition-count matrix (flat + `n` now), and a `RotamerLibrary` read as a dict |
| `structure_accessible_volumes` | stopped on a missing `ipyvolume`, which is a viewer and not a dependency |

### What holds them now

`test/expensive_test_docs_and_examples.py` runs **every** example and **every**
notebook
-- 42 of them -- with each notebook's own directory as the working directory,
which is what Jupyter does and what their relative paths assume. Opt-in, like
the other `expensive_test_` files. `test/label/test_dye_commands.py` runs six
commands for real, and asserts the walk *accepts* something: a sampler that
accepts nothing passes any "does it run" check.

### One decision left

`doc/manual/decays/decay_curves.ipynb`, `decay_forward_model.ipynb`,
`decay_objective_function.ipynb` and `doc/manual/programming/programming_imp_decorator.ipynb`
document `IMP.bff.DecayCurve`, `DecayConvolution`, `DecayLifetimeHandler`,
`DecayScale`, `DecayPattern` and `DecayLinearization` -- 50 references to an
API **deliberately deleted** by `93b7198` ("PRD-113 stage 0: delete the TCSPC
instrument layer"), and which exists nowhere in the stack now, tttrlib
included. They are not broken call sites; they are pages describing a layer
that left. Porting them to tttrlib's manual or dropping them is a call for the
maintainer, so they are skipped by name in the sweep, with the reason written
where the skip is.

## 2026-08-26 (night, last +5) — two landmines defused

**`.drot` could not read a payload that compressed well.** `DrotReader`'s
decompressor sized its output buffer as `compressed * 8` and grew it on
`NEEDS_MORE_OUTPUT`; brotli's one-shot decoder reports a too-small buffer as
`BROTLI_DECODER_RESULT_ERROR`, so the growth branch never fired and a good file
came back as "corrupt brotli stream". Every shipped `.drot` is under the ratio,
which is why nothing had noticed -- but a library of near-identical conformers
(a rigid dye on a short linker) is not. It streams now.

**Four `attach_probes` calls in `bin/imp_bff` were dead.** Three passed
`(hierarchy, chain, residue)` tuples and one passed four parallel lists; the
C++ takes `ProbeAttachment` values, so all four raised `TypeError` at the
binding. `dye sample-rotamer` then read `attached[0]["site"]` out of a list of
values. The command runs end to end now, and `resolve_probe_site` -- which
answers CA, N, C -- is where the site comes from.

Neither had a test, which is why neither was known to be dead. Both do now:
`test/io/test_drot.py::test_a_payload_that_compresses_hugely_still_reads` and
`test/label/test_attach_probes_cli.py` (3 tests, including the command itself
through `CliRunner`).

## 2026-08-23 — PRD-118: `.drot` v5 — the container becomes an unpackable archive

Final container shape: a `.drot` is now a **`tar.xz`** — `tar -tf` lists the
members, `tar -xf` unpacks them (`xz -t` still verifies). The members are
the library's parts as named files: `drot.json` (header), `template.cif`
(the mmCIF reference geometry), `rows.json` (the Z-matrix table),
`base.i16`/`theta.i16`/`phi.i16` (the grids), `weights.bin`. No bespoke
container code anywhere: Python `tarfile`, the CLIs, and C++ liblzma +
tar's 512-byte headers all read it; the loader still rebuilds and
cross-checks the `ZTree` from the embedded template. One tar payload
through one LZMA stream keeps the joint compression — the tar layering
costs ~0.4 KB/file (~40 KB over the corpus). **Sweep: 30.05 MB `.bcif` →
6.59 MB `.drot` = 4.56× smaller (78 % saved)**, 93/95 reconstruct exactly
(worst 0.018 Å), all 95 files pass `xz -t` *and* `tar -tf`, self-contained
load of all conformers ~8 s. Spec and commands in
`prototypes/drot_rotlib/CONVERTING.md`.

# Update Log

## 2026-08-26 (night, last +7) — a use-after-free that returned zeros

Chasing the last of the untested surface -- the two `bin/` programs nothing
covered, and the top-level commands -- turned up two real defects and one that
had been silently wrong for as long as the CIF writer has existed.

### A system CIF written anywhere else could not be read back

The writer records each component's MOL2 path *relative to the CIF*, and
computed the `..` count from the **separators** of the base directory rather
than its segments -- one short whenever the base has no trailing separator,
which is always:

```
base   /var/folders/cl/xxxx/T/tmp1        6 segments
mol2   /Users/me/dye.mol2
wrote  ../../../../../Users/me/dye.mol2   5 ups
read   /var/Users/me/dye.mol2             does not exist
```

Nothing noticed because every test writes its systems into the data tree,
where the common prefix is long and the miscount cancels. `imp_bff simulate`
was dead for any system built outside it. `relative_path` in
`TopologyBuild.cpp` counts segments now, and
`test/cgdye/test_system_paths.py` writes a system into a `tmp_path` -- which
is exactly the case that failed.

### A container read off a temporary was freed memory

Writing the test for that found something worse:

```python
IMP.bff.read_forcefield_cif(path).components.values()
```

returned `'\x00\x00\x00...'` where a path should be, and segfaulted under
pytest. An accessor returning `const std::map&` gives SWIG a **borrowed**
pointer into the owner, and the proxy does not keep the owner alive: the
system is a temporary, it dies as the expression unwinds, and the map proxy
outlives it. Silently, most of the time.

The accessors cannot return by value -- `Scoring.cpp` calls `get_bonds()` in
the *condition* of a loop over the bonds, so a copy per access is quadratic --
so the copy is made **at the language boundary**: `%owned_container_out` in
`IMP_bff.types.i` gives every `const container&` return an owning copy, which
costs one copy per Python access and hands the lifetime to Python. Eight
container types: the maps, and the vectors whose element is one of this
module's own values.

**Not** `std::vector<std::string>`, `<double>` or `<int>`. Those already come
back *converted* -- a Python list, a numpy array -- which is a copy, so they
were never borrowed; claiming them replaced the conversion with a raw proxy and
five tests stopped recognising their own results (`assert <RMF_HDF5.Strings> ==
('SD',)`). The first cut of this fix did claim them, and the suite said so.
`%naturalvar` is on beside it, for the same hazard in plain member variables.

`test/test_borrowed_containers.py` reads five kinds of container off a
temporary, and checks that the copy really is one.

### Two programs nothing had ever run

`imp_bff_dye_pdb2cif` passed `--dye-id`'s `None` default into a
`const std::string&`, so its documented default invocation raised a
`TypeError` before reading anything. `imp_bff_labelizer --show`, given a
structure instead of a container, failed with "does not begin with an EBML
header" -- true, and useless. Both fixed, both now in
`test/io/test_bin_programs.py`.

`imp_bff dock` runs (score 188.30 over 99 distances, 17 accessible volumes),
and so do `build-system`, `simulate` and the six `dye` commands.

## 2026-08-26 (night, last +6) — nine dead paths, found by running the docs

The unit suite was green at every step of the last three passes, and the CLI
and the example gallery were quietly rotting the whole time: **what breaks when
a C++ port changes a shape is the code nobody runs**, and nobody runs the
documentation.

Running every `bin/imp_bff` command and every example and notebook found nine.

### In the program

| command | what was wrong |
|---|---|
| `dye sample-dof-walk` | imported `LinkerSampler`, a name from before the linker sampler became C++; read `resolve_probe_site` as a dict; and then, once running, **accepted nothing** |
| `dye sample-langevin` | passed `timestep_fs=None` into a `double`, and asked `run()` to write an RMF, which it stopped doing when the sampler became C++ |
| `dye label-fusion` | `sys.exit` in a module that never imported `sys` |
| `flexfit` (FP branch) | `attach_probes` with four parallel lists, and another dict-style site |

`sample-dof-walk` is the interesting one. It rejected every proposal because
the probe is *bonded into* the site: its first atoms sit a bond length from
residues i-1 and i+1, and the walk counted those bonds as clashes. Excluding
the site and its two neighbours, and accepting a move that does not make the
clash count *worse* -- a walk that starts inside the protein can never leave,
otherwise -- gives 30/60 and 51/60 accepted with no clashes left. It also
built its geometry from a MOL2 while turning atoms read from a **PDB**, whose
order need not match; both come from the same file now.

### In the documentation

| page | what was wrong |
|---|---|
| `plot_convolution_routines`, `plot_pile_up` | `IMP.bff.decay_fconv*` -- decay convolution is **tttrlib's** by the placement rule; they call `tttrlib.fconv`/`fconv_per`/`fconv_simd` now |
| `hgbp1_label_and_sample` | `{"N", "CA", "C"} <= set(site)`, and `attached[0]["resnum"]` |
| `langevin_hgbp1_site481` | the `None` timestep and `out_rmf` |
| `plot_k2_uncertainty` | `k2_call` returned two values where its caller unpacked three -- it had *never* matched; `Kappa2Distribution` supplies the third |
| `t4l_pmi` | fetched AVs for a display helper that went with `representation/av.py` |
| `structure_cgdye` | `attach_probes` tuples, `SITE_KEEP_ATOM_NAMES` (a function now), a nested transition-count matrix (flat + `n` now), and a `RotamerLibrary` read as a dict |
| `structure_accessible_volumes` | stopped on a missing `ipyvolume`, which is a viewer and not a dependency |

### What holds them now

`test/expensive_test_docs_and_examples.py` runs **every** example and **every**
notebook
-- 42 of them -- with each notebook's own directory as the working directory,
which is what Jupyter does and what their relative paths assume. Opt-in, like
the other `expensive_test_` files. `test/label/test_dye_commands.py` runs six
commands for real, and asserts the walk *accepts* something: a sampler that
accepts nothing passes any "does it run" check.

### One decision left

`doc/manual/decays/decay_curves.ipynb`, `decay_forward_model.ipynb`,
`decay_objective_function.ipynb` and `doc/manual/programming/programming_imp_decorator.ipynb`
document `IMP.bff.DecayCurve`, `DecayConvolution`, `DecayLifetimeHandler`,
`DecayScale`, `DecayPattern` and `DecayLinearization` -- 50 references to an
API **deliberately deleted** by `93b7198` ("PRD-113 stage 0: delete the TCSPC
instrument layer"), and which exists nowhere in the stack now, tttrlib
included. They are not broken call sites; they are pages describing a layer
that left. Porting them to tttrlib's manual or dropping them is a call for the
maintainer, so they are skipped by name in the sweep, with the reason written
where the skip is.

## 2026-08-26 (night, last +5) — two landmines defused

**`.drot` could not read a payload that compressed well.** `DrotReader`'s
decompressor sized its output buffer as `compressed * 8` and grew it on
`NEEDS_MORE_OUTPUT`; brotli's one-shot decoder reports a too-small buffer as
`BROTLI_DECODER_RESULT_ERROR`, so the growth branch never fired and a good file
came back as "corrupt brotli stream". Every shipped `.drot` is under the ratio,
which is why nothing had noticed -- but a library of near-identical conformers
(a rigid dye on a short linker) is not. It streams now.

**Four `attach_probes` calls in `bin/imp_bff` were dead.** Three passed
`(hierarchy, chain, residue)` tuples and one passed four parallel lists; the
C++ takes `ProbeAttachment` values, so all four raised `TypeError` at the
binding. `dye sample-rotamer` then read `attached[0]["site"]` out of a list of
values. The command runs end to end now, and `resolve_probe_site` -- which
answers CA, N, C -- is where the site comes from.

Neither had a test, which is why neither was known to be dead. Both do now:
`test/io/test_drot.py::test_a_payload_that_compresses_hugely_still_reads` and
`test/label/test_attach_probes_cli.py` (3 tests, including the command itself
through `CliRunner`).

## 2026-08-23 — PRD-118: `.drot` v4 — one `.xz` container, self-contained

The ship format is settled: a `.drot` **is** a single `.xz` container (CRC32;
`xz -t` verifies every shipped file) whose payload carries everything —
the template as mmCIF (`_atom_site`), the Z-matrix row table, the per-rotamer
grids (int16: base 0.001 Å, θ/φ 0.1°) and varint weights. No sidecar `.pdb`:
`load_drot(path)` returns `ZTree` + coordinates + weights from the one file,
rebuilding the tree from the embedded template and cross-checking it against
the embedded rows (a built-in consistency check). Joint compression over the
whole payload beat v3's per-stream containers by 0.6 %, so the inner streams
went away; self-containment costs ~1 KB/file. **Sweep: 30.05 MB `.bcif` →
6.55 MB `.drot` = 4.59× smaller (78 % saved)**, 93/95 reconstruct exactly
(worst 0.018 Å), self-contained load of all ~30 k conformers in ~6 s; the
one sub-1× library (`A35_C1R_c30`) has a single rotamer — its 1.1 KB is
mostly the embedded template. Format and commands in
`prototypes/drot_rotlib/CONVERTING.md`.

# Update Log

## 2026-08-26 (night, last +7) — a use-after-free that returned zeros

Chasing the last of the untested surface -- the two `bin/` programs nothing
covered, and the top-level commands -- turned up two real defects and one that
had been silently wrong for as long as the CIF writer has existed.

### A system CIF written anywhere else could not be read back

The writer records each component's MOL2 path *relative to the CIF*, and
computed the `..` count from the **separators** of the base directory rather
than its segments -- one short whenever the base has no trailing separator,
which is always:

```
base   /var/folders/cl/xxxx/T/tmp1        6 segments
mol2   /Users/me/dye.mol2
wrote  ../../../../../Users/me/dye.mol2   5 ups
read   /var/Users/me/dye.mol2             does not exist
```

Nothing noticed because every test writes its systems into the data tree,
where the common prefix is long and the miscount cancels. `imp_bff simulate`
was dead for any system built outside it. `relative_path` in
`TopologyBuild.cpp` counts segments now, and
`test/cgdye/test_system_paths.py` writes a system into a `tmp_path` -- which
is exactly the case that failed.

### A container read off a temporary was freed memory

Writing the test for that found something worse:

```python
IMP.bff.read_forcefield_cif(path).components.values()
```

returned `'\x00\x00\x00...'` where a path should be, and segfaulted under
pytest. An accessor returning `const std::map&` gives SWIG a **borrowed**
pointer into the owner, and the proxy does not keep the owner alive: the
system is a temporary, it dies as the expression unwinds, and the map proxy
outlives it. Silently, most of the time.

The accessors cannot return by value -- `Scoring.cpp` calls `get_bonds()` in
the *condition* of a loop over the bonds, so a copy per access is quadratic --
so the copy is made **at the language boundary**: `%owned_container_out` in
`IMP_bff.types.i` gives every `const container&` return an owning copy, which
costs one copy per Python access and hands the lifetime to Python. Eight
container types: the maps, and the vectors whose element is one of this
module's own values.

**Not** `std::vector<std::string>`, `<double>` or `<int>`. Those already come
back *converted* -- a Python list, a numpy array -- which is a copy, so they
were never borrowed; claiming them replaced the conversion with a raw proxy and
five tests stopped recognising their own results (`assert <RMF_HDF5.Strings> ==
('SD',)`). The first cut of this fix did claim them, and the suite said so.
`%naturalvar` is on beside it, for the same hazard in plain member variables.

`test/test_borrowed_containers.py` reads five kinds of container off a
temporary, and checks that the copy really is one.

### Two programs nothing had ever run

`imp_bff_dye_pdb2cif` passed `--dye-id`'s `None` default into a
`const std::string&`, so its documented default invocation raised a
`TypeError` before reading anything. `imp_bff_labelizer --show`, given a
structure instead of a container, failed with "does not begin with an EBML
header" -- true, and useless. Both fixed, both now in
`test/io/test_bin_programs.py`.

`imp_bff dock` runs (score 188.30 over 99 distances, 17 accessible volumes),
and so do `build-system`, `simulate` and the six `dye` commands.

## 2026-08-26 (night, last +6) — nine dead paths, found by running the docs

The unit suite was green at every step of the last three passes, and the CLI
and the example gallery were quietly rotting the whole time: **what breaks when
a C++ port changes a shape is the code nobody runs**, and nobody runs the
documentation.

Running every `bin/imp_bff` command and every example and notebook found nine.

### In the program

| command | what was wrong |
|---|---|
| `dye sample-dof-walk` | imported `LinkerSampler`, a name from before the linker sampler became C++; read `resolve_probe_site` as a dict; and then, once running, **accepted nothing** |
| `dye sample-langevin` | passed `timestep_fs=None` into a `double`, and asked `run()` to write an RMF, which it stopped doing when the sampler became C++ |
| `dye label-fusion` | `sys.exit` in a module that never imported `sys` |
| `flexfit` (FP branch) | `attach_probes` with four parallel lists, and another dict-style site |

`sample-dof-walk` is the interesting one. It rejected every proposal because
the probe is *bonded into* the site: its first atoms sit a bond length from
residues i-1 and i+1, and the walk counted those bonds as clashes. Excluding
the site and its two neighbours, and accepting a move that does not make the
clash count *worse* -- a walk that starts inside the protein can never leave,
otherwise -- gives 30/60 and 51/60 accepted with no clashes left. It also
built its geometry from a MOL2 while turning atoms read from a **PDB**, whose
order need not match; both come from the same file now.

### In the documentation

| page | what was wrong |
|---|---|
| `plot_convolution_routines`, `plot_pile_up` | `IMP.bff.decay_fconv*` -- decay convolution is **tttrlib's** by the placement rule; they call `tttrlib.fconv`/`fconv_per`/`fconv_simd` now |
| `hgbp1_label_and_sample` | `{"N", "CA", "C"} <= set(site)`, and `attached[0]["resnum"]` |
| `langevin_hgbp1_site481` | the `None` timestep and `out_rmf` |
| `plot_k2_uncertainty` | `k2_call` returned two values where its caller unpacked three -- it had *never* matched; `Kappa2Distribution` supplies the third |
| `t4l_pmi` | fetched AVs for a display helper that went with `representation/av.py` |
| `structure_cgdye` | `attach_probes` tuples, `SITE_KEEP_ATOM_NAMES` (a function now), a nested transition-count matrix (flat + `n` now), and a `RotamerLibrary` read as a dict |
| `structure_accessible_volumes` | stopped on a missing `ipyvolume`, which is a viewer and not a dependency |

### What holds them now

`test/expensive_test_docs_and_examples.py` runs **every** example and **every**
notebook
-- 42 of them -- with each notebook's own directory as the working directory,
which is what Jupyter does and what their relative paths assume. Opt-in, like
the other `expensive_test_` files. `test/label/test_dye_commands.py` runs six
commands for real, and asserts the walk *accepts* something: a sampler that
accepts nothing passes any "does it run" check.

### One decision left

`doc/manual/decays/decay_curves.ipynb`, `decay_forward_model.ipynb`,
`decay_objective_function.ipynb` and `doc/manual/programming/programming_imp_decorator.ipynb`
document `IMP.bff.DecayCurve`, `DecayConvolution`, `DecayLifetimeHandler`,
`DecayScale`, `DecayPattern` and `DecayLinearization` -- 50 references to an
API **deliberately deleted** by `93b7198` ("PRD-113 stage 0: delete the TCSPC
instrument layer"), and which exists nowhere in the stack now, tttrlib
included. They are not broken call sites; they are pages describing a layer
that left. Porting them to tttrlib's manual or dropping them is a call for the
maintainer, so they are skipped by name in the sweep, with the reason written
where the skip is.

## 2026-08-26 (night, last +5) — two landmines defused

**`.drot` could not read a payload that compressed well.** `DrotReader`'s
decompressor sized its output buffer as `compressed * 8` and grew it on
`NEEDS_MORE_OUTPUT`; brotli's one-shot decoder reports a too-small buffer as
`BROTLI_DECODER_RESULT_ERROR`, so the growth branch never fired and a good file
came back as "corrupt brotli stream". Every shipped `.drot` is under the ratio,
which is why nothing had noticed -- but a library of near-identical conformers
(a rigid dye on a short linker) is not. It streams now.

**Four `attach_probes` calls in `bin/imp_bff` were dead.** Three passed
`(hierarchy, chain, residue)` tuples and one passed four parallel lists; the
C++ takes `ProbeAttachment` values, so all four raised `TypeError` at the
binding. `dye sample-rotamer` then read `attached[0]["site"]` out of a list of
values. The command runs end to end now, and `resolve_probe_site` -- which
answers CA, N, C -- is where the site comes from.

Neither had a test, which is why neither was known to be dead. Both do now:
`test/io/test_drot.py::test_a_payload_that_compresses_hugely_still_reads` and
`test/label/test_attach_probes_cli.py` (3 tests, including the command itself
through `CliRunner`).

## 2026-08-22 — PRD-117 batch 15: `pyext/src` deleted — zero .py files in pyext

The thirteen re-export shims (873 lines: `api.py`, `io/`, `io/cif.py`,
`label.py`, `scoring.py`, `representation/`, `restraints/`, `cgdye/*`) are
gone. Every public name is an attribute of `IMP.bff` itself — a SWIG wrapper
or a `%pythoncode` def in the one generated `__init__.py`:

- the two shims that held logic moved into `%pythoncode`: cif.i carries
  `forcefield_system_from_dict`/`as_forcefield_system`/`read|write_dye_
  forcefield_cif` (duck typing on Python objects has no C++ spelling; the
  write calls `_IMP_bff` directly so it cannot recurse into its own shadow),
  and swig.i-in's `__getattr__` is lazy-names-only (`_LAZY`), raising a real
  AttributeError for everything else — the api.py domain map died with the
  subpackages it indexed.
- 60+ files migrated to flat imports (`from IMP.bff import X`): 42 test
  files, `bin/imp_bff*`, the examples. Stale example surfaces fixed on the
  way (`IMP.bff.spectroscopy.kappa2`'s `s2delta`/`kappasq_all_delta` →
  `s2_delta_from_anisotropy`/`wobbling_kappa2_distribution_delta`;
  `representation.av.display_mean_av_positions` never survived the
  consolidation — the calls are commented with a note).
- The source-tree audits rewritten for the new shape: `test_public_api_names`
  (flat-name resolution, lazy build+cache, retired submodules stay gone,
  import survives click/RMF/IMP.pmi blocked), `test_import_discipline` (kept
  only its two source-independent guards — its own tripwire said "delete this
  file, its job is done"), `test_no_numba` (scans the `.i` `%pythoncode`
  blocks now), `test_shipped_files_compile` (examples + programs compile;
  nothing imports `IMP.bff.<submodule>`), `test_cgdye_integration` (dropped
  its vestigial `IMP.bff.cgdye` path helper; the program drives `bin/imp_bff`).

The build symlinks pyext/src into lib/IMP/bff at configure time; the stale
links were removed from the tree. `import IMP.bff` now loads RMF only behind
sim.i's try/except (unchanged) and nothing else optional.

Verification: `ninja IMP.bff` clean; non-medium **624 passed, 3 xfailed**
(count lower only because the deleted audits carried ~30 of their own tests;
sole failure is the pre-existing `av_reference_0.mrc` data file), cgdye
**162 passed**, expensive all-dyes **36 passed**, `imp_bff --help` runs.

## 2026-08-23 — PRD-118: `.drot` v3 — entropy-coded streams in `.xz` containers

Ship-size squeeze on top of the v2 exact store (the libraries must travel
with the package): measured the streams' entropy floor (order-0, 1.71×
headroom), then tested the candidates over the full corpus — per-stream LZMA
on the plain int16 grids wins (bit-packing *scrambles* the byte structure
the range coder exploits: 4.09 MB packed vs 3.59 MB plain on the probe);
weights are integer cluster populations, so they varint. The streams ship
inside **`.xz` containers (CRC32)**: +53 B/stream, +0.24 % over the raw
filter, for a built-in integrity check, `xz -t` debuggability and the simple
auto-detecting `lzma_stream_decoder` on the C++ side (filters self-describe
in the container header — no versioning hazard). **Final: 30.05 MB `.bcif`
→ 6.49 MB `.drot` = 4.63× smaller (78 % saved)**, per-library 2.09–13.83×
(median 4.33×); decode + full reconstruction of all ~30 k conformers in
~6 s; validation unchanged (93/95 exact, the 2 known-template-defect libs
apart). Numbers and rationale in `prototypes/drot_rotlib/CONVERTING.md`.

## 2026-08-23 — PRD-118: `.drot` space saving measured over all 95 libraries — and a shipped-template defect found

Encoded every recoverable DCD and compared against the shipped `.bcif` on
disk (`prototypes/drot_rotlib/results_space.json`):
**30.05 MB → 10.70 MB, 2.81× smaller (64 % saved)**, consistent per library
(2.67–3.15×); 3.02× against the DCD form. Validation in the same sweep:
**93/95 libraries reconstruct exactly** (worst mean RMSD 0.018 Å, all
≤ 0.05 Å). The two failures are a data defect, not a converter limit:
**`Cy3b_C2R.pdb` is byte-identical to `CF660R_C2R.pdb`** (residue `CF6` in
both) — a copy-paste template whose cutoff10 trajectory (real Cy3b
conformers) cannot match it, showing up as 10 Å template-vs-frames bond
deviations where every healthy library sits at ≤ 0.02 Å. Recorded on
PRD-116's registry-honesty list; numbers and caveats in
`prototypes/drot_rotlib/CONVERTING.md`.

# Update Log

## 2026-08-26 (night, last +7) — a use-after-free that returned zeros

Chasing the last of the untested surface -- the two `bin/` programs nothing
covered, and the top-level commands -- turned up two real defects and one that
had been silently wrong for as long as the CIF writer has existed.

### A system CIF written anywhere else could not be read back

The writer records each component's MOL2 path *relative to the CIF*, and
computed the `..` count from the **separators** of the base directory rather
than its segments -- one short whenever the base has no trailing separator,
which is always:

```
base   /var/folders/cl/xxxx/T/tmp1        6 segments
mol2   /Users/me/dye.mol2
wrote  ../../../../../Users/me/dye.mol2   5 ups
read   /var/Users/me/dye.mol2             does not exist
```

Nothing noticed because every test writes its systems into the data tree,
where the common prefix is long and the miscount cancels. `imp_bff simulate`
was dead for any system built outside it. `relative_path` in
`TopologyBuild.cpp` counts segments now, and
`test/cgdye/test_system_paths.py` writes a system into a `tmp_path` -- which
is exactly the case that failed.

### A container read off a temporary was freed memory

Writing the test for that found something worse:

```python
IMP.bff.read_forcefield_cif(path).components.values()
```

returned `'\x00\x00\x00...'` where a path should be, and segfaulted under
pytest. An accessor returning `const std::map&` gives SWIG a **borrowed**
pointer into the owner, and the proxy does not keep the owner alive: the
system is a temporary, it dies as the expression unwinds, and the map proxy
outlives it. Silently, most of the time.

The accessors cannot return by value -- `Scoring.cpp` calls `get_bonds()` in
the *condition* of a loop over the bonds, so a copy per access is quadratic --
so the copy is made **at the language boundary**: `%owned_container_out` in
`IMP_bff.types.i` gives every `const container&` return an owning copy, which
costs one copy per Python access and hands the lifetime to Python. Eight
container types: the maps, and the vectors whose element is one of this
module's own values.

**Not** `std::vector<std::string>`, `<double>` or `<int>`. Those already come
back *converted* -- a Python list, a numpy array -- which is a copy, so they
were never borrowed; claiming them replaced the conversion with a raw proxy and
five tests stopped recognising their own results (`assert <RMF_HDF5.Strings> ==
('SD',)`). The first cut of this fix did claim them, and the suite said so.
`%naturalvar` is on beside it, for the same hazard in plain member variables.

`test/test_borrowed_containers.py` reads five kinds of container off a
temporary, and checks that the copy really is one.

### Two programs nothing had ever run

`imp_bff_dye_pdb2cif` passed `--dye-id`'s `None` default into a
`const std::string&`, so its documented default invocation raised a
`TypeError` before reading anything. `imp_bff_labelizer --show`, given a
structure instead of a container, failed with "does not begin with an EBML
header" -- true, and useless. Both fixed, both now in
`test/io/test_bin_programs.py`.

`imp_bff dock` runs (score 188.30 over 99 distances, 17 accessible volumes),
and so do `build-system`, `simulate` and the six `dye` commands.

## 2026-08-26 (night, last +6) — nine dead paths, found by running the docs

The unit suite was green at every step of the last three passes, and the CLI
and the example gallery were quietly rotting the whole time: **what breaks when
a C++ port changes a shape is the code nobody runs**, and nobody runs the
documentation.

Running every `bin/imp_bff` command and every example and notebook found nine.

### In the program

| command | what was wrong |
|---|---|
| `dye sample-dof-walk` | imported `LinkerSampler`, a name from before the linker sampler became C++; read `resolve_probe_site` as a dict; and then, once running, **accepted nothing** |
| `dye sample-langevin` | passed `timestep_fs=None` into a `double`, and asked `run()` to write an RMF, which it stopped doing when the sampler became C++ |
| `dye label-fusion` | `sys.exit` in a module that never imported `sys` |
| `flexfit` (FP branch) | `attach_probes` with four parallel lists, and another dict-style site |

`sample-dof-walk` is the interesting one. It rejected every proposal because
the probe is *bonded into* the site: its first atoms sit a bond length from
residues i-1 and i+1, and the walk counted those bonds as clashes. Excluding
the site and its two neighbours, and accepting a move that does not make the
clash count *worse* -- a walk that starts inside the protein can never leave,
otherwise -- gives 30/60 and 51/60 accepted with no clashes left. It also
built its geometry from a MOL2 while turning atoms read from a **PDB**, whose
order need not match; both come from the same file now.

### In the documentation

| page | what was wrong |
|---|---|
| `plot_convolution_routines`, `plot_pile_up` | `IMP.bff.decay_fconv*` -- decay convolution is **tttrlib's** by the placement rule; they call `tttrlib.fconv`/`fconv_per`/`fconv_simd` now |
| `hgbp1_label_and_sample` | `{"N", "CA", "C"} <= set(site)`, and `attached[0]["resnum"]` |
| `langevin_hgbp1_site481` | the `None` timestep and `out_rmf` |
| `plot_k2_uncertainty` | `k2_call` returned two values where its caller unpacked three -- it had *never* matched; `Kappa2Distribution` supplies the third |
| `t4l_pmi` | fetched AVs for a display helper that went with `representation/av.py` |
| `structure_cgdye` | `attach_probes` tuples, `SITE_KEEP_ATOM_NAMES` (a function now), a nested transition-count matrix (flat + `n` now), and a `RotamerLibrary` read as a dict |
| `structure_accessible_volumes` | stopped on a missing `ipyvolume`, which is a viewer and not a dependency |

### What holds them now

`test/expensive_test_docs_and_examples.py` runs **every** example and **every**
notebook
-- 42 of them -- with each notebook's own directory as the working directory,
which is what Jupyter does and what their relative paths assume. Opt-in, like
the other `expensive_test_` files. `test/label/test_dye_commands.py` runs six
commands for real, and asserts the walk *accepts* something: a sampler that
accepts nothing passes any "does it run" check.

### One decision left

`doc/manual/decays/decay_curves.ipynb`, `decay_forward_model.ipynb`,
`decay_objective_function.ipynb` and `doc/manual/programming/programming_imp_decorator.ipynb`
document `IMP.bff.DecayCurve`, `DecayConvolution`, `DecayLifetimeHandler`,
`DecayScale`, `DecayPattern` and `DecayLinearization` -- 50 references to an
API **deliberately deleted** by `93b7198` ("PRD-113 stage 0: delete the TCSPC
instrument layer"), and which exists nowhere in the stack now, tttrlib
included. They are not broken call sites; they are pages describing a layer
that left. Porting them to tttrlib's manual or dropping them is a call for the
maintainer, so they are skipped by name in the sweep, with the reason written
where the skip is.

## 2026-08-26 (night, last +5) — two landmines defused

**`.drot` could not read a payload that compressed well.** `DrotReader`'s
decompressor sized its output buffer as `compressed * 8` and grew it on
`NEEDS_MORE_OUTPUT`; brotli's one-shot decoder reports a too-small buffer as
`BROTLI_DECODER_RESULT_ERROR`, so the growth branch never fired and a good file
came back as "corrupt brotli stream". Every shipped `.drot` is under the ratio,
which is why nothing had noticed -- but a library of near-identical conformers
(a rigid dye on a short linker) is not. It streams now.

**Four `attach_probes` calls in `bin/imp_bff` were dead.** Three passed
`(hierarchy, chain, residue)` tuples and one passed four parallel lists; the
C++ takes `ProbeAttachment` values, so all four raised `TypeError` at the
binding. `dye sample-rotamer` then read `attached[0]["site"]` out of a list of
values. The command runs end to end now, and `resolve_probe_site` -- which
answers CA, N, C -- is where the site comes from.

Neither had a test, which is why neither was known to be dead. Both do now:
`test/io/test_drot.py::test_a_payload_that_compresses_hugely_still_reads` and
`test/label/test_attach_probes_cli.py` (3 tests, including the command itself
through `CliRunner`).

## 2026-08-23 — PRD-118: the FRETpredict → `.drot` conversion runbook

`prototypes/drot_rotlib/CONVERTING.md` documents how a library entry becomes
`.drot`, every command verified end to end: DCD recovery from git `8fac573`,
encode (CLI weights now optional — omitting them leader-clusters a raw
trajectory, `cutoffN` → N/100 Å threshold), validate (acceptance: mean RMSD
≤ 0.015 Å, max ≤ 0.05 Å), the batch loop, and the raw-MD route. Two caveats
found while verifying are now written down rather than rediscovered: mdtraj
needs `top=` for XTC (the template PDB is it, and warns about the dye/linker
sharing resid 2 — harmless), and the one recoverable raw trajectory
(`A48_C1R/traj.xtc`, 28110 frames) has **92 atoms against the shipped 83** —
the upstream strip is undocumented and the upstream topology is not in the
repo, so the git-recovered DCDs are the supported conversion source. Linked
from the PRD's input-contract section and the prototype README.

## 2026-08-22 — PRD-118: the rotamer core becomes IMP-native and general — `IMP.bff.ZMatrix`

The FASPR port's engine stays vendored verbatim behind `faspr_pack` (that is
what the parity pin guards), but the *shared machinery* is now a proper IMP
citizen: `include/ZMatrix.h` + `src/ZMatrix.cpp` expose the IUPAC
`dihedral_deg`, `bond_angle_deg`, `internal2cartesian` (FASPR's
`Internal2Cartesian` construction, double precision, the exact inverse to
~1e-13), geometry-based `perceive_bonds`, and a `ZMatrix` class — a
molecule's spanning-tree internal coordinates with exact `encode`/`decode`
(`encode(decode(chi)) == chi` for arbitrary chi; the template round-trips to
itself to 1e-13).

**General over the molecule, by test**: the A48_C1R dye+linker template
perceives exactly its 87 MOL2 bonds and round-trips exactly; a protein
residue does the same through the same class; the chi the kernel reads off a
FASPR-repacked LYS equals the mdtraj-verified −65.9°. The dihedral
convention is pinned against an independent numpy implementation — the
standard+180 bug of the prototype survived every self round-trip, so the
test suite now makes the external convention check structural
(`test/representation/test_zmatrix.py`, 7 tests; faspr/representation/API
suites: 124 passed). SWIG: `VectorVector3D` instantiated for the
`std::vector<Vector3D>` surfaces.

What this buys: the remaining `.drot` half (store + `rotamer.i` loader +
pin re-derivation) is now assembly on a tested kernel rather than a port.

## 2026-08-22 — PRD-119: quenching anisotropy from the rotamers, decision-gated

Wrote [`okf/prds/prd-119.md`](prds/prd-119.md) — plan only, nothing
implemented. Two questions in the cheapest order:

- **Stage 0 is the continuum's trial.** Per-site reachability (17-parameter
  per-site overfit through the surrogate, after checking the surrogate's
  held-out 0.057 ns is not the limiter) plus the s33b/s33d floors, per site.
  If ≥ ~85–90 % of the 41 sites are reachable and the global fit beats the
  parameter-free s30 baseline (+0.571) by a real margin, the continuum is
  adequate and **the PRD closes with that finding** — the rotamer machinery
  is unnecessary. Honest prior: marginal (+0.573 vs +0.571).
- **Only if it fails:** the exact rotamer master equation
  `M = L − diag(k)` — exchange generator L **measured from the shipped
  A48_C1R MD** (reversible MSM on cutoff10 centroids; hard gate: reproduce
  the direct ⟨P2(μ(0)·μ(t))⟩), per-rotamer rates `k_i(d, χ)` on real
  conformer atoms with a ring-orientation factor, one global site-specific
  exchange-scaling exponent (contact rotamers slowed by
  `slow_factor^n_contact`). Every stage has a numeric kill criterion:
  frozen/Lipari–Szabo limits, emergent static-bound species (the 2026-08-22
  channel becomes a theorem of the formulation), `r∞/r0 = ⟨P2⟩` of the
  slowest mode, and stage 3 must beat both the continuum fit and the s30
  the slowest mode, and stage 3 must beat both the continuum fit and the s30
  baseline with the continuum-dead sites as headline. Then the learned
  components — **the point of the differentiable formulation, per the
  user's refinement: backprop to the interaction potential and to the
  quenching** (`s43` learned `k(d, χ)` surface, payoff = does it rediscover
  the Vaiana π-stacking criterion; `s43b` learned correction to the
  screening potential that sets the rotamer weights, payoff = does
  stickiness localise to Lys/Arg; both zero-regularised and
  Fisher-checked), with the deeper rungs (torsion PMF via a torch
  `internal2cart`, learned L correction) as conditional extensions.
  Surrogate/eigen/VV-VH machinery reused
  unchanged; a UDE in PRD-115's sense, not a collocation PINN. v1 of the
  plan built the rotamer model unconditionally; the user's correction —
  *if the continuum works, no rotamer is needed* — made stage 0 the gate.

## 2026-08-22 — PRD-117 phase 2, batch 14: `cif.i` to C++ — template CIF reader written, JSON system bridge

`pyext/IMP_bff.cif.i` lost all 558 lines of `%pythoncode` (the
`forcefield_system_from_dict` builder, the ihm.format template reader/writer
pair, the rotamer-library numpy IO and the `_parse_*_atom_site_rows` helpers).
The C++ surface in `CifIO.h` takes over every name:

- `read_component_template_cif`/`read_dye_template_cif` are **real now**: the
  shipped stub returned an empty template; this batch wrote the ihm C-reader
  integration (`_cgdye_template`/`_cgdye_feature`/`_cgdye_feature_atom`/
  `_cgdye_improper`/`_cgdye_metadata`) the same way ForceFieldCIF.cpp reads
  `_ff_*`. Returns the `ComponentTemplate` value (`features` a `FeatureMap`,
  rows flatnested `Feature`/`FeatureAtom`/`Improper`).
- The CifWriter quotes values at the writer, not per call site: a metadata
  value like `N1 and resname A48, N2 and resname A48` unquoted made the ihm
  reader split the row and the round trip failed short.
- `forcefield_system_from_dict` is a nlohmann bridge in C++
  (`forcefield_system_from_json` in DyeForceField.cpp): every field the dict
  literal carried, string-typed fields coerced with `str()` semantics (a
  numeric `"id": 0` is `"0"`), site references by `site_no` resolved. The
  writer lost its underscore (`write_dye_forcefield_cif`).
- `normalize_weights` returns the normalized copy (value semantics; the old
  in-place dict mutation has no spelling across the boundary).
- `as_forcefield_system`/`forcefield_system_from_dict`/the dict-accepting
  `write_dye_forcefield_cif` stay in the `IMP.bff.io.cif` package module --
  duck typing on a Python object has no C++ spelling.

Migrated consumers: `topology.i` (template value surface + the system bridge),
`sim.i` (`read_forcefield_cif` direct, sampling overrides on the typed value;
its `as_forcefield_system` indirection is deferred -- a module-scope io import
is circular, and shadowing the SWIG writer recursed), the io.cif shim, and the
tests (`test_template`, `test_dye_io`, `test_io`, `test_dye_topology`,
`test_rotamer_generation/sampling`, `test_physics_invariants`,
`test_combined_system`, the expensive all-dyes sweep -- whose `rng=` kwarg had
been dead since the Phase-1 dyesampling port).

Verification: `ninja IMP.bff` clean; the non-medium suite **656 passed,
3 xfailed** (only the pre-existing `test_access_av_feature` data failure);
cgdye **149 passed** and the expensive all-dyes sweep **36 passed**.

## 2026-08-22 — PRD-118: FASPR itself ported 1:1 to C++ — `IMP.bff.faspr_pack`

The first half of PRD-118's C++ stage: FASPR (version 20200309) is now in the
module, not just cloned. `utility/port_faspr.py` vendors
`junk/FASPR/src/*.cpp` into flat `src/Faspr*.cpp` + `src/faspr/*.h` inside
`namespace IMP::bff::faspr`, with **behaviour-identical, mechanical transforms
only**: namespace wrap, include-guard prefixes, `exit(0)` → throwing
`faspr_fail()`, MSVC pragmas dropped, one `std::swap` qualification (nested
namespace hiding, caught by the amalgamated `bff_all.cpp` build). The MIT
notice travels with every file. `src/Faspr.cpp` is FASPR's `main()` pipeline
minus the CLI, behind `IMP.bff.faspr_pack(pdb_in, pdb_out, rotlib, verbose)`
(`include/Faspr.h`, exposed through `swig.i-in`; FASPR's chatty stdout is
buf-swapped silent unless `verbose`).

**Parity is pinned, not asserted**: `test/faspr/test_faspr_port.py` repacks
T4L 3GUN with both the port and the *reference executable* built from the
same sources — atom sets, coordinates (≤ 0.002 Å) and every chi1 (≤ 0.5°)
must agree; backbone N/CA/C/O must be untouched; a missing library must
throw. All 4 tests pass, plus the surrounding suites (97 in the API/import/
vendored-headers/no-numba block). Example:
`examples/structure/faspr_sidechains.py`.

**The Dunbrack-2010 binary library (13.8 MB) is deliberately NOT
redistributed** — free for academic use, redistribution decision open in
PRD-118; callers pass its path, tests and the example skip politely without
it (`$IMP_BFF_FASPR_ROT_LIB` or the `junk/FASPR` clone).

Still open (second half): port the dye-rotamer `.drot` store
(`internal2cart` kernel, `rotamer.i` loader hook, pin re-derivation on
reconstructed frames).

## 2026-08-22 — PRD-118 A/B: side-chain packer vs the reference FASPR, and the 180° convention bug it caught

Built the FASPR executable from `junk/FASPR/` and ran `ab_faspr.py`
(`prototypes/drot_rotlib/`): repack 3GUN + mGBP2A with the reference, compare
every rotatable residue's chi against the prototype's top-prior Dunbrack
packer. **Assertion passes: 87/87 and 354/354** — on every residue where
FASPR's global-GMEC search lands on the top-probability rotamer, the
prototype reproduces it exactly (median per-chi deviation 0.0°); the
remaining 36 %/31 % is FASPR's pair-energy optimisation overriding the
library prior, an algorithmic difference by design.

**The first honest run failed 136/136 — every chi ~180° off — and that was a
real bug, not a comparison artefact:** the prototype's `dihedral_deg`
measured *standard + 180°*. It was exactly self-consistent with its own
placement primitive, so every internal round-trip check passed, while every
rotamer it built sat mirrored off the IUPAC/Dunbrack convention *and* every
(phi, psi) lookup read the mirrored bin. Two lessons now pinned in the
prototype docs: a convention that round-trips against itself proves nothing
about the world's convention, and an A/B against an independent reference
implementation is the check that catches it. After the fix (`rotlib.py`):
ZTree self-check 2.5e-14 Å, `.drot` v2 exactness unchanged (0.014 Å, 0.12°
dipole), side-chain round-trips exact, both A/B assertions pass.

## 2026-08-21 — PRD-117 phase 2, batch 13: `structureio.i` to C++

`pyext/IMP_bff.structureio.i` lost its reshape/dict/set wrappers (`read_dcd`
`read_dcd_header`, `read_trajectory`, `parse_pdb_atoms`, `parse_conect_bonds`,
`infer_bonds`, `write_mol2`, `read_score_series`, `write_pdb`, `compute_rmsd`,
`load_structure`, `structure_coordinates`, `read_xlink_table`, `apply_transform`)
and all twelve `_function` renames; the C++ functions are the public surface.
Flat coordinate buffers come back as numpy views; the `(n,3)` /
`(n_frames,n_atoms,3)` reshape is the caller's. Interfaces that changed:

- `read_dcd_header(path)` returns the `DCDHeader` value (`head.n_frames`,
  `head.n_atoms`, `head.endianness`, ...), not a dict; `read_dcd` returns a flat
  view (`max_frames` now defaults to -1) that callers reshape with the header's
  `n_atoms`.
- `parse_pdb_atoms`/`parse_conect_bonds`/`infer_bonds`/`read_xlink_table` are
  the C++ record lists/sets (iterate `AtomBond.a/.b` etc.) rather than dicts.
- `write_pdb(coords.ravel(), str(path), chain, res_name, transform, model_index)` /
  `compute_rmsd(a.ravel(), b.ravel(), mask, superpose)` / `apply_transform`.

Stays Python, deliberately: the RMF writers (`_LAZY`), which would pull
`IMP.rmf` into `required_modules`, and the IMP-object glue
(`load_structure_with_particles`, `read_angle_file`), which creates SWIG
objects (`read_pdb_hierarchy` is re-usable from Python since it is a C++
`Hierarchy` now).

Migrated: `test/io/test_formats.py`, `test/test_dcd_reader.py`,
`test/expensive_test_dcd_reader.py`, `bin/imp_bff_traj2bcif`.

Verification: `ninja IMP.bff` clean; the non-medium suite **652 passed,
3 xfailed** -- only the pre-existing `test_access_av_feature` data failure.

## 2026-08-21 — PRD-117 phase 2, batch 12: `fps.i` to C++

`pyext/IMP_bff.fps.i` lost its entire `%pythoncode` block (the module
constants, the `_fields_as_dict` tables, the `to_json_schema`/`validate_*`
dict bridges and the `_io.read_*`/`_write_*` wrappers) and all eleven
`_function` renames. The C++ `fps_json_schema`, `validate_position`,
`validate_distance`, `fps_schema_validate`, `read_old_lps_txt`,
`read_old_distances_txt`, `read_fps_json`, `write_fps_json`,
`fps_positions_for_docking`, `read_evaluators_json`, `write_evaluators_json`,
`fps_simulation_types`, `fps_av_simulation_types`, `fps_distance_types`,
`fps_position_fields` etc. are the public surface, speaking JSON text (a dict
is passed with `json.dumps(...)`; a returned `FPSDocument`'s sections are JSON
text the caller parses).

Interfaces that changed (prerelease; every caller migrated):

- `read_fps_json(str(path))` returns an `FPSDocument`; `positions`/`distances`/
  `score_sets`/`extra` are JSON strings, not dicts. No `p, d, s, e =` unpack.
- `write_fps_json(str(path), json.dumps(...), ..., validate=True)`.
- The module constants (`SIMULATION_TYPES`, `AV_SIMULATION_TYPES`,
  `DISTANCE_TYPES`, `SCHEMA_VERSION`, `POSITION_FIELDS`, ...) are the calls
  (`fps_simulation_types()`, `fps_av_simulation_types()`,
  `fps_distance_types()`, `fps_schema_version()`, ...); the io shim re-exports
  the calls.
- `validate_position`/`validate_distance`/`fps_schema_validate` take JSON text
  and return an `FPSValidation` (`.errors`/`.warnings` tuples, `.is_valid`).
- `histogram_rda(s1, s2)` gained defaults (empty axis spans the sample, its
  own `n_samples + 1` edges, normalized) so the "just a histogram" caller works
  without naming an axis.

Migrated consumers: `test/io/test_fps_schema.py`, `test/io/test_formats.py`,
`test/test_AVNetworkRestraint.py`, `pyext/IMP_bff.docking.i` (the
`_FPS_JSON_CACHE` conversion, now preserving `doc.*` JSON directly),
`pyext/IMP_bff.rotamer.i` (`read_rotamer_fps`, `write_rotamer_fps`), and
`pyext/IMP_bff.rotamer_ensemble.i` (`rotamer_ensembles_from_fps`). The
`RotamerEnsemble` subclass gained the `States` attribute surface (`points`,
`orientations`, `params`, `mean_position`, ...) over the C++ get_*/set_*
methods, and its pair-physics methods now adapt a plain `AccessibleVolume`
other as well. `load_rotamer_library_dcd` fixed at its single consumer to
pass `str(...)` paths (the old wrapper str'd them). cgdye/rotamer (24 tests)
green.

Verification: `ninja IMP.bff` clean; the non-medium suite **652 passed,
3 xfailed** -- only the pre-existing `test_access_av_feature` data failure.

## 2026-08-21 — PRD-117 phase 2, batch 11: `avmeandistance.i` to C++ (Phase-1 stragglers clear)

`pyext/IMP_bff.avmeandistance.i` lost the `%extend %pythoncode` property block
on `PositionUncertainty` (rmsf/mean_coords are the `get_*()` views, and
`rmsf_mean`/`rmsf_max`/`mobile_rmsf_mean` are `%attribute`s now) and the
`%pythoncode` wrapper around `estimate_position_uncertainty`, whose only job
was building a dict from the C++ record and writing the two output files. The
`_estimate_position_uncertainty` rename is gone; the C++ function and
`write_position_uncertainty_pdb`/`write_position_uncertainty_csv`/
`pdb_chain_ids` are the public surface. Nothing consumed the dict shape (no
test, shim or program), so there is nothing to migrate.

The `_LAZY["AVNetworkRestraintWrapper"]` builder stays Python: it subclasses
`IMP.pmi.restraints.RestraintBase`, and defining it at import time would make
`import IMP.bff` require IMP.pmi (the module's `required_modules` do not).

`avmodel.i` turns out to be already clean (the `%pythoncode` `rg` hit was a
comment); `av_pair_statistics`/`histogram_rda` are C++ in `StatesDistance.h`.

Verification: `ninja IMP.bff` clean; the non-medium suite **652 passed,
3 xfailed** — only the pre-existing `test_access_av_feature` data failure.

## 2026-08-21 — PRD-117 phase 1, batch 10: `quenching.i` to C++ — phase 1 complete

`pyext/IMP_bff.quenching.i` lost its one `%pythoncode` block (`_flat`/`_cube`/
`_kappa2` helpers and the ~dozen grid/map/trace reshapers) and all twelve
`_function` renames; the C++ `sphere_points`, `solvent_accessible_surface`,
`quenching_rate_per_frame`, `slow_factor_grid`, `quenching_rate_grid`,
`av_contact_mask`, `grid_axis`, `diffusion_coefficient_map`,
`quenching_rate_map`, `fret_rate_map`, `fret_rate_trace`,
`fret_rate_pair_trace` are the public flat surface. Interfaces/wrappers that
changed:

- The cube-reshaping (`_cube`), the `_flat` ravel and the `None`-kappa2
  resolution are caller-side now; `fret_rate_trace`/`fret_rate_pair_trace`
  take an explicit `kappa2` (pass `kappa2_isotropic()`).
- `radial_diffusion_map` is **C++** (QuenchingMap) and takes the profile
  pre-evaluated on an integer-Angstrom radius grid instead of a callable.

Two shared typemaps widened in `pyext/IMP_bff.types.i`: the
`std::vector<double>` `in` typemap now converts any array-like through
`PyArray_FROMANY(..., FORCECAST)` to a flat float64 copy (so an int-fold
density/mask reaches these kernels the way the Python `astype(float64)` did),
mirroring the `vector<int>` one.

Test callers migrated: `test_quenching_kernels` (helpers `_sas`/`_qgrid`/
`_sfgrid`/`_cmask`/`_ftrace`/`_fptrace`), `test_quenching_field`
(`_qrmap`/`_dcmap`/`_frmap` + `_axis`), `test_ported_kernels`, and
`test_photophysics_terms`.

**Phase 1 (the 17-file %pythoncode port) is done:** all of avdistance,
interactionterms, forcefield, observables, greedyolga, dye, distributions,
griddiffusion, avmodel, statesdistance, fret_pair_distribution, dyediffusion,
quenchingmodel, avbuilder, petquenching, dyesampling, photophysics and
quenching are %pythoncode-free. 14 `.i` files with `%pythoncode` remain (Phases
2-5): fps, cif, structureio, rotamer_ensemble, avmeandistance, topology,
scoring, label, sampling, rotamer, docking, sim.

Verification: `ninja IMP.bff` clean; the non-medium suite **652 passed,
3 xfailed** — only the pre-existing `test_access_av_feature` data failure.

## 2026-08-21 — PRD-117 phase 1, batch 9: `photophysics.i` to C++

`pyext/IMP_bff.photophysics.i` lost its one `%pythoncode` block (two dozen
reshape/tuple-split wrappers) and the `FRETRegimes` `%extend` dict surface;
the four `_function` renames are gone and the C++ kernels are the public flat
surface. Interfaces that changed (prerelease; callers migrated):

- `kappa(d0,d1,a0,a1)`/`kappa_distance(d1,d2,a1,a2)` are
  `dipole_kappa_distance(d1,d2,a1,a2)` -> two-value `VectorDouble`; `kappasq`
  is `wobbling_kappa2`.
- `kappa2_distribution_wobbling_in_cone`/`kappasq_all_delta` are
  `wobbling_kappa2_distribution`/`wobbling_kappa2_distribution_delta`: the
  `scale`/`hist` out-collections are now passed in as `VectorDouble()` and
  mutated in place; samples come back first (test helpers `_wobbling` /
  `_wobbling_delta` reshape and re-split).
- `kappa2_from_dipoles` -> `kappa2_dipole_matrix` (flat, caller reshapes to
  `(n_d, n_a)`).
- The convolution front door `convolve_distance_with_k2_ratio` was **moved
  into C++** (`LifetimeSpectrum.h/.cpp`): it runs `outer_product_histogram`,
  computes the product range, widens a zero-width (delta) range, and returns
  only the bins carrying more than `1e-10` of the peak, as two managed views.
- `FRETRegimes` no longer maps `["static"/"dynamic"/...]` (the dataclass dict
  surface); its members are attributes.

Verification: `ninja IMP.bff` clean; the non-medium suite **652 passed,
3 xfailed** — only the pre-existing `test_access_av_feature` data failure.

## 2026-08-21 — PRD-117 phase 1, batch 8: `dyesampling.i` to C++

`pyext/IMP_bff.dyesampling.i` lost its three `%pythoncode` blocks (the
`simulate_*` wrappers) and both `%extend` property blocks; the field-building,
the photon-trace split, and the in-place decay accumulation moved into C++.
Interfaces that changed (prerelease; callers migrated in `test/sampling/`,
`test/quenching/`):

- `simulate_dye_diffusion(density, slow_density=(), dg=0.5, t_max, t_step, D,
  slow_fact=(), random_seed=-1)` now builds the mobility **field in C++**: a
  length-`ng^3` `slow_fact` is the field, otherwise it is a scalar (applied
  where `slow_density` is nonzero) or empty (= uniform). A scalar `slow_fact`
  must be passed as `[x]`. The density is validated to be a perfect cube.
- `simulate_photon_trace(n_ph, k_quench, t_step=0.01, tau0=0.25,
  random_seed=-1)` returns `(delays, emitted)` as two views from one call.
- `simulate_quenched_decay(n_curves, decay, dt_tac, k_quench, t_step, tau0,
  random_seed=-1)` accumulates **in place** into the float64 histogram; a
  float32 histogram now raises rather than silently upcasts.
- `equilibrium_occupancy(diffusion_map, bounds, flux_form="smoluchowski")`
  returns a **flat** view; callers reshape to the cube (`closed_occupancy`).
- The trajectory's `.xyz`/`.accepted` are now `get_xyz()`/`get_accepted()`
  (flat; caller reshapes / casts to uint8); `.n_frames`, `.n_accepted`,
  `.n_rejected`, `.acceptance_ratio` are read-only `%attribute`s.
- `apply_rotamer_coordinates(hierarchy, coords)` moved to C++ (IMP::atom
  leaves -> IMP::core::XYZ), keeping the name `bin/imp_bff` imports.

Two shared typemaps were widened in `pyext/IMP_bff.types.i`: `const
std::vector<int>&` now converts any array-like through
`PyArray_FROMANY(..., FORCECAST)` to a flat int32 copy (float64 bounds, uint8
masks and F-contiguous float32 all land), and `swig.i-in` gained the
`INPLACE_ARRAY1` apply for the decay histogram before `PhotonSimulation.h`.

Verification: `ninja IMP.bff` clean; the non-medium suite **652 passed,
3 xfailed** — only the pre-existing `test_access_av_feature` data failure.

## 2026-08-21 — PRD-118: dihedral rotamer libraries (`.drot`) prototype + Dunbrack side-chain sampler

Wrote [PRD-118](prds/prd-118.md) and built the prototype in
`prototypes/drot_rotlib/` (python/numba; C++ port deferred — `src/` is in
flux under PRD-117). Two capabilities:

1. **`.drot` house format (driver a).** A library becomes a rigid template +
   the frame's own Z-matrix base/θ/φ + weights, with FASPR-style
   `Internal2Cartesian` reconstruction. Verified end to end on 5 recovered
   `.dcd`s (from git `8fac573`). **First pass stored only the moving dihedrals
   and looked like a hard 2-5 A / ~20 deg floor** — it was a format gap, not a
   limit: this MD also flexes bond angles (~3 deg/row), so the second pass
   stores base + θ + φ per row and is **lossless to the int grid** —
   reconstruction RMSD 0.009-0.015 A, transition-dipole error 0.09-0.13 deg
   median across five libs, `.drot` 0.33x the `.dcd` size. `.drot` is now a
   candidate drop-in for `.bcif`; remaining work is the C++ port + pin
   re-derivation. Table in `prototypes/drot_rotlib/results_validation.json`.
2. **Dunbrack side-chain sampler (driver b).** `dunbrack_sidechains.py` reads `junk/FASPR/dun2010bbdep.bin` directly and
   packs native side chains (`pack_protein_residue`) or continuous chi
   (`build_sidechain`) — PHE/GLU/TYR/ARG/ASN on 3GUN reproduce their chi to
   ~0.5°, PRO excluded (ring). Acknowledges FASPR (Bioinformatics 2020;
   36:3758-3765; github.com/tommyhuangthu/FASPR) in the PRD + prototype.

Found on the way: `import IMP.bff` is currently broken against the build tree
(`_IMP_bff.so` references a `compute_av_from_structure` overload the rebuilt
`libimp_bff` no longer exports) — the prototype therefore ships a standalone
DCD reader (`dcd.py`, bit-identical to IMP's on `A48_C1R_cutoff10`).

## 2026-08-21 — PRD-117 phase 1, batch 7: `petquenching.i` to C++

`pyext/IMP_bff.petquenching.i` lost its two `%pythoncode` blocks (the
`_pet_table`/`QUENCHER_ATOMS`/`PET_QUENCHING_REFERENCE`/
`STANDARD_AMINO_ACID_RESIDUES` constants-and-shims and the `ResidueSites`
`%extend` properties), the nine `_function` renames, and the two
`%attribute_py` (themselves `%pythoncode` macros — replaced with plain
`%attribute`). The C++ `quencher_atoms`, `pet_quenching_reference`,
`reference_*, normalize_*, *_for_residues`, `residue_sites` are now the public
SWIG surface, returning `std::map`/`std::vector` proxies and flat numpy views.

Interfaces that changed (prerelease; callers migrated):

- `QUENCHER_ATOMS`, `PET_QUENCHING_REFERENCE`, `STANDARD_AMINO_ACID_RESIDUES`
  module constants are gone — call the functions
  (`quencher_atoms()`, `pet_quenching_reference()`,
  `standard_amino_acid_residues()`). Their SWIG maps iterate and support `[]`,
  `in`, `len()`, `items()` but not `.get()`.
- `ResidueSites`: `.slow_centers`/`.quench_centers` properties are now
  `get_slow_centers()`/`get_quench_centers()` flat views (caller reshapes to
  `(n, 3)`); `len(found)` is `found.size()`.
- `residue_sites(..., table=None)` no longer accepts `None` for the table —
  the `%pythoncode` `_pet_table` shim that turned `None` into an empty map is
  gone; omit the argument instead.
- `Quencher.is_typed` / `PETParameters.is_transferred` are `%attribute`s now;
  `get_is_typed()`/`get_is_transferred()` methods are suppressed.

Migrated: `test/quenching/test_quenching_kernels.py`,
`test/quenching/test_quenching_model.py`, `test/test_label_system.py`
(`pet_quenching_reference` added to the `IMP.bff.label` shim).

Verification: `ninja IMP.bff` clean; non-medium suite **652 passed, 3 xfailed**,
only the pre-existing `test_access_av_feature` data failure.

## 2026-08-21 — PRD-117 phase 1, batch 6: `avbuilder.i` to C++

The last of the phase-1 five. `pyext/IMP_bff.avbuilder.i` lost all three
`%pythoncode` blocks and the `_compute_av`/`_compute_av_from_structure`
renames; `compute_av`, `compute_av_from_structure` and
`compute_avs_for_structure` are now the C++ functions SWIG exposes directly.

Two nlohmann JSON bridges were added to `AVBuilder.h`/`AVBuilder.cpp`
(vendored `internal/json.h`): the structure door takes one fps.json position
entry, and the batch door takes a whole `Positions` section, both serialised as
JSON text. Interfaces that changed (prerelease, callers migrated):

- `compute_av(atoms_xyzr, source_xyz, ...)` takes the `(N, 4)` x/y/z/vdW
  column stack directly instead of `(atoms_xyz, atoms_vdw, source_xyz, ...)` —
  the `%apply(double* IN_ARRAY2, …)` `(n_atoms, n_cols)` typemap collapses the
  three raw args into one 2-D view, and the pythoncode no longer
  `column_stack`s.
- `compute_av_from_structure(pdb, position_json, disc_step=-1.0)` takes a JSON
  **string**, not a dict; `disc_step <= 0` derives it from the position's
  `simulation_grid_resolution` and still refuses a declared one that disagrees
  (positions carrying `simulation_grid_resolution` without `disc_step`). The
  overload is SWIG-marshalled, so `disc_step` is positional (keywords vanish on
  overloads).
- `compute_avs_for_structure(positions_json, pdb_path, disc_step=-1.0)` returns
  `MapStringAccessibleVolume`; `pdb_path` may be a JSON array of paths indexed
  by `body_id`. Missing attachment sites come back as empty volumes, as before.

Callers migrated to `json.dumps(...)`: `test/representation/test_av_resolution.py`,
`test/representation/test_av_grid_registration.py`, and `bin/imp_bff:av_for_position`.
`pyext/src/representation/__init__.py` port-status note updated (15/17 covered).

Verification: `ninja IMP.bff` clean; the non-medium suite **652 passed,
3 xfailed** — only the pre-existing `test_access_av_feature` failure
(`references/av_reference_0.mrc` data missing, unrelated to the port) remains,
same as before this batch.

## 2026-08-21 — FASPR cloned into `junk/` as a possible rotamer / side-chain source

Cloned https://github.com/tommyhuangthu/FASPR into `junk/FASPR/` (Huang's fast
protein side-chain packing, Bioinformatics 2020; packs a backbone against the
Dunbrack 2010 rotamer library `dun2010bbdep.bin`, vendored). Two candidate uses
for imp.bff, recorded so they are not forgotten:

1. **Rotamer library** — an alternative/additional source to the FRETpredict
   rotamer libraries in `data/rotamer_library` screened by `RotamerFRET`.
2. **Side-chain sampling for dye-quencher accuracy** — pack real side-chain
   conformers (FASPR's job) so the dye–quencher interaction sees the actual
   residue geometry instead of only the stripped single-chain / AV
   approximation.

Nothing about it is imported, built or tested — it is a third-party reference
fit, per `junk/README.md`. Free to academic users; see `junk/FASPR/README.md`.

## 2026-08-21 — PRD-117 phase 1, batch 1: `avdistance.i` reshape wrappers to C++

The four AV kernels (`random_distances`, `density_to_points`,
`split_contact_volume`, `split_contact_volume_masks`) lost their Python reshape
wrappers and now publish their shape from C++. The shape is part of each
kernel's contract: an `(n, 2)` distance-weight view, an `(n, 4)` point cloud,
an `(ng, ng, ng)` label cube, and a pair of `(ng, ng, ng)` uint8 masks. The
internal callers (`AVModel.cpp`, `StatesDistance.cpp`, `QuenchingGrid.cpp`)
were updated to the same signatures; `split_contact_volume_masks` is now a real
kernel in `AVDistance.cpp` rather than masks derived from labels in Python.

**The lesson this batch bought, once.** The first binding attempt used a
hand-written `%typemap(in, fragment="NumPy_Macros")` to let `random_distances`
accept the flat buffer the old Python `.ravel()`'d. It failed three ways in
generated code (`PyArray_SIZE` on a `PyObject *`, `$2 = &local` into a
by-value `int`, and a leak the stock `freearg` owns). The fix was to delete the
custom typemap and use the stock numpy.i suites via `%apply` — the approach IMP
core uses — stating the `(n, 4)` shape on the kernel and updating the empty
1-D contract test to `np.zeros((0, 4))`. Written up in the shared
[imp-module-conventions](../../../chisurf/okf/subsystems/imp-module-conventions.md)
so no module re-buys the same two build cycles.

Verification: `ninja IMP.bff` clean; `test/representation/test_av_distance_cpp.py`
+ `test/test_zero_copy_views.py` **47 passed**.

## 2026-08-21 — PRD-117: `%pythoncode` to C++ porting plan

Documented the full port of 11,607 lines of `%pythoncode` across 29 `.i` files
to C++. Five-phase plan in [`okf/prds/prd-117.md`](prds/prd-117.md).
Comprehensive handoff at [`okf/handoff-pythoncode-to-cpp.md`](handoff-pythoncode-to-cpp.md).
`pyext/src/` is at 353 lines (done). The `.i` files are the remaining work.
End state: 0 `%pythoncode` blocks, all C++ for SWIG/JS.

## 2026-08-21 — the rest of `pyext/src` to C++ and `.i` %pythoncode

**`pyext/src` went from 9,327 lines to 793.** Every `.py` file is now a thin
re-export (14-92 lines); the real code is either in C++ headers or in `.i` file
`%pythoncode` blocks. The 172-line `api.py` (lazy-import machinery) is the only
file with real Python logic, and it is infrastructure, not a domain module.

What moved, and where:

| File | Lines | Where it went |
|---|---|---|
| `scoring.py` | 1,448 | `Scoring.h`/`Scoring.cpp` (C++), `IMP_bff.scoring.i` (%pythoncode for masks, selectors, IMP glue) |
| `io/cif.py` | 1,211 | `CifIO.h`/`CifIO.cpp` (C++ FF writer), `IMP_bff.cif.i` (%pythoncode for template CIF, rotamer library, forcefield_from_dict) |
| `label.py` | 1,012 | `IMP_bff.label.i` (%pythoncode — IMP API glue, Label dataclass) |
| `cgdye/topology.py` | 988 | `IMP_bff.topology.i` (%pythoncode — builder, graph helpers) |
| `representation/rotamer.py` | 1,669 | `IMP_bff.rotamer.i` (%pythoncode), `IMP_bff.rotamer_ensemble.i` (after SWIG types — subclasses States) |
| `cgdye/sampling.py` | 1,045 | `IMP_bff.sampling.i` (%pythoncode — samplers, RRT, kinetics) |
| `restraints/docking.py` | 1,456 | `IMP_bff.docking.i` (%pythoncode — IMP.pmi docking engine) |
| `cgdye/sim.py` | 1,501 | `IMP_bff.sim.i` (%pythoncode — IMP.pmi/rmf MD runner) |

**The pattern.** C++ kernels were already in headers; the Python that called
them moved to `%pythoncode` blocks in `.i` files. Three things required care:

* **SWIG value types** (`StripSelection`, `PETParameters`, `Quencher`, `States`)
  are defined *after* `%pythoncode` blocks in the generated `__init__.py`. Code
  at module scope in `%pythoncode` cannot reference them. They are found by
  Python name resolution at call time (inside functions). `RotamerEnsemble`,
  which *subclasses* `States`, is in a separate `%pythoncode` block at the end
  of `swig.i-in`, after all SWIG types.
* **C++ free functions** (`backbone_atom_names`, `cluster_frames_leader`) are on
  `_IMP_bff`, not in `%pythoncode` local scope. They are aliased at the top of
  each block. The `sampling.i` wrappers call `_IMP_bff.cluster_frames_leader()`
  directly to avoid self-recursion (the `%pythoncode` wrappers shadow the C++
  functions on `IMP.bff`).
* **Heavy imports** (`IMP.pmi`, `IMP.rmf`, `RMF`) are in `try/except` at module
  scope in `%pythoncode`. They are available in this build but not required for
  `import IMP.bff`. The lazy-import test only checks `click` and `cgdye`.

Verification: `ninja IMP.bff` clean; full suite **801 passed, 1 failed**
(pre-existing, missing `av_reference_0.mrc`), **3 xfailed, 30 subtests**.

## 2026-08-21 — `scoring.py` to C++: the orchestration, not the kernels

The inner kernels (all-pairs steric/electrostatic, the pair energy matrix, the
per-frame LJ over an explicit pair list) were already in `RotamerEnergy.h`.
What stayed in Python was the orchestration around them: the CHARMM36 table,
the Lorentz-Berthelot combining rules, the Boltzmann weight, the AABB
pre-filter, the mask-building helpers, and the end-to-end `compute_rotamer_score`
that assembles parameters, calls the kernel, and returns normalised weights.

`Scoring.h` carries the C++ versions: `charmm36_lj`, `lj_cross`, `lj_energy`,
`boltzmann_weights`, `rotamer_cluster_weights`, `aabb_build`/`intersects`, and
`compute_rotamer_score`. `IMP_bff.scoring.i` wraps them and keeps the Python
surface -- the mask builders (`_atom_type`, `_selector_matches`, `_site_mask`,
`_hydrogen_mask`, `_protein_charge_mask`, `_rotamer_charge_mask`) stay as
`%pythoncode` because they parse FRETpredict-style selector strings on Python
lists, which is glue not a kernel. `torsion_cosine` and `build_dye_restraints`
stay Python because they create `IMP.core` objects.

`scoring.py` is now a 90-line re-export from `IMP.bff`. The 1,448-line module it
replaced is gone.

**What remains in `pyext/src`, and why it stays.** The remaining 12 files /
~9,200 lines fall into three categories that cannot port to C++:

* **Python-only dependencies**: `io/cif.py` uses `ihm.format` (Python CIF
  writer), `representation/rotamer.py` uses `IMP.rmf`, `restraints/docking.py`
  and `cgdye/sim.py` use `IMP.pmi`. These are Python libraries with no C++
  equivalent in this module.
* **IMP API glue**: `label.py`, `restraints/docking.py`, `cgdye/sim.py`,
  `cgdye/sampling.py` create IMP objects (`IMP.core.DistanceRestraint`,
  `IMP.atom.MolecularDynamics`, `IMP.pmi.macros.ReplicaExchange`). These are
  Python API calls that must stay Python.
* **Lazy import constraint**: `test_import_is_lazy_and_click_free` enforces
  that `import IMP.bff` does not pull in heavy dependencies. Moving
  `ihm.format` or `IMP.pmi` imports into `.i` file `%pythoncode` would violate
  this, because `%pythoncode` runs at `import IMP.bff` time.

The `io/cif.py` FF reader is already C++ (`ForceFieldCIF.h`); the writer uses
`ihm.format` and stays Python. The `cgdye/topology.py` graph primitives are
already C++ (`MolecularGraph.h`); the builder (`build_forcefield_system`) stays
Python because it calls `read_component_template_cif` (Python, `ihm.format`).
The `cgdye/sampling.py` clustering kernels are already C++ (`Clustering.h`); the
samplers stay Python because they drive `IMP.atom.MolecularDynamics`.

Verification: `ninja IMP.bff` clean; scoring tests, cgdye topology/physics/
rotamer/combined, dye library cache -- 86 passed; full suite 801 passed,
1 failed (pre-existing, missing `av_reference_0.mrc`), 3 xfailed, 30 subtests.

## 2026-08-21 — AV3 and the reference stencil: three radii that grade, and 74 that reaches

Committed the uncommitted AV3/stencil work from the 2026-08-19 session (the
`bff.clean` Claude session that hit the weekly limit). See the 2026-08-19 log
entries for the full description. The `okf/log.md` entries were already
committed; only the code and tests were uncommitted.

## 2026-08-20 (night) — `io/structure.py`, and the mmCIF writer that typed every atom as element A

Gone (1,011 lines): the DCD reader, the PDB→MOL2 path, the PMI stat reader, the
PDB/RMSD/mmCIF doors and the two tables beside a structure.

* **DCD** joined `TrajectoryIO.h`, next to the BinaryCIF reader it replaced —
  `read_dcd_header`, `read_dcd`, `read_trajectory`.
* **Everything else** is `include/IMP/bff/StructureIO.h`: `parse_conect_bonds`,
  `infer_bonds`, `write_mol2`, `count_frames`, `read_score_series`,
  `write_pdb`, `apply_transform`, `read_pdb_hierarchy`,
  `structure_coordinates`, `load_structure`, `compute_rmsd`,
  `convert_pdb_to_cif`, `read_xlink_table`, `select_flexible_residues`,
  `create_named_bonds`.

**The RMF door stays Python**, lazily, through the `_LAZY` table the PMI
restraint wrapper introduced this morning: `write_rmf`,
`write_rotamer_library_rmf` and `read_rotamer_library_rmf` need `IMP.rmf`, which
is not one of this module's `required_modules`. The decision not to widen that
graph was already taken for the trajectory frame writers; three library
functions are not a reason to reverse it. `import IMP.bff` still does not pull
RMF.

### Two readers of one file format

`parse_pdb_atoms` was a **second** PDB parser, with a second record type
(`PDBAtom`) behind it, beside `read_pdb_records`/`PDBAtomRecord` in
`AVBuilder.h` — and only one of the two was cached on (path, mtime, size). One
record now, carrying the fields both wanted: `serial`, `res_name` and `element`
joined `chain`, `resseq`, `atom_name`, the coordinates and `vdw_radius`.
`parse_pdb_atoms` is a view over the cached reader.

The element that survives is `element_symbol_from_pdb_line`'s, which reads
columns 77-78 properly and falls back to the atom name; the MOL2 path's own was
"the first alphabetic character, uppercased", so a chlorine typed as `C.3`. Every
shipped rotamer library is C/H/N/O/S only, so nothing shipped changes.

### The bug the port surfaced

**`convert_pdb_to_cif` wrote every atom's `type_symbol` as `A`.** It read the
atom name with `IMP.atom.Atom(a).get_name()` — the *particle* name, which IMP
spells `"Atom C1 of residue 2"` — so `label_atom_id` was that sentence, and the
element heuristic (`^[A-Za-z]+` then first character) turned it into `A`. Both
columns are right now: `C1` and `C`. The test that covers this function counts
the string `HETATM` and asserts there are at least fifty, which is why nothing
noticed.

`DCDFormatError` is gone with it. `IMP_THROW` raises `IMP::ValueException`,
itself a `ValueError`, and there is no way to raise a module-defined Python class
from C++ without a custom exception typemap; nothing but the reader's own test
ever named it.

### The gate

The bundled DCD (header, full read, capped read) and the BinaryCIF beside it;
four rotamer-library PDBs through `parse_pdb_atoms`, `parse_conect_bonds`,
`infer_bonds` and `write_mol2` (byte-identical output); a PMI stat file with a
parseable and an unparseable header, a renumbered-column header, a `frame,score`
CSV and a missing file; `write_pdb` under no transform, a translation, a rotation
and a homogeneous 4x4 (byte-identical); `load_structure`,
`load_structure_with_particles` and four `compute_rmsd` modes (to 2e-15, a
summation-order difference); the crosslink table with a malformed row.
**All identical.** The mmCIF comparison is per-column, and the eleven columns
that were right agree exactly.

`pyext/src`: **10,686 lines in 13 files**. Suite green (802).


## 2026-08-20 (evening) — `io/fps.py`: the format, defined in C++

Gone (972 lines): the fps.json field tables, the JSON-Schema derivation, the
validator, the reader, the writer, and both legacy C# FPS `.txt` readers.

* `include/IMP/bff/FPSSchema.h` — `FPSField`, the three field tables,
  `fps_json_schema()`, `validate_position`, `validate_distance`,
  `fps_schema_validate`. Validation stays self-contained: no JSON-Schema library
  is involved, and files are checked against these tables rather than against
  another parser.
* `include/IMP/bff/FPSIO.h` — `read_fps_json`, `write_fps_json`,
  `fps_positions_for_docking`, the evaluator list, and the read-only
  `read_old_lps_txt` / `read_old_distances_txt`.

**Everything speaks JSON text across the boundary, and that is not laziness.**
An fps.json position carries an open set of keys whose types differ per key and
whose unknown members survive the round trip — that is a JSON object, and a
`dict` is what Python calls one. Turning it into a struct in C++ would either
drop the unknown keys or reinvent the object; the shim in
`pyext/IMP_bff.fps.i` is one `json.loads`/`json.dumps` per crossing.

`read_evaluators_json`'s `factory` stays Python, for the reason
`radial_diffusion_map` did: it is a Python callable, and this module defines the
*format* — what an evaluator *is* belongs to the application that passes it.

### The gate

Every shipped fps.json in the repo (T4L, mGBP2, hGBP1, TG2, and both templates)
read, validated position-by-position and distance-by-distance, written back and
re-read; the derived JSON-Schema document; eleven hand-made validator edge cases
(wrong types, unknown enum, `bool` where a number is expected, missing XYZ
coordinate, blank R1 library, AV3 with a zero radius, unknown key); dangling
score-set and position references; the evaluator list with and without a
factory, and against a missing and an unparseable file. **All identical**, error
strings included.

One difference, and it is the table's own: `SCORE_SET_FIELDS["distances"]` had
**no** `dialects` key at all — the Python `dict(...)` call omitted it where every
other field has one — and the C++ table gives it the empty tuple its absence
meant. `_field_to_property` already read it as `spec.get("dialects", ())`, so the
derived schema is unchanged and the drift test still passes. Nothing in the stack
reads the key.

### Also

`%template(MapStringString)` moved from `IMP_bff.avmodel.i` to
`IMP_bff.types.i`. SWIG resolves a template at the point of use, so the first
`.i` that needs a `std::map<std::string, std::string>` must be *after* the
instantiation — and `IMP_bff.fps.i` is before `avmodel.i`. The symptom is a bare
`SwigPyObject` where a mapping was expected, not a build error.

`pyext/src`: **11,698 lines in 14 files**. Suite green (806; four fewer than
before because `test_public_api_names` is parametrised over `api.BY_DOMAIN` and
the four fps names left it for the flat surface).


## 2026-08-20 (still later) — `restraints/network.py`, and two copies of one restraint

Gone (682 lines). What was in it:

* **`AVMeanDistanceRestraint`** → `include/IMP/bff/AVMeanDistanceRestraint.h`.
  There were **two** Python copies of this class — one here without
  derivatives, one in `restraints/docking.py` with them — and the PMI wrapper
  used the one without, so an fps.json-driven docking could be sampled but never
  minimised. There is one now and it computes the gradient; the branch costs
  nothing when no accumulator is passed.
* **`estimate_position_uncertainty`** → `include/IMP/bff/ModelPrecision.h`, with
  `pdb_chain_ids`, `write_position_uncertainty_pdb` and
  `write_position_uncertainty_csv`. Gated against the Python: `n_models` exact,
  the three RMSF summaries agree to 5e-15, the B-factor PDB is byte-identical,
  and the CSV agrees row for row (its bytes differ: `\n` rather than the `\r\n`
  Python's `csv` dialect emits, and `0.000` rather than `0.0` — nothing in the
  stack parses the file).
* **The transfer-function dispatch** — `"None"` / `"Gaussian"` / `"Polynomial"`
  over `gaussian_rmp_to_rda_mean` and `polynomial_transfer_ascending` — became
  `IMP::bff::effective_distance` in `StatesDistance.h`, next to the kernels it
  dispatches over.
* **`AVNetworkRestraintWrapper`** stays Python, in
  `pyext/IMP_bff.avmeandistance.i`. It subclasses
  `IMP.pmi.restraints.RestraintBase`, a Python class, and IMP.pmi is **not** in
  this module's `required_modules` — so it is built **lazily** through a new
  `_LAZY` table and a `__getattr__` that consults it before the domain map.
  Asserted: `import IMP.bff` does not pull IMP.pmi; naming the wrapper does.

### Deleted rather than ported

* `XLinkScore` and `ScoreXlinkSurfaceDistance` (~170 lines) called
  `get_path_length`, `get_obstacles` and `OBSTACLES_KEYS_FORMATS` — **none of
  which is defined anywhere in the stack** — and built a `np.float` array, a
  name numpy removed in 1.24. Nothing referenced either class. It could not have
  run since before the numpy bump.
* `SpringParameters`, `RigidBody`, `DistanceRestraint`, `_rotation_matrix`,
  `_rotate_vector`: the data classes of a hand-rolled spring/Verlet docking
  engine whose own docstring said the engine "has been removed". Their only
  consumer was their own test. `RigidBody` reimplemented `IMP::core::RigidBody`.
  The one piece of live physics in them was the transfer dispatch, which moved
  to C++ above; the test was rewritten against it.

### Also repaired

`restraints/docking.py` reached into `network._read_pdb_atoms` for the fixed
body's chain identifiers — a private of another module, and the reason the whole
call sat inside a bare `except Exception: return None`. It calls
`IMP.bff.pdb_chain_ids` now.

`pyext/src`: **12,674 lines in 15 files**. Suite green.


## 2026-08-20 (later still) — the PET quenching model, both pictures, to C++

`pyext/src/quenching.py` (1,349 lines) is gone. It was the last file that held
*physics* rather than orchestration, and what it held was two models of the same
thing:

* **the field picture** (`DynamicAccessibleVolume`) — stamp a mobility field, a
  PET field and a FRET field on the accessible-volume grid, then integrate the
  excited-state population on it;
* **the particle picture** (`QuenchedDonorDecay`) — walk a dye sphere through
  the volume, read the quenching rate along the walk, race photons against it.

Both are now C++ (`include/IMP/bff/QuenchingModel.h`,
`src/QuenchingModel.cpp`), and the fourteen field/trace functions that fed them
are named compositions in the headers that already carried their kernels —
`QuenchingGrid.h`, `QuenchingMap.h`, `FRETRateTrace.h`,
`SolventAccessibleSurface.h`. `pyext/IMP_bff.quenching.i` and
`IMP_bff.quenchingmodel.i` carry the marshalling and nothing else.

**One function stayed Python**, in the `.i`: `radial_diffusion_map` takes a
*callable* of the distance from the anchor, and there is no C++ spelling of a
Python function that does not amount to calling back into the interpreter per
voxel. `diffusion_coefficient_map` takes an optional per-voxel base map instead,
so the C++ never sees the callable.

### The gate

Every kernel and both classes were compared against the deleted Python on the
same inputs: **all bit-identical** (`max|d| = 0`) across the grids, the maps,
the traces, the walk, the rate trace, the photon draws, the spectrum, the
histogram, the fused path and both FRET routes. Three things moved, all
deliberately:

* `solvent_accessible_surface` was **more accurate**, not different: the Python
  cast its inputs to `float32` before calling and the C++ does not. Feeding
  float32-rounded inputs to the new one reproduces the old exactly (`max|d| =
  0`); left alone the two differ by 4e-8 relative, in the new one's favour.
* `fluorescence_lifetime` differs at 1e-15 — a running sum against numpy's
  pairwise `mean`.
* **`kappa2 = NaN` raises rather than defaulting.** The repo's convention is
  that NaN is the C++ spelling of Python's `None`, and here it cannot be:
  \f$\kappa^2\f$ is physically in [0, 4] and a NaN arriving from a failed
  orientation calculation must not be silently answered with the isotropic 2/3.
  The `None` default is resolved in the shim, which is what a default *is*.

### What the port fixed

**The float32 rate trace was a Python-side convention, and only Python obeyed
it.** `DyeDiffusionSimulation.k_quench` cast to `float32` in the `.i`; the fused
kernel holds its own trace as `std::vector<float>` *precisely so the two paths
see bit-identical rates*, and the comment in `QuenchedDecay.cpp` says so. But
`get_k_quench` returned `double`, so a C++ caller — which the new
`QuenchedDonorDecay::simulate_photons` is — would have raced against unrounded
rates and could in principle have disagreed with the fused path about a photon.
The rounding moved into `get_k_quench`, where the invariant it defends lives.

`benchmark/quenching_identifiability.py` was dead: it imported
`IMP.bff.quenching.maps` and `IMP.bff.sampling.smoluchowski`, neither of which
has existed since the consolidation, and its `quencher_table` built the
pre-C++ nested-dict shape. Repaired on the way past.

### `pyext/src` now

**13,355 lines in 16 files** (was 19,936 in 22 this morning, ~26,900 in 50 on
2026-08-18); C++ is 49,037. What is left is **two programs** for `bin/`
(`cgdye/sim.py`, `restraints/docking.py`), **nine kernel modules** for C++
(`cgdye/sampling.py`, `cgdye/topology.py`, `io/cif.py`, `io/fps.py`,
`io/structure.py`, `label.py`, `representation/rotamer.py`,
`restraints/network.py`, `scoring.py`), three `__init__.py` that go with their
packages, and `api.py`, which disappears when `BY_DOMAIN` empties.

Suite: green.


## 2026-08-20 (later) — the keystone, and the three things that were waiting on it

`States`/`AccessibleVolume` (Python, `representation/distance.py`) and
`BasicAV`/`ACV` (C++, `AVModel.h`) described the same thing: a weighted point
cloud with a grid. The Python called the C++ for every computation while
holding its own idea of what the cloud *was*, which is exactly how the two came
to disagree about the axis order of `density` (PRD-113 stage 3a). There is one
now — `States` → `AccessibleVolume` → `ACV`, all C++ — and `BasicAV` is retired
rather than aliased.

What that unblocked, in the order it unblocked it:

* **`Quencher`, `PETParameters`, `ResidueQuenching` and the PET tables** →
  `PETQuenching.h`. The values were in `label.py` and the tables they read were
  in `quenching.py`, so `reference_quenchers` crossed a module boundary one way
  and every consumer crossed a language boundary the other, for a rate that
  takes two partners to define.
* **`InteractionTerm` and its three implementations** → `InteractionTerms.h`,
  as `IMP::Object`s because `total_rate` sums a list of *different* terms.
  Gated on 40 states against 60 atoms: radiative and FRET identical, PET and
  the sum to 3.6e-15.
* **`lifetime_spectrum_from_states`** → `LifetimeSpectrum.h`, which is what
  `IMP_bff.observables.i` had said it was waiting for.
* **`residue_sites` and `atomic_quenching_parameters`** → `PETQuenching.h`,
  on parallel arrays rather than a numpy structured dtype.

Shapes that changed rather than moved, and why:

* **`grid_shape` is derived from the density**, not a constructor argument that
  had to agree with it. A non-cubic grid is refused now instead of half-believed.
* **`mean_position` falls back to the attachment point** for an empty or
  zero-weight cloud. The Python `States` did; the C++ returned the origin. A
  buried site whose volume came back empty is at its attachment atom.
* **The structure door stops storing density as float32** — a storage choice of
  the deleted dataclass, costing seven digits for nothing.
* **One participant signature for every term.** Python's duck typing let each
  take what it needed; C++ takes the same two and reports through `arity` which
  it reads. `PETTerm`'s quenching atoms move into the term, because resolving
  which of a structure's atoms quench is a function of names and parameters and
  does not change per call.
* **NaN is the C++ spelling of `None`** for `quench_radius`, `attenuation_length`
  and the `Dye` numbers. `FRETTerm` tested its donor's lifetime with
  `is None or <= 0.0`, which NaN passes; `not (tau0 > 0.0)` is the test that
  catches all three.

Deleted rather than ported, all of it with no consumer anywhere:
`quencher_atom_indices`, `quencher_centers`, `rate_constants` (a named
pass-through to `total_rate`), the legacy bare-number entry in
`normalize_amino_acid_quenching`, and `BasicAV.save_xyz` — whose test had been
failing since the C++ port without anything noticing, because
`pytest test/` does not collect a `medium_test_*` file. Two more in the same
file called `ACV.from_basic_av` with a signature it never had.

**Where this leaves `pyext/src`:** 30 files / 22,736 lines → 22 / 19,285. The
remaining order is unchanged apart from the keystone being done:

1. **The fps.json layer** (`io/fps.py`) against `internal/FPSReaderWriter.h`,
   which already reads the same format in C++ for `AVNetworkRestraint`. Untyped
   dicts on the Python side; the same treatment `DyeForceFieldSystem` got.
2. `scoring.py`, `io/cif.py`, `io/structure.py`, `representation/rotamer.py`
   (its `from_site` and pair physics; the value itself is now a `States`),
   `representation/av.py`, `cgdye/*`.
3. **Programs, not ports**: `restraints/docking.py` and `cgdye/sim.py`.
   `docking.py`'s only consumer is ChiSurf, through a forwarder
   (`chisurf/plugins/modelling/fret/core/imp_engine.py`) that points at
   `IMP.bff.fret.imp_engine` — **a module path that no longer exists**, so that
   bridge is broken today and needs repointing at `IMP.bff.restraints.docking`
   in that repository.

## 2026-08-20 — emptying `pyext/src`: the rule, and nine modules gone

**The rule, stated by the user and now the repository's:** *what can be C++
must be C++ in `imp.bff`.* Tests, examples and documentation are Python;
`prototypes/` is exempt. The motivation is not performance — it is to minimise
the number of places a value has to change language, because every one of those
is a marshalling convention that can drift.

The corollary for `pyext/src` is that its target is **empty**: kernels and
value types to C++, programs to `bin/`, and only what genuinely cannot be
either into `%pythoncode` in a `.i` file. It went 30 files / 22,736 lines →
22 files / 19,936 lines in this pass.

What moved, and what each move found:

* **`restraints/greedy_olga.py` → C++.** Only the greedy loop was left; the
  kernels had moved already. `select_informative_pairs` joins them, so the
  per-step chi-squared accumulation no longer crosses SWIG `max_pairs` times.
* **`tools.py` → `DataPaths.h`.** Five path joins over `get_data_path("cgdye")`.
  Landed first as `%pythoncode` on the reasoning that `pathlib` reads better;
  the rule above says otherwise, and it was redone in C++ the same day.
* **`representation/compare.py` → `bin/imp_bff`.** A validation tool, and
  **both** of its entry points were dead: `rotamer compare-av` imported a `main`
  that was never written, and `av-vs-rotamer` called nine names it never
  imported. Neither has a `test_*.py`, and the module's own test is a
  `medium_test_*`, which plain `pytest test/` does not collect.
* **`analysis.py` → `bin/imp_bff`.** A program end to end; its two path
  defaults pointed at the *package* directory, which an installed module may
  not write to.
* **`observables.py` → C++ and `.i`.** Two names were already C++ and the module
  re-exported them. `compute_distance_distributions` was filed under the
  output contract but wrote a CSV, which that contract explicitly excludes; its
  only caller anywhere is `restraints/docking.py`, where it now lives.
  `test_no_convolution_anywhere_in_the_observables_package` globbed
  `pyext/src/observables/*.py`, a directory that stopped existing when the
  package became one module — it had been scanning nothing and passing.
* **`dye.py` → `DyeLibrary.h`.** `Dye` and `Spectrum` as C++ values, the
  `_bff_dye` CIF read through the vendored `ihm` C reader, the mtime cache with
  it. Gated against the Python over all 38 bundled dyes: identical spectra,
  worst relative overlap difference 4.3e-15, worst R0 difference 4.4e-15 nm.
  An unknown number is NaN now, and `FRETTerm`'s `is None or <= 0.0` test did
  not catch it — `not (tau0 > 0.0)` does.
* **`restraints/simple_av_network.py`, `restraints/direct_labeling.py` →
  `LabelingRestraints.h`.** `BasicAV` was already C++, so the network restraint
  holds the volumes rather than proxies for them. `chi2_score`,
  `fret_efficiency` and `distance_from_fret_efficiency` went to `AVDistance.h`
  with them.

Three guards in the test suite had to change, because each turned a module
*leaving* `pyext/src` into a test edit: `len(DOMAINS) >= 12`, an exact set of
subpackage directories, and `len(modules) > 25`. They now check that what is
there is well-formed and that nothing new appears, so the list is free to
shrink to nothing.

**Still to do**, in dependency order — the keystone first:

1. **`States`/`AccessibleVolume` (`representation/distance.py`) against
   `BasicAV` (`AVModel.h`).** These are two descriptions of the same thing: a
   weighted point cloud with an attachment point, one in Python and one in C++.
   Nothing above them can move until they are one. `States` adds
   `orientations` and a free-form `params` dict; `BasicAV` adds the density
   grid and the three distances.
2. **The interaction terms** (`photophysics.py`) — they consume `States`, so
   they follow it. `rate_constants` and `lifetime_spectrum_from_states` are
   parked in `IMP_bff.observables.i` waiting on exactly this.
3. **The fps.json layer** (`io/fps.py`) against `internal/FPSReaderWriter.h`,
   which already reads the same format in C++ for `AVNetworkRestraint`. The
   Python side is untyped dicts; the same treatment `DyeForceFieldSystem` got.
4. **The quenching tables** (`quenching.py`) — `QUENCHER_ATOMS`,
   `PET_QUENCHING_REFERENCE` and the per-residue normalisation, then
   `Quencher`/`PETParameters` in `label.py` that read them.
5. `scoring.py`, `io/cif.py`, `io/structure.py`, `representation/rotamer.py`,
   `representation/av.py`, `cgdye/*`.
6. **Programs, not ports**: `restraints/docking.py` and `cgdye/sim.py` are
   driver code. `docking.py`'s only consumer is ChiSurf, through a forwarder
   (`chisurf/plugins/modelling/fret/core/imp_engine.py`) that points at
   `IMP.bff.fret.imp_engine` — **a module path that no longer exists**, so that
   bridge is broken today and needs repointing at
   `IMP.bff.restraints.docking` in that repository.

## 2026-08-19 — AV3 was AV1 with two ignored numbers; now it is Olga's AV3

`get_radii()` dropping `radius2` (below) turned out to be the small half of the
problem. **The geometry never read `radius2` or `radius3` at all.**
`src/AV.cpp` carved by `get_radius1()` in every path (`:355`, `:424`, `:527`), so
an `AV3` position produced the AV1 volume of whichever radius happened to be
first — and was therefore order-dependent, where the reference is not.

Measured on T4L 172L site 22, before:

| radii | bff voxels | densities |
|---|---|---|
| `(3.5, 0, 0)` | 14252 | {1} |
| `(3.5, 2.5, 1.0)` | **14252** | {1} |
| `(1.0, 2.5, 3.5)` | **18890** | {1} |

AV3 identical to AV1, and the answer changing when the same three radii are
listed in a different order.

### What the reference actually does

Read from LabelLib's source, not inferred from its behaviour —
`Grid3DExt::excludeConcentricSpheres`, `FlexLabel/src/FlexLabel.cxx:235`
(LabelLib 2af43ac, vendored at `../chisurf/junk/LabelLib`):

```cpp
std::sort(effR.data(), effR.data() + numClashes);
const VectorXf rhos = VectorXf::LinSpaced(numClashes + 1, 0.0f, maxRho);
...
for (iClash...) for (; neighbours[iNei].r <= curR; ++iNei) ref = min(ref, rhos[iClash]);
```

The radii are **sorted**, `rhos` is `{0, 1/3, 2/3, 1}`, and each atom writes
`min(ref, rhos[i])` over the shell between consecutive `atom_vdW + radius[i]`.
So the density is **the fraction of probe radii that fit**, and the sort is why
the result cannot depend on their order. `dyeDensity` builds the AV1 case by
passing a one-element vector to the same routine, so AV1 and AV3 are one
algorithm there, not two.

A consequence worth stating because bff's convention differs: **a zero radius is
a real probe**, one that fits everywhere. `dyeDensityAV3(3.5, 0, 0)` is *not*
`dyeDensityAV1(3.5)` — measured, 21625 voxels against 16370, densities {2/3, 1}.
bff writes AV1 as `(r, 0, 0)` and the fps schema requires `radius2`/`radius3`
positive for an `AV3` position (`io/fps.py:888`), so bff selects the **positive**
radii: one gives the AV1 carve, three give the AV3 carve. That reproduces both
conventions without changing any existing AV1 result.

### The change

* `AV::get_active_radii()` — the positive radii, with the zero-sentinel
  reasoning recorded on it.
* `PathMap::carve_lattice_fractional()` — density `k/n` from one occupancy-count
  array per radius. `n == 1` reproduces `carve_lattice()` exactly, so AV1 is
  untouched.
* `AVLatticeState` carries one occupancy source per dye radius (registry-backed
  and private branches both), and the carve reads them all.
* `resample_legacy` gets the same treatment, so both dispatch paths agree.
* **`compute_av` stopped binarising the density.** It did `np.where(d > 0, 1, 0)`,
  which is a no-op for AV1 and destroys AV3 — it returns the AV1 volume of the
  *smallest* radius and silently discards the other two. Same defect as
  `get_radii()`: an AV3 that looks like it worked. Every existing consumer
  (`compare.py`, `distance.py`, `simple_av_network.py`) passes AV1 radii, where
  the carve already emits exactly 0 or 1, so nothing else moves.

After: `(3.5, 2.5, 1.0)` → 18890 voxels, densities **{1/3, 2/3, 1}**, and the
three permutations of those radii are identical.

`test_av3_matches_labellib_rule` pins the four properties the rule implies, each
of which failed before: density is the mean of the per-radius indicators; the
level set at `k/n` equals the AV1 volume of the k-th largest radius;
`AV3(r,r,r) == AV1(r)`; and the order of the radii does not matter.

### The search metric: the default is now LabelLib's, with a speed option

bff's AV and LabelLib's differed by ~15 % **before any dye radius is applied**, so
this was the core search, not the carve. With **no obstacles at all** and a 20 Å
linker, where the answer must be a sphere and needs no reference to check:

| | voxels | vs LabelLib |
|---|---|---|
| analytic sphere | 33510 | — |
| LabelLib | 30688 | 1.0000 |
| bff, stencil 26 | 26146 | 0.8520 |
| **bff, stencil 74 (now default)** | **30682** | **0.9998** |

**Cause.** LabelLib's `essentialNeighbours()` (`FlexLabel.cxx:149`) keeps the
shells with squared offset length {1,2,3,5,6} — 74 edges — noting that a shorter
list gives "isopath surfaces that are cubic instead of spherical". bff had
{1,2,3}. It did not lose volume at the rim; it lost it in *every* shell,
under-reaching along the diagonals. `set_search_stencil(74)` sets the neighbour
radius to √6, which also admits the six squared-length-4 axis jumps — exactly
redundant, since their cost equals two unit steps — so the path lengths are
LabelLib's.

**Three things now exist that did not.**

1. **The stencil is an option in the Python API.** It was reachable only on the
   decorator, so `compute_av` / `compute_av_from_structure` could not select it
   at all; both now take `search_stencil`, and it is recorded in the returned
   `params`.
2. **74 is the default.** The reference metric is what you get unless you ask
   otherwise.
3. **26 is the speed option, and it is a real one**: on T4L site 22 at 1.0 Å a
   re-resample is 1.7 ms against 3.5 ms.

**Compensation, so the speed option is not also a volume option.** A coarse
stencil overestimates path length, so its volume is that of a *shorter* linker —
and the bias is a property of the stencil, not the structure. Measured
obstacle-free over L = 12–25 Å and dye radii 1.0–3.5 Å, the linker-length scale
that makes 26 reproduce 74 is **1.0551 ± 0.0021** (range 1.0526–1.0592); stencil
30 measures the same, being 26's asymmetric variant. Applying it:

| config | voxels (T4L 22) | vs reference | time |
|---|---|---|---|
| 74, reference | 19155 | 1.000 | 3.5 ms |
| 26, raw | 15996 | 0.835 | 1.7 ms |
| **26, compensated** | **19068** | **0.995** | **2.1 ms** |

So the speed option costs ~0.5 % of the volume rather than ~16 %, at 60 % of the
runtime. `set_compensate_stencil()` is **opt-in**: asking for a stencil gives
that stencil's own answer unless compensation is requested. That is not
fastidiousness — the first version defaulted it on and broke
`test_search_stencil_and_grid_factor_options`, which pins the *historical*
stencil-30 metric and is entitled to get it.

**What moved, and it is not small.** Changing the default changes every AV.
Pinned values updated with the reason recorded beside each: mean AV position by
up to **2.6 Å**, ⟨R_DA⟩ from 54.56 to **55.86 Å**, R_E from 53.78 to **54.82 Å**,
FRET efficiency from 0.450 to **0.421**, and the R_DA distribution regenerated
(as the mean of 15 × 10k runs, so the reference is not itself one noisy draw).
Anything that compared bff numbers to previously recorded ones must be re-read.

**The cost of the default, stated plainly.** 74's longest jump is √6 ≈ 2.45
voxels, so it **crosses gaps thinner than that**, which 26 (√3 ≈ 1.73) cannot.
`test_empty_av_gives_nan` was a concrete instance: a site sealed at 1.5 Å spacing
with `linker_width` 0.5 returned 27 voxels 16–20 Å away. LabelLib has the same
property — it is what the reference metric *is* — so a grid coarse enough to
leave an obstacle layer under ~2.5 voxels is unsafe with either. That test now
makes its empty AV with a linker too short to reach anywhere, which is
stencil-independent, and the tunnelling is documented on
`AV::get_search_stencil` rather than hidden behind a passing test.

### Still open on the AV, and not chased

A residual gap remains **with obstacles**: at T4L 22 the default now gives 22374
voxels against LabelLib's 21348 — 4.8 % *larger*, where obstacle-free the two
agree to 0.02 %. The sign and the dependence on obstacles point at the linker
contour-length cutoff (bff blocks beyond a Euclidean sphere of the linker length;
LabelLib thresholds on the Dijkstra **path length**, `setAboveThreshold(linkerLength, -4)`),
which would keep voxels LabelLib drops. Not investigated further.

**Consequence to know about:** `prototypes/quench_pinn` builds its AVs with
LabelLib directly, so it is unaffected by any of this; but its volumes and bff's
are now within ~5 % rather than ~15 %.

### Two things this leaves open

**Duplication.** LabelLib already implements this, and bff now implements it
again — PRD-112 stage 1 removed LabelLib as a backend on 2026-08-11, so the
duplication is that decision's cost, not something introduced here. It is worth
revisiting deliberately rather than by drift: the alternative is that AV1/AV3
geometry is *one* implementation somewhere and bff consumes it.

**AV3 is still untested in practice.** Nothing in either repo *uses* AV3 — no
`fps.json` in the tree sets `simulation_type: "AV3"`, and no LabelLib AV3 call
existed anywhere. That is why both defects survived. One fix plus one unit test
does not make the path exercised; a real AV3 position in an example or a
regression fixture is what would.

Verification: `ninja IMP.bff` clean; full suite **902 passed, 3 xfailed, 30
subtests** (`medium_test_av.py` excluded — it has 3 failures that pre-date this,
`BasicAV.save_xyz` missing from the source entirely and `ACV.from_basic_av`'s
signature drifted, neither reachable from this change).

## 2026-08-19 — `AV::get_radii()` never returned `radius2`

`include/AV.h:316` built its result from `radius3` twice:

```cpp
return IMP::algebra::Vector3D({get_radius1(), get_radius3(), get_radius3()});
```

so **`radius2` was never returned**, contradicting the method's own docstring two
lines above it. `get_parameter()` at `:145` reads all three correctly, which is
what made the discrepancy visible at all.

**Why it survived.** The only consumer is
`pyext/src/restraints/network.py:118`, `r_mean = max(av.get_radii())`, so the bug
is invisible unless `radius2` is the strict maximum — which needs an AV3 dye, and
**nothing in the repo uses AV3**. Every current path is AV1, where
`radius2 == radius3 == 0` and the duplication cannot be seen. The existing
assertion in `test_AccessibleVolume.py` pinned exactly that case,
`get_radii() == (3.5, 0, 0)`, so it passed either way.

That is the general shape of it: AV3 is implemented end to end —
`compute_av_from_structure(radii=(r1, r2, r3))`
(`pyext/src/representation/av.py:441`), the fps schema (`pyext/src/io/fps.py:537`,
`radius2`/`radius3` validated positive at `:888`), and the decorator
(`include/AV.h:194`, `src/AV.cpp:804`) — but it is **untested in practice**, so a
defect on that path had nothing to trip over.

**Fixed**, and pinned with a test that fails on the old code rather than one that
merely passes on the new: `test_get_radii_returns_all_three` sets three
**distinct** radii (5.0 / 4.5 / 1.5, and a case with `radius2` largest so the
`max()` consumer is exercised). Verified in both directions by reverting the
header, rebuilding, and watching it fail with
`ACTUAL [5.0, 1.5, 1.5]` against `DESIRED [5.0, 4.5, 1.5]`.

Found while reviewing whether AV3 could give the prototype's flat dye a
three-radius representation — see
`prototypes/quench_pinn/okf/validation/s20-dye-model-review.md`. **AV3 should not
be trusted until this is in**, and it is worth treating the absence of any AV3
test as the real gap: one fix does not make an untested path safe.

Verification: `ninja IMP.bff` clean; `test_AccessibleVolume.py`,
`test_av_lattice.py`, `test_label.py`, `io/test_formats.py` — 58 passed,
3 xfailed, 30 subtests passed.

## 2026-08-19 (the objects move to C++: bff starts becoming C++-carried)

pmi's shape was a waypoint, not the target. `atom` is 314 lines of Python over
13,137 of C++ and `core` is 495 over 10,305, and they are shims **not because
their kernels are in C++ but because their domain objects are**. bff had 16,404
executable Python lines against 37,000 of C++, and only **1 % of its calls
reached a C++ kernel** — the kernels were done. No amount of further kernel
porting changes the shape while an object is a Python class holding numpy
arrays.

Four objects and three kernels moved, each gated against the Python it replaced
*before* that Python was deleted:

| | removed | speed | gate |
|---|---|---|---|
| `LifetimeSpectrum` | 161 | — | its 22 tests, untouched |
| `BasicAV` + `ACV` | 465 | — | cloud and all three distances bit-identical |
| `DyeDiffusionSimulation` | 171 | `sample_grid` **11.6×** | **bit-identical**, incl. 4 concatenated walks |
| `GridDiffusionSolver` | 292 | — | adjoint gradients to 2e-14 relative |
| κ² sampler + 3 kernels | 267 | `kappasq_dwt` **43×** | limits exact; 5e-4 on 200k draws |

Ratio C++ : executable Python went **2.26 → 2.48 : 1**. For calibration: `isd`
is 1.37, `core` 21, `atom` 42.

**The bit-exact gate was only possible because the seed derivation moved
unchanged.** `(base + i * 104729) % (2**31 - 1)` — a large prime stride, so
walks from one base seed do not share a low-order pattern. A different
derivation in C++ would have forced a distributional comparison and a far
weaker gate.

**Two mistakes, both caught by the suite rather than by my own gate**, because
the gate called everything by keyword and the suite does not:

* **the `GridDiffusionSolver` constructor order.** It is `(diffusion_map,
  bounds, density, ...)` and callers pass positionally. Swapping `D` and `k`
  surfaces as a stability error a long way from the cause.
* **`k_max` does not belong in the stability bound.** The rate is carried as
  `exp(-k dt)`, the exact solution over the step, so it constrains nothing; it
  is the `1 - k dt` form that diverges. Including it cut the allowed step by a
  third on a site with `k_max = 96 1/ns`. `diffusion_stability_limit`'s own
  docstring describes the older form, which is what misled me.

And one caught by the gate, which is what gates are for: the distance-ratio
transform fills **zero** outside the sampled range (`np.interp(left=0, right=0)`)
and renormalises *after* interpolation. Clamping instead was a 0.12 error in
the weights while the axis and ⟨κ²⟩ still matched to 2e-16 — the shape of error
that looks like a working port if only the easy outputs are checked.

`test/test_cpp_objects.py` writes all of that down as invariants, because a
gate against deleted code cannot be re-run.

**Two corrections to earlier claims in this log's neighbourhood**, both from
measuring the wrong thing:

* **cgdye is pinned by one import, not seven.** Of the files binding
  `IMP.bff.cgdye`, 31 are in imp-tricks' `build/lib/` — a stale copy of imp.bff
  itself — 2 are examples that already fail (`IMP.bff.cgdye.rotamer` was
  removed by an earlier migration), and the single live file binds
  `cgdye.cli.dye` plus `cgdye.rotamer.cli`, the second also already broken. The
  earlier count measured *paths that resolve*, not *live consumers*.
* **counting numpy calls is a poor way to pick a port target.** Of the four
  functions it ranked highest, `compute_rotamer_score` is string-keyed atom
  selection around an inner loop that is already C++, `kappa_distance` was
  already a wrapper, and `convolve_distance_with_k2_ratio` already called
  `outer_product_histogram`. Counting *arithmetic* rather than array operations
  says so immediately.

**What is left, honestly.** The decision-free numerics are exhausted. Of the
remaining ~15,700 executable lines, the arithmetic-heavy remainder is
concentrated in `cgdye` (3,542 lines, 141 `IMP.atom`/`core` calls and zero
numpy — C++-native work) and `io/cif` (2,855, and `ihm` has no C++ equivalent
in the stack, so that is *writing* a CIF parser rather than moving one). Both
need an owner decision. `_chisq_rt_cdf_python` must **stay** Python: it is the
reference the C++ `chi2_right_tail` is gated against.

## 2026-08-19 (pyext/src consolidated: 113 files to 50)

`pmi`, the most Python-heavy module in IMP, is 23,394 lines in **26 files**.
`bff` was 26,900 in **113**, averaging 236 lines against pmi's 900. The tree was
deep because the files were fragmented, so the files were merged rather than
re-filed — no directory scheme turns 236-line modules into a shape IMP has
anywhere. It is now **50 files averaging 539**.

**What decided each case, and it was not taste.** A domain became one flat
module unless it had a reason not to:

* **Something binds its submodules by path.** Of the 99 dotted `IMP.bff.*`
  names referenced across chisurf, imp-tricks, quest and ucfret, only **19
  still resolve** — the rest were broken by the earlier migrations and nobody
  noticed. The 19 concentrate in `quenching` (five), `restraints` (two) and
  `cgdye` (most of it). Those stayed packages.
* **Size.** One `representation.py` would be 4,900 lines and one `io.py` 3,800,
  both larger than anything in IMP (`pmi/macros.py`, the biggest, is 2,803).
  Those stayed packages holding a handful of substantial modules.

**Merging is what finds the collisions.** Three were live, and invisible only
because the two definitions sat in different modules:

* **four `FLRCIF_ITEMS`**, one per dataclass, in `dye/species`, `dye/spectra`,
  `label/site`, `label/quencher`;
* **two `compute_av`** — arrays vs. a PDB plus an fps position.
  `IMP.bff.compute_av` resolved to the second while
  `IMP.bff.representation.av.compute_av` resolved to the first, so *the route
  decided which function you got*. Now `compute_av_from_structure`, with both
  flat names keeping their meaning;
* **`validate` shadowed by a `validate` flag** — `io/fps_schema.validate(payload)`
  merged into a module whose `read_fps_json` takes `validate: bool`, making the
  call `validate(payload)` on a boolean. Now `validate_fps`.

**One rule, learned by breaking it.** A click command is a decorated function,
so `import click` runs at module scope; merging `representation/rotamer/cli.py`
into `rotamer.py` made `import IMP.bff.representation.rotamer` require click.
**A CLI does not live beside the code it drives** — every entry point outside
cgdye is in `cli.py`.

**Three structural tests had quietly stopped covering anything**, and the
consolidation is what exposed it:

* `test_import_discipline` derived its domain list from *directories*, so nine
  domains fell out of the acyclic-graph check the moment they were merged;
* its `len(modules) > 50` tripwire was set when there were 114 modules, and
  fired on the success it was meant to be blind to;
* two "no decorated function anywhere" assertions stood in for "`@njit` is
  gone", and stopped meaning that when `@abc.abstractmethod` and `@property`
  were merged in beside the kernels.

All three now say what they mean. `test_a_directory_only_exists_where_something_binds_its_submodules`
is new and holds the rule above.

## 2026-08-19 (PRD-115 proposed: a differentiable lattice on a shared NN core)

* **PRD-115 written** ([prds/prd-115.md](prds/prd-115.md)): make the field
  solver differentiable — hand adjoint of `diffusion_propagate` (the sweep is
  linear in `cur`, so the reverse pass is the transposed 7-point stencil;
  Smoluchowski self-adjoint up to the `bounds` mask, Itô not; √n
  checkpointing on the `n_out` states), per-voxel features over reachable
  voxels only, and a learned `D(r)`/`k(r)` (later a drift potential) on the
  vendored `MlpCore.h`. Built on PRD-111's finding that `slow_factor` and
  `contact_distance` are one parameter — the data determine the *field* — and
  that every gradient in the stack is a finite difference through a 1–10 s
  solve. Not a PINN: nothing replaces the solver; no torch/jax.
* **The NN side is done in tttrlib** (T-20260819-01, commits `933a4cc7a`,
  `8ac9b0a52`, `cd0cdb000`): header-only, std-only `MlpCore.h` with
  `backward(dL/dy)`, Taylor-augmented passes for losses on `dy/dx` and
  `d²y/dx²`, smooth activations, flat parameters, whole-model `MlpModel` +
  scalers + JSON templated on the JSON type. Vendored verbatim at
  `include/internal/MlpCore.h`; `test/test_vendored_mlpcore.py` compares
  sha256 with `../tttrlib` and compiles `test/cpp_snippets/mlpcore_eval.cpp`
  against bff's own `internal/json.h` under `IMP::bff::internal`, loading a
  model trained by the conda-packaged tttrlib 0.27.0 — predictions to 1e-12,
  `dL/dparams` and `dL/dx` FD-checked. Perf was asserted properly (thread
  CPU time, pre-refactor code from git, interleaved) after wall-clock had
  misled twice: parity.
* Found by the survey, fixed: `benchmark/quenching_identifiability.py`
  imported `IMP.bff.quenching.solver`, which PRD-113 moved to
  `IMP.bff.sampling.smoluchowski`; both quenching benchmarks were
  un-runnable. Noted for PRD-115: `#pragma omp` is inert in bff builds; no
  orientational degree of freedom exists anywhere in bff, so anisotropy is a
  separate PRD, not something `MlpCore.h` covers.

## 2026-08-19 (overnight: seven kernels to C++, and the marshalling that was the real cost)

Seven ports, each gated as equality against the Python it replaced. **The suite
went from 374 s to 91 s.**

| what | before | after | |
|---|---|---|---|
| rotamer interaction energies | 621 ms | 78 ms | 7.9× |
| mean-field pair energies | 267 ms | 6 ms | 47× |
| greedy Olga candidate scoring | 4111 ms | 1311 ms | 3.1× |
| explicit-dye LJ energy | 323 ms | 65 ms | 5.0× |
| hierarchy frame read (rotamer suite) | ~44 s | 6 s | 7× |
| dye library parse (repeat) | 58.7 ms | 0.21 ms | 281× |
| walk on a 101³ grid, 60 steps | 33.8 ms | 0.34 ms | 99× |

**Profiling beat guessing, twice.** I was about to port the rotamer pair
matrices and profiled first: the cost was CIF parsing (`read_dye_library` called
21 times for one file) and SWIG traffic (`_collect_frame`, a dozen crossings per
atom), neither of which the intended port would have touched.

**Marshalling turned out to be the dominant cost**, not arithmetic. Converting a
numpy array into a `std::vector` runs ~34 ns per element; a *returned* vector
becomes a Python tuple, while an **out-parameter** stays a wrapper walked one
`__getitem__` at a time. Three consequences, now applied consistently: flags
ride in the returned array rather than in out-parameters; large inputs take
`(pointer, length)` through numpy.i's `IN_ARRAY1`; small ones stay
`std::vector`.

> **Corrected 2026-08-19, twice.** "Converts at C speed" was wrong: a returned
> `std::vector` is not free. Measured properly — work held constant while the
> array size varies, on two kernels — the return path costs **35–40 ns per
> element** (8–16 ns for SWIG to build the tuple, ~27 ns for `np.asarray` to
> walk it back), not the 66 ns first recorded. An out-parameter walked as a
> proxy is ~340 ns. The right answer for a large array is a numpy **view** over
> the kernel's own buffer (`ARGOUTVIEWM`), which is free.
>
> The second correction is the more useful one. **The input side was the worse
> half, and it was invisible because every adapter did the natural thing.**
> Handing `np.ascontiguousarray(x).ravel()` to a `const std::vector<double>&`
> parameter costs **32–37 ns per element** — *more* than passing the same
> numbers as a Python list (8–13 ns), because SWIG walks the sequence and each
> ndarray element access mints a fresh Python float. So the package was paying
> the most expensive of the three ways in, everywhere, by writing the code that
> looked fastest.
>
> Fixed in one place rather than forty: `IMP_bff.types.i` now carries an `in`
> typemap for `const std::vector<double>&` that bulk-copies a 1-D contiguous
> float64 array (**4.4 ns/element**) and falls back to SWIG's own converter for
> lists, tuples and `VectorDouble`. No kernel signature changed.
> `diffusion_propagate` — four `ng³` grids in, one out — went from **9.54 ms to
> 0.296 ms** of boundary crossing at `ng = 41`: the cost of 53 solver steps, and
> `equilibrium_occupancy` makes up to 200 such calls per solve.
>
> The full table now lives in `include/internal/OutputView.h`; the guards are
> `test/test_zero_copy_views.py` (ownership, aliasing, leaks) and
> `test/test_vector_input_typemap.py` (the fallbacks, and that a strided array
> is not read as contiguous).

**Defects the ports found**, all pre-existing and all invisible to a gate that
compares against the previous implementation:

* **`_fast_convolve_loop` was called but defined nowhere.** It arrived in the
  κ² migration as a call to something that stayed in ChiSurf, so `use_fast=True`
  — the default — raised `NameError` for every input over 1000 products.
* **The mean-field cross term used transposed parameters** whenever `d1 > d2`:
  an `(n₁, n₂)` distance matrix against `(n₂, n₁)` parameters. Flattened lengths
  match, so numpy broadcast it silently. **98.9 % wrong** on dyes whose elements
  differ in radius.
* **The dye-library cache exposed a test writing through a frozen dataclass**
  with `object.__setattr__`, which had been harmless only because every call
  re-parsed.
* **`compute_rotamer_score` validated its `potential` argument after an early
  return**, so a misspelling was accepted whenever the atom selection was empty.
* **No `#pragma omp` in this module has ever run.** IMP's CMake leaves
  `OpenMP_CXX_FLAGS` empty, so every figure above is single-threaded against
  single-threaded numpy — and a lower bound.
  (`okf/validation/openmp_is_not_enabled.md`.)

**And `pyext/` itself**: 21 entries to 6. Fourteen single-line `.i` files folded
into `swig.i-in`, which is IMP's own convention (`core` `%include`s 101 headers
directly); the three with real content renamed `IMP_bff.*.i` so they are
distinguishable from the vendored `numpy.i`. Two name collisions removed —
`fret/kappa2.py`, a 17-line re-export shim, and `distribution.py` beside
`distributions.py` in one directory meaning unrelated things.

## 2026-08-18 (PRD-113 stages 6-7: the output contract, io/, and the public surface)

* **`observables/`** — `LifetimeSpectrum`, `(amplitude, rate)` pairs. The
  contract is stated once and *tested*: no convolution, IRF, pileup, counting
  noise or binning anywhere in the package. C++ for the two real loops —
  evaluating the unconvolved decay, and coarse-graining 30 000 species to 128
  while preserving `sum(a)` and `sum(a·k)` exactly (population and initial
  slope) to 1e-12, max decay deviation 5.3e-5. The claimed error scaling is
  measured, not asserted in prose: halving the bin width quarters the error.
* The `exact` flag is a **checked** claim. A constant quenching rate — nothing
  to average over — must reproduce the spectrum's lifetime, and an alternating
  rate must visibly not. Without the second test the first proves nothing.
* **`io/`** — `fret/io.py` was 682 lines of three unrelated formats. Split on
  AST boundaries into `fps_schema`/`fps`/`fps_legacy`/`structure`. A first pass
  extracting only *function* spans silently dropped `AV_SIMULATION_TYPES`, a
  module-level constant in the gap between two functions; caught by diffing the
  set of top-level names before against after, which is now how the split is
  verified. A second check for unresolved names found the one real coupling —
  `read_fps_json` dispatches into the legacy reader by file extension.
* **`api.py` inverted** — `BY_DOMAIN` authored, `EXPORTS` derived, and the
  naming-families regex deleted. It matched export *names* against a token list
  that grew with every feature and could not tell a misfiled name from a
  well-filed one. Replaced by checks that do not grow: every export's module
  must lie inside its declared domain, no name in two domains, every domain
  importable.
* **The import-discipline test paid for itself on its first run**: a
  cross-domain relative import written to look local, an invalid escape
  sequence nothing surfaced, and a test-module **basename collision**
  (`test/io/test_io.py` vs `test/cgdye/test_io.py`) — pytest imports by
  basename, so one file was reported as a single ERROR line at the bottom of a
  780-test run and its assertions simply did not run. Basename uniqueness is
  now a test, because the failure mode is silent.
* **`dynamics/`** — the Brownian walk, the Smoluchowski solver and the
  excited-state Monte Carlo out of `quenching/`. They were filed under the first
  physics that used them, not under what they are. `pyext/src/` now matches
  PRD-113's target layout. The `api.py` domain check failed the instant they
  moved (eight exports still filed under `quenching`), which is the drift the
  old regex could not see.
* The runtime domain graph is acyclic and reads as a layering:
  `photophysics ← representation ← av ← restraints`, `dynamics` a leaf,
  `quenching → {dynamics, fret}`, `fret → {io, photophysics, representation,
  restraints}`, `cgdye → {dye, fret, io, representation}`.
* Suite **799 passing**, one pre-existing `IMP.em` MRC failure.

## 2026-08-18 (PRD-113 stages 3-5, tranches 7-10: **numba reaches zero**)

* `AVDistance.h`, `DistanceCalibration.h`, `OrientationFactor.h`,
  `BrownianWalk.h`, `PhotonSimulation.h` — the last **23** kernels.
  `_jit.py` deleted; `test/test_no_numba.py` parses every source file so it
  stays deleted. Suite **741 passed**, and it now runs in **134 s** against 614 s.
* **The gate had to be split three ways**, because "identical output to numba"
  is achievable for some kernels and impossible for others:
  * integer/index arithmetic — **bit-exact** (`split_contact_volume`, voxel for
    voxel over 8 random grids).
  * float reductions — **≤ 7.1e-15**, and the residual is *FMA contraction*:
    clang fuses `dg * i + r0`, numba's LLVM does not. Chasing it would be
    chasing codegen, so the tolerance says so.
  * RNG-carrying — **closed forms and distributions only**. No C++ generator
    reproduces numba's stream. Two `test_frozen_reference` tests asserted that
    stream (`n_accepted == 94567`; 6297 photons to 1e-9) and were replaced by
    the distributions they came from, with both spreads recorded in the
    docstrings.
* **Gating on closed forms found two defects that gating on the old numbers
  had preserved for years.**
  * **`kappasq_all` returned half the orientation factor.** `np.random.random(3)`
    normalised fills the cube's positive octant, not the sphere: `⟨κ²⟩ = 0.333`
    where the rigid isotropic limit is `2/3`. C++ gives 0.6663. **Fixed**;
    no callers, so nothing downstream carried it.
  * **The particle walk diffused at 3D.** Per-component step variance was
    `6 D dt` (the total 3-D MSD used as one component's width) where the
    convention is `2 D dt`. `GridDiffusionSolver` gives `⟨x²⟩ = 2Dt` exactly.
    **Fixed**, after tracing it to QuEst's 2019 docstring — a Berkeley teaching
    page's step *magnitude* transcribed as a per-component σ, unchanged through
    Cython → numba → C++ — and after correcting my own **10× unit error**: I
    wrote that `D = 40 Å²/ns` was implausibly high because "free Alexa488 is
    around 4 Å²/ns", and used that to argue the default had absorbed the factor
    of three. `1 Å²/ns = 10 µm²/s`, so 400 µm²/s **is** 40 Å²/ns — the default
    is the free-solution value entered correctly, and there is no calibration
    of `D` anywhere in the stack to have absorbed anything. Defaults unchanged.
    The fix does expose a real inconsistency the wrong width was masking: the
    particle model defaults `D = 40` (free) and the field model
    `free_diffusion = 8` (tethered), 5× apart under one name.
    `okf/validation/particle_vs_field_diffusion.md`.
  * **`gaussian_rmp_to_rda_mean` settled**: σ is the per-component width of the
    *separation vector*, so the correction is `σ²/Rmp` — the two transverse
    components give `⟨|ε⊥|²⟩ = 2σ²`, halved by the expansion. The canonical
    version had `σ²/(2Rmp)`, i.e. half; one function now, and `rmp ≤ 0` returns
    0 rather than the 3.6e11 Å the clamped denominator gave. Anything that
    fitted `sigma_rda` through the old canonical version has a σ that is √2 too
    large.
* **`import IMP.bff.av` as a first import had been raising ImportError** at
  every prior commit — `representation/__init__` → `distribution` → `av` →
  `representation`. Nothing in 676 tests caught it because every test imports
  `IMP.bff` first. Fixed by deferring the `av` import in `distribution.py`,
  which sits above the builder; all eight subpackages now import standalone.
* **De-duplication the port forced**: one `distance_sample_statistics` behind
  three reductions that disagreed at the limits; one Horner evaluator behind two
  `polynomial_transfer`s that read coefficients in opposite orders (the
  ascending caller now reverses and delegates); one `brownian_walk_in_volume`
  behind two near-identical walks, with mobility as a *field* rather than a
  scalar-and-mask pair. `chi2_score` is now literally one object.
* **Still an owner decision**: `gaussian_rmp_to_rda_mean` is `Rmp + s²/Rmp` in
  one module and `Rmp + s²/(2Rmp)` in the other, 45.8 against 45.4. A test
  asserts they still differ, so it cannot be closed by accident.
* **Two improvements the C++ made possible**: `photon_trace` is parallel *and*
  reproducible (each photon's generator seeded from `(seed, index)`, so the
  result does not depend on the scheduler — numba had to choose one or the
  other), and `quenched_decay` takes a seed at all, which the numba never did.

## 2026-08-18 (PRD-113 stage 5, tranches 5-6: maps, FRET traces, the solver — and a 40× regression, found and fixed)

* `QuenchingMap.h`, `FRETRateTrace.h`, `DiffusionSolver.h` — 7 more kernels.
* **Gates**: `slow_near_atoms`, `fret_map`, both grid stampers and the
  Smoluchowski step are **bit-exact**. `quenching_map` differs by 2.6× machine
  epsilon (`np.exp` vs `std::exp`), the trace kernels by 3.2×, the Itô step by
  4e-19. Each was isolated relatively rather than judged on an absolute number:
  3.6e-12 looked alarming until measured against values of order 141.
* **The port made the quenching suite 40× slower — 14 s to 558 s — and that is
  the interesting part.** numba's kernels were `parallel=True` over `prange`;
  the C++ was serial, allocated a grid per call, and crossed the SWIG boundary
  **once per step**. A solve is 1 000–10 000 steps.
* **The fix was to move the loop, not to tune the kernel.** `diffusion_propagate`
  runs the whole ping-pong loop in C++ with OpenMP over x-slabs: 1.33 ms/step →
  **0.023 ms**, 58×. Wiring `run()` alone changed nothing measurable — the hot
  path was `equilibrium()`, which takes tens of thousands of steps. With both
  inside: **488 s → 18.3 s**, now slightly better than the numba it replaced.
* The lesson generalises to the rest of the port: **a kernel called in a tight
  Python loop must take the loop with it**, or the binding cost dominates
  whatever the kernel gains.
* numba: **25 → 23**. Suite **676 passed**.

## 2026-08-18 (PRD-113 stage 5, tranche 3: ASA and the PET rate to C++)

* `SolventAccessibleSurface.h` / `.cpp` — `sphere_points`,
  `solvent_accessible_surface_area` (Shrake–Rupley), `quenching_rate_per_frame`.
* **Gate**: `quenching_rate_per_frame` bit-exact. The other two differ by up to
  2.4e-6 — **and that is the Python being wrong, not the port.** It built the
  golden spiral in `float32`; the C++ agrees with an independent `float64`
  construction to **1e-15**, where the old one was off by 2.4e-6. So the port is
  the more accurate of the two, and the difference is not a tolerance to relax
  later.
* **Another uncovered function, found the same way as the last one.** `pet.py`
  was missing `import IMP` after the delegation, and **139 quenching tests
  passed anyway** — nothing calls `quenching_rate_per_frame`. That is the second
  time in two tranches that a green suite has vouched for a function it never
  runs (`worm_like_chain_linker` was the first).
  `test/quenching/test_ported_kernels.py` now covers all three, including the
  neighbour-cutoff bug the ASA kernel was ported with: an occluder at 5 Å must
  be seen, which the old 2.45 Å cutoff would have missed.
* **The SWIG name collision recurs by construction.** Any C++ function SWIG
  binds into `IMP.bff` takes that flat name, so `api.py` must not also claim it
  — `sphere_points` hit this exactly as `normal_distribution` did. Worth
  expecting for every remaining tranche rather than rediscovering.
* numba: **35 → 33**. Suite **676 passed**.

## 2026-08-18 (PRD-113 stage 5, tranche 2: polymer chains to C++ — and a live breakage found)

* `PolymerChain.h` / `.cpp` — `gaussian_chain_ree`, `gaussian_chain`,
  `worm_like_chain` (the Becker–Rosa–Everaers multi-piece solution, both sides
  of the κ = 0.125 branch), `worm_like_chain_linker`.
* **Gate: agreement to ≤ 3 × 10⁻¹⁷** on every case. One line reported a relative
  failure at κ = 0.05 — on entries of order 1e-18, where a relative tolerance is
  the wrong instrument; the absolute difference is 6e-18.
* **The port found a function that had been raising for a whole commit.**
  `worm_like_chain_linker` was numba-jitted and called `normal_distribution`,
  which tranche 1 had just turned into a C++ delegation numba cannot type — so
  every call raised `TypingError`. **The suite did not notice: 656 tests passed
  with it broken**, because that function had no test at all.
* So its equality gate had no working reference. It was checked against an
  **independent numpy convolution** instead — agreement to 7e-18 — and
  `test/test_polymer_chain.py` now covers it, along with both branches of the κ
  expression and the contour-length cutoff.
* The lesson is about the gate, not the port: an equality gate only proves
  something when the thing being compared against still runs. Where it does not,
  say so and find another reference.
* **IMP concatenates a module's sources**, so a helper in an anonymous namespace
  in two `.cpp` files is a redefinition rather than two private copies. Shared
  as `include/internal/Normalize.h`.
* numba: **40 → 35**. Suite **665 passed**.

## 2026-08-18 (PRD-113 stage 5, tranche 1: distributions to C++)

* **numba must go entirely** (owner). 44 jitted functions across 13 files; this
  is the first four, and it establishes the pattern for the rest: C++ header +
  source, a SWIG `.i`, `Files.cmake`, a Python wrapper keeping the old
  signature, and **an equality gate against the numba before anything is
  deleted**.
* `Distributions.h` / `Distributions.cpp` — `poisson_0toN`,
  `normal_distribution`, `generalized_normal_distribution`,
  `distance_between_gaussian`.
* **Gate: 7 of 9 cases bit-for-bit identical**, the other two within
  **7 × 10⁻¹⁸** — the last bit, from a different summation order in the
  normalisation.
* **Two things the port caught by forcing the formulas to be read.**
  `generalized_normal_distribution` is *not* the exponential-power family its
  name suggests: it skews by transforming the axis,
  `z = −log(1 − κ(x−μ)/σ)/κ`, and evaluates the **standard** normal at `z` —
  `loc` and `scale` are already folded in. My first header documented the wrong
  distribution and the wrong default (`shape=2.0`, when the Python default is
  `0.0`, i.e. no skew). Both corrected against the body rather than the name.
* **A name collision, resolved deliberately.** SWIG binds C++ functions into
  `IMP.bff` directly, so `IMP.bff.normal_distribution` was already taken and
  `api.py`'s lazy hook could never fire — the test caught it. The C++ functions
  **are** the public surface; `IMP.bff.distributions` keeps numpy-returning
  wrappers for use inside the package, and `api.py` records why those names are
  deliberately absent.
* numba: **44 → 40**. Suite **656 passed**.

## 2026-08-18 (PRD-113 stage 4c: interaction terms)

* **The photophysics abstraction exists now**, shaped like a force field's
  terms: a functional form, an **arity**, and parameters looked up by type.
  1-body `RadiativeTerm`, 2-body `PETTerm` (dye × quencher atoms) and `FRETTerm`
  (dye × dye), with the N-body case shaped for but not implemented.
* **Rates add, and that is a first-class property.** `total_rate()` sums the
  channels because parallel deactivation channels add — which is exactly what
  `GridDiffusionSolver` already relies on when it sums a quenching map and a
  FRET map. Burying that inside each observable is how `fret_rate_trace` and
  `fret_rate_map` became two incompatible calls.
* **They are a consolidation, not a seventh layer, and the tests say so.**
  `PETTerm` reproduces `quenching_rate_map`'s law at the voxel centres to
  1e-9 (minus the `1/tau0` floor, which is the radiative term's job — that
  separation is the point of having terms). `FRETTerm` reproduces
  `fret_rate_trace` exactly. The kernels can move underneath the terms later
  without changing an answer.
* **`R0` is derived inside the term** from the two dyes, the medium's refractive
  index and κ² — not passed in. And the term reports
  `used_isotropic_kappa2`, so a representation that cannot resolve orientations
  (an accessible volume) is *told* it fell back to 2/3 rather than silently
  averaged.
* **Terms are representation-agnostic**: they consume
  `IMP.bff.representation.States`, so one implementation serves an AV, a rotamer
  library, a coarse-grained model and an MD trajectory. `needs_orientations`
  says which terms care.
* Suite **656 passed**.

## 2026-08-18 (PRD-113 stage 4b: the distance layers, and two that disagree)

* **`fret/distance.py` → `representation/distance.py`.** A distance between two
  labels is a property of what represents them, not of the FRET engine. It was
  the most complete of the implementations — pair statistics, the RDA histogram,
  the transfer polynomial, orientation-resolved pair geometry — so it is the
  canonical one. `fret/distance.py` is a re-export while the rest of `fret/`
  migrates.
* **I was wrong that there were six duplicates.** That claim came from matching
  *function names* (`comm -12` over `def` lines) without ever comparing
  signatures. Deleting on that basis broke 8 tests, which is how it was caught.
  The truth, checked properly:

  | | verdict |
  |---|---|
  | `fret_efficiency`, `distance_from_fret_efficiency` | identical — removed, re-exported |
  | `chi2_score` | identical but for a parameter name — kept, it is the only symbol `restraints/` imports |
  | `av_pair_statistics` | **different jobs**: array reduction here, sample-and-reduce there |
  | `gaussian_rmp_to_rda_mean`, `polynomial_transfer` | **different answers** |

* **Two same-named, same-signature functions give different numbers.**
  `gaussian_rmp_to_rda_mean(45, 6)` is `Rmp + σ²/Rmp` = **45.8** here and
  `Rmp + σ²/(2 Rmp)` = **45.4** there — a factor of two in the correction term.
  `polynomial_transfer(45, [0, 1, 0.02])` is **85.5** here (ascending
  coefficients) and **45.02** there (descending, `np.polyfit` order).
* **Which is live matters and is now known.** `fret/engine.py` calls the
  `representation.distance` versions, and its `fit_transfer_polynomial` produces
  `np.polyfit` coefficients — which only the descending evaluator reads
  correctly, so that pairing is consistent. The `distance_metrics` versions have
  **no consumers**. This is a trap that had not yet sprung, not a live defect.
* **Which `gaussian_rmp_to_rda_mean` is correct is a physics question** — whether
  `sigma` is the per-component width of an isotropic 3-D cloud or the width of
  the distance distribution — and is left for the owner rather than guessed.
  Both are documented side by side at the top of `distance_metrics.py`.
* Suite **641 passed**.

## 2026-08-18 (PRD-113 stage 4a: `photophysics/`, κ² unified, `spectroscopy/` gone)

* **933 lines of anisotropy modelling were unreachable.**
  `spectroscopy/kappa2.py` — κ² distributions, wobbling in a cone, order
  parameters, the κ²→distance-ratio conversion — had **no consumers at all**: it
  was never exported through `api.py`, and it imported `numba` *directly* rather
  than through `IMP.bff._jit`, so an installation without numba could not even
  load it. Both fixed; it is `photophysics/orientation.py` and is now exported.
* **κ² has one home.** The geometry (`kappa2_from_dipoles`, 44 lines, real
  consumers) and the distributions (933 lines, none) were in unrelated packages
  named the same thing. `photophysics/kappa2.py` and
  `photophysics/orientation.py` now sit side by side. `fret/kappa2.py` is a
  re-export while the rest of `fret/` migrates.
* **`spectroscopy/` is gone** — its two halves were an experiment-side decay
  wrapper (deleted in stage 0) and a forward model (moved here).
* Three docstrings became raw strings: LaTeX in them raised a `SyntaxWarning` on
  every import.
* Suite **641 passed**.

## 2026-08-18 (the sampler: the table transfers, the short test did not)

* PRD-111's posterior harness used the **stretch** move on a posterior whose
  condition number is ~10⁵. ChiSurf benchmarks exactly that case
  (`docs/development/benchmarks.md`) and reports 0.0036 → 0.0109 ESS per
  evaluation from stretch to an adaptive-covariance slice at κ=100, noting the
  gap grows with correlation. The machinery had been surveyed and then not used.
* **An 80-step comparison said the opposite** — slice 6× *worse* per evaluation —
  because a slice step costs 5.7 log-probability evaluations and neither chain
  had mixed (R̂ 2.3–2.5), so the ESS being compared was meaningless.
* **600 steps settles it**, one site, 12 walkers, 3.0 Å:

  | | evaluations | wall | ESS/eval | ESS/s | worst R̂ |
  |---|---|---|---|---|---|
  | stretch | 5 986 | 103 s | 0.00113 | 0.066 | **3.040** |
  | slice (adaptive covariance) | 35 428 | 753 s | **0.00223** | **0.105** | **1.137** |

  2× per evaluation and 1.6× per second — but the number that decides it is R̂:
  **stretch does not mix at all**, and got *worse* from 80 steps (2.333) to 600
  (3.040), so its ESS never meant anything. The slice sampler is the default in
  the harness now.
* The forward model still dominates: one six-site evaluation is ~0.36 s. A
  cheaper solve is a bigger lever than the proposal — the explicit scheme's
  `dt ≤ dg²/(6D)` sets the step count, and an implicit or exponential treatment
  of the *diffusion* term would lift it the way the exponential rate integration
  already lifted the reaction constraint.

## 2026-08-18 (PRD-113: stage 3e was mis-scoped — six distance layers, and they agree)

* **The C++ already has the distance API.** `AV.h` exports `av_distance`,
  `av_distance_quadrature`, `av_random_points`, `av_random_distances` and
  `av_distance_distribution`, operating on `AV` decorators. Stage 3e was written
  as "port `av/_kernels.py` to C++"; the six numba kernels are a **duplicate**
  of that, taking point arrays instead of decorators. It is de-duplication, not
  a port, and it belongs with the other distance layers in stage 4.
* **There are six implementations of dye-pair distances**: `av/_kernels.py`
  (numba), `av/basic.py` (wrapping it), `fret/distance.py`,
  `distance_metrics.py`, `representation/distribution.py` (once on each concrete
  class), and the C++ `av_distance` family.
* **They agree.** T4L A132×A65, 2180×3087 points, 2×10⁵ samples: `<R_DA>` 51.945
  against 51.940 Å, `<R_DA>_E` 51.696 against 51.692 Å, `Rmp` identical to all
  printed digits. The spread is Monte-Carlo sampling noise. **Unlike the
  transposed density, no defect is hiding in this one** — the merge is safe and
  the C++ is the implementation to keep.
* Remaining work for it: C++ overloads taking point arrays, plus SWIG typemaps
  for a *second* 2-D input array — `BFF.types.i` defines one
  (`double *input, int n_input1, int n_input2`) and a two-cloud signature needs
  two. Recorded rather than started, so the change lands in one piece.

## 2026-08-18 (PRD-113 stage 3d: the label/representation collision resolved)

* **`LabelDistribution*` moved to `representation/`.** `label/` means the
  *system* — which dye is attached where — and a label *distribution* is a
  representation of where it can be. One word, two questions; the collision was
  recorded in stage 2 and is now removed. `label/` exports only `Label`,
  `Quencher`, `PETParameters` and the reference tables.
* **They are a third implementation of the `States` surface** (`points`,
  `mean_position`, `n_points`) and a **fourth** copy of the distance layer —
  `dRmp`/`dRDA`/`dRDAE`/`pRDA` duplicated across both concrete classes, and
  again in `BasicAV`, and again in `fret/distance.py`. A `states` view bridges
  them for now; folding four distance layers into one is stage 4, and is a
  behaviour change that does not belong in a move.
* `representation/` now holds every way of describing where a dye is:
  `AccessibleVolume`, `LabelDistributionAV`, `DyeDistributionNormal`, and the
  rotamer library as a sibling — all over `States`.
* Gate 12/12. Suite **636 passed**.

## 2026-08-18 (PRD-113 stage 3c: States, and rotamers stop being an AV)

* **`States` is the abstraction that was missing.** An AV grid point, a rotamer,
  a coarse-grained conformer and an MD frame are the same kind of thing — a state
  the dye can occupy, with a weight, a position and possibly an orientation.
  Everything downstream consumes states, so distances and κ² are written once.
* **`RotamerEnsemble` no longer inherits `AccessibleVolume`.** That inheritance
  was *why* distance code worked for rotamers — by inheritance, not by design —
  and it forced a rotamer library to carry a grid it does not have:
  `density=zeros((0,0,0))`, `grid_step=0.0`, `grid_shape=(0,0,0)`. **No consumer
  ever read them**: every one used `points`, `mean_position`, `n_points` or
  `has_volume`, which is exactly the `States` surface. Both are now siblings.
* **One `AccessibleVolume`.** There were two, in `av/compute.py` and
  `fret/av.py`, with identical seven-field definitions in identical order —
  which is how they came to disagree about the axis order of `density` without
  anything noticing. Now `representation/types.py`, `AccessibleVolume(States)`.
* **`mu` and `orientations` are one array.** The rotamer code keeps its name;
  representation-agnostic code asks for `orientations`. Being the same object
  means no caller can set one and read a stale other — and it is what lets a
  consumer see that a rotamer library resolves dipoles where an AV does not.
* **A regex over-matched and took `_av_imp_bff` with the dataclass.** Caught by
  the suite; redone by locating the class and the next top-level definition
  rather than by pattern. Same lesson as the `api.py` rename map: structured
  code is not a place for a regex.
* Gate still 12/12 byte-identical. Suite **633 passed**.

## 2026-08-18 (PRD-113 stage 3b: one path-map core, two front doors)

* **`representation/pathmap.py` is now the one place `IMP.bff.AV` is driven.**
  Both builders set up the same decorator, call the same `resample()` and read
  the same `PathMap`; both had to learn the same lessons, and each had learned
  only some. The array door knew about the x-fastest voxel order; the structure
  door did not, and shipped a mirrored density (stage 3a). Sharing the core is
  what stops that happening a third time.
* **Gate: 12/12 byte-identical** — 3 sites × 2 resolutions × 2 doors, comparing
  `points`, `density`, `grid_origin`, `attachment_point`, `grid_step` and
  `grid_shape` by SHA-256. Held after each door was repointed and again after
  the dead code was removed. `av/compute.py` 405 → 350 lines.
* **What is deliberately *not* shared, and why.** The two doors disagree on
  conventions: the array door binarises the density to 0/1 float64 and keeps
  IMP's point weights; the structure door keeps raw float32 values and forces
  the weights to one. Those are **behaviour**, and this stage forbids behaviour
  changes — so `PathMapReading` returns the raw readings and each door applies
  its own conventions, with the disagreement recorded in the dataclass docstring
  rather than silently resolved. Choosing between them is a later, deliberate
  decision.
* The core keeps both hard-won details, one from each door: the AV is decorated
  onto its **own** particle with the source passed separately (otherwise the
  resampled map sits at the coordinate origin), and the build is serialised
  while `resample()` is not (which is ~96 % of the wall clock, so the split is
  what makes a threaded caller worth having).
* Suite **631 passed**, same single pre-existing `IMP.em` failure.

## 2026-08-18 (PRD-113 stage 3a: the AV density was transposed)

* **`IMP.bff.compute_av` returned `density` mirrored relative to its own
  `points`.** IMP numbers voxels with *x* fastest (`i = x + nx*y + nx*ny*z`), so
  `fret/av.py`'s C-order reshape into `(nx, ny, nz)` — which makes the *last*
  axis fastest — transposed the volume. Measured on T4L A132 at 1.5 Å: of 2180
  cloud points, **1548 (71 %)** landed on a voxel the density called occupied,
  against **100 %** after transposing. Confirmed at 1.5, 2.0 and 2.5 Å.
* **It survived because a mirrored volume looks right.** Same voxel count, same
  bounding box, same total volume — only the per-voxel relationship to the
  protein is wrong. Nothing but a voxel-by-voxel comparison catches it, and
  nothing was doing one.
* **The fix was already written, in the other builder.** `IMP.bff.av.compute`
  (the array path, from imp-tricks) does
  `reshape((nz, ny, nx), order="C").transpose(2, 1, 0)` and documents the hazard
  in a docstring — *"a transpose still overlaps the truth by ~83 % of its voxels,
  so the point cloud has to be compared against the grid"*. It also reads the
  point cloud from IMP rather than deriving it from the grid, **deliberately, so
  that two independent readings can disagree**. The structure path had neither
  the fix nor the test. This is the case for the merge in one example.
* **Downstream effect measured, not assumed**: on T4L A132 at 2.0 Å the mean
  donor lifetime is **3.5626 ns** correct against **3.3974 ns** as computed —
  4.9 % — with `k_max` 2.54 against 2.81.
* **PRD-110 and PRD-111 read `av.density` and built their maps from it**, so
  every number in them was computed on an accessible region mirrored relative to
  the protein. Their *conclusions* concern the model's sensitivity structure —
  θ unidentifiable from one decay, the multi-site rotation, `contact_distance`
  uninformative, Smoluchowski against Itô — and a mirrored AV is still a
  plausible AV-shaped region, so those are unlikely to hinge on the reflection.
  **That is an expectation, not a measurement**: re-running them is owed.
* **`test/fret/test_av_grid_registration.py`** pins it at three resolutions, plus
  a guard on the guard — asserting that a *transposed* density actually fails the
  test, so it cannot quietly become vacuous.
* No quenching test broke, because they pin against analytic results rather than
  recorded numbers. Suite **631 passed**.

## 2026-08-18 (PRD-113 stage 2: `label/`, PET as a pair property, CIF as the data format)

* **`Label` and `Quencher` replace two untyped dicts.** A labelling site was a
  `source_info` dict threaded through the AV builder, the quenching model and
  every benchmark — and it mixed *where the dye is attached* with *how its
  accessible volume is computed*. `Label` carries only the first half;
  `linker_length`, `radius1..3`, `allowed_sphere_radius` and
  `simulation_grid_resolution` are AV *representation* parameters, and a rotamer
  library has none of them. Pinned by a test.
* **PET is a pair property, not a quencher property** (owner, 2026-08-18). `kQ`
  depends on the redox potentials of *both* partners — a rhodamine, an oxazine
  and a cyanine see the same tryptophan differently — so a `Quencher` carries
  identity only (`comp_id`, the redox-active `atom_ids`, optional location) and
  `PETParameters` is keyed by `(dye, comp_id)`. My first cut put `rate_constant`
  on the `Quencher`, which would have hardened exactly that error.
* **The bundled table says so itself.** Its docstring reads *"reference PET
  parameters for a xanthene dye (Alexa488-like)"* — and the package applies it to
  every dye. `reference_pet_parameters(dye)` keeps that transfer but makes it
  visible: each entry records `measured_for`, so `is_transferred` reports when an
  assumption was made rather than a measurement used.
* **CIF is the data format** (owner). The dye library was 41 CSV files — one
  extinction/QY table and one curve file per dye — with no category names, no
  units and no provenance. Now one `dye_library.cif` (38 dyes, 26 169 spectrum
  points) read through `ihm.format`, the CIF layer the package already uses. The
  CSVs are gone.
* **Two CIF traps, both caught by checking rather than by assuming.**
  `ihm.format.CifReader` reads its keywords from the **handler's `__call__`
  signature**, so a `**kwargs` handler is silently handed nothing — the first
  reader returned zero dyes. And `CifWriter._repr` formats floats with `"%.3f"`,
  which rounded the spectra to three decimals and moved R0 by 0.0065 %
  (5.687415 vs 5.687781 nm); its docstring names the escape hatch — hand it a
  string. After that the CIF route reproduces the CSV R0 **bit for bit**
  (5.687781310385014).
* **`_bff_dye` and `_bff_dye_spectrum` are bff-native, and that is recorded.**
  No dictionary in the stack defines an item for a quantum yield, an extinction
  coefficient or a spectrum — checked across all ten `.dic` files in
  `../mmfdb/src/mmfdb/data`; the only matches are
  `_em_detector.detective_quantum_efficiency` and the NMR spectral categories.
  They should be proposed for `mmfdb_flr_ext.dic`. Identifiers still follow
  flrCIF (`chromophore_name` is `_flr_probe_list.chromophore_name`).
* **`_flr_fret_forster_radius.index_of_refraction` and `.kappa_squared` exist** —
  added by `mmfdb_flr_ext.dic`, not upstream IHM-FLR — so stage 1's new
  `refractive_index` parameter has a dictionary item after all, and a stored R0
  can be reproducible rather than a bare number.
* **A name collision found**: `label/` already meant label *distribution* (a
  representation concept) and now also means label = dye at a site (a system
  concept). `LabelDistribution*` is re-exported meanwhile and moves to
  `representation/` in stage 3.
* Suite **626 passed**, same single pre-existing `IMP.em` failure.

## 2026-08-18 (PRD-113 stage 1: `dye/` — the species)

* **A dye had two unrelated descriptions.** The *spectral* one in
  `fret/forster.py` (names, extinction coefficients, quantum yields, curves from
  the bundled tables); the *molecular* one in `cgdye/io/template_cif.py` (atoms,
  transition-dipole atoms, formal charges). Nothing connected them, so "what is
  AlexaFluor 488" had two answers. `IMP.bff.dye.Dye` is now the one answer, and
  it is deliberately **model-independent**: `linker_length` and
  `allowed_sphere_radius` are AV *representation* parameters, and a rotamer
  library has neither.
* **R0 is derived, not supplied.** `forster_radius(donor, acceptor, kappa2, n)`
  takes two `Dye`s and the medium. The old route hard-coded the refractive index
  as the literal `1.4**4` inside the calculation, so nothing could ask what R0
  would be in a different solvent — it is now a parameter (AlexaFluor488/594:
  5.6878 nm at n = 1.4, **5.8856 nm at n = 1.33**). The name-based
  `forster_radius_from_spectra` is kept for the rotamer code and agrees exactly.
* **Names aligned to flrCIF**, checked against `python-ihm`'s FLR model rather
  than guessed. `Dye.name` → `_flr_probe_list.chromophore_name`,
  `chromophore_center_atom` → `_flr_probe_descriptor.chromophore_center_atom`,
  `lifetime` → `_flr_reference_measurement_lifetime.lifetime`, plus
  `reactive_probe_name`/`probe_origin`/`probe_link_type`. flrCIF's word for a dye
  is **probe** (the reagent) and for the fluorescent moiety **chromophore**.
  Recorded as `FLRCIF_ITEMS`, same convention as `fret/fps_schema.py`.
* **flrCIF has no item for quantum yield, extinction coefficient or a spectrum.**
  The nearest, `_flr_fret_calibration_parameters.phi_acceptor`, is an analysis
  calibration value, not a species property. Those three are bff-native and
  marked `None` deliberately rather than by omission.
* **`_flr_reference_measurement_lifetime` is `(species_fraction, lifetime,
  species_name)`** — the lifetime-spectrum output contract for stage 6 already
  exists in the dictionary, named. Worth adopting verbatim there.
* **Two defects found and fixed on the way.** `api.py` had five names appearing
  twice; one (`read_component_template_cif`) was a genuine duplicate, and the
  other four were **not duplicates at all** — they live in a second dict, the
  rename map from export name to the module's real attribute. A regex dedup
  deleted that map and broke four exports; caught by the suite, restored. The
  lesson is the obvious one about regex edits on structured data.
* **A plan assumption corrected**: PRD-113 stage 1 said to move
  `cgdye/io/template_cif.py` into `dye/`. Having read it, it is a *generic
  component* template reader used by `topology/builder`, `analysis/density` and
  scripts — only `read_dye_template_cif` is dye-specific. It stays put and moves
  to `io/` at stage 7; `dye/library.py` calls into it meanwhile.
* Suite **605 passed** (from 599), same single pre-existing `IMP.em` failure.

## 2026-08-18 (PRD-113 stage 0: the instrument layer goes)

* **`IMP.bff` is a forward-model engine** — what `tttrlib` is to ChiSurf for
  photons, `IMP.bff` is for structure. It emits **experiment-neutral** quantities
  (lifetime spectra, rate constants, κ² distributions, distances) and leaves
  convolution, IRF, pileup, counting statistics and fitting against raw data
  outside. Scoring against *structural* data (distances, as `AVNetworkRestraint`
  does) stays. Recorded as [PRD-113](prds/prd-113.md); PRD-112 is absorbed into
  it, its LabelLib stage done and standing.
* **Deleted the TCSPC instrument layer**: `src/Decay*.cpp` (10) +
  `src/PhotonStatistics.cpp`, `include/Decay*.h` (11) +
  `include/internal/PhotonStatistics.h`, `pyext/Decay*.i` (11) and their
  `%include` lines, `test/test_Decay*.py` (9), and
  `pyext/src/spectroscopy/decay.py`. 44 files. It did convolution, pileup,
  linearisation and decay *fitting* — all experiment-side — and its only consumer
  was that one Python module. **Nothing in chisurf, quest, imp-tricks or tttrlib
  referenced any of it**, and no AV/PathMap source included it, so the cut was
  clean.
* **`spectroscopy/kappa2.py` deliberately survives.** It is the opposite kind of
  thing: a forward model producing orientation-factor distributions,
  wobbling-in-a-cone averages and order parameters from structure. It moves to
  `photophysics/` in stage 4, joining `fret/kappa2.py` so that "kappa squared"
  has one home instead of two. `spectroscopy/__init__.py` re-exported the deleted
  module and is now a docstring saying which half went and why —
  `test_tttrlib_is_optional.py` imports that package, so it had to keep
  importing.
* **Verified**: cmake reconfigured, module rebuilt, `DecayCurve` and friends
  absent from `IMP.bff`, `AV`/`PathMap`/`AVNetworkRestraint` intact. Suite **599
  passed**; the single failure, `test_AccessibleVolume::test_access_av_feature`,
  is the pre-existing unrelated `IMP.em` MRC error — the gate is that it neither
  disappears nor gains company, and it did neither.

## 2026-08-18 (solver stability; the rate term)

* **The explicit-step criterion ignored the rate term, which dominates it**
  (`okf/validation/quenching_solver_stability.md`). The update coefficient is
  `1 − 6D dt/dg² − k dt`; `diffusion_stability_limit` validated only the
  diffusion half. At T4L site 19, 2.5 Å, with a 25 Å²/ns ceiling on `D` setting
  the step: diffusion 0.16, **quenching 2.01**, sum 2.17 — past *divergence* —
  and the decay reached **7 × 10³⁶** while the solver called the step safe.
* **The divergence need not look like one.** Site 124 at a sum of 2.10 returned a
  smooth, finite, monotone, plausible decay that was **2.6 % wrong** at 25 ns.
  That is why it is pinned by a test rather than left to inspection.
* **Two thresholds, and they differ by two**: the coefficient goes negative at a
  sum of 1 and the scheme diverges at 2. `diffusion_stability_limit(d_max, dg,
  k_max)` returns the positivity bound, the stricter one, because a negative
  probability density is not an acceptable answer either.
* **The rate is now integrated exactly** — `exp(-k dt)` as a factor rather than
  `1 - k dt` subtracted — so it contributes no stability constraint at all and a
  strongly quenched site costs the same as a weak one. Pinned: a uniform rate
  reproduces `exp(-k t)` to 1e-10 whatever the step.
* **A published claim retracted.** On finding this I said PRD-111 stages 0–1 were
  contaminated. Wrong twice: I read the positivity bound as the divergence bound,
  and I computed the thresholds with the *new harness's* step, 1.6× larger than
  the one those stages used. They ran at sums of 0.92–0.98, and re-measuring
  reproduces them to ~1 % (eigenvalues 8.223e5/2.992e4/736/129/2.67 against
  8.209e5/2.994e4/735/128/2.66). **The divergence belongs to the new harness**:
  tightening the `D` ceiling to save compute raised the step and pushed the
  contact sites past 2.
* **A real hazard removed**: at 2.0 Å the fits could search `kQ_scale` to 10,
  where the sum reaches 8.6 and the old scheme would have diverged. They
  converged near 1.0 and never went there, but nothing was stopping them.

## 2026-08-18 (the flux form, and the estimator)

* **The flux discretisation was wrong, and it changes every earlier quenching
  finding** (`okf/validation/quenching_flux_form.md`). PRD-110 recorded it as
  *"a convention, not a derivation"*; that was too generous. Equilibrium is
  thermodynamics and mobility is kinetics, so a dye slowed by friction with no
  attraction must still be found uniformly across its accessible volume — while
  the inherited `d[i]p[i] − d[j]p[j]` flux gives `p ∝ 1/D`. The field's own
  canonical treatment agrees: the **Haas-Steinberg** equation writes the
  diffusion operator as `D ∂/∂r [ p(r) ∂/∂r ( N/p(r) ) ]` precisely so its
  stationary state is the *given* `p(r)` for any `D`. QuEst's notebook 04
  integrates it, but at constant `D` and uniform `p(r)` — the one case where both
  discretisations coincide.
* **`flux_form="smoluchowski"` is the default; `"ito"` reproduces the inherited
  behaviour.** Both closed forms verified against the kernel at a 16× mobility
  contrast: peak/min occupancy 1.0000 against 16.0000.
* **`slow_factor` was a disguised attraction.** Reducing `D` near an atom
  concentrates the dye there and the quenchers are exactly where `D` is reduced,
  so `slow_factor^n_contacts` made the attraction exponential. Mean lifetime on
  T4L A132 at `D = 0.5` going 0.985 → 0.90: **−32.2 % under `ito`, +0.2 % under
  the default** — a factor of ~160 in how much the parameter matters.
* **Four earlier findings are artifacts and are corrected in place**: PRD-110's
  closed-form `p ∝ 1/D` (right algebra, wrong operator — the practical half
  stands), its *"`slow_factor` only meaningful within a whisker of 1.0"*, its
  *"decay becomes exactly independent of `D`"*, and PRD-111 stage 0's
  *"`contact_distance` and `slow_factor` are one parameter"*. They traded because
  `slow_factor^n` set the strength of the spurious attraction; under the default
  `slow_factor` recovers to +1.1 % and `contact_distance` is simply
  **uninformative** (±61 %, its own eigenvector at λ = 2.66).
* **Where it counts the model got better**: the joint six-site fit now recovers
  `kQ_scale` to −0.0 % and `rC` to −0.1 %, in 312 forward solves against 516. The
  leading Fisher eigenvalue falls 34×, which is the size of the information that
  was really an attraction. The multi-site argument survives (rotation gain 396×).
* **PRD-111 stage 1 was not re-measured** — the free-dye nuisance, the 9.4 σ
  detectability and the benign geometry error all ran under `ito`.
* **The model now has no stickiness at all rather than the wrong kind.** Real
  dyes stick; the correct form takes a separate equilibrium,
  `flux_ij = D_ij p_eq,ij (p_i/p_eq,i − p_j/p_eq,j)`, reducing to what is
  implemented when `p_eq` is uniform. Future work, and its own identifiability
  question.
* **The estimator was `neyman_lsq`, hand-rolled, and it is the wrong one.**
  tttrlib ships the objective registry (`fit_objectives_json()`) and my
  benchmark reproduced `statistics::neyman` down to the `max(1, ·)` clamp — which
  tttrlib itself documents as *"biased low at small counts"*, against
  `poisson_mle` (2I\*) being *"what a TCSPC decay should normally be fitted
  with"*. These decays run to nothing inside the window, so most of the axis is
  that regime. One mitigation, and it flattered the study rather than the
  reverse: σ came from the noiseless model at the true θ, an oracle weighting no
  real fit has.
* **A point estimate is the wrong output**, given a 10⁵ condition number and a
  demonstrated second minimum — the answer is a posterior. Prior information is
  being discarded too: `PET_QUENCHING_REFERENCE` ships as *"starting values meant
  to be calibrated"*, which is a prior, and a bounded uniform search over
  `kQ_scale ∈ [0.05, 10]` throws it away.
* **All of it already exists, in the right layers.** Statistics in tttrlib
  (`twoIstar`, `statistics::{neyman,poisson,pearson,gauss}`); inference in
  ChiSurf (`FitGroup.run(local_first=True)` staged group fitting, `priors.py`
  MAP, `sample.py`/`ensemble.py` MCMC, `diagnostics.py` R̂/ESS/MCSE);
  forward model here. `benchmark/kq_sensitivity_analysis.py` does all three
  inside `imp.bff` and should not grow further — **PRD-111's inference stages
  move to ChiSurf.** The Fisher analysis stands either way, being a property of
  the model and the noise.
* **χ² was never going to show the problem**: 1.01–1.02 in every run, including
  those where `contact_distance` was 30–70 % wrong or pinned at its bound. A
  scalar cannot distinguish an under-determined fit from a misfitting one; the
  residual *shape* can, and nothing here looked at it.

## 2026-08-18 (PRD-111 stage 2 scoped; the acceptor observable)

* **Stage 2 rewritten around the acceptor decay, on the owner's physics
  (2026-08-18): quenching and FRET are spatially uncorrelated and
  non-homogeneous.** Quenching is set by contact with the residues immediately
  around the dye; FRET by `1/R⁶` to an acceptor elsewhere. Both are modulated by
  the same `r(t)` through unrelated geometry. So the acceptor decay is a second,
  **independent** projection of `p(r, t)` — which is what stage 0 showed the
  parameters need — and **nothing is separable**: the decay is not a quenching
  decay times a FRET decay, and neither rate may be averaged over the volume
  before the other is applied. Both must be evaluated along the *same*
  trajectory. That makes `fret_rate_trace` (acceptor averaged over its AV) the
  fast-acceptor approximation and `fret_rate_pair_trace` the general case; the
  known ~7 % gap between them on 148l E15→E90 is this effect, not a numerical
  detail.
* **The acceptor decay carries a rise**, its rate being the donor's *total*
  depopulation — radiative plus quenching plus transfer — seen on a channel free
  of the donor's own detection artefacts, with the amplitude ratio carrying the
  transfer efficiency.
* **It costs parameters as well as adding data.** The acceptor is a different
  dye at a different site with its own local quenchers, so it does **not** share
  the donor's θ. Whether it is a net gain for identifiability is measurable by
  the same Fisher construction as stage 0, and that is the stage-2 work; the gate
  is that it must reduce the joint condition number *after* its own parameters
  are counted.
* **A quenched sub-population is not the FRET-active sub-population**, so
  stage 1's free-fraction nuisance and the transfer efficiency interact — the
  free fraction cannot be assumed FRET-inactive.
* **No experimental data exists at the moment (owner), so "fit real decays" is
  blocked rather than pending.** Until data exists every result in PRD-111 is a
  statement about the model, not about a dye on a protein, and model-perturbation
  (stage 1) is the only validation available. Recorded as a standing caveat on
  the PRD, not a to-do.
* **`fret_rate_pair_trace`'s mismatch error named the wrong lever.** It asked for
  the same `t_max` and `t_step`, which are already shared; the number of **walks**
  is what has to agree. Unequal lengths are the normal case — the two dyes are
  quenched differently because they sit at different sites — so the message now
  says so, and says why truncating is refused (a trajectory concatenates one walk
  per excitation, so an arbitrary cut splits a walk). Matching the counts is
  QuEst-side work, recorded in `../quest/okf/log.md`.

## 2026-08-18 (PRD-111, quenching identifiability across sites)
* **PRD-111 stage 0 measured; the gate passes** (`okf/prds/prd-111.md`,
  `okf/validation/quenching_multisite.md`,
  `benchmark/kq_sensitivity_analysis.py`). PRD-110 stopped because θ is not
  recoverable from one decay; the successor question is *what has to be
  measured*, and the cheapest candidate is more sites. θ is **global** (one dye,
  one chemistry, one mobility) while the geometry is **per site**, so N decays
  share one θ and the joint Fisher information is the **sum** — which beats its
  terms only if the blind directions rotate with the geometry. Six T4L chain-A
  sites spanning quencher-in-contact (A124, TRP126 at 3.5 Å) to nearly-bare
  (A53): joint condition number **1.8 × 10⁵** against 1.96 × 10⁸ for the best
  single site. The rotation is the whole effect — joint `λ_min` **152.7** against
  **0.0456** for the sum of the per-site minima, a 3 350× gain. A joint fit
  confirms it: `free_diffusion` −52 % → **+2.5 %**, `kQ_scale` −69 % → **+4.3 %**,
  in **97 evaluations against 676**.
* **`contact_distance` and `slow_factor` are one parameter, not two.** Fisher
  predicts ±4.1 % and ±0.2 %; the fit misses by −33.4 % and −2.8 % in exactly
  compensating directions — at **reduced χ² = 1.016**. Not a failed fit but a
  second, exactly equivalent minimum, reproducible across resolutions (4.33 Å /
  0.957 at 1.5 Å, 4.51 Å / 0.960 at 2.0 Å). The two build the *same mobility
  field*: log D correlates at 0.927, equilibrium occupancies overlap **98.7 %**.
  `slow_factor` applies once per contacting atom and `contact_distance` sets the
  count, so `slow_factor^n(contact_distance)` is invariant along a valley. **No
  number of sites separates them.** `JᵀJ` error bars must not be quoted for
  either — a local linearisation cannot see a second minimum, and here it is
  confidently wrong about its *best*-determined parameter. The model needs
  re-parameterising on the quantity the data determines, not more data.
* **Repeated end to end at 2.0 Å**: joint condition 3.04 × 10⁵, rotation gain
  5 623×, largest principal angle 88.4°, same three parameters recovered and the
  same pair degenerate. Not a grid artefact.
* **`compute_av` silently discarded a declared grid resolution, and it had been
  discarding one all along.** It *writes* `simulation_grid_resolution` into
  `source_info` from its `disc_step` argument rather than reading it back, so
  both quenching benchmarks passed the resolution into a field that was
  overwritten: their `--resolution` flag was inert and **every run, PRD-110's
  included, was at the 1.5 Å default while the validation pages recorded 2.5 Å.**
  The numbers were always internally consistent — only the label was wrong, and
  both pages now carry the correction. `compute_av` raises on disagreement
  instead of discarding, the benchmarks assert the spacing they got, and
  `test/fret/test_av_resolution.py` (7 tests) pins it. A resolution is the one AV
  parameter whose being wrong is invisible in the result.
* **The stated hypothesis was wrong, and recorded as such.** Predicted before
  measuring: *diffusion only matters when the dye has to travel*. The ordering is
  the opposite — sensitivity to `free_diffusion` tracks **quenched fraction**,
  and correlates with median dye–quencher distance at **−0.56**. Transport
  becomes rate-limiting when quenching is fast enough to deplete. **To measure
  dye mobility, label where the dye is strongly quenched.**
* **A diagnostic that was pure arithmetic, caught and fixed.** The first
  subspace-angle comparison used rank-3 blind subspaces in a 5-parameter space;
  two subspaces of dimension `r` in `Rⁿ` must intersect once `2r > n`, so it
  reported 0.00° for all 15 site pairs regardless of geometry. Rank is now capped
  at `n/2` and the decisive number is the eigenvalue gain, which no choice of
  rank can flatter.
* **PRD-111 stage 1 measured; that gate passes too.** Two perturbations of the
  model, both realistic. A **per-site free-dye fraction** of 2–7 % (incomplete
  labelling, a stuck rotamer; it competes with `slow_factor` because both slow
  the decay) and a **5 % linker-length error** (data at 21 Å fitted at 20 Å,
  which grows the AV from 1 158 to 1 443 voxels).
* **Unmodelled free dye wipes out the multi-site gain**: `free_diffusion`
  −60.5 %, `kQ_scale` −29.5 % — back to PRD-110's single-decay errors. Six sites
  sharpen a correct model; they do not protect a wrong one. **But it is loudly
  detectable**: reduced χ² 1.683 on 379 dof is a **9.4 σ** rejection, so the
  failure mode is a fit that visibly fails, not a wrong number quietly
  published.
* **Modelling it works, and `free_diffusion` pays.** 11 parameters (5 global +
  6 per-site) restore χ² to 1.021 (0.3 σ) and bring `kQ_scale`, `rC` and
  `slow_factor` back to a few percent, at 2.3× the forward solves — but
  `free_diffusion` settles at −13 %, five times worse than stage 0, because a
  free fraction and a fast diffusion both raise the long-time intensity.
  **Stage 0's ±3 % on `free_diffusion` was a perfect-model number.** The fitted
  fractions themselves are badly recovered (0.040 → 0.009) and must never be
  quoted as labelling efficiencies.
* **The geometry error is benign — the opposite of the expected answer.** The
  5 % linker error was included expecting it to be the damaging one, since a
  joint fit converts per-site model error into global bias. With the nuisance
  fitted, `kQ_scale` returns to **−1.6 %** and `rC` to −2.1 % (both *better* than
  under the correct geometry) and χ² rises only to 1.110 — 1.5 σ, not a
  rejection. **Worry about the sample, not the structure.** Untested and
  different in kind: a *systematic* geometry error shared by every site, which a
  per-site nuisance cannot absorb.
* **All of stage 0 is synthetic data from the model being fitted**, so it
  measures the model's internal geometry under a perfect model and nothing about
  the world. Stage 1 (a per-site nuisance parameter, then real decays) is where
  that assumption is tested, and it is the assumption most likely to be false.

## 2026-08-17 (cgdye, PRD-107)
* **PRD-107 written and stage 0 landed** (`okf/prds/prd-107.md`; plan grilled
  with the user, decisions recorded in the PRD). Hygiene found three real
  defects hidden by stale-path skips: `scripts/rotamer_hgbp1_site481.py`
  shipped with a SyntaxError; `scripts/rrt_hgbp1_site481.py` called a fictional
  `IMP.kinematics.Sampler` API (rewritten on a torsion-space RRT in
  `sampling/rrt_imp.py`); `rotamer/io.load_protein_frames` never called
  `IMP.rmf.load_frame`, so every RMF frame returned the same coordinates.
  De-skipped `test_reference_attachments` (all 33 libraries, placement equals
  the reference transform to 1e-6), `test_rotamer_sampling`, the RRT script
  tests, the integration test (builds its system.cif itself), the RMF test.
  click made optional for the library (`utils.import_click`), the fake
  "Langevin" walk renamed `dof_walk`/`sample-dof-walk` (no shims — hard-rename
  policy), `RotamerFRET` no longer writes a `log` file into cwd. Cheap suite
  `test/cgdye test/fret`: 103 passed/14 skipped → **190 passed/5 skipped**.
* **PRD-107 stage 1:** one CHARMM36 table and one LJ kernel
  (`topology/dye.py::lj_energy`) behind all three cgdye scorers, `LJ_PARAMETERS`
  gone, scalar `evaluate` routed through `evaluate_batch`, unused
  `sampling/rrt.py` deleted (`rrt_imp.py` → `rrt.py`), a latent `_site_mask`
  list bug fixed. `test_lj_single_source.py`. Cheap suite 193 passed / 5 skipped.
* **PRD-107 stage 2 / PRD-106 cgdye half:** `pyext/src/fret/strip.py` is the
  one strip engine (grammar → `StripSelection`; PDB lines, obstacle arrays,
  IMP hierarchies; non-destructive by default); `fret/av.py` imports it,
  cgdye `strip_sidechain_at_site` delegates with its own keep-set (strips CB
  — the linker replaces the side chain from CA; AV keeps CB). PRD-106 status
  → in-progress with a record of what remains (full grammar, diagnostics,
  authored-mask fixtures). Cheap suite 199 passed / 5 skipped.
* **PRD-107 stage 3 — FRETpredict parity.** Scratch env with the local
  FRETpredict checkout: Hsp90 agrees to 5e-8 (the old "1.1e-3 gap" was a
  pin recorded at T=300 vs FRETpredict's 293). pp11 (cutoff10 libraries)
  found a real bug: every library name resolved to the cutoff-30 RMF
  template, so `cutoff10/20` requests silently scored 33 rotamers instead of
  711/123 (Es off by 0.18). Fixed in `rotamer/io.py` (cutoff substituted; the
  FRETpredict DCD sets in `data/rotamer_library` are canonical, RMF templates
  a cutoff-30 fallback); `load_protein_frames` reads multi-MODEL PDBs. Fixtures
  gzip'd under `test/cgdye/rotamer/data/`, pins with provenance in
  `test/references/cgdye_fretpredict_pins.json`, `test_fretpredict_pins.py`.
  `test/cgdye/rotamer` zero skips; cheap suite 209 passed / 0 skipped.
* **PRD-107 stage 4 — physics invariants.** `test_physics_invariants.py`
  (κ², R0, single-rotamer regimes, Boltzmann/mean-field/kinetic weights,
  master-equation limits, Kronecker pair). Found and fixed two physics bugs:
  `compute_exact_efficiency` solved the untransposed system (wrong for any
  non-symmetric transition matrix; fast-exchange 0.878 vs 0.764 correct), and
  the linker Metropolis sampler scored all pairs without 1-2/1-3/1-4
  exclusions while collapsing same-named MOL2 atoms into one site (self-pairs
  at r = 0, ~1e6–1e31 per frame). `dye_internal_system` fixes both; sampler
  pins recorded. Cheap suite 226 passed / 0 skipped.
* **PRD-108 stage 3 — Langevin/Brownian dye samplers; PRD-108 implemented.**
  `LangevinDyeSampler` (md: MD + Langevin thermostat; bd: Brownian dynamics,
  dt ≤ 0.5 fs with the stiff bonds), `dye sample-langevin`, hGBP1 script;
  thermodynamic pins (equipartition, MSD = 6Dt, torsion Boltzmann histogram,
  md ⇄ bd, attached dye at 481). Found and fixed a physics bug: torsion types
  in CHARMM convention were passed to `IMP.core.Cosine` (opposite sign) — planar
  π torsions were minimal at 90°, linkers eclipsed — in `dye_restraints.py` and
  `sim/runner.py` (`topology.dye.torsion_cosine`). Cheap suite 347 passed.
* **PRD-108 stage 2 — AV↔rotamer table.** `dye rotamer compare-av` /
  `cgdye/rotamer/compare_av.py` write `okf/validation/av_vs_rotamer.md` and
  drift pins; medium test with loose sanity bounds. Result: rotamer
  ensembles differ from AVs by 3–12 Å in ⟨R_DA⟩ and fit the 99 T4L distances
  worse (Σχ² 257 vs 155); FRETpredict's soft screening leaves most weight on
  rotamers 1.3–2.5 Å from the protein. C++ `R1` in AVNetworkRestraint is not
  motivated by the data; AV-parameter calibration is the follow-on. Also:
  the `dye` CLI's `main()` never existed (now defined; `rotamer` group attached).
* **PRD-108 stage 1 — fps.json `R1`.** Schema type + fields, C++ `IMP_WARN`
  on non-AV `simulation_type` (test with an R1 position in the T4L score set),
  `fps_positions_for_docking`, `write_rotamer_fps` /
  `distances_from_ensembles` / position payloads; JSON schema regenerated.
* **PRD-108 stage 0 — `RotamerEnsemble`.** `cgdye/rotamer/ensemble.py`
  (an `AccessibleVolume` with per-rotamer centre + dipole + weight, all atoms
  kept; fps `R1`), pair kernels `fret_pair_geometry/efficiencies/distribution`
  in `fret/distance.py`; `RotamerFRET._frame_fret` refactored on top with the
  FRETpredict pins unchanged. Tests + exports. PRD-108 page written.
* **PRD-107 stage 5 — flat `IMP.bff` API, docs; PRD-107 implemented.**
  `pyext/src/api.py` + a `%pythoncode` hook in `pyext/swig.i-in` expose 82
  cgdye/fret names lazily as `IMP.bff.<Name>` (the SWIG module never included
  `pyext/src/__init__.py` — dead file, deleted). Hard renames throughout
  (`System`→`DyeForceFieldSystem`, `calculate_r0`→`forster_radius_from_spectra`
  in `fret/forster.py`, `kappa2_from_vectors`→`kappa2_from_dipoles` in
  `fret/kappa2.py`, ... — full list in the PRD), no aliases;
  `test/test_public_api_names.py` asserts docstring + naming family per
  export. Manual notebook `structure_cgdye.ipynb`, README paragraph,
  `okf/cgdye.md`; notes appended to chisurf PRD-47/93/96. Cheap suite 310
  passed / 0 skipped.
* Note for the build: `pyext/src/Files.cmake` and `test/Files.cmake` are
  regenerated by the IMP configure step (`__init__.py` excluded by design);
  new/renamed Python files need a reconfigure (`cmake -S ../imp -B $B`) to be
  symlinked into the build tree, and stale dangling symlinks there must be
  pruned by hand.

## 2026-08-17
* **PRD-105 pass 7:** structure-of-arrays cloud shared by mean and
  quadrature (serial per-AV compute 2.02 → 1.84 ms/frame, bit-identical),
  scratch/capacity reuse, pair tasks queued in start-ready order.
* **Memory pass, PRD-105 code:** `leaks --atExit` over 14 restraint
  lifecycles in all modes (incl. async, pickle, forced recompute, coarse
  mode): no leaks in libimp_bff (only SWIG's one-time registration blocks).
  Guard Malloc over the full suite, the benchmark (all modes) and the
  lifecycle script: clean. Two ownership hardenings: an unfinished
  `evaluate_async()` job is joined by its own destructor and `job_` is the
  last member of the restraint (destroyed first, before the pool and the
  AVs it uses); the registry's coordinate snapshot is shared with its maps
  (`shared_ptr`) so a map handed to Python cannot dangle if it outlives the
  registry. ASan via `DYLD_INSERT_LIBRARIES` into the conda Python spins in
  `AsanInitFromRtl` on this macOS — unusable here; Guard Malloc +
  `MallocCheckHeap` + `leaks` are the tools that work.
* **Placement (user, 2026-08-17): "split at data level — bff owns structure,
  tttrlib owns data/algos."** The decay routines in imp.bff are the
  deprecated 2.25 copies (PRD-93: they leave in the next release; tttrlib is
  canonical and already carries the same `fconv_per_cs` bounds fix in its
  wrapper). The bff clamp is a stopgap in code that is going away, not an
  investment.
* **Memory: `decay_fconv_per_cs` wrote one past the fit array** — its `stop`
  and `conv_stop` are inclusive (the sibling routines' `stop` is exclusive)
  and the Python wrapper maps `stop=-1` to `n_fit`, so `fit[n_fit]` was
  written and the tail added there: heap corruption that surfaced as random
  aborts later in the same process (Guard Malloc pinned it to
  `test_DecayRoutines.test_fconv_per_cs`). Clamped to the last channel;
  results unchanged. Full suite now clean under Guard Malloc.
* **PRD-105: compilation/locality/arrangement pass.** PGO and branch hints:
  no gain; blocked layout not pursued (working set is cache-resident).
  `AVNetworkRestraint::evaluate_async()`/`wait_score()` added: the pool run
  proceeds on its own thread after the serial Model-reading phase, so the
  caller can load the next frame / evaluate other restraints meanwhile;
  bit-identical, tested; ~1 % on T4L (frame load is 0.04 ms), meant for
  callers with heavier per-frame work.
* **PRD-105 SIMD:** NEON row-form relaxation for interior tiles (~3 %),
  run-based vectorised cloud/penalty passes (noise). Serial per-AV sum ~2.0
  ms/frame; the kernel is bookkeeping/cache-bound, not arithmetic-bound.
* **PRD-105 pass 6 (search kernel):** obstacles encoded in the cost array
  (`BLOCKED_COST`), reused bucket/queue scratch, per-shape candidate list with
  shell/interior classification (occupancy-only tiles skip the sphere test),
  inline locations, `set_origin_fast` wired. Serial per-AV compute 4.1 → 2.25
  ms/frame; threaded frame ~0.7 ms best-of-8 on a loaded box, now with the
  correct 26-neighbour metric. All results bit-identical.
* **PRD-105 correction + memory pass:** the pass-3 "26-neighbour" stencil
  was 18 (`sqrt(3.0)` radius with `d² <= r²` excludes d²=3) — that, not
  tunnelling, caused the 15–25 % smaller AVs and the search halving; fixed
  (`sqrt(3.0)+1e-6`, test asserts 26/30 offsets), pins regenerated; the honest
  26-vs-30 change is 0.42 Å rms (max 3 Å). Speed is back to ~1.0 ms/frame
  quiet. `search_mode="euclidean"` (exact DDA visibility, straight linker;
  0.82 vs 0.77 ms — a model option, not a speed one) added; the cross-check
  it enabled found the stencil bug. Memory: `set_path_map_header` now rebuilds
  the location arrays on shape change (in-place origin refresh wrote the new
  size into old arrays); the full-suite 1-in-4 abort traced to an
  out-of-bounds read inside pip LabelLib (`Grid3DExt::excludeConcentricSpheres`,
  Guard Malloc) — `medium_test_av.py` skips the LabelLib backend unless
  `IMP_BFF_TEST_LABELLIB=1`; 0/8 crashes since. AV lattice suite clean under
  Guard Malloc / MallocCheckHeap.
* **PRD-105: `search_mode="euclidean"` added** (straight-linker AV: exact
  voxel visibility by DDA, Euclidean cost, no path search; option on
  `AVNetworkRestraint`/`AV`, default stays `dijkstra`). Model, not speed:
  ~40 % of the path-search voxels, 5.6 Å rms distance change, 0.82 vs 0.77
  ms/frame. Approximate visibility chains were tried and rejected (leak or
  over-conservative).
* **PRD-105: hybrid Euclidean/geodesic search and obstacle distance
  transform evaluated and rejected** (details in the PRD): visibility-seeded
  Dijkstra was slower with safe (face-step) visibility and tunnelled with
  26-step visibility; the obstacle D-field is worth ~0.03 ms wall and breaks
  local deltas. Tree left at 97cca77.
* **PRD-105 fourth pass ("all tricks", exact):** reach-only re-centring
  shared extents, one coordinate snapshot per frame, chord-run sphere raster,
  one pipelined pool run per evaluation (rasters → searches → carves →
  pairs with dependency waits; prepare reads the pending classification),
  spinning workers, integer occupancy into penalty/carve, interior flags,
  byte visited, in-place origin. Quiet-machine best 0.69–0.71 ms/frame
  (pre-PRD 37.7, ~53×); run-to-run variance on the shared box (0.7–1.0) now
  exceeds the remaining gains. Grid-spacing sweep recorded (2.5 Å: 0.51 ms at
  1.8 Å rms change; not applied — user parameter).
* **PRD-105 third perf pass (objective < 1 ms/frame): reached — 0.81 ms
  median quiet, 0.94–0.99 loaded (pre-PRD 37.7).** Exact: SoA tile arrays for
  the lattice path (`search_lattice`/`carve_lattice`, lazy `tiles` sync),
  bucket queue without duplicates or per-bucket sort (all tiles in an active
  unit bucket are final), fused source spheres, raw-location cloud, two pool
  tasks per AV + z-slab rasters, `timing_ms_total` diagnostics. **Model
  change:** the lattice path now uses a symmetric 26-neighbour stencil
  (`search_stencil=26`, `30` = historical): the historical offset loops ran
  `-2 ≤ d < 2`, so paths could cross a one-voxel wall towards −x/−y/−z only;
  the leak made T4L AVs 15–25 % larger and shifted means by up to 2.5 Å
  (distances 1.1 Å rms, max 10 Å) — pins regenerated, 30-stencil value kept
  as a regression pin, legacy anchoring unchanged. `quad_k` default 100 → 50
  (0.025 Å max error with moments). `search_grid_factor` (coarse search,
  fine carve) implemented and measured: 2.4 Å rms error for only 20 % speed
  → kept opt-in, default 1. Hardware note: 4P+4E cores make per-task CPU sum
  ~2× the serial CPU; compute wall is ~4–5× over serial.
* **PRD-105 second perf pass**: `set_origin` moved into the threaded compute
  phase, shared rasters refreshed on threads, persistent `internal::ThreadPool`
  (longest-first dynamic scheduling), and — the big one — the lattice path now
  runs `PathMap::find_path_dijkstra_bounded` (exact lazy Dijkstra that stops
  at cost ≥ ll/h; source tile left at default cost as before) on a monotone
  unit-bucket queue (edges ≥ 1 voxel ⇒ relaxations land in later buckets;
  each bucket sorted once by (cost, idx) = a heap's pop order). All results
  bit-identical (scores unchanged to the last digit); T4L default mode
  **4.2 → 1.8 ms/frame** (pre-PRD 37.7 → 21×). Fixed a 1-in-10 flake in
  `test_AccessibleVolume.test_distance_distributions` (MC histogram tolerance
  30000 was ~1.5× the expected 2N; now 60000).
* **PRD-105 perf pass** (same day, one commit): the AVs' compute phases run
  on threads (`AV::resample_prepare/compute/finish` split; restraint
  `set_number_of_threads`, default hardware concurrency; pair sums threaded
  too), neighbours enumerated inline in the search (no per-tile edge vectors,
  index heap, penalty mirror), `fill_sphere` without per-voxel Vector3D,
  cloud from reached tiles, locality-aware skip via per-generation change
  boxes on `AVOccupancyMap` (`get_changed_since`), `PathMap::set_exact_search`
  (lazy Dijkstra; agrees with the historical search on every reached tile
  except the source tile, not faster, off by default). T4L default mode
  18.9 → **4.2 ms/frame**, repeats 0.3 ms; all modes bit-identical to before
  (legacy pins unchanged; new thread-equivalence, local-move and exact-search
  tests). One race found and fixed on the way: with legacy anchoring the quad
  caches were not pre-built, so threaded pair sums wrote them concurrently —
  every AV now refreshes its cache in the compute pass before the pair loop.
* **PRD-105 implemented** (`okf/prds/prd-105.md`, status `implemented`;
  implementation record appended to the PRD). New: `AVOccupancyMap` /
  `AVOccupancyRegistry` (integer occupancy counts on the absolute lattice,
  exact subtract/add deltas, integer-voxel rolls with slab raster, grow-on-
  demand shared extent, `read_window` fill-zero-beyond), `AV::resample_lattice`
  (window `2·floor(ll/h+½)+1` centred on `round(s/h)`, skip when nothing that
  feeds the search changed, `force_full` for the exactness proof),
  `av_distance_quadrature` (≤K lattice blocks with centroid, weight and second
  moments; second-order corrected double sum — max error < 0.005 Å at K=100 on
  T4L @2 Å where centroids alone gave 0.26 Å), `AVNetworkRestraint(...,
  space_fixed, shared_map, distance, quad_k)` + `get_diagnostics_json` +
  `get_used_av`, SWIG kwargs shadow for the overloaded restraint constructor
  (kwargs had silently stopped working when the default ctor was added for
  serialization — `test_AVNetworkRestraint.py` was failing on it). Suite:
  `test/test_av_lattice.py` (6 modes, tier-1 bit-exact vs forced full recompute
  and shared vs private, legacy pins via sha256 in
  `test/references/prd105_legacy_pins.json`, order-independence with xfail on
  mc modes, quad K-curve vs the exact double-sum oracle, diagnostics,
  occupancy-map unit tests). Benchmark `benchmark/benchmark_av_screening.py`:
  T4L default 18.9 ms/frame vs corrected legacy 42.4 (pre-PRD code measured
  37.7 today), repeats 0.8 ms. **Two latent defects found in the legacy path
  and fixed for both paths:** `PathMap::update_tiles` had
  `edge_computed.resize(false, nvox)` (arguments swapped → flags resized to
  zero → `get_edges` read stale freed bits → edges never recomputed after the
  first evaluation, so a moved structure kept its old connectivity: scores
  were evaluation-order dependent beyond MC noise, and the old baseline was
  artificially cheap); `AV::resample` carved tiles with `density *= 0` and
  never restored them, so a moving structure's AVs shrank frame after frame
  (`inf` scores late in the T4L trajectory). Legacy single-evaluation results
  are byte-identical to the old code (verified against a stash-build of
  11cb3e7); repeated-evaluation legacy pins are from the fixed code. Existing
  pins in `test_AccessibleVolume.py` / `test_AVNetworkRestraint.py`
  regenerated once to lattice values (mean position shift 0.18 Å, restraint
  score 11.92 → 13.08 deterministic). Deviations from the locked text
  (recorded in the PRD): torus phase indexing replaced by an in-place shift
  of the retained block (same work bound, plain indexing everywhere); quad
  carries second moments. Pre-existing, untouched: `test_DecayConvolution.py`
  / `test_DecayScore.py` fail with the same overloaded-ctor kwargs TypeError.
* **PRD-105 design locked via grilling (15 decisions) — scope rewritten**
  (`okf/prds/prd-105.md`, status `specified`). The original fast-update-only
  scope is superseded: measurement showed all 159 T4L beads move every frame
  (incremental raster can't pay) and rolls fire ~every frame under any window
  policy (padded windows buy nothing). Locked v1: global absolute lattice as
  exactness reference, default on, legacy `space_fixed=False` opt-out
  (deprecated); minimal linker-defined rolling windows (PBC, integer-voxel
  rolls); one shared occupancy raster per (spacing, extra-radius) class —
  17 rasters → 4 classes on T4L, ~14.7 → ~2 ms/frame; cold Dijkstra
  unchanged (warm start = v2); deterministic lattice-quadrature distances
  K=100 replace MC (~14.5 → ~0.6 ms; MC-class accuracy 0.16 Å vs oracle;
  moments rejected at 2.2 Å bias; FFT oracle test-only; numpy timings
  implementation-unfair — C++ A/B required). Distances live in this PRD,
  independent of PRD-94. Two-tier exactness contract (tier 1 bit-exact maps,
  tier 2 deterministic distances); 6 flag modes (`space_fixed`, `shared_map`
  requires it, `distance=mc|quad`), all parameterized in CI with xfail on
  MC order-dependence (pre-existing defect, fixed by replacement); one merge,
  three staged commits; done = suite + tables, no speed gate. All grilling
  A/Bs ran on the real T4L restraint (`chi2_C1_33p`) in the arm64 env.
* **PRD-105 ideas extended — "semi space-fixed" grid** (idea 3 in the
  non-committed section): PBC (torus-indexed) window whose origin is only ever
  relocated by *integer multiples of the spacing* on a fixed global lattice.
  A roll is an O(1) index remap (cached integer occupancy counts stay valid —
  voxel centres stay on the lattice); between rolls the grid is exactly
  space-fixed so the bit-exact subtract/add deltas apply; per roll only the
  newly exposed face slab is rasterised. Wrap-around phantoms are provably
  masked if the window is padded by the reach (they land ≈2·ll from the source,
  inside the `fill_sphere` blocked-beyond-linker region). Conditions recorded:
  no torus-wrapping path edges, hysteresis on the roll trigger, fallback full
  recompute must quantise the origin with the same rule (one shared lattice),
  phase-aware coordinate/export variants. User asked for this variant
  ("pbc + relocation of origin, semi space fixed"); still collect-only, not
  implemented.

## 2026-08-16
* **PRD-106 authored in this bundle** (`okf/prds/prd-106.md`, feature PRD:
  strip mech — a selection-driven mechanism to strip side chains / whole
  residues / atom masks from the obstacle set). Grounding found in code: the
  fps.json `strip_mask` field is declared in `fret/fps_schema.py:163` and
  authored in the shipped T4L + TG2/flex examples (grammar `chain X and resid N
  and not name ...`, `or`, `resname`, parentheses) but has **zero consumers** —
  the AV build substitutes `allowed_sphere_radius` source clearance for
  stripping (`fret/av.py:124-132`); a whole-residue stripper
  `_strip_residue_atoms` (`fret/av.py:557`) is dead code; cgdye's
  `strip_sidechain_at_site` (`cgdye/labeling/attachment.py:73`) hard-codes a
  keep-set (`N CA C O OXT`, strips CB) that disagrees with the AV masks' keep
  (`CA CB C N O`). Scopes one engine in `fret` (grammar-faithful PyMOL subset,
  non-destructive, two outputs — stripped hierarchy + filtered obstacle array),
  live `strip_mask` in the AV build, cgdye delegation with behaviour pinned by
  the existing hGBP1 tests, dead-code removal. Cross-repo consumer: chimol's
  `_stripped_pdb_for` seam (T-20260815-01) hand-derives the same strip today.
* **Ideas collected into PRD-105 as a non-committed design section** (space-fixed
  grid, shared/integrative AV map, better-than-dense-cube data structures). User:
  "ideally make a space fixed grid ... points appear and disappear from grid",
  "move from simple grid to a better data struct"; collect, do not implement.
  Grounding found in code: grid is re-anchored to the source on every `resample`
  (`src/AV.cpp:161-165`, `src/PathMapHeader.cpp:55-74`) so voxel indices are not
  spatially stable across frames; every AV rasterises *all* root leaves into its
  own private cube (`src/AV.cpp:149`) — 17× duplicated on T4L. Ideas: world-
  anchored grid (envelope), rebase policy with integer-spacing-only exactness
  caveat, one shared occupancy grid per frame for all AVs, blocked/chunked grid +
  per-atom voxel footprint cache + hybrid sparse-hash, octree as stretch.
* **PRD-105 authored in this bundle** (`okf/prds/prd-105.md`, feature PRD:
  AV fast update — reuse the previous computation when local geometry barely
  changed). Every `AV::resample()` rebuilds the whole path map (obstacle
  rasterise + full Dijkstra; measured 96 % of an AV build on 148L at 0.5 Å
  grid, `pyext/src/av/compute.py`). The PRD scopes a fast-update path on
  `AV`/`PathMap`: short-circuit when nothing within the label's sphere of
  influence moved, locally re-rasterise + incrementally relabel only the
  changed neighbourhood otherwise, full fallback when the change is too
  large — under an exactness contract (fast == full recompute, or the fast
  path does not run). Consumers: trajectory screening
  (`plot_AVScreening.py`, `imp_engine.screen()`), full-AV sampling,
  `dock_minimize` refinement cycles. Rule (user, 2026-08-16): bff feature
  PRDs live in `okf/prds/` of this bundle, not in the shared chisurf bundle;
  `okf/index.md` scope-rule updated to match.

## 2026-08-11
* **PRD-97 stages 0–3: `IMP.bff.fret` created.** The FRET docking engine,
  the AV backend, six algorithm modules and the single fps.json reader moved
  here from ChiSurf's `fret/core` (details in
  `../chisurf/okf/prds/prd-97.md`). Rules applied during the move (user,
  2026-08-11): **no LabelLib** (IMP.bff's AV is the only backend) and **no
  numba** (the two jitted kernels are vectorised numpy). `fps_schema.py` is
  the authored fps.json definition — flrCIF item names where flrCIF has
  them — and `data/fps_json_schema.json` is derived from it with a drift
  test. `pyext/src/fps.py` deleted (both importers now use
  `IMP.bff.fret.io`). Fixed while landing: `src/AV.cpp` wrote radius1 into
  all three radii (`set_av_parameter`); the moved AV path had two crashes the
  old LabelLib fallback had hidden (ndarray truth test, argument-less
  `DensityHeader.get_origin()`); AV source clearance now scales with linker
  width (`allowed_sphere_radius >= lw/2 + grid/2` unless the position sets
  its own). `examples/structure/TG2/flex.fps.json` `S1_val_chi2` prefixes
  fixed; `hGBP1.fps.json` `577_577` dangling reference flagged, not fixed.
  Tests: `test/fret/` (19) plus cgdye all green — 93 passed, 14 skipped.
* **Bundle created.** `imp.bff` gets its own OKF bundle for repo-local
  knowledge, scoped deliberately narrow: cross-stack knowledge and all PRDs
  stay in `../chisurf/okf/`, photon-level concerns in `../tttrlib/okf/`, and
  the shared agent board is symlinked from `tttrlib/okf/agent-board.md`.
  Initial contents: `index.md`, `overview.md` (the sibling stack — tttrlib,
  imp.bff, imp, imp-tricks, chisurf, fpsimp, quest, ucfret — and the
  where-knowledge-lives map), and `references/index.md` pointing into the
  shared bundles. `AGENTS.md`/`CLAUDE.md` updated to match.

## 2026-08-19 — programs into `bin/`, and what the improper gap was hiding

- `scripts/` and `pyext/src/cgdye/scripts/` folded into `bin/`, IMP's place for
  installed programs: `imp_bff_traj2bcif`, `imp_bff_dye_pdb2cif`. The
  conversion half of the dye converter joined `io.structure`; the one-system
  Langevin driver became an example. Adding a program to `bin/` **disables the
  module** until it has a documented section in `README.md` — `setup_module.py`
  enforces it, and the symptom is a missing `_IMP_bff.so`, not a doc warning.
- A package-wide dead-code audit found essentially nothing (2 candidates in 570
  definitions, one a false positive). Its first three revisions each measured
  themselves rather than the code — see 1d1dfb5.
- Chasing a duplicated exclusion derivation surfaced
  [`impropers_are_dropped.md`](impropers_are_dropped.md): the combined-system
  builder never fills `impropers`, which left 81 of them unbuilt across the two
  shipped components, hid a `TypeError` in `scoring.build_dye_restraints`
  (fixed, 192a763), and is the only reason the two exclusion derivations agree.
