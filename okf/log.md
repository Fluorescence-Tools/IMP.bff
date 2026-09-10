# Update Log

## 2026-09-10

- **The Labelizer LS gets a notebook example** (`ipynb/example/labelizer_score.ipynb`): the label score via
  the native port (PRD-120) as an executed notebook -- `3GUN-processed.pdb` under the no-conservation model
  (T4L ships no ConSurf grades; `cs` dropped, the `--no-conservation` semantics), features + `ll_score_structure`
  in 0.1 s, all 162 residues scored, top site A135 (LS 1.89). One figure: LS along the sequence above the
  per-term heatmap; the structure embedded color-coded by LS (LS in the B-factor column, `3GUN_LS.pdb`,
  py3Dmol/3Dmol.js, red = high); a bridge cell runs `ll_pair_scores` (R0 52, n_refine 5, 4560 pairs) so the
  sites flow into `greedy_pair_selection.ipynb`; and a parity section replays the labelizer.org webserver
  (REST contract reverse-engineered from its SPA: `POST backend/load_pdb {chains,pdbID,payload}` →
  `POST backend/analysis` → poll `GET backend/analysis/{job}` → fetch `{pdbID}_LSlong.csv`; job
  `b6957587-94dd-4727-95bf-a5ca01a34264`, result cached beside the notebook). Measured parity: CR identical
  162/162, SS 156/162, SE 107/162 (same binned table, different SASA engines at bin edges), combined LS
  identical 104/162, median |diff| 0, max 0.60. The webserver evaluates neither tp nor ce ("coming soon"),
  so its LS is the corrected arithmetic -- the published model's zero-veto over weight-0 terms
  (`src/Labelizer.cpp:1404`) zeroes 34 further sites the webserver leaves standing. Its short `*_LS.csv`
  drops rows; only `*_LSlong.csv` is complete.

- **Greedy Olga gets a notebook example** (`ipynb/example/greedy_pair_selection.ipynb`): the pair-selection
  pipeline of `examples/labels/plot_pair_selection.py` as an executed Jupyter notebook, for Olga -- ensemble
  in (`t4l_docking.rmf3`, 100 frames, `fret.fps.json` `chi2_C1_33p`, 33 candidates), `pairwise_rmsd` +
  `select_informative_pairs` (err = 6%), and the requested plot of number of FRET pairs vs resulting
  resolution (expected mean RMSD). On this ensemble the curve answers "how many pairs": pair 1 (70-132_C1)
  buys 1.6 Å of the 4.7 Å prior, the knee is at ~8-10 pairs, and pairs past that re-measure separations
  already resolved. The executed PNG sits beside the notebook; a copy of the figure is in the notebook itself.

- **The command line, compiled** (owner: "all the python scripts in bin end up in a single compiled cpp
  file with corresponding subs ... distributed along with the pip wheel and/or used without the heavy IMP"):
  `include/Bin.h` + `src/Bin.cpp` are that file -- a CLI11 dispatcher (vendored verbatim as
  `include/internal/CLI11.h`, v2.7.2, BSD-3; the consensus header-only CLI, subcommands first-class) whose
  subs are the bin/ scripts that need nothing beyond the core. Two are ported end-to-end: `pdb2cif`
  (the one-liner it always was) and `traj2drot` (convert with weights/--cluster/--grid and the self-check,
  and --bundle; the bulk --all path stays Python for its directory walk). The wheel gains
  `imp_bff = IMP.bff:bin_cli` as a console script; `bin_main` is wrapped for both builds through
  `pyext/IMP_bff.core.i`. The IMP/pmi/RMF-bound modelling commands stay with the click program -- no
  dispatcher makes them lighter; this amends the "a program is Python because a click command is" rule for
  core-only subs, recorded here and in Bin.h. `test/test_bin_subs.py` drives the dispatcher as the console
  script does. First lesson of the port: CLI11's vector parse reads words **reversed** -- the (argc, argv)
  overload is the one that behaves.
- **The drot reads have been broken since the ptolib move -- this session's sub surfaced it, and the reads
  are fixed.** ptolib's `File::read` (so `PtoReader::data`) **decodes by the object's encoding**; the old
  hand-written walker returned raw bytes, and every reader kept its second decompression after the move:
  `read_drot` decompressed already-decoded JSON and failed with "corrupt brotli stream" on every shipped
  library. Fixed across `read_drot`, `drot_catalog`, `drot_provenance`, the Dunbrack records and
  `write_dunbrack_bin`, and the potentials manifest/table reads; `write_drot_bundle` had the inverse bug
  (wrote decoded bytes labelled `+brotli`) and now copies with `read_stored`, which is what its
  "bytes cross unread" comment always meant. The dyes/spinlabels containers read again, and the compiled
  `traj2drot` writes and verifies end-to-end. What remains of the fallout is on the board
  (T-20260910-01): byte-determinism of ptolib containers, the Dunbrack record-size arithmetic, container
  overhead vs the 8x ratio pin -- each with the test that holds it.

## 2026-09-09

- **The codec moved into the container layer** (imp-bff-ce + ptolib 0.4.0; owner ruling "it should be in
  ptolib", then "ptolib.h must change and support a set of standard codecs, all named"): ptolib embeds
  `zstd`, `brotli`, `lz4` and `deflate` in the vendored header's implementation TU — pristine trees under
  `thirdparty/`, spliced by `tools/embed_codecs.py` behind `PTOLIB_CODECS_REGION` fences, compiled hidden
  (the API macros' `visibility("default")` stands down via knob macros and one `[ptolib]` patch to brotli's
  `port.h`), registered under their format names, `PTOLIB_WITH_*`/`PTOLIB_NO_*` as overrides, and a new
  `codec_by_name` so a test can take a built-in out and put it back. The trigger was CI's gcc verdict on
  the private copy: clang had accepted the single-TU brotli's C/C++ linkage conflicts (`kIsBase64` declared
  inside one extern "C" and defined outside another's, ten file-local helpers defined by both fragment
  compressors) and gcc rejected them — repair work that would have been spent on a copy. imp.bff loses
  `src/brotli/` and `src/Brotli.cpp` (~40k lines), the registration shim in `Pto.cpp`, and the window knob
  (ptolib writes lgwin 24 at every quality; 24 was the shipped default) — `DrotEncoding.window` is gone,
  `quality` stays. Pinned and measured: the new codec's q11 bytes are **identical** to the old writer's
  (52/52 on the parity vector), ptolib's golden `.drot` fixtures decode through the embedded copy, and
  `test/io`'s drot suites fail identically before and after (local data incomplete — environment, not the
  codec). `sync_ptolib.sh` refreshed `include/internal/ptolib.h` (6.2 MB); `test_vendored_headers.py`
  unchanged. tttrlib's older copy keeps working; re-sync at its leisure (its board has the note).
- **CI narrowed to what it is for, and repaired** (owner: "ci.yml should ... 1. build bff as imp module,
  2. build bff as standalone pip pkg, 3. test; target platforms: osx, lnx, win"): ci.yml now has exactly
  `imp_module` (conda package, fast+medium tiers), `pip_wheel` (cibuildwheel: ubuntu x86_64, macos-14
  arm64, macos-15-intel x86_64 — macos-13 is retired by GitHub and sat queued forever — and **windows**
  for the first time: MSVC via msvc-dev-cmd, header-only Boost+Eigen from a cached vcpkg, cereal from
  `wheel_deps.sh`, delvewheel repair), `publish_wheels` on a release, and nothing else — the core conda
  package, the weekly expensive lane, the docs build and the anaconda deploy left ci.yml (the recipes and
  workflows to do them locally remain). Three build repairs rode along, each diagnosed from the fork's
  failing runs: the module build could not find `ihm_format.h` (IMP exports the header only inside its own
  tree — the module build now points at the vendored standalone copy); `find_package(Eigen3 3.3)` refused
  conda-forge's Eigen 5.0.1 (the floor became a breaker; the version is dropped, and the whole tree
  compiles clean against 5.0.1, measured locally); and Windows lost `rattler-build` between steps (its
  conda prefix bin dirs now go on PATH in the build step). Fork `tpeulen/IMP.bff` had never had a green
  run — everything was cancelled-by-concurrency or failed at one of the three.

- **A sampling run keeps its path.** `IMP::bff::dock` writes an RMF
  trajectory, a frame each time the walk improves -- not per proposal, which
  would be a file the size of the run. `DockingParameters::save_trajectory`
  now defaults to **true**, and the default is the point: a trajectory is what
  distinguishes a sampling run from a minimisation, and a walk whose path was
  thrown away cannot be judged for mixing, for the basin it settled in, or for
  whether it moved at all.
  Trap: `RmfStructureWriter` is a **value** that owns its file through a
  `shared_ptr`, not an `IMP::Object`, so `IMP::Pointer` will not hold it --
  the error is `no member named 'ref'`, from `PointerBase.h`, which does not
  name the class that is missing it.

- **`read_score_series` returned nothing for the format `count_frames`
  counted.** PMI writes its stat *values* quoted -- `{1: '10.0', 4: '0'}` --
  and `literal_int_map` sorted a quoted value into the *names* map, so every
  data line's numbers went where the column headings go and the value map was
  empty. The reader then skipped each line for want of a score key and
  answered `[]`, on a file whose frames `count_frames` had just counted
  correctly. A quoted value that reads as a number is now recorded as both.
  Found by chisurf's `test_read_score_series`, which had never run: the plugin
  that calls it has been unimportable since August.

- **One sampling interface, four backends -- and docking is now one of its
  callers** (owner, 2026-09-09: "i do not want too diverse sampling
  interfaces, make it such that there are different backends, zeus, emcee
  like, and plain normal mc").
  `Sampler` gains **`"slice"`** (aliases `zeus`, `ensemble_slice`): the
  ensemble slice sampler of Karamanis & Beutler. The halves and the
  differential direction of `"stretch"`, and a slice along that direction
  instead of a Metropolis proposal -- every move accepted, so there is no
  acceptance rate to tune. The header's "deliberately NOT here" note about it
  is retired. Backends are now `metropolis` / `stretch` / `slice` / `de`.
  **`IMP::bff::dock`** is the FRET-restrained sampler, in C++, and it does not
  write a walk of its own: the mobile bodies' poses are one vector -- six
  numbers each, a translation and a rotation vector -- and `Sampler` samples
  it against the network score. So `params.sampler` is a string and swapping
  emcee for slice for plain MC is not a rewrite. **`estimate_docking_errors`**
  repeats a run from independent random starts and reports the spread, which
  is what `bin/imp_bff`'s Python `estimate_errors` did.
  Two traps, both silent, both found by measurement rather than by reading:
  **the slice direction must be the difference of two OTHER walkers.** Built
  from the walker being moved (`partner - self`) the geometry depends on where
  that walker is, the line is not fixed, and the chain is not the target. It
  looks fine: acceptance is 100% by construction, the chain moves, the mean is
  right, and `slice_truncations` stays 0. What it loses is width -- **20% narrow
  on the correlated Gaussian**, and nothing else says so.
  And **the stepping-out must be allowed to overshoot**; a cap that binds
  truncates the tails the same silent way. `get_slice_truncations()` counts
  it, and `test_a_converged_slice_run_truncates_nothing` asserts it is zero.
  A rotation vector, not a quaternion, for the pose: a sampler that proposes a
  quaternion has to renormalise it, and the renormalisation is a projection
  the acceptance ratio does not know about.
  `test/sampler/test_sampler.py` runs its whole parametrised suite over the
  four backends, 48 pass. Its `test_an_unknown_algorithm_is_refused` used
  `"slice"` as the example of a name that is not a backend; that would now
  pass by accident, so it names one that really is not.

- **`write_points_mrc`: the voxelisation belongs with the format, not with the
  application** (owner, 2026-09-09). chisurf's AV export decided how a weighted
  point cloud becomes a grid -- the rounding, the origin, the extent and the
  `.mrc`/`.map`/`.ccp4` suffix rule -- and all four are format questions. They
  are here now, over `write_mrc_grid`, and chisurf names the file and nothing
  else. `rint((xyz - origin) / step)`, so a point sits at the *centre* of its
  voxel, which is the convention the AV raster uses: a cloud that came from one
  lands back on it. Weights sum per voxel; a cloud with no fourth column counts
  each point once. `test/io/test_points_mrc.py`, 17 cases.
  Trap: `%apply` reaches only what is wrapped *after* it, so a typemap for a
  function in `StructureIO.h` has to precede that `%include` -- placed after,
  the build succeeds and the Python signature silently keeps the raw
  `(pointer, n, cols)` triple.

- **`write_mrc_grid`: a dense grid to MRC2014, without `IMP::em`**
  (`include/StructureIO.h`, `src/StructureIO.cpp`). The free `write_mrc` takes
  a `GridHeader` and a `float*` and is `%ignore`d in Python, so the only way
  to write a map from a caller's own array was `DensityGrid` -- one wrapped
  call per voxel -- or `IMP::em`, which is what chisurf did. This is the
  sibling of `write_opendx`: same arguments, other format, C-order in and the
  x-fastest transposition done inside. It is a distinct **name** rather than
  an overload on purpose: SWIG turns keyword arguments off for overloads
  without saying so.

- **`RmfStructureWriter`: a structure as a PMI-shaped RMF trajectory**
  (`include/RmfIO.h`, `src/RmfIO.cpp`, `test/io/test_rmf_structure_writer.py`).
  chisurf built the same tree in Python through `IMP.rmf` -- a decorator call
  per atom per frame, and a dependency the connection layer's IMP does not
  even carry (`rmf` is not in `IMPBFF_IMP_MODULES`). This builds it through
  **RMF's own decorators**, the way `write_rotamer_library_rmf` already did,
  so it is **core**: it needs no IMP at all, only RMF. Per-frame scalars go to
  the `stat` category where `IMP.pmi.output.Output` puts them, as a JSON
  object, so the set of names need not be known in advance.
  It takes a `StructureTable`, which is why that record moved to the core
  (`src/StructureTable.cpp`) with the IMP-backed reader left behind in the
  layer as `StructureReader.h` / `src/imp/StructureReader.cpp`. The record
  also gained setters: a table is *built* by a caller with coordinates of its
  own, not only filled by the reader.
  Verified: chisurf writes a 9315-atom three-frame trajectory through it and
  `chimol` reads it back, in a Python where `IMP.atom`, `IMP.core` and
  `IMP.rmf` are all unimportable.
  Three defects the tests caught, all of which would have been silent:
  **`close()` flushed but did not close** -- RMF has no `close`, a file stays
  open until the last handle goes, and handles hide in the node list and the
  decorator factory as well as in the `FileHandle`, so reopening the path
  answered "This file is currently being written to";
  **a second constructor turned the real one into an overload**, and SWIG
  switches keyword arguments off for overloads without saying so, making the
  Python signature quietly positional-only;
  and **a `stat` key was looked up by type rather than by name** -- RMF
  namespaces keys by type, so a value that arrived as `4` after arriving as
  `4.75` created a *second* key of the same name and every reader asking for
  the first type got nothing back for that frame.

- **`read_structure_table`: IMP's PDB and mmCIF readers behind a flat table**
  (`include/StructureTable.h`, `src/imp/StructureTable.cpp`,
  `pyext/IMP_bff.structuretable.i`, `test/io/test_structure_table.py`).
  chisurf imported `IMP.atom` directly to read a structure, for two things the
  core cannot do: the **radius a docking score measures clashes against** (the
  particle's, after `read_pdb` -- `VdwRadii.h` says why it may not be swapped
  for the core's element-keyed table) and **mmCIF**, which the core does not
  parse at all. Both are here now, as columns -- xyz, radius, mass, bfactor,
  atom/residue id, chain, residue and atom name, element -- with no IMP type
  in the header, so the wheel that links IMP privately (PRD-139) offers them
  in a Python where `import IMP.atom` fails. Verified in that configuration:
  9315 atoms from a PDB and 4585 from an mmCIF, `get_build() == "core+imp"`,
  `IMP.atom` not importable; chisurf then reads structures and computes an
  accessible volume through it with **no IMP import anywhere**, and its
  field-for-field parity gate against the in-tree parser passes 36/36 on three
  structures.
  The walk out of the hierarchy is in C++ on purpose: doing it in Python cost
  a measured ~100 us per atom in SWIG traffic, 1.3 s of the 1.4 s it took to
  read 9315 atoms. Parsing was never the expensive part.
  `standalone/CMakeLists.txt` now names the layer headers that may be wrapped
  like core headers in one list, `_layer_array_doors` (`DyeDynamics.h`,
  `StructureTable.h`), instead of a hardcoded special case for the first one.
  Two traps, both silent: an int out-view **must** be spelled
  `out_view_i`/`n_out_view_i` or numpy's `%apply` does not match and the
  generated property calls a two-argument C++ method with none; and SWIG does
  not re-run when only a *header* changes -- deleting the generated
  `*_wrap.cxx` is not enough either, the entry `.i` has to be touched, or the
  shadow module keeps the previous signature and the build still reports
  success.


## 2026-09-07 (1) — the pair layer could not see a homodimer

Two defects, found while screening labelling sites on BmrA (a homodimeric ABC
transporter) for a two-colour FRET experiment. The screen wanted the only pairs
such an experiment can produce — residue *i* on one protomer to residue *i* on
the other — and the module could not express them.

* **`ll_pair_scores_two_states` matched sites on the residue number alone**
  (`src/LabelizerFret.cpp`, `std::map<int, std::size_t> by_seq`), which is the
  reference's convention (`fret_score.py:512`). On a multimer residue 100 of
  chains A, B and C is one key, the last one inserted wins, and an inter-chain
  pair therefore got the **same site on both sides of the second state**:
  `distance_2` came back 0, `E(0)` is 1, and the score was `jls·|E(d₁) − 1|`.
  Silently, for every such pair, with no status saying so. The map is keyed on
  `ll_residue_key(chain, seq_id)` now, with `LlFretOptions::chain_map` for the
  case the old comment was reaching for — two files that name a chain
  differently.
* **There was no way to say "donor in this chain, acceptor in that one."**
  `LlFretOptions::donor_chain` / `::acceptor_chain` (empty = today's behaviour)
  restore the reference's `chains_apo`/`chains_holo` argument
  (`fret_score.py:479`, filter at `:728`), which
  [`okf/labelizer-correspondence.md`](labelizer-correspondence.md) did not list
  as either ported or dropped. Without it a homodimer screen drowns in
  intra-protomer pairs, and on a file holding two copies of the assembly — 6R72
  has four chains, two dimers — it mixes in lattice contacts that produce
  "distances" of 43–121 Å.

`ll_cbeta_difference_map` has the same collision by construction: it is indexed
by residue number, so it describes one chain. It takes a `chain` argument now
(default empty, i.e. unchanged) and the header says why.

`bin/imp_bff_labelizer` grew `--second-structure`, `--donor-chain` and
`--acceptor-chain`; there was no two-state mode on the command line at all
before.

`test/label/test_labelizer_multimer.py`, seven tests on 3j0e (chains F/G/H).
Three of them fail against the old matching — the collision test reports
`distance_2 == 0.0` for the inter-chain pairs — which is what makes them
regression tests rather than decoration. The whole label suite is 323 green.

The study that found this is `prototypes/bmra/`: a notebook that vendors the
reference Python labelizer, patched, and runs it on cordeshub. Its
`vendor/labelizer/PATCHES.md` records 20 defects in the reference, of which the
one that matters most is `fret_score.py:220` — the joined label score is
`prod**0.5`, correct only for two scores, while the two-conformation case has
four. The port here has that one right already (`ll_joined_label_score` under
`LL_MODEL_CORRECTED`; `LL_MODEL_PUBLISHED` reproduces it on purpose).

## 2026-09-04 (28) — the excitation and emission crosstalk matrices moved here from chisurf

Owner: *"excitation and emission crosstalk matrix definition must be in bff
not in chisurf."* The compute/display line applied to calibration: what a
crosstalk matrix **is** — the label convention, how a light-path payload
becomes values, what forward and inverse mixing mean — is a definition, and a
second definition in Python would disagree with the first silently. ChiSurf
keeps the light-path calculator that *builds* a payload and the views that
display corrected signals; the algebra in between is now
`IMP/bff/CrosstalkMatrix.h` (+ `src/CrosstalkMatrix.cpp`, exposed through
`pyext/IMP_bff.crosstalk.i`), the same placement FcsMdf took.

The surface, in three pieces:

* **`CrosstalkMatrix`** — the labelled matrix both crosstalk matrices share
  (excitation: lasers × chromophores; emission: chromophores × detectors).
  `select(rows, columns)` re-orders against labels and **zero-fills** a
  requested label the matrix does not carry — a missing element of a light
  path is a dark element, not a broken one. `scale_row` / `scale_column` are
  how a per-chromophore quantum yield and a per-detector efficiency fold in.
* **`crosstalk_apply_mixing` / `crosstalk_invert_mixing`** — forward and
  inverse mixing over items (bursts/pixels), with the three solves of the
  Python module they replace: minimum-norm least squares (complete
  orthogonal decomposition, the pseudo-inverse a numpy caller means), the
  Tikhonov closed form, and **non-negative least squares** — Lawson–Hanson,
  ported because chisurf's stable unmixing leaned on scipy and the stable
  path is the one ill-conditioned emission matrices need. 30 random systems
  agree with `numpy.linalg.lstsq`, the closed ridge form and
  `scipy.optimize.nnls` to 1e-6.
* **`crosstalk_shuffle_unmix`** — the photon-reassignment unmixing: each
  detected photon attributed to exactly one source by a binomial chain,
  integer and exactly count-preserving, `std::mt19937_64`-seeded and
  reproducible, unbiased against the NNLS estimate on average.

Two traps met on the way, both silent:

1. **Eigen 5.0.1 mis-evaluates `a * b.col(d).asDiagonal()`.** With the
   diagonal wrapping a **Block**, the product degrades to `a * b(0, d)` — a
   scalar broadcast — so the shuffle's posterior read the abundance ratio and
   leaked a third of every donor photon into the acceptor *with the identity
   matrix*, where the answer is exact by inspection. A dense `VectorXd`
   diagonal is fine; only the Block wrapper breaks. The posterior is an array
   broadcast now, and `test_crosstalk_matrix.py` pins the identity case so
   the idiom cannot creep back. Anyone else in the stack reaching for
   `Block::asDiagonal()` on this Eigen should check their values.
2. **Zeroing an augmented system inside the loop wipes it.** The ridge-augmented
   NNLS zeroes its augmentation rows per item; when there is no augmentation
   (ridge 0) the "tail" *is* the measurement, and the first solve read a zero
   vector and answered zero for everything. Zero once at construction.

Deliberately left in chisurf: the scalar three-cube correction
(`correct_three_cube` / `three_cube_fret_efficiency`) — elementwise numpy on
display-scale arrays, with a batch twin already in tttrlib's
`SpectralCrosstalk`; delegating it is tttrlib's decision, not this move.
`chisurf/core/fluorescence/crosstalk.py` stays as the numpy adapter (payload →
`CrosstalkMatrix` → ndarray, reshape, int64 cast) and owns no arithmetic; its
former `rcond` argument is gone with the numpy pseudo-inverse it fed.

**Same-day follow-up — the tttrlib boundary, checked.** The three-cube
twin did not stay in chisurf long: it now forwards to tttrlib's
`SpectralCrosstalk.correct_three_cube_batch` (the registered,
A/B-validated owner), parity exact. The overlap that *does* concern this
repository: tttrlib's `invert_mixing_ridge` is the same normal-equation
math as this module's ridge path (≤1.5e-10 apart; the hand-rolled Gaussian
elimination is the cruder solver), with no caller outside tttrlib's own
tests. The general mixing algebra is here by the owner's placement, so
tttrlib's function is the superseded twin — deletion is tttrlib's decision
(its `spectral_crosstalk` registry entry names it). Verdict recorded in the
shared duplication register (chisurf `okf/prds/prd-105.md`, phase 4).

### State

`test/test_crosstalk_matrix.py` 15 passed (labels, the three solves against
references, the shuffle contract, and the Eigen-5 identity pin). On the
chisurf side the adapter's whole surface passes unchanged:
`test/fitting/` 1065 passed + `test_crosstalk`/`test_general_correction`/
`test_photon_shuffle`/`test_fret_species_decay`/light-path crosstalk 40+23,
and the backend pin `test/architecture/test_bff_is_the_backend.py::
test_the_crosstalk_matrix_definition_is_bffs` (11 passed) now asserts the
adapter is wired to this module and owns no second definition.

## 2026-09-02 (27) — a port's type is declared once: writes coerce, they do not convert

Owner: *"issue: type conversion. once port init, do not allow. but still
accept types of different kind."*

chinet inferred a port's dtype from every write, so an integer port became a
float port the first time anything stored 1.5 in it. That is what the type
system inherited, and it is the wrong half: a port whose type can change is
not typed. Whatever holds a reference to it has already decided what it is,
and one stray write should not answer differently on their behalf.

**The element type is now declared once** -- by the constructor or by
`set_value_type()` -- and a write of another kind is **accepted and coerced**,
never refused and never promoted. 1.5 into an int port stores 1; 3 into a
float port stores 3.0; 3.7 into a bool port stores `true`. Bool already worked
this way and was carved out as a special case; it is the general rule now.

**Vector-ness is not part of the declaration.** It is shape, not type, so an
array write is still allowed and gives an array port of the same element type.

### Three things that had to move with it

1. **A bare `Port()` now defaults to float, not int.** chinet could default to
   the integer code because the first write retyped it and the default never
   survived contact with data. With a declared type that default silently
   truncates the first float anybody stores -- `test_port_bounds` caught it
   immediately, writing a float array into a fresh port and getting integers.
   A float port is what almost every caller wants anyway.
2. **Both loaders restore the type *before* the value.** Restoring the value
   first would coerce it to the default type and lose it -- a saved 2.5 comes
   back 2, and the type restore that followed could not put the half back.
3. **The type is then applied a second time, after the value.** That is not
   redundancy: the first application makes the value coerce correctly, the
   second restores the document's exact code, including chinet's quirk that a
   constructor-built float *vector* reports the scalar code. A retype within
   the same element type converts nothing, so it costs no precision.

### What it buys

Exactness stops being fragile. Under inferred promotion one float write turned
an exact int64 port into a rounding one for the rest of the session; now the
port is still an integer port afterwards and the next wide value is still
exact. `test_a_wide_integer_port_stays_exact_across_a_float_write` says so.

### Divergence, stated plainly

This is a deliberate break from chinet's semantics, not a port of them. The
`reference_chinet` parity tests that pinned promotion were **inverted, not
deleted**, and each says in its own docstring what changed and why:
`test_port_value_type`, `test_an_int_port_accepts_a_float_and_coerces_it`,
`test_a_wide_integer_port_stays_exact_across_a_float_write`. `set_value_type()`
survives as the one declaration path, which is why loading a document uses it.

### State

`test/portnode/test_port_types.py` 56 tests, including the full accept/coerce
matrix over the three element types. bff 355 passed, chisurf 1468 passed (the
pre-existing pair).

## 2026-09-02 (26) — one slot buffer and a tag, replacing the double/int64 pair

Owner: *"why no storage buffer in bytes and type casting?"* -- a better
question than my answer to it deserved. (25) stored an int64 array **and** a
double mirror, and the answer was that `get_values_ref()` returns
`const std::vector<double>&` with no copy to ten hot-path call sites. That was
a real constraint and the wrong conclusion: the mirror created a bug class
(*which store is authoritative*) that had already cost three bugs in one
session.

### The shape that keeps both

Storage is now **one** `std::vector<double> buffer_` used as a slot array plus
the element type in `value_type_`:

* a float slot holds the double itself;
* an int or bool slot holds the **int64 bit pattern**, copied in and out with
  `memcpy`.

A double slot is already 8-byte aligned and exactly as wide as an `int64`, so
there is no alignment question, and `memcpy` keeps it clear of strict
aliasing -- no `reinterpret_cast` of one live object as another.

The zero-copy contract survives untouched. For a **float** port -- which every
one of those ten hot-path readers is -- `get_values_ref()` returns `buffer_`
itself. Nothing downstream changed, no span refactor, no signature churn.
For a non-float port it materialises `double_cache_`, which is **derived**:
never read back into `buffer_`, dropped on every write. A cache cannot become
a second authority; a mirror can, and did.

`bench_fit.py` after: graph TCSPC **2.53 ms**, FRET **6.18 ms**, against 2.26
and 6.21 recorded before the change. The hot path is where it was.

### What the single store fixed for free

`set_value_type()` used to need a rule for which array to trust -- the source
of (25)'s worst bug, where restoring a document destroyed the value it had
just restored. It now reads both forms out of the one buffer, retypes, and
writes back, so an int -> int retype is exact by construction rather than by a
special case.

`propagate_to_followers()` was pushing the double array into every follower;
for an integer port that rounded on the way into the follower's own storage --
the same mistake one layer further out, and nothing had caught it. It pushes
typed now.

### The one new risk, and its tests

A derived view can go stale. Five tests cover it: an integer port reads as its
*value* and not as its bit pattern (a reinterpreted slot would come back a
denormal near 1e-308, not 5), the view follows writes, follows a type change,
follows a link, and a wide integer survives being pushed to a follower.

Memory, incidentally: an integer port was 16 bytes an element (double plus
int64) and is 8, with the double view allocated only if something asks for it.

### State

`test/portnode/test_port_types.py` 45 tests. bff 339 passed, chisurf 1468
passed (the pre-existing pair). Adding a fourth element type is now a tag and
a conversion, not a fourth array to keep in step.

## 2026-09-02 (25) — Port types become real: an integer port now holds an integer

Owner: *"make ports have types."* -- said again after (24), which is the tell:
(24) gave the port a type *code* and a typed read, and left the storage a
`double` for every type. That is a label, not a type, and I had written the
caveat myself. This is the part that was missing.

### The two defects, both demonstrable

```
Port(value=2**53+1).value  ->  9007199254740992   type=float
```

Wrong twice in one line. **The type was wrong** because a Python int wider
than 32 bits does not convert to C++ `int`, so SWIG fell through to the
`double` overload and the port declared itself float. **The value was wrong**
because storage was `double` regardless.

### What it is now

`Port` carries `std::vector<long long> int_data_` as the exact store for the
int and bool element types, with `data_` kept as a **double mirror**. The
mirror is not laziness: `get_values_ref()` hands `data_` out with no copy and
the evaluation hot path -- `ChiSquared`, `JointChiSquared`, `Minimizer`,
`TcspcDecay`, the spectrum nodes, ten call sites -- reads it every iteration.
Those are all float ports, so they keep a branch-free zero-copy read and pay
nothing for a type system they do not use.

Every integer entry point widened from `int` to `long long` rather than
gaining an overload -- widening cannot be mis-dispatched, and (24) had just
been bitten by exactly that. `Node.cpp`'s two `static_cast<int>` writes went
with it, which also stops an expression yielding 3e9 wrapping through 32 bits.

Exact to `INT64_MAX` on construction, on assignment, through a link, through
`Port.get_json`, and through a saved `Session`.

### Three places the exactness leaked, each found by a test rather than by reading

1. **`set_value_type()` re-derived the integers from the double mirror.** So
   saving and reloading `2**60 + 7` came back `2**60`: the value survived the
   write and was destroyed by the type restore that followed it. Conversion
   now reads whichever store is authoritative.
2. **The kwargs constructor sent every list through `vector<double>`**, so a
   wide integer *vector* was rounded before it ever reached the int store,
   even though the scalar path was already exact.
3. **`Session` emitted through `get_value()` and read back through
   `get<double>()`** -- a path that never touches Python, and rounded twice.

### Two things worth keeping in mind

`Port(value=2**60+7)` then `p.value = 1.5` gives a **float** port holding
`1.5`, and that is correct: exactness is a property of the int type, not of
the value, and chinet's int -> float promotion is deliberate. What would be
wrong is promoting silently while still claiming to be an int port, and a test
says so.

A numpy `bool` array now builds a `PORT_BOOL_VECTOR` rather than an int
vector -- it was `dtype.kind in 'iub'` all in one branch before.

### State

`test/portnode/test_port_types.py` 40 tests. bff 339 passed, chisurf 1468
passed (the pre-existing pair). The `Port.h` file comment no longer says an
integer port is exact only to 2^53, because it no longer is.

## 2026-09-02 (24) — Port gets a type system: int, float, bool, each as a vector

Owner: *"improve bff Ports of Nodes, also vectors should be possible as input
and output, make Port have types, support int, float, bool."*

Vectors were already there -- constructor, `set_value_vector`,
`set_values_array`, `get_value_view`, and links propagate them -- so the real
gaps were **bool**, and that the "type" was a bare int (0..3) that a port did
not actually answer in: `Port(value=True).value` was `1.0`. A type nothing
reads back in is half a type.

### What the codes could and could not become

`Session` writes the code into the chinet document verbatim and
`test/session/chinet_fixture.jsonl` pins it, so **0-3 keep chinet's meaning**.
Bool took fresh codes (4, 5) rather than a renumbering, which leaves every
existing document readable. A named `PortValueType` plus
`port_value_type_element/_is_vector/_of/_name` replaces comparing magic
numbers by hand.

### Bool does not promote, and that is the point

int -> float is an *inferred* promotion: chinet takes numpy's dtype rules, so
writing 1.5 to an integer port makes it a float port. A bool port is a
**declared** type -- a flag -- so writing 3.7 to it stores `true` and leaves it
a bool port. Promoting a flag on assignment would silently destroy what the
port means, which is the opposite of what a type system is for. Only
`set_value_type()` takes a port out of bool. The int -> float rule is
untouched, and a test pins each separately.

### Three traps, all found by running rather than reading

1. **`bool` is a subclass of `int` in Python**, so `hasattr(True, '__index__')`
   is true and the integer branch swallowed every flag. `_set_value` tests
   `isinstance(v, bool)` first.
2. **The `std::vector<int>` overload made SWIG decide the element type, and it
   decided wrong.** A Python list of floats converts to `std::vector<int>` as
   happily as to `std::vector<double>` -- lossily -- so `Port(value=[1.5,2.5])`
   came back `[1, 2]`, and two existing tests caught it. The fix is to not
   overload across int/double on the binding surface at all: the vector
   constructor is gone and the method is named `set_value_vector_int`. The
   `value_type` argument already said what the caller meant.
3. **`true` is not a JSON number**, so the session loader dropped a saved flag
   and `set_value_type(4)` then read it back as `false`. The loader restores
   the value *before* the type, so that branch has to accept booleans itself.
   `emit_port_value` writes real `true`/`false` for bool ports; 0-3 render
   exactly as before.

### What did not change, deliberately

`Port(value=[1.5, 2.5])` still reports the *scalar* float code, not
`PORT_FLOAT_VECTOR` -- a chinet quirk, load-bearing because
`test_kwargs_ctor_chisurf_shape` pins it as the shape chisurf's `Parameter`
builds. Asserted rather than corrected. Likewise float stays sticky upward:
`set_value_vector_int` on an already-float port leaves it float.

Storage is still `double`, so an integer port is exact only to 2^53. Closing
that needs variant storage and was not part of this; the file comment now says
so instead of claiming more than it does.

### State

`test/portnode/test_port_types.py`, 26 tests. bff 325 passed, chisurf 1468
passed (the two failures are the pre-existing pair this session opened with).

## 2026-09-02 (23) — the Gaussian mixture in one call, and the third kernel hiding behind a flag

Owner: *"transfer to bff also distance btw gauss, so that all compute can
happen within bff and minimal cpp python transfer."*

The kernels were already C++. The *loop over them* was not: `Gaussians.
distribution` called into bff **once per component** and summed the returned
arrays in numpy, so a five-component model crossed the language boundary five
times to build one curve, each crossing returning a 96-point ndarray. That is
the half of the standing rule -- *the data stay where the computation is* --
that "it is already a C++ kernel" misses.

**`gaussian_distance_mixture` in `Distributions.h`** takes the axis and the
per-component arrays and returns the finished mixture. Five components ->
**one call**, measured. `GaussianDistances::evaluate` was rewritten to call the
same function, so the graph path and the Python path are now one
implementation and cannot drift; `test_the_node_and_the_free_function_agree`
pins that in both branches, at exactly zero.

### The graph stopped refusing

`minimizer.py` returned `None` for `is_distance_between_gaussians`, because
the node carried only the generalised normal — so `GaussianModel` had **no
graph at all** whenever that flag was on. The node has both kernels now and
the builder wires the flag through. `test_a_distance_between_gaussians_
refuses_the_graph` became `..._builds_the_graph`, and rather than asserting a
refusal it now checks the graph's curve against the model's (1e-8).

### Three kernels, not two — found by measuring

The obvious reading is that pda2c's plain Gaussian is the generalised normal
at `shape = 0`. It is not, and the mixture was briefly wrong because of it:

* `generalized_normal_distribution` evaluates the **standard** normal at
  \f$z=(x-\mu)/\sigma\f$, so it omits the \f$1/\sigma\f$ factor;
* `normal_distribution(norm=False)` keeps it.

With `norm=True` the constant divides out and the two agree — which is why
this hides. With components of **unequal width** left unnormalised it does
not: a **2e-2** error on a three-component mixture. So the kernel is chosen by
an explicit `GaussianMixtureKernel`, and a test asserts the two branches
*differ*, so nobody folds them back together.

The same care applies to `normalize_components`, which is a parameter rather
than a constant because **the reference is asymmetric between its own
branches**: the generalised-normal branch calls a function whose `norm`
defaults to true, the two-cloud branch one whose `normalize` defaults to
false. Hiding that in the kernel would silently change one of the two callers.

### Numbers

5 components -> 1 crossing per distribution (was 5), 53.5 us (generalised) and
43.3 us (two-cloud). Both branches agree with an independent numpy reference
built from the definition at **5.6e-17** and **0**. bff 125 passed; chisurf
1468 passed. Census unchanged.

Pre-existing and unrelated, confirmed not caused here: `mmfdb.repository` and
`IMP.bff.av` missing (import errors), a Qt abort in
`test_fret_structure_widget`, and the two failures this session started with.

## 2026-09-02 (22) — dedup: the last six copies, and a third `distance_between_gaussian` that was a different function

Owner: *"Dedup."*

Asked of the code rather than of memory: intersect every function name in
`chisurf.core.math` against `dir(IMP.bff)`, then check which of the shared
names still carry their own Python body. **14 shared names, 6 still
duplicated.** Each was compared numerically before being touched.

| | |
|---|---|
| `distributions.poisson_0toN` | forwarded, agreed to 1.2e-16 |
| `distributions.normal_distribution` | forwarded, 5.2e-16 |
| `distributions.generalized_normal_distribution` | forwarded, 5.3e-15 |
| `distributions.distance_between_gaussian` | forwarded, 0 |
| `special.i0_array` | forwarded, 3.8e-20 |
| `datatools.distance_between_gaussian` | **not a duplicate — see below** |

`sum_distribution` went too: `combine_distributions`' own docstring called the
two "functionally equivalent", and nothing called it.

### The one that was not a copy

`chisurf.core.math.datatools.distance_between_gaussian` shares its name with
two other functions in the same codebase and **is not the same function**. It
returns a plain Gaussian centred at the separation; the other two return the
distribution of the *distance between two Gaussian clouds*, which carries an
extra \f$r/s\f$ factor and an antisymmetric second term. On a typical axis
they are **0.93 apart in relative terms** — different functions, not copies
that drifted.

Forwarding it to bff on the strength of the matching name would have silently
replaced it. This is the same confusion that produced the worm-like-chain
linker bug one entry above, found a second time, in a second place, on the
same day. A name is not evidence.

**I removed it, and that was wrong**: the grep for callers covered `chisurf/`
and missed `test/math/test_datatools.py`, which tests it. It is restored, with
the collision written into its docstring. Renaming it to say what it is would
end the collision for good and is worth doing — it is left as a decision
rather than taken, because it is a public name.

### State

No function in `chisurf.core.math` now shares a name with an `IMP.bff` export
and still computes it in Python. Dead imports left behind by the forwarding
(`gamma`, `log`, `i0` in `rdf.py`; `math` in `distributions.py`) removed.
bff 95 passed, chisurf 244 passed, census unchanged. `Qd` and `linear_dist`
remain: dead, but not duplicates, so not this ticket's business.

## 2026-09-02 (21) — the worm-like chain, the distributions, and a bug the move uncovered

Owner: *"the compute of wlc and ising should be moved, also move the underlying
distributions to bff."*

### The move found bff computing the wrong function

Forwarding ChiSurf's `worm_like_chain` onto bff's meant comparing them first,
and they are **not the same function**. bff had

```cpp
pri *= std::exp(-d * kappa * a * (1.0 + b) * r_n / (1.0 - (b*r_n)*(b*r_n)));
```

where Becker-Rosa-Everaers -- and ChiSurf -- have \f$I_0\f$ of that same
argument. The argument is **negative** and \f$I_0\f$ is **even**, so
\f$I_0(-x)\f$ grows exactly where \f$\exp(-x)\f$ decays: a factor of 4e3 at
\f$x=-5\f$, 2e16 at \f$x=-20\f$.

It is not a precision bug, it is the wrong shape. Against the exact
Kratky-Porod \f$\langle R^2\rangle/L^2 = 2\kappa - 2\kappa^2(1-e^{-1/\kappa})\f$:

| \f$\kappa\f$ | exact | with `I0` | with `exp` |
|---|---|---|---|
| 0.05 | 0.095 | 0.077 | 0.070 |
| 0.50 | 0.568 | 0.421 | 0.090 |
| 2.00 | 0.852 | 0.770 | **0.007** |

A stiffer chain came out *more compact*. **Nothing outside this module's own
tests called it** -- ChiSurf ran its own Python copy -- so it was latent, and
the one change that would have made it live is the change that found it. That
is the "nobody compared them" hazard from `T-20260901-05`, a second time.

Two more divergences in the same family:

* `worm_like_chain_linker` convolved with `normal_density`. The broadening
  kernel is `distance_between_gaussian` -- the distance between two Gaussian
  *clouds*, which carries an extra \f$r/s\f$ factor and an antisymmetric
  second term -- not the density of an offset. It also fed the r^2-weighted
  chain into the convolution, which must see the bare one.
* The contour cut was a **mask**; the reference `break`s, which is a
  **prefix**. On an unsorted axis those differ, and nothing requires the axis
  to be sorted.

**The test that should have caught the kernel had confirmed it instead.**
`test_the_linker_convolution_matches_an_independent_one` built a numpy
reference and made the same two mistakes as the code. An independent reference
is only independent if it is derived from the definition rather than written
next to the implementation; it is now, and it says so.

### What moved

New `include/SpecialFunctions.h` / `src/SpecialFunctions.cpp` -- `i0`, the
Abramowitz & Stegun polynomial. **Including ChiSurf's transposed digit**:
A&S print 3.5156229, ChiSurf has 3.5156299, and the C++ carries 3.5156299 on
purpose. The contract of a port is "the same answer, faster", and every
published WLC fit in the stack was made with that constant. It costs ~1e-6
near |x|=3.75. If it is ever corrected it should be corrected once,
deliberately, with the refits that implies. **Flagged, not fixed.**

`saw_nu` added to `PolymerChain.h`. Everything now forwards from
`chisurf/core/math/functions/{rdf,special}.py`: `i0`, `gaussian_chain`,
`gaussian_chain_ree`, `saw_nu`, `worm_like_chain`, `worm_like_chain_linker`,
`distance_between_gaussian`. Every one is parity-exact against the pre-port
code taken from git -- worst case 3.4e-15, most at 0.

The distributions that already had bff twins were checked rather than assumed,
and all agree to machine precision (`normal_distribution`,
`generalized_normal_distribution`, `distance_between_gaussian`, `poisson_0toN`,
`gaussian_chain`): 0 to 3e-17. Only the worm-like chain had drifted.

### What did not move, and why

* `combine_distributions` takes an arbitrary **callback**. Porting it needs a
  SWIG director per element, which `okf/log.md` 2026-09-01 already measured as
  a regression (1.54 ms against scipy's 1.33). The specialised case that
  matters is already a node -- `GaussianDistances`.
* `Qd`, `linear_dist`, `sum_distribution`, `i0_array` (the ChiSurf one) have
  **no callers anywhere in the stack**. Porting dead code buys nothing.
* `worm_like_chain`'s `distance` flag has **never done anything** in ChiSurf --
  accepted and dropped. bff implements it properly, so the forwarder passes
  `False` explicitly to preserve the existing answer. Honouring it would change
  every WLC fit. **Flagged, not fixed** -- and the Kratky-Porod table above says
  the r^2 form is the better-behaved one, so this is worth a decision.

### Numbers

`IsingChainModel` 55215 -> 1090 us (**50.7x**). `WormLikeChainModel` 266 -> 210
us (1.3x), `SawNuModel` 224 -> 211 us (1.1x) -- small, and
`test/minimizer/bench_models.py` said in advance that they would be: those two
were 2.0 and 1.5 ms of movable compute between them. They were moved for the
placement rule and for the bug, not for speed, which is the honest way to
describe it.

bff suites 239 passed (27 in `test_polymer_chain.py`, including a regression
test that pins \f$\langle R^2\rangle\f$ *rising* with \f$\kappa\f$ rather
than any particular number). Census unchanged: no model builds a graph that
disagrees with its own curve.

## 2026-09-01 (20) — where the model compute actually is, and the Ising chain to C++

Owner: *"benchmark chisurf mdls and move compute to bff."*

### The census was ranking the wrong thing

`test/minimizer/census_models.py` answers *whether* a model becomes a C++
graph. Both open tickets in this family (`T-20260901-08`, `-09`) were scoped
off it, and it cannot rank: a model nobody fits scores the same as the one
every session runs. So the second half got written —
**`test/minimizer/bench_models.py`** — which times one `update_model()` per
model and profiles it into three buckets: `native` (inside a `tttrlib`/`IMP`
extension), `numpy`, `py`. It ranks by `us/LM-iter` — cost times
`(n_free + 1)`, because LM differences the objective once per free parameter
— restricted to the share that is **not** already C++.

Three things fell out, and two of them contradict the tickets:

1. **The 7 models that build a graph are the cheap ones**: 4.1 ms of 2174 ms
   of total per-LM-iteration compute, **0.2%**. The graph work so far bought
   correctness and structure, not speed. That is worth knowing before the
   next session promises a speed-up from finishing it.
2. **`MaxEntLifetimeModel` is the most expensive model in ChiSurf by a factor
   of three (~190 ms a curve) and has nothing to move** — 97% of it is
   already inside `tttrlib.solve_tcspc_mem_lifetime`, which is where the
   layering rule puts it. A benchmark that only measured *time* would have
   sent someone to port it. This is the entire reason the split exists.
3. **The distributions `T-20260901-08` scopes are 0.2% of movable compute**
   between them. Worth doing for the graph; not worth doing for speed.

### The one that mattered

`ising_chain` — 619 ms of the 1192 ms of movable compute, **52% of it, in one
function**. `chisurf/core/math/functions/rdf.py` looped over 2000 k-points and
stepped a 2-vector through `n` 2x2 matrix multiplies per point: ~80k
interpreter iterations of arithmetic numpy has no way to vectorise. It went to
`include/PolymerChain.h` / `src/PolymerChain.cpp` beside the worm-like chain,
and `rdf.ising_chain` is now a thin forwarder on the
`kappa2_to_distance_ratio` pattern.

| | before | after |
|---|---|---|
| kernel | 55.6 ms | 0.92 ms (**61x**) |
| `IsingChainModel.update_model` | 55.2 ms | 1.14 ms (**48x**) |
| its rank in the table | 1 | 10 |
| **all movable model compute** | **1192 ms** | **587 ms** |

Parity **5.6e-17** against the numpy original over five parameter sets, and
**6.9e-17** on the model's own axis against the pre-port code taken from git —
which is the comparison that actually guards the port, since it is the one a
user would notice. 7 tests added; polymer + minimizer suites 74 passed.

**The transfer product is formed unscaled, deliberately**, and the header says
why: `phi` is used only as the ratio `phi(k)/phi(0)`, so a rescaling cancels
only if it is identical at every `k` — a per-`k` rescaling silently changes the
answer. The cost is that a long, strongly coupled chain can overflow; those
entries end up zero rather than NaN, and a test pins that, because a NaN there
propagates into the FRET rate spectrum downstream.

### Two cautions the benchmark carries in its own docstring

**It over-reports `py` for linear algebra.** Functions are bucketed by the file
they are defined in, so `np.linalg.inv` — a Python wrapper over LAPACK — is
charged to `py`. `DeerTikhonovModel` reads 99% `py` and a real part of that is
BLAS. A high `py%` on a linear-algebra model is a question, not a verdict. The
`native` column has no equivalent failure: extension calls are always
attributed correctly, so a row reading `native` really is done.

**`DeerTikhonovModel` is the noisiest row** — 46 to 232 ms a curve on an idle
machine, without growing within a process. Its rank is real; its absolute
number is not, and a single-run before/after on it proves nothing.

Both are why `T-20260901-15` says to fix that model's arithmetic before
porting it: `_gcv_score` recomputes `A.T @ A` and `L.T @ L` inside the loop
over all 24 alphas though neither depends on alpha, and forms a full
`A @ inv @ A.T` to take its trace. That is a cheaper experiment than a port.

## 2026-09-01 (19) — the transport nobody counted: the graph was rebuilt once per consumer

Owner: *"minimize transport btw bff and chisurf."*

### What was crossing

Counted rather than guessed — `graph_objective` calls,
`ChiSquared.set_data_arrays` calls and `Port.set_values_array` calls during
one `fit.run()`, and then during one `covariance_matrix()`:

```
tcspc   run():  graph builds 1 | set_data_arrays 1 | update_model 2
        + one covariance_matrix(): graph builds 1, set_data_arrays 1
```

`run()` is clean. **`covariance_matrix()` was not**: every call rebuilt the
entire graph and re-copied `x`/`y`/`ey` into `ChiSquared`. Timed:

| | graph build | whole call | build's share |
|---|---|---|---|
| tcspc | 175 µs | 249 µs | **70%** |
| FRET | 523 µs | 978 µs | **53%** |

And there are **six** consumers of that curvature (the posterior view,
derived quantities, both sampler preconditioners, `Fit.grad`, the property),
so the data were marshalled once per consumer. This was introduced by (17),
which routed all six through the graph and made each of them build one — the
arithmetic moved to C++ and the transport moved with it, which is precisely
the half of the rule that gets forgotten.

### The fix, and the bug inside the fix

`_cached_graph` builds once and keys the result on what must not have
changed. **The first key was wrong in an instructive way**: it included
`id(fit.data.x)` and `id(fit.data.y)`, and `DataCurve` builds a *fresh array
on every access*, so those ids are the address of a temporary CPython
immediately reuses. Two reads in one expression compared **equal**; two reads
either side of a call compared **unequal**. So the naive stability check said
the key was fine and the cache never hit once — visible only because the
measurement did not move.

The key is now structure and window (`id(model)`, `id(fit.data)`, the
free-parameter identities, `xmin`, `xmax`), and `Fit.run()` drops the cache,
so a graph cannot outlive the run it was built for. **Uncovered, and written
into the docstring rather than left implicit**: mutating a data buffer *in
place* between two analyses of one finished fit. Fingerprinting the contents
would mean reading the whole array per call, which is the transport being
removed.

### What it is worth

`covariance_matrix()`: tcspc **249 → 72.8 µs (3.4x)**, FRET **978 → 448 µs
(2.2x)**. Per fit that is one graph build instead of up to six, and one copy
of the data instead of up to six.

`test_the_graph_is_not_rebuilt_for_every_curvature` pins the count *and* that
the cached answer is bit-identical to the first — a cache that returns a
different answer is worse than no cache. `test_a_rerun_does_not_reuse_the_previous_graph`
pins the invalidation.

### Where transport now stands

Per `fit.run()`: **one** graph build, **one** data copy, **two** Python
`update_model()` calls. The two remaining known crossings are unchanged and
both already ticketed: the post-fit curve is recomputed for display rather
than read off the port (`T-20260901-11`), and `DataCurve` still holds
`x`/`y`/`ey` in numpy beside the node's copy — one copy per fit, so "the data
stay local" is close but still not literally true.

### Tests

`chisurf test/fitting` **1003 passed** plus the same two pre-existing
failures, and the two new tests above.

## 2026-09-01 (18) — the reconvolution stops waiting for itself: 5.9x on the kernel, 3x on a FRET fit

Owner: *"bff and tttrlib are the engine, chisurf is the car. move more to
bff cpp to make the engine go fast."* Ticket `T-20260901-12`, opened by (16).

### The one-line version

`fconv_per_cs_ad` computed one species at a time, and each species is a
**serial dependency chain** over the channels. The species are independent of
one another, so eight of their recursions now run interleaved and the
pipeline is full instead of empty. **A FRET fit is 18.6 -> 6.18 ms.**

### Measured before it was written, which is why it was written

A standalone A/B of the serial kernel against a blocked one, min-of-many,
interleaved in one process, 512 channels:

```
species    serial     B=2      B=4      B=8
      1    2.0 us    1.9      2.0      2.5
      2    3.9 us    1.9      2.0      2.5
      4    7.8 us    3.8      2.0      2.5
     53  103.2 us   50.3     27.6     17.4     <- 5.9x
     97  188.2 us   91.2     49.3     32.4     <- 5.8x
```

The prediction in (16) was ~62 µs of the 106 was FMA latency; recovering
86 µs of it says that was right. **Plain `-O3` gives the same figures as
`-mcpu=native`**, so this is ILP and not autovectorisation — it should hold
on x86, which matters because x86 has no AVX kernel for this variant and
takes the scalar path for the plain-`double` case too.

### What it cost, stated plainly

Interleaving changes the **order** the per-species contributions are summed
into `fit[i]`: a block is summed and added once rather than each species
being added in turn. Measured deviation from the serial body: **5e-16**
relative to the curve peak, against the 1e-10 to 1e-14 the callers pin
curves at. That is a few ULP, and it is a real change rather than no change.

**Below `FCONV_AD_BLOCK_MIN = 2` the serial body is kept**, for two reasons
that happen to agree: a one-species spectrum has nothing to interleave and
B=8 is *slower* there (2.5 vs 2.0 µs), and keeping it means **every
single-exponential result in the library is bit-identical to what it was**.
`DecayFitNExp` calls this with `numexp = 1`.

### Where it landed, and where it did not

**tttrlib owns the original and the change is there**, not here: the kernel
is `tttrlib/modules/spectroscopy/decay/include/DecayConvolution.h`, its test
is `tttrlib/test/cpp/test_fconv_interleave.cpp`, its registration is
`tttrlib/CMakeLists.txt`, and its own write-up is `tttrlib/okf/log.md`
2026-09-01. This entry records the *downstream* effect. `fconv_per_cs_ad` is
now a dispatcher over `fconv_per_cs_ad_serial` (the old body, kept and named)
and the blocked one, with the test beside it — bit-identity below the
threshold, a few ULP above it, and a check that appending zero-amplitude
species changes nothing, which is the failure a padding lane would cause.
Vendored into bff by `cp`; `test_decay_convolution_copy_is_identical.py`
still passes.

tttrlib's *compiled* library was rebuilt too, which was not free to skip:
its editable install is `editable.rebuild=false`, so the extension was stale
by design and `DecayFit23/24/NExp` would have kept the old kernel
indefinitely. Rebuilt through the supported hook
(`tttrlib.__loader__.rebuild()`); all four decay translation units
recompiled, and `tttrlib test/python/decayfit` is **172 passed, 1 skipped**.
`test_ad_gradient` passes against the new header — the `Dual` instantiation
and DecayFit23/24's objective gradients, which is what (12) asked to be
checked — and `test_fconv_interleave` passes beside it.

Because that install is **shared**, both suites were then re-run against the
rebuilt library: imp.bff 1830 passed, chisurf `test/fitting` 1003 passed plus
the same two. A rebuild under a sibling checkout is exactly what breaks it
silently, so testing only the repository the change was made in is not
enough.

### The fit, end to end

`fit.run()`, best-of-many, graph path:

| fit | before (17) | now | |
|---|---|---|---|
| parse | 0.477 ms | 0.501 ms | 1 species — serial body, unchanged |
| group | 2.128 ms | 2.237 ms | unchanged |
| TCSPC | 2.193 ms | 2.262 ms | 1 species — **bit-identical** |
| TCSPC VV | 2.993 ms | **2.639 ms** | 1.13x |
| FRET | 18.135 ms | **6.179 ms** | **2.94x** |

And against the numpy path the graph replaced: TCSPC **5.8x**, VV **6.5x**,
FRET **4.4x** — the FRET row was 1.2x when the chain was built in (14)/(15),
and 1.5x after the amplitude threshold in (16).

The decay node alone, which is what actually changed: **106 µs -> 20.0 µs**
at 53 species. A single-species TCSPC decay is 4.2 µs, exactly what (16)
recorded for it, because it takes the same code it always did.

### The shape of this, worth remembering

Three sessions attacked a FRET fit and only the third one moved it:
(14)/(15) removed the crossing — **1.2x**, because the crossing was not the
cost; (16) stopped computing species that weigh nothing — 1.2x more; (17)
took the error estimate out of Python — and this one made the arithmetic that
remained wait less. Each step was worth measuring the node before writing
code, and (16)'s per-node table is what pointed at this one.

### Tests

`imp.bff` **1830 passed**, none failed — and again after tttrlib was rebuilt
underneath it. `chisurf test/fitting` **1003 passed** plus the same two
pre-existing failures (`test_fit_state`, `test_pcf_experiment`), likewise
twice. `tttrlib test/python/decayfit` **172 passed, 1 skipped**. The curve-parity tests that pin the graph against
`update_model()` at `rtol=1e-10` pass untouched, which is the check that
mattered: a 5e-16 shift is four orders inside them.

## 2026-09-01 (17) — the optimiser, and everything it needs, is bff's

Owner: *"lmdif etc must be also in bff."* Ticket `T-20260901-10`; the plan is
`okf/handover-error-estimate-2026-09-01.md`.

### What was wrong, in one line

A fit optimised entirely in C++ and then called the Python model **five more
times** to rebuild a Jacobian for its error bars -- 32% of a TCSPC fit. And
underneath that: **three** implementations of one bounded Levenberg-Marquardt
lived in this stack, and **two** of its finite-difference Jacobian, *with
different step rules*.

### The result

`fit.run()`, best-of-6, with Python `update_model()` calls counted inside it:

```
                 before                                after
fit        total  optimise errors rest  calls   total  optimise errors rest  calls
parse      0.50ms   60%      4%   36%     2     0.46ms   60%      4%   36%     2
tcspc      3.26ms   41%     32%   27%     7     2.25ms   59%      2%   39%     2
tcspc VV   4.51ms   44%     34%   22%     8     3.00ms   68%      1%   31%     2
FRET      24.3 ms   80%     12%    8%     9    18.6 ms   90%      0%   10%     2
```

**Two Python model evaluations per `run()`, whatever the model is**, and none
of them for the covariance -- which is what a parse fit always did and what
the ticket asked for. The clock follows: a TCSPC fit is 1.45x, VV 1.50x, FRET
1.31x, on top of the 4-6x the graph had already bought.

### Why MINPACK's own covariance had to be refused, and what replaced it

`lmdif` differences at `h_j = sqrt(epsfcn)*|x_j|` -- a step **relative to the
parameter**, floored only at exactly zero. A parameter that converges near
zero is differenced at a step near zero and its Jacobian column is round-off;
`covar` then reports an enormous variance for it. Nothing about that is wrong
in the optimiser, and *scipy's own* covariance for the same fit is no better.
Every real TCSPC model hits it: a scatter fraction of 1.5e-5 beside a
lifetime of 3.15.

The fix is not a better matrix in numpy but the **same differences taken over
the graph**, at chisurf's `approx_grad` step rule -- `eps * max(|x|, 1)`, an
**absolute floor**, which is exactly what the optimiser's rule lacks.
`Minimizer::compute_covariance` does that, and `epsfcn` was not touched: the
optimiser's step is tuned for convergence (1e-6 is worth 23 more converged
fits in 88) and the covariance's for resolution. Two jobs, two steps.

Pinned against numpy before it was wired, on all four fixtures:

```
parse    max rel on sigma 2.1e-14      tcspc VV  3.5e-07
tcspc    max rel on sigma 2.8e-07      FRET      1.4e-07   (same column dropped)
```

The FRET row is the one that matters most: `E_FRET` is a free parameter
nothing can move, so its column is exactly zero **by design**, and both paths
must *drop* it rather than invert a singular `J'J`. They do.

### Three implementations became one

| where | lines | what happened |
|---|---|---|
| `imp.bff/src/Minimizer.cpp` | 1041 | the one that survives |
| `chisurf/core/math/optimization/leastsqbound.py` | 794 | **deleted** |
| `chisurf/plugins/.../lltf/core/optimization/leastsqbound.py` | 365 | **deleted** |

`minimize`'s fallback is now `ResidualNode` + `bff.Minimizer` -- the same C++
optimiser driving a Python residual through a director -- so the algorithm is
bff's whatever the model is.

**The measurement that once rejected that fallback: I could not reproduce the
handover's "it is a wash".** Timing `minimize` alone, interleaved,
min-of-many, with the covariance switched off so both do the same work:

```
                scipy      director    ratio
parse           0.863 ms   0.961 ms    1.11x
tcspc          11.74  ms  12.47  ms    1.06x
```

So the director is **6-11% slower**, not equal. It was 1.29x / 1.22x before I
cached `ResidualNode`'s port lookups -- `get_input_port` was being called once
per parameter per residual evaluation, for an answer settled in the
constructor. The remaining few per cent is paid deliberately: it is a few per
cent of the path taken by models that build *no* graph, against two
implementations of one algorithm that disagree quietly. This stack has
already paid for that once -- the permuted covariance in (9) was found only
because two of them were finally run against each other.

The reference did not die with the implementation: `leastsqbound.py` is
frozen at **`test/minimizer/reference_leastsqbound.py`**, imported by nothing,
and the parity tests still assert the port is 1:1 against it. A skipped parity
test on a 1:1 port is worth nothing.

### The scope surprise: six other consumers

`covariance_matrix` is not only the error estimate. It is also the posterior
view, error propagation onto derived quantities, both sampler
preconditioners, and `Fit.grad`. Rather than edit six call sites, the choice
was made **inside `covariance_matrix` itself**: it tries
`curvature_over_the_graph` first and falls through to numpy for a model that
builds no graph. Measured, per call:

```
              update_model() calls        covariance   grad
              graph   numpy
parse           0       5                 2.1e-14      bit-identical
tcspc           0       6                 2.8e-07      7.5e-06 abs
tcspc VV        0       7                 3.5e-07      7.5e-06 abs
FRET            0       8                 2.2e-07      1.1e-05 abs
```

### Two things found on the way, neither of them mine

* **A Python `Node` held by a `Minimizer` can be collected, and it segfaults.**
  `Minimizer` holds `std::shared_ptr<Node>`, but a SWIG director keeps only a
  *weak* reference to its Python proxy, so a `ResidualNode` that falls out of
  scope while its minimiser is alive is a use-after-free. Reproduced on plain
  `run()`, so it predates this work; advertised as `T-20260901-12`. It is the
  same defect `IMP_SWIG_DIRECTOR` fixes for `MinimizerObserver`. Binding the
  node to `_` is enough to lose it, because the next tuple unpacking rebinds
  `_` -- which is how it was found.

* **A sampler is chaotic in the last digits of the error estimates that seed
  it.** `test_a_skewed_ratio_...` in chisurf compares a sampled quantile
  against a hand-picked factor of 1.5. The error estimates moved by **3e-8**
  when they started coming from `Minimizer` rather than numpy -- and the
  upper arm of the sampled interval moved by **7%**, from 0.1535 to 0.1648,
  taking the assertion from 1.51 to 1.48. The property it defends is intact
  (the symmetric width is still 48% wider than one arm and 43% narrower than
  the other); the factor was tuned to three digits against a chain that
  cannot honour three digits. Loosened to 1.4, with the measurement written
  into the test.

### Where the tests stand

* `imp.bff`: **1830 passed**, none failed, up from 1816 -- the new ones
  cover `compute_covariance`: agreement with numpy, the step rule *read off
  the objective it is handed* rather than inferred from the matrix, the zero
  column dropped, external coordinates, the ports restored, `p + 1`
  evaluations.
* `chisurf test/fitting`: **1003 passed** with the same two pre-existing
  failures (`test_fit_state`, `test_pcf_experiment`), plus
  `test_optimization_progress.py` (7) rerouted at bff's
  `minimizer_expected_evaluations` / `minimizer_reported_total`.
* `test_the_error_estimates_do_not_move` **was not testing anything** for
  TCSPC: both paths fell back to the same numpy `covariance_matrix`, so it
  compared numpy against numpy. It now has a sibling in each of the three
  graph-fit files that compares the matrix actually handed over against
  `covariance_matrix`'s -- verified by breaking the C++ step rule on purpose
  and watching them fail.
* `census_models.py` unchanged: 7 of 42, 3 of 14 under VV, no model builds a
  graph that disagrees with its own curve.

## 2026-09-01 (16) — the engine's own cost: a decay pays per species, and most of them were free of charge

Owner: *"bff and tttrlib are the engine, chisurf is the car. move more to
bff cpp to make the engine go fast."*

### Where the time actually is, measured per node

Profiling the graph a node at a time, rather than guessing, put the whole
question in one table (min-of-many, one `evaluate()` per row):

```
FRET objective = 195 us      donor   0.4 us    species 1
                             fret    1.0 us    species 97
                             dist    1.0 us
                             decay 192.3 us    species 97   <-- 98.6%
                             chi2    1.3 us
TCSPC objective =  5.2 us    decay   4.2 us    species 1
```

The producers this session added cost **2.4 us of 195**. Essentially all of a
FRET evaluation is `fconv_per_cs_ad`, and it is linear in the species count:
1 species is 4.2 us, 97 species is 192 us.

### The species count is not a physical quantity, it is a discretisation

A model whose spectrum comes from a *distribution* has as many species as
that distribution has bins, **whatever their weight**. A FRET decay over the
96-point `rda_axis` is 97 species. Asked how many of them carry weight, on a
Gaussian at 45 A with sigma 6:

| threshold (relative to the largest amplitude) | species | max abs. curve change / peak |
|---|---|---|
| `0` (exact) | 97 | 0 |
| `1e-14` | **53** | **1.0e-15** |
| `1e-12` | 44 | 5.2e-14 |
| `1e-10` | 37 | 2.0e-11 |
| `1e-6` | 26 | 2.1e-7 |

`TcspcDecay` gained `set_amplitude_threshold(relative)`, applied after
`absolute_amplitudes`/`normalize_amplitudes` because those describe what an
amplitude *means* and this describes which of them are worth computing.

**The node defaults to `0`** -- which drops only exact zeros and is therefore
bit-exact, and still fires, because a donor-only fraction of zero contributes
a whole block of them. ChiSurf asks for `1e-14`
(`minimizer.AMPLITUDE_THRESHOLD`). That split is deliberate and it is the
same lesson as (15): a *library* that silently drops terms is the failure
this module spent the morning fixing, so the application that knows what its
data are worth makes the choice, in the open, with the table above beside it.

`1e-14` is chosen because **the change it makes to the curve is smaller than
the change from summing the same species in a different order.** It is not a
trade of accuracy for speed; it is declining to compute terms whose
contribution double precision cannot represent. The looser rows are real
trades and are deliberately not taken.

### What it is worth, A/B'd in one session rather than across two

| | threshold off | threshold `1e-14` |
|---|---|---|
| FRET, whole `run()` | 22.40 ms | **18.30 ms** (1.22x) |
| FRET, decay node alone | 194 us | **106 us** (1.83x) |
| TCSPC, whole `run()` | 2.189 ms | 2.204 ms (unchanged) |

The TCSPC row is the one worth keeping. Across sessions it *looked* like this
change had made a plain lifetime fit 1.46x faster -- and it had not: a
single-species model has nothing to prune, and the apparent gain was the
machine being in a different state than it had been an hour earlier. The A/B
in one process says so plainly. **A before/after taken across two sessions on
this laptop measures the laptop** (the same lesson as 2026-09-01 (5), learned
again).

### The remaining cost is a serial recursion, and it is tttrlib's

At 53 species the decay node is 106 us. `fconv_per_cs_ad` is two loops per
species over 512 channels, each a *serial dependency chain*
(`fitcurr = (fitcurr + c1) * expcurr + c2`), so it is latency-bound rather
than throughput-bound: at ~4 cycles per dependent FMA, 53 x 512 x 2 x 4 is
~62 us of the 106 at 3.5 GHz. **The species are independent of one another**,
so interleaving several recursions would fill the pipeline -- but that
changes the order the per-species contributions are summed, and the header is
a byte-identical vendored copy of tttrlib's kept in step by a test. So it is
a tttrlib change vendored back, the same route `LatticeDiffusion.h` took, and
it is advertised as `T-20260901-12` rather than forked here.

Also removed: a heap allocation and a 512-double copy *per objective
evaluation* (the response was copied into a local before being normalised;
it is a member buffer now, walked twice instead of three times, same
arithmetic in the same order). **Honestly, it does not measure** -- 108 vs
106 us, inside the noise. It is kept because it removes an allocation from
the per-evaluation path, not because a benchmark asked for it.

### Tests

`imp.bff` **1830 passed** (+14: the threshold's own tests, plus one existing
composition test rewritten -- it compared the decay's spectrum against the
whole of what the anisotropy published, and at `l1 = 0` a VV channel appends
VH components scaled by exactly zero, which the node now drops. The property
it defends is better stated against the *contributing* part). `chisurf`
**1023 passed** plus the four known pre-existing failures.

One defect caught while writing it, before it could ship: pruning first
shrank `n_lifetimes_` in place, and nothing on the scalar-port path restores
that count -- so a component whose amplitude an optimiser walked *through*
zero would have disappeared permanently rather than for one curve. The kernel
now takes a separate `n_active_`, and `test_pruning_does_not_consume_the_ports`
walks an amplitude to zero and back.


## 2026-09-01 (15) — the photophysics joins the instrument: polarisation and FRET as nodes

Owner: *"improve all chisurf mdls by making bff the compute kernel."*

Closes the first half of `T-20260901-04`. Three things landed, and the first
one is a bug the other two would have inherited.

### A graph that builds is not the same as a graph that is the model

The census in (14) asked which models build a bff graph. It did not ask
whether the graph they build *is* them. Asking that found two:
`MaxEntLifetimeModel` and `LifetimeMixtureModel` both **override**
`lifetime_spectrum` -- one from a maximum-entropy inversion, the other from
other fits' spectra -- while `_lifetime_objective` builds the spectrum from
the `lifetimes` group's own amplitude and lifetime parameters. Measured on a
512-channel decay: the graph's curve and the model's were **793.8 counts
apart** at identical parameters, and nothing said so.

Refusing on the free parameters alone does not catch it, which is why it
survived. `MaxEntLifetimeModel`'s extra parameters -- its grid, its entropy
weight -- ship **fixed**, so it offers exactly the four free parameters the
plain model does, every port is placed, and the graph builds. The check is
now on the property itself: a model whose `lifetime_spectrum` is not
`LifetimeModel`'s does not get the plain multi-exponential graph.

`test/minimizer/census_models.py` is the census, checked in, and it now
compares each built graph's curve against `update_model()` and prints
**MISMATCH** rather than a pass. *Being representable is not the same as
being represented*, and only the second is allowed to run.

### The instrument was only half the statement

`TcspcDecay` holds what every TCSPC model shares -- reconvolve, add scatter,
scale to the data, add a background. What differs between models is entirely
*upstream* of it, and (14) left the `lifetime_spectrum` port as a mechanism
with no user. It has users now, in `SpectrumNode.h`/`.cpp`:

| node | what it produces |
|---|---|
| `LifetimeSpectrumNode` | the interleaved spectrum from scalar `a`/`t` ports -- the brick everything else starts from |
| `AnisotropySpectrum` | the polarised spectrum a VV, VH or VV/VH detector sees |
| `GaussianDistances` | a distance distribution as a weighted sum of generalised normals |
| `FretSpectrum` | a donor spectrum quenched over a distance distribution, with the donor-only fraction mixed in |

so a model is a chain and the instrument never learns what happened
upstream:

    [GaussianDistances] -> [FretSpectrum] -> [AnisotropySpectrum] ->
    TcspcDecay -> ChiSquared -> Minimizer

**Why this composes at all is the physics, not a convenience.** A
polarisation multiplies the decay by \(r(t)\), and the product of two sums
of exponentials is a sum of exponentials over the Cartesian product with the
harmonic mean of the time constants. A FRET quenching adds a transfer rate,
and rates add, so the quenched spectrum is the product of the donor's rates
and the transfer rates. Both are therefore *spectrum* transforms, which is
what lets each be one node with no opinion about the others.

### What it measures, and one number that is honestly small

```
TCSPC lifetime  14.5 -> 3.2 ms    4.58x
VV, one rotation 18.6 -> 4.4 ms   4.18x
FRET (96-point Gaussian)  29.8 -> 24.3 ms   1.23x
```

The FRET row is small and the reason is worth writing down rather than
explaining away: a 96-point distance distribution against a one-species
donor is a **97-species reconvolution**, and that is genuinely most of the
fit in either path. The crossing was already a small share of it. What the
graph buys there is not speed, it is that the model is now *in* the graph --
so a group of FRET fits, or a sampler over one, no longer has a Python leaf.

### The parity that licenses all of it

Every new node was checked against ChiSurf's own kernel before anything was
wired: `AnisotropySpectrum` against `calculcate_spectrum` over 160 random
(spectrum, rotation, g, l1, l2, polarisation) draws, `FretSpectrum` against
`distribution2rates` + `rates2lifetimes` over 30, `GaussianDistances`
against `combine_distributions` over 30 -- all to 1e-12 or better, entry for
entry **and in the same order**. Ordering matters here in a way it does not
elsewhere: a mixed polarisation channel is the *union* of two scaled spectra
rather than their sum, so its length and its layout are observable, and a
fit built on a differently-ordered one is wrong in a way no chi-square would
report.

Then the models: the graph's curve equals `update_model()` to ~1e-15
relative for VM/VV/VH/VV-VH at one and two rotations, with and without
amplitude normalisation, and for `FRETModel` and `GaussianModel` at
`xDOnly` 0 and 0.25 under every polarisation.

### Two things learned about *testing* a fit, not about fitting

**A degenerate fixture measures conditioning, not correctness.** The first
polarised parity test failed with the two optimisers a whole sigma apart --
because a rotational correlation time fitted to magic-angle data is a flat
direction, and a rotational *amplitude* is redundant against `r0` by
construction (its error estimate comes back `nan`). Same again for the
Gaussian weight of a single-component distribution. The fixtures now
generate their data from the model itself and fix what is redundant; the
objective is pinned separately, and to machine precision, by a test that
does no fitting at all.

**`E_FRET` is free and neither path can move it.** Its value comes from a
callable and its setter is ignored, so a write is a no-op in numpy *and* in
the graph -- where its port is deliberately wired to nothing. Both see a
column of zeros and agree. Refusing instead would refuse every FRET model,
since the parameter ships free.

### What still falls back, and why each one is its own kernel

Every remaining `FRETModel` subclass **is** a distance distribution:
`SingleDistanceModel` (a histogrammed discrete set), `WormLikeChainModel`,
`SawNuModel`, `IsingChainModel`, `FRETrateModel` (a rate spectrum given
outright), `MaxEntFRETModel`, `FRETStructure`, and `PDDEMModel` (which
replaces the lifetime spectrum entirely). Also refused: a kappa-squared
*spectrum* (`orientation_mode = "slow"`), which convolves the distances
before they become rates, and `bin_lifetime`, which coarse-grains the
spectrum afterwards. Outside the decay families, DEER, ICS, PCH/FIDA,
PDA2C/3C, MFD and the structure model are untouched. Advertised as
`T-20260901-06`.

### The rule this was converging on, now written down

The owner stated it in general on the same day: *"keep it all in bff and
tttrlib, and only things that get displayed should be moved to chisurf /
Python, so that the data stays local."* Written up as a cross-stack page --
`../chisurf/okf/architecture/compute-display-line.md`, linked from the
architecture index and summarised in `AGENTS.md` -- with a decision procedure
(*who consumes the value: a human, or another computation?*), what
legitimately stays in chisurf, and where it is enforced.

Measuring the distance from it found the largest remaining violation, and it
is **not** in the optimisation, which crosses zero times per iteration. It is
the **error estimate**: `fit.run()` timed and Python `update_model()` calls
counted gives 0.50 ms / 2 calls for a parse fit but 3.26 ms / **7 calls** for
a TCSPC one, five of them rebuilding a Jacobian in numpy *after* the
optimiser stopped -- 32% of the fit. The cause is that `lmdif` differences at
a step **relative to the parameter**, so a scatter fraction of 1.5e-5 beside a
lifetime of 3.15 gets a round-off column, and `_forward_differences_resolved`
correctly refuses the C++ covariance. The fix is to difference the *graph* at
a step that resolves (chisurf's own `approx_grad` already floors the step at
1.0, which is the specification), not to compute a better matrix in numpy.
Handed over as **`T-20260901-10`** with
[`okf/handover-error-estimate-2026-09-01.md`](handover-error-estimate-2026-09-01.md);
`T-20260901-11` covers the model curve being recomputed for display rather
than read off the graph.

### Tests

`imp.bff` **1816 passed** (+17: `test/spectrum/`). `chisurf test/fitting`
**1000 passed** plus the two known pre-existing failures (`test_fit_state`,
`test_pcf_experiment`, both failing on `HEAD`); `test/architecture` gained
two provenance tests and keeps its own two pre-existing failures
(`test_guarded_imports`, `test_optional_backends`, verified against `HEAD`).
`test/minimizer/bench_fit.py` grew a VV table and a FRET table.


## 2026-09-01 (14) — which models actually run in bff, and the port that composes the rest

Owner: *"make sure that the fits, and models, live and work in bff and
tttrlib. chisurf is just the glue and interface."*

### The census, because the directive was not actionable without one

52 model classes live under `chisurf.core.models`. Asking `graph_objective`
which of them build a bff graph, rather than assuming, moved the picture a
long way:

**Already entirely in bff** -- and this was the surprise -- every
*equation-based* family, because they are `Expression` nodes and have been
since that landed. Verified by construction, not by reading:

```
parse.ParseModel                   graph=YES   free=[a, t, b]
fcs.ParseFCSModel                  graph=YES   free=[b, N, td, s]
pcf.ParsePCFModel                  graph=YES   free=[A, sigma, Mode]
stopped_flow.ParseStoppedFlowModel graph=YES   free=[a1, k1, c]
```

So **FCS, PCF and stopped-flow are not "waiting for a port"** -- they are
compiled equations over vectors in C++ already. Add the plain TCSPC
`LifetimeModel` (2026-09-01 (12)) and any `FitGroup` over any of these
(2026-09-01 (11)), and the fitted-in-bff set is much larger than the earlier
entries implied. `T-20260901-04` was written on the assumption that "FCS and
PDA still fall back"; half of that is wrong and the ticket is corrected.

**Still falling back**: the `LifetimeModel` *subclasses* (FRET, Gaussian,
SawNu, WormLikeChain, IsingChain, SingleDistance, PDDEM, the mixture and
max-ent models), plus DEER, ICS, PCH/FIDA, PDA2C/PDA3C, MFD and the structure
model.

### The one thing that unblocks the whole subclass family

Reading them together makes the shape obvious: **every one of those TCSPC
subclasses has the same instrument model** -- reconvolve, add scatter, scale
to the data, add a background -- and differs *only* in how the (amplitude,
lifetime) pairs are arrived at. A FRET model derives them from a distance, a
distribution model from a distance distribution, a mixture from two spectra.

`TcspcDecay` took its pairs only from scalar `a0`/`t0`/... ports, so each
subclass would have needed the deriving done in the caller, once per
iteration -- the arrangement that measured as a *regression* in (9). It now
also accepts the interleaved spectrum on an input **port**
(`set_spectrum_from_port(true)`, key `lifetime_spectrum`), so the graph
composes:

    <whatever computes a spectrum> -> TcspcDecay -> ChiSquared -> Minimizer

The node keeps its opinion about the instrument and holds none about the
photophysics. `absolute_amplitudes` and `normalize_amplitudes` still apply,
because they describe what the *model* means by an amplitude rather than
where the number came from; an odd-length spectrum is refused rather than
rounded down, which would silently drop the last amplitude's lifetime; and
the port does not sanitise, so a NaN lifetime reaches `ChiSquared` and makes
the misfit infinite instead of being floored to `tiny` and reading as a good
fit.

The test that matters drives it from a real upstream node (an `Expression`
standing in for a photophysics model) through a linked vector port, and gets
the same curve as the scalar path to 0.0 -- what is being checked is the
wiring, not the arithmetic.

### What is left, honestly

The port is the *mechanism*; no subclass uses it yet. Each still needs its
spectrum producer to become a node -- and for the FRET family much of that
physics is already on this side (`LifetimeSpectrum`,
`lifetime_spectrum_from_rates`, `convolve_distance_with_k2_ratio`,
`RotamerFret`), so those are wiring jobs rather than ports. The anisotropy
case is different and better: VV/VH is `dfa_vv_vh_convolved` in tttrlib
already.

### Tests

`imp.bff` **1799 passed** (+5, `SpectrumPortTests`). `chisurf test/fitting`
**963 passed** plus the two known pre-existing failures. The benchmark is
unmoved -- TCSPC 13.76 -> 3.17 ms, **4.34x** -- which is the point: the
scalar path had to keep costing exactly what it did.

### Also this session

Verified the previous entry's open question. The three setup suites that
exercise the module moved out of the view in (13) -- `test_fcs_setup`,
`test_setup_user_migration`, `test_setup_prerequisites` -- **34 passed** once
`../mmfdb/src` is on the path. The failures reported against them earlier
were a missing sibling checkout on `sys.path`, not the move. Worth writing
down because it cost a turn: **`mmfdb` is a sibling checkout, not installed**,
so `PYTHONPATH=.:../mmfdb/src` is required for `test/fio`.

## 2026-09-01 (13) — the factor graph stops being two implementations, and one of them was wrong

Owner: *"bff must become the core engine of chisurf. in chisurf there is a lot
of stuff that can probably be simplified once bff fully landed."*

`chisurf.core.fitting.factorgraph.FactorGraph` now answers every structural
query through `IMP.bff.FactorGraph`. Moralisation, greedy elimination
(min_fill and min_degree), the maximal cliques, the junction tree, the
separators, the sampling blocks, the treewidth and the relevance queries
(`affected_factors`, `affected_fits`, `unexplained_variables`) were **two
implementations of the same algorithms** -- one in C++ here since the PRD-68
port, one in ~200 lines of Python there -- and are now one. The engine is a
lazily built mirror on a new `FactorGraph.engine` property, dropped by
`invalidate()`.

**What stays in chisurf, and why it is not laziness.** *What a variable is* --
which parameter, at which position of which model's flat vector, belonging to
which local fit -- and *what a factor is*. Discovery is the application's
knowledge: it walks models, resolves parameter links, decides which prior is
more than its bounds. Triangulating a graph is not. Also staying: `describe()`,
because it renders parameter *names* rather than keys, and `markov_graph()`,
because it returns a `chisurf.core.graph.Graph` the posterior views draw.

### The bug this found

**bff's `is_complete()` double-counted, and a star hit it exactly.** The moral
adjacency is stored symmetrically, so summing the neighbour-set sizes gives
`2|E|`; it was compared against `n(n-1)/2`. That is true whenever
`|E| = n(n-1)/4` -- and **three datasets around one shared parameter is four
variables with three edges**, `2*3 == 4*3/2`. The graph was then taken for a
clique and the completeness shortcut fired: one clique over everything, a
treewidth of **3** where the answer is **1**, and a sampler offered one
4-dimensional block where there are three 2-dimensional ones.

It is the shape a global fit *has*. It survived because the existing bff star
fixture has ten variables and sixteen edges, which does not hit the
coincidence, and because nobody was comparing the two implementations -- which
is precisely what the duplication cost. Fixed in `FactorGraph::is_complete()`;
regression tests on both sides, the minimal four-variable star and a genuinely
complete four-variable graph so the guard cannot overshoot.

That is the argument for this whole direction, stated as a measurement: two
implementations of one algorithm do not average out, they disagree, and the
disagreement is silent until something forces them to answer the same
question.

### Two canonical orders, deliberately kept

A clique is a set, so its member order is a convention. bff orders members by
**flat-vector position**, which is what a C++ caller lining a clique up
against a parameter array wants; chisurf's documented contract is **sorted by
key**, which is what a caller using a clique as a dict key or a graph node
relies on. chisurf sorts on receipt. Sampling *blocks* keep the engine's
order, because a block is consumed as a slice of the parameter vector.

### Tests

`chisurf test/fitting` **963 passed** plus the two known pre-existing
failures; `imp.bff` **1794 passed**. Two chisurf tests were rewritten rather
than repaired: they counted copies of the Python moral graph, which no longer
exists. The property they defended -- computed once per graph, dropped by
`invalidate()` -- is now stated against the engine's identity, which is the
observable that replaced it.

### Left open

The line count barely moved (-14) because the delegation kept the docstrings;
the substance is that ~200 lines of algorithm became zero. What is genuinely
next, and is advertised as **T-20260901-05**:

- `chisurf/core/fitting/graphview.py` (880 lines) is *layout* over this graph.
  It should shrink now that the structure comes from one place, and the parts
  that are graph algorithms rather than presentation belong beside the engine.
- `chisurf/core/graph/` (1966 lines: a graph type, algorithms, layout,
  GraphML) is a general-purpose graph library living in an application. Some
  of it is already in bff.
- Sessions are **already** on bff (`bff.Session.load` in
  `chisurf/core/project/project.py`), so that half of the owner's example is
  landed; what is left there is the `.csp` archive plumbing around it.

## 2026-09-01 (12) — the decay on the graph: bff builds the network, tttrlib computes the curve

Owner: *"all computation should be offloaded to bff or tttrlib in chisurf.
bff for network building, tttrlib for some scoring, bff tttrlib must interact
without chisurf in the middle, minimize friction."*

Closes `T-20260901-03`. `TcspcDecay.h`/`.cpp`: **a TCSPC decay as a node.**
The model curve of a multi-exponential decay through a real instrument --
reconvolved with a measured response, shifted against it, scaled to the data,
plus scatter and a background -- with the amplitudes, lifetimes and nuisance
scalars on ports. `TcspcDecay -> ChiSquared -> Minimizer` is then a whole
lifetime fit in one C++ graph, exactly what `Expression` gives a parse model.

**The arithmetic is not in bff, and that is the point.** The periodic
reconvolution is tttrlib's `fconv_per_cs_ad<double>` and the timeshift is its
`shift_lamp_ad<double>`, both from the vendored byte-identical copy of
`DecayConvolution.h`. Same decision as the expression engine, for the same
reason: a decay is a *curve*, `AGENTS.md` puts curves in tttrlib, and a second
convolution is how two libraries end up disagreeing about what a lifetime is.
What bff contributes is the network -- ports, ordering, invalidation -- and
the composition of the instrument model around those kernels. chisurf is no
longer between them: it hands over the response and the data once and reads
the answer back once.

**One small change in tttrlib made the whole thing header-only.**
`shift_lamp` had its body in `DecayConvolution.cpp`, which a header-only
consumer cannot reach; it now lives in the header as `shift_lamp_ad<T>` and
the exported function calls it. `fconv_per_cs_ad` was already that shape and
is the precedent. So bff needs one vendored header and links nothing --
which is what "minimize friction" bought here, against the alternative of
vendoring a 983-line .cpp and its `Registry.h`/`Verbose.h`/`info.h` chain.

**tttrlib's `shift_lamp` *is* chisurf's `shift_array`, sign-flipped.** Found
while looking for a third implementation to avoid writing.
`shift_lamp(out, v, -s)` equals `shift_array(v, s)` to 0.0 for every shift
tested -- integer, fractional, both signs. The two index in opposite
directions *and* interpolate toward opposite neighbours, and the flips
cancel. They part company at exactly `s = 0`, where `shift_lamp` zeroes the
last sample (`out_right = tsint + 1` is 1 when `tsint` is 0) and
`shift_array` returns the input untouched; chisurf short-circuits `shift ==
0.0` before calling, so the node does too and is exact. **chisurf's numpy
`shift_array` is a duplicate of a tttrlib kernel** and is now a candidate for
deletion rather than a maintained twin.

**chisurf's `rescale_w_bg` is *not* tttrlib's, despite the name.** tttrlib
guards `decay[i] > 0` and adds `1e-12` to the squared weight; chisurf guards
`e > 0` *and* a finite weight, adds no epsilon, and returns the factor
instead of rescaling in place. The node reproduces chisurf's, deliberately
and with the divergence written down: unifying them moves fitted amplitudes,
which is a decision about the fit and not about where code lives.

**What the node refuses**, each because it is a term the node does not carry
and a term left in Python would put the per-iteration crossing back: pile-up
(Coates), a DNL linearisation table, a measured background curve, VV/VH
polarisation (which rotates the spectrum before convolving), a non-periodic
convolution mode, the convolution switched off, a response whose length
differs from the data's, and any free parameter the node cannot place --
which is what actually catches a `LifetimeModel` subclass that adds a term.

### A wrong error bar, found by the new model and true of the old path too

The first version reported a standard deviation on the scatter fraction
**140x** the finite-difference one. It is not the port's doing:
`lmdif` differences forward at `h_j = sqrt(epsfcn)*|x_j|`, a step *relative
to the parameter*, and this fit's free vector spans `1.5e-5` (a scatter
fraction) to `3.15` (a lifetime) -- so the scatter is probed at ~1e-8, its
Jacobian column is round-off, and `covar` reports nonsense for it.
**scipy's own `leastsqbound` covariance for the same fit is no better**:
0.094 for the same parameter and exactly *zero* for two others. Both
implementations fail together, which is how it is known to be the estimator.

`minimize` now refuses to stash the optimiser's covariance when any
`|x_j| < sqrt(epsfcn) * max|x|` -- with that step, such a parameter is being
differenced below the resolution the objective is computed to. The whole
matrix goes, not one column: the columns are not independent and dropping one
silently changes what the others mean. Refusing costs a Jacobian rebuild
(about a third of a *parse* fit, much less of a decay fit) and gives exactly
the numpy path's answer, which is what `test_the_error_estimates_do_not_move`
now pins at 5e-3. A well-scaled parse fit still keeps the matrix, and there
is a test for that too -- the guard must not throw away the saving it was
built to allow.

### Measured

`test/minimizer/bench_fit.py`, now three tables (interleaved, min-of-many):

```
parse fit, 512 points, 3 free      1.40 / 1.33 / 0.53 ms    2.64x
FitGroup, 4 members, 5 free        3.84 -> 2.34 ms          1.64x
TCSPC lifetime, 512 ch, 4 free    13.92 -> 3.20 ms          4.35x
```

The decay gains most because it has the most to gain: its per-iteration
Python was a reconvolution, a shift, a rescale and four numpy copies, not a
three-term expression.

### Tests

`imp.bff` **1792 passed** (+21: `test/decay/test_tcspc_decay.py`, plus the
drift guard on the vendored header). The ones that matter pin the node to
something that is *not* the node: the convolution against tttrlib's exported
`fconv_per_cs`, the timeshift against chisurf's `shift_array`, and the whole
curve against a numpy transcription of `LifetimeModel.update_model`.
`chisurf test/fitting` **962 passed** (+16,
`test/fitting/test_graph_fit_tcspc.py`) plus the two known pre-existing
failures. tttrlib's changed translation unit compiles clean; the refactor is
a body move and the bff tests pin its behaviour transitively.

### Left open

FCS, PDA and the `LifetimeModel` *subclasses* (FRET, anisotropy, PDDEM) still
fall back. Each adds terms to the same decay, so each is a node beside this
one rather than a rewrite of it -- the anisotropy one is the obvious next,
since VV/VH is `dfa_vv_vh_convolved` in tttrlib already.

## 2026-09-01 (11) — the group on the graph: 1.01x becomes 1.74x, and the link that never leaves C++

Closes `T-20260901-02` (handover: `okf/handover-fit-graph-2026-09-01.md`).
All the code is in chisurf; the bff half shipped in (10).

**The gap this closes.** A plain `Fit` on a parse model has run entirely in
C++ since (9) and is 2.85x. A `FitGroup` was **1.01x** -- and a `FitGroup` is
what the GUI builds, `global_optimize_local_first` shipping `false`, so a
group is exactly one optimisation over `GlobalFitModel` and
`graph_objective` refused it. The speed-up reached almost nobody.

**The mapping is exact, which is the whole reason this is a graph rather than
an approximation of one.** `GlobalFitModel.weighted_residuals` is its
members' residuals concatenated; `JointChiSquared`'s residual is its
members' blocks end to end, in the order they were added.
`GlobalFitModel.parameters` is each member's free parameters followed by the
group's globals; that is the order of the minimiser's ports. Nothing had to
be reinterpreted.

**The sharing never leaves C++, and that took the most care.** A member
parameter that is linked is *not free*, so it is not in
`GlobalFitModel.parameters` and gets no port of the optimiser's. Its equation
port is instead `set_link`ed to the port of whatever it follows -- the same
relation ChiSurf's `Parameter.link` is, one level down. Two details that are
not obvious from the code:

- The chain is **walked**, not stepped once. A master may itself be linked,
  and a master that is *fixed* stops nothing: the port chain resolves
  through it exactly as `Parameter.value` does. Stopping at the first master
  would freeze a follower at its start value whenever the link went through
  a fixed parameter.
- A **global parameter appears in no equation.** Members reach it only by
  linking, so it gets a standalone `Port` that the optimiser writes and the
  followers follow. There is nothing else for it to be.

**A group is refused whole.** A member the graph cannot represent cannot be
left in Python and joined to the others: `JointChiSquared` has one objective,
and half a group crossing the boundary once per iteration measures like the
director path -- which (9) recorded as a *regression*, not a speed-up. (It
is 1.06-1.11x, not 1.16x; see the revision on (9) and entry (17). It became
the fallback in (17) all the same.) The
lesson from (9) carried straight over: the crossing is the cost, and only
removing it entirely helps.

**Masks refuse the group, and that is a statement about the Python path, not
the C++.** `GlobalFitModel.weighted_residuals` concatenates its members'
residuals *unmasked* -- no member's mask is consulted -- and `_apply_fit_mask`
then applies the group's own mask, which is the *selected* member's, only when
its window happens to be as long as the whole concatenation. That is not an
objective worth reproducing in C++, so a mask anywhere in the group sends it
back to scipy.

**Three things fell out of the refactor.**

- `_member_objective` builds one dataset's half -- equation, ports, axis,
  data, window, noise model -- and deliberately does *not* decide what drives
  those ports. Whether a variable is the optimiser's, a follower, or a
  constant is a question about the whole fit. The single-fit and group paths
  share it.
- **A bug in the single-fit path, fixed by the same code.** A parse parameter
  linked to another *in the same model* was not free, so the old builder
  treated it as a constant at its start value while the numpy path had it
  follow its master. It now links.
- An equation variable that is `redundant` or callable-driven refuses the
  graph. Such a value is *derived* from other parameters, so it is neither a
  constant nor something the optimiser writes, and freezing it would be
  silent.

**A segfault in `Node::update()`, found by that third item.** Wiring the
single-fit path's links made the intra-model case reachable, and it crashed
the interpreter -- not an exception, a stack overflow. `Node::update()` walked
each linked input to its source *node* and updated it if invalid; for an input
linked to another port of the **same** node, that node is `this`, so it
re-entered itself forever. `inputs_valid()` had carried the guard
(`output_node.get() == this` -> `continue`) since it was written; `update()`
never got it. One line in `src/Node.cpp`, and a regression test in
`test/portnode/test_port_node.py` -- the crash is worth a test of its own
because two ports of one node linked together is *legal* and already had one
(`test_port_link_same_node_allowed`), which passed because it never updated
the node.

**Measured** (`test/minimizer/bench_fit.py`, which now carries a `FitGroup`
table beside the single-fit one; interleaved and min-of-many, as (5) taught):
four members of 512 points sharing one lifetime, five free parameters.

```
scipy     3.818 ms
graph     2.196 ms      <- 1.74x   (was 1.01x)
```

The single-fit table is unchanged at 1.30 / 1.31 / 0.46 ms, 2.82x.

**Why 1.74x and not 2.9x, which is not a disappointment.** Of the graph's
2.20 ms, 0.26 ms builds the graph and 1.13 ms is the LM loop -- so the
optimisation itself is about 2.2x. The remaining ~0.8 ms is the rest of
`FitGroup.run`: every member updated, the error estimates, the result
snapshot. That work is identical on both paths and dilutes the whole-run
ratio. A group also has more free parameters than a member, so it takes more
iterations of a loop that is already the fast one.

**Error estimates are unchanged, deliberately.** `FitGroup.model` is the
*selected member's* model, so `update_error_estimates` has always described
that member. The optimiser's covariance is stashed against the identity of
the parameters it was computed for, so a multi-member group's stash simply
does not match and the finite-difference path runs as before; a one-member
group's does match, and there the group objective *is* the member's.
`FitGroup.run` now drops a stale `_cpp_covariance` at the start, as
`Fit.run` already did.

**Tests.** `chisurf/test/fitting/test_graph_fit.py` 14 -> 28. The one that
earns its place is `test_the_group_is_not_two_separate_fits`: two datasets
with *different* true lifetimes, so the shared parameter must land strictly
between the two each recovers alone. A group that quietly optimised its
members one at a time would pass every other test in the file. It mirrors the
C++-level claim in `test/minimizer/test_joint.py`.
`chisurf test/fitting`: **946 passed**, plus the two known pre-existing
failures (`test_fit_state`, `test_pcf_experiment`, both failing on `HEAD`).
`imp.bff`: 285 passed (the self-link regression is the new one).

Still open: **T-20260901-03** -- `graph_objective` only builds for a model
whose curve is one compiled equation, so every TCSPC/FCS/PDA model still
falls back to scipy and is exactly as fast as before.

## 2026-09-01 (10) — the grouping in bff, and the sanitiser that made a broken model look perfect

Owner: *"bff also should reflect the grouping, do not necessarily call fit
group."*

`JointChiSquared.h`/`.cpp`: **one misfit over several datasets.** Its residual
is the members' residuals laid end to end and its chi-square is their sum,
which is what minimising a joint fit means -- one Levenberg-Marquardt step
moves the shared parameters using every dataset's curvature at once, rather
than each dataset in turn hoping they agree.

**Half of a joint fit was already here.** Sharing a parameter between datasets
is `Port::set_link` -- the members' ports follow one master, cycles refused --
so the coupling that makes a joint fit *joint* has been on this side of the
boundary since the Port runtime landed. Only the joint objective was missing,
and that is all this class is.

**Not named after ChiSurf's `FitGroup`, deliberately.** A `FitGroup` is a
container that also owns a selection, a result history, plots and a run
policy. This is the arithmetic underneath it, so the name says what the object
*is*. It is also not restricted to `ChiSquared` members: anything presenting a
residual vector on an output port qualifies, including a Python `Node`
director, so a group may mix representable and unrepresentable members and
still take one step.

The members are reached the ordinary way -- each member's residual port is
*linked* to one of this node's inputs -- so `Node::update()` walks the whole
tree from one call and nothing new had to learn how deep a member's graph
goes. Members keep their own data, window, mask and noise model, which is the
normal case rather than an awkward one.

### The defect this found, which is the important part of the entry

Writing a test for "a NaN in one member makes the group infinite" found that
**it did not** -- and then that the same was true of a *single* fit. chinet's
ports floor a NaN to `np.finfo(float).tiny` and clamp infinities, so a stored
value is always JSON-serialisable. For a document that is right. For the
numeric transport of a fit it is exactly backwards:

| | numpy path | graph path, before |
|---|---|---|
| model curve goes NaN | residuals NaN | curve floored to ~0 |
| chi-square | NaN | **14.0** -- finite and plausible |
| MINPACK | rejects the step | **accepts it** |

A model that blows up therefore looked like a model that is zero everywhere,
which for data near zero is a *good* fit -- so the optimiser was attracted to
the parameters that break the model, silently, and the graph path disagreed
with the numpy path in the dangerous direction. This had been true of every
`Expression -> ChiSquared` graph since entry (9), and no test would have found
it except one that deliberately broke a model.

Fixed with `Port::set_sanitize()`, **default unchanged**: an ordinary port
still behaves as chinet's, and only the fit's own transport -- an
`Expression`'s curve, a `ChiSquared`'s model input and residual output, a
`JointChiSquared`'s blocks -- opts out. `ChiSquared` has to clear it in an
`update()` override rather than in `evaluate()`, because `Node::update()`
writes the model input *before* `evaluate()` runs and the NaN would already
have been floored.

Two smaller things the same test found:

* `ChiSquared::set_data_arrays` and `set_mask_array` -- the **numpy** doors,
  the ones a per-iteration caller is told to use -- did not `set_valid(false)`
  where their `std::vector` twins did. Data set through them left the node
  claiming to be valid, so the next `update()` skipped `evaluate()` and served
  residuals against the *previous* data.
* `Node::update()` copied every linked input's value per call. It now reads
  the reference. A joint objective puts one full residual vector per member
  through that line on every iteration of a fit.

`std::vector<std::shared_ptr<Node> >` is not a template this module may name,
and wrapping it anyway leaks (SWIG finds no destructor), so the members are
reached by `get_member(i)` and `get_member_names()`.

`test/minimizer/test_joint.py`: 13 tests. The one that matters is
`test_the_group_is_not_two_separate_fits` -- two datasets with *different*
true lifetimes, where a shared parameter must land between the two separate
answers. A group that quietly optimised its members one at a time would land
on one of them and pass every other test in the file.

imp.bff `test/minimizer`, `test/chi2`, `test/expression`, `test/sampler`,
`test/portnode`, `test/factorgraph`, `test/session`: **283 passed.**

**Handed over:** [`okf/handover-fit-graph-2026-09-01.md`](handover-fit-graph-2026-09-01.md)
writes the remaining work up to be picked up cold -- the exact ChiSurf-to-bff
mapping, the eight traps that cost this session time, the build and test
commands, and `test/minimizer/bench_fit.py` to run before and after.

**Still open:** chisurf's `GlobalFitModel` is not yet wired to this, so a
`FitGroup` remains the slow path -- ticket `T-20260901-02`, now unblocked and
with a much shorter description, because the mapping is exact:
`GlobalFitModel.weighted_residuals` is a concatenation of its members' and
`GlobalFitModel.parameters` is a concatenation of theirs plus the global
ones, which is precisely this node's shape.

## 2026-09-01 (9) — the whole fit in C++: 2.85x, and three measurements that changed the design

Owner: *"chisurf must use imp.bff as compute engine... minimize swig boundary
crossing. Keep data on bff, only move data that is displayed. also port
minimizer, optimizer that are used in chisurf to bff."* Then: *"super
important is speed."*

`Minimizer.h`/`Minimizer.cpp`: MINPACK's `lmdif` -- `enorm`, `fdjac2`,
`qrfac`, `qrsolv`, `lmpar`, `covar` -- transcribed line by line, plus
chisurf's `leastsqbound` bounds transform, over free-parameter `Port`s and a
`Node` objective. Shaped like `Sampler`, which is the stochastic half of the
same idea.

**A plain `chisurf` parse fit is 2.85x**: 1.35 ms -> 0.47 ms at 512 points and
three free parameters, to the same answer to six decimals. Interleaved and
min-of-many, because load on this box drifts by more than most of the effects
below.

### Three measurements, each of which changed what got built

**1. Replacing the optimiser alone is worth nothing (1.02x).** Entry (6)
predicted it and it held. MINPACK's own arithmetic is ~4% of a fit; the
callback is 57%. The port is a **precondition**, not a speed-up: an objective
in C++ under an optimiser in Python still re-enters the interpreter every
iteration.

**2. The obvious fallback is a REGRESSION, so it was thrown away.** The
natural shape for a model bff cannot represent is `ResidualNode` -- the C++
optimiser driving the Python residual through a `Node` director. Measured:
**1.54 ms against scipy's 1.33 ms.** Wrapping a Python callback in a C++ loop
buys nothing over wrapping it in a Fortran one and costs the director
dispatch on top. A refused graph therefore falls back to `leastsqbound`.
(The director is 1.0 ms of that 1.54 before a separate fix: it was handing
back its residual through `set_value_vector`, i.e. converting 512 doubles
into a Python list per iteration. `Port::set_values_array` -- a numpy
typemap -- took it from 2.54 ms to 1.54. Still not enough to justify it.)

> **Revised twice; read the second revision.**
>
> *(15) said the director is "no longer a regression -- a wash or marginally
> faster",* from three interleaved min-of-many runs of whole `fit.run()`
> calls: scipy 1.307 / 1.293 / 1.299 ms against director 1.301 / 1.274 /
> 1.289. **That comparison was confounded and the conclusion was too
> strong.** A whole `run()` includes the error estimate, and by then a scipy
> run's error estimate had already moved into C++ while a director run's had
> not -- so the two rows were not doing the same work.
>
> **(17) re-measured `minimize` alone, with the covariance switched off so
> both do the same work: the director is 1.11x scipy on the parse fit and
> 1.06x on a TCSPC decay.** Still slower, by a few per cent rather than by a
> fifth. Half of what remained was `ResidualNode` calling `get_input_port`
> once per parameter per residual evaluation, for an answer settled in its
> constructor; caching that took it from 1.29x / 1.22x to the figures above.
>
> So the original conclusion -- *only removing the crossing helps* -- stands,
> and is still why the graph exists. What changed is the *size* of the
> penalty, from one that made the director an obviously wrong fallback to
> one worth paying: (17) deleted ChiSurf's second optimiser and a plugin's
> third, and a few per cent on the models that build no graph is a cheap
> price for having one implementation of the algorithm rather than three.
> See [`handover-error-estimate-2026-09-01.md`](handover-error-estimate-2026-09-01.md)
> and entry (17).

**3. The graph is the only thing that helps.** `Expression -> ChiSquared ->
Minimizer`: parameters are ports the optimiser writes in C++, the curve is
computed in C++, the data live in the node, and **nothing crosses per
iteration**. 0.16 ms for the optimisation alone against 1.33 ms for the whole
scipy fit.

### The blocker was avoided rather than solved, and that is the design

Driving the model's *own* ports from C++ is the obvious wiring and is wrong:
`chisurf.core.parameter.Parameter.value` keeps a `_frozen_value` cache inside
`frozen_structure()`, so a port written from C++ is **not seen** by the next
Python read -- the model would report its pre-fit values silently. So the
graph is **private to the fit**: built from the equation, the data and the
current parameter values, run, and then the answer written back through the
ordinary setters with one `update_model()`. The only place the numbers change
language is once per fit, which is also exactly "keep the data on bff and
move only what is displayed".

`graph_objective` refuses -- and the fit falls back -- on **semantics**, never
on C++ capability: a smooth prior (it adds rows the graph does not produce),
an equation the engine will not compile, a free parameter the equation does
not carry, a duplicate variable name. Bounds are *not* a refusal: they
synthesise a uniform prior that contributes no rows, and treating that as a
prior would have refused every bounded fit.

### Error estimation: a third of a fit, deleted

`update_error_estimates` rebuilt the Jacobian by finite differences -- `p + 1`
model evaluations, **32% of every fit** -- for a matrix `lmdif` already has in
its final `R`. It now takes the optimiser's, when there is one.

Entry (6) called this out as needing to be done *deliberately* rather than as
part of a speed pass, so it was: agreement with the finite-difference matrix
was measured **before** wiring, at **5e-4 relative** on the parameter standard
deviations across three models, against the four significant figures the
parameter table displays. It is offered rather than imposed -- stashed only on
convergence and only at full rank (a zero on the diagonal means a parameter
the data do not constrain, which chisurf reports by *dropping* the column,
and that is not what MINPACK's rank-deficient `covar` produces) -- and the
reader checks the stash against the *identity* of the current free
parameters, because a re-parse rebuilds them and a covariance for the
previous set would line up by length and mean nothing.

### What the parity bar could be, measured

Unbounded, the transform is the identity and the *whole trajectory* is
reproduced: every evaluation point, the same evaluation count, `info` equal,
x to 1e-12. Bounded it cannot be -- `leastsqbound` maps coordinates with
numpy's `sin`/`arcsin` and this maps them with libm's, `i2e(e2i(x))`
round-trips to **2.2e-16**, and twenty-odd LM iterations amplify one ulp to
~1e-7, inside `xtol = 1.49012e-8` *relative*. scipy's own answer moves by the
same amount when the bounds are merely respelled.

A fixture lesson worth keeping: a sum of exponentials started at (1, 1, 1, 1)
sits on its own permutation symmetry, where two pairs of Jacobian columns are
*identical* and which of two equal norms the pivoted QR takes is decided by
the last bit. Both optimisers find the same minimum there and label its terms
differently. The first version of these tests was asking a question that has
no answer.

### Two defects found on the way, one in each repository

* `chisurf`'s `leastsqbound(full_output=1)` returned a **permuted
  covariance**. MINPACK's `ipvt` is 1-based in Fortran, so the code did
  `ipvt - 1`; SciPy made it 0-based (1.18 `_minpack_py.py:488` builds its own
  `cov_x` from a bare `perm = ipvt`), and subtracting one produces `-1`,
  which numpy reads as *the last element* -- so the un-pivot was not a
  permutation at all. Against `scipy.optimize.leastsq`'s own `cov_x`: every
  element wrong, by up to **9x**. Latent only because error estimation took
  its own Jacobian -- and entry (6) had named reusing this matrix as the way
  to remove 32% of a fit, which is what this entry then did. Fixed by
  *detecting* the convention, so both keep working.
* `chisurf`'s `_internal2external_grad` tests `is None` while its transform
  lambdas test `_is_unbounded`, so `(-inf, inf)` yields `(inf - -inf) *
  cos(v) / 2` = inf. The two tests agree in the C++.

### Smaller things that are now true

* `Port::get_values_ref()` -- the values with no copy, C++-only. `Expression`
  reading its axis, `ChiSquared` reading the model curve and `Minimizer`
  reading the residuals each copied a full-length vector per iteration.
  Worth ~4%; kept because it is strictly less work, not because it showed up.
* `Expression::evaluate()` no longer copies every operand into a
  `vector<vector<double>>` per call, and allocates nothing in the loop.
* `ChiSquared` gained a `residuals` output port, written only when the node
  has one -- a `Sampler` wants the scalar and would otherwise pay for a copy
  of the whole residual vector per move.
* The observer is `IMP_SWIG_DIRECTOR`, not a bare `%feature("director")`,
  **and that was a segfault first**: C++ holds it across a whole run while
  SWIG's director keeps only a weak pointer to the Python proxy, so
  `m.set_observer(MyObserver())` -- how anyone would write it -- lost its
  proxy to the collector. Cancellation is a *return value*, so a cancelled
  fit stops at the last accepted point and no Python exception unwinds a C++
  loop holding raw buffers.

### Tests

`test/minimizer/test_minimizer.py` 31 tests (parity included);
`test/chi2`, `test/expression`, `test/sampler`, `test/portnode`,
`test/factorgraph`, `test/session`: **271 passed**. chisurf `test/fitting`:
**932 passed** (was 909; +14 `test_graph_fit.py`, +9
`test_leastsqbound_covariance.py`), the two known pre-existing failures
(`test_fit_state`, `test_pcf_experiment`), both confirmed against `HEAD`.

### What is NOT faster, said plainly

**A `FitGroup` is unchanged.** `global_optimize_local_first` ships as
`false`, so a group runs exactly one optimisation -- over `GlobalFitModel`,
whose residual is the *concatenation* of its members' -- and the graph
refuses it. Since a `FitGroup` is what the GUI builds, the 2.85x above is
today a plain-`Fit` result. Ticket `T-20260901-02`; the shape is a
concatenating residual node in C++ over one `ChiSquared` per member, and the
cross-member parameter *links* are already `Port::set_link`, i.e. already on
the right side of the boundary.

## 2026-09-01 (8) — the residual path in C++, for every model, and the design that had to be thrown away first

Owner: *"wire residual path in bff, for all mdls."* Done, **2.65x** on that
path. The route there is worth recording, because the obvious design was
slower than what it replaced.

**Why it covers all models.** `calculate_weighted_residuals` takes the model's
*curve*, not the model. Whatever computed that curve -- this library, numpy,
arbitrary Python -- the residual afterwards is the same operation, so a single
kernel serves the whole model tree and there is no per-model work.

**The first design lost, and measurably.** A `bff::ChiSquared` holding the data
and reused across iterations, with numpy setters and a per-data cache: it
measured **0.72x**, i.e. slower than numpy, interleaved so load was not the
explanation. The profile said the cache lookup was **4.45 us** of a 5.64 us
call while the C++ call it guarded was **0.61 us**. The cause is worth
knowing: **`data.y is data.y` is False** -- `Curve` stores `x` and `y` as rows
of one 2xN array and each access builds a fresh view -- so the identity check
never matched and a `ChiSquared` was rebuilt on *every* call. The stable
identity is the storage array or `.base`, not the view.

Rather than repair the cache, the design went. Caching the data in C++ also
means holding a *copy*, so a caller editing its curve through `unlocked()`
would be served residuals against stale data. The kernel is now a **stateless
free function** reading all three buffers where they lie: no node, no cache,
no copy, no staleness, and one boundary crossing. **2.65x at 512 points and
2.64x at 4096**, against the numpy expression *including* the slicing around
it, and `get_wres` as a whole 1.15x.

**A real divergence the edge cases caught.** With `xmin = -5`, Python slicing
means "five from the end" and returns *nothing*; the C++ window clamps to zero
and returned 100 residuals. `ChiSquared::resolve_window` has always had that
clamp, so this was latent in bff rather than introduced here. Negative `xmin`
now takes the numpy path -- reproducing Python's slice arithmetic in C++ would
only invite the next edge case, and the rare path may as well take the one
that defines the behaviour.

Verified on **80 combinations**: mismatched data/model lengths, empty,
inverted, past-the-end, zero-width and negative windows, both noise models,
**bit-identical every time**. Pinned by `ResidualKernelTests`, which compares
against the previous implementation rather than against itself, checks that the
C++ kernel is actually resolved (or the comparison passes by testing numpy
twice), and checks the numpy fallback still works, since it is the definition
of the behaviour.

`test/fitting` **906 passed** plus the two known pre-existing failures;
`test/expression` and `test/chi2` green in bff.

One measured caveat on the ambition. The residual is 16% of the callback, so
2.65x on it is 1.15x on the callback and less on the whole fit. **Crossing the
boundary is not free** -- a SWIG call is ~0.6 us before any work -- so moving a
*small* kernel across adds a crossing rather than removing one, which is
exactly what the first design demonstrated. The remaining large item is the
opposite trade: crossing **once** for a whole iteration instead of once per
part.

## 2026-09-01 (7) — the error estimates were wrong, by 22x, whenever the data had no sigma

Owner: *"the error estimates were/are all wrong."* Correct, and the defect is
specific enough to state exactly.

`covariance_matrix` returns ``(J'J)^-1`` for ``J = d(weighted residuals)/dp``.
That is the parameter covariance **only when the weights are real standard
deviations** -- weighted least squares, which is what
`calculate_weighted_residuals` documents for its ``"default"`` noise model and
what ``curve_fit(absolute_sigma=True)`` computes. Checked against scipy on a
400-point exponential with a known sigma: chisurf agrees to 0.3%, which is
finite-difference noise. **That path was never wrong.**

**`DataCurve` sets `ey` to ones whenever none is supplied** -- five sites in
`core/data.py`, including the ordinary two-column ``x, y`` file. The residuals
are then unweighted, the covariance comes out in units of "a residual of 1",
and the reported errors are too large by ``1/sqrt(chi2r)``:

| | b | a1 | t1 |
|---|---|---|---|
| reported, no `ey` | 0.078 | 0.291 | 0.375 |
| correct (sigma from residuals) | 0.0036 | 0.0133 | 0.0171 |

**A factor of 22.** Not a rounding disagreement -- a number a reader would act
on. And silent: nothing about a two-column file announces that its error bars
have become meaningless.

Fixed by estimating sigma from the residuals when the data supplies none,
which is the standard treatment and what ``curve_fit`` does by default. The
scale is applied *only* when ``ey`` is exactly one everywhere, so **a dataset
with genuine uncertainties is untouched and its error bars do not move** --
pinned by its own test, because a detector that fired too eagerly would
silently shift every published error bar in the other direction. Both regimes
now match scipy: 0.00391/0.01454/0.01875 against `absolute_sigma=True`, and
0.00357/0.01329/0.01714 against `absolute_sigma=False`.

The display layer is **not** at fault and was checked: `parameter_table.py`
formats `error_estimate` to 4 significant figures and blanks a non-finite one,
and `FittingParameter.error_estimate` returns NaN when nothing has been
computed, so an unknown error shows empty rather than stale. Two mapping
hazards were tested and are sound: a **fixed** parameter is absent from
`model.parameters` and reports no error while the free ones keep their correct
values, and a parameter with an identically zero derivative reports NaN rather
than inheriting a neighbour's column. A stale estimate surviving a re-fit is
reachable in principle -- `_error_estimate` persists and only
`important_parameters` are reassigned -- but not through a change of equation,
which rebuilds the parameter objects.

`test/fitting/test_error_estimates.py`, 7 tests, checks both regimes against
scipy rather than against chisurf's own output, and one of them guards the
*size* of the old defect so a future regression cannot pass by being 1% wrong.
`test/fitting`: **906 passed**, the two known pre-existing failures.

## 2026-09-01 (6) — where a fit's time actually goes, and why porting leastsq is the wrong lever

Owner: *"did you port leastsq to bff? that could add more speed."* No, and
having measured it: **it would not add much.** The optimiser is not the cost.

`fit.run()` at 512 points, three free parameters, measured without a profiler
(cProfile inflates this path about 2.3x):

| | ms | share |
|---|---|---|
| **residual callbacks (12x)** | **0.441** | **57%** |
| -- of which `update_model` | 0.186 | 24% |
| -- of which the rest of `get_wres` | 0.255 | 33% |
| **post-fit error estimation** | **0.247** | **32%** |
| MINPACK plus the `leastsqbound` wrapper | ~0.084 | ~11% |

MINPACK's `_lmdif` is compiled Fortran already: its *own* time is
`tottime 0.003` of 0.085 in the cumulative profile, about **4%**. Replacing it
with a C++ Levenberg-Marquardt would swap one compiled optimiser for another
and leave the 57% untouched, because that 57% is the *callback* -- the loop
returning to Python on every iteration to write parameters, evaluate the model
and form residuals.

**The productive port is therefore the residual, not the optimiser**, and it is
the thing the original handover already named: `Expression -> ChiSquared ->
Sampler` as one C++ graph, so an iteration never re-enters the interpreter.
bff already has `ChiSquared.h` and `Sampler.h`; what is missing is the wiring
from a chisurf model to them. That only helps models bff can represent -- which
parse models now are, and arbitrary Python models are not.

**Error estimation is a third of every fit**, and it recomputes a Jacobian.
Two observations, one acted on:

- `approx_grad` costs `1 (f0) + p (steps) + 1 (restore)` evaluations. The
  restore is necessary and its reason is documented in the source. **`f0` was
  not**: `covariance_matrix` takes the gradient around
  `xk = model.parameter_values`, and both callers of
  `update_error_estimates()` run `self.update()` immediately before, so the
  model already holds the residuals at `xk`. Threading `f0` through takes a
  fit from **13 model evaluations to 12**, and the covariance comes out
  **bit-identical** -- `np.array_equal` true, max difference exactly `0.0`,
  checked on a 3- and a 5-parameter model.
- **Honestly, that is worth about 1%, not the 8% I first estimated.** I had
  assumed model evaluations dominated the fit; they are 24% of it. One
  evaluation of thirteen is ~9 us of ~772 us, which is below this machine's
  noise -- 0.810 ms against 0.813 ms, indistinguishable. The change is kept
  because it is strictly less work for a provably identical answer, not
  because it is measurable.

Not done, and recorded as the alternative that was rejected on inspection:
MINPACK can return `cov_x` built from the `fjac`/`ipvt` it already computed,
which is what `scipy.optimize.curve_fit` uses, and it would remove the whole
32%. **`leastsqbound` transforms parameters when bounds are present**, so that
covariance is in internal coordinates and would need the transform's Jacobian
applied before it means anything externally. Silently wrong error bars are a
worse outcome than slow ones, so this needs doing deliberately rather than as
part of a speed pass.

`test/fitting`: **899 passed**, the two known pre-existing failures.

## 2026-09-01 (5) — the fit's cost was never the arithmetic

Owner: *"improve fit mdl in chisurf to be faster."* Profiled before touching
anything, and the profile said something worth recording: after the previous
entries, **evaluating the equation was 21% of a fit iteration and chisurf's
own machinery was the rest.**

First correction to my own method. A fit runs inside
`factorgraph.frozen_structure()`, which collapses parameter reads to a dict
lookup; measuring `update_model()` outside it measures a path no fit takes.
Unfrozen 17.81 us, frozen 9.57 us. Everything below is frozen.

Component breakdown at 512 points, FCS:

| | us | share |
|---|---|---|
| **`self.y = <array>`** | **4.15** | **45%** |
| `compute_curve_bound` | 1.92 | 21% |
| `numpy.array([p.value ...])` | 1.18 | 13% |
| `fit.data.x`, `ascontiguousarray`, `super()` | 0.50 | 5% |

So the single largest cost of a model evaluation was **assigning the result**,
at more than double the computation.

`Curve.y`'s setter wraps its write in `self.unlocked('d')`. Three things were
wrong with that, found by measuring each:

1. It was a `@contextlib.contextmanager`, so every assignment built and tore
   down a generator frame. Now a small class with `__enter__`/`__exit__`.
   Semantics unchanged, re-entrancy included.
2. **The real cost, and not where I first looked.** `self._unlock_depth = ...`
   went through `Curve.__setattr__`, which exists to lock any per-sample array
   on its way in -- with two further `__setattr__` overrides above it in the
   MRO. Assigning that private int cost **0.760 us** against **0.136 us** for
   the dict write, and it happens twice per unlock. It is bookkeeping, never
   an array, so none of that machinery applies to it.
3. Resolving `'d'` to an attribute name is a property of the class, asked
   again on every iteration. Cached per (class, arguments).

`ParseModel` also stopped rebuilding its parameter array per step (filled in
place; the size cannot change) and stopped re-checking the axis for contiguity
(the same object for a fit's life).

**Measured honestly, which took two attempts.** A plain before/after A/B was
useless: load on this machine drifted from 11.9 to 14.8 between runs, which is
larger than the effect. Re-done with both implementations **interleaved in one
process**, alternating and min-of-many:

> `self.y = v` -- **new 2.34 us, old 4.54 us, 1.94x faster**

That is per model evaluation for **every curve model in chisurf**, not only
parse models, since they all assign `y` through the same setter.

Tests: `test/core` + `test/fitting` **2053 passed, 5 failed**, and all five are
pre-existing. The two known ones from the previous entry, plus three in
`test_burst_manifest.py` and `test_tttrlib_registry.py` that I had not run
before. Established, not assumed: none of those files mentions `unlocked`,
`Curve`, `Expression` or `compute_curve`, and both changed files were stashed
back to HEAD and the three re-run -- they fail identically without them. The
registry two look downstream of another session's live tttrlib work, which is
also what `test_registry_completeness.py` has been failing on over there.
`test/core/test_curve.py`, `test_curve_locking.py` and the curve doctests are
green.

## 2026-09-01 (4) — a fit re-passed its variable names on every iteration

Owner: *"no repassing of str needed on fit."* Correct, and it was costing more
than it looks.

`compute_curve()` takes the variable names on every call, so SWIG rebuilt a
`std::vector<std::string>` -- allocating and copying each name -- per fit
iteration, and the engine then looked each one up again. Measured by holding
the axis at 8 points so the arithmetic vanishes and only call overhead
remains: **1 name 1.28 us, 3 names 1.44, 6 names 1.66, 10 names 2.14** --
about **0.1 us per name**, none of which can change while a fit runs.

`bind_parameters(names, axis)` now resolves the mapping once into a small
`vector<int>` -- one slot per engine variable, a sentinel for the axis -- and
`compute_curve_bound(values, axis)` evaluates against it with **no string
crossing the boundary at all**:

| equation | 8 pts | 512 | 4096 |
|---|---|---|---|
| `b+a1*exp(-x/t1)` | **57%** faster | **26%** | 6% |
| FCS | **59%** | **30%** | 7% |

The 512-point column is the one that matters: that is the length an FCS curve
actually has. Results are bit-identical to the named call (`rtol=1e-15`), as
an optimisation's should be.

End to end through `ParseModel.update_model()`, C++ against `eval()`:
**1.84x** and **2.41x** at 512 points, 1.39x and 1.70x at 4096.

Three things the binding has to refuse rather than guess, each with a test: a
parameter count that disagrees with the binding (it would read the right
equation from the wrong slots), evaluating before binding, and a binding left
over from a previous equation -- `set_expression()` drops it, because keeping
it would evaluate the new equation through the old slots. chisurf binds
against `_parameters_equation` rather than `_keys`, since the values array is
built from that same list, so the two cannot drift.

`test/expression` **79 passed**; `test_parse_uses_bff.py` 26 tests over 69
subtests; chisurf `test/fitting` **899 passed, 2 failed** -- the same two
pre-existing failures established as unrelated in the previous entry.

**The bff suite cannot be reported clean, and it is not ours.** It ran 1677
passed / 12 failed, and after rebuilding against current sources 9 remain, in
`test_av_lattice.py`, `test_docking_values.py`, `test_fps_screening_ab.py`,
`test_fps_point_positions.py` and `test_ProbeNetworkRestraint.py`. Grounds for
saying they are somebody else's: none of those files mentions `bff.Expression`,
`compute_curve` or `bind_parameters`; `test/expression` is 79/79; the earlier
run this session -- before this change -- failed **four** tests and all four
were FPS *export*, a different set entirely; and `src/AV.cpp` was modified at
**07:11**, during the run, with `Docking.h` at 00:43. The failing areas track
another session's live edits, not this one's.

## 2026-09-01 (3) — the integration test, and proof the guard bites

The gap left by the previous entry is closed. `update_model()` is now driven
through a real `Fit`, which is what it needed: a `ModelCurve` is constructed
*by* a `Fit` and reads its axis from `fit.data`, so a stub cannot assign
`self.y` -- the pattern already used by `test_models_regression.py`.

Six tests, and the distinction from the `parse_code()` sweep matters: that one
proves an equation **compiled** for the engine, this one proves the fit
**uses** it. A model can compile and still fall back on every iteration if
`compute_curve` throws, because the fallback catches and logs rather than
failing. Nothing else would report it.

- one update evaluates in C++ and not in Python -- counters `(1, 0)`
- **fifty consecutive iterations stay in C++**, none drifting onto the
  interpreter, which is the shape a real fit has
- the curve equals the one numpy computes, because fast is worthless if wrong
- the `xD` model that used to raise `NameError` on every evaluation now
  produces a finite curve end to end
- **all 67 shipped equations evaluate in C++ through a real fit**
- and the guard bites: dropping the compiled expression the way a regression
  would makes the counters report `(0, 1)` and `evaluates_in_cpp` False,
  rather than the model carrying on looking healthy

That last one is there because the others pass vacuously without it -- a
counter that never increments would satisfy every assertion above.

`test/fitting/test_parse_uses_bff.py` is now 20 tests over 67 subtests; with
the two existing parse test files, 32 pass. `test_models_regression.py` passes
again as well (3) -- it had failed to collect earlier only because of another
session's in-flight FPS work, which has since been rebuilt.

Whole of `test/fitting`: **893 passed, 2 failed**. The two --
`test_fit_state.py::test_fret_gaussian_model_get_set_state_preserves_gaussians`
and `test_pcf_experiment.py::test_pcf_config_block_present` -- are
**pre-existing and unrelated**. Established rather than asserted: neither file
mentions `ParseModel`, `parse.` or the attributes touched here, and
`parse.py` is the only chisurf file changed, so it was stashed back to HEAD
and both tests were re-run. They fail identically without the change. Stash
popped and the file verified byte-identical to the working version afterwards.

## 2026-09-01 (2) — asserting chisurf parses with bff, and the two bugs that assertion found

Owner: *"assert that chisurf parse uses bff parsing, so that compute stays in
bff for max speed."* The assertion is
`chisurf/test/fitting/test_parse_uses_bff.py`, 14 tests over all **67** shipped
equations. Writing it was worth more than the assurance: it found two
pre-existing defects, neither caused by the move to C++.

**Why an assertion was needed at all.** `ParseModel` keeps `eval()` as a
fallback on purpose -- an equation may use a numpy or scipy function the
engine lacks. That fallback is silent by design, so anything which stopped the
C++ path compiling (a renamed method, a missing `IMP.bff`, an engine that
rejects a spelling) would break nothing: every fit would just go back to an
interpreter round trip per iteration and nobody would find out except by
profiling. So `evaluates_in_cpp` and `evaluation_counts` were added to make
the path taken observable, and the suite pins that **every shipped equation
parses onto the engine**, driving the real `parse_code()` rather than reaching
past it.

The sharpest test is not "does it compile" but **"do the scanner and the
engine find the same free names"** -- a disagreement there means the two
evaluators bind different things and quietly return different curves. It
failed on three equations.

**Defect 1: `pi` was a fit parameter, initialised to 1.0.** ChiSurf's scanner
rewrote every free name to `a[i]`, `pi` included, so three shipped
dye-diffusion models computing `...-4*pi*Rdye*Ddye*Nq**2/Vav*x*...` were
running with **pi = 1** -- a factor of pi wrong in the quenching term --
unless a user noticed the stray "pi" in the parameter table and typed
3.14159 into it. `models.yaml` gives it no `initial:`, so nothing corrected
it. `pi` and `e` are now left alone: `eval()` picks them up from the
`from numpy import *` already at the top of that module, which is also how the
engine resolves them, so both paths now agree **and are right**.

**Defect 2: a name beginning with `x` was split.** The scanner's first rule
was the bare `r"x"`, which matched the leading character of any such name, so
`xD` became the axis `x` followed by a parameter `D` and the generated code
read `xa[1]`. The shipped two-state quenching model
`p0*((1-xD)*(...)+xD*(...))` therefore failed with
**`NameError: name 'xa' is not defined` on every evaluation** -- it has been
broken outright, not subtly wrong. Narrowed to `r"x(?!\w)"`. Confirmed
against the old generated string before fixing, and the repaired model now
produces a finite curve end to end.

Both are pinned by their own regression tests. Existing parse tests
(`test_parse_widget_expression_editor.py`, `test_parse_latex.py`) still pass,
12 of them.

`update_model()` is now driven too, through a real `Fit` -- see the following
entry; the gap noted here is closed.

## 2026-09-01 — Table retired, and the last of ExprTk with it

Owner: *"can retire table, i guess."* Done. Independently re-verified the
premise before deleting anything: `bff.Table`/`bff::Table` across imp.bff,
chisurf, imp-tricks, tttrlib, quest, ucfret and fpsimp returns **three hits,
all in its own test file**.

Deleted: `include/Table.h`, `src/standalone/Table.cpp`, `test/table/`, and
**`include/internal/exprtk.h` -- 1.6 MB**, which Table was the last thing
keeping alive. `Expression::compute_pointers()` went too: Table was its only
caller and it was `%ignore`d from Python, so nothing outside could reach it.
`grep -rl exprtk` over every header, source, interface and CMake file in the
repository now returns **nothing**.

Verified after: `bff.Table` gone from Python, `compute_pointers` gone,
`Expression` intact; `test/expression` 72 passed; full suite **1656 passed**.

Four failures in that run are **not ours** and are actively in flight:
`test/io/test_fps_export_formats.py` (3) and
`test/representation/test_fps_error_estimation.py` (1), from the session
writing the FPS export work -- `src/FPSExport.cpp` was modified at 00:18 and
the second test file at 00:20, mid-run. Neither test touches `bff.Table` or
`Expression`; the one "Table" match in them is `IMP.bff.FPSResultTable`, a
different class.

**A build-system trap worth recording, because it cost most of the time here
and had nothing to do with the deletion.** The module went from building to
`ninja: error: unknown target 'IMP.bff-lib'`, which reads like the deletion
broke the build. It did not. Another session had added
`bin/imp_bff_fps_export` without a `README.md` section, and this repository
*disables the whole module* when a program in `bin/` is undocumented --
`setup_module.py` enforces it, and the symptom is a missing target, not a doc
warning. That is already recorded in the 2026-08-19 entry; it is recorded
again because it presented as somebody else's regression. Written a section
for the tool from its own help text, which re-enabled the module for everyone.

A second, smaller one: after editing `src/Files.cmake`, `build.ninja` kept a
stale `standalone/Table.cpp` rule through two full `cmake .` runs and only
dropped it after touching `Files.cmake` and `src/CMakeLists.txt`. If a deleted
source is still being compiled, the generate step has not noticed --
re-touching the file that lists it is the fix.

## 2026-08-31 (21) — chisurf's parse models evaluate in C++; ndx filters measured end to end

Owner: *"make sure that ndx uses tttrlib for filter and is fast, chisurf eval,
parser, must use bff. check speed."* Both done, and the chisurf half needed a
new entry point before it was worth doing at all.

### chisurf: `eval()` -> `IMP.bff.Expression`

`ParseModel.update_model` ran `eval(self.code)` on every fit iteration, over a
string a `re.Scanner` had rewritten to replace each free name with `a[0]`,
`a[1]`... It now compiles the **original** equation with `bff.Expression` --
bff binds by name, so the rewrite is unnecessary and would only hide the names
from it -- and evaluates in C++. `eval()` stays as a fallback, deliberately:
an equation may use a numpy or scipy function the engine does not implement,
and refusing it would break a model a user already has.

**The obvious wiring was slower than the interpreter, which is why this is not
a one-line change.** Measured, for the shape a fit actually has (N *scalar*
parameters and one vector axis):

| entry point | 2exp, 512 | 2exp, 4096 | FCS, 512 | FCS, 4096 |
|---|---|---|---|---|
| python `eval` | 7.8 us | 30.5 | 6.6 | 13.0 |
| `compute()` | 19.1 | **134.8** | 16.6 | **117.0** |
| `compute_columns()` | 5.6 | 34.5 | 3.5 | 18.1 |
| **`compute_curve()`** | **5.3** | **31.1** | **2.9** | **13.2** |

`compute()` takes `std::vector`, so from Python every call converts the axis
to a list -- four times *slower* than the interpreter it replaces.
`compute_columns()` reads numpy directly but wants every operand the same
length, so a scalar has to be materialised as a full column and the memcpy per
operand costs more than the arithmetic. Both lose at 4096 points.

So `Expression::compute_curve()` was added: scalars stay scalars, the axis is
read where it lies, nothing is copied. The engine already broadcasts a scalar
-- its `is_vector` flag is exactly that -- so this stops lying to it about the
shape of the data. It is the fastest of the four at every size measured.

**Across the whole shipped catalogue** -- 67 unique equations -- against
`eval()`: **zero mismatched, zero refused**, 1.38 ms -> 0.58 ms at 512 points
(**2.40x**) and 4.60 -> 3.66 ms at 4096 (**1.26x**). The margin is larger at
the short curves, which is where FCS and most fits live.

### ndx: already tttrlib, now measured through its own API

`DataSource.query_mask` was already delegating to
`tttrlib.DataStore.select_expression`; the `np.where` calls elsewhere in
`data_source.py` are log-scaling and histogram binning, not filtering. What
was missing was a measurement through the API the CLI and GUI actually call
rather than through `DataStore` directly:

| rows | `query_mask` | vs pandas |
|---|---|---|
| 100k | 0.116-0.164 ms | **5.4-7.0x** |
| 1M | 1.26-2.14 ms | **2.0-4.9x** |

Every query checked against `frame.eval` first; no mismatches. These are
slower than the raw `count_expression` figures (0.076-0.092 ms at 100k)
because `query_mask` also saves and restores the store's own selection and
copies the mask out to numpy. That is the honest cost of the call a user
makes, and it is the number worth quoting.

bff **1642 passed, 4 xfailed, 0 failed**.

Two cross-session tears hit during this and are worth recording as evidence
the board's hazard section is not theoretical: `fps_bootstrap` and then
`FPSMolecule_path_get` were undefined for a while because another session was
editing `Docking.h` and `FPSExport.h` (the latter at 23:12, mid-run). Both
cleared on a rebuild, neither was ours.

## 2026-08-31 (20) — the engine is a copied header; imp.bff depends on nothing

Owner: *"simply copy the header over."* Done, and it is the right shape for
this component.

`include/internal/ExpressionEngine.h` is now a **byte-identical copy** of
tttrlib's `modules/core/include/ExpressionEngine.h`, included as
`<IMP/bff/internal/ExpressionEngine.h>` — the same convention the vendored
ExprTk already used. imp.bff's build now needs **nothing outside IMP**:

| | before today | after |
|---|---|---|
| own evaluator | 1,891 lines + 1.6 MB ExprTk | none |
| link edge to tttrlib | — | **none** (`otool -L` → 0) |
| build dependency | — | **none** (`dependencies.py` clean, `dependency/` gone) |

The intermediate designs were both heavier than the problem. Linking
`libtttrlib` coupled two build systems, an ABI and an install contract to a
computation that opens no file and calls nothing outside libm. Taking the
header from an installed tttrlib was lighter but still made this repository
need tttrlib *present* to build a file with no dependencies of its own.

**The hazard of a copy is drift, so it is a test rather than a hope.**
`test/expression/test_engine_copy_is_identical.py` compares the two files by
SHA-256 against the sibling checkout. tttrlib owns the original; changes go
there and are copied here, never the reverse. Verified the guard actually
bites: a two-line edit to the copy fails it, restoring passes. It *skips*
rather than fails when the sibling checkout is absent, so a release tarball or
a bff-only CI job is not reported as broken for a comparison it cannot make.

That leaves one implementation of the language in the ecosystem, in tttrlib,
with a mechanically-checked copy in imp.bff — which is what "not maintaining
two split code bases" actually asks for. A second *file* is not a second code
base; a second *implementation* is, and there is no longer one.

bff **1641 passed, 4 xfailed, 0 failed** (72 in `test/expression`, including
the two drift guards). No performance change: FCS 2.86x / 2.31x / 1.69x /
1.28x / 0.98x, 12.5 us at 4096 points.

Still standing, unchanged: `include/internal/exprtk.h`, 1.6 MB, kept alive
only by `src/standalone/Table.cpp`.

## 2026-08-31 (19) — header-only: one engine, and no link edge at all

Owner, on the previous entry: *"but i do not get, why bff should depend on
tttrlib, they should just use the same header files. make the stuff header
only!"* Right, and the earlier design was heavier than the problem.

`ExpressionEngine` is pure arithmetic — a tokeniser, a parser and SIMD kernels
over arrays. It opens no file, holds no resource and calls nothing outside
libm. Making imp.bff **link** `libtttrlib` for that coupled two build systems,
an ABI and an install contract to a computation that needs none of them.

So the engine is now **header-only**: `modules/core/src/ExpressionEngine.cpp`
is gone, its 1,593 lines merged into `ExpressionEngine.h` (1,900 lines), the
eight non-template member definitions marked `inline` (`run` is a template and
must not be). Both callers include the one file —
`tttrlib::data::DataStore` gates burst columns with it, `IMP::bff::Expression`
evaluates model equations with it — and neither depends on the other's
libraries.

**Verified there is no link edge**: `otool -L libimp_bff.dylib | grep -c
tttrlib` returns **0**. imp.bff's `dependency/tttrlib.description` now
declares `libraries=""` with only a header.

What remains is a build-time need for tttrlib's *headers* to be present, which
is a much lighter thing than a library dependency but is not nothing — bff
includes `<tttrlib/ExpressionEngine.h>` from the conda prefix.

**An honest note on T-20260831-14.** The `install(EXPORT)` /
`tttrlibConfig.cmake` work done an hour earlier is **not needed for this**. The
header install already existed and is what is actually used; `find_package`
support was built for a link dependency that no longer exists. It is kept
because it stands on its own — there was previously no supported way to be a
C++ consumer of tttrlib at all — but it was not on the critical path, and
saying otherwise would be dressing up a detour as a plan.

No regression, and if anything a shade faster now that the compiler can inline
across what used to be a translation-unit boundary: FCS 3.05x / 2.29x / 1.66x
/ 1.21x / 1.00x against numpy, 12.4 us at 4096 points (12.5 linked, 13.0 with
bff's own engine). bff **1639 passed, 4 xfailed, 0 failed**; tttrlib's
DataStore expression tests 49 passed; 12,000 fuzzed queries agreeing with
numpy, `known-divergence=0`.

Still standing: `include/internal/exprtk.h`, 1.6 MB, kept alive only by
`src/standalone/Table.cpp`. Retiring Table removes the last of it.

## 2026-08-31 (18) — one engine: imp.bff now evaluates with tttrlib's

**T-20260831-14 then -12.** `src/standalone/Expression.cpp` went **1891 lines
-> 300**. What is left is the part that is genuinely imp.bff's: presenting the
evaluator as a `Node` with ports, and the numpy-facing entry points SWIG
wraps. The tokeniser, parser, RPN compiler, folding, CSE and the block
evaluator are gone from this repo; `tttrlib::data::ExpressionEngine` is the
only copy now.

**The install contract (T-14).** A correction first: this log said tttrlib
"installs no C++ headers". It does — `BUILD_LIBRARY` and `INSTALL` default ON
and a plain install lays down 172 headers. The *wheel* sets
`BUILD_LIBRARY=OFF`, and I generalised from the wheel. The real gap was that
`find_package(tttrlib)` failed: no `install(EXPORT)`, no config file. Fixed
with `EXPORT`/`INCLUDES DESTINATION`, a generated `tttrlibConfig.cmake` and a
`SameMajorVersion` version file. Two snags: the vendored `tiff` is a
build-tree target in no export set and blocked the whole export (tttrlib links
third-party libraries at *directory* scope, so they land in every target's
interface) — cleared the interface, which is accurate rather than a dodge
since a consumer of the dylib does not re-link what it already records; and
`EXPORT_NAME`, or the targets export under internal names. Proven with an
out-of-tree consumer that `find_package`d, linked and ran. Installed into the
arm64 conda env at the owner's direction, after checking nothing collided.

**Wiring (T-12).** `dependency/tttrlib.description` plus `tttrlib` in
`dependencies.py`; IMP's own machinery does the rest, and `build_info/tttrlib`
reports `ok=True`.

**The one real regression, and its cure.** Dropping ExprTk cost the functions
only it implemented — `floor(x/2)` was the case that proved it. Rather than
accept the loss, the missing set went into the engine: **18 unary and 2 binary
functions**, numpy-spelled. So `hypot(x,y)` now returns
`[1.118 2.5 3.905 5.315]` — the values it was *silently getting wrong*
before — and T-20260831-10 is cured rather than guarded. Two details worth
keeping: `round` is `std::nearbyint`, half-to-even as numpy rounds, not
`std::round`'s half-away-from-zero; and `sign` returns numpy's -1/0/1 with NaN
for NaN, not `copysign`, which has no zero case.

Four wrapper bugs found by the existing tests and fixed, each of which the
tests were right to catch: the curve goes to the output port keyed by the
node's **own name**, not a fixed `"value"`; a missing input port throws rather
than returning quietly; `get_number_of_compilations()` counts trips to the
parser *for the equation currently held*, so it is 1 after a new equation, not
cumulative; and a **failed** `set_expression` must leave the previous program
intact — compiling in place would let a mistyped equation destroy a running
fit, so it compiles into a probe and commits only on success.

Verified: bff **1639 passed, 4 xfailed, 0 failed**; tttrlib 3457 passed (its
one failure is another session's image-kernel WIP); 22,000 fuzzed queries all
agreeing with numpy with **known-divergence=0**, so T-20260831-09's min/max
NaN bug is gone from bff too, by construction. **No performance regression on
the model-curve path**: FCS 2.95x / 2.29x / 1.73x / 1.25x / 0.99x against
numpy, 12.5 us at 4096 points where bff's own engine measured 13.0 — the port
carried the constant-folding work.

Not finished: `include/internal/exprtk.h` is **still there, 1.6 MB**, because
`src/standalone/Table.cpp` still includes it. Table has zero consumers outside
its own test (T-20260831-07 recommended retiring it), so deleting both is the
last step — but that removes a public class, which is the owner's call.

## 2026-08-31 (17) — the wrong curve now refuses, and what actually blocks one codebase

Two things, one delivered and one a correction.

**T-20260831-10 fixed.** The ExprTk fallback returned a silently constant
curve for any multi-argument function — `hypot(x,y)` flat where numpy varies,
`if(x>2,1,0)` all zeros — with no exception and no warning. It now **refuses**,
which is what `Expression.h` promised all along: an equation that cannot be
compiled is refused *"so a caller can fall back rather than get a wrong
curve"*. Returning the wrong curve was the one outcome the design ruled out.

The guard keys on **arity, not on falling back**, which is the distinction
that makes it safe: unary functions vectorise correctly in ExprTk, so
`floor(x/2)` still works. `min`/`max`/`pow` still answer because the vector
engine implements them — but they are on the guard list anyway, since
`floor(min(x,y))` drags the whole expression into the fallback and takes `min`
with it. Identifier boundaries are checked, so `summary + xmin` is untouched.
No shipped equation is affected: all 86 take the vector path. Five regression
tests; `test/expression` 70 passed, full suite **1605 passed, 0 failed**.

**Correction on T-20260831-12, the ticket that ends the duplicated engine.**
An earlier note called its feasibility "settled" because an IMP module may
carry its own `CMakeModules/Find*.cmake`. That was half the question and the
wrong half: `find_package` needs something installed to find, and **tttrlib
installs no C++ headers and exports no CMake package config**.
`tttrlib_install_modules` (`cmake/TTTRLibModule.cmake:420`) installs
`LIBRARY`/`ARCHIVE`/`RUNTIME` only, into the Python wheel's "bindings"
component; the conda prefix has no tttrlib headers, and the only C++ artefact
anywhere is a dylib inside a build tree. GSL, FFTW3 and HDF5 are *system*
libraries with installed headers — the precedent covers the mechanism, not the
availability.

So imp.bff cannot link tttrlib today, and that is not an imp.bff problem:
**there is currently no supported way to be a C++ consumer of tttrlib at
all.** Filed as **T-20260831-14**, the prerequisite for -12, with the wheel's
component split flagged as the thing not to break. Which route to take —
a real install contract versus compiling tttrlib's one source file from the
sibling checkout — is left as an owner decision rather than an agent's.

## 2026-08-31 (16) — the engine moves down into tttrlib, and gating halves

User directive: *"must work on tttrlib, not going to maintain two split code
bases; must be fast on pto and ndx."*

The direction was forced rather than chosen. Layering is tttrlib → imp.bff and
`imp.bff/dependencies.py` names IMP modules only, so `DataStore` **cannot**
call `bff::Expression`; and `AGENTS.md`'s placement rule — *"photons/curves →
tttrlib"* — points the same way for both consumers, since burst columns are
photon-derived and a model curve is a curve. So the evaluator core goes down
into tttrlib as `modules/core/{include,src}/ExpressionEngine.{h,cpp}`, and
`bff::Expression` is left to become a thin `Node` wrapper.

**Gating is about twice as fast, and the comparison is conservative** — the
baseline was taken at load average 9.7 and the new numbers at 14.3:

| rows | before | after | vs pandas after |
|---|---|---|---|
| 100k | 0.154-0.187 ms | 0.076-0.092 ms | 10-16x |
| 1M | 1.576-1.956 ms | 0.755-0.910 ms | 3.3-4.2x |
| 5M | 7.905-9.622 ms | 3.831-4.543 ms | 2.5-4.4x |

Two things bought that. The engine is **templated on the working type**, so
float32 columns — which is how ndxplorer stores burst parameters — run in
float32 rather than being widened, giving twice the SIMD lanes and an answer
bit-identical to pandas over the same columns. And `compute_mask` packs
`BitMask` words **straight out of the block loop**, where the old code wrote a
float per row and packed bits in a second pass. One exception, recorded rather
than buried: `g != 3` on float64 went 0.802 → 0.835 ms, the only cell in the
table that did not improve.

**pto is now measured at all**, which it was not before: a store written to a
PTO file and read back gates at 0.759-0.895 ms at 1M rows, indistinguishable
from the in-memory store. Coming off disk costs the gate nothing.

Correctness: 49 new DataStore expression tests; the grammar fuzzer ported
across and run for 15,000 random valid queries with **every one agreeing with
numpy**; full tttrlib suite 3,457 passed. Two failures were run down and are
**both other people's**: `test_registry_completeness` lists only image-kernel
symbols from another session's untracked WIP, and ndxplorer's three
`_bff_table` failures reference an attribute that occurs zero times in
`data_source.py` — stale since the `bff::Table` migration. The three ndxplorer
tests that actually check queries against pandas pass.

**T-20260831-09 was fixed in the port rather than carried over**: `min`/`max`
propagate NaN in both operand orders now.

Two things this did *not* finish, both filed. ExprTk survives in tttrlib as a
fallback and **carries the same silent multi-argument bug** found in bff this
afternoon — `hypot(g,r) > 5` keeps 8 rows where numpy keeps 6, `atan2(g,r) >
1.0` keeps 0 where numpy keeps 5, while `pow`/`min`/`max` are right precisely
because the new engine handles them (**T-20260831-13**). And bff still has its
own copy of the engine, so right now there are *two* — **T-20260831-12** is
the ticket that ends that, and until it lands the user's actual request is not
yet met.

The subagent doing this work hit its session limit and died mid-task; this
session took over its build lock, finished, and verified every number above
independently rather than accepting the report.

## 2026-08-31 (15) — the FCS loss was never the thing the ticket said it was

T-20260831-05 asked for common subexpression elimination, on the handover's
reading that FCS recomputes `x/1.2`. **It does not.** In
`0.3+1/1.7*(1+x/1.2)**(-1)/sqrt(1+1/2.1**2*x/1.2)` the second occurrence is
`(1/2.1**2*x)/1.2` -- a structurally different subtree -- and the only shared
thing is the leaf `x`, which is not worth caching. The whole of that row's
loss was **`**(-1)`**: it tokenises as a *negated* constant (`OP_CONST 1`,
`OP_FUN F_NEG`), and the constant-power fold only recognised a bare
`OP_CONST`, so the exponent stayed a buffer and every element paid a
`std::pow()`. Isolated: `x**(-1)` cost 39.7 us at 4096 points against 5.1 us
for `1/x`. numpy never pays it because it rewrites `arr ** -1` into
`np.reciprocal` before the ufunc ever runs.

Fixed by folding **every** constant subtree before the power fold, using the
same `apply_scalar_fun`/`apply_binary_scalar` the evaluator's own scalar path
uses -- so the arithmetic is bit-identical and only its timing moves.

FCS: **1.76x / 0.88x / 0.56x / 0.34x / 0.25x** -> **3.10x / 2.36x / 1.68x /
1.22x / 0.95x**, or 46.5 us -> 13.0 us at 4096 points. Those are this
session's independent re-run at load 9.8; the implementing session measured
3.2 / 2.4 / 1.75 / 1.29 / 1.01 at load 14. The first four cells agree, the
last does not: **4096 points is parity, not a win**, and both sessions
independently flagged that same cell as the doubtful one. The other three
curve rows did not move, and for `0.3+2.0*x` that is provable rather than
statistical -- its compiled program is instruction-identical.

CSE landed too and **earns its place on its own merits, none of them this
row**: A/B on the same binary with the cache disabled, at 4096 points,
`exp(-x/1.5)*exp(-x/1.5)+exp(-x/1.5)` 41.8 -> 17.5 us and a 3-state FRET-FCS
with 29 repeated subtrees 81.4 -> 46.6 us, while equations without real
repeats are neutral. It only shares a subtree where recomputing costs more
than the two block copies sharing costs, and `OP_SAVE`/`OP_LOADC` carry the
slot's *type*, which is what keeps it correct against the typed-stack
invariants from entry (12).

`test/expression` is **65 passed** (was 55). The grammar fuzzer returns
identical known-divergence counts before and after the change (21 and 17 on
seeds 1 and 2), which is the strongest evidence that nothing about the
semantics moved.

**The full suite is 1580 passed, 4 xfailed, and 1 failed** --
`test_docking_values.py::test_screening_ranks_a_library_and_never_reports_a_silent_nan`,
which asserts a screening CSV header of `["pdb", "score"]` and now gets one
ending `sigma3`. Not ours: `include/Docking.h` was modified at 19:01, mid-run,
by another live session adding columns. The same suite was 1569/0-failed at
18:05. Recorded rather than waved away, because "another agent did it" is the
easiest excuse in a shared tree and it should carry a timestamp.

Follow-on **T-20260831-11**: let a stack slot alias a column instead of
memcpy'ing it -- the largest remaining lead, worth two of FCS's ~11 block
passes and most of `0.3+2.0*x`. Build lock released; **T-20260831-09** and
**-10** are unblocked.

## 2026-08-31 (14) — two open questions closed, and a wrong curve nobody could see

Both read-only tickets came back, and one of them found a bug worse than the
one the fuzzer found this morning.

**T-20260831-06, is ExprTk worth keeping: no, after four named gaps.**
(`okf/validation/exprtk_fate.md`.) All 86 shipped equations take the vector
path and none falls back, so ExprTk is a validation gate and a dead branch.
It costs **45.3% of the shipped dylib's `__TEXT`** and takes `Expression.cpp`
from ~2 s to ~23 s to compile, twice over counting `Table.cpp`'s
`exprtk<float>`. Runtime is not an argument in either direction. The agent
also ported `tokenize` + `compile_vector_program` to Python to observe the
vector parser's verdict on input `compile()` rejects before it can be seen —
1,753 comparable inputs, 1,753 agreements.

**T-20260831-10, and this is the real find: the fallback returns a silently
constant curve.** ExprTk's *multi-argument* functions collapse to element 0
and broadcast it. `hypot(x,y)` gives `[1.118 1.118 1.118 1.118]` where numpy
gives `[1.118 2.5 3.905 5.315]`; `atan2` likewise; `if(x>2,1,0)` returns all
zeros. Unary functions are fine, so the split is arity, not vectorisation.
Reproduced independently before filing. No exception, no warning — a
plausible flat curve, which is exactly the outcome `Expression.h` designs
against when it says an uncompilable equation is *refused* "so a caller can
fall back rather than get a wrong curve". Unreachable from the shipped
catalogue, which is why 1,569 tests never saw it, and reachable by anyone who
types `hypot`, `atan2`, `if`, `clamp`, `root` or `inrange`.

**T-20260831-07, `bff::Table`: retire it.**
(`okf/validation/table_vs_datastore.md`.) Zero consumers outside its own
test; ndxplorer's `query_mask` already calls tttrlib's `select_expression`.
The layering rule agrees on both halves, so the tie-break never fires.
Deletion is confined to imp.bff. Two defects found in passing, both the
signature of code nobody runs: ten unreachable lines in `count()`, and a test
class defined after `unittest.main()`.

Both notes correct this repo's own documentation: the claim that ExprTk still
handles "comparisons, booleans" is stale in the handover *and* in the comment
at `Expression.cpp:1514` — the vector engine has done both since the typed
boolean stack landed.

T-20260831-05 (CSE) is still running and holds the build lock; it has taken
`Expression.cpp` from 1,561 to 1,891 lines adding `OP_SAVE`/`OP_LOADC`.
Tickets 09 and 10 are both blocked on it, since all three touch that file.

## 2026-08-31 (13) — a fuzzer that reaches the evaluator, and the bug it found in ten minutes

The handover's fuzz harness lived in `/tmp/bang.py` and threw random
*characters* at the parser. 1,000,000 inputs, zero crashes — but re-run and
measured, **193,951 of 200,000 inputs are refused and only 2.8% reach the
evaluator at all.** It was testing `is_supported`, not the arithmetic behind
it. Now `test/expression/fuzz_expression.py`, with the structural mode kept
and a second one added: build a random expression *tree*, render it twice —
once in the engine's syntax, once as the numpy expression that means the same
thing — and demand the two agree. That reaches the evaluator on 92% of cases.

It found a real bug on the first serious run. **`min`/`max` are not
commutative under NaN**: `min(y, nan)` returns `y`, `min(nan, y)` returns
`nan`. The kernel is a plain ternary `(a > b) ? b : a`, which is neither
numpy's rule (propagate, either order) nor C's `fmin` (ignore, either order),
while `Expression.h` states the contract as "semantics follow numpy".
Confirmed on both the SIMD body and the scalar tail. Filed as
**T-20260831-09**, not fixed: T-20260831-05 (CSE) owns `Expression.cpp` and
the build lock, and editing under it is the exact collision the board's
hazards section documents.

Two things had to be built before the fuzzer could be trusted, and both are
the interesting part:

- **The oracle was wrong first.** Rendering constants as bare Python literals
  makes `(-3.7) ** 2.1` a *complex number* — Python's semantics — where C's
  `pow` and numpy both give NaN. That manufactured failures that looked like
  engine bugs. Constants now render as `np.float64(...)`.
- **Ill-conditioned cases get no verdict.** `sin((z**z) ** (z/2.192))` reaches
  an argument near 8e7, where one ulp in becomes 1e-7 out, so the engine and
  numpy differ in the seventh digit while both are right. Rather than loosen
  the tolerance for trig — which would hide real bugs everywhere else — the
  harness nudges every input by one ulp and re-evaluates the twin: if numpy
  disagrees with *itself* by more than the tolerance the engine is held to,
  the case is unanswerable and is counted, not reported.

Result: **75,000 valid expressions across five seeds, zero unexplained
failures**, with the two explainable categories counted separately
(`known-divergence`, `ill-conditioned`) so a genuinely new failure is visible
the moment it appears. Every case also cross-checks the byte-mask evaluator
against the double one — the invariant from entry (12), and the one most
exposed to a stack optimisation like CSE.

Work handed out on the board: **T-20260831-05** (CSE for FCS-shaped
equations, holds the build lock), **T-20260831-06** (whether ExprTk still
earns its 1.6 MB), **T-20260831-07** (`bff::Table` vs `tttrlib::DataStore` —
**came back "retire": the class has zero consumers outside its own test**, and
the migration it was written for already went to DataStore), plus
**T-20260831-08** advertised unowned for the tttrlib gating port, which needs
coordination with the live tttrlib sessions rather than a second agent in the
same build.

## 2026-08-31 (12) — the gate answers in bytes, and three bugs in the typed stack

`Expression::compute_mask()` was the one item the expression-engine handover
marked **known-bad**: it allocated a full `std::vector<double>`, ran the
ordinary evaluator into it and narrowed afterwards -- strictly more work than
the plain evaluation it existed to beat, so every mask measurement taken
before today measured a stub.

It is now what the handover predicted it would be. The block loop already
carried a typed boolean stack, so `evaluate_vector_program` and
`evaluate_mask_program` collapsed into one `evaluate_program(..., double* out,
unsigned char* out_mask)` differing only in the block's closing store.
Measured at **1.3-2.2x the double path** and 1.3-9.1x pandas
(`benchmark/expression_mask.py`), with the margin over doubles *growing* with
row count -- which is what a saving in output bytes ought to do.

Making that path live exposed **three bugs of one shape**: the typed stack's
slot type was never reconciled at its boundaries. `is_bool` survived a push
that reused the slot, so `a>0 and b<1 or c>2` loaded `c` into a slot still
claiming to hold `b<1`; `and`/`or`/`not` read the boolean stack even when an
operand was a column of doubles; and a comparison used as a number read the
double stack, so **`(x>2)*3` computed `x*3`** -- wrong in shipped behaviour,
not merely latent. Fixed with a `booleanise()`/`numerify()` pair at every
boundary, and truthiness pinned to numpy's rule (nonzero is true, so a
negative and a NaN are both true; the stub's `> 0.5` threshold had both
wrong).

`test/expression/test_expression_mask.py`, 21 tests, includes a cross-check
that the mask and double evaluators can never again disagree. Suite: **1569
passed, 4 xfailed**. Load average was 12-18 throughout from concurrent builds
and test runs, so the timings are a floor, not a measurement.

## 2026-08-31 (11) — the labelizer port re-checked against the reference

Owner: *"check again the status against the labelizer."* After a week of
renames around it -- `DyeLibrary` -> `ProbeLibrary`, `dye_library.cif` ->
`probe_library.cif`, `find_dye` -> `find_probe`, and the whole `cgdye` ->
`cgprobe` move -- the question is whether the port still reproduces what it
claimed on 2026-08-24.

**It does, number for number.** `benchmark/labelizer_ab.py` on 1DDB-39 prints
what `okf/validation/labelizer_ab.md` recorded: `cr` 195/195 exact, `ss`
195/195 exact, `se` 146/195 bin-exact with bias +0.015 A, r = 0.9763, 92.3 %
within a bin, and `cs` off by the same constant 0.72310 that the shipped
example's own feedback loop produces. `test/label/` is 316 passed. The port's
four translation units include only `AV*.h`, `SolventAccessibleSurface.h`,
`StatesDistance.h`, `Pto.h` and `PtoProfile.h`, which is why the restructure
went past it without touching it. The reference itself has not moved either:
its scoring modules are unchanged since 2025-09-20, and the commits since are
nginx and docker-compose.

**Two things had drifted, both in the trim.** The docs still named
`dye_library.cif` in three places that describe the present tense
(`okf/labelizer-correspondence.md`, `prd-120.md` twice); fixed. And both
labelizer examples wrote their containers into the **current directory**, so
running them from the checkout left `1DDB-39.mmfdb.pto`, `dyes.mmfdb.pto`,
`MalE_apo_holo.mmfdb.pto` and two recovered PDBs sitting in the repo root,
where a derived file is indistinguishable from an input. They write into a
`tempfile.mkdtemp()` now, like `plot_fret_restrained_md.py` already did, and
the five leftovers are gone.

That is the same defect this port keeps meeting from a new angle -- a derived
artifact somewhere nothing checks it -- and the root of the tree still holds
about thirty more of them from other examples (`mGBP2_*.mrc/.txt`,
`A48_C1R_dry.prmtop`, `hgbp1_rotamer.fps.json`). Those belong to their authors,
but the cure is the same three lines.

Left alone deliberately: `data/dyes.mmfdb.pto` keeps its name though the CIF it
derives from was renamed, and `../labelizer-backend` has a staged deletion of
`terms.dic` that is not ours to resolve.

## 2026-09-01 (6) — the hard sphere becomes configurable, and a silent 147 A walk

Owner: *"make the hardsphere configurable... make it possible to scale such that
uses hard sphere."* `DockingParameters::clash_radii_source` (`"imp"` default,
`"olga"` available) and `clash_radii_scale` (default **1.0**, so nothing moves
unless asked). Source picks the table, scale multiplies it; both reuse
`VdwRadii.h` and the same spellings as `AV::set_radii_source`, so there is one
mechanism, not two. On `build_docking_assembly`, `score_structures`,
`dock_minimize`, `fps_bootstrap`, the CLI and `FPSProject` (schema **1.7** —
which also carries `coarse_clash`, because a project with the first two and not
the third resumes into an exception). **1771 pass.**

**No radii are written on the structure.** A substituted set rides *shadow
spheres* — new `XYZR` particles at the atoms' coordinates, members of the same
rigid bodies, not hierarchy leaves — so the atoms keep what `read_pdb` gave them
and an `AV` under `"imp"` still sees its own. FPS's D12 (mutating shared
molecules and never reverting) is not reproduced. Coarse beads plus a substituted
table is **refused**: a 2.5 A bead stands for a residue, a table keyed by atom
name has no entry for one, and a scale on it is `bead_radius` spelled twice.

**A silent catastrophic defect, found on the way, and it is IMP-level not
ours.** Adding rigid-body members after the model has computed dependencies once
leaves the body's position constraint reporting a **stale output list**; a
scoring function built afterwards then drops it from
`get_required_score_states()` entirely. Measured: 3 required states -> 1,
HIV-RT body-1 gradient exactly `(0,0,0)`, and `dock_minimize` walked the DNA
**147 A to the origin while reporting a better score** (31.69), because the
frozen proxies rode along. `RigidBody::add_member` calls `Model::clear_caches()`,
which does not cover it; `set_has_dependencies(false)` on the score states does.
Reproduced in 25 lines of plain IMP with no `IMP.bff` involved.

**One global scale is not defensible, measured rather than argued.** Bisected
against the HIV-RT interface, **0.8206** reproduces Olga's pair count, **0.8386**
its overlap, **0.8476** its energy — three answers, and at the first the energy
is off by 2x. IMP's ratio to Bondi is per element (C 0.8083, N 0.838, P 0.841,
O 0.889, S 0.900) and IMP's carbon alone spans 1.70-2.275 A where Bondi's is a
single 1.70. A scale matches one statistic of the contact distribution, never
the distribution. Selecting the source is the honest route; the scale is a
convenience and is documented as one.

**The bootstrap: resolved, not solved.** On HIV-RT `resolved`, clash goes from
**74 % of the score to 5 %** and the spread stops being identically zero:

| clash radii | parent | clash/score | RMSD spread |
|---|---|---|---|
| `imp` x1.00 (default) | 150.69 | 0.737 | **0.0000 +/- 0.0000** |
| `olga` x1.00 | 35.88 | 0.054 | **0.0689 +/- 0.0155** |
| `imp` x0.80 | 39.96 | 0.022 | 0.4040 +/- 0.0953 |
| `ev_weight = 0` | 24.34 | 0.000 | 2.1297 +/- 1.4584 |

0.069 A is still **thirty times** below the clash-free 2.13 A, and ~0.4 A only
arrives once clash is under 2 % of the score — well past any radii set FPS would
recognise. The radii artefact was real and is now controllable; **it was not the
whole story.** Whatever else pins this pose is not the size of an atom, and that
is the next thing to find rather than something to explain away.

Not done: `refine_docking` and `screen_structures` take no clash-radii arguments
(they take no `clash_tolerance` either) and run at the defaults.

## 2026-09-01 (5) — the default goes back to IMP's radii, and a test that could not see it

Owner, reversing entry (4): *"no use the IMP radii, charmm, as otherwise wont be
consistent with IMP docking."* Right, and for the reason now written into the
code in four places: the excluded-volume term is `clash_container`, which reads
`IMP::core::XYZR` — the **particles'** radii. A volume on Olga's table and a
clash term on IMP's are two halves of one score disagreeing about how big an
atom is. Internal consistency beats matching Olga's numbers.

Default `"olga"` -> **`"imp"`** (`AV_RADII_IMP = 0`), schema **1.5 -> 1.6**.
Olga's table stays vendored and selectable — it is what reproduces Olga-era
references. `"model"`, the previous spelling, is **refused rather than aliased**:
two names for one quantity is the failure being avoided, and that spelling was
one day old and had never left this repo. `"imp"` rather than `"charmm"` because
the value means *the radius on the particle*, which is CHARMM-derived only after
`read_pdb` — on a bead model `"charmm"` would be a lie.

**1741 pass**, exactly the baseline. Every pin returned to its pre-Olga value
except two AV mean positions, which land ~0.002 A away because the separate,
still-uncommitted attachment-atom change also moved them; those were re-pinned to
the measured values rather than forced.

**A test that could not detect the default it was running under.**
`test_fps_av_parity.py` writes FPS's Bondi radii onto the particles, so the whole
A/B is premised on the volume reading *those*. Under the Olga default that
substitution became **dead code and the test kept passing anyway** — because both
tables are Bondi-derived and carbon is 1.70 in each. It agreed for the wrong
reason. It now names its radii source explicitly and says why in the comment.
That is the shape to watch for: a test whose premise is silently satisfied by
something other than the code under test.

**The overturned finding survives the revert as a checkable pair.** Entry (4)
showed that "the authored clearance buries the structure" was a *radii* artefact,
not a clearance one. Reverting the default would have quietly re-established the
wrong conclusion, so the test now asserts **both** behaviours against their own
sources — default 9 authored / 24 derived, Olga 33 / 33 — and its docstring says
it is a statement about the two tables, not about which is default.

**The price of the default is pinned, not hidden.** With IMP's radii, nine of the
33 published <R_DA> have **no model value at all**, and the Zenodo agreement is
worse: **+0.223 / 0.908 / 0.99602** over 24 pairs against Olga's **-0.023 /
0.713 / 0.99536** over all 33. Against the FPS oracle, IMP **-0.315 / 1.204 /
0.99259** over 24; Olga +0.020 / 1.624 / 0.97772 over 33, or -0.616 / 1.131 /
0.99399 over the same 24. `fps_screening_ab.md` leads with the default, labels
Olga's as the alternative, states the cost in its own paragraph and tells anyone
reproducing Olga-era numbers to set `"radii_source": "olga"` — one sentence, so
the default does not read as free.

Consequence worth recording: this makes the HIV-RT bootstrap collapse a plain
**calibration** question rather than an inconsistency. FPS's `ClashTolerance`
constants were fitted against Bondi-scale radii; on IMP's larger ones they
over-penalise. A number to recalibrate, not a mismatch to reconcile.

## 2026-09-01 (4) — Olga's radii, vendored; and a finding of ours that was a radii artefact

Owner: *"the clash radii are different in fps, must use consistent set. use what
is in olga, include in bff and use that."* Done. `data/olga_vdw_radii.csv`
(derived; `src/VdwRadii.cpp` is the definition and `olga_vdw_radii_csv()`
regenerates it, the same contract as `data/fps_json_schema.json`),
`include/VdwRadii.h`, `AV::set_radii_source("olga"|"model")` **defaulting to
olga**, and a per-position `radii_source` field at schema **1.5**.
**1741 pass.**

**Three schemes were in play and they are not interchangeable**: FPS keys by
*element* (Bondi, Angstrom); IMP uses *united-atom* radii carrying implicit
hydrogens (carbon 1.85-2.275); Olga keys by **atom name** — 128 entries, `CA`,
`CB`, `CD1` each listed — with a flat fallback. Not a numeric swap, a different
lookup.

**The unit was established, not assumed.** `Olga/src/AV/Position.cpp:118-124`
scales the coordinate *and* the radius by the same `10.0f` on one line, because
pteros stores nm and `calculateAV` works in Angstrom. So C 0.17 nm =
**1.70 Angstrom** — Bondi's carbon, and FPS's, which is independent
corroboration. Pinned. A 10x error here would have been silent.

The fallback is **1.50 Angstrom flat, no element lookup** — within 0.01 of
Olga's own oxygen, so not transparent, but below every heavy atom, so it only
bites carbons. Misses on the shipped fixtures: T4L **none**; HIV-RT DNA 85 of
1018 atoms, all old spellings (`C1*`, `O1P/O2P`, `O3*`-`O5*`). Reproduced, not
corrected. Olga's JSON also overrides its own built-in P = 0.1 with 0.186, so
**1.86 is Olga's phosphorus**.

**Both A/Bs improve, and one widens.** Against the Zenodo table with the ACV on:
+0.22 / 0.91 / 0.9960 on 24 pairs becomes **-0.02 / 0.71 / 0.9954 on all 33** —
every published <R_DA> now has a model value where nine did not. Against the FPS
oracle the headline rmsd got *worse* (1.20 -> 1.62) **because the set grew**;
over the same 24 pairs Olga's radii win on both (rmsd 1.20 -> 1.13, r 0.9926 ->
0.9940). Both figures are asserted so neither can be quoted alone.

**A finding of ours was overturned.** Entry (10) recorded "the authored clearance
buries the structure" — 8 of 17 sites empty at `allowed_sphere_radius: 2`. That
was a **radii** artefact, not a clearance one: with Olga's table the same file
resolves all 33 pairs. The test's assertion is inverted and the old behaviour is
still asserted under `radii_source="model"`, so the attribution stays checkable.
Volumes grow accordingly on T4L — 132 CB +16 %, 55 CB +5 %, and **99 CB from 710
to 29 586 voxels, x42**, which is the site the whole empty-volume thread of
2026-08-31 was about.

**The HIV-RT bootstrap collapse is not resolved, and the reason is structural.**
Clash 110.98 of 150.69 -> 110.85 of 140.58; spread 0.000 +/- 0.000 either way.
The radii source is the **volume's**; the clash term is `clash_container`, which
reads `IMP::core::XYZR` — the *particles'* radii — so it never saw the change.
That the collapse *is* a radii artefact is nonetheless measurable: statically
over the protein-DNA interface at k = 8, united-atom radii give 268 overlapping
pairs / 90.5 Angstrom overlap / energy 198.4 against Olga's 30 / 7.6 / 11.1 — a
factor of **17.9**. Giving `clash_container` its own radii source is the
remaining half and moves every docking number, so it is a protocol decision and
was left for the owner.

Concurrency note: another agent was editing this tree at the same time, dropped
`include/Minimizer.h` with no `src/Minimizer.cpp` and broke the shared build,
and its `test/minimizer` segfaulted two full-suite runs. Both had resolved by the
verification run above.

## 2026-09-01 (3) — G9: the contact volume was real, and it was the +2 A

The last gap, and it **overturns entry (10)**. That entry concluded the +2 A
offset against the Zenodo table "is in the reference, not the code". Wrong — it
was the accessible-contact volume, silently inert.

**Wired, not invented.** The definition was authoritative and in reach: Olga's
`path2points()` (`../ucfret/thirdparty/olga/src/AV/fretAV.cpp:236`), the program
the shipped files were authored for. `PathMap::apply_contact_weighting` marks a
cloud voxel as *in contact* when the obstacle raster inflated by the **dye**
radius — Olga's `occupancyVdWDye`, the same array the carve uses — lies within
`thickness` of it, then scales those voxels to carry `trapped_fraction` of the
cloud's total weight. Two departures documented: shares rather than a fixed
ratio (identical for uniform base weights, i.e. Olga's own case; better defined
for AV3 and chain weighting), and no wrap-around at the grid faces. Olga's
whole-voxel **truncation** of the layer *is* reproduced, and the measurement
below is why. `space_fixed=False` refuses it and warns: that path exists to
reproduce pre-PRD-105 numbers byte-for-byte.

**The defect behind the defect.** `AV::set_av_parameter` read
`j.value("contact_volume_trapped_fraction", -1)` — an **`int`** default, so
nlohmann deduced `int` and **every fps.json trapped fraction was truncated to 0
or 1**. Invisible while the value was unused. The first round of measurements
was taken with it still in place and said the ACV *lengthened* <R_DA>; that was
the bug, not the feature. One character of type, and it inverted the sign of the
answer.

**Against the Zenodo table**, 3GUN, 24 pairs: ACV **off** +2.06 A / rmsd 2.54 /
r 0.9901; ACV **on** **+0.22 A / 0.91 / 0.9960**. Under the authored clearance,
+1.70 -> **+0.06 A**. So the file **is** the parameter set the table was computed
with, and the program was **Olga**, which has an ACV — not FPS, which has none.
The owner's "FPS has no ACV" was right; the inference drawn from it was not —
and the reason the two facts sit together is chronological: **the ACV did not
exist yet when FPS was written** (owner, 2026-09-01). It is an Olga-era
addition. The two references are therefore from two eras, and there is no single
answer to match: the FPS oracle is pre-ACV and the comparison against it drops
`contact_volume_*` as a key FPS never had, while the Zenodo table asks for an
ACV and honouring the ask is what closes it. A file gets an ACV if and only if it
asks for one.

The discretisation is decided empirically too: a true 3 A sphere gives
+0.49 / 1.04 where Olga's truncated shell gives +0.22 / 0.91. The truncation is
part of what produced the reference.

**The FPS oracle A/B is unchanged at -0.31 A / 1.20 A / r 0.9926**, because the
comparison now drops `contact_volume_*` for the same reason it drops
`allowed_sphere_radius`: FPS has no such key. Two references, two conventions,
and the port now matches each on its own terms.

On the shipped T4L file the ACV moves the mean 2.3 A on average and **closer to
the nearest atom at 16 of 17 sites**; mean <R_DA> 46.91 -> **43.88 A** (all 33
pairs shorter), the model-vs-experiment bias +3.93 -> **+0.84 A**, and the
network score 22.5277 -> **11.3268**. Five pinned numbers moved and were re-pinned
with the reason inline. **1698 pass**, one xfail flipped to a passing test.

Unverified and named: Olga was **read, not run** — there is no Olga oracle here
(it is a Qt GUI), so the parity rests on the source plus the agreement with the
table. That the table was produced **by Olga, with an ACV**, is no longer an
inference from the 24-pair fit — the owner stated it (2026-09-01), so the fit is
corroboration of a known provenance rather than the evidence for it. And AV3
classifies contact against
the first radius's raster where Olga uses each radius's own; no in-tree AV3 file
requests an ACV, so it is unmeasured.

## 2026-09-01 (2) — G1: a project, and a selection that was being silently mis-named

`FPSProject.h`/`.cpp` land: structure paths in body order, the labelling source
(fps.json **or** the legacy pair), the selected distances, the five per-mode
parameter blocks at FPS's shipped values, the conversion and AV globals, and the
poses — so `capture_poses`/`apply_poses` become store-and-continue. Both doors
per D3: the legacy `OptionsManager.cs::Export` **text twin** is read, the `.bin`
stays a non-goal, and the extended fps.json is read *and* written. Schema
**1.3 -> 1.4**, additive, with a new `FPS_OBJECT` field type for the nested
blocks. `imp_bff_fps project init/show/convert`, and `-P/--project` on every run
mode with `--save-project`/`--resume`. **1681 pass.**

**A silent mis-naming, found and fixed.** FPS's `SelectedDistances` is a
`Boolean[]` parallel to the **line order** of the legacy distances file, but
`read_old_distances_txt` returns a JSON object whose keys come out **sorted** —
so naming a positional selection through it pairs every flag with the wrong
distance. On HIV-RT the two `p66_K287C` entries at positions 8-9 became
`p51_K173C_p_19bp` and `p51_E194C_p_1bp`. That is a wrong answer with no
symptom: the run scores a different distance set than the file asked for.
`read_old_distances_order` now preserves file order, and without the distances
file the flags survive as unnamed `selected_flags` rather than being guessed at.

Three `Export` traps handled and commented, all of which produce a wrong read
rather than an error: array elements are space-separated and **unquoted on the
following line**, so a path containing a space is unrecoverable (FPS loses it
too); a **zero-length** array still writes that line, blank, so it must be
consumed unconditionally; and `Double.ToString()` follows the writer's culture,
so `0,0005` reads as **0** under `strtod` unless retried.

Kept deliberately separate: FPS's integrator settings are *stored and not
obeyed*, and `max_iterations` sits beside a distinct `iterations`, because
200 000 Verlet steps is not 500 conjugate-gradient iterations and silently
equating them would be the same category of error as the rest of this log.

**Measured, and not a bug — the size of an approximation.** A dock freezes its
volumes at the input pose and rides their means on the bodies; a **resume**
re-samples them in the pose it resumes from. HIV-RT, `-c resolved --shuffle 0
--seed 1 -n 30`: **41.8586** at the input pose, **30.6265** at the end,
**39.0207** at that same end pose with the volumes recomputed. The gap between
the last two is what the frozen-volume approximation is worth; `dock --help` now
says so with the numbers.

Also: every overridable CLI flag now defaults to `None`, because click cannot
otherwise tell "typed the default" from "typed nothing" — and that distinction
is exactly what decides whether a project's stored setting survives.

## 2026-09-01 — G2 and G5: the bootstrap, the tables, and three FPS defects made into choices

`fps_bootstrap` and `FPSExport.h`/`.cpp` land, with `bin/imp_bff_fps_export`
(`errors` / `table` / `screen`). **1661 pass.**

**The three FPS defects are implemented as *choices*, not guesses** — correct by
default, FPS-compatible pinned beside it, difference measured.

* **The perturbation.** `BFF_SPLIT_NORMAL` (two-piece, continuous, half-masses
  `sigma±/(sigma+ + sigma-)`) is the default; `FPS_SIGN_SPLIT_NORMAL` reproduces
  FPS. Over 400 000 draws at 3/10 A: mean **+2.7926** (FPS) against **+5.5852 A**
  (split, *exactly twice*), P(X>0) 0.500 against 0.769, density ratio at the
  target 0.30 against 1.00. **Worth knowing before choosing: fixing the
  discontinuity doubles the outward drift.** Neither is unbiased — an asymmetric
  density has its mean away from its mode. What the split normal buys is
  quantiles that mean what the error bars claim.
* **The zero-noise pins.** `perturb_deselected` defaults **true**, so every
  scored distance gets noise; false reproduces FPS's downward bias. On HIV-RT
  `resolved`, 20/0 against 18/2.
* **The RMSD sign error.** `pose_rmsd(..., fps_sign_convention)`, pinned to
  machine precision through the identity
  `rmsd_correct^2 - rmsd_fps^2 = (4/N) sum (U r_i) . t`. On HIV-RT the effect is
  ~1 % of the RMSD.
* **`BestFitRotation`** is not reproduced: the superposition is computed against
  the stated reference where it is used, so fitting a row onto itself gives the
  identity rather than FPS's zero-matrix 180 degrees.

**A silent failure found in our own code, not FPS's.** IMP's `create_rigid_body`
starts a body in its **principal-axis frame**, so `capture_poses` is *not* the
transform a PyMOL script needs — it has to be composed with the inverse of the
input frame. HIV-RT body 0 is a 55 degree rotation at the input pose. The
failure mode is why it matters: every structure rotates by the same amount, so
the overlay still *looks* right. `FPSMolecule::frame` and a test now pin it.

**And a physical one worth carrying into G7's defaults.** With FPS's
`ClashTolerance = 0.5` (k = 8) and IMP's united-atom radii, HIV-RT's
protein-DNA interface scores **102.9 of a parent score of 134.8 as clash**, and
the bootstrap spread collapses to **0.000 +/- 0.000 A** — the restraints cannot
move anything. At `ev_weight = 0` the same run gives **1.935 +/- 1.133 A**.
`parent_e_clash` is reported, the program warns, and it is pinned. This is the
radii question from 2026-08-31 (10) arriving with consequences: FPS's clash
constants assume FPS's Bondi radii.

Also carried: `chi2` raw with `chi2_r` beside it (FPS's `dof = max(N-6(M-1),1)`);
`chi2_bond` a subset; Filter `Number` 1-based with `--fps-filter-number` to
restore FPS's 0-based one; `_tmp.pdb` absolute and removed; `camera=0`; all *N*
transforms written; **one** R table rather than FPS's two, because this port's
model distance is simulated where FPS's second Dock-mode file is a polynomial
estimate wearing the name its Filter mode gives a simulated one. `results.json`
replaces the unreadable `.bin`.

Smaller finds: `-0.000` reached exported files and made identical exports compare
unequal; the unity build means anonymous namespaces do not isolate helpers, so
`open_for_write`/`stem_of`/`basename_of` collided with `StructureIO.cpp`.

**Not A/B'd against FPS.** No FPS reference output exists for error estimation or
for any export file, so all of this is pinned against the reference pages and
against closed forms — never against FPS output. PRD-121's open item 1 already
records that there is no FPS docked pose anywhere.

Left open and named: `estimate_errors` in `bin/imp_bff` measures optimiser
reproducibility, not measurement uncertainty; the two now point at each other but
it keeps its name. `fps_label_positions` reads a `dye` key no current fps.json
writes, so the pseudoatom colouring path is exercised but never coloured.

## 2026-08-31 (12) — fixed positions become scorable; bonds and the protocol land

**A position with no volume can now be scored.** `ProbeNetworkRestraint` sent
every position through `search_labeling_site`, and an `XYZ` position has no
chain/residue/atom, so the selection matched the whole structure and threw
*ambiguous* — a network containing one was **unbuildable**. Such a position is
now carried as a **point**. Two spellings: `XYZ`, a fixed coordinate, optionally
transported onto the structure being scored by a Kabsch fit of its
`reference_atoms` (which is where G3's `fit_reference_atoms` finally reaches
*scoring* rather than only reporting); and `ATOM`, a plain atom of the structure
— FPS's fifth LP type (`LabelingPositions.cs:203`), schema **1.3**.

FPS's rule (`FilterEngine.cs:305-307`) — a distance with a point end is scored
as **R_mp** whatever the file's `distance_type` — is applied to the
*measurement*, once, so the restraint, the transfer function, the CSV column and
the table agree by construction and the column reports what was scored rather
than what was asked for. Pinned: a fixed point written in a rotated frame comes
back to its atom to **2.2e-14 A**, and three files differing only in
`distance_type` give **identical** model distances while the same two positions
as *volumes* differ by more than 1 A — so the equality is the rule, not the
geometry.

**G4 (bonds).** Both ends `ATOM` makes a bond (`SpringEngine.cs:111-116`; an
`XYZ` end does not qualify, FPS requires `AtomID > 0`); both anchors drop to
0.4 A and leave clash detection; `e_bond` is reported as a **subset** of the
score, never added. FPS's D12 — the radius write leaking into shared `Molecule`
objects forever — is **not** reproduced, and a test proves a second assembly
sees the reader's radius. Measured: excluded volume exactly **0.0** with the
anchors dropped against **0.19531** without.

**G7 (protocol).** `optimize_selected`, `max_force`, `clash_tolerance` on
`DockingParameters`, with FPS's Huber form `chi2_score_capped`. Clashes are
**never** gated, as in FPS. The linear tail is verified as linear: slope
**400.0** to 1e-9 past `drmax`, at three different errors. `k_clash =
2/ClashTolerance^2` confirmed across four values (CT 0.5 is exactly 4x CT 1.0).
The scoring door keeps IMP's historical `k = 1` — its numbers are pinned — so
FPS's constant is opt-in through `DockingParameters`.

**A silent wrong answer, found and fixed here.**
`score_structures(...).pairs[0]` raised `IndexError` while
`r = score_structures(...); r.pairs[0]` worked: SWIG returns a bare
`std::vector` **member of a temporary** as a pointer into freed memory, and a
freed vector reports size **0** — so the expression read as *"this structure has
no distances"*. A wrong answer, not an error, and it looked intermittent because
binding to a name happened to work. Typing the member as `PairDistances` (the
`IMP_VALUES` vector, which IMP wraps properly) fixes it; the three assignment
sites construct from the iterators. Regression test in `test_docking_values.py`.
**1637 pass.**

Two more found and not fixed, both recorded rather than papered over:
`ProbeNetworkRestraint::get_model_distance` called **before the first
`evaluate()`** on a registry-backed network **segfaults** in
`AVOccupancyMap::atom_reach` — the shared registry's coordinate snapshot is only
taken in `begin_evaluation`, so this is PRD-105 registry-lifecycle work. And
`read_old_lps_txt` was silently dropping FPS's `ATOM` LP lines as an unknown
dialect, so a converted legacy file had **fewer restraints than the original**
— that one *is* fixed.

Still missing from G7: FPS's clash-free random shuffle (±5 A per axis, uniform
orientation, accumulating until no clash, molecule 0 fixed). `shuffle_bodies` is
still a single displacement with no clash test.

## 2026-08-31 (11) — G3 lands, and two FPS reference pages

**G3 (screening diagnostics), by subagent, verified here.** `ScreenedStructure`
grew `chi2_r`, `sigma1/2/3`, `invalid_r`, `ref_rmsd` and the per-structure pair
table; the sigma counting follows `FilterEngine.CalculateChi2` exactly
(`error_pos` when the model is long, `error_neg` when short, counts nesting).
New value types `ReferenceAtom`/`ReferenceFit` and the kernels
`fps_reference_atoms`, `fit_reference_atoms`, `fit_reference_positions`,
`reference_rmsd` — the Kabsch fit **reuses**
`IMP::algebra::get_transformation_aligning_first_to_second` rather than becoming
a second copy of the same SVD-with-determinant-flip. `reference_atoms` is in the
schema at **1.2**, and `imp_bff_fps screen` reports the columns and gained
`--pairs-csv`. **1595 pass** (was 1547).

Measured on T4L `chi2_C2_33p`: `3GUN_faspr_port.pdb` 22.3165 / chi2_r 1.3525 /
sigma 12-3-0; `3GUN.pdb` 22.5277 / 1.3653 / 10-3-0. `2*score == chi2_r*33` on
both, because `score_model` is half the chi-square — now documented on the field
rather than left to be rediscovered.

**The gap that keeps it from being finished.** `ProbeNetworkRestraint` cannot
build a network containing an `XYZ` position at all:
`create_av_decorated_particles` calls `search_labeling_site` for every used
position, and an `XYZ` position has no chain/residue/atom, so the selection
matches everything and throws *"ambiguous labelling site"*. The transported
coordinate is therefore computed and reported but **has nowhere to be scored** —
a distance touching a fixed position still cannot enter chi-square. FPS's rule is
that either end of type `None` makes the model distance R_mp; wiring that belongs
with the scoring core. Also open: `read_old_lps_txt` drops the `ATOM` lines
trailing an LP line (FPS's `ReadRefAtoms`), so reference frames arrive only
through `.json` — and no in-tree fixture contains a reference atom or an `XYZ`
position at all.

Pre-existing and confirmed, not introduced: `data/fps_json_schema.json` **at
HEAD is at `x-schema-version` 1.0**, stale against the authored tables by the
whole 1.1 step. The drift test would fail at HEAD. The working tree already
carried the fix.

**Two FPS reference pages**, from read-only research agents, in
`okf/references/`. They cite `File.cs:line` throughout and each carries a defect
index; reading them is cheaper than rediscovering what they found.

* [`fps-export-formats.md`](references/fps-export-formats.md) — the six outputs,
  byte by byte. **`SimulationResult.RMSD` has a sign error**: it accumulates
  `|U*r - t|^2` where the displacement is `U*r + t`, so every RMSD FPS has ever
  printed is wrong unless the two translations coincide. And `BestFitRotation`
  is **never computed by any `Save*` method** — it is a side effect of the GUI's
  RMSD refresh mutating the live array, so an exported "Overlay" is a chain of
  *pairwise* fits, not a common superposition, and a zero matrix (spurious 180
  degree rotation) if best-fit was off.
* [`fps-sampling-and-protocol.md`](references/fps-sampling-and-protocol.md) —
  error estimation is a parametric bootstrap on a **sign-split** normal (each
  half mass 1/2 regardless of width, so the density jumps by `err-/err+` at the
  target and the perturbation carries a mean shift of
  `(err+ - err-)/sqrt(2*pi)` — +2.79 A for 10/3 errors). `rkT` is `1/kT`, so
  larger is **colder**; the default 10 accepts a 1-chi-square uphill move with
  probability 4.5e-5. `Molecule.Selected` means **not** randomised, the opposite
  of what the UI comment says.

**And the answer to a question this PRD had open.** FPS computes docking AVs on
**isolated subunits** and `Refinement.RedoAV` rebuilds them **inside the docked
complex** (`new Molecule(sr)` merges every subunit in its docked pose, the AV is
computed there, and the mean position is mapped back to the home frame). So a
partner subunit does occlude a dye, and refinement is the step that accounts for
it — a first-order self-consistent loop, not iterated. That is the design to
follow, and it explains `p66_K287C` in the HIV-RT example.

Also landed: `bin/imp_bff_fps_av` and `bin/imp_bff_fps_distance` (G6), FPS's AV
dialog and distance calculator, with FPS's own dye presets and grid rule. The
distance tool reads FPS's duplicate-expanded `.xyz` and reports the unique voxel
count separately, so 5473 lines are never mistaken for 5473 voxels.

## 2026-08-31 (10) — the +2 A was the reference, not the code

Owner: *"FPS has no ACV. find issue with 2Ang."* Both halves right, and the
second one took the oracle to settle.

**Ask FPS.** `prototypes/fps_oracle/` was built to reproduce a shipped cloud; it
also *generates* one, so it can be pointed at the same T4L structure with the
same file and asked what FPS computes. Over the 33 pairs:

| | bias | rmsd | r |
|---|---|---|---|
| **this module vs FPS's own routine** | **-0.31 A** | **1.20 A** | **0.9926** |
| this module vs the published table | +2.05 A | 2.54 A | 0.9901 |
| **FPS's own routine vs the published table** | **+2.37 A** | 2.73 A | 0.9735 |

The two implementations agree. **Neither reproduces the table, and FPS misses
it by more than this module does.** So there was never a bias in `IMP.bff` to
explain, and the ACV story in entry (9) was wrong twice over: FPS has no ACV,
and the premise it was explaining did not exist.

**Not a mislabelled distance type** either -- on FPS's own clouds <R_DA> is
+2.03, <R_DA>_E +3.11, R_mp -1.40 against the table, so no statistic fits. **Not
AV1-vs-AV3**: recomputing with FPS's own AV3 radii from `Fps/data/linker.txt`
(donor 5/4.5/1.5, acceptor 11/3/1.5) makes it *worse*, +6.36 A.

What is left is the file. `FRET_screening.fps.json` is **not the parameter set
the Zenodo table was computed with**, and says so on its face: seven of its
seventeen sites carry a linker width of **2.5 or 3.5 A** where FPS's own
`linker.txt` gives **4.5 A** for both dyes, and it carries
`allowed_sphere_radius`, `contact_volume_thickness` and
`contact_volume_trapped_fraction`, three keys FPS has no concept of. It is an
`IMP.bff`/Olga file sitting beside an FPS table; pairing them measures the
pairing.

**The end-to-end pin moved accordingly.** It is now the oracle
(`test/references/fps_screening_oracle_pins.json` -- FPS, reproducible,
parameter-matched by construction), and the published table stays in the suite
only as the record of a mismatch that would otherwise be rediscovered by the
next person. [`okf/validation/fps_screening_ab.md`](validation/fps_screening_ab.md)
rewritten; G9 keeps the defect and loses the causal claim.

The lesson is the one this log keeps relearning in a new costume: **an
agreement measured against the wrong reference is not a measurement.** Three
times now on this port -- a `//4`, a duplicate-expanded cloud, and a file that
was never the table's.

## 2026-08-31 (9) — PRD-121 phase 0 closed: two A/B pins, and the contact volume that does nothing

**The single-volume pin.** `test/input/fps/` gets the FPS frame and
`p_1bp_D.xyz` -- FPS's own cloud, verbatim -- with `p66`'s numbers in
`test/references/fps_av_pins.json` rather than its 3 MB cloud, because
`prototypes/fps_oracle/` regenerates that one exactly. 1.3 MB, not 53.
`test_fps_av_parity.py` pins it: `p_1bp` at **0.985** of FPS's volume and
**0.121 A** of its mean position, **99.04 %** of this module's voxels within
one voxel of an FPS voxel and **100 %** within two, furthest 1.130 A -- and the
reverse, 99.0 % of FPS's within one of ours. The two are one volume on lattices
that do not line up.

Two traps the test now states rather than leaving for the next reader: FPS's
`.xyz` line count is the **sum of densities**, not the voxel count (5473 lines
over 3187 voxels), and its `Dmp` is the **density-weighted** mean. Getting
either wrong makes agreement look like disagreement -- it did, twice, earlier
in this log. A hand-rolled neighbour search in the first cut of the test made
the same class of error a third time: it searched a fixed +-1 cell
neighbourhood at a two-cell tolerance and under-reported coverage as 99.7 %
when it is 100 %.

**The end-to-end pin, and what it found.** Zenodo 3376527 holds what the real
FPS produced for 421 T4L structures. One of them is **3GUN**, which this module
already ships -- so the whole path is measurable against FPS's own answer with
**no new structure data at all**. `test_fps_screening_ab.py`;
[`okf/validation/fps_screening_ab.md`](validation/fps_screening_ab.md).

* **The authored clearance buries the structure.** `allowed_sphere_radius: 2`
  leaves **8 of 17** sites empty and only 9 of 33 pairs with a model value;
  deriving it leaves **one** and gives 24. That field is an `IMP.bff` key
  somebody chose, not an FPS parameter -- FPS seeds from
  `LinkerInitialSphere x linker_width` and has no such concept, and its seed is
  *unconditional* where this module's is carved out of an already-inflated
  obstacle map. Any fps.json with a small explicit clearance is asking for
  something other than it appears to.
* **A systematic +2.05 A.** Not scatter: r = **0.9901** over 24 pairs, rmsd
  2.54 A, and this module's <R_DA> consistently *longer* than FPS's.
* **The contact volume is inert.** `contact_volume_thickness: 3`,
  `trapped_fraction: 0.55`: turning it off changes **nothing** -- on T4L
  132/CB the volume, the per-voxel weights and the mean position are
  bit-identical. Accepted and ignored, and pinned as a **strict xfail**. It was
  first written up here as the likely cause of the +2 A; entry (10) shows that
  was wrong.

Phase 0 is closed. 1546 pass, 4 xfail.

## 2026-08-31 (8) — `imp_bff_fps`, a docking example, and a segfault the example found

**A program for the FPS run modes.** `bin/imp_bff_fps` with `score`, `dock`,
`refine`, `screen` and `convert`. Every mode takes either an fps.json or the
**legacy C# FPS pair** directly (`--positions` / `--distances`), so a decade of
files runs without a conversion step; `convert` exists only for when the
fps.json should be kept. Each mode's `--help` carries worked examples on
shipped data rather than placeholders. `README.md` gains the section that
`setup_module.py` requires of any new `bin/` program.

**Example data that is actually a docking problem.** T4L is one rigid body and
mGBP2 has a single distance, so neither demonstrates docking. Shipped
`examples/structure/HIV_RT/` instead: HIV-1 reverse transcriptase (p66/p51,
1R0A) with its DNA primer/template -- FPS's own docking test case, two bodies,
11 positions, 20 distances -- with the original `LabelingPositions.txt` and
`Distances.txt` beside the converted `hiv_rt.fps.json`, so the legacy readers
have a shipped example too. `examples/structure/fret_docking.py` walks
convert -> score -> dock -> refine -> screen on it. Docking from a random start
reaches chi2 38.4, 31.0 and 27.5 on three seeds against 59.0 for the deposited
complex, which is also the honest lesson the example makes: one run is one
local minimum.

Two things that file is honest about rather than hiding. Its `all` score set
carries all 20 distances and scores **inf**, because the two involving
`p66_K287C` have no model value -- that site is buried at the protein-DNA
interface and its volume is **empty**; `resolved` is the other 18. And the
reason it is buried here at all: **FPS builds each volume on its own subunit in
isolation, `IMP.bff` builds them in the assembled complex.** A modelling
difference worth a PRD-121 open item, because it decides whether a docking
score is separable from the pose being searched.

**Running the example segfaulted, and the cause was worth the trip.**
`search_labeling_site` read `p_residue[0]` and *then* asked whether
`p_residue` held exactly one element -- and `IMP_USAGE_CHECK` compiles out of a
release build, so a structure that does not contain the chain a position names
indexed element zero of an empty vector and the process died. That is the
**ordinary** case for a screen: candidates differ, and some do not carry every
site. `screen_structures` already caught `IMP::Exception` to record such a
structure with a NaN; there was simply no exception to catch. Now there is one,
naming the site. Two regression tests in `test_docking_values.py`.

Found in the same function and removed: a bare
`std::clog << p_chain << p_residue << p_atom` that printed
`["Chain B", ...]["GLU"]["Atom CB of residue 194"]` to stderr for **every
labelling site of every evaluation** -- the noise that has been burying the
output of every run in this log. It is an `IMP_LOG_VERBOSE` now.

1535 pass.

## 2026-08-31 (7) — one mechanism for the attachment atom, and the invariant that caught two

Owner: *"do what FPS does with the attachment atom."* Two things came out of
making that literally true.

**The subtraction is a true exclusion, measured.** Running the same volume two
ways -- the attachment atom present with `drop_source_obstruction()` doing the
work, and the attachment atom given radius zero so it is genuinely not an
obstacle -- gives **identical point sets** on both FPS reference sites, 3331
and 50 920 voxels, zero difference either way. So the count subtraction is not
an approximation of FPS's `if (i == atom_i) continue;`, it is that.

**Trying to do it in two places broke an invariant, and the invariant said so.**
The first attempt at extending the drop beyond the default path made the anchor
transparent in the *obstacle radii override* as well. Both mechanisms are
correct alone; together they disagree, because the shared raster cannot honour
a per-volume override (it is one obstacle set for every volume in its class)
while a private raster is built from exactly that override. `shared` and
`private` then computed different volumes for the same position -- 141.02
against 143.09 on T4L frame 0, 14 of 17 AVs differing --  and
`TestTier1Exactness::test_shared_equals_private` failed. It failed only on a
**moving** structure; on a static one the two agreed, which is what a test over
five trajectory frames is for. Collapsed back to one mechanism: the anchor is a
normal particle in the override, and `drop_source_obstruction()` runs on every
lattice path, shared or private. shared == private on all five frames again.

**`space_fixed=False` is deliberately left obstructing itself.** It is the
pre-PRD-105 opt-out and `references/prd105_legacy_pins.json` exists to prove it
still reproduces the old numbers; giving it new physics would remove the only
thing it is for. The consequence is stated where the two paths part: a volume
computed with `space_fixed=False` is smaller than the same volume on the
default path. `p66` is 22 182 voxels there against 50 920 on the default path.

Parity against FPS's own cloud, final for this pass -- FPS's atoms and radii,
`p_1bp`: **0.99** of its voxel count at **0.121 A**. `p66`: **0.92** on heavy
atoms, and still empty on explicit hydrogens, where FPS's three-voxel link hop
tunnels through channels narrower than a voxel and this module's search will
not. 1520 pass.

Noted in passing: `test/expression/test_expression.py`'s four "compiled once"
tests flaked once in a full run and passed alone and on a re-run. They count
compilations across a shared plan cache, so they are sensitive to state left by
another run. Unrelated to this work, but they are a false-alarm generator.

## 2026-08-31 (6) — the attachment atom does not obstruct its own linker

Owner, on reading the parity table: *"you must drop the attachment atom,
otherwise it won't work."* Correct, and it was the last structural difference
between this module's AV and FPS's. FPS drops it before it rasterises anything
(`av_routines.cpp:50`, `if (i == atom_i) continue;`); `IMP.bff` kept it, so
every volume grew out of a source sitting inside its own inflated sphere, and
the `allowed_sphere_radius` had been standing in for the missing rule.

**How, without losing PRD-105.** The obvious implementation -- give the source
radius zero in the obstacle override -- costs the shared occupancy raster,
because that raster is one obstacle set for every volume in its (spacing,
extra-radius) class and each volume drops a *different* atom. But the raster
stores a per-voxel **atom count**, so one atom's contribution subtracts
exactly: `drop_source_obstruction()` decrements the window inside
`r_source + extra`, a voxel only that atom covered falls to zero and opens, a
voxel any other atom covers stays blocked. Sharing is untouched. Applied to
the linker pass and to every dye radius of the carve, with the source's radius
read in `prepare()` because the compute phase may not touch the Model.

**Measured against FPS's own cloud** (`p_1bp`, its parameters, its radii):
0.88 -> **0.99** of FPS's voxel count, mean position 0.179 -> **0.121 A**.
`p66`, which was **empty**, computes at **0.92** on heavy atoms. On explicit
hydrogens `p66` is still empty: FPS's link search hops `linknodes = 3` voxels
and tunnels through channels narrower than a voxel, which is exactly the leak
PRD-105 deliberately closed. That difference is real, is FPS's, and is not
being reproduced.

**A bug in the first cut, worth recording.** The subtraction fired even when
the attachment atom was not in the obstacle set at all -- the array door builds
its obstacles from a caller's list and the source is not among them -- so it
opened voxels no rule had opened. `test_av3_matches_labellib_rule` caught it:
two voxels out of 10 683, on the level-set property rather than the
mean-of-indicators one. The guard is `source_radius <= 0` means nothing to
subtract.

**And a second copy of the phantom weight**, found while editing next to it:
`st.last_mean` on the lattice fast path had its own `sum = 2.0` -- the source
in the numerator once, in the denominator twice -- and *that* is the value
`resample()` writes into the AV's coordinates. Yesterday's fix to
`get_mean_position()` had left it. Two places, one quantity, exactly what
AGENTS.md is about.

**Blast radius, all repinned deliberately:** the T4L quadrature score
22.5933 -> 22.5277 (`test_docking_values`, `test_ProbeNetworkRestraint`), the
stencil-30/26 values 22.1322 -> 22.1203 and 22.7055 -> 22.6933, two AV mean
pins, and the Olga A/B decay tolerance 1e-2 -> 1.5e-2 -- where the pair
**selection and its order are unchanged** and only the RMSD the selection
reaches moved, by 0.3 %. Two fixtures that deliberately built an *inaccessible*
volume at T4L residue 99 stopped being empty (437 voxels) and were moved to
linker width 4.5, where the site is inaccessible at every clearance from 1 to
6 A. 1520 pass.

## 2026-08-31 (5) — PRD-121 phase 0: two AV defects, and a pin that had recorded one

Both found by measuring against the rebuilt FPS AV, and both produced a
plausible-looking number for a volume that did not exist.

**One clearance derivation, not two.** `max(1.5, 0.5*linker_width + 0.5*grid)`
lived in `compute_av_from_structure()` alone, so the two doors onto the same
volume disagreed about a parameter neither caller passes: the fps.json door
derived it, the **decorator** door took a flat 1.5 and returned **nothing** at
FPS's standard linker width of 4.5 A. It reached users --
`imp_bff av-export -c A -r 132 -a CB --linker-width 4.5` wrote an `.xyz` whose
first line was `0`, printed a mean position, and exited 0. The rule is now
`AV::get_effective_allowed_sphere_radius()`, negative means *derive*, and
`AVBuilder` delegates rather than keeping its own copy. `data/fps_json_schema.json`
said the default was `1.5`, which was never the value an omitted key took --
now the sentinel, because "no default" in that table means *required*.

**A phantom unit of weight.** `AV::get_mean_position()` started its weight sum
at `1.0` while the source contributed `1.0` more, so the numerator held the
source once and the denominator counted it twice. Every mean position was
pulled toward the origin by `(1+W)/(2+W)` -- negligible on a big cloud,
**0.175 A** on T4L residue 99 at 81 points, which is exactly where a
constrained site needs it most -- and an **empty** volume returned precisely
half the source coordinate. That is where the `source/2` signature came from.
An empty volume now reports its anchor, `get_mean_position(False)` returns NaN
when there is nothing to average, and a resample that finds no voxel **warns**,
naming the clearance, because in every case seen so far that is the cause.

**The pins had recorded it.** Four of the twelve `prd105_legacy_pins.json`
cases are empty volumes (`n == 0`), and each pinned "mean" was half its source
atom -- a fiction asserted as reference data. The pin file now carries an
`empty` flag and the test asserts emptiness explicitly for those cases instead
of comparing a centre. PRD-105's guarantee is untouched: the *maps* are still
byte-identical, only the derived mean moved, which is what the fix is.
`restraint_mp` (33) and `traj_mp` (165) regenerated for the same reason.

**What is deliberately not fixed.** `IMP.bff` keeps the attachment atom in the
obstacle set; FPS drops it (`av_routines.cpp:50`). So the derived clearance
clears the *inflation* but not the source atom's **own** radius, and at a fine
grid the source is still walled in: T4L 132/55/19/86 all come back empty at
width 4.5 / grid 0.5, and all are non-empty under
`r_source + width/2 + grid/2`. That rule would also grow volumes that already
compute -- T4L 132 at width 0.5 goes 66 805 -> 85 456 points -- so it is a
change to a physics default across the module, not a bug fix, and it is the
owner's. `test/representation/test_av_source_clearance.py` pins today's answer
so that changing it is deliberate; 1520 pass.

Measured on the way, and worth carrying into the parity table: reading FPS's
atoms and radii, the `p_1bp` volume is **0.99** of FPS's at 0.121 A once the
attachment atom is made transparent as FPS does, against 0.88 with it opaque.

## 2026-08-31 (4) — FPS surveyed, its AV rebuilt and run, and PRD-121 authored

Read the C# FPS toolkit (`../chisurf/junk/fps`, `Fluorescence-Tools/fps` at
`eb3489f`) against this module to plan a port. The survey is
[`okf/prds/prd-121.md`](prds/prd-121.md); the summary is that the docking half
is mostly here already and the **product** is not.

Already covered: AV1/AV3 and mean positions, R_mp/⟨R_DA⟩/⟨R_DA⟩_E/σ_DA,
asymmetric ±χ², both legacy `.txt` readers, rigid-body assembly with excluded
volume, `dock_minimize` (whose own docstring calls itself "the FPS approach"),
`refine_docking`, `screen_structures`, pose capture/resume, and the polynomial
distance conversion. Three of those — `screen_structures`, `refine_docking`,
`score_structures` — exist in C++ with **no door**: no CLI reaches them.

Seven gaps, two of which carry the PRD:

* **Error estimation is not what `imp_bff dock-errors` measures.** FPS takes
  the *model* distances as truth, perturbs each by a Gaussian scaled by that
  distance's own asymmetric ±error, and re-docks — a parametric bootstrap of
  the coordinate uncertainty. `dock-errors` re-docks from random starts, which
  measures optimiser reproducibility. Different quantity, same name.
* **`screen_structures` returns `(path, score)`.** FPS also returns 1/2/3σ
  violation counts, `InvalidR`, the per-structure R table in both R_mp and the
  file's own distance type, and RefRMSD — a Kabsch fit of per-position
  *reference atoms*, which is how an `XYZ` position is transplanted onto each
  library structure.

The other five: no project object; no per-repetition results table and none of
the exports (PyMOL, overlays, Rtable, chi2table); no bond restraints (in FPS a
distance between two plain atoms drops both vdW radii to 0.4 Å and reports
`Ebond` separately); no doors for the three tool dialogs; and the protocol
details that change the answer (`OptimizeSelected` Selected/All/SelectedThenAll,
the clash-free shuffle, the force cap).

Found while reading, worth recording: `AVS` (real dye structure) and `EDF`
(external density file) are **enum values FPS never implemented** —
`LabelingPositions.cs` parses them, `AVEngine` has only `Calculate1R` and
`Calculate3R`. They are not being ported because there is nothing to port. And
the legacy project `.bin` is a gzipped .NET `BinaryFormatter` payload
(`OptionsManager.cs:287`), but `om.Export` writes a readable text twin beside
it on **every** save — so the project is readable without implementing MS-NRBF.

**The fixture is FPS's own test bundle**, and testing it found a bug.
`../chisurf/examples/4w_junction/fps_test_data/` is FPS's shipped test data,
2014-05-19, still carrying its author's `C:\Users\doroshen\Desktop\` paths.
The directory name is wrong — the contents are HIV-1 reverse transcriptase
(p66/p51, 1R0A) with its DNA primer/template, not a four-way junction. It
carries docking inputs, 20 screening frames, and — the find — **two AV point
clouds FPS itself exported**, `p66(D).xyz` (130 531 points) and `p_1bp(D).xyz`
(5 473 points), each with its mean position and its parameters in the header.
Those parameters *disagree with the LP file beside them*: they were typed into
the AV interface by hand, so the header is the truth and the `.txt` is not.
The grid follows from them with nothing guessed — FPS's
`dg = max(min(0.2·L, 0.2·W, 0.4·Rᵢ), 0.4)` gives 0.6 Å for both, and the grid
extent `IMP.bff` builds matches exactly.

Trying the pin found **an AV that fails silently**. When the clearance is too
small the volume comes back empty, `get_is_valid()` still answers `true`,
nothing is logged, and `get_mean_position()` returns **exactly half the source
atom's coordinate**. It reaches users:
`imp_bff av-export -p 3GUN.pdb -c A -r 132 -a CB --linker-width 4.5` writes an
`.xyz` whose first line is `0`, prints a mean position of `(0.68, -8.60, -2.73)`
— half of CB at `(1.353, -17.196, -5.469)` — and exits 0; at the default width
1.5 the same command writes 2734 points. **4.5 Å is FPS's standard linker
width**, used by every position in both fixtures, so this is not an exotic
input. With `linker_width = 4.5` and `dg = 0.6` the derived
`allowed_sphere_radius` is 1.5 and the volume is empty; 2.5 is still empty; 3.0
produces points.

**Then FPS's own AV was rebuilt and run**, which is the only way this was
going to be settled. `Fps.Native/av_routines.cpp` is x86 SSE and this machine
is arm64, so it was compiled `-arch x86_64` and run under Rosetta behind a
driver replicating `Molecule.cs` (PDB parsing), `StaticData.cs` (the element
ladder, the Bondi table) and `AVEngine.Calculate3R` (the grid). It reproduces
both shipped clouds **exactly** — `p_1bp` 3187 unique voxels and `p66` 55 514,
identical point sets, zero difference either way. It lives in
`prototypes/fps_oracle/`; FPS is LGPL-2.1 and is **not vendored**, the script
compiles it in place from the sibling checkout. That is an **oracle**, not a fixture: FPS's AV can
now be run on any structure, site and parameters, so Phase 0 is not limited to
the two clouds that happen to be on disk.

**The apparent disagreement was a reading error, twice over.** FPS's `.xyz` is
a **duplicate-expanded, density-weighted cloud** and its line count is not a
voxel count: `calculate3R` returns `n += dn`, the *sum* of densities — one per
(voxel, dye radius) that fits — not the number of voxels
(`av_routines.cpp:414`), and `AVEngine` sizes the point array by that sum and
emits a voxel once per radius that fits. `p_1bp(D).xyz` is 5473 lines over
**3187** unique voxels (901 once, 2286 twice — R₁ = 11 Å never fits);
`p66(D).xyz` is 130 531 lines over **55 514**. FPS's `Dmp` is the mean of all
5473 lines, i.e. the density-weighted mean, reproduced here to the last digit.
`IMP.bff` carries the same information as unique voxels plus a weight column
(`0.3333` = one radius of three) — same content, different packing, which is
precisely the marshalling difference AGENTS.md is about.

**Parity, measured properly.** `p_1bp` at L = 8.5, W = 4.5, R = 11/3/1.5,
dg = 0.6, clearance 3.0, against FPS's 3187 voxels: heavy atoms with IMP's
united-atom radii give 2816 (0.88) at 0.324 Å; heavy atoms with FPS's Bondi
radii give 3178 (1.00) at 0.222 Å; **all 17 733 atoms with Bondi radii — what
FPS actually reads — give 2991 (0.94) at 0.104 Å**. And they are the same
volume, not merely the same size: **100 %** of bff's voxels lie within one
voxel of an FPS voxel, 99.8 % within 0.52 Å (the largest offset two
same-spacing lattices can have), 95.3 % within half a voxel. bff finds no
volume FPS does not; it is conservative by 6 %. An exact voxel-set comparison
reports **zero** overlap however well the two agree, because bff anchors on the
global absolute lattice (PRD-105) and FPS on the source atom — the wrong metric
to reach for.

What is left of G8 is one documented knob: FPS's `LinkerInitialSphere` is 0.5
and `allowed_sphere_radius` is not the same quantity, and the mapping has never
been written down.

**One trap, recorded so it is not repeated:** `AV::get_map().get_xyz_density()`
returns **one entry per point**, not a flat `x, y, z, w` array. Dividing its
length by four under-reports the volume by exactly 4×, and was the first of the
two reading errors above. `IMP.bff.write_av()` is the honest count.

Owner decisions (2026-08-31): IMP backend, so **scores are pinned and
coordinates explicitly are not**; five programs in `bin/` mirroring FPS's own
split; both project formats in (legacy text export, extended `.fps.json`) and
only fps.json out; both input doors everywhere (the legacy positions+distances
pair, or one fps.json), which puts `reference_atoms` into the schema at 1.2.

Two open items closed the same day (owner): there is **no FPS docked
reference anywhere**, so docking is pinned self-consistently and against the
crystal structure, never against FPS — D1 had given up coordinate parity and
this makes it unconditional; and the fixtures ship as a **reduced set in
`test/input/` with `.bcif` permitted**, except the two FPS `.xyz` clouds, which
ship verbatim because re-encoding a reference answer is how a pin stops being
one.

The second (end-to-end) screening pin was already in the tree and unused:
`prototypes/fast_label_score/cache/anisotropy/zenodo3376527_FRET_screening/` —
421 T4L structures × 33 pairs of ⟨R_DA⟩ plus χ² for three states, produced by
the real FPS, with its own `FRET_screening.fps.json` and experimental distances
beside it. Phase 0 is G8 first — localise the volume deficit on the two clouds,
where one site gives one answer — then wire the Zenodo set, where 421 answers
say whether the fix generalises. No engine work before an AV matches.

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

## 2026-08-30 — three new dyes parameterized the AmberTools way; exhaustive queue running

The original 12-library queue is **finished** (T51_C3R completed at 2 fs,
723 ns/day — the C16/water 4 fs diagnosis was right); all trajectories are
on cordeshub awaiting clustering. The three high-priority dyes — **ATTO 655
(T65), Cy3B (C3B), Lumiprobe AF488 (AF4)**, each on C2R/L1R/B1R — went
through `prototypes/dye_library/parameterize_dye.py`: the standard
AmberTools route applied to the **whole label** (ACE–linker–NME–dye as ONE
molecule): antechamber AM1-BCC/GAFF2 → parmchk2 → `loadmol2`. One fit, one
frcmod, no dye/linker junction to splice. All nine labels validate (charges
±0.005 of the conjugated net: 0/0/−2) and produce `labels/<stem>.mol2 +
.frcmod` in the repo and on cordeshub. `make_library.build()` routes through
the label automatically when one exists, seeds inject by (residue, name).

Hard-won specifics, all recorded in the script: no AMBER-format
parameterization webserver covers these dyes (CHARMM-DYES has Cy3B in
GROMACS/CHARMM36 with its own linker — wrong force field family, kept only
as a structure cross-check); sqm needs mol2 input with bond orders (PDB
bond-perception fails on the conjugated dyes) and a clash-free geometry
(the dye is Kabsch-placed onto a shipped template's C99/O99/d1 frame, then
MMFF-relaxed via RDKit); RDKit needs the formal charges derived
structurally (bond-order-sum 4 on N = the cation; one O⁻ per
delocalized S/C group) and carried to the label by coordinate match; and
tleap reads mol2 fixed-width — **four-character residue names from mol2
substructures lose their first character in the prmtop** ("Cy3b"→"y3b"),
hence the three-character codes. **Exhaustive libraries** per decision:
27 runs (9 stems × 3 donor-rank replicas × 140 ns) now grinding serially on
cordeshub (`run_queue_newdyes.sh`), then `cluster_dihedral.py --cutoff 1`
over the concatenated replicas.

## 2026-08-28 — cordeshub dye-library queue: L2R batch unblocked, T51_C3R NaN outstanding

State of the `prototypes/dye_library` MD queue on `cordeshub` (resumed from the
"rotamer" session): the overnight pass finished 8 of 12 seeded libraries
(C5W/C3W/C5N/C7N/C3N `_L1R`, C5N_C2R, T51_C3R — see below, A64_B1R); three
`L2R` dyes (T51/T48/Tth) died in tleap with
`1-4: cannot add bond 33 75`. Root cause: `compact()` writes the compact PDB
with OpenMM's `PDBFile.writeFile`, which emits **CONECT records**, and the
solvated `build()` then re-adds the N99–C99 bond explicitly — tleap refuses a
duplicate. Every finished library was seed-grafted (tleap-written compact PDB,
no CONECT), so the compaction path had never actually reached the solvated
build before. Fixed in `compact()` by dropping CONECT (lossless: all other
bonds come from residue templates), applied identically to the local
`prototypes/dye_library/make_library.py` and the cordeshub copy (backup at
`make_library.py.bak-20260828`); smoke-tested end-to-end on T51_L2R at
`--ns 0.02`. Also generated the missing `A53_C1R_seed.pdb` locally via
`cheap_library.py A48 C1R --donor A53_C1R` (the `A48_C1R_gate` entry had in
fact already completed on 08-26; its log was below my first listing's fold).
Queue relaunched detached (`~/dye_library/queue.log`): the skip guard re-runs
only the unfinished — the three L2R dyes, ~4 h each serial on the RTX 4000 Ada.

**Outstanding:** T51_C3R's 01:57 "finish" was itself a NaN —
`Particle coordinate is NaN` at the *first 4 fs block* of `produce()`'s ramp
(minimize + 160 ps at 1–2 fs survive; the queue retried it and it reproduced
deterministically, same velocity seed). **Diagnosed same day**, by replaying
the ramp instrumented and stepping 4 fs one at a time: the ramp's PE rise
(−198.6 → −159.2 MJ/mol) is exculpatory — controls C5N_C2R (+40.6 MJ) and
A64_B1R (+48.7 MJ) rise the same and still hold 4 fs for 150 ns. The failure
is a genuine 4 fs detonation: C16 of the C3R linker is flung 42 Å in a single
step as a neighbouring water is pressed to near-zero distance (a hydrophobic
carbon in a tight dye–linker pocket); three runs at 1–2 fs survive, three at
4 fs die within ~100 ps. Ship path: `--dt-fs 2.0` added to `make_library.py`
(same HMR/physics, honest `md_record.timestep_fs`, ~2× wall time), queued via
a cordeshub watcher that waits for the L2R queue to drain and then runs
T51_C3R at 2 fs into the usual `~/dye_out/T51_C3R.log`.

## 2026-08-27 (late) — four headers folded into the families they belong to

`DunbrackLibrary.h` and `Faspr.h` are gone into `RotamerLibrary.h`, and
`RotamerStatistics.h` and `DistanceCalibration.h` into `RotamerEnsemble.h` and
`StatesDistance.h`. Four files that each declared one or two functions and a
page of prose.

**Is the Dunbrack code needed, given `IMP::rotamer::RotamerLibrary`?** Yes,
and the reason is the file format rather than the question. IMP's reads the
Dunbrack *text* library and answers `get_rotamers(residue, phi, psi, thr)`;
this module's reads the **binary** `dun2010bbdep.bin` that FASPR ships, and
the `.drot.pto` container that carries it. Neither can read the other's, and
both are used here: `get_anchor_cb_position` goes through IMP's, and
`faspr_pack` needs ours -- it unpacks the shipped container back to a scratch
`dun2010bbdep.bin` because the vendored engine seeks in that file by byte
offset. `write_dunbrack_bin` reproduces the input byte for byte, which is both
the correctness proof and what keeps the parity pin meaningful. That the two
coexist is now said in the header rather than left to be rediscovered.

`Faspr.h` was 74 lines of documentation for **one function**. What a caller
needs -- what it does, the parameters, what it throws, that it is
deterministic -- is eight lines; the rest was porting narrative and a
redistribution question that PRD-118 owns. The citation and the licence
pointer stay, now aimed at `RotamerLibrary.h`, which is where the vendored
sources' notices point too.

1208 passed, 3 xfailed, 30 subtests; the docs sweep 36. 77 public headers,
down from 81.


## 2026-08-27 (night) — the selection language, and PRD-106 closed

**A `strip_mask` is a selection expression now**, not a four-term dialect.
`SelectionExpression.h`/`.cpp`: a tokenizer, a precedence parser, an evaluator
over a flat atom table, and a compiler to **#IMP::atom::Selection** -- which is
the point. IMP already has the algebra (`set_intersection`, `set_union`,
`set_difference`) and the predicates (`set_chain_ids`, `set_residue_indexes`,
`set_residue_types`, `set_atom_types`, `set_element`), so `chain A and resi
10-20 and not name CA+CB` compiles term by term and composes with everything
that takes a Selection. Four kinds of term have no IMP predicate -- wildcards,
the fields IMP's hierarchy does not carry (`segi`, `alt`, `id`), `index`, and
the whole-structure operators -- and enter the algebra as a Selection over
exactly the particles they matched.

**Two spellings, one parser, detected.** The files in this stack are written
both ways: `resi`/`resid`, `name CA+CB`/`name CA CB`, `50-60`/`50 to 60`,
`segi`/`segname`, `id`/`serial`, `byres`/`same residue as`, `solvent`/`water`.
Both are accepted, including mixed -- every shipped T4L mask *is* mixed. The
spelling is detected from the expression's own markers, because the two
disagree about one thing: whether `index` counts from 0 or from 1. A caller
may declare it instead. `within`/`beyond`/`around` are answered against
coordinates; `pbc` is refused rather than ignored, there being no cell.

**A strip is a radius of zero, not a shorter list** (maintainer, 2026-08-27).
The first cut handed each volume its own filtered particle list. Giving the
stripped atoms **no size** instead is simpler and better: a rasteriser asks
`distance < radius`, so an atom of radius zero blocks nothing, and every
volume indexes the same particles whatever it strips -- which is what the
incremental window machinery wants, and it mutates nothing on a model that
another thread is reading (AVs resample in parallel).

It needed one rule in each rasteriser, because both inflate every radius by
the probe radius before testing it: **a non-positive radius is not inflated
and is not sampled**. Without that, a stripped atom still carves a
probe-sized hole. Three pieces of delta bookkeeping came with it, each a
silent-wrong if missed: `full_raster` records a transparent atom's `last_`
entry anyway, so a later frame has a truthful "before"; `apply_local`
subtracts only what was added; and the moved-set classification does *not*
skip a transparent atom, because `r != last[3]` is exactly what catches one
that has just become transparent.

The rule is one rule everywhere: `strip_obstacles`, the array form, zeroes a
radius rather than dropping a row, so the array keeps row-for-atom
correspondence with the structure. A test proves the equivalence through the
public API -- a volume built from zeroed rows equals one built from the
shorter array those rows were dropped from.

**PRD-106 is done**: six of its seven acceptance criteria as written, and the
seventh's substance by a different rule. Asked whether it really was, the
check found two holes in my own tests -- worth recording, because both would
have passed silently:

* The authored-mask fixture ran **every** shipped mask against T4L's 3GUN,
  including TG2's 26. Those name TG2 residues, so they selected nothing, and
  the assertion was `set() == set()`: 23 of 43 cases proved nothing. Each file
  is checked against its own structure now, and an empty result has to be a
  *checked* claim -- the residue really has no atom outside the keep set.
  Which is how TG2 position A52 turned out to be an **alanine**: its only
  side-chain atom is CB, which the mask keeps, so it legitimately strips
  nothing. That is also why acceptance 7 ("a mask selecting zero atoms
  raises") is not met and should not be.
* The AV parity test compared the number of density points, where the
  criterion says *bit-identical*. Two clouds of equal size can be different
  clouds; it compares the densities now.

Asked a second time, the **Requirements** section (which I had not checked --
only the acceptance criteria) gave two more:

* `strip_obstacles`, requirement 1's second output -- the `(M, 4)` obstacle
  array the flat kernels take -- existed before the C++ port and was lost in
  it. Restored, with the attachment atom kept even when the mask names it.
* Requirement 3 asks for the `("N","CA","C","O","OXT")` keep-set in exactly
  one place, and there were two: the labelling site's, and the selection
  language's answer to `backbone`. They are one list now
  (`protein_backbone_atom_names`), ordered N-CA-C-O-OXT rather than sorted,
  because a mask written from it should read the way a person writes one.

Two findings on the way:

* The decorator path read no `strip_mask` at all. `AV::set_av_parameter` took
  `linker_length`, the radii, `linker_width` and `allowed_sphere_radius`; the
  field was absent from it, so `probe_network_restraint_set`, docking and the
  lattice suite built **unstripped** volumes from files that ask to be
  stripped. It reads it now, and the T4L quadrature score moved from
  22.69935116408601 to 22.593263882270314 -- the number the file had always
  been asking for. `set_av_parameter` takes JSON *text* now, like every other
  reader in this module, rather than an `nlohmann::json` no Python caller can
  build.
* A masked volume cannot join a **shared** occupancy raster: a shared grid is
  one obstacle set for every volume in it, and each position strips its own
  residue. Masked volumes take a private map. `test_av_lattice.py` measures
  the lattice on a mask-free copy of that fixture, so its pre-PRD-105 pins
  still mean what they were recorded to mean, and three new tests pin the
  interaction itself.

**fps.json 1.1** adds a document-level `strip_mask` -- what the *structure*
carries that no volume should be blocked by (waters, ions, an additive) --
beside the per-position one, which is about the labelling site. The two are a
**union**, folded in by the readers so a position taken alone is
self-contained. Additive, so every 1.0 file is a valid 1.1 file.

1208 passed, 3 xfailed, 30 subtests; 37 in the expensive sweeps. PRD-106's
own document, its `Consumers` and `Requirements` sections and its index entry
were rewritten to describe what was built rather than what was wrong: a closed
PRD that still states the defect in the present tense is a trap for the next
reader.


## 2026-08-27 (evening) — TCSPC out of the docs, and every PRD says where it stands

**TCSPC is tttrlib's, and now the documentation agrees.** PRD-113 stage 0
deleted the instrument layer from this package; the manual kept describing it.
`doc/manual/decays/` is gone -- three of its pages called `IMP.bff.DecayCurve`
and friends, two were empty stubs whose toctree entries (`_decay_convolution`)
matched no file, and the sixth was a page about ChiSurf. With them went
`concepts/decay.ipynb` (phasor theory) and `programming/programming_imp_decorator.ipynb`,
whose subject was `DecayLifetimeHandler` -- the decorator pattern is taught by
`structure_accessible_volumes.ipynb` and `getting_started.ipynb` on live API.
`examples/spectroscopy/plot_convolution_routines.py` and `plot_pile_up.py` were
literal duplicates of tttrlib's own
`examples/fluorescence_decay/`, and are deleted rather than moved.

`test/expensive_test_docs_and_examples.py` skipped four pages by name for this
reason; the skip list is gone and **every** notebook under `doc/` now runs. The
manual index carries a note saying where decays live: `tttrlib`'s
`modules/spectroscopy/decay` -- the `fconv` family, `BlindIRF`,
`DecayFit23`--`26`, `MaxEntTcspc`, `DecayStatistics`. `test/standards_exceptions`
lost 45 entries naming deleted `Decay*` classes and `decay_*` functions.

**Every PRD states four things under its title**: Status, Done, Missing, Next.
Sixteen PRDs, one block each, and `okf/prds/index.md` opens with a table of the
status and next move of all of them. Two were verified rather than transcribed:

* **PRD-113 is done.** All eight stages, checked against the tree. Of its three
  recorded todos, one is open (the κ² surrogate, which it never found), one is
  superseded (PRD-110/111 have their own reasons now) and one was already done
  (`PET_QUENCHING_REFERENCE` no longer exists; `QuenchingModel` takes
  `PETParameters`). Its Layout and stage-7 sections describe `pyext/src` domain
  packages that PRD-117 later deleted -- the header says to read them as
  history.
* **PRD-117 is done.** `pyext/src` is empty and no `.i` file holds library
  logic. What remains is `types.i`'s `%attribute_np*` macros -- numpy
  `.reshape` sugar with no C++ spelling -- and one `import numpy as np`.

The index was also missing PRD-109 entirely, and claimed PRD-117 still had
11,607 lines to port.

1099 passed, 3 xfailed, 30 subtests; docs/examples sweep 38 passed with nothing
skipped.


## 2026-08-27 (later still) — a cleanup pass: one of each thing, and the port narrative out of the code

Consolidation, naming and comment scrub across the module. Pre-release, so
nothing was kept for compatibility's sake.

**One of each thing.** The merges had left duplicates that agreed by luck:

| was | now |
|---|---|
| `has`/`txt`/`dbl` in `ForceFieldCIF.cpp`, again in `ProbeLibrary.cpp`, and as `tc_has`/`tc_str`/`tc_int`/`tc_bool` in `ComponentTemplate.cpp` | `internal/CifReader.h`, mirroring the writer header |
| `upper()` in `AVBuilder.cpp`, `Scoring.cpp`, `StripMask.cpp` | `internal/Text.h`, with the other shared string helpers |
| `normalize_weights` (leaves a zero-sum library at zero) and `rs_normalize` (makes it uniform) | `normalize_weights_in_place`, uniform fallback -- a library whose weights were never filled in is one where every conformer counts the same, and leaving them at zero silently drops it out of any average |
| `effective_distance` in `Docking.cpp` shadowing the public one | `modelled_distance`, which says which converter it goes through |

**Three writers with no caller are gone**: `write_component_template_cif`,
`write_probe_template_cif`, `write_rotamer_library`. Nothing in `bin/`,
`examples/`, `doc/` or the siblings called them -- templates are authored by
hand and rotamer libraries are written as `.drot` or RMF -- so their only job
was round-tripping their own readers. Those tests now read the **shipped**
templates and a fixture numpy writes, which is the only arrangement in which a
reader test fails because the reader is wrong about the format rather than
merely self-consistent with our writer.

**Names against flrCIF and mmfdb.** `../mmfdb/src/mmfdb/data/mmfdb_flr_ext.dic`
is the authority. `_mmfdb_probe.probe_type` enumerates `dye`,
`fluorescent_protein`, `spin_label`, `unspecified` -- which `ProbeType` already
matched -- so the work was the surface around it: `DYE_PAIR_*` ->
`PROBE_PAIR_*`, `LlDyeModel`/`LL_DYE_*` -> `ProbeModel`/`PROBE_MODEL_*`,
`DEFAULT_DYE_RADIUS` -> `DEFAULT_PROBE_RADIUS`,
`pair_distribution_from_dyes` -> `pair_distribution_from_probes`,
`data/rotamer_library/R0/dye_library.cif` -> `probe_library.cif` (its constant
already read `PROBE_LIBRARY_CIF`), and the `_mmfdb_operation.algorithm` value
`dye_library` -> `probe_library`. `REFERENCE_DYE` stays: it names an actual
dye, which is a probe *type* in the dictionary. Four include guards still said
`IMPBFF_DYE…`/`IMPBFF_LABELDYNAMICS_H`, and six `\file` tags in `src/` named
the header instead of the source.

**The port narrative is out of the code.** The PRD-113 banner (8 files), ~35
inline "which is what the Python did" justifications, the narrative blocks
atop ten SWIG interface files, and the same in `bin/`, `test/` and one
notebook. Rationale that explains a live decision was kept and moved to the
present tense; what a reader needs is the rule, not its provenance. The
history stays here, which is what this file is for.

**Deduplication against IMP** turned up less than expected: no second PDB or
MOL2 parser, no local RMSD or superposition, no hand-rolled rotation math.
Two overlaps are deliberate and now say so in their headers --
`IMP::saxs::SolventAccessibleSurface` returns a *fraction* of an isolated
surface for a form factor where this returns area in square Angstrom, and
`IMP::algebra::get_uniform_surface_cover` is a different point set from
`sphere_points`, which every sampled quantity in the module is pinned against.

1101 passed, 3 xfailed, 30 subtests; docs/examples sweep 42 passed.


## 2026-08-27 (later) — the same defect a fourth time: three probe fields the container dropped

Yesterday's note said any derived artifact wants a test that reads *it*, not a
copy of it, and that only the probe container had that gap closed. Going
looking closed it properly and turned up a fourth instance of the same defect.

`data/potentials.pto` and the two `template_av_*.fps.json` turned out to be
guarded already -- `test_potential_tables.py` reads the shipped container and
checks its manifest, shapes, NaN handling and proline sentinel;
`test_fps_schema.py` validates the shipped templates against the schema. So the
files were fine. **The fields were not.**

### `dipole_atoms`, `positive_atoms`, `negative_atoms` never survived a round trip

`probe_write_pto` wrote five flrCIF columns and the optical-property triples,
and silently dropped those three. They had no flrCIF item, so they were left
out; the reader had nothing to look for.

`dipole_atoms` is the one that does damage without looking like it. A `Probe`
with no dipole is **legal** -- it means the probe is modelled as an isotropic
emitter -- so the loss produces a valid probe, not an error. An anisotropic
probe comes back isotropic, kappa^2 falls back to 2/3, and the R0 that follows
is wrong by a factor nothing in the output reveals. The charges decide how a
probe behaves beside a charged surface, which is not decoration either.

Fixed: `_mmfdb_probe.{dipole_atoms,positive_atoms,negative_atoms}` in the
dictionary (arrays, following `_mmfdb_spectrum.wavelengths`), written and read.
All fourteen `Probe` fields now round-trip.

### The test enumerates rather than checks three

`test_no_probe_field_is_silently_dropped` walks `dir(probe)` and compares every
field before and after, so a field added later is covered without anyone
remembering. Checking the three that were broken would have left the next one
to be found the same way this one was.

### And the column check was checking a prefix

`test_every_column_name_is_a_dictionary_item` asserted that each column name
*started with* `_flr_probe_list.` / `_flr_probe_descriptor.` / `_mmfdb_probe.`.
That passes for `_mmfdb_probe.dipole_atoms` whether or not the item was ever
declared -- which is precisely the failure the container exists to prevent, a
name that looks like a dictionary item and resolves to nothing. It now looks
each name up in `mmfdb_flr_ext.dic`. Verified the check discriminates:
`_mmfdb_probe.not_a_real_item` and `_flr_probe_list.invented` both fail it. A
guard that cannot fail is the same defect wearing a test's clothes.

1100 passed, 3 xfailed. `data/dyes.mmfdb.pto` regenerated so its column
metadata matches.

## 2026-08-27 — the labelizer port, re-aligned to the restructured package

The package moved a long way under the labelizer port: `cgdye` became
`cgprobe`, `AVNetworkRestraint` became `ProbeNetworkRestraint`, and
`DyeSampling.h`, `DyeDiffusion.h`, `DyeForceField.h`, `CifIO.h` and
`LabelingRestraints.h` were deleted outright along with most of `test/cgdye/`.
This pass checked what that cost the port and fixed what it found.

**The C++ needed nothing.** It builds clean and 1099 pass. The labelizer
includes only `AV*.h`, `StatesDistance.h`, `SolventAccessibleSurface.h`,
`ProbeLibrary.h`, `Pto.h` and its own headers -- no deleted header, no renamed
symbol. Surviving a restructuring of that size untouched is the useful evidence
that "code like AV must not be duplicated" was actually followed, rather than
asserted.

### The shipped container had gone stale, and nothing could see it

`data/dyes.mmfdb.pto` was written 2026-08-25, before `probe_type` and `vendor`
existed. Every one of its 40 probes read back `unspecified` with no vendor
while the CIF beside it had both.

The reason it went unnoticed is worth more than the fix: **every test in
`test_probe_container.py` writes its own container into `tmp_path` and reads
that back.** That proves the round trip perfectly and says nothing whatever
about the file that ships. Regenerated, and
`test_the_shipped_container_still_matches_the_library` now compares the shipped
artifact against `read_probe_library()` field by field, so the next schema
change fails a test instead of shipping.

Same shape as the `probe_type`-holding-a-vendor bug and the reference's
`-1`/`0` sentinel collision: a value that is wrong because nothing reads it.
Three instances in three days is a pattern, not a coincidence -- a derived
artifact needs a test that reads *it*, not a copy of it.

### A rename that reached a section title

`examples/labels/README.txt` had been swept to "ProbePosition simulations" over
a 17-character underline, which is both wrong (the gallery is AV decorators,
screening, path maps and label-site scoring -- not probe positions) and a
docutils error. Now "Labelling simulations", underlined to match. This is the
`dye_name` lesson from the sweep two days ago in a different costume: a renamed
identifier is free, a renamed *literal* is a break, and prose in an `.rst`
title is a literal.

### Docs re-pointed

PRD-120 and `labelizer-correspondence.md` carried `DyeContainer.h`,
`dye_write_pto`, `resolve_dye_name`, `test_dye_container.py`,
`cgdye_fretpredict_pins.json` and "38 dyes". All corrected, and PRD-120 gained
a **What changed after it shipped** section, because a reader who finds a
symbol that no longer exists cannot otherwise tell a rename from a deletion.
The PRD index entry still described the 1.8 dictionary pass as the single term
it started as; it now lists what it actually became.

Verified end to end rather than by inspection: `imp_bff_labelizer` both ways
(195 sites; 20 pairs against the regenerated container, R0 53.8 A), both
examples, and the notebook's whole `bff.*` surface resolved.

## 2026-08-27 (later) — `CifIO.h` is gone: IMP parses, this module only writes

"IMP has CIF, why need another CIF?" -- it doesn't, and it never did. What
looked like a second CIF implementation was one file named after a format.

**Nothing here parses CIF.** Every read in this module -- the force-field
system, the `_cgprobe_*` template, the probe library, trajectories -- drives
`ihm_reader` from `ihm_format.h`, the C parser IMP vendors at
`modules/core/dependency/python-ihm/src/` and uses for its own
`IMP::atom::read_mmcif`. That was already true before this change.

**Writing has nothing to borrow.** The vendored ihm library is read-only (the
string `write` does not occur in its header), `IMP::atom` has `read_mmcif`,
`read_bcif` and no writer, and `IMP.mmcif` is a Python module over python-ihm
that builds deposition systems -- unreachable from C++, and with no spelling
for `_ff_*`, `_cgprobe_*` or `_flr_*`. So ~60 lines of category emission stay,
and the useful question is how many copies of them there are. There were
three: `CifWriter` in `CifIO.cpp`, a weaker `cif_val` in the same file, and a
hand-rolled `_atom_site` loop in `StructureIO.cpp`.

That second one was a live defect. `CifWriter::quote` carries a comment
explaining that quoting must happen in the writer, because the reader splits
rows on whitespace and an unquoted multi-word value makes the loop come up
short. Forty lines below it, `cif_val` quoted on space, `#` and a leading `_`
only -- not commas, not tabs. The force-field writer used `cif_val` for all 58
of its values.

**And consolidating them surfaced a bug in the good one.** `CifWriter::quote`
escaped an embedded `"` by doubling it. That is the CSV convention; CIF has no
escape at all. A quoted value ends at its delimiter *followed by whitespace or
end of line*, so an embedded quote is ordinary text unless whitespace follows
it -- and `""` reads back as two literal quotes. `a "quoted" word` was written
as `"a ""quoted"" word"` and came back `a ""quoted"" word`. The delimiter is
now **chosen** rather than escaped: `"` when it closes the value, else `'`,
else a semicolon text field (whose delimiters each start a line, so nothing
inside can end it early -- `write_loop` puts the next column on its own line
after one). Four values pin the three branches in
`test_an_awkward_value_survives_the_forcefield_round_trip`.

**The fold.** `CifIO.h`/`.cpp` held four unrelated things that shared a file
format. Each went to the header that owns its data:

| was in `CifIO.h` | now in |
|---|---|
| `write_probe_forcefield_cif`, site-id and range utilities | `ForceFieldCIF.h` -- beside `read_forcefield_cif` |
| `ComponentTemplate` + its reader, writer, `region_features` | `ComponentTemplate.h` -- named for the thing |
| `read_rotamer_library` / `write_rotamer_library` | `RotamerLibrary.h` -- with the value they return |
| `CifWriter`, `cif_val` | `internal/CifWriter.h` -- one copy, internal |

`pyext/IMP_bff.cif.i` is `IMP_bff.componenttemplate.i`. `StructureIO.cpp`'s
`convert_pdb_to_cif` writes through the shared class, so a component or atom
name with a space in it is quoted the way every other category's is.
`RotamerLibrary.h` documented itself against `CifIO.h`; it now points at its
own bottom half.


## 2026-08-27 — Dye -> Probe, and cgdye -> cgprobe

The last of the generalisation. A probe is a dye, a spin label or a quencher,
and the structural layer never needed to know which; the naming now says so.
`Probe` rather than `Label` because the tree was already five headers to one
that way -- `ProbeAttachment`, `ProbeContainer`, `ProbeLibrary`,
`ProbeNetworkRestraint`, `ProbeRestraints` -- and a second session was adding
more while this one ran.

| was | is |
|---|---|
| `DyeDiffusion.h/.cpp` | `ProbeDiffusion.h/.cpp` |
| `DyeForceField.h/.cpp` | `ProbeForceField.h/.cpp` |
| `DyeSampling.h/.cpp` | `ProbeSampling.h/.cpp` |
| `LabelDynamics.h/.cpp` | `ProbeDynamics.h/.cpp` -- the one `Label*` holdout |
| `DyeForceFieldSystem` (26 uses) | `ProbeForceFieldSystem` |
| `DyeDiffusionSimulation`, `DyeDiffusionTrajectory`, `DyePairMeasures`, `DyeDistributionNormal` | `Probe...` |
| `build_dye_restraints`, `write_dye_forcefield_cif`, `simulate_dye_diffusion`, `dye_radius`, `find_dye_structure`, ... | `..._probe_...` |
| `bin/imp_bff_dye_pdb2cif`, `--dye-id` | `bin/imp_bff_probe_pdb2cif`, `--probe-id` |
| `data/cgdye/`, `test/cgdye/`, `okf/cgdye.md` | `data/cgprobe/`, `test/cgprobe/`, `okf/cgprobe.md` |
| `get_cgdye_data_dir` | `get_cgprobe_data_dir` |

### What is **not** renamed, and why

**Wire format is not naming.** These are keys in files this package did not
write and cannot rewrite:

- `dye_name` in the fps.json schema -- every fps.json ever written carries it;
- `dye_template`, the default component-template name.

The template CIF's categories **were** renamed, on the maintainer's call:
`_cgdye_template`, `_cgdye_feature`, `_cgdye_feature_atom`, `_cgdye_improper`,
`_cgdye_metadata` and `_cgdye_integration` are `_cgprobe_*`, in the reader, the
writer and the four shipped templates together.

That is a format break, so it fails **loudly**: a template in the old spelling
matches no category, every handler stays silent, and the reader used to return
an empty value -- not even obviously empty, because an unnamed template takes
the file's stem for a name, so `name` was never empty and a guard on it never
fired. The context carries a `saw_category` flag now, and an old file gets the
sed line that migrates it.
`test_a_template_in_the_old_spelling_fails_loudly` pins that.

The sweep renamed all of them. `test_shipped_json_schema_matches_authored_definition`
caught `dye_name` on the first run; the rest came out by reading the diff for
changed *string literals*, which is the check worth keeping for a rename of
this size: **a renamed identifier is free, a renamed literal is a format
break.**

`labelizer` keeps its name (a program), and so do `dyes.mmfdb.pto` and
`dyes.drot.pto` (shipped binaries whose rename buys nothing).

### A self-inflicted one, worth remembering

Restoring the `_cgdye_*` categories with a substring replace of `_cgprobe_`
also hit **`get_cgprobe_data_dir`**, which contains it -- so the header
reverted while the `.cpp` kept the new name. The library exported one symbol
and the wrapper called the other, and every import died with
`symbol not found: get_cgdye_data_dir`. The same over-broad-replace mistake as
the container typemap two rounds ago: a substring is not a token.

### Verification

`ninja` clean; **1092 passed, 3 xfailed, 30 subtests** -- the same count as
before the rename, which is what a rename should do -- plus the 53 of
`expensive_test_docs_and_examples.py`, `test_bin_programs.py` and
`test_dye_commands.py`.

## 2026-08-26 — `DyeLibrary` is `ProbeLibrary`, and a column that lied

`DyeLibrary` modelled a dye, and that was already untrue of the package:
`ProbeAttachment.h` spoke of probes while holding a `Dye`, the rotamer stores
are partitioned `dyes` / `spinlabels` / `sidechains`, and the flrCIF vocabulary
the fields follow calls all of them probes. So `Probe` is the type, with a
`ProbeType` of `dye` / `fluorescent_protein` / `spin_label` / `unspecified`.

**No compatibility shim.** The owner's call, and the right one pre-release:
`DyeLibrary.h` and `DyeContainer.h` are deleted, every call site is renamed,
and there is one name for the thing.

### The column called `probe_type` held the vendor

`_bff_dye.probe_type` was the second column of the bundled library CIF, and its
values are `ATTO`, `Lumiprobe`, `AlexaFluor`. That is a vendor, not a type --
and the reader never read it at all, so nothing noticed. Adding a real
`probe_type` would have collided with a name already used for something else.

The categories are now `_bff_probe` / `_bff_probe_spectrum`, the vendor column
is called `vendor` and is read into `Probe::vendor` instead of dropped, and
`probe_type` means what it says. Every one of the 40 bundled rows is `dye`.

### `forster_radius` refuses a spin label by name

It used to say "both dyes need a spectrum", which reads as a missing file. A
nitroxide has no spectrum *by construction*, so it now says so and names the
probe: a category error and a missing file should not look the same.
`PROBE_UNSPECIFIED` still counts as fluorescent -- a library that does not
state a type is a dye library, and refusing it would break callers over a
column they never had.

### Dictionary 1.8 gains `_mmfdb_probe`

`probe_type` and `vendor` have no flrCIF item, and squatting on
`_flr_probe_list.` for them would collide with a future IHM-FLR one. So
`_mmfdb_probe.{probe_id,probe_type,vendor}` in the mmfdb namespace, bound by
`_mmfdb_schema.table_name` onto the `probes` table the databases already have:
new columns, not a new table. `PtoProfile` enumerates the four type values and
`probe_write_pto` checks against them at the write.

`_mmfdb_probe.probe_type` is a different question from
`_flr_probe_list.probe_origin`, which is *how the probe got onto the molecule*
(intrinsic / extrinsic). A fluorescent protein is intrinsic and a spin label
extrinsic; neither fact tells you which one fluoresces.

### `clustering` added to `_mmfdb_operation.operation_type`

Asked for by the rotamer session, and it was a real gap rather than a
convenience: `mfdb_operation_tags` throws on a non-term, so without it there
was no honest way to record which protocol built a `.drot` library. `analysis`
would have fitted the signature and said nothing -- the same failure mode as
the `probe_type` column above, a wrong value that survives because nothing
looks at it.

Not `ndxplorer_clustering`, which is one tool's burst-selection step. This one
is grouping conformations: a rotamer library built from a trajectory, MSM
microstate construction, `cluster_frames_leader`. More than one caller.

It matters right now because the regenerated rotamer libraries are
**deliberately not comparable** to the shipped ones -- dihedral k-means with a
periodic (cos, sin) metric rather than leader RMSD, which stops the trans well
at 180 deg from splitting in two (0.948 self-fidelity against 0.864; the
FRETpredict libraries they descend from sit at 0.815). Two libraries of the
same probe with different conformer counts have to be able to say why.

### Not renamed, deliberately

`imp_bff_dye_pdb2cif`, `--dye`, `dye_name` in the FPS schema, `DyeDiffusion`,
`cgdye`. The FPS field is literally `dye_name` in an external format, and the
cgdye subsystem simulates dyes. Renaming the description of a field called
`dye_name` to "Probe name" makes it worse, not more general.

1083 passed, 3 xfailed on the full suite.

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
besides PMI's bookkeeping is `probe_network_restraint_set` in
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
`pair_distribution_from_probes` answering with the typed pair values instead of
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

## 2026-09-06 — PRD-135, av.pinn

- [PRD-135](prds/prd-135.md): a voxel-encoder PINN replacing `qpinn/forward.py::simulate`,
  built and accepted on the forward criteria at the frozen 300k checkpoint (τ_x 0.102 ns,
  r_inf 0.0199, calibrated) and on reverse recovery at 88 % coverage. Findings that changed
  the design on the way: a leaking validation split, a silent OOM death of a trainer that
  loaded every chunk, a precision head pinned at its bound (the emitted covariance was a
  constant), a non-finite Hessian behind every "identified, not covered" verdict, and a kink
  in the likelihood behind every "not converged". The identifiability-vs-sites table
  (`prototypes/quench_pinn/figures/s90_identifiability_set2.png`) is the number the
  literature search needed.
- `qpinn/measurements.py`: reads the table's `__error` columns (none populated for the four
  fitted proteins); pools every Alexa Fluor 488 spelling; registers p27 (1JSU) — which loads
  nothing usable, and why is in the PRD.
- The evidence bundle was synced to cordeshub and `pinn_table.csv` recompiled from the
  current `measurements.csv` (14,889 rows, 370 proteins; 36 conflicts and 697 disagreeing
  overlaps flagged for review, previous table kept as `pinn_table_2026-08-27.csv.bak`).

## 2026-09-07 — PRD-136, dye-timewarp

- [PRD-136](prds/prd-136.md): a site-conditioned, time-coarsened generative model of dye+linker
  dynamics learned from MD (Timewarp line; TITO design). Prior-art survey: no published learned
  dye-dynamics surrogate; the machinery (Timewarp, ITO/TITO, UniSim, FBM) exists and TITO's code is
  MIT (`junk/tito`, with pretrained checkpoints). Decisions by interview: Cartesian equivariant state,
  hybrid MD campaign on heinzehub, gates from PRD-119/135 + PRD-99 + PRD-116, Alexa488-C1R + Alexa647-C2R.
- `prototypes/dye_timewarp/s00_inventory.py`: the 19 free-dye runs on SD1TB are byte-identical to
  `MANIFEST.tsv` (2660 ns, 532,000 frames at 5 ps); the lag ladder 5 ps … 5 ns is supported; the
  AF4_C2R and C3B_C2R replica sets are the noise floor for every later gate.
- Stage 0/1 the same day: TITO (`junk/tito`, MIT) is reused as the model core; the free-dye corpus grew from 19 to
  **65 runs / 12.9 µs** by pulling the 53 finished heinzehub runs (SHA-256 verified, 7 duplicates dropped). `s02` on the
  full corpus: dipole rotation converged and replica-tight (C(t) spread ≤ 0.026), but the slowest torsion processes run
  5–190 ns, so a single 140 ns run is not converged in its marginals -- the FRETpredict-protocol libraries inherit that,
  and every marginal gate here is scored against pooled replicas. Stage-1 gate PASS. Full connectivity OOMs a 20 GB GPU
  at ~100 atoms; a 0.6 nm radius graph trains at batch 32. `free_v1` launched on heinzehub GPU 1.
- Arm B started the same day: `s06_build_site.py` builds and runs a labelled site (graft + fragment + tleap + OpenMM);
  `s07_campaign.sh gold` queued on heinzehub GPU 0 behind the open-linker MD. Gold set = the 13 cysteine-maleimide
  Alexa488 sites (5 T4L, 8 hGBP1); the ten T4L pAcF sites are a different linker chemistry and are parked.
- Start conformers for arm B now come from the shipped rotamer-ensemble placement (`s06a`, user's suggestion), which
  also yields each site's effective sample size — the gold sites span ESS 5.7 (T4L 127) to 115 (T4L 86).
- Owner decision: parameterise the Alexa488–pAcF ketoxime label (`P1R`) so the ten T4L anisotropy sites become
  simulable; filed as `T-20260907-05` for the rotamer-library builder (PRD-116), consumer PRD-136.
- `T-20260907-05` resolved without work: AMBER-DYES `B1R` is the Alexa488–pAcF ketoxime label (found by
  rotamer-simulation, confirmed by substructure match); the ten T4L anisotropy sites enter the campaign as `B1R`.

## 2026-09-07 — an independent core, steps 1–3 (PRD-137, cross-stack; owner decisions in the chisurf bundle)

- **Why**: IMP is hard to install and slow to release; chisurf needs bff to fit fluorescence data and needs no IMP.
  Measured: chisurf uses 38 `IMP.bff` symbols and 36 need no IMP. Architecture chosen: an IMP-free core plus a
  *connection layer* (`src/imp/`, later) so bff keeps working inside IMP. Not a replacement for IMP.
- **Step 1, `include/Base.h`**: one door for `IMP_THROW` / `IMP_SHOWABLE_INLINE` / `IMP_VALUES` / `IMP_OBJECTS`.
  95 files rewired, 0 of ~600 call sites touched. `libimp_bff` exports 2182 symbols before and after. The
  `IMPBFF_STANDALONE` branch compiles and runs with no IMP on the include path (`test/test_base_header.py`).
- **Step 3, `src/RmfIO.cpp`**: written against RMF's own API, no `IMP.rmf`. 33/34 shipped `.rmf3` bit-identical; the
  T4L docking trajectory matches IMP to 3.5e-6 Å (float32 storage). Two bugs the oracle caught: reference frames are
  *frame* data (compose per frame, not once — 43 Å otherwise), and `ReferenceFrameFactory::get_is` answers from the
  current frame (position on frame 0 before the structural walk). Dropping `rmf` from `dependencies.py` is deferred:
  it works out-of-tree but silently disables IMP.bff in-tree (configure ordering writes `build_info/RMF` after bff).
- **Step 2, `include/DensityGrid.h`**: `PathMap` no longer derives from `IMP::em::SampledDensityMap`; `IMP::em` is
  confined to `write_map_feature` and the new explicit `PathMap::create_density_map()`. AV output bit-identical
  (`test/test_density_grid.py`, 9 golden records). Five regressions found and fixed, none visible to the golden
  records: `set_origin` must rebuild the coordinate caches; SWIG must *see* the new base or `PathMap` silently loses
  every inherited method; the header is `float`-backed like IMP's (double moved mean positions in the 8th figure);
  index arithmetic promotes to double after float inputs, as IMP's does; and obstacles are re-read from the *live*
  particles at every sample — IMP held decorators, a value copy sampled frame 0 forever (15 Å on trajectory tests).
  `PathMapHeader`'s serialized form changed (11 fields, was ~40). Suite: 1780 passed, 0 failed.
- **Benchmarks**: `benchmark/benchmark_av.py` (the AV solver timing gate, with `--noise`) and `benchmark/paired_ab.py`
  (A/B two self-contained builds alternately, median of per-round ratios; judge on median *and* minimum; 11 rounds —
  seven flipped verdict at load 13–18). Provisional baseline in `okf/validation/av_solver_baseline.md`.
- **Fork**: `tpeulen/IMP.bff` created (public, parent `Fluorescence-Tools/IMP.bff`), remote `fork` added. Owner
  decision: commit locally on `independent-core`, **do not push yet**. Layout stays flat (IMP-compatible); merges only.
- Owner narrowed PRD-136 to Alexa488 maleimide (`A48_C1R`) for a conference: `s07 conf` tranche (13 C1R sites,
  fragments first), probe loop on A48_C1R, pAcF/B1R and A64_C2R runs deferred; GPU 0 requested from rotamer-simulation.
- **PRD-138 merges done** (imp-bff-ce): 14 commits on `independent-core`, 97 → 75 headers, 106 → 76 sources,
  40 → 35 topic `.i`; suite 1780/0 after each. Lessons: SWIG `%include` is first-wins, so a merge that folds a
  late-wrapped header into an early-wrapped one must move the survivor *late* and carry the value declarations
  with it (Rotamer, FRET); typemaps are global, so headers can move freely past them (Quenching); a header-only
  merge does not reconfigure — run `cmake .` (Cif). Narrowed rows: QuenchingModel (cycle via InteractionTerms),
  FPSExport/Project, ProbeDynamics, ProbeNetworkRestraint (IMP-side, PRD-137), and `AVModel` → `States` +
  `AVModel` after the owner objected that an AV is one representation ("there are also rotamers"). Pre-existing,
  uncollected failure noted on the board: `medium_test_av.py` / `States.pRDA()`. Build lock released.
- **PRD-138 interface pass done** (imp-bff-ce): 19 `compute_/build_/find_` free functions → `get_`/`create_`
  (hard rename, no aliases; 79 files; three genuine IMP.bff references in `../chisurf` patched, uncommitted);
  raw-pointer audit: no change (IMP-object returns in layer code, one hidden buffer, one non-owning association);
  one `%typemap(out) std::vector<double>` → ndarray owning the moved vector, 91 functions, zero copy. Suite 1780/0.
  Downstream to expect: bff vector returns are ndarrays, not tuples (`assert not x` → `len(x) == 0`).
- **PTO container now ptolib** (fable-5.1/b19b35ca, T-20260907-07 on the shared board): the
  hand-written EBML walker in `Pto.cpp` (PtoWriter/PtoReader, ~570 lines) is replaced by a verbatim
  copy of ptolib (https://github.com/tpeulen/ptolib, private; `include/internal/ptolib.h`,
  `utility/sync_ptolib.sh`, byte-compared by `test/test_vendored_headers.py`), the same header tttrlib
  now builds on. `Pto.h` keeps `PtoWriter`/`PtoReader`/`PtoObject` as a thin face so the 25 call sites
  compile unchanged, and exposes `.file()` for tags, annotations and `pto::DataStore` tables. Files written
  from now on carry two indexes, an Info element and random uids; every shipped `.drot.pto`,
  `potentials.pto` and `.mmfdb.pto` still opens (the reader accepts the old single-index layout, read-only).
  `Pto.cpp` is the module's `PTOLIB_IMPLEMENTATION` TU with `PTOLIB_JSON_INCLUDE` pointed at the vendored
  nlohmann copy. The header compiles as C++14 (ptolib shims std::filesystem) so the module's standard and
  its macOS 10.13 floor are untouched.
- **PRD-137 steps 4+5 (main cut)** (imp-bff-ce): connection layer = `src/imp/` behind the flat bridge
  `src/ImpLayer.cpp` (IMP's tooling: only `src/*.cpp`+`src/internal/*.cpp` compiled, only `include/*.h`+
  `include/internal/*.h` linked -- a public `include/imp/` cannot exist); Scoring's five restraint factories →
  Potentials.h; kb_kcal → Scoring.h; `get_av` Model-free via `OccupancyGrid` (core raster; `AVOccupancyMap` is
  the layer's particle view) + `get_av_lattice` (the decorator's lattice path over spheres), 18/18 array-door
  records byte-identical; `AVPairDistanceMeasurement` → FPS.h; the Model doors (`resample_av`,
  `get_av_from_structure`, `get_avs_for_structure`) → AV.h. No core file includes a layer header. Suite 1780/0.
  Shared-tree note: another agent's uncommitted ptolib vendoring does not link under the unity build
  (PTOLIB_IMPLEMENTATION skipped by the header guard) -- built against HEAD's Pto files, restored after.
- **PRD-137 step 5, interface complete** (imp-bff-ce, 813a418): 21 Hierarchy/Particle overloads → `HierarchyBridge.h`
  (layer), ProbeAttachment → layer, PathMap takes spheres from a pluggable source (`set_path_map_particles`; Python
  `PathMap.set_particles` kept via %extend). No core header names a particle/hierarchy/model. Suite 1780/0, oracle
  18/18. Residue 5c (implementation): `load_protein_frames`, `load_structure`, the trajectory loader's PDB branch,
  `selection_from_expression` still read PDBs / compile the AST through IMP::atom -- marked at their definitions.
- **PRD-137 step 5c** (imp-bff-ce): the core reads PDB (`internal/PdbFrames.h`, IMP's rules incl. float-rounded
  coordinates) and MOL2 itself; `selection_from_expression` built from the core's evaluation. Gated by equality
  tests against the IMP road on every shipped .pdb/.mol2, 240 selection records, AV oracle 18/18. Suite 1783/0.
  No core source names a particle/hierarchy/model. Left: PathMap's IMP::em interop (core MRC writer), and the
  IMP::Object/algebra/Pointer vocabulary in nine headers (Base.h shims, step 6).
- **PRD-137 em cut** (imp-bff-ce): `DensityGrid::write_mrc` writes a correct MRC2014 map (real statistics,
  n*spacing cells, NVERSION, machine stamp; validated by `mrcfile` and read back by IMP.em); `write_map_feature`
  uses it; `create_density_map` is the layer's `EmBridge.h` (Python spelling kept). Owner ruling: correct MRC,
  no byte-copy of IMP's (which wrote NaN statistics and, for features, cell lengths ignoring the spacing).
  The core names IMP::em nowhere. Suite 1786/0 (one failure = another agent's in-progress ptolib test).
- **PRD-137 step 6a** (imp-bff-ce): `standalone/include/` -- the IMP/ shim tree (config, Object, Pointer,
  constants, algebra::VectorD) -- and the core compiles with no IMP on the path (expensive test). RRT's
  Transformation3D conversions → `AlgebraBridge.h` (layer). The Labelizer has a replaceable AV door: core road
  (own PDB reader + vdW table) vs IMP's road installed by the layer at load; radii table choice deferred to the
  owner. Suite 1784/0. Next: 6b, the standalone CMake build of libimp_bff (Eigen, cereal, Boost headers, RMF,
  vendored python-ihm C parser), then the standalone SWIG entry, then packaging.
- **PRD-137 step 6b+6c** (imp-bff-ce): `standalone/CMakeLists.txt` builds `libimp_bff` with no IMP (Eigen,
  cereal, Boost headers, RMF, vendored python-ihm C parser; links RMF/libomp/libc++ only), and now the Python
  module: `standalone/pyext/IMP_bff_standalone.i` wraps the same `IMP_bff.core.i` the IMP build wraps, under
  `IMP_bff_standalone.macros.i` (IMP_SWIG_* equivalents, VectorD tuples, IMP's exception family as Python classes
  with IMP's std:: mapping, IMP's director registry, %implicitconv). `import IMP.bff` with no IMP loaded; the
  surface is the IMP build's minus IMP's config constants and the layer, checked name by name. One lane for both
  builds: `test/conftest.py` deselects files naming another IMP module or a layer name when `IMP.atom` is absent
  -- standalone 718/0, IMP build 1790/0. Found on the way: a value-vector `in` typemap must copy the item before
  dropping its reference (a wrapped std::vector's `__getitem__` hands out an owning temporary); the two
  `AVPairDistanceMeasurement` definitions had stayed in the layer. Shared-tree note: another agent's ptolib /
  ExpressionEngine move (uncommitted) does not link under the unity build in either state, so the IMP gate was
  built with their six files at HEAD (plus `cmake .`, since IMP links headers by symlink at configure) and they
  were put back byte-for-byte. Next: 6d packaging (pyproject on scikit-build-core, wheels, the second recipe).
- **PRD-137 step 6d** (imp-bff-ce): two packages, one import. Package `bff` = the core (root `pyproject.toml`,
  scikit-build-core over `standalone/`; `conda-recipe/recipe-core.yaml`), package `imp.bff` = the IMP module
  build as before; both import as `IMP.bff`, both own `site-packages/IMP/bff/`, so the recipes constrain each
  other and `IMP.bff.get_build()` tells them apart (owner's naming, 2026-09-08). The wheel carries ~2 MB of
  `data/` and fetches `rotamer_library/` + `cgprobe/` (62 MB, 284 files) through `data/registry.json` (sha256,
  `utility/data_registry.py`) with pooch from `https://www.peulen.xyz/downloads/imp-bff-data/`
  (`get_data_path()`, `fetch_data()`, `imp_bff_fetch_data`; `IMP_BFF_DATA` is now a PATH-like list). No RMF on
  PyPI: `utility/wheel_deps.sh` builds RMF 1.7.1 (avro only) + cereal on a wheel image. CI gains `build_core`,
  `wheels` (Linux x86_64, macOS arm64/x86_64; Windows deferred) and `publish_wheels` (trusted publishing on
  release). Gates: wheel built and installed into a venv with no IMP -- import, shipped data, registry files
  fetched from a local server and verified; the same after `delocate-wheel` against an RMF built by
  `utility/wheel_deps.sh` (the script patches RMF's CMakeLists, which runs FindHDF5 directly, so the library
  links no HDF5); standalone lane 724/0; IMP build 1794/0; the core conda recipe built with rattler-build,
  passed its recipe tests, installed from the local channel and ran a slice of the suite (3 failures = tests
  that assume tttrlib is installed, a pre-existing suite property). Found on the way: `%exception` must precede the `bff_config.h` include
  (get_data_path's wrapper had no handler and an IOException terminated the process); `install(TARGETS)`
  COMPONENT binds per artifact group; install_name_tool + strip leave a stale signature on arm64 macOS
  (SIGKILL) -- no strip, codesign after install; conda-forge's clang refused `static const int kAxisSlot`
  ODR-used in Expression.cpp (now constexpr). Owner's steps: upload the two data directories to the host,
  register `bff` on pypi.org with the workflow as trusted publisher. PRD-137 is implemented end to end.
- **PRD-137 6d, the data leaves git** (imp-bff-ce): `data/rotamer_library/` + `data/cgprobe/` (62 MB) are no
  longer tracked -- `data/registry.json` is, the files are on the download host (`/Volumes/Downloads/imp.bff`,
  `https://www.peulen.xyz/downloads/imp.bff/`), `utility/data_registry.py --fetch` restores a checkout, CI caches
  the download keyed by the registry and downloads only on a miss (tttrlib's pattern), conda builds fetch before
  packaging (`pooch` in build requirements), and one shared `pyext/IMP_bff.data.i` gives both builds
  `get_data_path()` that fetches a registry file on first use (it wraps the extension's function, which both
  shadow modules call -- IMP's tooling defines the shadow's `get_data_path` after any module code). Root
  scratch is gitignored and the recipes honour it (`use_gitignore`).
- **PRD-137 6d, LabelLib's interface** (imp-bff-ce): `IMP.bff.labellib` is LabelLib's API name for name --
  `dyeDensityAV1/AV3`, `minLinkerLength`, `addWeights`, `meanDistance`, `meanEfficiency`,
  `sampleDistanceDistInv`, `Grid3D` and the `_arr` twins, same argument order, (4, N) or (N, 4) atoms, grids
  flat with x fastest -- so a LabelLib script moves over by changing its import. New C++ behind it:
  `get_av_from_pdb` (the core's file door, the road the Labelizer's AV mode already took) and
  `get_linker_path_lengths` (the lattice search read out as path lengths, a `DensityGrid`; unreachable voxels
  negative, as LabelLib marks them). The lattice search and the obstacle preparation are now shared helpers,
  so a volume and its path lengths are two read-outs of one search. Gated by `test/test_labellib_interface.py`
  (17 cases; the parity ones run against the real LabelLib where it is installed -- mean distance within 0.5 A,
  mean efficiency within 0.02). Owner's ruling: the Langevin/attachment dye road stays IMP-only.
  README gains install instructions for both packages, the LabelLib migration and a first-volume example.
- **PRD-137, where dye simulation stops without IMP** (imp-bff-ce, owner ruling 2026-09-08): the core covers
  rotamer-library dyes, grid diffusion, RRT linker sampling, probe force fields and FASPR packing; attaching a
  dye onto a hierarchy and Langevin dynamics on it stay IMP-only. Linking IMP into the core was considered and
  is possible but does not remove the shims: no translation unit can hold both the shim tree and IMP's headers
  (they redefine `IMP::Exception`, `IMP::Object`, `IMP::Pointer`, `IMP::algebra::VectorD` -- measured), while
  the same sources built against real IMP with no shims do compile (83 of 84 in a hand-rolled IMP-mode build;
  the one failure is a logging-macro flag). So IMP-as-a-library is a second configuration of these sources, and
  the shims stay because they are what lets the core build with no IMP present at all (the pip package). PRD-139
  carries the follow-up question: whether the IMP-mode build should leave IMP's in-tree tooling for this
  project's own CMake (`find_package(IMP)`), which would retire the flat-header rule, `src/imp/Headers.cmake`,
  the `src/ImpLayer.cpp` bridge and `dependencies.py`. The core conda package (`bff`) builds and passes its recipe tests on this machine
  with the data fetched at build time; `build_core.bat` gives it Windows, so CI runs it on all three platforms.
- **What a core Langevin would need** (imp-bff-ce, measured 2026-09-08, recorded in PRD-137): the topology is
  already core (`ProbeForceFieldSystem` and its builders); the force-field energies are array-based and
  IMP-free in their signatures even though they live in the layer's `Potentials.h` (`clash_energy`, `go_energy`,
  `lennard_jones_bead_energy`, `generalized_born_energy`, `ramachandran_energy`, `residue_asa` -- all over
  `const std::vector<double>& xyz`), so they would move like `States.h` did. Missing everywhere in bff:
  **gradients**. IMP's restraints supply the derivatives its integrator steps on. So a bff-only Langevin is
  analytic gradients plus an integrator, not attachment plumbing. Stays IMP-only by the owner's ruling.
- **A wheel can carry IMP** (imp-bff-ce, measured 2026-09-08, PRD-139): "IMP has no wheel" was the wrong
  objection -- a wheel needs IMP's C++ libraries, not its Python package, and can vendor them exactly as this
  wheel already vendors librmf, Boost and libicu. The closure bff needs is 13.9 MB of shared libs,
  `libimp_kernel` links no Python, and IMP's CMake guards every pyext directory with `if(NOT IMP_STATIC)`, so
  IMP's wrappers can be left out. One hard condition: the Python surface must stay IMP-free, or SWIG needs
  IMP's .i files and the extension imports `_IMP_kernel` at init. If every shipped configuration links IMP,
  the shim tree, Base.h's second branch and the two-lane test hook can all go -- one build, and dye Langevin
  works from `pip install bff`. The price: every user carries IMP's code, CI builds IMP per platform, and
  PRD-137's "links no IMP" premise is given up. Owner's decision, recorded as three options in PRD-139.
- **Minimales, eingebettetes IMP: gemessen und gebaut** (imp-bff-ce, owner's choice 2026-09-08, PRD-139): the
  layer names only atom, core, algebra, container plus one `IMP::em` and one `IMP::rotamer` function; the core
  mentions IMP in comments only. Required libraries kernel/algebra/display/score_functor/core/container/atom =
  10.1 MB, the two narrow spots (em, statistics, rotamer) 1.4 MB more. A minimal IMP -- ten modules, library
  targets only, no SWIG, nothing linked against Python -- builds in **57 s**; `libimp_bff` with the connection
  layer builds against it in **71 s** and is 6.4 MB. `standalone/CMakeLists.txt` gained `IMPBFF_WITH_IMP`
  (exclusive with the shims by nature), `utility/wheel_deps.sh` builds that IMP behind `WITH_IMP=1`, `Base.h`'s
  IMP branch now includes `<IMP/log_macros.h>` (IMP_WARN was undefined and its stream argument got compiled),
  and `bff_config.h` self-defines `IMPBFF_STANDALONE` only when `IMPBFF_WITH_IMP` is absent. Layer symbols
  verified present; the IMP-free lane still 724/0. Open: the Python surface -- the layer's IMP-typed signatures
  would drag `_IMP_kernel` into the extension, so array and file doors are what keep a wheel IMP-free.
- **Das Wheel von 20,7 auf 5,4 MB** (imp-bff-ce, 2026-09-08, PRD-139): the fat was not IMP but RMF's transitive
  chain -- librmf links Boost.Iostreams, which as packaged (conda-forge and Homebrew alike) links Boost.Regex
  and all of ICU: 34.8 MB of the 45 MB the repaired wheel bundled. RMF itself is one source file
  (`src/RmfIO.cpp`, 93 uses) behind four functions whose only callers are two `bin/imp_bff` commands; the whole
  lane passes without them (724/0). So `IMPBFF_WITH_RMF`, opt-*out* (`IMPBFF_NO_RMF`) so IMP's tooling keeps
  building what it always built; the wheel turns it off, both conda packages keep it. Also: `install.strip` back
  on (the install re-signs afterwards), and `wheel_deps.sh` now fixes RMF's install name, which pointed into the
  build directory -- anything linked against a hand-built RMF failed at import, and only delocate hid it.
  Hidden visibility was tried and reverted: 0.23 MB for four broken tests, because the exception classes' typeinfo
  stops being exported and the catch that translates them sits in the extension while the throw is in the library.
- **Dye-Langevin ohne IMP-Modul** (imp-bff-ce, 2026-09-08, PRD-139): `attach_dye_to_pdb` and
  `run_dye_langevin` (`include/DyeDynamics.h`, `src/imp/DyeDynamics.cpp`) are the connection layer's dye roads
  behind file paths and arrays -- model, hierarchies and particle indexes are made and destroyed inside. On a
  build with the embedded minimal IMP a real trajectory came back (5 frames of 83 atoms with its energies) from
  a Python where `IMP.atom` is not importable. The rule that makes it work: nothing wrapped may name an IMP
  type, or SWIG needs IMP's .i files and the extension imports `_IMP_kernel` at load. Cost of the plumbing:
  SWIG gets IMP's macro spellings (`IMP_VALUES` -> `IMP::Vector`, both `VectorD` spellings), `IMP::Vector` gets
  the value typemaps, the vector `out` typemap binds a reference (SwigValueWrapper), the wrapper includes
  `<boost/range/distance.hpp>` (IMP's VectorD calls boost::distance), the director registry no longer asks for
  `get_ref_count` before holding a subclass (IMP::Object is unwrapped here -- it held nothing and a cancelled
  minimizer run segfaulted), and the wheel carries IMP's atom data (3.1 MB; score_functor's 18 MB stays out).
  `get_build()` now answers core / core+imp / imp. `test/test_dye_dynamics_doors.py` asserts in a subprocess
  that running a dye imports no IMP module. Lanes: 724/0 IMP-free; IMP-linked differs only by the new tests.
- **Eine Form für jede Simulation** (imp-bff-ce, owner's review 2026-09-08, PRD-139): the first cut of the dye
  doors was two free functions taking PDB paths; the owner refused both halves -- a path where an array belongs,
  and a one-shot call where a simulation object belongs ("like in OpenMM"), general enough for Langevin, MD,
  diffusion, lattice diffusion and a learned model, "maybe json style". `include/Simulation.h` is that shape:
  `SimulationTrajectory` (the module's one trajectory record; `LangevinTrajectory` is a typedef) and
  `ProbeSimulation` (`get_type`, `get_parameters`/`set_parameters`, `get_positions`/`set_positions`, `minimize`,
  `step`, `run`, `has_energy`). Parameters are JSON because that is what lets one interface carry kinds with
  disjoint knobs; state stays arrays. `DyeSimulation` takes two `ProteinFrame`s and the dye's MOL2, keeps IMP
  behind a pimpl, and refuses parameter changes after the system is built. `ProbeDiffusionSimulation` joins the
  interface: its old `run(D, ...)` is `simulate(...)` and `run(n_steps, write_every)` now means what it means
  everywhere else. New bridge `hierarchy_from_protein_frame` (the inverse of `protein_frame_from_hierarchy`)
  is what lets arrays reach IMP's machinery with no file in between. The same six lines drive both simulations;
  lanes 724/0 IMP-free, IMP-linked equal in the same environment.
- **OpenMP läuft zum ersten Mal** (imp-bff-ce, 2026-09-08, PRD-140 Schritt 0): the standalone build now gets
  `-Xclang -fopenmp` -- CMake finds it on a fresh configure, and where it does not, the CMake looks for `libomp`
  itself. `built_with_openmp()` went false -> true and `parallel_threads()` 1 -> 8; the suite passes with the 19
  pragmas executing for the first time (724/0). Two traps recorded in
  `okf/validation/openmp_is_not_enabled.md`: the empty flags live in the *cache* of ~/dev/imp's build, not in
  this compiler, and the flags string has to be `separate_arguments`-split or `-Xclang -fopenmp` reaches the
  compiler as one unknown option. Also measured, and against expectation: `diffusion_propagate` at ng=41 is
  *slower* with eight threads (0.138 vs 0.115 ms/step, machine at load 18.9 so not a verdict) -- the pragma sits
  inside the per-step loop over x-slabs, and a 41-voxel slab does not pay for four thousand forks and joins.
- **Die Rechen-Tür, und das f32-Tor** (imp-bff-ce, 2026-09-08, PRD-140 Schritt 1): `include/Compute.h` is where a
  kernel can be sent somewhere else -- a C-ABI plugin struct (`ImpBffComputeBackend`, versioned), a `dlopen`
  loader that declines rather than crashes, and `get_compute_backend_name()` so nobody measures a CPU and
  reports a GPU. `diffusion_propagate` consults it and hands over the *whole* propagation or none of it; a
  backend may decline per call. Python does the finding, because that is where `import wgpu` is natural:
  `enable_gpu()` honours `IMP_BFF_GPU=off` and `IMP_BFF_GPU_LIBRARY`, looks for the plugin beside the module and
  for wgpu-native inside a wgpu-py installation, and warns once if it found a plugin but no library. No OS ships
  WebGPU; wgpu-py does, for every platform this library targets, and it is already installed here.
  **The gate that had to come first**: WGSL has no f64. Measured (benchmark/f32_diffusion_gate.py,
  okf/validation/f32_holds_the_diffusion_stencil.md): the same sweep in f32 against f64 deviates 1.2e-6 relative
  in the fluorescence trace and 1.5e-6 in the final density at ng=41 over 4000 steps -- the scheme is
  contractive, so error is damped rather than accumulated, and longer runs come out *closer*. The adjoint is not
  covered and needs its own measurement. Lane 732/0.
- **Der Stencil auf der GPU: 5-12x** (imp-bff-ce, 2026-09-08, PRD-140): `benchmark/gpu_diffusion_wgsl.py` is the
  WGSL the backend will run, held against the CPU kernel. Against the eight-threaded CPU it is 5.3-5.6x at every
  grid size; against the single thread that shipped until today, 5.9x at ng=41 growing to 12.3x at ng=81 -- and
  that growth is the point, because `dg^-5` puts the science exactly where the CPU is worst. Deviation 1.2e-6 to
  1.1e-5, matching what the numpy f32 gate predicted from the other direction. Two traps recorded in the file
  and worth knowing before writing the C plugin: `queue.write_buffer` runs immediately while dispatches run at
  submit, so a uniform cannot carry a per-dispatch value; and a dynamic storage-buffer offset did not take
  effect in this wgpu at all (every reduction wrote slot 0, with a constant as well as with the sum) -- a
  `copy_buffer_to_buffer` between passes is the way around both, and costs nothing. Also measured on the way:
  OpenMP is worth 1.0x at ng=41, 1.8x at 61 and 2.3x at 81, which is the same small-slab story the OpenMP note
  tells. Prototyped through wgpu-py because it already carries wgpu-native; the C plugin is next.
- **Drei Optimierungen am Stencil** (imp-bff-ce, 2026-09-08, PRD-140): dispatching only the active voxels
  (a third to two fifths of the cube) is worth 1.3-1.5x; a workgroup of 256 another 4%; and precomputing the
  stencil weights 1.6x. The last is the one worth remembering: `d`, `decay` and `bounds` are constant over a
  propagation, so writing the step as `[decay*(1-sum a)]*p0 + sum [decay*a]*cur[m]` and computing the six `a`
  once turns eighteen scattered reads per voxel into six, with the weights in a coalesced layout -- and it is
  exact, deviation against the f64 CPU unchanged. Together 7.7-9.6x against the eight-threaded CPU where the
  naive kernel managed 4.8-5.9x. Untried and listed in the note: f16 weights, tiling with workgroup memory,
  renumbering the density into the compacted layout, hoisting the CPU parallel region out of the step loop, and
  -- the real prize -- an integrator that does not need four thousand steps.
