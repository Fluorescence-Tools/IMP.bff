# Handover — PRD-117 `%pythoncode` → C++ port

**Date:** 2026-08-21 (evening), updated 2026-08-22
**Repo:** `/Users/tpeulen/dev/imp.bff`
**Build:** `cd /Users/tpeulen/dev/imp/cmake-build-arm64 && ninja IMP.bff`
**Tests:** `cd /Users/tpeulen/dev/imp.bff && /Users/tpeulen/mambaforge/envs/arm64/bin/python -m pytest test/ -q -p no:cacheprovider --ignore=test/medium_test_av.py --ignore=test/cgprobe -k "not medium_test"`
**Suite state after last clean build:** 656 passed, 3 xfailed, ONLY pre-existing failure is `test_access_av_feature` (missing `references/av_reference_0.mrc` data file — unrelated to the port). cgprobe 149 passed, expensive all-dyes 36 passed.

## STATUS UPDATE (2026-08-22) — Phase 2 four files in

**Phase 1 complete; Phase 2 four of five files done.** Batches 6-14 all
committed green on `dev` (last `4e25e38`): avbuilder (nlohmann JSON doors),
petquenching, dyesampling, photophysics, quenching (Phase 1), then
avmeandistance (`2b0f558`), fps (`c1e2d5f`), structureio (`64a9261`),
cif (`4e25e38`). cif.i's port included writing the real ihm C-reader template
CIF integration (the C++ stub returned empty) and the
`forcefield_system_from_json` nlohmann bridge in DyeForceField.cpp.

**NEXT = `scoring.i` (671 lines)** — the first file that is a C++ *authoring*
effort, not shim removal:

1. `charmm36_lj`/`lj_cross`/`lj_parameter_arrays` exist in Scoring.h but are
   `%ignore`d (they return `std::array<double,2>`); the %pythoncode re-implements
   them plus `CHARMM36_LJ` (a dict tests index with `["rmin_half"]`). Plan:
   change the C++ returns to out-views/wrappable types, un-ignore, delete the
   Python duplicates, migrate the dict-indexing test callers.
2. Port `lj_energy`, `boltzmann_weights`, `rotamer_cluster_weights`,
   `BoundingBoxFilter`, `rotamer_mean_field_weights(_multi_dye)`,
   `compute_rotamer_score` orchestration, `compute_lj_pair_sites`/
   `site_element_map`/`DyeInternalEnergyEvaluator` (typed-system walkers) to C++.
3. The FRETpredict **selector string mini-language** (`_selector_matches`,
   `selector_resnames`, `_selector_atom_names`, `_site_mask`, charge masks) --
   `selector_atom_indices` already exists in C++ (Scoring.h); the mask builders
   around it need the same treatment or a documented Python landing.
4. `torsion_cosine`/`build_dye_restraints` are IMP-object glue (Cosine,
   DistanceRestraint) -- sanctioned Python, like the RMF lazy builders.

Then the four 1000+ line files (`label.i` 987, `topology.i` 987, `sampling.i`
1047, `rotamer.i` 1404, `docking.i` 1429, `sim.i` 1488) -- heavy IMP.pmi/IMP
API glue; large genuinely-Python parts there will need the same
sanctioned-exception treatment (lazy builders, object glue) with the kernels
still going to C++.

## The command from the user

> "continue, MUST be done! no python."
> "you can leave the attributes in swig you can map methods to attributes"
> "use nlohman json, header only can add internal"

So: port ALL `%pythoncode` in `pyext/*.i` to C++. **No Python in any `.i` file.** No hand-written `PyObject*` typemaps. Use SWIG stock machinery + numpy.i suites + `nlohmann::json` (vendored at `include/internal/json.h`) for dict/JSON bridges.

## THE ONE RULE that cost the most time

You may change **anything** (interfaces, tests, fixtures) — this is prerelease. I change C++ signatures freely and update the bff tests to the new surface. **The non-medium bff test suite is the source of truth.**

## Patterns that are locked in (do not re-litigate)

1. **Array fields = `get_*()`/`set_*()` methods returning flat numpy views**
   (`ARGOUTVIEWM_ARRAY1/2/3`). Caller reshapes. NOT `.density` attributes.
   Reason: SWIG `%attribute` compiles only for **value-returning** getters; it
   FAILS for `std::string` getters AND for `void get_x(double**,int*)` getters.
   There is no `%attribute_np` in this repo.
2. **Scalar fields = `%attribute`** (native, no Python). e.g.
   `%attribute(IMP::bff::States, int, n_points, get_n_points);`
   Caveat: `%attribute` **ignores** the matching `GetMethod`, so the method is
   replaced by the attribute (tests must use `obj.n_frames`, not `get_n_frames()`).
3. **Input grids = shared typemaps in `pyext/IMP_bff.types.i`** now accept ANY
   contiguous numpy array flattened + `None`→empty:
   - `const std::vector<double>&` (already existed; was 1-D-only, I widened it,
     added `None` handling and the any-dim flatten)
   - `const std::vector<int>&` (I added it, mirroring the double one; handles
     int8/uint8/int32/uint32/int64 contiguous + None)
   These replaced every `.ravel()` the old pythoncode did.
4. **Duck-typed fixtures are gone.** `QuenchedDonorDecay(av, atoms)` and
   `DynamicAccessibleVolume` now take **real** `AccessibleVolume`/`ObstacleAtoms`
   values (no `FakeAV`/structured-numpy). Tests build
   `AccessibleVolume(density=..., grid_step=..., attachment_point=...)` and
   `ObstacleAtoms()` default-constructed with fields assigned
   (`.chains=`, `.res_ids=`, `.res_names=`, `.atom_names=`, `.coords=`) — SWIG
   value structs do NOT accept keyword ctor args (only default ctor).
5. **dtype on the caller**: density grids are uint8; the C++ vector<int>/vector<double>
   typemaps accept them flattened, but a mixed call like `AccessibleVolume(density=sp.density.astype(np.float64), ...)` is the pattern (cube → flat is the typemap's).

## Completed files (Phase 1 = 17 files): 14 done, committed

Commit trail (all on `dev`):
- `1cf23b2` — avdistance, interactionterms, forcefield, observables, greedyolga, dye, distributions, griddiffusion
- `7aa112a` — avmodel.i → C++ (scalar `%attribute`, array `get_*`) + direct AV test migrations
- `a1e2c76` — statesdistance.i → C++ (FRETPairGeometry/Efficiencies lose dict/`[...]` surfaces → `get_R()`/`get_E()`/members; LabelDistribution passthroughs; LabelDistributionAV takes atoms_xyz/atoms_vdw)
- `bb51a88` — fret_pair_distribution as C++ free function
- `51d4a7c` — quenchingmodel.i + dyediffusion.i → C++

**Phase-1 remaining (5 files — this is THE next work):**

| File | `%pythoncode` blocks | Notes |
|---|---|---|
| `avbuilder.i` | 1 | **IN PROGRESS RIGHT NOW** — see below |
| `petquenching.i` | 2 | `_pet_table`/QUENCHER_ATOMS dict constants, `residue_sites`/`quencher_atoms` etc. free functions |
| `dyesampling.i` | 3 | `simulate_dye_diffusion`, `simulate_photon_trace`, etc. + ProbeRotamerLibrary/DyeDiffusionTrajectory properties |
| `photophysics.i` | 2 | kappa2/s2 wrappers, `kappa2_from_dipoles` etc. — many reshape-only wrappers |
| `quenching.i` | 1 | `_cube`/`_flat`/`_kappa2` helpers + ~15 kernel wrappers + `radial_diffusion_map` (see handoff) |

Then **Phases 2–5** are far larger: `fps.i`, `cif.i`, `structureio.i`, `rotamer_ensemble.i`, `avmeandistance.i`, `topology.i`, `scoring.i`, `label.i`, `sampling.i`, `rotamer.i`, `docking.i`, `sim.i` (~9,000 lines of dict/JSON/IMP.pmi glue — the heavy lifting is IMP.pmi subclasses and JSON bridges, user wants them in C++/nlohmann).

## IN PROGRESS: avbuilder.i (uncommitted!)

Files currently modified (uncommitted):
- `include/AVBuilder.h`
- `src/AVBuilder.cpp`
- `okf/` log/index + untracked `okf/prds/prd-118.md` touches (not the port)

What I have done for avbuilder:
- Added two **nlohmann-based overloads** in AVBuilder.h, taking the position/paths as **JSON strings** (SWIG-marshallable as `std::string`), parsed internally:
  - `compute_av_from_structure(const std::string& pdb, const std::string& positions_json, const std::string& disc_step="")`
  - `compute_avs_for_structure(const std::string& positions_json, const std::string& pdb_path_or_json, const std::string& disc_step="")`
- Implemented in AVBuilder.cpp: `av_from_position(pdb, nlohmann::json, disc_step)` helper (reads chain/resseq/atom_name/linker_length/linker_width/radius{1,2,3}/strip_mask/allowed_sphere_radius/contact_volume_*; validates `simulation_grid_resolution` vs `disc_step`, throws on disagreement per the old pythoncode), and the two public functions that `nlohmann::json::parse` their strings and call the typed `compute_av_from_structure`/`find_attachment_point`.
- Vendored header `internal/JSON.h` is included in AVBuilder.cpp. The C++ built clean (18/18 linking) BEFORE the throw fix — re-verify after.

## Current state — build status: green

- `ninja IMP.bff` completed clean (turning), but the last run was BEFORE the last `throw` edit. MUST rerun.
- The `.i` (`pyext/IMP_bff.avbuilder.i`) STILL has the old `%rename(_compute_av)` `%rename(_compute_av_from_structure)` and the full `%pythoncode` block — those must be **removed** so SWIG exposes the public C++ names directly. The `%apply(double* IN_ARRAY2,... (double* atoms_xyzr, int n_atoms, int n_cuds))` line stays (needed by the typed `compute_av`).
- After that: test the new surface — THE tests use `compute_av_from_structure(pdb_path, dict(...), disc_step=...)` passing a **dict**, which will now get the JSON-string variant. The user said interfaces change freely, so either (a) migrate the callers/prints the position via JSON text at the call site, or (b) expose BOTH. The existing dict→JSON conversion is the caller's job now (json.dumps) OR add a `%pythoncode`-free map. Recommend: migrate callers/tests to `json.dumps(position)`. Note `compute_avs_for_structure(pdb_path)` in the old pythoncode accepted `(positions, pdb_path, disc_step)` where positions is a dict keyed by name — now `(positions_json, pdb_or_json, disc_step)`.

## Existing commit history to honor
Same style: message describes the move, mentions the interfaces that changed, includes "non-medium suite ... N passed". I have been committing WITHOUT user prompt but the user said "write handoff/log commit" earlier and the change-tracking workflow (okg/log.md dated bullets, commit locally) is sanctioned.

## Log/OKF discipline
- `okf/log.md` has dated entries per change; add one for the avbuilder finish.
- `okf/prds/prd-117.md` is the PRD. There is an untracked `okf/prds/prd-118.md` (not mine).
- The shared stack has `../../chisurf/okf/subsystems/imp-module-convent`. Do not touch cross-repo without checking the change-tracking workflow.

## Memory (mnemosyne)
Private memory id `927802c301c67c3b` has the milestone. A newer memory may exist; recall first.

## Own the tests

The full non-medium suite currently passes EXCEPT the 5 remaining avbuilder-induced failures will appear once avbuilder.i is wrapped. Then proceed to petquenching.i, dyesampling.i, photophysics.i, quenching.i — same pattern: strip `%pythoncode`, move conversion/reshape/validation to C++ with the locked typemaps + `%attribute` scalars + `get_*` arrays, update the small set of failing tests to the flat/method/attribute surface, full suite to green, cut a commit.

**Do not** port the DEFORMATION tests / medium_test_av (ignored/unchanged): they carry the missing-data failure.