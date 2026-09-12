# Approved flat taxonomy: completion evidence

The owner reaffirmed `.omx/plans/prd-bff-file-taxonomy.md` on 2026-09-12 as
controlling over the later reader-ticket hold. The complete mapping is now
implemented. `test/api_taxonomy.json` freezes 121 retired paths, 91 retired
public identifiers, and the decisions for conditional investigation rows.

## Final structure

- Graph, Fit, Inference, Bayesian, MCMC and Optimization retain the completed
  stage-1 families. Labelizer has independent Features, Score, FRET and IO
  files; generic measurement selection is ProbePairSelection.
- Probe representations and simulations have their canonical files. Probe
  rotamer representation/placement is `ProbeRotamer`; its transfer driver
  is `FRETRotamer`; FPS payloads/workflows are `FPSRotamer`. Probe `.drot`
  libraries and protein-sidechain/Dunbrack libraries have separate headers,
  sources and bindings. The packer is `pack_protein_sidechains`.
- Photophysics owns its spectrum, interaction, quenching, photon and crosstalk
  files. FCS and TCSPC use their acronym casing. The five spectrum-producing
  nodes have separate files; generic Gaussian/Polymer distance names remain.
- `MMFDBProfile` contains only the domain vocabulary. There is no bff PTO
  reader/writer/object facade; consumers call `pto::File` directly. Its
  implementation, FASPR, legacy FPS I/O and Functions are source-private.
- Bridge, command-line, compute-backend, OpenMP and IMP compatibility files
  use the approved names. The compute plugin C ABI stays unchanged.
- States remains generic: its weighted point-cloud/statistical surface does
  not establish probe-only ownership. RmfIO casing remains as explicitly
  deferred in the plan. FftLength and Rebin are internal.

## Verification

- Full IMP suite: **2530 passed, 19 skipped, 3 xfailed, 276 subtests passed**.
- Full standalone suite: **1225 passed, 23 skipped, 203 subtests passed**
  with the complete path/export manifest.
- ChiSurf graph fitting, FCS/TCSPC, crosstalk and MFD consumer slice:
  **170 passed**. imp-tricks rotamer consumers: **3 passed**.
- All **108 public headers** compile independently using actual IMP flags;
  the generated umbrella also compiles. Standalone builds succeed both as
  a unity translation unit and with `IMPBFF_UNITY=OFF`, including RMF.
- All 87 rotamer implementation units and 16 photophysics implementation
  comparisons preserve tokens after approved name/path substitutions and
  relocation. Header dependency repairs add includes only.
- Refreshed DecayConvolution and its Dual companion are verbatim tttrlib
  copies. Old/new double convolution and shifting produced **23,040
  bit-identical values** in the compiled comparison harness.
- All four shipped potential tables read, and a 10,000-value compressed
  table roundtrips. Existing data and serialized input fixtures are unchanged.
- Python syntax/static checks pass on changed Python files and CLI; focused
  mypy checks pass on the new API manifest tests. Nine Python CLI tests pass.

## Corrections needed to verify the complete migration

The old potential reader sent a domain-specific length prefix into a generic
Brotli decoder. Reading stored bytes and using its existing decoder restores
the unchanged format. The frame-to-IMP bridge now updates the Mass decorator
already installed by Atom, preserving the element-table mass instead of
raising on duplicate setup. A regression test pins explicit element overrides.

External chinet reference adapters use chinet's own Node/Port names again;
capability checks and tests now refer to current bff exports. Deterministic
library content is compared independently of ptolib's random container UIDs;
the high-compression fixture is large enough that fixed index overhead cannot
mask the original output-buffer regression. Four existing CLI undefined names
are repaired with standard-library imports and existing bff functions.

## Reproduce

Use the arm64 environment and rebuild the IMP module with `IMP.bff-python`.
Run `python -m pytest test -q -o addopts=''` from this checkout. For standalone,
set `PYTHONPATH` to the build's `python` directory and set `IMP_BFF_DATA` /
`IMP_BFF_EXAMPLES` to this checkout's data/examples (GPU can be disabled with
`IMP_BFF_GPU=off`). The same test directory selects the appropriate lane.

The independent source build is:

```sh
cmake -S standalone -B /private/tmp/prd141-nonunity -G Ninja \
  -DCMAKE_PREFIX_PATH=/Users/tpeulen/mambaforge/envs/arm64 \
  -DIMPBFF_UNITY=OFF -DIMPBFF_PYTHON=OFF
ninja -C /private/tmp/prd141-nonunity -j3 imp_bff
```

## Existing worktree changes

The pre-existing reader correction is preserved as an unstaged hunk in
`src/ProteinSidechainRotamerLibrary.cpp`; the corresponding Dunbrack test/data
registry edits remain separate. Verification uses the working tree with those
existing corrections intact. The taxonomy commit excludes those user-owned
hunks. Unrelated sibling edits and the untracked backend test remain separate.
Full legacy plotting runs and optional licensed datasets are not claimed.

## Concurrent vendor-package work

After the taxonomy builds passed, T-20260912-ptolib-modular replaced the
single-header vendor layout with a separately built source package. That
owner's CMake, forwarding header, implementation shim, vendor tests and
thirdparty package changes are preserved as separate working-tree changes.
The taxonomy commit contains the independently tested private implementation
TU and does not absorb the modularization. The isolated package-integration build also passes: 375 API/container/rotamer
tests with 6 expected skips, plus all four shipped potential-table reads.
This check uses `/private/tmp/prd141-modular-integration` and does not modify
the other task's vendor/build edits. Taxonomy commit: `54b36b4`; downstream
commits: ChiSurf `fafb002e4`, imp-tricks `9cce724`.
