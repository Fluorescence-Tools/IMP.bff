# Handover: BFF-native model search for TCSPC and FCS

Date: 2026-09-12
Branch: `independent-core`
Status: safe verified stopping point; native core and first TCSPC/FCS factories build and pass focused tests. ChiSurf migration is deferred.

## Owner decisions

- ChiSurf is the application/car; BFF is the engine.
- Parameters, links, model topology, objectives, fitting decisions, diagnostics and MCTS live in BFF.
- Required model-search scope is TCSPC and FCS, including future global composition.
- No numerical fallback is allowed. If BFF cannot represent a model/objective completely, return unsupported and do nothing.
- No duplicated Python evaluator, optimizer, reward or topology declaration.
- Candidate graphs do not exchange parameter values. Every topology links to one BFF-owned canonical parameter registry.

## Implemented in BFF

### Generic search

`include/ModelSearch.h`, `src/ModelSearch.cpp`

- `ModelSearch`: native PUCT traversal, lazy expansion, deterministic seed, Dirichlet exploration, cycle/collapse handling, cancellation and result bookkeeping.
- `ModelSearchProblem`: model-independent C++ problem boundary with no SWIG director.
- `TabularModelSearchProblem`: callback-free finite reference/test problem.
- `FittingModelSearchProblem`: optimized one-graph fitting search with declarative masks/seeds, native `FitMinimizer`, cached snapshots and rollback.
- `MultiStructureModelSearchProblem`: one canonical parameter registry shared across complete structure-specific objective graphs. It validates exact owner identity, rejects linked followers/vector ports/duplicate IDs/conflicting owners, owns upstream graph nodes, caches canonical snapshots, runs native minimizers and activates the winning topology.
- Initial structures are optimized before root scoring. This matters for a one-topology problem whose only child is terminal.

Public names were cleaned before release: `ModelSearch*`, `FittingModelSearchProblem`, `MultiStructureModelSearchProblem`; no `FitModelSearch*` aliases remain.

### TCSPC factory

`include/TCSPCModelSearch.h`, `src/TCSPCModelSearch.cpp`

- `TCSPCLifetimeSearchFactory` creates canonical lifetime/instrument owner ports.
- Builds complete `TCSPCDecay -> FitChiSquared` graphs for each requested component count.
- Topologies link directly to the same owners; native BIC scoring and add/remove/stop actions.
- `TCSPCLifetimeSearchSpace` retains the problem and graph objects and exposes topology/parameter inspection.
- Current scope is plain multi-exponential TCSPC over a supplied response. FRET, anisotropy, parsed decay, PDDEM, distributed acceptor, mixture, structure and MaxEnt are not implemented here.

### FCS factory

`include/FCSModelSearch.h`, `src/FCSModelSearch.cpp`

- `FCSModelSearchFactory.create_analytical` returns a ready `MultiStructureModelSearchProblem`.
- Complete native graphs cover 2D/3D analytical Gaussian diffusion, one/two diffusion components, and zero/one relaxation term.
- Canonical IDs cover N, baseline, structure parameter, diffusion times/fraction and relaxation amplitude/time.
- Native weighted/Poisson objective selection, fit range/mask validation, BIC metadata and topology actions.
- FCS seeds use the last-lag baseline and half-amplitude crossing; this materially improves root convergence.
- MDF and saturation/kinetics nodes exist elsewhere in BFF but are not yet registered as search topologies.

### Build integration

- `src/Files.cmake`: includes `ModelSearch.cpp`, `TCSPCModelSearch.cpp`, `FCSModelSearch.cpp`.
- `pyext/include/IMP_bff.core.i`: generic ModelSearch API.
- `pyext/include/IMP_bff.tcspcdecay.i`: TCSPC factory.
- `pyext/include/IMP_bff.fcs.i`: FCS factory/config.

## Verification evidence

Environment: `/Users/tpeulen/mambaforge/envs/arm64`
Build tree: `/Users/tpeulen/dev/imp/cmake-build-arm64`

Final build completed successfully:

```sh
ninja -C /Users/tpeulen/dev/imp/cmake-build-arm64 IMP.bff-python -j4
```

Final focused/API gate:

```sh
/Users/tpeulen/mambaforge/envs/arm64/bin/python -m pytest -q \
  test/mcts/test_model_search.py \
  test/mcts/test_tcspc_model_search.py \
  test/mcts/test_fcs_model_search.py \
  test/test_public_api_names.py \
  test/test_base_header.py
```

Result: **85 passed**. `git diff --check` is clean.

Tests cover deterministic/lazy search, cancellation, failed-candidate rollback, same-graph and multi-graph fitting, canonical owner validation, graph lifetime retention, joint objective/shared link handling, TCSPC component topology, FCS topology family, exact native fitting and invalid configurations.

## Current limits and next BFF work

1. Add FCS MDF and full kinetics/saturation topologies using existing native nodes. Define which modes are comparable under one selection score.
2. Extend TCSPC factory with instrument nuisance topology as native structure choices where scientifically valid.
3. Add BFF FRET/anisotropy topology factories only if they remain in the required TCSPC scope; donor/reference links must use canonical owners.
4. Implement a BFF joint problem factory that composes TCSPC and FCS member factories over one canonical registry. It must support local/shared parameters, member masks/ranges/noise models, count each shared parameter/prior once and expose joint plus per-member diagnostics.
5. Add structured BFF capability/refusal records so ChiSurf asks BFF what a model can search instead of declaring topology.
6. Add versioned search-result/session serialization for `.cs.pto` and history: schema, canonical IDs, chosen topology, seed/budget, termination, score terms and validation status.
7. Broaden numerical fixtures and compare MCTS with ordinary fit/multistart at equal evaluation/time budgets.

Do not implement dozens of ChiSurf parameter-group declarations. The registered-model census found that the prior one-static-graph approach could only describe current topology and could not add components or switch equations. The multi-structure canonical registry is the required foundation.

## Deferred ChiSurf migration

Tracked at:

`../chisurf/.omx/plans/chisurf-changes-after-bff-tcspc-fcs-mcts.md`

The ChiSurf working tree contains transitional native bridge/dispatcher/GUI work created before the ownership boundary was finalized. Treat Python topology/group declarations as temporary fixtures. Once the BFF factories and capability API cover the required cases, replace them with thin request/result routing and delete the legacy Python MCTS environment/tree/reward code.

## Resume point

Start with the BFF joint factory or FCS MDF topology. Read this file and the shared agent board, claim a new BFF ticket, inspect the dirty tree before editing, and rerun the 85-test command above before and after the next change. Do not launch concurrent Ninja builds against the same IMP build tree.
