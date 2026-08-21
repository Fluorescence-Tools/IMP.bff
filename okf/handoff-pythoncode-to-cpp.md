# Handoff: %pythoncode to C++ port

**PRD**: [`okf/prds/prd-117.md`](prds/prd-117.md)
**Status**: ready to start. `pyext/src/` is done (353 lines of thin re-exports).
11,607 lines of `%pythoncode` across 29 `.i` files remain.

## Quick start

```bash
# Build
cd /Users/tpeulen/dev/imp/cmake-build-arm64 && ninja IMP.bff

# Test
cd /Users/tpeulen/dev/imp.bff && pytest test/ --ignore=test/medium_test_av.py

# Configure (only if CMakeLists changed)
cd /Users/tpeulen/dev/imp/cmake-build-arm64 && cmake -S ../imp -B .
```

Source files go in `src/` (listed in `src/Files.cmake`, git-tracked, `git add -f`).
Headers go in `include/`. SWIG wrappers go in `pyext/`. The `.i` file loses its
`%pythoncode` blocks but keeps its `%include`, `%template`, `%extend`, `%rename`,
`%feature("shadow")`, and `%attribute_*` directives.

## Architecture: what the %pythoncode blocks do

Every `%pythoncode` block falls into one of these categories. None requires Python.

### Category 1: numpy shims (~2,000 lines)

Reshape flat `std::vector<double>` to expected shape. C++ returns the right
shape natively.

```python
# BEFORE (%pythoncode):
def sphere_points(n):
    return np.asarray(_IMP_bff._sphere_points(int(n)),
                      dtype=np.float64).reshape(-1, 3)

# AFTER (C++): already returns std::vector<std::array<double, 3>>
# Just %include the header, drop the %pythoncode.
```

Files: `petquenching.i`, `quenching.i`, `dyediffusion.i`, `avmodel.i`,
`dyesampling.i`, `griddiffusion.i`, `photophysics.i`, `statesdistance.i`,
`quenchingmodel.i`, `avbuilder.i`, `avdistance.i`, `distributions.i`

### Category 2: IMP API glue (~4,200 lines)

Calls IMP C++ objects through SWIG. C++ can call them directly.

```python
# BEFORE (%pythoncode):
import IMP.atom
import IMP.core
model = IMP.Model()
hier = IMP.atom.read_pdb(str(path), model)
coords = _IMP_bff._structure_coordinates(hier)

# AFTER (C++): directly holds IMP::atom::Hierarchy, calls API
```

Files: `label.i`, `topology.i`, `scoring.i`, `sampling.i`, `rotamer.i`,
`docking.i`, `sim.i`, `structureio.i`, `avmeandistance.i`

### Category 3: dict/JSON bridges (~1,500 lines)

Python dicts to C++ typed objects.

```python
# BEFORE (%pythoncode):
def forcefield_system_from_dict(d):
    sys_ = IMP.bff.DyeForceFieldSystem(str(d.get("name") or ""))
    # ... build from dict ...

# AFTER (C++): DyeForceFieldSystem already exists, parse JSON or take typed args
```

Files: `cif.i` (forcefield_system_from_dict ~150 lines), `fps.i`,
`topology.i` (build_forcefield_system ~500 lines)

### Category 4: Lazy class definitions (~200 lines)

Python classes defined at first access via `_LAZY` dict.

```python
# BEFORE (%pythoncode):
class AVNetworkRestraintWrapper(IMP.pmi.restraints.RestraintBase):
    ...
_LAZY["AVNetworkRestraintWrapper"] = _build_av_network_restraint_wrapper

# AFTER (C++): C++ class, SWIG wraps it, or keep lazy for PMI dependency
```

Files: `avmeandistance.i` (AVNetworkRestraintWrapper),
`structureio.i` (write_rmf, read/write_rotamer_library_rmf),
`rotamer_ensemble.i` (RotamerEnsemble subclasses States)

### Category 5: IMP.pmi subclasses

Two classes subclass `IMP.pmi.restraints.RestraintBase`, which is a Python
class. For JS/SWIG, these either:
- Use the C++ IMP API directly (IMP::core::Restraint)
- Reimplement without PMI
- Stay as lazy %pythoncode if PMI dependency is acceptable

Files: `avmeandistance.i` (AVNetworkRestraintWrapper),
`docking.i` (multiple PMI subclasses), `sim.i` (PMI macros/tools)

## File-by-file inventory

### Phase 1: Trivial numpy shims (17 files, ~1,500 lines)

| # | File | Lines | What to do |
|---|---|---|---|
| 1 | `forcefield.i` | 12 | One `%extend` with `exclusions()` method returning set. Move to C++ or keep as thin shim. |
| 2 | `interactionterms.i` | 17 | Two properties (kQ, rC) and a `_term_rates` helper. Properties already on C++ getters; just expose. |
| 3 | `observables.i` | 18 | `__len__` and `decay()` reshape. `decay()` is the reshape; C++ can return flat, SWIG reshapes. |
| 4 | `greedyolga.i` | 38 | `select_informative_pairs` wrapper: converts numpy arrays, calls C++ rename, returns tuple. |
| 5 | `avdistance.i` | 44 | Three functions: `random_distances`, `density_to_points`, `split_contact_volume_masks`. All reshape. |
| 6 | `dye.i` | 46 | Two dict constants (`DYE_FLRCIF_ITEMS`, `FORSTER_RADIUS_FLRCIF_ITEMS`). Move to C++ static data or keep as Python dict. |
| 7 | `distributions.i` | 98 | ~10 functions, all `np.asarray(_IMP_bff._foo(...))`. Pure reshape. |
| 8 | `griddiffusion.i` | 98 | `_cube` reshape, properties, `run`/`gradient`/`equilibrium` wrappers. |
| 9 | `avbuilder.i` | 109 | `compute_av`, `compute_av_from_structure`, `compute_avs_for_structure`. Dict-to-typed conversion. |
| 10 | `petquenching.i` | 121 | Dict-returning wrappers, `_pet_table` helper, properties. |
| 11 | `dyediffusion.i` | 140 | `_as_cube` reshape, properties, `run` wrapper. |
| 12 | `avmodel.i` | 157 | `_av_flat`/`_av_params` helpers, properties on States/AV/ACV. |
| 13 | `dyesampling.i` | 183 | Properties on RotamerLibrary/DyeDiffusionTrajectory, `simulate_dye_diffusion` wrapper. |
| 14 | `photophysics.i` | 223 | ~15 functions, all reshape + call C++. |
| 15 | `statesdistance.i` | 195 | Properties on FRETPairGeometry/FRETPairEfficiencies, dict surface. |
| 16 | `quenchingmodel.i` | 313 | `_as_obstacles`/`_as_volume` helpers, properties, update methods. |
| 17 | `quenching.i` | 246 | `_flat`/`_cube`/`_kappa2` helpers, ~15 functions. |

### Phase 2: File I/O + dict bridges (2 files, ~1,100 lines)

| # | File | Lines | What to do |
|---|---|---|---|
| 18 | `fps.i` | 185 | JSON round-trip wrappers. Keep as thin Python shims (JSON is Python-native). |
| 19 | `cif.i` | 560 | `forcefield_system_from_dict` (~150 lines, dict-to-typed), template CIF reader/writer (~250 lines, uses `ihm.format`), rotamer library IO (~100 lines). |

### Phase 3: IMP API glue — medium (3 files, ~800 lines)

| # | File | Lines | What to do |
|---|---|---|---|
| 20 | `structureio.i` | 373 | ~15 reshape wrappers + 3 lazy RMF builders. RMF stays lazy Python. |
| 21 | `rotamer_ensemble.i` | 248 | `RotamerEnsemble` subclasses `States`. Must become C++ class. |
| 22 | `avmeandistance.i` | 156 | `estimate_position_uncertainty` wrapper + `AVNetworkRestraintWrapper` (IMP.pmi subclass, stays lazy Python). |

### Phase 4: IMP API glue — large (4 files, ~4,500 lines)

| # | File | Lines | What to do |
|---|---|---|---|
| 23 | `topology.i` | 989 | `build_forcefield_system` and helpers. Pure math + dict assembly. |
| 24 | `scoring.i` | 673 | CHARMM36 table, mask builders, `compute_rotamer_score`, IMP restraints. |
| 25 | `label.i` | 989 | `Label` dataclass, `backbone_frame`, `attach_dyes`, FP domain detection. |
| 26 | `sampling.i` | 1,049 | `LangevinDyeSampler`, `LinkerSampler`, RRT algorithms. |

### Phase 5: IMP API glue — largest (4 files, ~5,300 lines)

| # | File | Lines | What to do |
|---|---|---|---|
| 27 | `rotamer.i` | 1,406 | `load_rotamer_library`, `RotamerFRET`, file I/O + orchestration. |
| 28 | `docking.i` | 1,431 | IMP.pmi docking engine. Heavy PMI dependency. |
| 29 | `sim.i` | 1,490 | IMP.pmi MD/MC runner. Heavy PMI dependency. |

## Build system

### Adding a new C++ file

1. Create `include/ModuleName.h` and `src/ModuleName.cpp`
2. Add `ModuleName.cpp` to `src/Files.cmake` (the `cppfiles` variable)
3. `src/bff_all.cpp` is auto-generated — it includes all `.cpp` files from `Files.cmake`
4. Build: `cd /Users/tpeulen/dev/imp/cmake-build-arm64 && ninja IMP.bff`

### SWIG registration

New C++ headers are `%include`d in `pyext/swig.i-in`. The include order matters:
headers defined earlier are available to headers defined later.

Current order (from `swig.i-in`):
```
types.i → DataPaths.h → pathmap.i → av.i → HierarchyFrame.h →
petquenching.i → quenching.i → fps.i → dye.i → avdistance.i →
BrownianWalk.h → DiffusionSolver.h → PhotonSimulation.h → QuenchedDecay.h →
RotamerEnergy.h → scoring.i → greedyolga.i →
LifetimeSpectrum.h → Mol2IO.h → SequenceAlignment.h → LinkerGeometry.h →
RotamerStatistics.h → Clustering.h → MolecularGraph.h →
forcefield.i → ForceFieldCIF.h → cif.i →
label.i → topology.i → rotamer.i → sampling.i → sim.i → docking.i →
avmodel.i → avbuilder.i → structureio.i →
photophysics.i → interactionterms.i → observables.i →
statesdistance.i → labelingrestraints.i → avmeandistance.i →
dyediffusion.i → dyesampling.i → griddiffusion.i → quenchingmodel.i →
rotamer_ensemble.i → __getattr__
```

### Pattern for porting one file

For each `.i` file being ported:

1. Read the `%pythoncode` blocks. Identify what they do.
2. Check if the C++ header already has the functions. Most do — the `%pythoncode`
   is thin wrappers around C++ calls.
3. If the C++ header already has it: just `%include` the header, remove the
   `%pythoncode`. The SWIG typemaps and `%extend` blocks stay.
4. If the C++ header is missing something: add it to the header and
   implementation, then `%include` it.
5. Build and test.

### What stays in .i files

Even after the port, `.i` files still have:
- `%include` of C++ headers
- `%template` for STL containers
- `%rename` for name conflicts
- `%ignore` for unwrappable methods
- `%feature("shadow")` for keyword argument constructors (SWIG limitation)
- `%attribute_py` / `%attribute_np` for property access
- `%apply` for numpy typemaps
- `%extend` for thin Python-side additions (like `__len__`, `__getitem__`)

What disappears:
- `%pythoncode % { ... }` blocks containing implementation logic

## Special cases

### radial_diffusion_map (quenching.i)

Takes a Python callable `f(distance)`. For C++/JS, redesign to take a
pre-evaluated array:

```cpp
// BEFORE: f is a Python callable
// AFTER: take std::vector<double> values pre-evaluated on the grid
std::vector<double> radial_diffusion_map(
    const std::vector<double>& density, int ng, double dg,
    const std::vector<double>& values);
```

### read_evaluators_json factory (fps.i)

Takes a Python callable `factory(dict)`. For C++/JS, return raw JSON or
use `std::function`:

```cpp
// BEFORE: factory is a Python callable
// AFTER: return the raw JSON string, let caller parse
std::string read_evaluators_json(const std::string& path);
```

### RotamerEnsemble (rotamer_ensemble.i)

Subclasses `States` (a SWIG value type). In C++ this is straightforward
inheritance:

```cpp
class RotamerEnsemble : public States {
    std::vector<double> atoms_;  // (n, n_atoms, 3) flat
    std::vector<std::string> atom_names_;
    // ...
};
```

The Python `from_site` classmethod is complex (~70 lines of orchestration).
Move to a C++ static factory or a free function.

### AVNetworkRestraintWrapper (avmeandistance.i)

Subclasses `IMP.pmi.restraints.RestraintBase`. For JS/SWIG:
- Option A: reimplement as C++ class extending IMP::core::Restraint
- Option B: keep as lazy %pythoncode (acceptable if PMI is available)
- Option C: drop PMI dependency, use AVNetworkRestraint directly

### docking.i and sim.i (Phase 5)

Heavy IMP.pmi dependency. These are the hardest files. They subclass
PMI macros, use PMI tools, and drive MD/MC through PMI. For JS:
- PMI is Python-only. These must be reimplemented without PMI.
- Use IMP C++ API directly (IMP::atom::MolecularDynamics,
  IMP::core::MonteCarlo, etc.)
- Or: keep as %pythoncode if PMI will never be needed in JS

### Lazy RMF builders (structureio.i)

Three functions (`write_rmf`, `write_rotamer_library_rmf`,
`read_rotamer_library_rmf`) are lazy because they import RMF. RMF is
available as C++ headers. Move to C++ and `#include <RMF/*.h>`.

## Testing

### What tests cover

- `test/` — full test suite (801 passed, 1 failed pre-existing)
- `test/test_public_api_names.py` — verifies all public names resolve
- `test/test_lazy_import.py` — verifies `import IMP.bff` does not pull
  optional modules (IMP.pmi, RMF)

### After each phase

1. Build: `ninja IMP.bff`
2. Run: `pytest test/ --ignore=test/medium_test_av.py`
3. Verify: `grep -c '%pythoncode' pyext/IMP_bff.*.i` should decrease

### End state

Zero `%pythoncode` blocks across all `.i` files. All logic in C++ headers
and implementations. The `.i` files contain only SWIG directives.

## Key C++ types available

Headers in `include/` (55 headers):

| Header | What it covers |
|---|---|
| `AVModel.h` | States, AccessibleVolume, ACV |
| `AVNetworkRestraint.h` | AV scoring restraint |
| `AVMeanDistanceRestraint.h` | Mean-position FRET restraint |
| `DyeForceField.h` | DyeForceFieldSystem, FFSite, FFBond, etc. |
| `ForceFieldCIF.h` | mmCIF reader |
| `CifIO.h` | mmCIF writer, template/rotamer IO |
| `Scoring.h` | CHARMM36 LJ, boltzmann, AABB, rotamer score |
| `RotamerEnergy.h` | Pair energy matrix, all-pairs kernels |
| `FRETPair.h` | Pair geometry, efficiencies |
| `OrientationFactor.h` | kappa^2 distributions |
| `ExchangeFRET.h` | Master equation, averaging limits |
| `InteractionTerms.h` | PETTerm, FRETTerm, RadiativeTerm |
| `LifetimeSpectrum.h` | (amplitude, rate) pairs |
| `DyeLibrary.h` | Dye, Spectrum |
| `PETQuenching.h` | Quencher, PETParameters, ResidueQuenching |
| `QuenchingModel.h` | DynamicAV, QuenchedDonorDecay |
| `QuenchingGrid.h` | slow_factor_grid, quenching_rate_grid |
| `QuenchingMap.h` | quenching_rate_map, fret_rate_map |
| `FRETRateTrace.h` | fret_rate_trace |
| `DyeDiffusion.h` | DyeDiffusionSimulation |
| `DyeSampling.h` | simulate_dye_diffusion, photon_trace |
| `BrownianWalk.h` | Brownian walker |
| `DiffusionSolver.h` | Field-picture solver |
| `GridDiffusionSolver.h` | Grid diffusion solver |
| `PhotonSimulation.h` | Photon Monte Carlo |
| `StatesDistance.h` | fret_pair_geometry, histogram_rda |
| `DistanceCalibration.h` | fit_transfer_polynomial |
| `AVDistance.h` | random_distances, density_to_points |
| `AVBuilder.h` | compute_av |
| `StripMask.h` | StripSelection, parse_strip_mask |
| `StructureIO.h` | PDB/MOL2/DCD IO |
| `TrajectoryIO.h` | Trajectory reading |
| `MolecularGraph.h` | Bond graph, angles, torsions, rings |
| `GreedyOlga.h` | Informative pair selection |
| `Distributions.h` | WLC, Gaussian, Poisson |
| `PolymerChain.h` | Chain models |
| `LinkerGeometry.h` | Rotatable DOF |
| `RotamerStatistics.h` | Weighted averaging |
| `Clustering.h` | Leader clustering |
| `Mol2IO.h` | MOL2 reading |
| `SequenceAlignment.h` | Local alignment |
| `FPSIO.h` | fps.json read/write |
| `FPSSchema.h` | fps.json validation |
| `SolventAccessibleSurface.h` | SASA |
| `PathMap.h` | AV solver internals |
| `AVOccupancyMap.h` | Shared obstacle raster |
| `LabelingRestraints.h` | Chi-squared scoring |
| `ModelPrecision.h` | Position uncertainty |
| `DyeDiffusion.h` | Particle-picture walk |
| `DataPaths.h` | Shipped data location |
