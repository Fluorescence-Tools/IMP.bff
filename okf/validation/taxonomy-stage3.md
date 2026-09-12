# Probe taxonomy migration: verification and existing failures

**Historical checkpoint:** superseded by [the completed migration](taxonomy-completion.md). The owner overrode the reader hold; the recorded failures are resolved or corrected in that continuation.

The simulation, accessible-volume, template/topology/force-field, scoring and
probe-data file migration follows PRD-141. No numerical algorithm or data
fixture changed. The simulation and AV bodies compare identically after the
approved symbol substitutions and extraction; independent review checked
SWIG macros, plurals, map templates and serialized strings.

## Before changes

- AV: 66 passed, 3 xfailed, 30 subtests (AccessibleVolume, AV lattice,
  AV distance and contact-volume tests).
- Template/topology/IO/data-path/scoring/potentials: 86 passed, 11 failed,
  5 errors. Every failure/error comes from corrupt-brotli exceptions while
  reading potential-table objects. These are present before file renaming.
- Simulation/objects/quenching: 97 passed, 1 skipped, 9 failed, 6 errors.
  Five failures and six fixture errors are `Particle N already set up as
  Mass` in molecular simulation setup. Four failures are stale caller tests
  using `run(t_max=...)`; the existing duration-based method is `simulate`.
  Those four callers were repaired without changing their assertions.

## Static and build checks

The new simulation headers/core sources and the molecular IMP source compile
independently. AV/core file checks exposed missing includes supplied by the
unity build: OutputView for accessible-volume views, fstream for force-field
CIF output, numeric for rotamer scoring, and States for the FPS distance enum.
The changes add those includes. The IMP public all-header and standalone
aggregate (with its normal compatibility preamble) compile under C++17.

Standalone configuration now prunes copied public headers no longer present
in the source tree, preserving generated bff_config.h. Otherwise a reused
build directory would install retired API headers; the copied public tree
was checked to contain no obsolete header after reconfiguration.

## Runtime results

Both IMP and standalone library/binding builds complete. The combined IMP
slice has **643 passed, 1 skipped, 3 xfailed, 30 subtests passed**, plus the
same **16 failures and 11 errors** recorded before editing: 16 potential-table
corrupt-brotli cases and 11 molecular Mass-setup cases. No failing/error case
falls outside those two baseline groups. The four stale `run(t_max=...)`
caller failures are gone. The standalone slice has **242 passed, 6 skipped**.

A new portable zero-diffusion fixture checks both builds' shared interface,
trajectory type, coordinate shape, absent energy arrays and femtosecond time
axis. IMP-only modules are tested under IMP; their direct selection under
standalone bypasses the repository collection filter and is not a supported
standalone test run. Ruff and mypy pass on the added public/simulation guards;
all changed Python files parse. The imp-tricks follow-up is consumer-only: its source tree has no remaining
rotamer implementation. Legacy example imports/value assumptions and dead
CLI registrations are migrated to the existing bff API in the companion
change (`imp-tricks` commit `6ee5385`). Its three durable numerical tests,
CLI help and both example help commands pass on the applied files. Synthetic
placement/weight parity is within 1e-13; a real T4L backbone/ensemble smoke
check also passes. Complete scientific plotting runs are not claimed here.

Reproduce the IMP comparison (arm64 environment, repository root):

```sh
python -m pytest test/label test/test_public_api_names.py \
  test/test_probe_simulation_contract.py test/test_AccessibleVolume.py \
  test/test_av_lattice.py test/representation/test_av_distance_cpp.py \
  test/representation/test_av_contact_volume.py test/cgprobe/test_template.py \
  test/cgprobe/test_topology.py test/cgprobe/test_dye_topology.py \
  test/cgprobe/test_io.py test/cgprobe/test_system_paths.py test/potentials \
  test/scoring test/test_dye_dynamics_doors.py test/test_cpp_objects.py \
  test/quenching/test_quenching_model.py test/quenching/test_quenching_kernels.py \
  -q -o addopts=''
```

For standalone, set `PYTHONPATH` to the standalone build's `python` directory,
`IMP_BFF_DATA` / `IMP_BFF_EXAMPLES` to the checkout directories, and
`IMP_BFF_GPU=off`. Run `test/label`, `test/scoring`,
`test/test_probe_simulation_contract.py`, the probe-pair numerical contract,
and the public-header/retired-name guards.

## Remaining scope

The separate rotamer reader changes remain untouched. The shared ticket
T-20260910-01 must release that work before the probe/side-chain library
split. PTO/compatibility migration also requires editing that library, so
it cannot be declared complete while the ownership gate remains in force.
