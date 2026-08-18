# `pyext/src` — what each directory is

Fifteen packages sit here and they are not all the same kind of thing. This says
which is which, because the tree cannot.

## The domains

These are the public surface. Every name in `api.py` resolves into one of them,
and `test_public_api_names.py` enforces that a name filed under a domain really
lives there.

| | |
|---|---|
| `dye/` | the species: spectra, photophysics parameters, topology, the library |
| `label/` | a dye attached at a site, its linker and frame, and the quenchers it sees |
| `representation/` | where the dye can be — accessible volume, rotamer library, Gaussian, polymer |
| `scoring/` | is a configuration allowed, and how heavily does it count |
| `dynamics/` | how the model is realised in time — Brownian, Smoluchowski, excited-state kMC |
| `photophysics/` | interaction terms and the orientation factor; rate constants |
| `observables/` | the output contract: `(amplitude, rate)` pairs, unconvolved |
| `io/` | the formats — fps.json and its legacy ancestors, structures, templates |
| `restraints/` | scoring against *experimental* structural data (distinct from `scoring/`) |

The conceptual model behind them — representation, scoring, sampling, analysis,
each with a variant per representation — is in
[`okf/architecture.md`](../../okf/architecture.md). It is deliberately not the
directory layout: the stages cut across the representations, so a tree has to
pick one axis and betray the other.

## Older organisation, still load-bearing

Named for where the code came from rather than what it is. They work, they are
tested, and they are on the public surface; they have not been reconciled with
the domains above.

| | |
|---|---|
| `av/` | the accessible-volume builders. Two front doors — arrays, and structure+fps — that PRD-113 stage 3 set out to merge and did not |
| `fret/` | the FRET engine, docking, P(R_DA), the strip tool |
| `quenching/` | the PET model: chemistry, grids, and the site-level model objects |

## Off the domain layout

| | |
|---|---|
| `cgdye/` | explicit all-atom dye modelling under a force field. Molecular mechanics, which is IMP's own territory. Importable by module path, **not** a domain, **no** flat names. Nothing in the package imports it |

## Utility

| | |
|---|---|
| `tools/` | paths to shipped data, and IMP helpers |
| `cli/` | what `bin/imp_bff` imports. Not the whole command line — see its docstring |

## The one loose file

`api.py` — the public surface. `BY_DOMAIN` is authored, `EXPORTS` is derived
from it. Nothing else lives at this level, on purpose.
