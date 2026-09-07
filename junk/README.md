# junk/

Code that is **out of the package** but kept in the repository for reference.

`FASPR/` and `DEERpredict/` are the two things here that were never in `IMP.bff`: third-party clones
(Huang's side-chain packer, `dun2010bbdep.bin` Dunbrack rotamer library) kept
as external reference. See `okf/log.md` 2026-08-21 for why it is here and the
two ways it could feed cgdye (rotamer library source; real-side-chain sampling
for dye–quencher accuracy).

Nothing here is importable, built, installed or tested. It is not a staging
area and not a deprecation shelf — it is the record of what something used to
look like, kept because reading the old version is occasionally the fastest way
to answer "why is it like this".

The rule that governs what lands here, from the PRD-113 cleanup (2026-08-18):

> **Capability never moves to `junk/`.** A file may only be moved here once
> everything it could still do lives somewhere in `IMP.bff`. If a function is
> the only implementation of something, it gets folded into the domain it
> belongs to *first*, and only the husk is moved.

**Amended 2026-08-19.** That rule covers code being *relocated*. It does not
cover code being *retired*, and the NMR restraints are the first of those: the
capability is gone from `IMP.bff` on purpose, by the owner's decision, and
nothing replaced it. Retirement is the one case where the capability does move
here, and it is marked as such in the table below so the two are never confused
— a husk means "look elsewhere in the package", a retirement means "this is the
only copy left outside git history".

What is here, and where its capability went:

| here | what it was | where the capability lives now |
|---|---|---|
| `ported-to-cpp/` | Python whose job moved into C++ | the C++ named in each file's docstring. **Not a retirement**: the capability is in the module, gated against exactly this code. Kept because the edge cases are the part worth re-reading -- what `_calculate_ws` does when every partition function is zero, what `_weighted_average_sd_se` does with a NaN frame |
| `distance_metrics.py` | six functions sharing names with `representation.distance`, three of them genuinely different questions | `IMP.bff.representation.distance` — `polynomial_transfer_ascending`, `distance_sample_statistics`, `mean_position_distance`; the rest were verified identical and deleted |
| `cgdye-scripts/` | 24 one-off driver and analysis scripts, mostly for a single system (`hgbp1_site481`) | the library functions they called are all still in `IMP.bff.cgdye`. **One exception**: writing mol2 existed *only* in `pdb_to_mol2.py`, so it was lifted into `IMP.bff.cgdye.io.mol2` (`write_mol2`, `parse_pdb_atoms`, `parse_conect_bonds`, `infer_bonds`) before the script moved |
| `fret_distance_shim.py` | a 28-line re-export of `representation.distance` left behind by PRD-113 stage 4b | `IMP.bff.representation.distance`, which is where it pointed |
| `cgdye-scripts/tests/` | the two tests that drove junked scripts | nothing — a test whose subject is gone should not stay in the suite skipping itself green |

Three scripts stayed in the package because tests exercise them
(`compare_av_rotamer`, `convert_dye_pdb`, `langevin_hgbp1_site481`).

Two things this cleanup taught, worth repeating before the next one:

* **"Nothing imports it" is not "dead".** `sampling/segments.py` had no
  importers by module path and is reached by the `dye label-fusion` CLI
  command; `write_mol2` had no importers and was the only mol2 writer in the
  repository. Both would have been lost by a reachability scan alone.
* **Grep by import misses tests that locate a file by path.**
  `test_rrt_hgbp1_script.py` found its script through the filesystem, so the
  first scan reported the script as unreferenced. It then *skipped itself*
  green when the script vanished, which is worse than failing.

If you need something from here, the move is to lift it into a domain module
with a test, not to import it.

| `nmr_restraints/` | the NMR-restraint CIF dialect, its reader and writer, the cgdye MD runner's distance-restraint path, and its test | **retired** — no replacement, by decision on 2026-08-19 |
