# Labelizer → IMP.bff: what became what

The map from the reference Python package
(`../labelizer-backend/src_labelizer/labelizer/`, 17 modules, ~4.2 kLOC) to the
native port (PRD-120). Read this when you have a line of the reference in front
of you and want the C++ that replaced it, or the reverse.

Conventions: everything ported carries the `labelizer_` / `Labelizer` prefix so the ported
model is distinguishable from the module's own physics. Reference line numbers
are as of the checkout at `../labelizer-backend` on 2026-08-24.

## Public files

The C++ API has four headers, each with a matching source file:
`LabelizerFeatures.h` for structure properties, `LabelizerScore.h` for the
labelability model, `LabelizerFRET.h` for pair scoring and its accessible-volume
provider, and `LabelizerIO.h` for scored containers. Include the header for the
operation you call; the Python API remains flat under `IMP.bff`.

## Module by module

| reference module | lines | became | note |
|---|---|---|---|
| `labelizer.py` | 466 | `labelizer_score_structure`, `bin/imp_bff_labelizer` | orchestration; the `Labelizer` class does not survive as a class |
| `labeling_parameter.py` | 267 | `LabelizerTable`, `labelizer_load_table`, `labelizer_lookup`, `labelizer_lookup_key`, `labelizer_parameter_scores` | the base class becomes a table plus a dispatch |
| `labeling_score.py` | 286 | `labelizer_labeling_score`, `LabelizerParameter`, `labelizer_model_paper` | the weighted geometric mean |
| `cysteine_resemblance.py` | 85 | `labelizer_parameter_scores`, tag `cr` | pure lookup |
| `secondary_structure.py` | 158 | **`labelizer_dssp`** + tag `ss` | the DSSP *binary* became C++ |
| `solvent_exposure.py` | 163 | **`labelizer_residue_depth`**, `labelizer_relative_solvent_accessibility`, `labelizer_half_sphere_exposure` + tag `se` | the MSMS *binary* became C++ |
| `conservation_score.py` | 155 | `labelizer_read_consurf` + tag `cs` | import only, as in the reference |
| `methionin_exclusion.py` | 100 | `labelizer_parameter_scores`, tag `me` | |
| `tryptophan_proximity.py` | 234 | `labelizer_parameter_scores`, tag `tp` | ported *working*; raises in the reference |
| `charge_environment.py` | 351 | `labelizer_parameter_scores`, tag `ce` | ported *working*; raises in the reference |
| `fret_score.py` | 978 | `LabelizerFRET.h` in full | the biggest module |
| `measurement_score.py` | 196 | `LabelizerPairScore`, `labelizer_pair_scores*` | the CSV base class is gone |
| `fluorophore.py` | 218 | **not ported** — `ProbeLibrary.h` already had it | see *Not ported*, below |
| `label_lib_functions.py` | 148 | **not ported** — `AVBuilder.h` / `StatesDistance.h` | LabelLib is banned here |
| `pdbhelper.py` | 108 | `labelizer_read_structure` | over `read_pdb_records` |
| `auxiliary_functions.py` | 82 | **not ported** | MSMS binary discovery, plus dead code |
| `config.py` | 147 | `LabelizerOptions`, `LabelizerFRETOptions` | constants become defaults on a struct |
| `__init__.py` | 65 | `labelizer_model_paper`, `labelizer_available_tables` | |

## Function by function, where it is not obvious

| reference | imp.bff |
|---|---|
| `Labelizer.calc_parameter_score()` (`labelizer.py:216`) | `labelizer_parameter_scores` |
| `Labelizer.calc_labeling_score()` (`:246`) | `labelizer_labeling_score` |
| `Labelizer.calc_fret_score()` (`:261`) | `labelizer_fret_pair_scores` / `labelizer_pair_scores_two_states` |
| `FRETScore.calc_measurement_scores(chains_apo, chains_holo)` (`fret_score.py:479`) and the chain filter at `:728` | `LabelizerFRETOptions::donor_chain` / `::acceptor_chain`, plus `::chain_map` for two files that name a chain differently. **Ported 2026-09-07**; it was missing until then, and without it a homodimer cannot be screened at all — see `okf/log.md` 2026-09-07 (1). |
| `Labelizer.zip_files()` (`:423`) | `labelizer_write_pto` — one container, not a zip |
| `LabelingParameter.factory(tag)` (`labeling_parameter.py:56`) | a tag string on `LabelizerParameter`; no class hierarchy |
| `LabelingParameter.calc_score_frequency(v)` (`:184`) | `labelizer_lookup` (numeric) / `labelizer_lookup_key` (categorical) |
| `LabelingParameter._load_logic(m)` (`:198`) | `labelizer_load_table`, cached on name |
| `LabelingParameter._save_pdb/_save_csv` (`:~230`) | **gone** — see *What did not come across* |
| `LabelingScore._labeling_score(...)` (`labeling_score.py:187`) | `labelizer_labeling_score`, weight as repeat count |
| `label_score_model_paper` (`labelizer.py:40`) | `labelizer_model_paper()` |
| `Labelizer.parameter_score_model` (`:72`) | the `table` field of each `LabelizerParameter` |
| `FRETScore.joined_label_score(...)` (`fret_score.py:192`) | `labelizer_joined_label_score` |
| `FRETScore.calc_measurement_score(...)` (`:225`) | `labelizer_pair_score_double` |
| `FRETScore.calc_measurement_score_single(...)` (`:263`) | `labelizer_pair_score_single` |
| the negative control (`:243`) | `labelizer_pair_score_negative_control` — the reference computes and discards it |
| `FRETScore._fret_efficiency(r)` (`:260`) | `IMP::bff::fret_efficiency`, which already existed |
| `FRETScore._get_cb_position(res)` (`:317`) | `labelizer_cbeta_position` |
| `FRETScore._get_gly_cb_vector(...)` (`:293`) | folded into `labelizer_cbeta_position` |
| `FRETScore._calc_simple_mean_position(...)` (`:326`) | `labelizer_alpha_cone_mean_position` |
| `FRETScore._calc_single_av(...)` (`:404`) | `compute_av_from_structure` — the module's own AV |
| `FRETScore._calc_distances_single(...)` (`:641`) | `labelizer_fret_pair_scores`, with the AVs cached per site |
| `FRETScore._calc_heat_map(...)` (`:920`) | `labelizer_cbeta_difference_map` |
| `Fluorophore.calc_foerster_radius(a)` (`fluorophore.py:130`) | `IMP::bff::forster_radius` — already existed, derived not supplied |
| `llf.effDistance(av1, av2, R0, n)` (`label_lib_functions.py:138`) | `model_distance(s1, s2, "RDAMeanE")` |
| `llf.meanAV(av)` (`:122`) | `AccessibleVolume::get_mean_position` |
| `pdbhelper.remove_hetatoms(model)` (`pdbhelper.py:16`) | the `protein_only` argument of `labelizer_read_structure` |
| `ChargeEnvironment.calc_global_charge()` (`charge_environment.py:160`) | `labelizer_global_charge` — net, positive and negative formal charge. Offered by the reference and consumed by nothing there. Its PQR variant, which reads a partial charge out of the **occupancy** column, is not ported: see the header. |
| `llf.saveXYZ/savePqr/saveAV/loadAV` (`label_lib_functions.py`) | `write_path_map` and `AccessibleVolume::get_points` / `PathMap::get_xyz_density`, plus IMP.em's map writers — the capability is here and is not a text format of our own |

## Constants

All of `config.py` moved onto two option structs, so a caller sets them rather
than editing a module-level singleton.

| reference | imp.bff |
|---|---|
| `LS_THRESHOLD = 0.5` | `LabelizerFRETOptions::label_score_threshold` |
| `AV_CALC_RADIUS = 20.0` | `LabelizerFRETOptions::alpha_cone_radius` |
| `COLLISION_RADIUS = 1.7` | `LF_COLLISION_RADIUS` (file-local in `LabelizerFRET.cpp`) |
| `GS_PRECISE = 0.8` / `GS_FAST = 1.2` | `LabelizerFRETOptions::grid_resolution` |
| `N_KALININ_PRECISE = 100000` | no equivalent — the AV distance is deterministic quadrature, not Monte Carlo |
| `SPEED = "precise"` | gone; the caller sets the numbers directly |
| `max_bad_aa_distance = 6`, `min_bad_aa_solvent_exposure = 0.4` | `LabelizerOptions::exclusion_distance`, `::exclusion_exposure` |
| `HSE_RADIUS = 13.0` | `LabelizerOptions::hse_radius` |
| `RADII`/`WEIGHTS` (`tp`, `ce`) | file-local in `LabelizerScore.cpp`, at their use site |
| the `IMG_PIXELSIZE` / `ENDCAP_RANGE` / `R_DIST_*` block | **not ported** — dead in the reference, copied in from another project |

## Not ported, and why

* **`fluorophore.py`** — `ProbeLibrary.h` already models a dye and derives
  \(R_0\) from spectra (`forster_radius`). Porting a second one would have
  given the package two descriptions of a dye, which PRD-113 exists to prevent.
  The reference's linker length, width and three radii are *not* dye
  properties: they are AV parameters, and they live on `LabelizerFRETOptions` here and
  on the `AV` decorator in general.
* **`label_lib_functions.py`** — a wrapper over LabelLib, **banned in imp.bff**
  (owner rule, 2026-08-11). Every call has a native counterpart above.
* **`auxiliary_functions.py`** — MSMS binary discovery (there is no MSMS now),
  plus two classes nothing calls that contain a Python-2 bug.
* **`_save_pdb`** — writing a dimensionless score into a PDB B-factor column,
  a field that means "isotropic displacement parameter, in Å²". The container
  replaces it. This is also the mechanism behind the shipped example's
  conservation feedback loop; see `okf/validation/labelizer_ab.md`.
* The six CSV writers, the heat-map JSON writer and `zip_files` — one
  `.mmfdb.pto` replaces all of them.

## Behaviour that deliberately differs

| | reference | here |
|---|---|---|
| DSSP | external binary; **fails on macOS** (`secondary_structure.py:126`) | `labelizer_dssp`, native. Reproduces the binary exactly on 1DDB (195/195) |
| residue depth | MSMS binary; **does not execute** (32-bit ppc/i386) | `labelizer_residue_depth`, native. No bias, r = 0.976 vs the reference's output |
| accessible volume | LabelLib | `compute_av_from_structure` |
| AV distance | `ll.meanEfficiency`, 100 000 MC samples | deterministic lattice quadrature |
| AV rebuild cost | inner site rebuilt per outer site (`fret_score.py:653`) | placed once per site |
| a score that was not computed | `-1` (excluded) and `0` (no contribution), in the score column | **absent**, with a `status` saying why |
| a sequence gap | array-indexed turns bridge it silently | `labelizer_dssp` refuses to bridge; see `test_dssp_does_not_infer_a_peptide_bond_across_a_sequence_gap` |
| matching sites across two conformations | by residue number (`fret_score.py:512`), so chains collide on a multimer | by `(chain, seq_id)`, with `LabelizerFRETOptions::chain_map` when the names differ. The reference's behaviour is not selectable: it returns `distance_2 = 0` for every inter-chain pair, which is not a modelling choice but a defect |
| the Cbeta difference map on a multimer | indexed by residue number, every chain on the same row | `labelizer_cbeta_difference_map` takes a `chain` |
| a β-bulge | (the reference's DSSP joins the ladder) | joined by union-find over bridges; without it sheets fragment into isolated `B` |
| a `G` shorter than three residues | (the reference's DSSP does not emit one) | suppressed — the span must be free, or the residues fall through to `T` |
| output | 6 CSV + 4 PDB + 1 JSON + zip | one `.mmfdb.pto` |

Three arithmetic differences are **selectable, not fixed** — `LabelizerModel`,
defaulting to the published behaviour. See PRD-120.

## Where a reference defect is reproduced on purpose

| defect | reference | reproduced under |
|---|---|---|
| joined label score uses `prod ** 0.5`, not `prod ** (1/N)` | `fret_score.py:220` | `LABELIZER_MODEL_PUBLISHED` |
| a weight-0 term reading exactly `0.0` zeroes the whole score | `labeling_score.py:166` | `LABELIZER_MODEL_PUBLISHED` |
| the excluded amino acid is always MET whatever is configured | `labelizer.py:233` | `LABELIZER_MODEL_PUBLISHED` |

And two that are **fixed unconditionally**, because they make the reference
raise rather than compute: `tp` and `ce` call their base `set_up` without a
table name (`tryptophan_proximity.py:40`, `charge_environment.py:63`), and the
LabelLib guard tests its import the wrong way round (`fret_score.py:121`, moot
here).

## Data

| reference | here |
|---|---|
| `resources/probabilities/*.json` (29 tables) | `data/labelizer/probabilities/` — verbatim, they are the fitted model |
| `resources/fluorophores/*.fluo`, `*.csv` | **all 15 covered.** 13 were already in `data/rotamer_library/R0/probe_library.cif` under vendor spellings with identical QY and ε (both sets descend from the same upstream data); `Atto532` and `Atto643` were imported from Labelizer's own tables by `utility/import_labelizer_dyes.py`. |
| the fluorophore database as a whole | `data/dyes.mmfdb.pto` (`ProbeContainer.h`), in `_mmfdb_optical_property` / `_mmfdb_spectrum` terms. ChiSurf's `spectra.db` is the large one — 2 288 probes in the same schema. |
| `examples/1DDB/*` | `test/input/labelizer/` — the A/B fixtures |
| `tests/test_methionin_exclusion.py` | the golden 3j0e counts, in `test/label/` |
| Supplementary Data 1, sheet "MalE - LS" | `test/input/labelizer/malE_*_reference.csv` — the β-rich A/B, and the only published output the two-state layer can be checked against |

Fixture provenance and the two caveats that go with it are in
[`test/input/labelizer/README.md`](../test/input/labelizer/README.md).
