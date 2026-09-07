---
okf_version: "0.2"
title: "IMP.bff against FPS's published screening numbers (T4L, 3GUN)"
description: "End-to-end A/B of the screening path against the real FPS, and the accessible-contact volume that explained the +2 A offset against the published table. Under the default (IMP/CHARMM) radii: 24 of 33 pairs at +0.22 A bias, 0.91 A rmsd, r = 0.9960. Olga's radii, selectable, reach all 33 at -0.02 A / 0.71 A."
tags: [validation, fps, screening, av, acv, radii, prd-121]
---

# IMP.bff against FPS's published screening numbers

**Date:** 2026-09-01 (first written 2026-08-31) · **Drivers:**
[`test/representation/test_fps_screening_ab.py`](../../test/representation/test_fps_screening_ab.py),
[`test/representation/test_av_contact_volume.py`](../../test/representation/test_av_contact_volume.py),
[`test/representation/test_olga_vdw_radii.py`](../../test/representation/test_olga_vdw_radii.py)
· **PRD:** [121](../prds/prd-121.md)

## What is being compared

Zenodo 3376527 (Sanabria *et al.*, *Nat. Commun.* **11**, 1231, 2020) holds what
was published for **421 T4L structures**: 33 donor–acceptor pairs of ⟨R_DA⟩
each, plus χ² for three states, beside the `FRET_screening.fps.json` that names
the labelling sites. One of the 421 is **3GUN**, which this module already
ships — so the whole path can be measured against that answer without shipping a
new structure.

This is the complement to the single-volume pin in
[`test_fps_av_parity.py`](../../test/representation/test_fps_av_parity.py). That
one localises a disagreement; this one is where a **systematic** one cannot
hide. It is the one that found a bias, and then found what caused it.

## The radii, and why the default is not the one that fits best

Three radii sets are in circulation here and they are not interchangeable:
FPS's element-keyed **Bondi** (C 1.70, N 1.55, O 1.52, P 1.80, S 1.80 Å),
IMP's **CHARMM-derived united-atom** set from `read_pdb` (carbon 1.85–2.275 Å,
carrying implicit hydrogens), and **Olga's**, keyed by *atom name* rather than
element. An accessible volume must use exactly one of them.

**The default is IMP's** — the radii already on the particles
(`AV::set_radii_source`, spelled `"imp"`). The owner's decision, 2026-09-01:
*"use the IMP radii, charmm, as otherwise wont be consistent with IMP
docking."* The reason is the other half of a docking score. Excluded volume is
`clash_container`, which measures overlap through `IMP::core::XYZR` — the
**particles'** radii. Inflate the accessible volumes by Olga's table while the
clash term is built on IMP's, and the two halves of one score disagree about
how big an atom is, silently, in a quantity neither half reports. That is the
"one quantity, two places" failure this repository is organised against, and
consistency with IMP docking was judged worth more than agreement with Olga's
published numbers.

**The default is not free, and this page will not pretend otherwise.** With
IMP's radii:

* **nine of the 33 published ⟨R_DA⟩ have no model value at all** — the volumes
  at those sites (all the pairs involving site 132) come back empty, so every
  row below marked "24 pairs" is over a subset of the published table;
* the agreement with the Zenodo table is **worse**: +0.22 Å / 0.91 Å rmsd over
  24 pairs, where Olga's radii give −0.02 Å / 0.71 Å over all 33.

That is the price of internal consistency. **Anyone reproducing Olga-era
numbers should set `"radii_source": "olga"` on the positions** (schema 1.6,
per position); Olga's table stays vendored (`data/olga_vdw_radii.csv`,
`IMP::bff::olga_vdw_radius`) for exactly that.

Olga's radii were the default for part of 2026-09-01, before this reversal.
Every number below was re-measured under the current default rather than
carried over.

## Result

Two things changed on 2026-09-01 and the table has to be read with both of them
named: the accessible-contact volume the file asks for stopped being ignored
(PRD-121 **G9**), and Olga's radii became available as a per-position choice.
Every row is labelled with its radii set; the rows marked **default** are what
this module reports out of the box.

| | radii | pairs | bias | rmsd | r |
|---|---|---|---|---|---|
| **this module, no ACV, vs FPS's own routine** | **IMP (default)** | **24** | **−0.31 Å** | **1.20 Å** | **0.9926** |
| this module, no ACV, vs FPS's own routine | Olga | 33 | +0.02 Å | 1.62 Å | 0.9777 |
| *…the same 24 pairs as the default row* | Olga | 24 | −0.62 Å | 1.13 Å | 0.9940 |
| **this module, no ACV, vs the published table** | **IMP (default)** | **24** | **+2.05 Å** | **2.54 Å** | **0.9901** |
| this module, no ACV, vs the published table | Olga | 33 | +2.05 Å | 2.40 Å | 0.9890 |
| FPS's own routine vs the published table | FPS Bondi | 33 | +2.37 Å | 2.73 Å | 0.9735 |
| **this module, ACV on, vs the published table** | **IMP (default)** | **24** | **+0.22 Å** | **0.91 Å** | **0.9960** |
| this module, ACV on, vs the published table | Olga | 33 | −0.02 Å | **0.71 Å** | 0.9954 |
| *…the same 24 pairs* | Olga | 24 | −0.23 Å | 0.72 Å | **0.9966** |
| **this module, ACV on, vs FPS's own routine** | **IMP (default)** | **24** | **−2.15 Å** | **2.52 Å** | **0.9868** |

`prototypes/fps_oracle/` runs FPS's own routine on the same structure with the
same file, pinned in `test/references/fps_screening_oracle_pins.json`, so both
implementations can be measured against the `.dat` and against each other.

**Three things are established.** Without a contact volume this module and FPS
agree with each other and *neither* reproduces the table. With the contact
volume the file asks for, this module reproduces the table — 24 of its pairs to
0.91 Å rmsd by default, all 33 to 0.71 Å if Olga's radii are selected — with no
fitted parameter anywhere in the path. And the radii set is worth **more than
the search convention**: switching it moves the agreement with the published
table by more than the residual against FPS's own routine.

**The comparison against FPS deserves reading twice.** Olga's radii give a
*smaller* bias against the FPS oracle (+0.02 Å against −0.31 Å) and a *worse*
rmsd and correlation (1.62 Å / 0.9777 against 1.20 Å / 0.9926) — but that is a
different, larger set. Over the **same** 24 pairs Olga's radii are better on
both (−0.62 Å, 1.13 Å, r = 0.9940). The nine pairs they add are all pairs with
site **132**, which united-atom radii bury, and they are the hard ones. Both
figures are pinned in
`test_fps_screening_ab.py::test_rda_agrees_with_fps_itself` so neither can be
quoted without the other, and so that the default is not mistaken for the more
accurate choice. It is the more *consistent* one.

## The unit that had to be right

**Olga's numbers are nanometres upstream and Angstrom here, and the factor is
Olga's own, not an inference from the magnitudes.**
`Olga/src/AV/Position.cpp:118-124` (`coordsVdW`) scales the coordinate *and*
the radius by the same `10.0f` on one line, because pteros stores coordinates
in nm and `calculateAV()` works in Å:

```cpp
xyzw.emplace_back(frame.coord.at(i)[0] * 10.0f,
                  frame.coord.at(i)[1] * 10.0f,
                  frame.coord.at(i)[2] * 10.0f,
                  pterosVDW(system, i) * 10.0f);
```

So carbon 0.17 nm is **1.70 Å**, which is Bondi's carbon and FPS's. Pinned in
`test/representation/test_olga_vdw_radii.py::test_the_unit_is_nanometre_times_ten`.

**The fallback for an unknown atom name is flat.** `pterosVDW` is
`vdWRMap.value(name, 0.15)` (`Position.cpp:110`) — 0.15 nm = **1.50 Å** for any
name the table misses, with no element lookup behind it. It is not
near-transparent: at 1.50 Å it is within 0.01 Å of Olga's own oxygen. But it is
smaller than every heavy atom in the table and much smaller than the
united-atom carbon an unrecognised carbon would otherwise have taken, so it
matters where the name is a carbon. On the shipped fixtures: T4L 3GUN, **no**
misses; the HIV-RT DNA, 85 atoms of 1018 (`C1*`–`C5*`, `O1P`, `O2P`, `O3*`–`O5*`,
`C7` — the old spellings of names Olga writes `OP1`/`OP2`/`C1'`). That is
reproduced rather than corrected: Olga reading the same file would do the same.
It is also why the HIV-RT `resolved` score differs between the two sets
(59.0404 by default, 59.2547 under Olga's table).

**What the two sets do to the volumes.** T4L 3GUN, grid 0.5 Å, clearance 2.0:

| site | IMP united-atom (default) | Olga | mean moves |
|---|---|---|---|
| 132 CB | 82 418 voxels | 95 836 (+16.3 %) | 1.05 Å |
| 55 CB | 136 823 | 144 032 (+5.3 %) | 0.38 Å |
| 99 CB | 710 | 29 586 (×42) | 17.9 Å |

Site 99 is the one to read: at 710 voxels it is a numerical accident, not a
volume. Under the default it stays one.

## The authored clearance buries the structure — and that is the radii

Two statements, both true, and the second is the one that was nearly lost:

**1. Under the default radii, the authored clearance buries the structure.**
With `allowed_sphere_radius: 2` as `FRET_screening.fps.json` writes it, eight
of seventeen sites come back empty and only **9 of 33** pairs have a model
value, where deriving the clearance leaves one empty and gives **24**. FPS
seeds from `LinkerInitialSphere × linker_width` **unconditionally**, where this
module's clearance is a sphere carved out of an obstacle map already inflated
by half the linker width. Any fps.json carrying a small explicit clearance is
asking for something other than it appears to.

**2. What empties those volumes is the radii, not the clearance.** Select
Olga's table and the *same* authored clearance resolves **all 33** pairs, the
same as deriving it, and the two agree to 0.03 Å of bias against the published
table. A carbon of 2.275 Å walls a source in where one of 1.70 Å does not.

The two clearance conventions really are different quantities — that has not
been retracted. What was retracted is the claim that the difference is what
empties the volumes. Both behaviours are asserted, each against its own radii
source, in
`test_fps_screening_ab.py::test_the_authored_clearance_buries_the_structure_and_it_is_the_radii`,
so the attribution stays checkable rather than remembered. It survived the
default moving to Olga's radii and back on the same day, because it is a
statement about the two sets and not about which of them is default.

## The +2 Å was the accessible-contact volume after all

An earlier version of this page said the offset was "in the reference, not the
code", on the grounds that FPS has no ACV and so the ACV could not explain a
disagreement between FPS and the table. The first half is right — **FPS has no
contact volume at all** — and the second half is a non sequitur, because it
assumes the table was produced by FPS. It was not.

`FRET_screening.fps.json` asks for `contact_volume_thickness: 3` and a fitted
`contact_volume_trapped_fraction` at every one of its seventeen sites. Those
fields were accepted by `IMP::bff::AV` and did nothing (PRD-121 **G9**). Wiring
them up (2026-09-01) moves this module from +2.05 Å to **+0.22 Å** of the
published ⟨R_DA⟩ and cuts the rmsd from 2.54 Å to 0.91 Å. The conclusion is the
opposite of the earlier one: the file **is** the parameter set the table was
computed with, and the program that computed it was ACV-capable — Olga, not
FPS. Selecting Olga's radii on top of that closes the rest and widens it to all
33 pairs (−0.02 Å / 0.71 Å), which is corroboration of the same conclusion from
a second direction.

**The two references are from two eras, and the port matches each on its own
terms.** The accessible-contact volume did not exist yet when FPS was written
(owner, 2026-09-01); it is an Olga-era addition. So there is no single "correct"
answer to compare against here, and no contradiction between the two rows of the
table above:

* the **FPS oracle** is pre-ACV by construction, and the comparison against it
  drops `contact_volume_*` exactly as it drops `allowed_sphere_radius` — keys
  FPS has no concept of — giving −0.31 Å;
* the **Zenodo table** is Olga-era and asks for an ACV, and honouring that ask
  gives +0.22 Å.

A file gets an ACV if and only if it asks for one, which is what keeps both
comparisons honest. Reading the offset against one reference as evidence about
the other is the mistake this page made in its first version.

It holds under the other clearance convention too. With the authored
`allowed_sphere_radius: 2`, which leaves only nine pairs resolvable under the
default radii, the ACV moves those nine from +1.70 Å / 2.35 Å rmsd to
**+0.05 Å / 1.06 Å**.

The alternatives that were tried and do not fit are still worth recording,
because they are what makes the ACV explanation load-bearing rather than
merely available:

| statistic (FPS's own clouds, no ACV) | bias vs the table |
|---|---|
| ⟨R_DA⟩ | +2.03 Å |
| ⟨R_DA⟩_E | +3.11 Å |
| R_mp | −1.40 Å |

and recomputing with FPS's own AV3 dye radii from `Fps/data/linker.txt` (donor
5/4.5/1.5, acceptor 11/3/1.5) makes it **worse**, +6.36 Å.

## What the contact volume is, and one number that pins its discretisation

The rule is Olga's, taken from `Olga/src/AV/fretAV.cpp` (`path2points`): a cloud
voxel is *in contact* when excluded volume — the obstacle raster inflated by the
dye radius — lies within `contact_volume_thickness` of it, and the contact
voxels are then weighted so they carry `contact_volume_trapped_fraction` of the
cloud's total weight. Olga rounds the layer down to whole voxels first
(`deltaIlist` takes an `int`), which at `thickness = 3` and a 2 Å grid is a
one-voxel shell rather than a 3 Å sphere.

That truncation is reproduced rather than corrected, and the table says why: a
true 3 Å sphere gives **+0.49 Å bias, 1.04 Å rmsd** against the published
numbers where Olga's shell gives **+0.22 Å, 0.91 Å**. The trapped fraction is a
*fitted* number, and it is only meaningful against the region it was fitted
for. The consequence to know is that the contact layer is quantised by the
grid: a thickness below `simulation_grid_resolution` is no layer at all.

Note that the fractions were fitted against **Olga's** obstacle set, and the
default here is not that set. That is the sharpest form of the price above: a
fitted constant is being applied over a slightly different geometry. It is also
why `radii_source` lives on the *position* rather than on the run — the file
that carries the fitted fraction is the file that should name the radii it was
fitted with.

## What it does to the shipped T4L file

`examples/structure/T4L/fret.fps.json` asks for the same thing — thickness 3 at
every site, trapped fractions 0.33–0.72 — so every number this module reports
for it moved on 2026-09-01 when the ACV was honoured.

| | before | after (default radii) |
|---|---|---|
| mean ⟨R_DA⟩, 33 pairs of `chi2_C1_33p` | 46.91 Å | **43.88 Å** — −3.04 Å, all 33 shorter (−5.46 … −1.40) |
| model vs the file's experimental distances (`chi2_C2_33p`) | +3.93 Å bias, 5.78 Å rmsd | **+0.84 Å**, 3.93 Å |
| `ProbeNetworkRestraint` score, `chi2_C2_33p` | 22.5277 | **11.3268** |

Selecting Olga's radii moves the same three to 43.48 Å, +0.44 Å / 3.91 Å and
11.2866 — a further improvement, and not the default, for the reason at the top
of this page.

The mechanism is visible one site at a time: the fitted trapped fractions
(0.33–0.72) are all well above the *geometric* contact share of the cloud
(0.13–0.47, mean 0.20), so the surface layer is up-weighted. The mean position
moves 2.3 Å on average (0.04–4.6 Å) and **closer to the nearest protein atom at
sixteen of the seventeen sites**, by 1.81 Å.

## A second defect behind the first

Wiring the fields up immediately found that `AV::set_av_parameter()` read the
fraction as `j.value("contact_volume_trapped_fraction", -1)` — an **`int`**
default, so nlohmann deduced `int` and every fps.json fraction was truncated to
0 or 1. It could not be seen while the value was unused. The first measurements
made after the wiring were taken with it, and said the ACV *lengthened* ⟨R_DA⟩
by +1.47 Å; with the truncation fixed it shortens them by 3.04 Å. Both are in
this page's history for the reason such things usually are.

## What survives from the first pass

**1. The authored clearance buries the structure**, under the default radii —
see the section above, which also says what actually causes it.

**2. FPS is not the ACV.** The A/B against FPS's own routine must be run with
the contact volume **off**, because FPS has none. `_scored()` in the test drops
`contact_volume_*` for the same reason it drops `allowed_sphere_radius`: they
are keys FPS has no concept of, and leaving them in turns the A/B into a
comparison of two different models. With them off the −0.31 Å / 1.20 Å /
r = 0.9926 agreement with FPS is unchanged by any of this.

## The other half of the score: the clash term's own radii

Everything above is about the **accessible volume's** radii. The other half of
a docking score is the excluded volume, `clash_container`, and until 2026-09-01
it had no radii control at all — it read `IMP::core::XYZR`, the radii the
particles carry. That is why moving the *volume* to Olga's table left the clash
term exactly where it was (PRD-121 G2: 110.98 of 150.69 before, 110.85 of
140.58 after).

FPS's clash constant is `k_clash = 2/ClashTolerance²` and FPS applies it to
**Bondi** radii. Applied to IMP's united-atom radii, the same interface reads
as a far worse clash. Measured statically at the input pose over HIV-RT's
protein–DNA interface (8016 protein atoms against 1018 DNA atoms) at k = 8,
which is FPS's `ClashTolerance = 0.5` for refinement and error estimation:

| clash radii | overlapping pairs | total overlap | energy |
|---|---|---|---|
| `imp` (united-atom, the default) | **268** | **90.49 Å** | **198.40** |
| `olga` (Bondi scale) | **30** | **7.59 Å** | **11.11** |

A factor of **17.9** in energy, from the radii alone — the coordinates, the
constant and the pose are identical.

`DockingParameters::clash_radii_source` and `clash_radii_scale` (and
`--clash-radii-source` / `--clash-radii-scale` on `imp_bff_fps dock` and
`score`) are that control. The source picks the table, the scale multiplies it;
the defaults are `"imp"` and `1.0`, so nothing moves unless asked. The radii are
carried by **shadow spheres** — new particles at the atoms' coordinates, in the
same rigid bodies — rather than written onto the structure, so a volume built
afterwards still sees the radii its own `radii_source` asked for. FPS's D12
(`SpringEngine` mutating session-global `Molecule` objects and never reverting)
is not reproduced.

### One global scale is not defensible, and here is the measurement

The obvious cheap alternative is a single multiplier on IMP's radii. It does
not work, and the reason is that IMP's ratio to Bondi is **per element**:

| element | n | IMP min/mean/max (Å) | Bondi | Bondi / IMP mean |
|---|---|---|---|---|
| C | 5687 | 1.700 / 2.1033 / 2.275 | 1.70 | **0.8083** |
| O | 1764 | 1.700 / 1.7106 / 1.770 | 1.52 | 0.8886 |
| N | 1520 | 1.700 / 1.8495 / 1.850 | 1.55 | 0.8381 |
| P | 48 | 1.700 / 2.1406 / 2.150 | 1.80 | 0.8409 |
| S | 15 | 2.000 / 2.0000 / 2.000 | 1.80 | 0.9000 |

**The carbon answer is 0.8083** — 1.70 / 2.1033, the population mean over this
structure. It is not the answer for oxygen (0.889) or sulfur (0.900), and IMP's
carbon is not one number: it spans 1.70–2.275 Å where Bondi's is a single 1.70,
so no factor reproduces even carbon exactly.

Bisecting the scale against each of the three interface statistics gives three
different answers:

| statistic matched | scale | pairs | overlap | energy |
|---|---|---|---|---|
| Olga's 30 pairs | **0.8206** | 30 | 5.05 Å | 5.45 |
| Olga's 7.595 Å overlap | **0.8386** | 43 | 7.59 Å | 8.83 |
| Olga's energy 11.114 | **0.8476** | 50 | 9.27 Å | **11.11** |

At the pair-count match the energy is off by a factor of two; at the energy
match the pair count is off by 1.7×. A scale can reproduce **one** number of
the contact distribution and not the distribution. So: **use the scale to
loosen a clash term; use the source to reproduce FPS.** The two compose, and a
caller who wants the closest single-knob approximation of FPS's hard sphere on
this structure should use the source, not 0.83.

### What it does to the two doors

Scoring (`imp_bff_fps score`, k = 1, HIV-RT `resolved`, full AVs). The
restraint half is byte-identical in every row, which is the check that this
knob is the clash term's alone and does not reach the volumes:

| clash radii | score | of which clash | score − clash |
|---|---|---|---|
| `imp` ×1.0 (default) | **59.0404** | 24.8005 | 34.2399 |
| `imp` ×0.90 | 38.8090 | 4.5691 | 34.2399 |
| `imp` ×0.8476 | 35.6303 | 1.3904 | 34.2399 |
| `imp` ×0.80 | 34.6113 | 0.3714 | 34.2399 |
| `olga` ×1.0 | 35.6291 | 1.3893 | 34.2399 |

(`olga` at k = 1 gives 1.3893 = 11.1141/8, the static census divided by the
constant — the two measurements are the same number reached two ways.)

### The G2 bootstrap, before and after

This is the pass/fail the control was built for. PRD-121 G2 records that
HIV-RT's bootstrap collapses to **0.000 ± 0.000 Å** because the clash term
carries most of the parent score. Re-measured on the same cheap configuration
the test runs (`fps_error_estimation_parameters`, score set `resolved`,
`n_frames = 60`, `coarse_clash` off, 3 replicas, seed 3), varying only the
clash radii:

| clash radii | parent score | of which clash | clash / score | RMSD spread (Å) |
|---|---|---|---|---|
| `imp` ×1.00 (**default**) | 150.6884 | 110.9848 | **0.737** | **0.0000 ± 0.0000** |
| `olga` ×1.00 (FPS-like) | 35.8845 | 1.9295 | 0.054 | **0.0689 ± 0.0155** |
| `imp` ×0.85 | 43.2916 | 2.7068 | 0.063 | 0.0983 ± 0.0365 |
| `imp` ×0.80 | 39.9571 | 0.8630 | 0.022 | 0.4040 ± 0.0953 |
| `olga` ×0.90 | 39.7326 | 0.5171 | 0.013 | 0.3539 ± 0.2444 |
| `ev_weight = 0` | 24.3351 | 0.0000 | 0.000 | 2.1297 ± 1.4584 |

**Read this carefully; it is a partial result.** The collapse *is* resolved in
the strict sense — at FPS's own hard sphere the spread is no longer identically
zero, and the clash term stops dominating (73.7 % of the parent score → 5.4 %).
But 0.069 Å is still **thirty times** smaller than the 2.13 Å the same run
gives with the excluded volume switched off. Moving the clash term to Bondi
radii removes 98 % of the clash *energy* and does not remove the pinning. The
spread rises monotonically with the remaining clash fraction and only reaches
~0.4 Å once the clash term is down to 1–2 % of the score, which is well past
any radii set FPS would recognise.

So the honest statement is: **the radii artefact was real and is now
controllable, and it was not the whole story.** Whatever else holds this pose
is not the size of an atom. The `ev_weight = 0` row remains the number to
compare a bootstrap against before believing it, and `parent_e_clash` remains
the number to read first.

## Provenance, stated rather than inferred

The Zenodo data was produced **with Olga, using an accessible-contact volume**
(owner, 2026-09-01). This page originally reasoned the other way round — the fit
improved when the ACV was honoured, therefore the producer was ACV-capable — so
it is worth being explicit that the fit is now **corroboration of a known
provenance**, not the evidence for it. The two facts that looked contradictory,
"FPS has no ACV" and "the ACV closes the offset", are chronological: the ACV
postdates FPS.

## What this does *not* say

The 3GUN row is one of 421; the other structures are not in this checkout, so
"reproduces the published table" is one structure's pairs, not the study — and
by default it is 24 of 33 of them, not all. The −0.31 Å / 1.20 Å residual
against FPS without the ACV is not zero either: it is the two search
conventions (this module's tight stencil against FPS's three-voxel hop) and the
clearance mapping, both open in PRD-121. Nothing here validates the trapped
fractions themselves — they are fitted per site, and this measures that applying
them reproduces the numbers they were fitted to produce, not that they are
right. And nothing here says the default radii are the *better* ones for
reproducing Olga; the table says plainly that they are not. They are the ones
the rest of an IMP docking score is computed with.

## Reproducing

```bash
python -m pytest test/representation/test_fps_screening_ab.py \
                 test/representation/test_av_contact_volume.py \
                 test/representation/test_olga_vdw_radii.py -q
```
