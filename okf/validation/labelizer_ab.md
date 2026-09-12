# The native Labelizer against the reference implementation

**Date: 2026-08-24. Two structures.** 1DDB model 39 chain A (mouse BID, 195
residues, all-helical) — the worked example the reference ships, at full double
precision. And **maltose-binding protein in two conformations**, 1OMP (apo) and
1ANF (holo), 370 residues each, from the paper's Supplementary Data 1 sheet
"MalE - LS" (five decimals). MBP is α/β, so it is the case that reaches the
strand assignment; 1DDB cannot.

The port replaces two external programs and one banned library. Whether that
changes a number is a measurement, not an assertion, and this is it. Every
tolerance below was stated before the comparison was run
(`test/label/test_labelizer_ab.py` holds them).

## Result

| term | what changed | agreement |
|---|---|---|
| `cr` cysteine resemblance | nothing — pure table lookup | **195/195 exact** |
| `ss` secondary structure | DSSP binary → native Kabsch–Sander | **195/195 exact**, zero score delta |
| `cs` conservation | nothing — file import + table lookup | **exact to 16 digits**; the shipped reference is not reproducible, see below |
| `se` solvent exposure | MSMS residue depth → native surface | **146/195 bin-exact**; depth agrees to 0.14 Å rms with **no bias** |

## Secondary structure: exact, after one version difference

The reference shells out to DSSP and raises `RuntimeError("Unknown platform")`
on macOS (`secondary_structure.py:126`), so the term could not be evaluated at
all on this machine. `labelizer_dssp` is a native Kabsch–Sander implementation and
reproduces the reference binary's eight-state assignment on **every one of the
195 residues** — the state histogram matches term for term
(`- 64, H 96, S 21, G 3, T 6, I 5`) and the parameter score delta is exactly
zero.

Getting there took two corrections, both real:

1. **Turn extent.** DSSP marks the residues *between* an n-turn's
   hydrogen-bonded pair — `k+1 … k+n-1`, not the donor `k`. Including `k`
   over-assigned `T` and cost nine residues (93.3 % → 97.9 %).
2. **π-helix precedence.** The remaining four disagreements were all `I → H`,
   in one contiguous stretch. DSSP 2 and earlier let a 4-helix win wherever a
   4-turn and a 5-turn overlap; **DSSP 3.0 reversed it** (Touw 2015). Assigning
   `I` before `H` took it to 195/195. The reference's CSVs were therefore
   produced against DSSP ≥ 3, and the port matches that.

The second is worth keeping in mind if these numbers are ever compared against
an older DSSP: the disagreement would reappear, and it would not be a bug here.

## Solvent exposure: no bias, and the residual is quantisation

The published default is MSMS residue depth (`N_SE11_MEAN_SURFACE_DIST`), and
the MSMS binaries the reference ships are 32-bit ppc/i386 Mach-O that do not
execute. `labelizer_residue_depth` builds the solvent-excluded surface natively (the
probe is rolled on the golden-spiral point set already in
`SolventAccessibleSurface.h`, and the unoccluded positions are pulled back by
one probe radius).

Against the reference, at probe 1.4 Å and 590 points per atom:

| | |
|---|---|
| bias | **+0.015 Å** mean, +0.023 Å median |
| scatter | 0.163 Å sd |
| correlation | **r = 0.976**, slope 0.953 |
| within one table bin (0.268 Å) | **92.3 %** |
| bin-exact score | 146/195 = 74.9 % |

There is **nothing to correct**: the bias is a tenth of the scatter and the
scatter is close to the quantisation floor. The reference depth is only
recoverable to half a bin (the published score is binned), which alone
contributes 0.077 Å sd; removing it in quadrature leaves ≈ **0.14 Å** of genuine
disagreement between an analytic and a sampled surface.

Sampling matters more than the probe: 590 points gives 74.9 % bin-exact against
66.7 % at 200. Known limitation: the native surface is the **contact** surface
only and omits the reentrant (toroidal) patches MSMS includes, which is where
the remaining scatter should live — in crevices. Adding them is the obvious
next step if this term is ever tightened.

## Conservation: exact, on an input that is not

`cs` looked wrong — a constant 0.7231 delta on every residue — and it is not.

The shipped `1DDB-39_cs.csv` has **two distinct values across 195 residues**
where the table has ten bins (`cr` has 20, `se` has 10). The example passes
`prot1_cs=[".../1DDB-39_cs.pdb"]`, and labelizer's `_save_pdb` writes its
*output* to that same path — so re-running the example feeds the previous run's
scores back in as conservation grades. The lookup has a **two-cycle**:

```
grade 1.63  →  score 2.3553071957924936  →  written to the B-factor as 2.36
grade 2.36  →  score 1.6322095472510827  →  written to the B-factor as 1.63
```

Each score rounds to the other's grade, so the pair is a fixed point of the
whole read-score-write loop and the example is stuck on it. Those two numbers
are exactly what the shipped CSV contains. Our
lookup reproduces both values **to all sixteen digits**
(`test_the_conservation_lookup_is_exact_to_the_last_digit`).

So the conservation machinery is verified and the shipped reference cannot be
reproduced from any shipped file — the input needed is one generation back and
no longer exists in the distribution. The real ConSurf grades for this entry
(`1DDB-conservationscore-39-A.pdb`, 187 distinct B-factors) correlate with the
reference's implied bins at **r = +0.07**, i.e. not at all; that file is not
what produced the CSVs either.

`test_the_shipped_conservation_reference_is_degenerate` guards this claim so it
cannot rot silently.

## Accessible volumes: characterised, not pinned

The reference calls LabelLib, which is banned in this package (owner rule,
2026-08-11), so the dye cloud comes from `compute_av_from_structure`. There is
no reference output for the pair layer, so this is described rather than
pinned.

The two-tier design is justified quantitatively. Screening all 15 051 pairs of
the 174 labelable sites with the analytic alpha cone takes **0.01 s**; rebuilding
five of them with real accessible volumes changes them substantially:

| pair | cone | AV | Δd |
|---|---|---|---|
| A3–A30 | 1.8931 @ 52.4 Å | **1.3606 @ 47.0 Å** | −5.4 Å |
| A30–A66 | 1.8672 @ 53.0 Å | 1.5119 @ 48.0 Å | −5.0 Å |
| A66–A80 | 1.8292 @ 50.7 Å | 1.5548 @ 48.4 Å | −2.3 Å |

The cone **systematically over-reaches** by 2–5 Å, and because the pair score
peaks sharply at `R = R0` a 10 % distance error is a 30 % score error.

The blunt version: **of the top five pairs by cone score, zero are still in the
top five after refinement.** The cheap screen is a way of choosing which pairs
are worth building volumes for; it is not a way of choosing between them, and
`n_refine` should not be left at zero for real work. Rebuilding those five cost
0.20 s against 0.021 s for screening all 15 051.

One cost difference, deliberate and numerically neutral: the reference rebuilds
the inner site's volume inside the inner loop (`fret_score.py:653`), making its
exact mode quadratic in AV builds. Sites are placed once here.

## Maltose-binding protein: the strand assignment, and two conformations

The gap this page used to name as its biggest — no β-rich reference case — is
closed. The paper's Supplementary Data 1 carries a per-residue table for MBP in
both conformations, and the reference's own `SS` column there holds **72 `E`
and 5 `B`**.

| | 1ANF (holo) | 1OMP (apo) |
|---|---|---|
| `cr` cysteine resemblance | **370/370** | — |
| `ss` secondary structure | **359/369 = 97.3 %** | **360/370 = 97.3 %** |

`cr` is exact on all 370 residues to the last digit the supplement prints
(max \|Δ\| 4.9 × 10⁻⁶, which is the file's own rounding). Two independent
conformations giving the same DSSP agreement is worth more than either number
alone: an implementation tuned to one fold would not do that.

Getting from 94.6 % to 97.3 % found two real defects, both invisible on 1DDB:

1. **β-bulges were fragmenting sheets.** Ladders were joined only between
   strictly adjacent bridges, so a bulge — a step of one residue on one strand
   against several on the other — broke a strand into isolated bridges. Cost:
   seven `E → B` and several `E → S`. Fixed with a union-find over bridges
   joining any two of the same type that advance together with one side
   stepping exactly one (up to five on the other).
2. **`G` was being emitted as a leftover stub.** A 3-10 helix is three
   residues; assigning whatever remained after `H` had taken the rest produced
   one- and two-residue `G` runs at helix C-termini where the reference has
   `T`. Six of them. Fixed by requiring the whole span to be free.

**Ten residues remain**, all interior to strands, and they are *not* an energy
threshold artifact — the published −0.5 kcal/mol cutoff is optimal here
(94.3 % at −0.4, 94.9 % at −0.6, 89.2 % at −0.3). What is left is DSSP's
sheet-level bookkeeping (its two-best-partner rule and sheet grouping) rather
than its energetics, and closing it means reimplementing that bookkeeping.

### The two-state layer does what it is for

No reference output exists for the pair score, so this is a physical check
rather than a comparison — and MBP is the textbook case, closing around its
ligand by a hinge bend. Ranking all 41 041 pairs by the two-state score:

| pair | d(apo) | d(holo) | Δd |
|---|---|---|---|
| A53–A352 | 52.4 Å | 40.4 Å | **−12.0 Å** |
| A29–A219 | 56.9 Å | 43.8 Å | **−13.0 Å** |
| A55–A352 | 54.9 Å | 43.2 Å | −11.7 Å |
| A3–A352 | 62.9 Å | 50.2 Å | −12.7 Å |

Every one of the top twenty moves by more than 5 Å and at least eighteen of
them *contract* — the score is finding the domain-spanning pairs, which is
exactly what it should rank first. **A29 is one of the sites the paper itself
labelled** (Figure 4's `MalE29`). The Cβ difference map peaks at 13.8 Å over
the same motion, is symmetric, and is zero on its diagonal.

## What is not covered

* **Ten interior strand residues on MBP** — see above; sheet-level bookkeeping.
* **The pair layer has no published reference.** The checks above are physical
  (does the closure rank first, is the map symmetric), not comparisons.
* `tp` and `ce` cannot be compared: both crash on load in the reference
  (`tryptophan_proximity.py:40`, `charge_environment.py:63`) and carry weight 0
  in the published model. They are ported, working, and unvalidated.
* Conservation on MBP is not checked — no ConSurf grades ship for it, and the
  1DDB conservation reference is the degenerate one described above.

## See also

* [`okf/labelizer-correspondence.md`](../labelizer-correspondence.md) — the
  module-for-module map from the reference to this port.
* [`okf/prds/prd-120.md`](../prds/prd-120.md) — the decisions.

## Reproducing

```bash
pytest test/label/test_labelizer_ab.py
python examples/labels/plot_labelizer_score.py
```
