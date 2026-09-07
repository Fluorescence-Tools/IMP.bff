# Labelizer A/B fixtures

Where each of these came from, so a later reader can tell a reference from a
convenience and can replace one without guessing what it was for.

## Structures

| file | what | source |
|---|---|---|
| `1DDB-39.pdb` | mouse BID, NMR model 39, chain A, 195 residues, **all-helical** | the reference package's own worked example, `labelizer/examples/1DDB/`. Trimmed to model 1 — the reference reads `structure[0]` unconditionally, so the other 19 models are unused and cost 4.6 MB. |
| `1anf.pdb` | maltose-binding protein, **holo** (maltose-bound, closed), 370 residues, α/β | the reference package's `test_data/` |
| `1OMP.pdb` | maltose-binding protein, **apo** (open), 370 residues | RCSB |
| `3j0e.pdb` | three chains (F, G, H), 1010 residues | the reference package's `test_data/` |
| `1DDB-conservationscore-39-A.pdb` | ConSurf grades in the B-factor column, 187 distinct values | the reference package's example directory |
| `1DDB-39_cs.pdb` | **not** a conservation input despite the name — the reference's *output*, whose B-factors are its own conservation scores | see the caveat below |

## Reference output

| file | what | source |
|---|---|---|
| `1DDB-39_{cr,cs,se,ss}.csv` | per-residue parameter scores, full double precision | produced by the reference implementation; shipped in its example directory |
| `malE_1ANF_reference.csv` | per-residue CS/SE/CR/SS/LS for the holo conformation, 370 rows | Gebhardt *et al.*, *Nat. Commun.* **16**, 3305 (2025), Supplementary Data 1, sheet "MalE - LS", columns `*.1` |
| `malE_1OMP_reference.csv` | the same for the apo conformation | same sheet, unsuffixed columns |
| `malE_two_state_reference.csv` | the published combined two-state score and its delta | same sheet, `combined` and `Labeling Score` columns |

The MalE tables are printed to **five decimals**, so comparisons against them
are at `1e-5`. That is the file's precision, not a tolerance on the port: `cr`
agrees on all 370 residues to the last digit printed.

## Two caveats worth knowing before using these

**`1DDB-39_cs.pdb` is contaminated, and deliberately kept anyway.** The
reference's example passes it as the conservation *input* while the reference's
`_save_pdb` writes its conservation *output* to the same path, so every run
feeds the previous run's scores back in as grades. The lookup has a two-cycle,
which is why `1DDB-39_cs.csv` holds exactly two distinct values across 195
residues where the table has ten bins. It is kept because reproducing the
reference means feeding it what the reference was fed; the honest conservation
input for that entry is `1DDB-conservationscore-39-A.pdb`, whose grades
correlate with the shipped CSV at r = +0.07, i.e. not at all.

**1DDB has no β structure.** 96 H and not one strand, so nothing in it reaches
the DSSP bridge, ladder or bulge logic. That is what the MalE fixtures are for,
and it is why they were added: two real defects in the sheet assignment were
invisible until a β-rich reference existed.

See [`okf/validation/labelizer_ab.md`](../../../okf/validation/labelizer_ab.md).
