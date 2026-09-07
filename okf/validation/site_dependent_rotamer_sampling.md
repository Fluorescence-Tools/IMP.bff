---
type: validation
title: A dye library is a fixed sample, and a site decides how much of it survives
description: Measured effective sample size per labelling site across the shipped cutoff ladder, and what the shallow default costs the FRET observable. At a tight site the cutoff-30 library carries ESS 1.7 and lands 0.068 in E away from the deep library; at open sites the same comparison moves E by 0.008.
resource: /Users/tpeulen/dev/imp.bff
tags: [validation, imp.bff, rotamer, sampling, cgprobe, fret, kappa2, prd-116, prd-118]
timestamp: '2026-08-24T00:00:00Z'
---

# The thing that is easy to miss

A dye rotamer library is **one fixed sample of the free dye**, generated once by
MD with no protein in sight. At a labelling site it is *reweighted* — every
conformer gets a Boltzmann weight against the local steric and electrostatic
environment — and it is never *resampled*. So the quality of a site's answer is
not set by the library's size. It is set by the overlap between that fixed
sample and the conformations the site actually allows, and that overlap is a
property of the site, not of the library.

The count of surviving rotamers hides this: a site can keep twenty conformers
and put three quarters of the weight on one of them. What does not hide it is
Kish's effective sample size,

    ESS = (sum w)^2 / sum w^2

which says how many independent conformers the weighted ensemble is really
carrying. It is on every ensemble as
`RotamerEnsemble.effective_sample_size`.

# Measured, on the bundled parity fixtures

Two structures, both sites of each, across the three shipped cutoffs
(`test/cgprobe/rotamer/data`, the FRETpredict pins' own fixtures):

| system | site | cutoff | library | w > 0 | ESS | ESS/N | top-1 weight |
|---|---|---|---|---|---|---|---|
| Hsp90 | 452 | 10 | 575 | 421 | 186.2 | 0.32 | 0.03 |
| Hsp90 | 452 | 30 | 37 | 20 | 10.4 | 0.28 | 0.20 |
| Hsp90 | **637** | 10 | 446 | 293 | **71.7** | 0.13 | 0.10 |
| Hsp90 | **637** | 20 | 40 | 20 | **9.5** | 0.22 | 0.23 |
| Hsp90 | **637** | 30 | **7** | **2** | **1.7** | 0.24 | **0.72** |
| pp11 | 0 | 10 | 711 | 644 | 310.1 | 0.44 | 0.02 |
| pp11 | 0 | 30 | 33 | 32 | 18.4 | 0.56 | 0.10 |
| pp11 | 12 | 30 | 37 | 34 | 15.3 | 0.41 | 0.13 |

**Hsp90 residue 637 with the default library is an ensemble of 1.7 conformers.**
Seven rotamers exist at cutoff 30, two survive the site, and one of them holds
72 % of the weight. Every orientational quantity that site reports — its
`kappa2`, and through it `R0` and the efficiencies — is that one conformer's.

# What it costs the answer

Running the same case at each cutoff, changing nothing else:

| | cutoff30 | cutoff20 | cutoff10 | default vs deepest |
|---|---|---|---|---|
| Hsp90 `Es` | 0.4677 | 0.4036 | 0.3999 | **ΔEs 0.068** |
| Hsp90 `<kappa^2>` | 0.9693 | 0.6253 | 0.6459 | **Δ 0.32** |
| pp11 `Es` | 0.7464 | 0.7540 | 0.7386 | ΔEs 0.008 |
| pp11 `<kappa^2>` | 0.7639 | 0.7245 | 0.6844 | Δ 0.08 |

At the tight site the default is **0.068 in E** away from the deep library —
several times a good experiment's precision, and four thousand times the
tolerance the FRETpredict parity pins are held to. At the open sites the same
change is 0.008. That gap between the two rows *is* the site dependence: the
same library, the same code, and an error that varies by an order of magnitude
with where the dye is attached.

`<kappa^2>` shows the mechanism. At cutoff 30 the tight site reports 0.97,
far from the isotropic 2/3, not because the dye is genuinely ordered there but
because two conformers cannot represent a distribution of orientations.

# What to do about it

1. **Read the ESS.** It is free and it is the only number that says whether a
   site's answer is supported. Below roughly ten, do not trust the third
   decimal of anything that site reports.
2. **Deepen where it is low.** `cutoff10` ships for every library in
   `dyes.drot.pto` and is 10–100x larger; at Hsp90 637 it takes the ESS from
   1.7 to 71.7. Loading it is cheap since the container work of 2026-08-24
   (the whole 95-library corpus loads in 0.51 s), so depth is now a decision
   about accuracy rather than about time.
3. **The real fix is resampling, not reweighting.** A `.drot` conformer *is* a
   dihedral vector (PRD-118), so conformers can be generated near the survivors
   and rebuilt through the same `internal2cartesian` kernel — importance
   sampling conditioned on the site instead of filtering a free-dye sample.
   That changes the numbers by design and so is a decision, not a default; the
   store was built to make it possible and nothing here does it yet.

Note the shipped default is FRETpredict's default (cutoff 30), and the parity
pins are recorded against it. That is correct for parity and wrong for
accuracy at a tight site — the pins prove this package reproduces FRETpredict,
including where FRETpredict is under-sampled.
