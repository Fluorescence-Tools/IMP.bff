---
type: reference
title: "The four stages: how to think about IMP.bff, independent of where the files are"
description: A dye model in this package is a representation, a scoring of it, a way of sampling it, and an analysis that projects it onto an experiment. Each stage has a variant per representation -- accessible volume, rotamer library, coarse-grained dye. This is a conceptual map, deliberately not a directory layout.
resource: /Users/tpeulen/dev/imp.bff
tags: [reference, imp.bff, architecture, conceptual]
timestamp: '2026-08-18T00:00:00Z'
---

# The four stages

**This is a way of thinking, not a directory layout.** The repository does not
mirror it and does not need to — the point is that any question about the dye
model can be placed in one of four boxes, and that each box has a variant per
representation. Reorganising the tree to match would be churn for no gain.

|  | accessible volume | rotamer library | coarse-grained dye |
|---|---|---|---|
| **1. representation** | the allowed region, uniform on it | listed conformers with weights | explicit beads under a simplified force field |
| **2. scoring** | occupancy: reachable or not | clash energy per conformer | the force field's internal energy |
| **3. sampling** | pathfinding on the grid | library screening | MD / Langevin |
| **4. analysis** | dye density, and the projection to an experiment | | |

## What each stage answers

**Representation — where can the dye be?** Three answers, and the choice is a
modelling decision, not an implementation detail. An accessible volume says
"anywhere in this region, equally"; a rotamer library says "one of these
listed states, with these weights"; a coarse-grained dye says "wherever the
force field allows, weighted by Boltzmann".

**Scoring — is this configuration allowed, and how heavily does it count?**
This is what turns candidates into an ensemble. Note it is *not*
`IMP.bff.restraints`, which asks a different question — does the model agree
with a *measurement*. Scoring never sees an experiment.

**Sampling — how do we enumerate or draw configurations?** Pathfinding for a
volume, screening for a library, dynamics for an explicit dye. The Brownian
walk and the Smoluchowski solver are the same dynamics sampled or integrated,
so they belong to the same box.

**Analysis — what do we report?** The dye density, and the projection onto
what an experiment would see. The projection has its own contract
(`IMP.bff.observables`): experiment-*neutral* quantities only — lifetime
spectra, rate constants — with convolution, pileup, binning and counting noise
left outside the package.

## Why it is not the directory layout

The stages cut across the representations, so a directory tree has to pick one
axis and betray the other. The tree is cut by **domain** — `dye`, `label`,
`representation`, `dynamics`, `photophysics`, `observables`, `io`, `scoring`,
`restraints` — because that is what an import has to name. The table above is
what a reader should hold in their head while reading it.

Where the two happen to agree, they agree: `representation/` really is stage 1,
`scoring/` really is stage 2, `observables/` really is the projection half of
stage 4. Where they do not, the table is the authority on *meaning* and the
tree is the authority on *imports*.
