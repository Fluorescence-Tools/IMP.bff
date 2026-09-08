---
type: validation
title: "Can chisurf use bff's Dataset, and why it does not already"
description: Asked whether chisurf can use bff for the magic-angle error propagation, and why bff curves are not used throughout. Half the answer is that it already does — chisurf's model side calls bff nodes today, including AnisotropySpectrum and ChiSquared. The data side does not, and the reasons are real: DataCurve carries an axis, file provenance and a GUI object model across 61 files, and until today bff's Dataset had no axis at all. The strongest argument for moving is one chisurf already concedes in its own code: it flattens 2-D data to 1-D and carries the shape in a metadata dictionary, used by two readers and the fit controller. That is a shape, badly.
resource: /Users/tpeulen/dev/imp.bff
tags: [validation, imp.bff, chisurf, dataset, prd-140]
timestamp: '2026-09-08T00:00:00Z'
---
# Can chisurf use it, and why not all along

## It already does, on the model side

`chisurf/core/fluorescence/mfd/patterns.py` builds an anisotropy spectrum by
constructing `bff.AnisotropySpectrum`, wiring its ports and evaluating it.
`chisurf/core/fitting/minimizer.py` sets `ChiSquared`'s noise model per fit.
So the question is not whether chisurf can use bff nodes — it does, in the
place where the model is computed.

What has never moved is the **data**.

## What chisurf's DataCurve is, and bff's Dataset is not

`chisurf/core/data.py`'s `DataCurve` is `Curve` plus `ExperimentalData`:

| | DataCurve | Dataset |
|---|---|---|
| values | `y` | any rank |
| axis | one `x` array, and `ex` for its uncertainty | **added 2026-09-08**: any number of coordinates, one value per point; `ex` still absent |
| uncertainty | a stored `ey` array | a noise *family*, and propagation over sources |
| mask | yes | yes |
| N-D | flattened to 1-D plus `meta_data['grid']` = `{ndim, shape, order}` | a real shape |
| provenance | filename, metadata, readers, save/load | none |
| object model | `DataGroup`, `DataCurveGroup`, GUI bindings, write-locked arrays | none |

Two of those rows are the whole answer.

**Why not already: the axis, the provenance and the object model.** A dataset
with no `x` cannot be a curve, which until today was disqualifying on its own.
Note what the coordinate had to become, though, on the owner's correction:
**every value carries its own coordinate**, because an axis can be any shape
and samples need not lie on a lattice — and **how many coordinates a dataset
has is independent of its rank**, because a list of bursts is rank 1 and
carries an efficiency and a stoichiometry. `DataCurve`'s single `x` is the
rank-1 separable case of that, which is the direction a migration would have
to go rather than the reverse.
The rest — where the file came from, how it is grouped, how a widget binds to
it — is chisurf's business and should stay there. `DataCurve` is named in 61
files.

**Why it should move anyway: the N-D row.** chisurf already carries 2-D data
as a flattened 1-D curve with a dictionary saying `ndim`, `shape` and `order`,
and it is not decoration — `experiments/pda2c/reader.py` and
`experiments/mfd/reader.py` write it and `gui/widgets/fitting/fit_controller.py`
reads it to build selection masks. That is a shape, reimplemented in a
metadata dict because the curve class could not hold one. bff's Dataset holds
one properly, and scores a rank-3 dataset exactly as it scores the flattened
equivalent.

## And the thing that prompted the question

The magic angle. `VV + 2G·VH` is not Poisson — its variance is
`VV + 4G²·VH`, the coefficient squares — and anisotropy is worse, because
`(VV − G·VH)/(VV + 2G·VH)` shares both channels above and below the line, so
numerator and denominator are correlated. bff's `Dataset` now propagates over
the independent sources and gets both right; `DataCurve` carries a stored `ey`
and has nowhere to record what the number was built from.

So a chisurf curve can be *constructed* through bff today: build the channels
as sources, compose the observable, take `sqrt(variance)` into `ey`. That is a
small, useful integration and it needs no migration at all.

## The shape of a migration, if one is wanted

1. **Now, no migration**: chisurf constructs combined observables through
   `Dataset` and keeps the result in `DataCurve`'s `ey`. It fixes the
   propagation and touches one call site per observable.
2. **Then**: `DataCurve` holds a `Dataset` rather than raw arrays, keeping its
   own `x`, provenance and grouping. The 61 files see no change.
3. **Only if it earns it**: the `meta_data['grid']` path retires in favour of
   the dataset's shape, which is the one place the current design is doing
   something the abstraction does properly.

What should not move: readers, file provenance, groups, GUI binding. bff
processes data; it should not learn where a file came from.

## Still missing before step 2

- `ex`, a coordinate's own uncertainty. `DataCurve` has it for its single `x`;
  `Dataset` has coordinates but no uncertainty on them. It should probably be
  the same propagation machinery the values use rather than a second stored
  array, which is a design question and not a gap to fill blindly.
- No consumer yet uses `Dataset`'s rank > 1 in anger — the rank-3 case is
  tested, not used, and a migration is the wrong time to discover what an
  image needs.
