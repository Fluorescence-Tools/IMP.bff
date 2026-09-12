/*
 * A FitDataset, and the one helper that lets an existing curve class inherit its
 * probability calculus without rewriting how it stores anything.
 */

%pythoncode %{
def sync_dataset(dataset, y, ey=None, x=None, mask=None, shape=None,
                 coordinate_name="x"):
    """Fill a `FitDataset` from the arrays a measured curve already holds.

    Written so a class that stores its data as numpy arrays -- chisurf's
    `DataCurve`, for instance -- can add `IMP.bff.FitDataset` to its bases and
    call this whenever those arrays change. It then *is* a FitDataset: it scores,
    it propagates, and it can be handed to `FitChiSquared.set_dataset`, without
    anything about its own storage changing.

    The arrays go through the ndarray setters and not the list ones, because
    a curve's values change on every write and the list path walks them one
    Python object at a time -- measured at 3.9 ms against 316 us for 117k
    points. A sync that costs more than the arithmetic gets skipped, and then
    the two halves drift.

    The noise family is chosen rather than demanded: given errors, the data
    carries a stored variance of `ey**2`; given none, it is counts. Say so
    afterwards with `set_noise_family` if that is wrong -- and see PRD-140 for
    why storing `sqrt(y)` in `ey` and calling it a variance is Neyman, and
    biased.
    """
    import numpy as _np
    y = _np.ascontiguousarray(_np.asarray(y, dtype=float).ravel())
    if shape is not None:
        dataset.set_values(list(y), list(int(s) for s in shape))
    else:
        dataset.set_values_array(y)
    if ey is not None:
        ey = _np.ascontiguousarray(_np.asarray(ey, dtype=float).ravel())
        if ey.size == y.size and _np.any(ey != 0.0):
            dataset.set_noise_family(FIT_NOISE_FAMILY_STORED)
            dataset.set_stored_variance_array(
                _np.ascontiguousarray(ey.astype(float) ** 2))
        else:
            dataset.set_noise_family(FIT_NOISE_FAMILY_POISSON)
    else:
        dataset.set_noise_family(FIT_NOISE_FAMILY_POISSON)
    if x is not None:
        x = _np.ascontiguousarray(_np.asarray(x, dtype=float).ravel())
        if x.size == y.size:
            dataset.set_coordinate_array(0, coordinate_name, x)
    if mask is not None:
        mask = _np.ascontiguousarray(_np.asarray(mask, dtype=float).ravel())
        if mask.size == y.size:
            dataset.set_mask_array(mask)
    return dataset
%}
