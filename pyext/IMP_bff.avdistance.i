/*
 * The array kernels an accessible volume is measured with.
 *
 * `representation/av.py` carried a `_kernels` section that was six functions
 * doing nothing but `.ravel()` on the way in and `.reshape()` on the way out of
 * these. The shape belongs to the caller's picture of the data, so it is stated
 * once, here, on the kernel itself -- rather than in a module whose only
 * content was stating it a second time.
 */

%include "IMP/bff/AVDistance.h"

%pythoncode %{
def random_distances(p1, p2, n_samples, seed=0):
    """``(n_samples, 2)`` of distance and weight product, from two clouds.

    Each sample picks one point from each cloud independently, so the pairs
    sample the joint distribution.
    """
    out = _IMP_bff.random_distances(
        np.ascontiguousarray(np.asarray(p1, dtype=np.float64)).ravel(),
        np.ascontiguousarray(np.asarray(p2, dtype=np.float64)).ravel(),
        int(n_samples), int(seed))
    return np.asarray(out, dtype=np.float64).reshape(int(n_samples), 2)


def density_to_points(density, dg, r0, threshold=0.0):
    """A density grid as an ``(n, 4)`` weighted point cloud."""
    d = np.ascontiguousarray(np.asarray(density, dtype=np.float64))
    nx, ny, nz = d.shape
    out = _IMP_bff.density_to_points(
        d.ravel(), int(nx), int(ny), int(nz), float(dg),
        np.asarray(r0, dtype=np.float64).ravel(), float(threshold))
    return np.asarray(out, dtype=np.float64).reshape(-1, 4)


def split_contact_volume_masks(density, dg, radius, rs, r0):
    """``(contact, free)`` uint8 masks over the density grid.

    ``uint8``, not ``float64``: these are binary masks, and at ``ng = 92`` a
    float64 pair costs 12.5 MB per labelling site against 1.6 MB. A residue
    scan builds one per site.
    """
    d = np.ascontiguousarray(np.asarray(density, dtype=np.float64))
    ng = int(d.shape[0])
    centres = np.ascontiguousarray(np.asarray(rs, dtype=np.float64)).reshape(-1, 3)
    r = np.asarray(radius, dtype=np.float64).ravel()
    if r.size != centres.shape[0]:
        r = np.full(centres.shape[0], float(r[0]))
    label = np.asarray(_IMP_bff.split_contact_volume(
        d.ravel(), ng, float(dg), r, centres.ravel(),
        np.asarray(r0, dtype=np.float64).ravel()), dtype=np.int32
    ).reshape(ng, ng, ng)
    return ((label == AV_VOXEL_CONTACT).astype(np.uint8),
            (label == AV_VOXEL_FREE).astype(np.uint8))
%}
