/*
 * Distances between two labels, whatever represents them.
 *
 * `representation/distance.py` asked all of these of a `States`, and every one
 * of them was a `.ravel()`, a call, and a `.reshape()`. The pair geometry and
 * its efficiencies came back as dictionaries of matrices, which is what a C++
 * value is for.
 *
 * The dict surface survives on the two results, because that is what a caller
 * indexes -- `geometry["R"]`, `out["static"]` -- and there is no reason to make
 * every one of them change.
 */

IMP_SWIG_VALUE(IMP::bff, FRETPairGeometry, FRETPairGeometries);
IMP_SWIG_VALUE(IMP::bff, FRETDistanceConverter, FRETDistanceConverters);

// The obstacle array, straight from numpy.
%apply(double* IN_ARRAY2, int DIM1, int DIM2) {
    (double* atoms_xyzr, int n_atoms, int n_cols)
};

%attribute_py(IMP::bff::LabelDistribution, std::string, position_name,
              get_position_name);
%attribute_py(IMP::bff::LabelDistribution, std::string, simulation_type,
              get_simulation_type);
%attribute_py(IMP::bff::LabelDistribution, double, simulation_grid_resolution,
              get_simulation_grid_resolution);
%attribute_py(IMP::bff::LabelDistribution, int, n_points, get_n_points);
%attribute_py(IMP::bff::DyeDistributionNormal, double, width, get_width);
%attribute_py(IMP::bff::FRETDistanceConverter, double, forster_radius,
              get_forster_radius, set_forster_radius);
%attribute_py(IMP::bff::FRETDistanceConverter, double, sigma, get_sigma,
              set_sigma);
IMP_SWIG_VALUE(IMP::bff, FRETPairEfficiencies, FRETPairEfficienciesList);

%rename(_fret_pair_geometry) IMP::bff::fret_pair_geometry;
%rename(_fret_pair_efficiencies) IMP::bff::fret_pair_efficiencies;
%rename(_histogram_rda) IMP::bff::histogram_rda;
%rename(_fit_transfer_polynomial) IMP::bff::fit_transfer_polynomial;
%rename(_av_pair_statistics) IMP::bff::av_pair_statistics;

// The distributions and the polynomial evaluators live here too: their renames
// have to precede the headers, and the headers are included below.
%include "IMP_bff.distributions.i"

%feature("shadow") IMP::bff::LabelDistributionAV::LabelDistributionAV %{
def __init__(self, atoms_xyz, atoms_vdw, linker_length=20.0, linker_width=0.5,
             dye_radii=(3.5, 0.0, 0.0), source_xyz=None,
             simulation_grid_resolution=1.5, position_name=""):
    xyz = np.ascontiguousarray(np.asarray(atoms_xyz, dtype=np.float64)).reshape(-1, 3)
    vdw = np.asarray(atoms_vdw, dtype=np.float64).ravel()
    r = (list(dye_radii) + [0.0, 0.0, 0.0])[:3]
    # No `source_xyz` means the first obstacle. The Python this replaces took a
    # residue number and an atom name and then ignored both, returning index 0
    # from a loop whose body was a comment saying so.
    src = xyz[0] if source_xyz is None else np.asarray(source_xyz, dtype=np.float64)
    _IMP_bff.LabelDistributionAV_swiginit(self, _IMP_bff.new_LabelDistributionAV(
        np.ascontiguousarray(np.column_stack([xyz, vdw])),
        np.asarray(src, dtype=np.float64).ravel(),
        float(linker_length), float(linker_width),
        float(r[0]), float(r[1]), float(r[2]),
        float(simulation_grid_resolution), str(position_name)))
%}

%feature("shadow") IMP::bff::DyeDistributionNormal::DyeDistributionNormal %{
def __init__(self, origin, width=6.0, n_points=50000, seed=0,
             position_name=""):
    _IMP_bff.DyeDistributionNormal_swiginit(
        self, _IMP_bff.new_DyeDistributionNormal(
            np.asarray(origin, dtype=np.float64).ravel(), float(width),
            int(n_points), int(seed), str(position_name)))
%}

// The renames above have to precede both headers, so `FRETPair.h` is included
// here rather than beside its kernels: the value types it now carries are what
// this file gives a Python surface.
%include "IMP/bff/DistanceCalibration.h"
%include "IMP/bff/FRETPair.h"
%include "IMP/bff/StatesDistance.h"

%extend IMP::bff::LabelDistribution {
    %pythoncode %{
        @property
        def points(self):
            """The cloud as ``(n, 4)``, computing it on first use."""
            return self.get_accessible_volume().points

        @property
        def mean_position(self):
            return self.get_accessible_volume().mean_position

        @property
        def origin(self):
            o = _IMP_bff.LabelDistribution_get_origin(self)
            return None if o.size == 0 else o

        @property
        def states(self):
            """The representation-agnostic view: whatever produced the cloud,
            a consumer that wants positions and weights asks for this."""
            return self.get_accessible_volume()

        def dRmp(self, other):
            return self.get_accessible_volume().dRmp(
                other.get_accessible_volume())

        def dRDA(self, other, n_samples=50000):
            return self.get_accessible_volume().dRDA(
                other.get_accessible_volume(), n_samples)

        def dRDAE(self, other, forster_radius=52.0, n_samples=50000):
            return self.get_accessible_volume().dRDAE(
                other.get_accessible_volume(), forster_radius, n_samples)

        def pRDA(self, other, axis=None, n_samples=50000):
            return self.get_accessible_volume().pRDA(
                other.get_accessible_volume(), axis, n_samples)
    %}
}

%extend IMP::bff::FRETDistanceConverter {
    %pythoncode %{
        def __call__(self, value, distance_type):
            return self.get_effective_distance(float(value), int(distance_type))
    %}
}

%extend IMP::bff::FRETPairGeometry {
    %pythoncode %{
        @property
        def R(self):
            """``(n1, n2)`` separations, A."""
            return _IMP_bff.FRETPairGeometry_get_R(self).reshape(self.n1, self.n2)

        @property
        def kappa2(self):
            """``(n1, n2)`` orientation factors."""
            return _IMP_bff.FRETPairGeometry_get_kappa2(self).reshape(
                self.n1, self.n2)

        @property
        def weight(self):
            """``(n1, n2)`` pair weights, normalised to sum 1."""
            return _IMP_bff.FRETPairGeometry_get_weight(self).reshape(
                self.n1, self.n2)

        def __getitem__(self, key):
            """The dictionary surface the Python result had."""
            return {"R": self.R, "kappa2": self.kappa2, "weight": self.weight,
                    "kappa2_avg": self.kappa2_avg}[key]
    %}
}

%extend IMP::bff::FRETPairEfficiencies {
    %pythoncode %{
        @property
        def E(self):
            """``(n1, n2)`` per-pair efficiencies."""
            return _IMP_bff.FRETPairEfficiencies_get_E(self).reshape(
                self.n1, self.n2)

        @property
        def rate_ratio(self):
            """``(n1, n2)`` of k_FRET / k_rad. **Not sanitised** -- an infinite
            rate at zero separation is true, and the caller decides."""
            return _IMP_bff.FRETPairEfficiencies_get_rate_ratio(self).reshape(
                self.n1, self.n2)

        @property
        def k_fret(self):
            """``(n1, n2)`` FRET rates, 1/ns, or ``None`` without a lifetime."""
            k = _IMP_bff.FRETPairEfficiencies_get_k_fret(self)
            return None if k.size == 0 else k.reshape(self.n1, self.n2)

        @property
        def R(self):
            """``(n1, n2)`` separations, A -- from the geometry."""
            return _IMP_bff.FRETPairEfficiencies_get_R(self).reshape(
                self.n1, self.n2)

        @property
        def kappa2(self):
            """``(n1, n2)`` orientation factors -- from the geometry."""
            return _IMP_bff.FRETPairEfficiencies_get_kappa2(self).reshape(
                self.n1, self.n2)

        @property
        def weight(self):
            """``(n1, n2)`` pair weights -- from the geometry."""
            return _IMP_bff.FRETPairEfficiencies_get_weight(self).reshape(
                self.n1, self.n2)

        @property
        def forster_radius_nm(self):
            """R0 in nm, which is the unit the spectra derive it in."""
            return self.forster_radius / 10.0

        def __getitem__(self, key):
            """The dictionary surface the Python result had."""
            out = {"static": self.static_efficiency, "dynamic1": self.dynamic1,
                   "dynamic2": self.dynamic2, "kappa2_avg": self.kappa2_avg,
                   "E": self.E, "rate_ratio": self.rate_ratio,
                   "R": self.R, "kappa2": self.kappa2, "weight": self.weight,
                   "forster_radius": self.forster_radius,
                   "forster_radius_nm": self.forster_radius_nm}
            if self.k_fret is not None:
                out["k_fret"] = self.k_fret
            return out[key]

        def __contains__(self, key):
            try:
                self[key]
            except KeyError:
                return False
            return True
    %}
}

%pythoncode %{
def _xyz(points):
    """The first three columns of a cloud, flat. Accepts `(n, 3)` or `(n, 4)`."""
    p = np.ascontiguousarray(np.asarray(points, dtype=np.float64))
    if p.ndim == 2 and p.shape[1] >= 3:
        p = np.ascontiguousarray(p[:, :3])
    return p.ravel()


def fret_pair_geometry(points1, weights1, points2, weights2, mu1=None, mu2=None):
    """Distances, kappa^2 and pair weights over all (i, j) of two state sets.

    Without dipoles kappa^2 is the isotropic 2/3 everywhere, which is what an
    accessible volume can say: a point cloud carries no orientation.
    """
    empty = np.zeros(0)
    return _IMP_bff._fret_pair_geometry(
        _xyz(points1), np.asarray(weights1, dtype=np.float64).ravel(),
        _xyz(points2), np.asarray(weights2, dtype=np.float64).ravel(),
        empty if mu1 is None else _xyz(mu1),
        empty if mu2 is None else _xyz(mu2))


def fret_pair_efficiencies(geometry, forster_radius, tau0=None):
    """The efficiencies of a pair geometry, in the three averaging limits."""
    return _IMP_bff._fret_pair_efficiencies(
        geometry, float(forster_radius),
        -1.0 if tau0 is None else float(tau0))


def fret_pair_distribution(points1, weights1, points2, weights2, *,
                           forster_radius, mu1=None, mu2=None, tau0=None):
    """:func:`fret_pair_geometry` followed by :func:`fret_pair_efficiencies`."""
    geometry = fret_pair_geometry(points1, weights1, points2, weights2, mu1, mu2)
    return fret_pair_efficiencies(geometry, forster_radius, tau0)


def histogram_rda(s1, s2, rda_axis=None, rda_min=1.0, rda_max=200.0,
                  n_rda_bins=100, use_log=False, n_samples=50000,
                  normalize=False):
    """``(histogram, bin_edges)`` of the pair-distance distribution."""
    if rda_axis is None:
        if use_log:
            rda_axis = np.logspace(np.log10(rda_min), np.log10(rda_max),
                                   n_rda_bins, dtype=np.float64)
        else:
            rda_axis = np.linspace(rda_min, rda_max, n_rda_bins,
                                   dtype=np.float64)
    axis = np.ascontiguousarray(np.asarray(rda_axis, dtype=np.float64))
    p = _IMP_bff._histogram_rda(s1, s2, axis, int(n_samples), bool(normalize))
    return p, axis


def fit_transfer_polynomial(s1, s2, distance_type, forster_radius=52.0,
                            degree=3, n_samples=10000):
    """Fit ``Rmp -> RDAMean`` (or ``RDAMeanE``) by translating the two clouds.

    Coefficients highest power first, as ``numpy.polyfit`` returns them.
    """
    return np.asarray(_IMP_bff._fit_transfer_polynomial(
        s1, s2, str(distance_type), float(forster_radius), int(degree),
        int(n_samples)), dtype=np.float64)


def av_pair_statistics(s1, s2, forster_radius=52.0, n_samples=50000):
    """``(Rmp, <R_DA>, <R_DA>_E, sigma_R)`` in one pass."""
    return tuple(_IMP_bff._av_pair_statistics(
        s1, s2, float(forster_radius), int(n_samples)))
%}
