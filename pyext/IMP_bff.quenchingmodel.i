/*
 * A labelled site's donor decay: the field picture and the particle picture,
 * both as C++ objects.
 *
 * These were the two orchestrating classes of `quenching.py`, and they were
 * orchestrating across the boundary: every field is `ng^3` doubles, and the
 * Python held all of them, handed them down to a kernel and took them back on
 * every call. At `ng = 41` one field is 69 000 doubles and a donor decay
 * touches four of them per step.
 *
 * The shims below keep two Python conveniences that have no C++ spelling and
 * that every caller uses. First, `atoms` may be a numpy **structured array**
 * -- `chain`, `res_id`, `res_name`, `atom_name`, `coord` -- which is how a
 * structure reader hands its atoms over; a structured dtype is a numpy idea, so
 * the C++ takes the fields as parallel arrays and the shim unpacks. Second,
 * `av` may be anything carrying `density`, `grid_step` and `attachment_point`,
 * because that is all either model reads of a volume and a test has no reason
 * to build a whole one.
 */

IMP_SWIG_VALUE(IMP::bff, ObstacleAtoms, ObstacleAtomsList);
IMP_SWIG_VALUE(IMP::bff, DynamicAccessibleVolume, DynamicAccessibleVolumes);
IMP_SWIG_VALUE(IMP::bff, QuenchedDonorDecay, QuenchedDonorDecays);

// `run()` on the solver is not this; `donor_decay` returns one buffer whose
// first value is the resolved time step, and the shim splits it.
%rename(_donor_decay) IMP::bff::DynamicAccessibleVolume::donor_decay;
%rename(_decay_histogram) IMP::bff::QuenchedDonorDecay::decay_histogram;
%rename(_fret_rate_trace_paired) IMP::bff::QuenchedDonorDecay::fret_rate_trace_paired;
%rename(_fret_rate_trace_cloud) IMP::bff::QuenchedDonorDecay::fret_rate_trace_cloud;
%rename(_simulate_photons) IMP::bff::QuenchedDonorDecay::simulate_photons;
%rename(_photons_fused) IMP::bff::QuenchedDonorDecay::photons_fused;
%rename(_lifetime_spectrum) IMP::bff::QuenchedDonorDecay::lifetime_spectrum;

%feature("shadow") IMP::bff::ObstacleAtoms::ObstacleAtoms %{
def __init__(self, chains=None, res_ids=None, res_names=None, atom_names=None,
             coords=None):
    """The obstacles a quenching model needs of a structure.

    Every field optional: the field picture reads only `res_name`, `atom_name`
    and `coord`, and only the particle picture groups atoms into residues.
    """
    _IMP_bff.ObstacleAtoms_swiginit(self, _IMP_bff.new_ObstacleAtoms())
    self.chains = [] if chains is None else [str(c) for c in chains]
    self.res_names = [] if res_names is None else [str(c) for c in res_names]
    self.atom_names = [] if atom_names is None else [str(c) for c in atom_names]
    self.res_ids = [] if res_ids is None else [int(i) for i in res_ids]
    self.coords = ([] if coords is None else
                   [float(v) for v in np.asarray(coords, dtype=np.float64).ravel()])
%}

%feature("shadow") IMP::bff::DynamicAccessibleVolume::DynamicAccessibleVolume %{
def __init__(self, av=None, atoms=None, tau0=4.0, dye_radius=3.5,
             free_diffusion=8.0, contact_distance=6.5, slow_factor=0.985,
             flux_form="smoluchowski"):
    """An accessible volume with diffusion, quenching and FRET rate fields.

    :param av: the :class:`IMP.bff.AccessibleVolume` to decorate, or anything
        carrying ``density``, ``grid_step`` and ``attachment_point``.
    :param atoms: the obstacles, as :class:`IMP.bff.ObstacleAtoms` or a numpy
        structured array with ``res_name``, ``atom_name`` and ``coord``.
    :param tau0: unquenched donor lifetime in ns.
    :param dye_radius: dye sphere radius in Angstrom.
    :param free_diffusion: unhindered diffusion coefficient in A^2/ns.
    :param contact_distance: dye-to-atom distance counted as contact, for the
        mobility field.
    :param slow_factor: mobility multiplier applied per contacting atom.
    :param flux_form: where the dye sits at equilibrium -- see
        :func:`IMP.bff.equilibrium_occupancy`.
    """
    _IMP_bff.DynamicAccessibleVolume_swiginit(
        self, _IMP_bff.new_DynamicAccessibleVolume(
            _as_volume(av), _as_obstacles(atoms), float(tau0),
            float(dye_radius), float(free_diffusion), float(contact_distance),
            float(slow_factor), str(flux_form)))
%}

%feature("shadow") IMP::bff::QuenchedDonorDecay::QuenchedDonorDecay %{
def __init__(self, av=None, atoms=None, tau0=4.0, quenching_table=None,
             critical_distance=7.0, slow_radius=10.0, dye_radius=None,
             diffusion_coefficient=40.0, slow_fact=0.05, t_step=0.004,
             t_max=10000.0, n_photons=100000, n_trajectories=-1,
             random_seed=None):
    """A labelled site's donor decay, from the structure and the PET chemistry.

    :param av: the donor's :class:`IMP.bff.AccessibleVolume`, or anything
        carrying ``density``, ``grid_step`` and ``attachment_point``.
    :param atoms: the obstacles, as :class:`IMP.bff.ObstacleAtoms` or a numpy
        structured array with ``chain``, ``res_id``, ``res_name``,
        ``atom_name`` and ``coord``.
    :param tau0: unquenched donor lifetime in ns.
    :param quenching_table: per-residue interaction table; the defaults of
        :func:`IMP.bff.amino_acid_quenching_defaults` when omitted.
    :param critical_distance: contact radius inherited by residues whose table
        entry leaves ``quench_radius`` unset.
    :param slow_radius: radius of each residue's sticky sphere.
    :param n_photons: excitation events for the photon Monte-Carlo.
    """
    _IMP_bff.QuenchedDonorDecay_swiginit(
        self, _IMP_bff.new_QuenchedDonorDecay(
            _as_volume(av), _as_obstacles(atoms), float(tau0),
            _pet_table(quenching_table), float(critical_distance),
            float(slow_radius),
            DEFAULT_DYE_RADIUS if dye_radius is None else float(dye_radius),
            float(diffusion_coefficient), float(slow_fact), float(t_step),
            float(t_max), int(n_photons), int(n_trajectories),
            -1 if random_seed is None else int(random_seed)))
%}

%attribute_py(IMP::bff::DynamicAccessibleVolume, double, tau0, get_tau0);
%attribute_py(IMP::bff::DynamicAccessibleVolume, double, dg, get_dg);
%attribute_py(IMP::bff::DynamicAccessibleVolume, std::string, flux_form,
              get_flux_form);
%attribute_py(IMP::bff::QuenchedDonorDecay, double, tau0, get_tau0);
%attribute_py(IMP::bff::QuenchedDonorDecay, double, dg, get_dg);
%attribute_py(IMP::bff::QuenchedDonorDecay, int, n_photons, get_n_photons);

%include "IMP/bff/QuenchingModel.h"

%pythoncode %{
_OBSTACLE_FIELDS = ("chain", "res_id", "res_name", "atom_name", "coord")


def _as_obstacles(atoms):
    """`ObstacleAtoms` from one, or from a numpy structured array."""
    if atoms is None:
        return _IMP_bff.ObstacleAtoms()
    if isinstance(atoms, ObstacleAtoms):
        return atoms
    names = getattr(getattr(atoms, "dtype", None), "names", None) or ()
    def _field(key):
        return atoms[key] if key in names else None
    return ObstacleAtoms(
        chains=_field("chain"), res_ids=_field("res_id"),
        res_names=_field("res_name"), atom_names=_field("atom_name"),
        coords=_field("coord"))


def _as_volume(av):
    """An `AccessibleVolume`, or one built from whatever carries its three parts.

    A model reads exactly `density`, `grid_step` and `attachment_point` of a
    volume, so anything carrying those three is one as far as it is concerned.
    """
    if av is None:
        return AccessibleVolume()
    if isinstance(av, AccessibleVolume):
        return av
    return AccessibleVolume(
        density=getattr(av, "density", None),
        grid_step=float(getattr(av, "grid_step", 1.0)),
        attachment_point=getattr(av, "attachment_point", None))
%}

%extend IMP::bff::ObstacleAtoms {
    %pythoncode %{
        def __len__(self):
            return int(self.size())

        @property
        def coord(self):
            """``(n, 3)`` atom coordinates."""
            return _IMP_bff.ObstacleAtoms_get_coords(self).reshape(-1, 3)
    %}
}

%extend IMP::bff::DynamicAccessibleVolume {
    %pythoncode %{
        @property
        def density(self):
            """Binary occupancy of the accessible volume."""
            return _cube(_IMP_bff.DynamicAccessibleVolume_get_density(self))

        @property
        def x0(self):
            """The grid anchor -- the attachment point."""
            return _IMP_bff.DynamicAccessibleVolume_get_x0(self)

        @property
        def bounds(self):
            """Where the dye may be: the domain mask for the solver."""
            return _cube(_IMP_bff.DynamicAccessibleVolume_get_bounds(self))

        @property
        def diffusion_map(self):
            return _cube(_IMP_bff.DynamicAccessibleVolume_get_diffusion_map(self))

        @property
        def quenching_rate_map(self):
            return _cube(
                _IMP_bff.DynamicAccessibleVolume_get_quenching_rate_map(self))

        @property
        def fret_rate_map(self):
            """The FRET field, or ``None`` until an acceptor is set."""
            flat = _IMP_bff.DynamicAccessibleVolume_get_fret_rate_map(self)
            return None if flat.size == 0 else _cube(flat)

        @property
        def rate_map(self):
            """Total decay rate per voxel: quenching, plus FRET if set."""
            return _cube(_IMP_bff.DynamicAccessibleVolume_get_rate_map(self))

        @property
        def occupancy(self):
            """The equilibrium distribution of the dye, normalised.

            Under the default ``flux_form="smoluchowski"`` this **is** the AV
            density, uniform over the accessible voxels, and does not depend on
            the mobility -- equilibrium is thermodynamics, mobility is kinetics,
            and a dye slowed by friction with no attraction is still found
            everywhere it can reach.
            """
            return _cube(_IMP_bff.DynamicAccessibleVolume_get_occupancy(self))

        def update_diffusion_map(self, radial_profile=None):
            """Build the mobility field: free diffusion, slowed by nearby atoms."""
            base = (np.zeros(0) if radial_profile is None else
                    radial_diffusion_map(self.density, self.dg, radial_profile))
            _IMP_bff.DynamicAccessibleVolume_update_diffusion_map(
                self, np.ascontiguousarray(base, dtype=np.float64).ravel())
            return self.diffusion_map

        def update_quenching_map(self, quencher, rC=None):
            """Build the PET field from the dye's pair parameters.

            :param quencher: ``{comp_id: PETParameters}`` for one dye.
            :param rC: overrides every per-atom characteristic distance with one
                electron-transfer length, which is how ChiSurf drove it.
            """
            _IMP_bff.DynamicAccessibleVolume_update_quenching_map(
                self, ({} if not quencher else dict(quencher)),
                float("nan") if rC is None else float(rC))
            return self.quenching_rate_map

        def update_fret_map(self, acceptor, forster_radius=52.0,
                            acceptor_step=2):
            """Build the FRET field against an acceptor's accessible volume."""
            _IMP_bff.DynamicAccessibleVolume_update_fret_map(
                self, acceptor, float(forster_radius), int(acceptor_step))
            return self.fret_rate_map

        def update_occupancy(self, t_step=None, iterate=False, **kwargs):
            """The equilibrium occupancy, in closed form.

            Propagating to it instead is possible but slow and, in the ``"ito"``
            form on a real site where the compounding slow factor makes ``D``
            span orders of magnitude, may not converge at all. Pass
            ``iterate=True`` to do it the long way anyway.
            """
            if iterate:
                _IMP_bff.DynamicAccessibleVolume_update_occupancy_by_iteration(
                    self, -1.0 if t_step is None else float(t_step),
                    int(kwargs.pop("n_steps", 20000)),
                    float(kwargs.pop("tolerance", 1e-8)),
                    int(kwargs.pop("n_check", 100)))
            else:
                _IMP_bff.DynamicAccessibleVolume_update_occupancy(self)
            return self.occupancy

        def donor_decay(self, t_max=50.0, t_step=None, n_out=10):
            """Integrate the donor decay on the grid, from the equilibrium start.

            :returns: a :class:`IMP.bff.GridDiffusionResult` -- time axis in ns,
                surviving excited-state fraction, final density.
            """
            n_out = max(1, int(n_out))
            flat = _IMP_bff.DynamicAccessibleVolume__donor_decay(
                self, float(t_max), -1.0 if t_step is None else float(t_step),
                n_out)
            step = float(flat[0])
            ng = self.get_ng()
            n_reports = flat.size - 1 - ng ** 3
            fluorescence = flat[1:1 + n_reports]
            final = flat[1 + n_reports:].reshape(ng, ng, ng)
            time = np.arange(n_reports, dtype=np.float64) * step * n_out
            return GridDiffusionResult(time, fluorescence, final)
    %}
}

%extend IMP::bff::QuenchedDonorDecay {
    %pythoncode %{
        @property
        def density(self):
            return _cube(_IMP_bff.QuenchedDonorDecay_get_density(self))

        @property
        def x0(self):
            return _IMP_bff.QuenchedDonorDecay_get_x0(self)

        @property
        def sites(self):
            """The slow and quench centres of every residue in the structure.

            Memoised on the proxy as well as in C++: `get_sites` returns a
            reference, but SWIG copies a value type on the way out, so reading
            this in a loop would copy the whole table each time.
            """
            cached = self.__dict__.get("_sites")
            if cached is None:
                cached = _IMP_bff.QuenchedDonorDecay_get_sites(self)
                self.__dict__["_sites"] = cached
            return cached

        @property
        def has_walk(self):
            """Whether a walk has been run, by either path."""
            return _IMP_bff.QuenchedDonorDecay_get_has_walk(self)

        @property
        def n_frames(self):
            """Frames the walk produced -- reported by the fused path too."""
            return _IMP_bff.QuenchedDonorDecay_get_n_frames(self)

        @property
        def mean_k_quench(self):
            """Mean quenching rate along the walk, 1/ns."""
            return _IMP_bff.QuenchedDonorDecay_get_mean_k_quench(self)

        @property
        def collision_fraction(self):
            """Fraction of frames inside a quenching sphere."""
            return _IMP_bff.QuenchedDonorDecay_get_collision_fraction(self)

        @property
        def quenching_table(self):
            return dict(_IMP_bff.QuenchedDonorDecay_get_quenching_table(self))

        @property
        def quenching_rate_map(self):
            return _cube(
                _IMP_bff.QuenchedDonorDecay_get_quenching_rate_map(self))

        @property
        def slow_factor_map(self):
            return _cube(_IMP_bff.QuenchedDonorDecay_get_slow_factor_map(self))

        @property
        def diffusion(self):
            """The Brownian walk, run on first access."""
            return _IMP_bff.QuenchedDonorDecay_get_diffusion(self)

        @property
        def k_quench(self):
            """The quenching rate the dye sees, frame by frame, 1/ns."""
            return _IMP_bff.QuenchedDonorDecay_get_k_quench(self)

        @property
        def photon_trace(self):
            """``(delay_times, emitted)``, sampled on first access."""
            return (_IMP_bff.QuenchedDonorDecay_get_delays(self),
                    _IMP_bff.QuenchedDonorDecay_get_emitted(self).astype(np.uint8))

        @property
        def quantum_yield(self):
            """Emitted photons over excitation events."""
            return _IMP_bff.QuenchedDonorDecay_get_quantum_yield(self)

        @property
        def fluorescence_lifetime(self):
            """The species-averaged lifetime of the emitted photons."""
            return _IMP_bff.QuenchedDonorDecay_get_fluorescence_lifetime(self)

        def update_grids(self):
            """Stamp the quenching-rate and stickiness fields onto the AV grid."""
            _IMP_bff.QuenchedDonorDecay_update_grids(self)
            return self.quenching_rate_map, self.slow_factor_map

        def simulate_photons(self, k_quench=None):
            """Race photons against the quenching rate along the trajectory."""
            _IMP_bff.QuenchedDonorDecay__simulate_photons(
                self, np.zeros(0) if k_quench is None else
                np.ascontiguousarray(np.asarray(k_quench, dtype=np.float64)).ravel())
            return self.photon_trace

        def photons_fused(self):
            """The walk, the rate along it, and the photon race -- in one call.

            Identical to :meth:`simulate_diffusion` then
            :meth:`simulate_photons`, photon for photon. What it skips is the
            **trajectory**: at the default ``t_max`` that is 20 million doubles
            to produce a few thousand photons, and nothing downstream of the
            decay wants the coordinates. Use the split calls when the trajectory
            *is* the point; use this one inside a fitting loop.
            """
            _IMP_bff.QuenchedDonorDecay__photons_fused(self)
            return self.photon_trace

        def lifetime_spectrum(self, n_species=128):
            """The decay as ``(amplitude, rate)`` pairs -- the neutral output.

            Prefer this to :meth:`decay_histogram`: a spectrum carries no bin
            width and no time range, so whatever owns the instrument can
            convolve, bin and add noise on its own terms.
            """
            return _IMP_bff.QuenchedDonorDecay__lifetime_spectrum(
                self, int(n_species))

        def decay_histogram(self, n_bins=4096, time_range=(0.0, 50.0)):
            """Histogram of the photons that were actually emitted.

            **Only emitted photons.** A quenched excitation comes back with
            ``dt = 0``, so histogramming the whole trace piles every non-emitted
            photon into the first bin -- a spike of photons that never existed.
            """
            n_bins = int(n_bins)
            flat = _IMP_bff.QuenchedDonorDecay__decay_histogram(
                self, n_bins, float(time_range[0]), float(time_range[1]))
            return flat[:n_bins + 1], flat[n_bins + 1:].astype(np.int64)

        def fret_rate_trace(self, acceptor, forster_radius=52.0, kappa2=None,
                            r_min=7.0):
            """Per-frame FRET rate against an acceptor.

            *acceptor* may be a :class:`QuenchedDonorDecay` whose walk has been
            run -- then both dyes are resolved in time and the frames are paired
            -- or anything with a ``points`` cloud, which is the fast-acceptor
            limit.
            """
            k2 = _kappa2(kappa2)
            if isinstance(acceptor, QuenchedDonorDecay):
                return _IMP_bff.QuenchedDonorDecay__fret_rate_trace_paired(
                    self, acceptor, float(forster_radius), k2, float(r_min))
            points = np.asarray(getattr(acceptor, "points", acceptor),
                                dtype=np.float64)
            return _IMP_bff.QuenchedDonorDecay__fret_rate_trace_cloud(
                self, np.ascontiguousarray(points[:, :3]).ravel(),
                float(forster_radius), k2, float(r_min))

        def fret_efficiency(self, acceptor, forster_radius=52.0, kappa2=None):
            """``1 - QY_DA / QY_D``, both from the photon Monte-Carlo.

            Taking the ratio of two simulated quantum yields rather than an
            analytic formula keeps the quenching in: the donor is quenched by
            PET in both terms, so what is left is the FRET.
            """
            rates = self.fret_rate_trace(acceptor, forster_radius, kappa2)
            return _IMP_bff.QuenchedDonorDecay_fret_efficiency(
                self, np.ascontiguousarray(np.asarray(rates, dtype=np.float64)).ravel())
    %}
}
