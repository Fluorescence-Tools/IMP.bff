/*
 * Sampling a tethered dye: the walk, the photons, and the equilibrium.
 *
 * `sampling.py` was the shapes these kernels' callers want -- a named tuple for
 * the trajectory, a mobility field built from a scalar and a mask, an
 * in-place-accumulated decay curve. The trajectory is a C++ value now and the
 * mobility field is built where the walk is, which is what collapsed two
 * near-identical kernels into one.
 */

IMP_SWIG_VALUE(IMP::bff, DyeDiffusionTrajectory, DyeDiffusionTrajectories);
IMP_SWIG_VALUE(IMP::bff, RotamerLibrary, RotamerLibraries);

%apply(int* IN_ARRAY1, int DIM1) {(int* density, int n_density)};
%apply(double* IN_ARRAY1, int DIM1) {(double* mobility, int n_mobility)};
%apply(double* IN_ARRAY1, int DIM1) {(double* diffusion_map, int n_diffusion_map)};
%apply(int* IN_ARRAY1, int DIM1) {(int* bounds, int n_bounds)};

%rename(_simulate_dye_diffusion) IMP::bff::simulate_dye_diffusion;
%rename(_equilibrium_occupancy) IMP::bff::equilibrium_occupancy;
%rename(_find_reference_rotamer_files) IMP::bff::find_reference_rotamer_files;
// `load_rotamer_library` is taken: `representation/rotamer.py` has one that
// resolves a *registry name* to a file and reads it. This one takes a PDB and a
// trajectory, and the public spelling of it is `load_rotamer_library_dcd`.
%rename(_load_rotamer_library) IMP::bff::load_rotamer_library;

%include "IMP/bff/DyeSampling.h"

%extend IMP::bff::RotamerLibrary {
    %pythoncode %{
        @property
        def coords(self):
            """``(n_frames, n_atoms, 3)`` conformer coordinates, A."""
            return _IMP_bff.RotamerLibrary_get_coords(self).reshape(
                self.n_frames, self.n_atoms, 3)

        @property
        def weights(self):
            """``(n_frames,)``, normalised to sum 1."""
            return _IMP_bff.RotamerLibrary_get_weights(self)

        def __getitem__(self, key):
            """The dictionary surface the loader's result had."""
            return {"coords": self.coords, "weights": self.weights,
                    "atom_names": list(self.atom_names)}[key]

        def get(self, key, default=None):
            try:
                return self[key]
            except KeyError:
                return default
    %}
}

%extend IMP::bff::DyeDiffusionTrajectory {
    %pythoncode %{
        @property
        def xyz(self):
            """``(n_frames, 3)`` dye-centre positions, relative to the grid
            anchor: add the attachment point to place them in the structure."""
            return _IMP_bff.DyeDiffusionTrajectory_get_xyz(self).reshape(-1, 3)

        @property
        def accepted(self):
            """``(n_frames,)`` uint8 -- whether each step was taken."""
            return _IMP_bff.DyeDiffusionTrajectory_get_accepted(self).astype(
                np.uint8)

        @property
        def n_frames(self):
            return self.get_n_frames()

        @property
        def n_accepted(self):
            return self.get_n_accepted()

        @property
        def n_rejected(self):
            return self.get_n_rejected()

        @property
        def acceptance_ratio(self):
            return self.get_acceptance_ratio()

        def __len__(self):
            return self.get_n_frames()
    %}
}

%pythoncode %{
def simulate_dye_diffusion(density, slow_density=None, dg=0.5, t_max=10000.0,
                           t_step=0.002, D=40.0, slow_fact=1.0,
                           random_seed=None):
    """One Brownian trajectory of the dye centre in its accessible volume.

    :param density: ``(ng, ng, ng)`` binary occupancy of the volume.
    :param slow_density: binary contact grid, used with a **scalar**
        ``slow_fact``; ignored when ``slow_fact`` is itself a grid.
    :param dg: voxel edge, A.
    :param t_max, t_step: total time and step, ns.
    :param D: diffusion coefficient, A^2/ns, in the standard convention
        (``<dx^2> = 2 D t`` per component) -- the same D
        :class:`GridDiffusionSolver` takes. The default 40 is free Alexa488 in
        water and is an uncalibrated transferable starting value, as every
        document in this stack that mentions it says.
    :param slow_fact: a scalar factor applied inside ``slow_density``, or a
        ``(ng, ng, ng)`` per-voxel factor grid from :func:`slow_factor_grid`.
    :param random_seed: ``None`` draws freely.
    """
    occupancy = np.ascontiguousarray(np.asarray(density, dtype=np.int32))
    ng = int(occupancy.shape[0])
    slow = np.asarray(slow_fact, dtype=np.float64)

    # Mobility is a *field*: one scaling per voxel. The scalar-plus-mask form is
    # one way to build it.
    if slow.ndim == 3:
        mobility = np.ascontiguousarray(slow, dtype=np.float64)
    elif float(slow) == 1.0 or slow_density is None:
        mobility = np.empty(0, dtype=np.float64)     # uniform medium
    else:
        mask = np.asarray(slow_density, dtype=bool)
        mobility = np.where(mask, float(slow), 1.0)

    return _IMP_bff._simulate_dye_diffusion(
        occupancy.ravel(), np.ascontiguousarray(mobility).ravel(), ng,
        float(dg), float(t_max), float(t_step), float(D),
        -1 if random_seed is None else int(random_seed))


def simulate_photon_trace(n_ph, k_quench, t_step=0.01, tau0=0.25,
                          random_seed=None):
    """``(delay_times, emitted)`` for ``n_ph`` excitation events.

    ``emitted`` is a uint8 flag; the delay time is 0 where the dye was quenched
    before emitting. The flag comes back inside the returned array because a
    SWIG out-parameter is a wrapper numpy walks one element at a time -- 13.7 ms
    of a 152 ms trace, measured.
    """
    flat = photon_trace(
        int(n_ph),
        np.ascontiguousarray(np.asarray(k_quench, dtype=np.float64).ravel()),
        float(t_step), float(tau0),
        -1 if random_seed is None else int(random_seed))
    packed = np.asarray(flat, dtype=np.float64).reshape(-1, 2)
    return np.ascontiguousarray(packed[:, 0]), packed[:, 1].astype(np.uint8)


def simulate_quenched_decay(n_curves, decay, dt_tac, k_quench, t_step, tau0,
                            random_seed=None):
    """Add a trajectory-driven fluorescence decay curve into ``decay``.

    The smooth decay of the same model :func:`simulate_photon_trace` samples
    photon by photon -- no shot noise, and far cheaper when only the curve is
    wanted.
    """
    decay += np.asarray(quenched_decay(
        int(n_curves), int(np.asarray(decay).shape[0]), float(dt_tac),
        np.ascontiguousarray(np.asarray(k_quench, dtype=np.float64).ravel()),
        float(t_step), float(tau0),
        -1 if random_seed is None else int(random_seed)), dtype=np.float64)


def equilibrium_occupancy(diffusion_map, bounds, flux_form="smoluchowski"):
    """The stationary occupancy of the volume, in closed form.

    ``"smoluchowski"`` (the default, and the physics) gives ``p_eq`` uniform on
    the accessible domain, independent of D; ``"ito"`` gives ``p_eq ~ 1/D``, so
    the dye piles up wherever it moves slowly. Equilibrium is thermodynamics and
    mobility is kinetics: a dye slowed by friction with no attractive
    interaction is still found uniformly across its volume.
    """
    d = np.ascontiguousarray(np.asarray(diffusion_map, dtype=np.float64))
    b = np.ascontiguousarray(np.asarray(bounds)).astype(np.int32)
    out = _IMP_bff._equilibrium_occupancy(d.ravel(), b.ravel(), str(flux_form))
    return np.asarray(out, dtype=np.float64).reshape(d.shape)


def find_reference_rotamer_files(lib_dir, dye_name, cutoff=30):
    """``(pdb, trajectory, weights)`` for a reference dye+linker name.

    ``weights`` is ``None`` when there is no weights file, which is how the
    loader spells "weight the frames uniformly".
    """
    pdb, traj, weights = _IMP_bff._find_reference_rotamer_files(
        str(lib_dir), str(dye_name), int(cutoff))
    return pdb, traj, (weights or None)


def load_rotamer_library_dcd(pdb_path, dcd_path, weights_path=None,
                             max_frames=None):
    """A reference rotamer library from a PDB plus a trajectory.

    The name says DCD for the format it was written for; BinaryCIF is what this
    package stores, and what the shipped libraries are.
    """
    return _IMP_bff._load_rotamer_library(
        str(pdb_path), str(dcd_path),
        "" if weights_path is None else str(weights_path),
        -1 if max_frames is None else int(max_frames))


def apply_rotamer_coordinates(dye_hier, coords):
    """Put one rotamer's coordinates onto a dye hierarchy, in place."""
    import IMP.algebra
    import IMP.atom
    import IMP.core

    atoms = list(IMP.atom.get_by_type(dye_hier, IMP.atom.ATOM_TYPE))
    arr = np.asarray(coords, dtype=np.float64)
    if arr.shape[0] != len(atoms):
        raise ValueError(
            f"Atom count mismatch: coords={arr.shape[0]} "
            f"hierarchy={len(atoms)}")
    for a, c in zip(atoms, arr):
        IMP.core.XYZ(a).set_coordinates(
            IMP.algebra.Vector3D(float(c[0]), float(c[1]), float(c[2])))
%}
