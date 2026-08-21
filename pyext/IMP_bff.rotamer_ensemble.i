/*
 * RotamerEnsemble -- defined here, at the end of swig.i-in, because it
 * subclasses States (a SWIG value type that is defined after all %pythoncode
 * blocks in the main body).
 */

%pythoncode %{
class RotamerEnsemble(States):
    """A rotamer library placed and screened at one labelling site (fps ``R1``).

    A **sibling** of :class:`IMP.bff.AccessibleVolume`, not a subclass of it. It
    used to inherit from the concrete AV, which is why distance code worked for
    rotamers -- by inheritance rather than by design -- and it then had to carry
    a grid it does not have, filling ``density``, ``grid_step`` and
    ``grid_shape`` with empty placeholders. No consumer ever read them: they all
    used ``points``, ``mean_position``, ``n_points`` or ``has_volume``, which is
    the :class:`IMP.bff.States` surface both representations share.

    From ``States``, which is C++: ``points`` (N, 4) chromophore centre +
    weight, ``attachment_point`` (CA), ``orientations`` (the per-rotamer
    transition dipoles), ``position_name``, ``params`` (carries
    ``simulation_type='R1'``, library, temperature, ...). Adds ``atoms``
    (N, n_atoms, 3) in the protein frame, ``atom_names``, ``resnames``,
    ``energies``, ``partition`` (Z) and the ``library`` name.

    ``mu`` is the name the rotamer code uses for the dipoles, and it *is*
    ``States.orientations`` -- a property over the one array, so a caller cannot
    set one and read a stale other. It is why this is a plain class rather than
    a dataclass: the states are a C++ value, and its constructor is the one that
    has to run.
    """

    def __init__(
        self,
        *,
        points,
        attachment_point,
        position_name: str = "",
        params: Optional[dict] = None,
        mu=None,
        atoms=None,
        atom_names: tuple = (),
        resnames: tuple = (),
        energies=None,
        partition: float = 0.0,
        library: str = "",
        chain: str = "",
        residue: int = 0,
    ):
        States.__init__(
            self, points=points, attachment_point=attachment_point,
            orientations=mu, position_name=position_name, params=params or {})
        self.atoms = (np.zeros((0, 0, 3)) if atoms is None
                      else np.asarray(atoms, dtype=np.float64))
        self.atom_names = tuple(atom_names)
        self.resnames = tuple(resnames)
        self.energies = (np.zeros(0) if energies is None
                         else np.asarray(energies, dtype=np.float64))
        self.partition = float(partition)
        self.library = str(library)
        self.chain = str(chain)
        self.residue = int(residue)

    @property
    def mu(self) -> np.ndarray:
        """(N, 3) per-rotamer transition dipoles -- ``States.orientations``."""
        return self.orientations

    @mu.setter
    def mu(self, value):
        self.orientations = value

    def __repr__(self) -> str:
        return (f"RotamerEnsemble({self.n_rotamers} rotamers, "
                f"{self.library!r}, Z={self.partition:.3g})")

    # -- construction ------------------------------------------------------
    @classmethod
    def from_site(
        cls,
        structure,
        chain: Optional[str],
        residue: int,
        library,
        *,
        temperature: float = 298.15,
        electrostatic: bool = False,
        potential: str = "lj",
        ignore_h: bool = True,
        sigma_scaling: float = 0.5,
        epsilon_scaling: float = 1.0,
        frame_index: int = 0,
        position_name: str = "",
    ) -> "RotamerEnsemble":
        """Place ``library`` on ``(chain, residue)`` of ``structure`` and screen it.

        ``structure`` is a PDB / multi-MODEL PDB / RMF path or a frame dict
        from ``load_protein_frames``; ``library`` a registry name
        (``'AlexaFluor 488 C1R cutoff30'``), a path, or a loaded library dict.
        Scoring is FRETpredict's (LJ or Gauss, optional Debye–Hückel; the
        labelled residue and hydrogens are not obstacles).
        """
        frame = _frame(structure, frame_index)
        lib = _library(library)
        ca, n, c = resolve_backbone_site(frame, chain, residue)
        rotamers = transform_library_to_site(lib["coords"], ca, n, c)
        score = compute_rotamer_score(
            rotamers,
            frame["coords"],
            frame["atom_names"],
            frame["resnames"],
            lib["atom_names"],
            lib.get("metadata", {}),
            lib.get("resnames"),
            protein_residue_indices=frame.get("residue_indices"),
            protein_chain_ids=frame.get("chain_ids"),
            site_residue=residue,
            site_chain=chain,
            rotamer_weights=lib.get("weights"),
            temperature=temperature,
            ignore_h=ignore_h,
            electrostatic=electrostatic,
            potential=potential,
            sigma_scaling=sigma_scaling,
            epsilon_scaling=epsilon_scaling,
        )
        metadata = dict(lib.get("metadata", {}) or {})
        names = list(lib["atom_names"])
        lib_res = lib.get("resnames")
        centre_idx = selector_atom_indices(names, metadata.get("r", []), lib_res)[0]
        mu_idx = selector_atom_indices(names, metadata.get("mu", []), lib_res)
        if len(mu_idx) >= 2:
            mu = rotamers[:, mu_idx[1], :] - rotamers[:, mu_idx[0], :]
        else:
            mu = rotamers[:, 1, :] - rotamers[:, 0, :]
        mu = mu / np.linalg.norm(mu, axis=1, keepdims=True)
        centres = rotamers[:, centre_idx, :]
        points = np.hstack([centres, score.weights[:, None]]).astype(np.float64)
        ca_v = np.asarray(ca, dtype=np.float64)
        return cls(
            points=points,
            attachment_point=ca_v.copy(),
            position_name=position_name or f"{chain or ''}{residue}",
            params={
                "simulation_type": SIMULATION_TYPE_R1,
                "library": metadata.get("library_name", metadata.get("name", str(library))),
                "chain": chain or "",
                "residue": int(residue),
                "temperature": float(temperature),
                "electrostatic": bool(electrostatic),
                "potential": potential,
                "ignore_h": bool(ignore_h),
                "sigma_scaling": float(sigma_scaling),
                "epsilon_scaling": float(epsilon_scaling),
                "partition": float(score.partition),
            },
            mu=mu,
            atoms=rotamers,
            atom_names=tuple(names),
            resnames=tuple(lib.get("resnames") or ()),
            energies=np.asarray(score.energies, dtype=np.float64),
            partition=float(score.partition),
            library=str(metadata.get("library_name", metadata.get("name", str(library)))),
            chain=chain or "",
            residue=int(residue),
        )

    # -- accessors ---------------------------------------------------------
    @property
    def centres(self) -> np.ndarray:
        """(N, 3) chromophore centres in the protein frame (Å)."""
        return self.points[:, :3]

    @property
    def weights(self) -> np.ndarray:
        """(N,) normalised Boltzmann × library weights."""
        return self.points[:, 3]

    @property
    def n_rotamers(self) -> int:
        return int(self.points.shape[0])

    # -- pair physics ------------------------------------------------------
    def pair_geometry(self, other: "RotamerEnsemble") -> dict:
        """R_ij, κ²_ij and w_i·w_j against another ensemble (or any AV: κ² = 2/3)."""
        mu_other = getattr(other, "mu", None)
        if mu_other is not None and np.asarray(mu_other).shape[0] != other.points.shape[0]:
            mu_other = None
        return fret_pair_geometry(self.centres, self.weights, other.points[:, :3], other.points[:, 3], self.mu, mu_other)

    def pair_distribution(self, other: "RotamerEnsemble", forster_radius: float, tau0: Optional[float] = None) -> dict:
        """The pair's FRET rate distribution and averages (see ``fret.distance.fret_pair_efficiencies``).

        ``forster_radius`` in Å for κ² = 2/3.
        """
        geometry = self.pair_geometry(other)
        return fret_pair_efficiencies(geometry, forster_radius, tau0)

    def fret_efficiencies(
        self,
        other: "RotamerEnsemble",
        forster_radius: Optional[float] = None,
        *,
        donor: Optional[str] = None,
        acceptor: Optional[str] = None,
    ) -> dict:
        """E_static, E_dynamic1, E_dynamic2 and ⟨κ²⟩ for this (donor) and ``other`` (acceptor).

        Give ``forster_radius`` (Å, κ² = 2/3) or the dye names -- then R0 is
        computed from the spectra at the pair's ⟨κ²⟩, as FRETpredict does.
        """
        geometry = self.pair_geometry(other)
        if forster_radius is None:
            if donor is None or acceptor is None:
                raise ValueError("give forster_radius (A) or donor and acceptor names")
            r0_nm = forster_radius_from_spectra(donor, acceptor, geometry["kappa2_avg"])
            # FRETpredict applies its k2-dependent R0 with the isotropic formula
            # 1/(1 + (2/3/k2)(r/R0)^6); fret_pair_efficiencies expects R0 at k2 = 2/3
            forster_radius = r0_nm * 10.0
            return fret_pair_efficiencies(geometry, forster_radius)
        return fret_pair_efficiencies(geometry, forster_radius)


def rotamer_ensembles_from_fps(
    fps_json,
    structure,
    library_map: Optional[Dict[str, str]] = None,
    *,
    frame_index: int = 0,
    **kwargs,
) -> Dict[str, RotamerEnsemble]:
    """One :class:`RotamerEnsemble` per fps.json position that names a library.

    A position's ``rotamer_library`` / ``library`` field selects the library;
    ``library_map`` (position name → library name) overrides or supplies it
    for AV-only files. Positions without a library are skipped. ``kwargs`` go
    to :meth:`RotamerEnsemble.from_site`.
    """
        # (was: from .fps import ...) -- now in this module

    positions, _distances, _score_sets, _extra = read_fps_json(fps_json)
    frame = _frame(structure, frame_index)
    library_map = dict(library_map or {})
    out: Dict[str, RotamerEnsemble] = {}
    for name, payload in positions.items():
        pos = RotamerPosition.from_payload(name, payload)
        lib_name = library_map.get(name) or pos.library
        if not lib_name:
            continue
        out[name] = RotamerEnsemble.from_site(
            frame, pos.chain, pos.residue, lib_name, position_name=name, **kwargs)
    return out

%}
