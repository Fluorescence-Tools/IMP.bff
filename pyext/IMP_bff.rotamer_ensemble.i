/*
 * RotamerEnsemble -- defined here, at the end of swig.i-in, because it
 * subclasses States (a SWIG value type that is defined after all %pythoncode
 * blocks in the main body).
 */

%pythoncode %{
import json as _json

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
        p = IMP.bff.MapStringString()
        for k, v in (params or {}).items():
            p[k] = v if isinstance(v, str) else str(v)
        States.__init__(
            self, points=points, attachment_point=attachment_point,
            orientations=mu, position_name=position_name, params=p)
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

    # -- the States surface, as the attribute names callers have always used ---
    # The base `States` is a C++ value behind get_*()/set_*()`; these properties
    # restore the dataclass attribute surface this subclass used (points,
    # orientations, params, ...) on top of it.
    @property
    def points(self) -> np.ndarray:
        return np.asarray(self.get_points(), dtype=np.float64).reshape(-1, 4)

    @points.setter
    def points(self, value):
        self.set_points(np.ascontiguousarray(
            np.asarray(value, dtype=np.float64)).ravel())

    @property
    def attachment_point(self) -> np.ndarray:
        return np.asarray(self.get_attachment_point(), dtype=np.float64)

    @attachment_point.setter
    def attachment_point(self, value):
        self.set_attachment_point(np.asarray(value, dtype=np.float64).ravel())

    @property
    def orientations(self) -> np.ndarray:
        o = np.asarray(self.get_orientations(), dtype=np.float64)
        return o.reshape(-1, 3) if o.size else o

    @orientations.setter
    def orientations(self, value):
        self.set_orientations(np.ascontiguousarray(
            np.asarray(value, dtype=np.float64)).ravel())

    @property
    def params(self) -> dict:
        return dict(self.get_params())

    @params.setter
    def params(self, value):
        m = IMP.bff.MapStringString()
        for k, v in (value or {}).items():
            m[k] = v if isinstance(v, str) else str(v)
        self.set_params(m)

    @property
    def mean_position(self) -> np.ndarray:
        return np.asarray(self.get_mean_position(), dtype=np.float64)

    @property
    def position_name(self) -> str:
        return self.get_position_name()

    @position_name.setter
    def position_name(self, value):
        self.set_position_name(str(value))

    @property
    def n_points(self) -> int:
        return self.get_n_points()

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
    @staticmethod
    def _pts4(obj):
        """The `(N, 4)` points of any `States` (ensemble or AV value)."""
        if isinstance(obj, RotamerEnsemble):
            return obj.points
        p = np.asarray(obj.get_points(), dtype=np.float64)
        return p.reshape(-1, 4)

    @staticmethod
    def _mu_of(obj):
        """The `(N, 3)` dipoles of any `States`, or None."""
        mu = getattr(obj, "mu", None)
        if mu is not None:
            return mu
        o = np.asarray(obj.get_orientations(), dtype=np.float64)
        return o.reshape(-1, 3) if o.size else None

    def pair_geometry(self, other: "RotamerEnsemble") -> dict:
        """R_ij, κ²_ij and w_i·w_j against another ensemble (or any AV: κ² = 2/3)."""
        op = self._pts4(other)
        om = self._mu_of(other)
        if om is not None and om.shape[0] != op.shape[0]:
            om = None
        g = fret_pair_geometry(self.centres, self.weights, op[:, :3], op[:, 3], self.mu, om)
        n1, n2 = g.n1, g.n2
        return {
            "R": np.asarray(g.R, dtype=np.float64).reshape(n1, n2),
            "kappa2": np.asarray(g.kappa2, dtype=np.float64).reshape(n1, n2),
            "weight": np.asarray(g.weight, dtype=np.float64).reshape(n1, n2),
            "kappa2_avg": g.kappa2_avg,
        }

    def pair_distribution(self, other: "RotamerEnsemble", forster_radius: float, tau0: Optional[float] = None) -> dict:
        """The pair's FRET rate distribution and averages (see ``fret.distance.fret_pair_efficiencies``).

        ``forster_radius`` in Å for κ² = 2/3.
        """
        g = self.pair_geometry(other)
        op = self._pts4(other)
        om = self._mu_of(other)
        if om is not None and om.shape[0] != op.shape[0]:
            om = None
        eff = fret_pair_efficiencies(
            fret_pair_geometry(
                self.centres, self.weights, op[:, :3], op[:, 3], self.mu, om),
            forster_radius, tau0)
        n1, n2 = eff.n1, eff.n2
        return {
            "R": np.asarray(eff.R, dtype=np.float64).reshape(n1, n2),
            "kappa2": np.asarray(eff.kappa2, dtype=np.float64).reshape(n1, n2),
            "weight": np.asarray(eff.weight, dtype=np.float64).reshape(n1, n2),
            "k_fret": np.asarray(eff.get_k_fret(), dtype=np.float64).reshape(n1, n2),
            "E": np.asarray(eff.get_E(), dtype=np.float64).reshape(n1, n2),
            "rate_ratio": np.asarray(eff.get_rate_ratio(), dtype=np.float64).reshape(n1, n2),
            "static": eff.static_efficiency,
            "dynamic1": eff.dynamic1,
            "dynamic2": eff.dynamic2,
            "kappa2_avg": eff.kappa2_avg,
        }

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
        op = self._pts4(other)
        om = self._mu_of(other)
        if om is not None and om.shape[0] != op.shape[0]:
            om = None
        eff = fret_pair_efficiencies(
            fret_pair_geometry(
                self.centres, self.weights, op[:, :3], op[:, 3], self.mu, om),
            forster_radius)
        return {"static": eff.static_efficiency, "dynamic1": eff.dynamic1,
                "dynamic2": eff.dynamic2, "kappa2_avg": eff.kappa2_avg,
                "forster_radius_nm": forster_radius / 10.0}


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
    positions = _json.loads(read_fps_json(str(fps_json)).positions)
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
