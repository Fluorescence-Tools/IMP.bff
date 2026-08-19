"""@namespace IMP.bff
Restraints for handling distances between accessible volumes.
"""

from __future__ import annotations, print_function

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple
import json
import os
import pathlib
import typing

import numpy as np

import IMP
import IMP.algebra
import IMP.atom
import IMP.bayesianem
import IMP.bff
import IMP.bff.tools
import IMP.container
import IMP.core
import IMP.isd
import IMP.pmi.restraints
import IMP.pmi.tools
import IMP.rmf

__all__ = [
    'estimate_position_uncertainty',
]

# --------------------------------------------------------------------------
# AVNetworkRestraint
# --------------------------------------------------------------------------
"""@namespace IMP.bff
Restraints for handling distances between accessible volumes.
"""

class AVMeanDistanceRestraint(IMP.Restraint):

    def __init__(
            self,
            m: IMP.Model,
            av1: IMP.bff.AV, av2: IMP.bff.AV,
            dist: IMP.bff.AVPairDistanceMeasurement,
            sigma: float,
            weight: float = 1.0
    ):
        IMP.Restraint.__init__(self, m, "BiStableDistanceRestraint %1%")
        self.dist = dist
        self.av1 = av1
        self.av2 = av2

        self.sigma = sigma
        self.weight = weight
        self.d1 = IMP.core.XYZ(av1)
        self.d2 = IMP.core.XYZ(av2)
        forster_radius = dist.forster_radius
        distance_range = (1, 2.5 * forster_radius)
        self.dc = IMP.bff.tools.FRETDistanceConverter(
            forster_radius=forster_radius,
            sigma=sigma,
            distance_range=distance_range
        )
        self.particle_list = [av1.get_particle(), av2.get_particle()]

    def unprotected_evaluate(self, da):
        d_exp = self.dist
        av1 = self.av1
        av2 = self.av2
        r: IMP.algebra.Vector3D = \
            IMP.core.XYZ(av1).get_coordinates() - \
            IMP.core.XYZ(av2).get_coordinates()
        d_mp = r.get_magnitude()
        d_mod = self.dc(d_mp, d_exp.distance_type)
        score = d_exp.score_model(d_mod)
        return score

    def do_get_inputs(self):
        return self.particle_list


class AVNetworkRestraintWrapper(IMP.pmi.restraints.RestraintBase):

    @staticmethod
    def add_used_dyes_to_rb(used_avs: typing.List[IMP.bff.AV]):
        for dk in used_avs:
            dye = used_avs[dk]
            p_dye = dye.get_particle()
            # The coordinates of an AV are the mean AV the density map. Thus, the position of the
            # AV changes when the AV is resampled.
            dye.resample()

            # dye_xyz = IMP.core.XYZ(p_dye)
            # print("Adding dye to RB:")
            # print("-- Position key:", dk)
            # print("-- Label parameter:", dye)
            # print("-- Mean dye position:", dye_xyz)

            p_att = dye.get_source()
            if IMP.core.RigidBodyMember.get_is_setup(p_att):
                rbm = IMP.core.RigidBodyMember(p_att)
                rb: IMP.core.RigidBody = rbm.get_rigid_body()
                rb.add_member(p_dye)

    def add_xyz_mass_to_avs(self):
        """
        IMP_NEW(Particle,p2,(m,"p2"));
        IMP::core::XYZR d2=IMP::core::XYZR::setup_particle(
        m,p2->get_index(),IMP::algebra::Sphere3D(
        IMP::algebra::Vector3D(1.0,4.0,6.0),1.0));
        atom::Mass mm2 = atom::Mass::setup_particle(p2, 30.0);
        rbps.push_back(d2);
        """
        for ak in self.used_avs:
            av: IMP.bff.AV = self.used_avs[ak]
            r_mean = max(av.get_radii())
            # add radius
            av_d = IMP.core.XYZR.setup_particle(av)
            av_d.set_radius(r_mean * 1.0)
            # add Mass
            av_m = IMP.atom.Mass.setup_particle(av, 0.1)
            av_m.set_mass(r_mean * 2.0)

    """Restraint for Accessible Volume (AV) decorated particles

    The AVs of the decorated particles are recomputed when the
    score is evaluated. Computing an AV (searching for the points
    that are accessible is computationally costly (expensive
    restraint).
    """
    def __init__(
            self,
            hier: IMP.atom.Hierarchy,
            fps_json_fn: str,
            score_set: str = "",
            weight: float = 1.0,
            mean_position_restraint: bool = False,
            sigma_DA: float = 6.0,
            label: str = "AVNetworkRestraint",
            occupy_volume: bool = True
    ):
        """

        :param hier:
        :param fps_json_fn:
        :param score_set:
        :param weight:
        :param mean_position_restraint:
        :param sigma_DA:
        :param label:
        :param occupy_volume:
        """
        # some parameters
        m = hier.get_model()
        self.mdl: IMP.Model = m
        self.hier: IMP.atom.Hierarchy = hier
        super(AVNetworkRestraintWrapper, self).__init__(m, label=label, weight=weight)

        self.model_ps = []
        self.model_ps += [k.get_particle() for k in IMP.atom.get_leaves(hier)]

        name = self.name
        self.mean_position_restraint = mean_position_restraint
        if pathlib.Path(fps_json_fn).is_file():
            self.av_network_restraint = IMP.bff.AVNetworkRestraint(
                hier,
                fps_json_fn,
                name,
                score_set
            )
        else:
            raise FileNotFoundError("{}".format(fps_json_fn))
        self.rs = IMP.RestraintSet(m, 'AVNetworkRestraint')
        self.used_avs = dict([(v.get_name(), v) for v in self.av_network_restraint.get_used_avs()])
        if not self.mean_position_restraint:
            self.rs.add_restraint(self.av_network_restraint)
        else:
            self.used_distances = self.av_network_restraint.get_used_distances()
            self.add_used_dyes_to_rb(self.used_avs)
            for dk in self.used_distances:
                d_exp = self.used_distances[dk]
                av1 = self.used_avs[d_exp.position_1]
                av2 = self.used_avs[d_exp.position_2]
                r = AVMeanDistanceRestraint(m, av1, av2, d_exp, sigma=sigma_DA)
                self.rs.add_restraint(r)
        if occupy_volume:
            self.add_xyz_mass_to_avs()
        self.set_weight(weight)

    def evaluate(self):
        """Evaluate the score of the restraint."""
        return self.rs.unprotected_evaluate(None) * self.weight

    def add_to_model(self, add_to_rmf=True):
        IMP.pmi.tools.add_restraint_to_model(self.mdl, self.rs,
                                             add_to_rmf=add_to_rmf)


# --------------------------------------------------------------------------
# uncertainty
# --------------------------------------------------------------------------
"""Model precision from repeated docking (FPS-style positional uncertainty).

Given the best-scoring structures of several independent docking runs, superpose
them on the fixed (reference) body and measure how much each atom of the mobile
body wanders — the per-atom RMSF. This is the FPS "precision of the model":
a mean structure whose B-factor column carries the positional uncertainty, plus
a per-atom CSV.
"""

def _read_pdb_atoms(path: str):
    """Return ``(lines, xyz)`` for ATOM/HETATM records of a PDB file."""
    lines: List[str] = []
    xyz: List[tuple] = []
    chains: List[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                lines.append(line.rstrip("\n"))
                xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                chains.append(line[21])
    return lines, np.asarray(xyz, dtype=float), np.asarray(chains)


def _kabsch(mobile: np.ndarray, target: np.ndarray):
    """Rigid transform (R, t) minimising ``||R·mobile + t - target||``."""
    mc, tc = mobile.mean(0), target.mean(0)
    h = (mobile - mc).T @ (target - tc)
    u, _s, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return r, tc - r @ mc


def estimate_position_uncertainty(
    pdb_paths: Sequence[str],
    fixed_chains: Sequence[str],
    out_pdb: Optional[str] = None,
    out_csv: Optional[str] = None,
) -> Dict:
    """Superpose docked models on the fixed body and report per-atom RMSF.

    Parameters
    ----------
    pdb_paths : sequence of str
        Best-scoring PDB per docking run (same atom ordering — they are the same
        structure docked from different starts).
    fixed_chains : sequence of str
        Chain identifiers of the fixed/reference body used for superposition.
    out_pdb : str, optional
        Write the mean structure with B-factor = per-atom RMSF here.
    out_csv : str, optional
        Write a per-atom ``index,chain,rmsf`` table here.

    Returns
    -------
    dict
        ``{n_models, rmsf_mean, rmsf_max, mobile_rmsf_mean, uncertainty_pdb,
        uncertainty_csv}``.
    """
    paths = [p for p in pdb_paths if p and os.path.exists(p)]
    if len(paths) < 2:
        return {"n_models": len(paths), "rmsf_mean": float("nan"),
                "rmsf_max": float("nan"), "mobile_rmsf_mean": float("nan"),
                "uncertainty_pdb": None, "uncertainty_csv": None}

    lines0, xyz0, chains = _read_pdb_atoms(paths[0])
    fixed_mask = np.isin(chains, list(fixed_chains))
    if not fixed_mask.any():
        fixed_mask = np.ones(len(chains), dtype=bool)  # no fixed body: align on all

    aligned = [xyz0]
    for p in paths[1:]:
        _lines, xyz, _ch = _read_pdb_atoms(p)
        if xyz.shape != xyz0.shape:
            continue  # skip models with a different atom count
        r, t = _kabsch(xyz[fixed_mask], xyz0[fixed_mask])
        aligned.append((r @ xyz.T).T + t)

    stack = np.stack(aligned)                  # (n_models, n_atoms, 3)
    mean = stack.mean(0)
    rmsf = np.sqrt(((stack - mean) ** 2).sum(-1).mean(0))  # per-atom

    if out_pdb:
        _write_bfactor_pdb(lines0, mean, rmsf, out_pdb)
    if out_csv:
        import csv
        with open(out_csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["atom_index", "chain", "rmsf"])
            for i, (c, v) in enumerate(zip(chains, rmsf)):
                w.writerow([i, c, round(float(v), 3)])

    mobile = rmsf[~fixed_mask] if (~fixed_mask).any() else rmsf
    return {
        "n_models": len(aligned),
        "rmsf_mean": float(rmsf.mean()),
        "rmsf_max": float(rmsf.max()),
        "mobile_rmsf_mean": float(mobile.mean()),
        "uncertainty_pdb": out_pdb if out_pdb else None,
        "uncertainty_csv": out_csv if out_csv else None,
    }


def _write_bfactor_pdb(lines, coords, bfactors, out_pdb) -> None:
    """Rewrite ATOM lines with mean coordinates and RMSF in the B-factor column."""
    with open(out_pdb, "w") as fh:
        for line, (x, y, z), b in zip(lines, coords, bfactors):
            b = min(999.99, float(b))
            fh.write(f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:60]}{b:6.2f}{line[66:]}\n")
        fh.write("END\n")


# --------------------------------------------------------------------------
# engine
# --------------------------------------------------------------------------
"""Lightweight rigid-body data structures (legacy).

The hand-rolled spring/Verlet docking engine that once lived here has been
removed in favour of :mod:`.imp_engine` (IMP + IMP.bff). Only the plain data
classes used by OLGA-style evaluators remain: :class:`RigidBody`,
:class:`DistanceRestraint`, :class:`SpringParameters`. ``SpringParameters``
mirrors the legacy C# FPS simulation knobs, which flrCIF records as
``_flr_FPS_global_parameter.sim_*``.
"""

# ---------------------------------------------------------------------------
# Parameter container
# ---------------------------------------------------------------------------

@dataclass
class SpringParameters:
    """Simulation parameters mirroring FPS ``FPSParameters``."""
    viscosity_factor: float = 0.85
    time_step_factor: float = 0.005
    max_iterations: int = 50000
    max_force: float = 100.0
    clash_tolerance: float = 0.0
    k_clash: float = 10.0
    rkT: float = 1.0
    E_tolerance: float = 1e-4
    K_tolerance: float = 1e-4
    F_tolerance: float = 1e-4
    T_tolerance: float = 1e-4
    optimize_selected: int = 0  # 0=all, 1=selected, 2=selected-then-all


# ---------------------------------------------------------------------------
# Rigid Body
# ---------------------------------------------------------------------------

@dataclass
class RigidBody:
    """A rigid body that can be translated and rotated as a unit.

    Attributes stored in the *body* (local) frame.
    """
    name: str
    atoms_local: np.ndarray  # (N, 4) xyzr in the local frame (com = 0)
    com: np.ndarray  # (3,) center of mass in global frame
    rotation: np.ndarray  # (3, 3) rotation from local → global
    translation: np.ndarray  # (3,) alias for com (for clarity)
    velocity: np.ndarray = field(default_factory=lambda: np.zeros(3))
    angular_velocity: np.ndarray = field(default_factory=lambda: np.zeros(3))
    mass: float = 1.0
    inertia: np.ndarray = field(default_factory=lambda: np.eye(3))

    def global_coords(self) -> np.ndarray:
        """Return (N, 3) atom coordinates in the global frame."""
        return self.atoms_local[:, :3] @ self.rotation.T + self.com

    def global_xyzr(self) -> np.ndarray:
        """Return (N, 4) xyzr in the global frame."""
        coords = self.global_coords()
        return np.column_stack([coords, self.atoms_local[:, 3]])

    def apply_force_at_point(
        self,
        force: np.ndarray,
        point: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Apply a force at a point on the body.

        Returns (linear_force_increment, torque_increment).
        """
        r = point - self.com
        torque = np.cross(r, force)
        return force, torque

    def random_shake(self, translation_scale: float = 5.0, rotation_scale: float = 0.5):
        """Randomly perturb position and orientation (for initial clash removal)."""
        self.com += np.random.uniform(-translation_scale, translation_scale, size=3)
        angle = np.random.uniform(-rotation_scale, rotation_scale)
        axis = np.random.randn(3)
        axis /= np.linalg.norm(axis) + 1e-12
        rot = _rotation_matrix(axis, angle)
        self.rotation = self.rotation @ rot


# ---------------------------------------------------------------------------
# Distance restraint (harmonic spring)
# ---------------------------------------------------------------------------

@dataclass
class DistanceRestraint:
    """A single FRET distance restraint between two AV positions on two bodies.

    Positions are stored as offsets from each body's COM so they move
    with the body during docking.
    """
    name: str
    body_a: int
    offset_a: np.ndarray  # AV mean position relative to body A's COM
    body_b: int
    offset_b: np.ndarray  # AV mean position relative to body B's COM
    distance_exp: float
    error_neg: float
    error_pos: float
    distance_type: str = "RDAMean"
    forster_radius: float = 52.0
    active: bool = True
    position_name_a: str = ""
    position_name_b: str = ""
    sigma_rda: float = 0.0
    convfun: Optional[np.ndarray] = None
    transfer_function_type: str = "Polynomial"

    def global_position_a(self, bodies) -> np.ndarray:
        """Return AV mean position on body A in the global frame."""
        ba = bodies[self.body_a]
        return ba.com + ba.rotation @ self.offset_a

    def global_position_b(self, bodies) -> np.ndarray:
        """Return AV mean position on body B in the global frame."""
        bb = bodies[self.body_b]
        return bb.com + bb.rotation @ self.offset_b

    def get_effective_distance(self, rmp: float) -> float:
        """Apply the transfer function to the Rmp distance.

        Parameters
        ----------
        rmp : float
            Distance between mean positions.

        Returns
        -------
        float
            Effective distance.
        """
        import IMP.bff.representation.distance as _dist
        if self.distance_type == "Rmp" or self.transfer_function_type == "None":
            return rmp
        elif self.transfer_function_type == "Gaussian":
            if self.sigma_rda > 0.0:
                return float(_dist.gaussian_rmp_to_rda_mean(rmp, self.sigma_rda))
            return rmp
        elif self.transfer_function_type == "Polynomial":
            if self.convfun is not None:
                return float(_dist.polynomial_transfer(rmp, self.convfun))
            if self.sigma_rda > 0.0:
                return float(_dist.gaussian_rmp_to_rda_mean(rmp, self.sigma_rda))
            return rmp
        return rmp


# ---------------------------------------------------------------------------
# Rotation helpers
# ---------------------------------------------------------------------------

def _rotation_matrix(axis: np.ndarray, angle: float) -> np.ndarray:
    """Build a (3,3) rotation matrix from an axis-angle pair."""
    c = np.cos(angle)
    s = np.sin(angle)
    t = 1.0 - c
    x, y, z = axis
    return np.array([
        [t * x * x + c, t * x * y - s * z, t * x * z + s * y],
        [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
        [t * x * z - s * y, t * y * z + s * x, t * z * z + c],
    ])


def _rotate_vector(v: np.ndarray, axis: np.ndarray, angle: float) -> np.ndarray:
    """Rotate a vector by an axis-angle (Rodrigues' formula)."""
    c = np.cos(angle)
    s = np.sin(angle)
    return v * c + np.cross(axis, v) * s + axis * np.dot(axis, v) * (1.0 - c)


# --------------------------------------------------------------------------
# xlink
# --------------------------------------------------------------------------
class XLinkScore(typing.TypedDict):
    """
    Attributes
    ----------

    individual: list
        The scores to all potential cross linking partners. The score is the
        fraction of distances of all path that are shorter or equal to the
        cross-linker length

    total: np.ndarray
        The sum of scores to all potential cross linking partners. The score is the
        fraction of distances of all path that are shorter or equal to the
        cross-linker length

    """
    individual: list
    total: np.ndarray


class ScoreXlinkSurfaceDistance(IMP.pmi.restraints.RestraintBase):
    """
    :param linker_length:
    :param linker_width:
    :param radius:
    :param simulation_grid_spacing:
    :param min_points: the minimum number of points the the computed volume. If
        there are less points the site is not accessible and the score will be
        zero
    :param verbose:
    :return:
    """

    model: IMP.Model
    obstacles: np.ndarray
    xlinks: typing.Dict[int, typing.Dict]
    linker_length: float = 20.0
    linker_width: float = 2.0
    radius: float = 1.0
    simulation_grid_spacing: float = 3.5
    min_points: int = 100
    verbose: bool = False
    name_map: dict = None

    def compute_scores(self, **kwargs) -> XLinkScore:
        linker_length = kwargs.get('linker_length', self.linker_length)
        linker_width = kwargs.get('linker_width', self.linker_width)
        radius = kwargs.get('radius', self.radius)
        simulation_grid_spacing = kwargs.get('simulation_grid_spacing', self.simulation_grid_spacing)
        verbose = kwargs.get('verbose', self.verbose)
        min_points = kwargs.get('min_points', self.min_points)
        obstacles = self.obstacles
        scores = list()
        total_scores = list()
        for xlink_key in self.xlinks:
            xlink = self.xlinks[xlink_key]
            protein_1 = xlink['protein_1']
            protein_2 = xlink['protein_2']
            residue_1 = xlink['residue_1']
            residue_2 = xlink['residue_2']
            attachment_idx_1 = np.where((obstacles['res_id'] == residue_1) & (obstacles['protein_name'] == protein_1))[0]
            attachment_idx_2 = np.where((obstacles['res_id'] == residue_2) & (obstacles['protein_name'] == protein_2))[0]
            x_links_scores = []
            for idx_1 in attachment_idx_1:
                origin = obstacles['xyz'][idx_1]
                for idx_2 in attachment_idx_2:
                    target = obstacles['xyz'][idx_2]
                    if np.linalg.norm(origin - target) > linker_length:
                        if verbose:
                            print("Eucledian distance > linker length")
                        x_links_scores.append(0.0)
                    else:
                        av = get_path_length(
                            obstacles=self.obstacles,
                            idx=idx_1,
                            linker_length=linker_length,
                            linker_width=linker_width,
                            radius=radius,
                            simulation_grid_spacing=simulation_grid_spacing
                        )
                        if len(np.array(av.points()).T) < min_points:
                            if verbose:
                                print(protein_1, residue_1, "is not accessible")
                            x_links_scores.append(0.0)
                        else:
                            if verbose:
                                print("Eucledian distance < linker length.. Testing surface distance")
                            points = av.points()
                            xyz = points[0:3]
                            dist = points[3]
                            dist_eq = np.linalg.norm(xyz.T - target, axis=1) + dist
                            shorter = np.sum(dist_eq < (linker_length + radius))
                            score = shorter / len(dist_eq)
                            x_links_scores.append(score)
            indiviudal_scores = np.array(x_links_scores)
            scores.append(indiviudal_scores)
            total_score = indiviudal_scores.sum()
            total_scores.append(total_score)
            if verbose:
                print(xlink_key, ":", xlink, "score:", total_score)
        re: XLinkScore = {
            'total': np.array(total_scores, dtype=np.float),
            'individual': scores
        }
        return re

    def __init__(
            self,
            root_hier: IMP.atom.Hierarchy = None,
            verbose: bool = False,
            xlink_settings_file: str = '',
            rmf_file: str = '',
            weight: float = 1.0,
            label: str = 'XLinkSurfaceDistance'
    ):
        # create a new ScoreXlinkRMF object with predefined settings
        base_dir = pathlib.Path(xlink_settings_file).parent
        with open(xlink_settings_file, 'r') as fp:
            xlink_settings = json.load(fp)
        if root_hier is None:
            model = IMP.Model()
        else:
            model = root_hier.get_model()
        super().__init__(
            model,
            weight=weight,
            label=label
        )
        self.model = model
        keys, formats = list(zip(*OBSTACLES_KEYS_FORMATS))
        obstacles = np.zeros(0, dtype={
                'names': keys,
                'formats': formats
            }
        )
        self.obstacles = obstacles
        self.xlinks = dict()
        self.linker_length = xlink_settings['linker_length']
        self.linker_width = xlink_settings['linker_width']
        self.radius = xlink_settings['radius']
        self.simulation_grid_spacing = xlink_settings['simulation_grid_spacing']
        self.min_points = xlink_settings['min_points']
        self.verbose = verbose
        xlink_file = xlink_settings['xlink_file']
        if pathlib.Path(xlink_file).is_file():
            self.xlinks = IMP.bff.tools.read_xlink_table(
                fn=xlink_file
            )
        elif (base_dir / xlink_file).is_file():
            self.xlinks = IMP.bff.tools.read_xlink_table(
                fn=str(base_dir / xlink_file)
            )
        else:
            raise FileNotFoundError("Could not find XL file %s" % xlink_file)
        if pathlib.Path(rmf_file).is_file():
            self.obstacles = get_obstacles(
                rmf_file=rmf_file,
                model=self.model,
                name_map=self.name_map
            )
        # Add custom metadata (will be saved in RMF output)
        self.rs.filename = xlink_settings_file

    def __call__(
            self,
            rmf_file: str = '',
            xlink_file: str = '',
            frame_index: int = 0,
            hier: IMP.atom.Hierarchy = None,
            *args, **kwargs
    ):
        if pathlib.Path(xlink_file).is_file():
            self.xlinks = IMP.bff.tools.read_xlink_table(
                fn=xlink_file
            )
        if pathlib.Path(rmf_file).is_file():
            self.obstacles = get_obstacles(
                rmf_file=rmf_file,
                hier=hier,
                model=self.model,
                frame_index=frame_index,
                name_map=self.name_map
            )
        return self.compute_scores(**kwargs)

    def get_score(self) -> float:
        total_score = float(np.sum(self()['total']))
        return total_score
