import IMP
import IMP.core
import IMP.atom
import IMP.bff

import typing
import numpy as np
import scipy.stats


def display_mean_av_positions(
        used_avs: typing.List[IMP.bff.AV],
        dye_radius: float = 3.5
):
    """Decorate the accessible volume particles
    with a mass and radii, so that they appear
    in the rmf file"""
    for dye in used_avs:
        # Add radius and mass to dye
        p_dye = dye.get_particle()
        if not IMP.core.XYZR.get_is_setup(p_dye):
            p_dye = IMP.core.XYZR.setup_particle(p_dye)
            p_dye.set_radius(dye_radius)
        else:
            p_dye = IMP.core.XYZR(p_dye)
        if not IMP.atom.Mass.get_is_setup(p_dye):
            IMP.atom.Mass.setup_particle(p_dye, 100)

        # add dye to atom hier
        p_att = dye.get_source()
        h_att = IMP.atom.Hierarchy(p_att)
        IMP.atom.Hierarchy.setup_particle(p_dye)
        h_dye = IMP.atom.Hierarchy(p_dye)
        h_att.add_child(h_dye)
        # create bond between dye and source
        if not IMP.atom.Bonded.get_is_setup(p_dye):
            IMP.atom.Bonded.setup_particle(p_dye)
        if not IMP.atom.Bonded.get_is_setup(p_att):
            IMP.atom.Bonded.setup_particle(p_att)
        if not IMP.atom.get_bond(IMP.atom.Bonded(p_dye), IMP.atom.Bonded(p_att)):
            IMP.atom.create_bond(
                IMP.atom.Bonded(p_dye), IMP.atom.Bonded(p_att), 1)


class FRETDistanceConverter(object):

    _forster_radius_: float = 52.0

    # current center-to-center distance
    _dist_center_center_: float = 50.0

    # Variables for distance distribution computation
    _distances_: np.ndarray = None  # distance array for distance distribution
    _fret_efficiencies_: np.ndarray = None  # FRET efficiency array corre
    _density_: np.ndarray = None

    # lookup table to cache the computation
    _d_mean_: np.ndarray = None
    _d_mean_fret_: np.ndarray = None
    _e_mean_: np.ndarray = None

    def _update_fret_efficiencies_(self, forster_radius: float):
        RDA = self._distances_
        R0 = forster_radius
        self._fret_efficiencies_ = 1. / (1 + (RDA / R0)**6.0)

    def _update_lookup(self):
        self._d_mean_ = np.empty_like(self._distances_)
        self._d_mean_fret_ = np.empty_like(self._distances_)
        self._e_mean_ = np.empty_like(self._distances_)
        for i, d in enumerate(self._distances_):
            self._update_density(d, self._sigma_)
            self._d_mean_[i] = self.calc_distance_mean()
            self._d_mean_fret_[i] = self.calc_distance_mean_fret()
            self._e_mean_[i] = self.calc_fret_efficiency_mean()

    def _update_density(self, dcc: float, w: float):
        """Update the distance distribution

        @param dcc center-to-center distance
        @param w parameter controlling the distribution width
        """
        d = self._distances_
        n1 = scipy.stats.norm(-dcc, w)
        n2 = scipy.stats.norm(dcc, w)
        p = n1.pdf(d) + n2.pdf(d)
        p /= p.sum()
        self._density_ = p

    @property
    def dist_center_center(self):
        return self._dist_center_center_

    @dist_center_center.setter
    def dist_center_center(self, v):
        self._dist_center_center_ = v

    @property
    def sigma(self):
        return self._sigma_

    @sigma.setter
    def sigma(self, v: float):
        self._sigma_ = v
        self._update_lookup()

    @property
    def forster_radius(self):
        return self._forster_radius_

    @forster_radius.setter
    def forster_radius(self, v: float):
        self._forster_radius_ = v
        self._update_fret_efficiencies_(v)
        self._update_lookup()

    def calc_distance_mean(self):
        p = self._density_
        v = self._distances_
        return p @ v

    def calc_fret_efficiency_mean(self):
        p = self._density_
        v = self._fret_efficiencies_
        return p @ v

    def calc_distance_mean_fret(self):
        E = self.calc_fret_efficiency_mean()
        R0 = self.forster_radius
        return R0 * (1. / E - 1.)**(1. / 6.)

    @property
    def fret_efficiency_mean(self):
        return np.interp(self.dist_center_center, self._distances_, self._e_mean_)

    @property
    def distance_mean(self):
        return np.interp(self.dist_center_center, self._distances_, self._d_mean_)

    @property
    def distance_mean_fret(self):
        return np.interp(self.dist_center_center, self._distances_, self._d_mean_fret_)

    def __call__(self, value: float, distance_type: int):
        self.dist_center_center = value
        if distance_type == IMP.bff.DYE_PAIR_DISTANCE_MEAN:
            return self.distance_mean
        elif distance_type == IMP.bff.DYE_PAIR_DISTANCE_E:
            return self.distance_mean_fret
        elif distance_type == IMP.bff.DYE_PAIR_DISTANCE_MP:
            return value

    def __init__(
            self,
            forster_radius: float = 52.0,
            sigma: float = 6.0,
            distance_range: typing.Tuple[float, float] = (1.0, 100.0),
            n_distances: int = 128
    ):
        self._forster_radius_ = forster_radius
        self._sigma_ = sigma
        self._distances_ = np.linspace(*distance_range, n_distances, dtype=np.float64)
        self._density_ = np.zeros_like(self._distances_)
        self._update_fret_efficiencies_(self._forster_radius_)
        self._update_lookup()


def read_xlink_table(fn: str) -> typing.Dict[int, typing.Dict]:
    """Read a xlink table

    :param fn: filename of the xlink table
    :return:
    """
    # Read the xlink table
    xlinks = {}  # a dict of the xlinks, the keys are used to address the xlinks
    xlink_idx = 0
    with open(fn, 'r') as fp:
        lines = fp.readlines()
        for line in lines[1:]:
            try:
                line = line.strip()
                protein_1, residue_1, protein_2, residue_2 = line.split(',')
                residue_1 = int(residue_1)
                residue_2 = int(residue_2)
                d = {
                    'protein_1': protein_1,
                    'protein_2': protein_2,
                    'residue_1': residue_1,
                    'residue_2': residue_2
                }
                xlinks[xlink_idx] = d
                xlink_idx += 1
            except ValueError:
                pass
    return xlinks

