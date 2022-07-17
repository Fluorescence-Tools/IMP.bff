"""@namespace IMP.bff
Restraints for handling distances between accessible volumes.
"""

from __future__ import print_function

import typing

import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container
import IMP.bayesianem
import IMP.isd
import IMP.pmi.tools
import IMP.pmi.restraints

import IMP.bff
import IMP.bff.tools


class AVMeanDistanceRestraint(IMP.Restraint):

    def __init__(
            self,
            m: IMP.Model,
            av1: IMP.bff.AV, av2: IMP.bff.AV,
            dist: IMP.bff.AVPairDistanceMeasurement,
            sigma: float,
            weight: float = 1.0
    ):
        """
        input two AV particles, the two equilibrium distances, their amplitudes,
        and their weights (populations)
        sigma is the width of the av1,av2 distance distribution

        :param m:
        :param av1:
        :param av2:
        :param dist:
        :param sigma:
        :param weight:

        Examples
        --------
        import IMP
        import IMP.core
        import IMP.atom
        import IMP.bff

        >>> rmf_file = './test/data/screen/h1_run1_8031.rmf3'
        >>> get_obstacles(rmf_file)

        >>> rmf_file = './test/data/screen_kala/h1_run3_4059_.rmf'
        >>> get_obstacles(rmf_file)

        Returns
        -------
        obstacles: np.ndarray
            An array

        """
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
            dye.resample()
            print("Adding", dk, "with parameters", dye, "to RB at location", IMP.core.XYZ(dye))
            p_att = dye.get_source()
            if IMP.core.RigidBodyMember.get_is_setup(p_att):
                rbm = IMP.core.RigidBodyMember(p_att)
                rb: IMP.core.RigidBody = rbm.get_rigid_body()
                rb.add_member(p_dye)

    def add_xyz_mass_to_avs(self):
        for ak in self.used_avs:
            av: IMP.bff.AV = self.used_avs[ak]
            r_mean = max(av.get_radii())
            IMP.core.XYZR.setup_particle(av, r_mean)
            IMP.atom.Mass.setup_particle(av, 1)
            av_h = IMP.atom.Hierarchy(av)
            self.hier.add_child(av_h)

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
        self.mdl : IMP.Model = m
        self.hier : IMP.atom.Hierarchy = hier
        super(AVNetworkRestraintWrapper, self).__init__(m, label=label, weight=weight)

        self.model_ps = []
        self.model_ps += [k.get_particle() for k in IMP.atom.get_leaves(hier)]

        name = self.name
        self.mean_position_restraint = mean_position_restraint

        self.av_network_restraint = IMP.bff.AVNetworkRestraint(hier, fps_json_fn, name, score_set)
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
