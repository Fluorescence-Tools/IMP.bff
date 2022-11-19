"""@namespace IMP.bff
Restraints for handling distances between accessible volumes.
"""

from __future__ import print_function

import pathlib

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
