/*
 * The mean-position FRET restraint, the precision of a docked model, and the
 * one restraint in this package that has to stay Python.
 *
 * `AVMeanDistanceRestraint` had two Python copies -- one in
 * `restraints/network.py` without derivatives and one in `restraints/docking.py`
 * with them -- and the wrapper used the one without, so an fps.json-driven
 * docking could only be sampled, never minimised. There is one now, and it
 * computes the gradient; the branch costs nothing when no accumulator is passed.
 *
 * `AVNetworkRestraintWrapper` stays Python because it subclasses
 * `IMP.pmi.restraints.RestraintBase`, which is a Python class in a module
 * `IMP.bff` does not depend on. It is built **lazily**: defining it at import
 * time would make `import IMP.bff` require `IMP.pmi`, and IMP.pmi is not in this
 * module's `required_modules`.
 */

IMP_SWIG_OBJECT(IMP::bff, AVMeanDistanceRestraint, AVMeanDistanceRestraints);
IMP_SWIG_VALUE(IMP::bff, PositionUncertainty, PositionUncertainties);

%feature("shadow") IMP::bff::AVMeanDistanceRestraint::AVMeanDistanceRestraint %{
def __init__(self, m, av1, av2, dist, sigma=6.0, weight=1.0):
    """Score the separation of two mean dye positions against a measurement.

    :param m: the model the two particles belong to.
    :param av1: the donor particle or :class:`IMP.bff.AV` decorator.
    :param av2: the acceptor particle or decorator.
    :param dist: an :class:`IMP.bff.AVPairDistanceMeasurement`.
    :param sigma: per-component width of the separation vector, for the
        Rmp to <RDA> conversion.
    """
    def _index(p):
        p = p.get_particle() if hasattr(p, "get_particle") else p
        return p.get_index() if hasattr(p, "get_index") else p
    this = _IMP_bff.new_AVMeanDistanceRestraint(
        m, _index(av1), _index(av2), dist, float(sigma), float(weight))
    try:
        self.this.append(this)
    except __builtin__.Exception:
        self.this = this
%}

%include "IMP/bff/AVMeanDistanceRestraint.h"
%include "IMP/bff/ModelPrecision.h"

// The scalar RMSF summaries as read-only attributes; the array fields stay
// `get_rmsf()`/`get_mean_coords()` numpy views, reshaped by the caller.
%attribute(IMP::bff::PositionUncertainty, double, rmsf_mean, get_rmsf_mean);
%attribute(IMP::bff::PositionUncertainty, double, rmsf_max, get_rmsf_max);
%attribute(IMP::bff::PositionUncertainty, double, mobile_rmsf_mean,
           get_mobile_rmsf_mean);

%pythoncode %{
def _build_av_network_restraint_wrapper():
    """Define `AVNetworkRestraintWrapper` on first use.

    At import time this would make ``import IMP.bff`` require ``IMP.pmi``, which
    is not one of this module's ``required_modules`` -- so the class is built
    when it is first named and cached on the module by `__getattr__`.
    """
    import pathlib
    import IMP.atom
    import IMP.core
    import IMP.pmi.restraints
    import IMP.pmi.tools

    class AVNetworkRestraintWrapper(IMP.pmi.restraints.RestraintBase):
        """The fps.json AV network restraint, wrapped for IMP.pmi.

        The AVs of the decorated particles are recomputed when the score is
        evaluated, which is the whole cost of this restraint -- pass
        ``mean_position_restraint=True`` to score the separation of the volumes'
        mean positions instead, which is an approximation whose content is that
        the shape of a volume does not change when the structure moves.
        """

        @staticmethod
        def add_used_dyes_to_rb(used_avs):
            """Make each dye a member of its attachment atom's rigid body.

            The coordinates of an AV are the mean of its density, so the AV has
            to be resampled before it is added: the position it is added *at* is
            the one the rigid body will carry it by.
            """
            for key in used_avs:
                dye = used_avs[key]
                dye.resample()
                source = dye.get_source()
                if IMP.core.RigidBodyMember.get_is_setup(source):
                    body = IMP.core.RigidBodyMember(source).get_rigid_body()
                    body.add_member(dye.get_particle())

        def add_xyz_mass_to_avs(self):
            """Give each AV a radius and a mass, so it occupies volume."""
            for key in self.used_avs:
                av = self.used_avs[key]
                r_mean = max(av.get_radii())
                IMP.core.XYZR.setup_particle(av).set_radius(r_mean)
                IMP.atom.Mass.setup_particle(av, 0.1).set_mass(r_mean * 2.0)

        def __init__(self, hier, fps_json_fn, score_set="", weight=1.0,
                     mean_position_restraint=False, sigma_DA=6.0,
                     label="AVNetworkRestraint", occupy_volume=True):
            m = hier.get_model()
            self.mdl = m
            self.hier = hier
            super().__init__(m, label=label, weight=weight)

            self.model_ps = [k.get_particle() for k in IMP.atom.get_leaves(hier)]
            self.mean_position_restraint = mean_position_restraint
            if not pathlib.Path(fps_json_fn).is_file():
                raise FileNotFoundError("{}".format(fps_json_fn))
            self.av_network_restraint = AVNetworkRestraint(
                hier, fps_json_fn, self.name, score_set)
            self.rs = IMP.RestraintSet(m, "AVNetworkRestraint")
            self.used_avs = dict(
                (v.get_name(), v)
                for v in self.av_network_restraint.get_used_avs())
            if not mean_position_restraint:
                self.rs.add_restraint(self.av_network_restraint)
            else:
                self.used_distances = \
                    self.av_network_restraint.get_used_distances()
                self.add_used_dyes_to_rb(self.used_avs)
                for dk in self.used_distances:
                    d_exp = self.used_distances[dk]
                    self.rs.add_restraint(AVMeanDistanceRestraint(
                        m, self.used_avs[d_exp.position_1],
                        self.used_avs[d_exp.position_2], d_exp, sigma=sigma_DA))
            if occupy_volume:
                self.add_xyz_mass_to_avs()
            self.set_weight(weight)

        def evaluate(self):
            """Evaluate the score of the restraint."""
            return self.rs.unprotected_evaluate(None) * self.weight

        def add_to_model(self, add_to_rmf=True):
            IMP.pmi.tools.add_restraint_to_model(
                self.mdl, self.rs, add_to_rmf=add_to_rmf)

    return AVNetworkRestraintWrapper


_LAZY["AVNetworkRestraintWrapper"] = _build_av_network_restraint_wrapper
%}
