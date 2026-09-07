/**
 * \file IMP/bff/LinkerGeometry.h
 * \brief Applying a torsion/angle configuration to a linker.
 *
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_LINKERGEOMETRY_H
#define IMPBFF_LINKERGEOMETRY_H

#include <IMP/bff/bff_config.h>

#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! A linker's rotatable degrees of freedom, and how to apply a setting of them.
/** The topology -- which atoms move when a torsion or bond angle turns -- is
    fixed for the whole of a sampling run, so it is set once here and only the
    configuration changes per call.

    Rebuilding an `IMP.core.XYZ` decorator per moving atom per rotation would
    cost a boundary crossing each time, which for a probe linker is thousands
    per configuration.

    Rotations are applied **sequentially, in the order given**, each about the
    axis as it stands after the previous one. That is not the same as applying
    them to the reference geometry independently, and it is what the sampler
    means: a torsion further down the chain turns what the ones before it
    already moved. */
class IMPBFFEXPORT LinkerGeometry {
    std::vector<double> base_;                       //!< (n_atoms, 3), the reference
    std::vector<int> torsion_fixed_, torsion_moving_;
    std::vector<std::vector<int> > torsion_sets_;
    std::vector<int> angle_b_, angle_c_, angle_a_;
    std::vector<std::vector<int> > angle_sets_;

public:
    LinkerGeometry() {}

    //! Reference coordinates, `(n_atoms, 3)` flattened.
    void set_coordinates(const std::vector<double>& xyz) { base_ = xyz; }
    const std::vector<double>& get_coordinates() const { return base_; }

    //! One rotatable bond: the axis runs `fixed` -> `moving`, and `moves` turn.
    void add_torsion(int fixed, int moving, const std::vector<int>& moves);

    //! One rotatable angle at `b`, in the plane of `a`-`b`-`c`; `moves` turn.
    void add_angle(int b, int c, int a, const std::vector<int>& moves);

    unsigned int get_number_of_torsions() const { return (unsigned int) torsion_fixed_.size(); }
    unsigned int get_number_of_angles() const { return (unsigned int) angle_b_.size(); }

    //! Which atoms each degree of freedom is defined by, and what it moves.
    /*! Rows into the coordinate array, in the order the degrees of freedom
        were added -- which is the order a configuration vector is in. A
        caller that wants to check this class against its own rotation needs
        to know what it claims to turn. */
    const std::vector<int>& get_torsion_fixed() const { return torsion_fixed_; }
    const std::vector<int>& get_torsion_moving() const { return torsion_moving_; }
    std::vector<int> get_torsion_set(unsigned int i) const {
        return i < torsion_sets_.size() ? torsion_sets_[i] : std::vector<int>();
    }
    const std::vector<int>& get_angle_b() const { return angle_b_; }
    const std::vector<int>& get_angle_c() const { return angle_c_; }
    const std::vector<int>& get_angle_a() const { return angle_a_; }
    std::vector<int> get_angle_set(unsigned int i) const {
        return i < angle_sets_.size() ? angle_sets_[i] : std::vector<int>();
    }

    //! The reference geometry with `config` applied, `(n_atoms, 3)` flattened.
    /** `config` is the torsion settings followed by the angle settings, in
        radians. An axis shorter than 1e-8 is skipped rather than normalised:
        two atoms at the same point define no rotation. */
    std::vector<double> apply(const std::vector<double>& config) const;
};

IMPBFF_END_NAMESPACE

#endif //IMPBFF_LINKERGEOMETRY_H
