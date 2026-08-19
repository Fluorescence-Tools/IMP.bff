/**
 *  \file IMP/bff/DyeForceField.h
 *  \brief A coarse-grained dye/protein system: sites, bonds, and their terms.
 *
 * The topology and parameters of one simulation system — what
 * #IMP::bff::cgdye builds, writes to mmCIF, reads back, and hands to IMP to
 * become particles and restraints.
 *
 * It was a nested Python dictionary: fifteen keys, read by six modules,
 * holding lists whose elements were positional — a bond was
 * `[site_a, site_b, length, type_id]` and the third slot was a length by
 * convention only. Nothing could state what a field meant, nothing checked a
 * field was present, and every consumer wrote `system.get("bonds", [])`
 * because nothing guaranteed it was. Typed, none of that is expressible.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_DYEFORCEFIELD_H
#define IMPBFF_DYEFORCEFIELD_H

#include <IMP/bff/bff_config.h>

#include <IMP/value_macros.h>
#include <IMP/showable_macros.h>

#include <map>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! One molecule taking part in the system.
struct IMPBFFEXPORT FFComponent {
    std::string mol2_path;
    std::string pdb_path;
    //! `"fixed"` or `"mobile"` — whether the sampler may move it.
    std::string role;
    IMP_SHOWABLE_INLINE(FFComponent, out << "FFComponent(" << role << ")");
};

//! One coarse-grained interaction site.
struct IMPBFFEXPORT FFSite {
    std::string id;
    std::string component;
    std::string atom_name;
    //! Chemical element, which is what picks the site's LJ type.
    //!
    //! Present because the two producers of this structure disagreed about it:
    //! the mmCIF reader set `site_no` and no element, the topology builder set
    //! an element and no `site_no`, and `read_dye_forcefield_cif` guarded with
    //! `s.get("site_no")` because it might not be there. One type, one shape.
    std::string element;
    //! Index within the system, and within its component's own numbering.
    int site_no = 0;
    int site_serial = 0;
    //! Angstrom and Dalton. The defaults are carbon-ish, and are what the
    //! reader falls back to when the file omits them.
    double radius = 1.7;
    double mass = 12.0;
    IMP_SHOWABLE_INLINE(FFSite, out << "FFSite(" << id << ")");
};

//! A harmonic bond between two sites.
struct IMPBFFEXPORT FFBond {
    std::string site_a, site_b;
    //! Equilibrium length, Angstrom.
    double length = 0.0;
    //! Key into DyeForceFieldSystem::get_bond_types().
    std::string type_id;
    IMP_SHOWABLE_INLINE(FFBond, out << "FFBond(" << site_a << "-" << site_b << ")");
};

//! A harmonic angle over three sites.
struct IMPBFFEXPORT FFAngle {
    std::string site_a, site_b, site_c;
    //! Equilibrium angle, radians.
    double theta = 0.0;
    std::string type_id;
    IMP_SHOWABLE_INLINE(FFAngle, out << "FFAngle(" << site_b << ")");
};

//! A torsion or improper over four sites. Its parameters live in its type.
struct IMPBFFEXPORT FFTorsion {
    std::string site_a, site_b, site_c, site_d;
    std::string type_id;
    IMP_SHOWABLE_INLINE(FFTorsion, out << "FFTorsion(" << type_id << ")");
};

//! Periodic torsion parameters.
struct IMPBFFEXPORT FFTorsionType {
    int periodicity = 1;
    //! Phase, radians.
    double phase = 0.0;
    //! Force constant, kcal/mol.
    double k = 0.0;
    IMP_SHOWABLE_INLINE(FFTorsionType, out << "FFTorsionType(n=" << periodicity << ")");
};

//! Lennard-Jones parameters for one site type.
struct IMPBFFEXPORT FFLJType {
    std::string element;
    //! \f$R_{min}/2\f$ in Angstrom, and \f$\epsilon\f$ in kcal/mol — the
    //! CHARMM convention, not \f$\sigma\f$.
    double rmin_half = 0.0;
    double epsilon = 0.0;
    IMP_SHOWABLE_INLINE(FFLJType, out << "FFLJType(" << element << ")");
};

//! A fluorophore in the system, for the flrCIF probe list.
/*!
    One per *mobile* component: in this model a component that the sampler may
    move is a dye, and a fixed one is what it is attached to.
*/
struct IMPBFFEXPORT FFProbe {
    int id = 0;
    std::string name;
    //! `"extrinsic"` for a dye added to the structure, `"intrinsic"` for a
    //! native fluorophore such as a tryptophan.
    std::string origin;
    std::string link_type;
    IMP_SHOWABLE_INLINE(FFProbe, out << "FFProbe(" << name << ")");
};

//! The soft-sphere non-bonded term.
struct IMPBFFEXPORT FFNonbonded {
    bool enabled = true;
    double k = 5.0;
    //! Angstrom.
    double cutoff = 6.0;
    IMP_SHOWABLE_INLINE(FFNonbonded, out << "FFNonbonded(k=" << k << ")");
};

//! How the system is to be sampled. Defaults are the ones the reader supplies.
struct IMPBFFEXPORT FFSampling {
    double temperature_K = 300.0;
    double friction_ps = 10.0;
    double timestep_fs = 0.25;
    int n_steps = 500000;
    int write_every = 1000;
    int minimize_steps = 200;
    IMP_SHOWABLE_INLINE(FFSampling, out << "FFSampling(" << temperature_K << " K)");
};

IMP_VALUES(FFComponent, FFComponents);
IMP_VALUES(FFSite, FFSites);
IMP_VALUES(FFBond, FFBonds);
IMP_VALUES(FFAngle, FFAngles);
IMP_VALUES(FFTorsion, FFTorsions);
IMP_VALUES(FFTorsionType, FFTorsionTypes);
IMP_VALUES(FFLJType, FFLJTypes);
IMP_VALUES(FFProbe, FFProbes);
IMP_VALUES(FFNonbonded, FFNonbondeds);
IMP_VALUES(FFSampling, FFSamplings);

//! A coarse-grained dye/protein system: what to simulate, and with what terms.
/*!
    Groups are named sets of sites, and there are three kinds because they
    answer three different questions: `groups` is the general naming mechanism,
    `rb_groups` says which sites move as one rigid body, and `md_fixed_groups`
    says which are held still during molecular dynamics. They were three keys
    of the same dictionary and nothing said they were different in kind.

    Members are **site id tokens**, as strings. The mmCIF schema spells a
    member two ways -- by index (`n`, or `n_start`/`n_end` for a run) or by id
    (`site_id`, or `site_id_start`/`site_id_end`) -- and the Python dictionary
    this replaces carried whichever the file happened to use, so a consumer
    could receive `4` or `"CX4/S1"` for the same site and had no way to know
    which. One canonical form is the point of having a type at all; resolving
    an index to its id is the reader's job.
*/
class IMPBFFEXPORT DyeForceFieldSystem {
    std::string name_;
    std::map<std::string, FFComponent> components_;
    std::vector<FFSite> sites_;
    std::map<std::string, std::vector<std::string> > groups_;
    std::map<std::string, std::vector<std::string> > rb_groups_;
    std::map<std::string, std::vector<std::string> > md_fixed_groups_;
    std::vector<std::string> fixed_groups_;
    std::map<std::string, double> bond_types_;
    std::map<std::string, double> angle_types_;
    std::map<std::string, FFTorsionType> torsion_types_;
    std::map<std::string, FFTorsionType> improper_types_;
    std::map<std::string, FFLJType> lj_types_;
    std::vector<FFBond> bonds_;
    std::vector<FFAngle> angles_;
    std::vector<FFTorsion> dihedrals_;
    std::vector<FFTorsion> impropers_;
    std::vector<FFProbe> probes_;
    FFNonbonded nonbonded_;
    FFSampling sampling_;

public:
    DyeForceFieldSystem(const std::string& name = "") : name_(name) {}

    std::string get_name() const { return name_; }
    void set_name(const std::string& n) { name_ = n; }

    const std::map<std::string, FFComponent>& get_components() const { return components_; }
    const std::vector<FFSite>& get_sites() const { return sites_; }
    const std::map<std::string, std::vector<std::string> >& get_groups() const { return groups_; }
    const std::map<std::string, std::vector<std::string> >& get_rb_groups() const { return rb_groups_; }
    const std::map<std::string, std::vector<std::string> >& get_md_fixed_groups() const { return md_fixed_groups_; }
    const std::vector<std::string>& get_fixed_groups() const { return fixed_groups_; }
    const std::map<std::string, double>& get_bond_types() const { return bond_types_; }
    const std::map<std::string, double>& get_angle_types() const { return angle_types_; }
    const std::map<std::string, FFTorsionType>& get_torsion_types() const { return torsion_types_; }
    const std::map<std::string, FFTorsionType>& get_improper_types() const { return improper_types_; }
    const std::map<std::string, FFLJType>& get_lj_types() const { return lj_types_; }
    const std::vector<FFBond>& get_bonds() const { return bonds_; }
    const std::vector<FFAngle>& get_angles() const { return angles_; }
    const std::vector<FFTorsion>& get_dihedrals() const { return dihedrals_; }
    const std::vector<FFTorsion>& get_impropers() const { return impropers_; }
    const std::vector<FFProbe>& get_probes() const { return probes_; }
    FFNonbonded get_nonbonded() const { return nonbonded_; }
    FFSampling get_sampling() const { return sampling_; }

    void set_components(const std::map<std::string, FFComponent>& v) { components_ = v; }
    void set_sites(const std::vector<FFSite>& v) { sites_ = v; }
    void set_groups(const std::map<std::string, std::vector<std::string> >& v) { groups_ = v; }
    void set_rb_groups(const std::map<std::string, std::vector<std::string> >& v) { rb_groups_ = v; }
    void set_md_fixed_groups(const std::map<std::string, std::vector<std::string> >& v) { md_fixed_groups_ = v; }
    void set_fixed_groups(const std::vector<std::string>& v) { fixed_groups_ = v; }
    void set_bond_types(const std::map<std::string, double>& v) { bond_types_ = v; }
    void set_angle_types(const std::map<std::string, double>& v) { angle_types_ = v; }
    void set_torsion_types(const std::map<std::string, FFTorsionType>& v) { torsion_types_ = v; }
    void set_improper_types(const std::map<std::string, FFTorsionType>& v) { improper_types_ = v; }
    void set_lj_types(const std::map<std::string, FFLJType>& v) { lj_types_ = v; }
    void set_bonds(const std::vector<FFBond>& v) { bonds_ = v; }
    void set_angles(const std::vector<FFAngle>& v) { angles_ = v; }
    void set_dihedrals(const std::vector<FFTorsion>& v) { dihedrals_ = v; }
    void set_impropers(const std::vector<FFTorsion>& v) { impropers_ = v; }
    void set_probes(const std::vector<FFProbe>& v) { probes_ = v; }
    void set_nonbonded(const FFNonbonded& v) { nonbonded_ = v; }
    void set_sampling(const FFSampling& v) { sampling_ = v; }

    //! Sites belonging to a component, by name.
    std::vector<int> get_component_sites(const std::string& component) const;

    IMP_SHOWABLE_INLINE(DyeForceFieldSystem,
                        out << "DyeForceFieldSystem(\"" << name_ << "\", "
                            << sites_.size() << " sites, " << bonds_.size()
                            << " bonds, " << components_.size() << " components)");
};
IMP_VALUES(DyeForceFieldSystem, DyeForceFieldSystems);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_DYEFORCEFIELD_H
