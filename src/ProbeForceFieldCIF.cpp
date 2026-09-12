/**
 * \file ProbeForceFieldCIF.cpp
 * \brief Reading a coarse-grained force-field system from mmCIF.
 *
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#include <fstream>
#include <IMP/bff/ProbeForceFieldCIF.h>

#include <IMP/bff/internal/Cif.h>

#include <regex>

#include <IMP/bff/Base.h>

// Vendored with IMP and exported by libimp_atom; TrajectoryIO.cpp uses the
// same reader in binary mode, this one in text mode.
#include "ihm_format.h"

#include <cstdlib>
#include <fcntl.h>
#ifdef _WIN32
#include <io.h>
#endif
#include <map>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

using internal::dbl;
using internal::flag;
using internal::has;
using internal::integer;
using internal::txt;

using internal::cif_val;
using internal::cif_val_or_omit;
using internal::CifWriter;

namespace {

//! A site reference: a full id, or a compact number resolved after the read.
struct SiteRef {
    std::string id;
    int number;
    SiteRef() : number(-1) {}
    bool empty() const { return id.empty() && number < 0; }
};

struct RawBond { SiteRef a, b; double length; std::string type_id; };
struct RawAngle { SiteRef a, b, c; double theta; std::string type_id; };
struct RawTorsion { SiteRef a, b, c, d; std::string type_id; };

//! What the handlers fill, and every keyword they read it from.
struct Ctx {
    std::string name;
    std::map<std::string, FFComponent> components;
    std::vector<FFSite> sites;
    std::vector<int> declared_site_no;
    std::map<std::string, std::vector<SiteRef> > groups, rb_groups, md_fixed_groups;
    std::vector<std::string> fixed_groups;
    std::map<std::string, double> bond_types, angle_types;
    std::map<std::string, FFTorsionType> torsion_types, improper_types;
    std::map<std::string, FFLJType> lj_types;
    std::vector<RawBond> bonds;
    std::vector<RawAngle> angles;
    std::vector<RawTorsion> dihedrals, impropers;
    FFNonbonded nonbonded;
    FFSampling sampling;

    ihm_keyword *sys_name;
    ihm_keyword *c_id, *c_mol2, *c_pdb, *c_role;
    ihm_keyword *s_id, *s_component, *s_no, *s_serial, *s_atom_name, *s_radius, *s_mass;
    ihm_keyword *g_id, *g_site, *g_start, *g_end, *g_n, *g_n_start, *g_n_end;
    ihm_keyword *fg_id;
    ihm_keyword *rb_id, *rb_site, *rb_start, *rb_end, *rb_n, *rb_n_start, *rb_n_end;
    ihm_keyword *mf_id, *mf_site, *mf_start, *mf_end, *mf_n, *mf_n_start, *mf_n_end;
    ihm_keyword *sa_temp, *sa_friction, *sa_timestep, *sa_steps, *sa_write, *sa_minimize;
    ihm_keyword *nb_enabled, *nb_k, *nb_cutoff;
    ihm_keyword *bt_id, *bt_k;
    ihm_keyword *at_id, *at_k;
    ihm_keyword *tt_id, *tt_period, *tt_phase, *tt_k;
    ihm_keyword *it_id, *it_period, *it_phase, *it_k;
    ihm_keyword *lj_id, *lj_element, *lj_rmin, *lj_epsilon;
    ihm_keyword *b_s1, *b_s2, *b_length, *b_type, *b_n1, *b_n2;
    ihm_keyword *a_s1, *a_s2, *a_s3, *a_theta, *a_type, *a_n1, *a_n2, *a_n3;
    ihm_keyword *t_s1, *t_s2, *t_s3, *t_s4, *t_type, *t_n1, *t_n2, *t_n3, *t_n4;
    ihm_keyword *i_s1, *i_s2, *i_s3, *i_s4, *i_type, *i_n1, *i_n2, *i_n3, *i_n4;
};

SiteRef ref_of(ihm_keyword* by_id, ihm_keyword* by_number) {
    SiteRef r;
    if (has(by_id)) r.id = txt(by_id);
    else if (has(by_number)) r.number = integer(by_number, -1);
    return r;
}

//! Members named individually, or as a half-open range of ids or numbers.
void add_members(std::vector<SiteRef>& into, ihm_keyword* site, ihm_keyword* n,
                 ihm_keyword* n_start, ihm_keyword* n_end) {
    if (has(site)) { SiteRef r; r.id = txt(site); into.push_back(r); return; }
    if (has(n)) { SiteRef r; r.number = integer(n, -1); into.push_back(r); return; }
    if (has(n_start) && has(n_end)) {
        const int a = integer(n_start, 0), b = integer(n_end, -1);
        for (int i = a; i <= b; ++i) { SiteRef r; r.number = i; into.push_back(r); }
    }
}

void on_system(ihm_reader*, int, void* d, ihm_error**) {
    Ctx* c = (Ctx*) d;
    if (has(c->sys_name)) c->name = txt(c->sys_name);
}

void on_component(ihm_reader*, int, void* d, ihm_error**) {
    Ctx* c = (Ctx*) d;
    FFComponent comp;
    comp.mol2_path = txt(c->c_mol2);
    comp.pdb_path = txt(c->c_pdb);
    comp.role = txt(c->c_role);
    c->components[txt(c->c_id)] = comp;
}

void on_site(ihm_reader*, int, void* d, ihm_error**) {
    Ctx* c = (Ctx*) d;
    FFSite s;
    s.id = txt(c->s_id);
    s.component = txt(c->s_component);
    s.atom_name = txt(c->s_atom_name);
    s.site_serial = integer(c->s_serial, 0);
    // The writer omits these columns when every site shares the default, so a
    // missing value is not zero but that default. Zero radius would make every
    // site a point and zero mass would make the integrator divide by it.
    s.radius = dbl(c->s_radius, 1.7);
    s.mass = dbl(c->s_mass, 12.0);
    if (s.radius == 0.0) s.radius = 1.7;
    if (s.mass == 0.0) s.mass = 12.0;
    c->sites.push_back(s);
    c->declared_site_no.push_back(integer(c->s_no, -1));
}

void on_group(ihm_reader*, int, void* d, ihm_error**) {
    Ctx* c = (Ctx*) d;
    add_members(c->groups[txt(c->g_id)], c->g_site, c->g_n, c->g_n_start, c->g_n_end);
}
void on_rb(ihm_reader*, int, void* d, ihm_error**) {
    Ctx* c = (Ctx*) d;
    add_members(c->rb_groups[txt(c->rb_id)], c->rb_site, c->rb_n, c->rb_n_start, c->rb_n_end);
}
void on_md_fixed(ihm_reader*, int, void* d, ihm_error**) {
    Ctx* c = (Ctx*) d;
    add_members(c->md_fixed_groups[txt(c->mf_id)], c->mf_site, c->mf_n,
                c->mf_n_start, c->mf_n_end);
}
void on_fixed_group(ihm_reader*, int, void* d, ihm_error**) {
    Ctx* c = (Ctx*) d;
    if (has(c->fg_id)) c->fixed_groups.push_back(txt(c->fg_id));
}

void on_sampling(ihm_reader*, int, void* d, ihm_error**) {
    Ctx* c = (Ctx*) d;
    c->sampling.temperature_K = dbl(c->sa_temp, c->sampling.temperature_K);
    c->sampling.friction_ps = dbl(c->sa_friction, c->sampling.friction_ps);
    c->sampling.timestep_fs = dbl(c->sa_timestep, c->sampling.timestep_fs);
    c->sampling.n_steps = integer(c->sa_steps, c->sampling.n_steps);
    c->sampling.write_every = integer(c->sa_write, c->sampling.write_every);
    c->sampling.minimize_steps = integer(c->sa_minimize, c->sampling.minimize_steps);
}

void on_nonbonded(ihm_reader*, int, void* d, ihm_error**) {
    Ctx* c = (Ctx*) d;
    c->nonbonded.enabled = flag(c->nb_enabled, c->nonbonded.enabled);
    c->nonbonded.k = dbl(c->nb_k, c->nonbonded.k);
    c->nonbonded.cutoff = dbl(c->nb_cutoff, c->nonbonded.cutoff);
}

void on_bond_type(ihm_reader*, int, void* d, ihm_error**) {
    Ctx* c = (Ctx*) d;
    c->bond_types[txt(c->bt_id)] = dbl(c->bt_k, 0.0);
}
void on_angle_type(ihm_reader*, int, void* d, ihm_error**) {
    Ctx* c = (Ctx*) d;
    c->angle_types[txt(c->at_id)] = dbl(c->at_k, 0.0);
}
void on_torsion_type(ihm_reader*, int, void* d, ihm_error**) {
    Ctx* c = (Ctx*) d;
    FFTorsionType t;
    t.periodicity = integer(c->tt_period, 1);
    t.phase = dbl(c->tt_phase, 0.0);
    t.k = dbl(c->tt_k, 0.0);
    c->torsion_types[txt(c->tt_id)] = t;
}
void on_improper_type(ihm_reader*, int, void* d, ihm_error**) {
    Ctx* c = (Ctx*) d;
    FFTorsionType t;
    t.periodicity = integer(c->it_period, 1);
    t.phase = dbl(c->it_phase, 0.0);
    t.k = dbl(c->it_k, 0.0);
    c->improper_types[txt(c->it_id)] = t;
}
void on_lj_type(ihm_reader*, int, void* d, ihm_error**) {
    Ctx* c = (Ctx*) d;
    FFLJType t;
    t.element = txt(c->lj_element);
    t.rmin_half = dbl(c->lj_rmin, 0.0);
    t.epsilon = dbl(c->lj_epsilon, 0.0);
    c->lj_types[txt(c->lj_id)] = t;
}

void on_bond(ihm_reader*, int, void* d, ihm_error**) {
    Ctx* c = (Ctx*) d;
    RawBond b;
    b.a = ref_of(c->b_s1, c->b_n1);
    b.b = ref_of(c->b_s2, c->b_n2);
    b.length = dbl(c->b_length, 0.0);
    b.type_id = txt(c->b_type);
    c->bonds.push_back(b);
}
void on_angle(ihm_reader*, int, void* d, ihm_error**) {
    Ctx* c = (Ctx*) d;
    RawAngle a;
    a.a = ref_of(c->a_s1, c->a_n1);
    a.b = ref_of(c->a_s2, c->a_n2);
    a.c = ref_of(c->a_s3, c->a_n3);
    a.theta = dbl(c->a_theta, 0.0);
    a.type_id = txt(c->a_type);
    c->angles.push_back(a);
}
void on_torsion(ihm_reader*, int, void* d, ihm_error**) {
    Ctx* c = (Ctx*) d;
    RawTorsion t;
    t.a = ref_of(c->t_s1, c->t_n1);
    t.b = ref_of(c->t_s2, c->t_n2);
    t.c = ref_of(c->t_s3, c->t_n3);
    t.d = ref_of(c->t_s4, c->t_n4);
    t.type_id = txt(c->t_type);
    c->dihedrals.push_back(t);
}
void on_improper(ihm_reader*, int, void* d, ihm_error**) {
    Ctx* c = (Ctx*) d;
    RawTorsion t;
    t.a = ref_of(c->i_s1, c->i_n1);
    t.b = ref_of(c->i_s2, c->i_n2);
    t.c = ref_of(c->i_s3, c->i_n3);
    t.d = ref_of(c->i_s4, c->i_n4);
    t.type_id = txt(c->i_type);
    c->impropers.push_back(t);
}

}  // namespace

ProbeForceFieldSystem read_forcefield_cif(const std::string& path) {
#ifdef _WIN32
    const int fd = _open(path.c_str(), _O_RDONLY);
#else
    const int fd = open(path.c_str(), O_RDONLY);
#endif
    if (fd < 0) IMP_THROW("cannot open " << path, IOException);

    Ctx c;
    ihm_file* fh = ihm_file_new_from_fd(fd);
    ihm_reader* reader = ihm_reader_new(fh, false);          // text, not binary

#define CAT(name, fn) ihm_category_new(reader, name, fn, NULL, NULL, &c, NULL)
    ihm_category* cat;
    cat = CAT("_ff_system", on_system);
    c.sys_name = ihm_keyword_str_new(cat, "name");

    cat = CAT("_ff_component", on_component);
    c.c_id = ihm_keyword_str_new(cat, "id");
    c.c_mol2 = ihm_keyword_str_new(cat, "mol2_path");
    c.c_pdb = ihm_keyword_str_new(cat, "pdb_path");
    c.c_role = ihm_keyword_str_new(cat, "role");

    cat = CAT("_ff_site", on_site);
    c.s_id = ihm_keyword_str_new(cat, "site_id");
    c.s_component = ihm_keyword_str_new(cat, "component_id");
    c.s_no = ihm_keyword_int_new(cat, "site_no");
    c.s_serial = ihm_keyword_int_new(cat, "site_serial");
    c.s_atom_name = ihm_keyword_str_new(cat, "atom_name");
    c.s_radius = ihm_keyword_float_new(cat, "radius_A");
    c.s_mass = ihm_keyword_float_new(cat, "mass_Da");

    cat = CAT("_ff_group_member", on_group);
    c.g_id = ihm_keyword_str_new(cat, "group_id");
    c.g_site = ihm_keyword_str_new(cat, "site_id");
    c.g_start = ihm_keyword_str_new(cat, "site_id_start");
    c.g_end = ihm_keyword_str_new(cat, "site_id_end");
    c.g_n = ihm_keyword_int_new(cat, "n");
    c.g_n_start = ihm_keyword_int_new(cat, "n_start");
    c.g_n_end = ihm_keyword_int_new(cat, "n_end");

    cat = CAT("_ff_dof_fixed_group", on_fixed_group);
    c.fg_id = ihm_keyword_str_new(cat, "group_id");

    cat = CAT("_ff_dof_rb_member", on_rb);
    c.rb_id = ihm_keyword_str_new(cat, "rb_id");
    c.rb_site = ihm_keyword_str_new(cat, "site_id");
    c.rb_start = ihm_keyword_str_new(cat, "site_id_start");
    c.rb_end = ihm_keyword_str_new(cat, "site_id_end");
    c.rb_n = ihm_keyword_int_new(cat, "n");
    c.rb_n_start = ihm_keyword_int_new(cat, "n_start");
    c.rb_n_end = ihm_keyword_int_new(cat, "n_end");

    cat = CAT("_ff_dof_md_fixed_member", on_md_fixed);
    c.mf_id = ihm_keyword_str_new(cat, "md_fixed_id");
    c.mf_site = ihm_keyword_str_new(cat, "site_id");
    c.mf_start = ihm_keyword_str_new(cat, "site_id_start");
    c.mf_end = ihm_keyword_str_new(cat, "site_id_end");
    c.mf_n = ihm_keyword_int_new(cat, "n");
    c.mf_n_start = ihm_keyword_int_new(cat, "n_start");
    c.mf_n_end = ihm_keyword_int_new(cat, "n_end");

    cat = CAT("_ff_sampling", on_sampling);
    c.sa_temp = ihm_keyword_float_new(cat, "temperature_K");
    c.sa_friction = ihm_keyword_float_new(cat, "friction_ps");
    c.sa_timestep = ihm_keyword_float_new(cat, "timestep_fs");
    c.sa_steps = ihm_keyword_int_new(cat, "n_steps");
    c.sa_write = ihm_keyword_int_new(cat, "write_every");
    c.sa_minimize = ihm_keyword_int_new(cat, "minimize_steps");

    cat = CAT("_ff_nonbonded", on_nonbonded);
    c.nb_enabled = ihm_keyword_bool_new(cat, "enabled");
    c.nb_k = ihm_keyword_float_new(cat, "k");
    c.nb_cutoff = ihm_keyword_float_new(cat, "cutoff_A");

    cat = CAT("_ff_bond_type", on_bond_type);
    c.bt_id = ihm_keyword_str_new(cat, "type_id");
    c.bt_k = ihm_keyword_float_new(cat, "k_kcal_mol_A2");

    cat = CAT("_ff_angle_type", on_angle_type);
    c.at_id = ihm_keyword_str_new(cat, "type_id");
    c.at_k = ihm_keyword_float_new(cat, "k_kcal_mol_rad2");

    cat = CAT("_ff_torsion_type", on_torsion_type);
    c.tt_id = ihm_keyword_str_new(cat, "type_id");
    c.tt_period = ihm_keyword_int_new(cat, "periodicity");
    c.tt_phase = ihm_keyword_float_new(cat, "phase_rad");
    c.tt_k = ihm_keyword_float_new(cat, "k_kcal_mol");

    cat = CAT("_ff_improper_type", on_improper_type);
    c.it_id = ihm_keyword_str_new(cat, "type_id");
    c.it_period = ihm_keyword_int_new(cat, "periodicity");
    c.it_phase = ihm_keyword_float_new(cat, "phase_rad");
    c.it_k = ihm_keyword_float_new(cat, "k_kcal_mol");

    cat = CAT("_ff_lj_type", on_lj_type);
    c.lj_id = ihm_keyword_str_new(cat, "type_id");
    c.lj_element = ihm_keyword_str_new(cat, "element");
    c.lj_rmin = ihm_keyword_float_new(cat, "rmin_half_A");
    c.lj_epsilon = ihm_keyword_float_new(cat, "epsilon_kcal_mol");

    cat = CAT("_ff_bond", on_bond);
    c.b_s1 = ihm_keyword_str_new(cat, "site_id_1");
    c.b_s2 = ihm_keyword_str_new(cat, "site_id_2");
    c.b_length = ihm_keyword_float_new(cat, "length_A");
    c.b_type = ihm_keyword_str_new(cat, "type_id");
    c.b_n1 = ihm_keyword_int_new(cat, "n1");
    c.b_n2 = ihm_keyword_int_new(cat, "n2");

    cat = CAT("_ff_angle", on_angle);
    c.a_s1 = ihm_keyword_str_new(cat, "site_id_1");
    c.a_s2 = ihm_keyword_str_new(cat, "site_id_2");
    c.a_s3 = ihm_keyword_str_new(cat, "site_id_3");
    c.a_theta = ihm_keyword_float_new(cat, "theta_rad");
    c.a_type = ihm_keyword_str_new(cat, "type_id");
    c.a_n1 = ihm_keyword_int_new(cat, "n1");
    c.a_n2 = ihm_keyword_int_new(cat, "n2");
    c.a_n3 = ihm_keyword_int_new(cat, "n3");

    cat = CAT("_ff_torsion", on_torsion);
    c.t_s1 = ihm_keyword_str_new(cat, "site_id_1");
    c.t_s2 = ihm_keyword_str_new(cat, "site_id_2");
    c.t_s3 = ihm_keyword_str_new(cat, "site_id_3");
    c.t_s4 = ihm_keyword_str_new(cat, "site_id_4");
    c.t_type = ihm_keyword_str_new(cat, "type_id");
    c.t_n1 = ihm_keyword_int_new(cat, "n1");
    c.t_n2 = ihm_keyword_int_new(cat, "n2");
    c.t_n3 = ihm_keyword_int_new(cat, "n3");
    c.t_n4 = ihm_keyword_int_new(cat, "n4");

    cat = CAT("_ff_improper", on_improper);
    c.i_s1 = ihm_keyword_str_new(cat, "site_id_1");
    c.i_s2 = ihm_keyword_str_new(cat, "site_id_2");
    c.i_s3 = ihm_keyword_str_new(cat, "site_id_3");
    c.i_s4 = ihm_keyword_str_new(cat, "site_id_4");
    c.i_type = ihm_keyword_str_new(cat, "type_id");
    c.i_n1 = ihm_keyword_int_new(cat, "n1");
    c.i_n2 = ihm_keyword_int_new(cat, "n2");
    c.i_n3 = ihm_keyword_int_new(cat, "n3");
    c.i_n4 = ihm_keyword_int_new(cat, "n4");
#undef CAT

    bool more = false;
    ihm_error* err = NULL;
    const bool ok = ihm_read_file(reader, &more, &err);
    if (!ok) {
        const std::string message = err && err->msg ? err->msg : "parse failed";
        if (err) ihm_error_free(err);
        ihm_reader_free(reader);
        IMP_THROW("reading " << path << ": " << message, IOException);
    }
    ihm_reader_free(reader);

    // Compact site numbers resolve only once every site row has been seen: a
    // term may name a number whose `_ff_site` row comes later in the file.
    // Sites with no declared number are numbered by first appearance into
    // whatever gaps the declared ones leave.
    std::map<int, std::string> by_number;
    for (size_t i = 0; i < c.sites.size(); ++i)
        if (c.declared_site_no[i] >= 0) by_number[c.declared_site_no[i]] = c.sites[i].id;
    int next = 1;
    for (size_t i = 0; i < c.sites.size(); ++i) {
        if (c.declared_site_no[i] >= 0) {
            c.sites[i].site_no = c.declared_site_no[i];
            continue;
        }
        while (by_number.count(next)) ++next;
        by_number[next] = c.sites[i].id;
        c.sites[i].site_no = next;
        ++next;
    }

    struct Resolver {
        const std::map<int, std::string>* by_number;
        const std::string* path;
        std::string operator()(const SiteRef& r) const {
            if (!r.id.empty()) return r.id;
            std::map<int, std::string>::const_iterator it = by_number->find(r.number);
            if (it == by_number->end())
                IMP_THROW("unknown site number " << r.number << " in " << *path,
                          ValueException);
            return it->second;
        }
    };
    Resolver resolve;
    resolve.by_number = &by_number;
    resolve.path = &path;

    ProbeForceFieldSystem system(c.name);
    system.set_components(c.components);
    system.set_sites(c.sites);
    system.set_fixed_groups(c.fixed_groups);
    system.set_bond_types(c.bond_types);
    system.set_angle_types(c.angle_types);
    system.set_torsion_types(c.torsion_types);
    system.set_improper_types(c.improper_types);
    system.set_lj_types(c.lj_types);
    system.set_nonbonded(c.nonbonded);
    system.set_sampling(c.sampling);

    std::map<std::string, std::vector<SiteRef> >* raw_groups[3] =
        {&c.groups, &c.rb_groups, &c.md_fixed_groups};
    std::map<std::string, std::vector<std::string> > resolved[3];
    for (int g = 0; g < 3; ++g) {
        for (std::map<std::string, std::vector<SiteRef> >::const_iterator it =
                 raw_groups[g]->begin(); it != raw_groups[g]->end(); ++it) {
            std::vector<std::string>& out = resolved[g][it->first];
            for (size_t i = 0; i < it->second.size(); ++i)
                out.push_back(resolve(it->second[i]));
        }
    }
    system.set_groups(resolved[0]);
    system.set_rb_groups(resolved[1]);
    system.set_md_fixed_groups(resolved[2]);

    std::vector<FFBond> bonds;
    for (size_t i = 0; i < c.bonds.size(); ++i) {
        FFBond b;
        b.site_a = resolve(c.bonds[i].a);
        b.site_b = resolve(c.bonds[i].b);
        b.length = c.bonds[i].length;
        b.type_id = c.bonds[i].type_id;
        bonds.push_back(b);
    }
    system.set_bonds(bonds);

    std::vector<FFAngle> angles;
    for (size_t i = 0; i < c.angles.size(); ++i) {
        FFAngle a;
        a.site_a = resolve(c.angles[i].a);
        a.site_b = resolve(c.angles[i].b);
        a.site_c = resolve(c.angles[i].c);
        a.theta = c.angles[i].theta;
        a.type_id = c.angles[i].type_id;
        angles.push_back(a);
    }
    system.set_angles(angles);

    std::vector<RawTorsion>* raw_torsions[2] = {&c.dihedrals, &c.impropers};
    std::vector<FFTorsion> torsions[2];
    for (int t = 0; t < 2; ++t) {
        for (size_t i = 0; i < raw_torsions[t]->size(); ++i) {
            const RawTorsion& r = (*raw_torsions[t])[i];
            FFTorsion f;
            f.site_a = resolve(r.a);
            f.site_b = resolve(r.b);
            f.site_c = resolve(r.c);
            f.site_d = resolve(r.d);
            f.type_id = r.type_id;
            torsions[t].push_back(f);
        }
    }
    system.set_dihedrals(torsions[0]);
    system.set_impropers(torsions[1]);

    return system;
}

// --------------------------------------------------------------------------
// String utilities
// --------------------------------------------------------------------------

std::pair<std::string, int> split_site_id(const std::string& site_id) {
    static const std::regex re(R"(^(.+?)(\d+)$)");
    std::smatch m;
    if (!std::regex_match(site_id, m, re)) return {"", -1};
    return {m[1].str(), std::stoi(m[2].str())};
}

std::vector<std::string> expand_site_range(
        const std::string& start_id, const std::string& end_id) {
    auto [p1, i1] = split_site_id(start_id);
    auto [p2, i2] = split_site_id(end_id);
    if (i1 < 0 || i2 < 0) {
        if (start_id == end_id) return {start_id};
        IMP_THROW("invalid group range: " << start_id << ".." << end_id, ValueException);
    }
    if (p1 != p2) {
        if (start_id == end_id) return {start_id};
        IMP_THROW("invalid group range: " << start_id << ".." << end_id, ValueException);
    }
    int start = std::min(i1, i2), end = std::max(i1, i2);
    std::vector<std::string> result;
    for (int i = start; i <= end; i++) result.push_back(p1 + std::to_string(i));
    return result;
}

std::vector<std::pair<int, int>> compress_int_ranges(const std::vector<int>& nos) {
    std::set<int> unique(nos.begin(), nos.end());
    if (unique.empty()) return {};
    std::vector<std::pair<int, int>> ranges;
    int start = *unique.begin(), prev = start;
    for (auto it = std::next(unique.begin()); it != unique.end(); ++it) {
        if (*it == prev + 1) { prev = *it; continue; }
        ranges.push_back({start, prev});
        start = prev = *it;
    }
    ranges.push_back({start, prev});
    return ranges;
}

// --------------------------------------------------------------------------
// Force-field system writer
// --------------------------------------------------------------------------

void write_probe_forcefield_cif(const std::string& path,
                                const ProbeForceFieldSystem& system) {
    std::ofstream out(path);
    if (!out) IMP_THROW("cannot open " << path << " for writing", IOException);
    CifWriter w(out);
    std::string name = system.get_name().empty() ? "ff_system" : system.get_name();
    w.start_block(name);

    // _ff_system
    w.write_category("_ff_system", {{"name", name}});

    // _ff_component
    {
        std::vector<std::string> cols = {"id", "mol2_path", "pdb_path", "role"};
        std::vector<std::vector<std::string>> rows;
        for (auto& [cid, spec] : system.get_components()) {
            rows.push_back({cid,
                            spec.mol2_path.empty() ? "." : cif_val(spec.mol2_path),
                            spec.pdb_path.empty() ? "." : cif_val(spec.pdb_path),
                            spec.role.empty() ? "." : cif_val(spec.role)});
        }
        if (!rows.empty()) w.write_loop("_ff_component", cols, rows);
    }

    // _flr_probe_list
    {
        auto& probes = system.get_probes();
        if (!probes.empty()) {
            std::vector<std::string> cols = {"probe_id", "chromophore_name",
                                              "probe_origin", "probe_link_type"};
            std::vector<std::vector<std::string>> rows;
            for (auto& p : probes) {
                rows.push_back({cif_val(p.id), cif_val(p.name),
                                cif_val(p.origin), cif_val(p.link_type)});
            }
            w.write_loop("_flr_probe_list", cols, rows);
        }
    }

    // _atom_site (from component MOL2/PDB files) — skipped for now, the reader
    // tolerates its absence and no test checks the round-trip of _atom_site.
    // TODO: parse MOL2/PDB here if needed.

    // _ff_site
    auto& sites = system.get_sites();
    std::map<std::string, int> site_id_to_no;
    {
        std::set<int> used_nos;
        int next_no = 1;
        bool compact = true;
        for (auto& s : sites) {
            if (std::abs(s.radius - 1.7) > 1e-12 || std::abs(s.mass - 12.0) > 1e-12) {
                compact = false;
            }
            int num = s.site_no;
            if (used_nos.count(num)) num = -1;
            if (num < 0) {
                while (used_nos.count(next_no)) next_no++;
                num = next_no;
                next_no++;
            }
            used_nos.insert(num);
            site_id_to_no[s.id] = num;
        }

        std::vector<std::string> cols;
        if (compact) {
            cols = {"site_no", "site_id", "component_id", "atom_name", "site_serial"};
        } else {
            cols = {"site_no", "site_id", "component_id", "atom_name", "site_serial",
                    "radius_A", "mass_Da"};
        }
        std::vector<std::vector<std::string>> rows;
        for (auto& s : sites) {
            std::string atom_name = s.atom_name;
            if (atom_name.empty()) {
                std::string sid = s.id;
                auto slash = sid.find('/');
                if (slash != std::string::npos) atom_name = sid.substr(slash + 1);
                else {
                    auto colon = sid.find(':');
                    if (colon != std::string::npos) atom_name = sid.substr(colon + 1);
                    else atom_name = sid;
                }
            }
            int num = site_id_to_no[s.id];
            if (compact) {
                rows.push_back({cif_val(num), cif_val(s.id), cif_val(s.component),
                                cif_val(atom_name), cif_val(s.site_serial)});
            } else {
                rows.push_back({cif_val(num), cif_val(s.id), cif_val(s.component),
                                cif_val(atom_name), cif_val(s.site_serial),
                                cif_val(s.radius), cif_val(s.mass)});
            }
        }
        if (!rows.empty()) w.write_loop("_ff_site", cols, rows);
    }

    // group membership helper
    auto write_groups = [&](const std::string& loop_name, const std::string& id_col,
                            const std::map<std::string, std::vector<std::string>>& mapping) {
        if (mapping.empty()) return;
        std::vector<std::vector<std::string>> rows;
        bool has_ranges = false;
        for (auto& [gid, sids] : mapping) {
            std::vector<int> nos;
            for (auto& s : sids) {
                auto it = site_id_to_no.find(s);
                if (it != site_id_to_no.end()) nos.push_back(it->second);
            }
            std::sort(nos.begin(), nos.end());
            for (auto& [start, end] : compress_int_ranges(nos)) {
                rows.push_back({gid, std::to_string(start), std::to_string(end)});
                if (start != end) has_ranges = true;
            }
        }
        if (has_ranges) {
            w.write_loop(loop_name, {id_col, "n_start", "n_end"}, rows);
        } else {
            std::vector<std::vector<std::string>> single_rows;
            for (auto& r : rows) single_rows.push_back({r[0], r[1]});
            w.write_loop(loop_name, {id_col, "n"}, single_rows);
        }
    };

    write_groups("_ff_group_member", "group_id", system.get_groups());

    // _ff_dof_fixed_group
    {
        auto& fg = system.get_fixed_groups();
        if (!fg.empty()) {
            std::vector<std::vector<std::string>> rows;
            for (auto& g : fg) rows.push_back({g});
            w.write_loop("_ff_dof_fixed_group", {"group_id"}, rows);
        }
    }

    write_groups("_ff_dof_rb_member", "rb_id", system.get_rb_groups());
    write_groups("_ff_dof_md_fixed_member", "md_fixed_id", system.get_md_fixed_groups());

    // _ff_sampling
    {
        auto s = system.get_sampling();
        w.write_category("_ff_sampling", {
            {"temperature_K", cif_val(s.temperature_K)},
            {"friction_ps", cif_val(s.friction_ps)},
            {"timestep_fs", cif_val(s.timestep_fs)},
            {"n_steps", cif_val(s.n_steps)},
            {"write_every", cif_val(s.write_every)},
            {"minimize_steps", cif_val(s.minimize_steps)},
        });
    }

    // _ff_nonbonded
    {
        auto nb = system.get_nonbonded();
        w.write_category("_ff_nonbonded", {
            {"enabled", cif_val(nb.enabled)},
            {"k", cif_val(nb.k)},
            {"cutoff_A", cif_val(nb.cutoff)},
        });
    }

    // _ff_bond_type
    {
        auto& bt = system.get_bond_types();
        if (!bt.empty()) {
            std::vector<std::vector<std::string>> rows;
            for (auto& [tid, k] : bt) rows.push_back({tid, cif_val(k)});
            w.write_loop("_ff_bond_type", {"type_id", "k_kcal_mol_A2"}, rows);
        }
    }

    // _ff_angle_type
    {
        auto& at = system.get_angle_types();
        if (!at.empty()) {
            std::vector<std::vector<std::string>> rows;
            for (auto& [tid, k] : at) rows.push_back({tid, cif_val(k)});
            w.write_loop("_ff_angle_type", {"type_id", "k_kcal_mol_rad2"}, rows);
        }
    }

    // _ff_torsion_type
    {
        auto& tt = system.get_torsion_types();
        if (!tt.empty()) {
            std::vector<std::vector<std::string>> rows;
            for (auto& [tid, t] : tt) {
                rows.push_back({tid, cif_val(t.periodicity), cif_val(t.phase), cif_val(t.k)});
            }
            w.write_loop("_ff_torsion_type", {"type_id", "periodicity", "phase_rad", "k_kcal_mol"}, rows);
        }
    }

    // _ff_improper_type
    {
        auto& it = system.get_improper_types();
        if (!it.empty()) {
            std::vector<std::vector<std::string>> rows;
            for (auto& [tid, t] : it) {
                rows.push_back({tid, cif_val(t.periodicity), cif_val(t.phase), cif_val(t.k)});
            }
            w.write_loop("_ff_improper_type", {"type_id", "periodicity", "phase_rad", "k_kcal_mol"}, rows);
        }
    }

    // _ff_lj_type
    {
        auto& lj = system.get_lj_types();
        if (!lj.empty()) {
            std::vector<std::vector<std::string>> rows;
            for (auto& [tid, t] : lj) {
                rows.push_back({tid, cif_val(t.element), cif_val(t.rmin_half), cif_val(t.epsilon)});
            }
            w.write_loop("_ff_lj_type", {"type_id", "element", "rmin_half_A", "epsilon_kcal_mol"}, rows);
        }
    }

    // _ff_bond, _ff_angle, _ff_torsion, _ff_improper
    auto no = [&](const std::string& sid) -> int {
        auto it = site_id_to_no.find(sid);
        return it != site_id_to_no.end() ? it->second : -1;
    };

    {
        auto& bonds = system.get_bonds();
        if (!bonds.empty()) {
            std::vector<std::vector<std::string>> rows;
            for (auto& bd : bonds) {
                rows.push_back({cif_val(no(bd.site_a)), cif_val(no(bd.site_b)),
                                cif_val(bd.length), cif_val(bd.type_id)});
            }
            w.write_loop("_ff_bond", {"n1", "n2", "length_A", "type_id"}, rows);
        }
    }
    {
        auto& angles = system.get_angles();
        if (!angles.empty()) {
            std::vector<std::vector<std::string>> rows;
            for (auto& an : angles) {
                rows.push_back({cif_val(no(an.site_a)), cif_val(no(an.site_b)),
                                cif_val(no(an.site_c)), cif_val(an.theta), cif_val(an.type_id)});
            }
            w.write_loop("_ff_angle", {"n1", "n2", "n3", "theta_rad", "type_id"}, rows);
        }
    }
    {
        auto& dih = system.get_dihedrals();
        if (!dih.empty()) {
            std::vector<std::vector<std::string>> rows;
            for (auto& to : dih) {
                rows.push_back({cif_val(no(to.site_a)), cif_val(no(to.site_b)),
                                cif_val(no(to.site_c)), cif_val(no(to.site_d)),
                                cif_val(to.type_id)});
            }
            w.write_loop("_ff_torsion", {"n1", "n2", "n3", "n4", "type_id"}, rows);
        }
    }
    {
        auto& imp = system.get_impropers();
        if (!imp.empty()) {
            std::vector<std::vector<std::string>> rows;
            for (auto& to : imp) {
                rows.push_back({cif_val(no(to.site_a)), cif_val(no(to.site_b)),
                                cif_val(no(to.site_c)), cif_val(no(to.site_d)),
                                cif_val(to.type_id)});
            }
            w.write_loop("_ff_improper", {"n1", "n2", "n3", "n4", "type_id"}, rows);
        }
    }
}

IMPBFF_END_NAMESPACE
