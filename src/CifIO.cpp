/**
 *  \file CifIO.cpp
 *  \brief Writing compact FF topology mmCIF tables, template CIF, and rotamer
 *         library IO — the remaining io/cif.py surface, in C++.
 */

#include <IMP/bff/CifIO.h>
#include <IMP/bff/DyeForceField.h>
#include <IMP/bff/Mol2IO.h>

#include <IMP/exception.h>

#include "ihm_format.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <ios>
#include <regex>
#include <set>
#include <sstream>
#include <string>

IMPBFF_BEGIN_NAMESPACE

// --------------------------------------------------------------------------
// mmCIF writer (replaces ihm.format.CifWriter)
// --------------------------------------------------------------------------

namespace {

class CifWriter {
public:
    explicit CifWriter(std::ostream& out) : out_(out) {}
    void start_block(const std::string& name) {
        out_ << "data_" << name << "\n";
    }
    void write_category(const std::string& name,
                         const std::vector<std::pair<std::string, std::string>>& kv) {
        out_ << "#" << "\n";
        for (auto& [k, v] : kv) {
            out_ << name << "." << k << " " << quote(v) << "\n";
        }
        out_ << "#" << "\n";
    }
    void write_loop(const std::string& name,
                    const std::vector<std::string>& cols,
                    const std::vector<std::vector<std::string>>& rows) {
        out_ << "#" << "\n" << "loop_" << "\n";
        for (auto& c : cols) out_ << name << "." << c << "\n";
        for (auto& row : rows) {
            for (size_t i = 0; i < row.size(); i++) {
                if (i) out_ << " ";
                out_ << quote(row[i]);
            }
            out_ << "\n";
        }
        out_ << "#" << "\n";
    }
private:
    //! A value the ihm reader will tokenize as one token. Quoting here (and
    //! not at each call site) is what makes a metadata value like
    //! `N1 and resname A48, N2 and resname A48` survive a round trip -- the
    //! reader splits rows on whitespace, so an unquoted multi-word value is
    //! read as several and the loop comes up short. Already-quoted and the CIF
    //! nulls (`.`/`?`) pass through untouched.
    static std::string quote(const std::string& s) {
        if (s.empty() || s == "." || s == "?") return s;
        if (s[0] == '"' || s[0] == '\'') return s;
        const bool plain = s.find_first_of(" \t,#'\"") == std::string::npos
                && s[0] != '_' && s[0] != ';';
        if (plain) return s;
        std::string escaped;
        for (char ch : s) {
            if (ch == '"') escaped += "\"\"";
            else escaped += ch;
        }
        return "\"" + escaped + "\"";
    }
    std::ostream& out_;
};

std::string cif_val(const std::string& s) {
    if (s.empty()) return ".";
    // Quote if contains space, #, or starts with _
    if (s.find(' ') != std::string::npos || s.find('#') != std::string::npos
        || (!s.empty() && s[0] == '_')) {
        return "\"" + s + "\"";
    }
    return s;
}

std::string cif_val(double v) {
    std::ostringstream ss;
    ss << v;
    return ss.str();
}

std::string cif_val(int v) {
    return std::to_string(v);
}

std::string cif_val(bool v) {
    return v ? "YES" : "NO";
}

std::string cif_val_or_omit(const std::string& s) {
    if (s.empty()) return ".";
    return cif_val(s);
}

}  // namespace

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

void write_dye_forcefield_cif(const std::string& path,
                                const DyeForceFieldSystem& system) {
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

// --------------------------------------------------------------------------
// Component template CIF reader (through the ihm C reader)
// --------------------------------------------------------------------------

namespace {

// The template CIF reader: same ihm callback pattern ForceFieldCIF.cpp uses.
bool tc_has(ihm_keyword* k) {
    return k && k->in_file && !k->omitted && !k->unknown;
}
std::string tc_str(ihm_keyword* k) {
    return tc_has(k) && k->data.str ? std::string(k->data.str) : std::string();
}
int tc_int(ihm_keyword* k, int fallback) {
    return tc_has(k) ? k->data.ival : fallback;
}
bool tc_bool(ihm_keyword* k) { return tc_has(k) && k->data.bval; }

struct TemplateCtx {
    ComponentTemplate* t;
    ihm_keyword *name;
    ihm_keyword *f_id, *f_type, *f_rb, *f_md, *f_color;
    ihm_keyword *fa_id, *fa_atom, *fa_occ;
    ihm_keyword *i_center, *i_type;
    ihm_keyword *m_key, *m_value;

    ComponentTemplate::Feature& feature(const std::string& fid) {
        ComponentTemplate::Feature& f = t->features[fid];
        if (f.feature_type.empty()) f.feature_type = "dof";
        return f;
    }
};

void tc_on_template(ihm_reader*, int, void* d, ihm_error**) {
    TemplateCtx* c = (TemplateCtx*) d;
    const std::string n = tc_str(c->name);
    if (!n.empty()) c->t->name = n;
}

void tc_on_feature(ihm_reader*, int, void* d, ihm_error**) {
    TemplateCtx* c = (TemplateCtx*) d;
    const std::string fid = tc_str(c->f_id);
    if (fid.empty()) return;
    ComponentTemplate::Feature& f = c->feature(fid);
    const std::string ft = tc_str(c->f_type);
    if (!ft.empty()) f.feature_type = ft;
    f.rb = tc_bool(c->f_rb);
    f.md_fixed = tc_bool(c->f_md);
    const std::string rc = tc_str(c->f_color);
    if (!rc.empty() && rc != ".") f.region_color = rc;
}

void tc_on_feature_atom(ihm_reader*, int, void* d, ihm_error**) {
    TemplateCtx* c = (TemplateCtx*) d;
    const std::string fid = tc_str(c->fa_id);
    const std::string anm = tc_str(c->fa_atom);
    if (fid.empty() || anm.empty()) return;
    ComponentTemplate::FeatureAtom a;
    a.name = anm;
    a.occurrence = tc_int(c->fa_occ, 1);
    c->feature(fid).atoms.push_back(a);
}

void tc_on_improper(ihm_reader*, int, void* d, ihm_error**) {
    TemplateCtx* c = (TemplateCtx*) d;
    const std::string ca = tc_str(c->i_center);
    const std::string it = tc_str(c->i_type);
    if (ca.empty() || it.empty()) return;
    ComponentTemplate::Improper imp;
    imp.center_atom = ca;
    imp.type = it;
    c->t->impropers.push_back(imp);
}

void tc_on_metadata(ihm_reader*, int, void* d, ihm_error**) {
    TemplateCtx* c = (TemplateCtx*) d;
    const std::string k = tc_str(c->m_key);
    const std::string v = tc_str(c->m_value);
    if (k.empty() || v.empty()) return;
    ComponentTemplate* t = c->t;
    if (k == "center_atom") {
        t->center_atom = v;
    } else if (k == "dipole_atom_1") {
        t->dipole_atom_1 = v;
    } else if (k == "dipole_atom_2") {
        t->dipole_atom_2 = v;
    } else if (k == "positive_atoms" || k == "negative_atoms") {
        // Comma-separated on the wire; split on any whitespace around commas.
        std::vector<std::string>& out =
                k == "positive_atoms" ? t->positive_atoms : t->negative_atoms;
        std::string tok;
        std::istringstream in(v);
        while (std::getline(in, tok, ',')) {
            const std::size_t a = tok.find_first_not_of(" \t");
            const std::size_t b = tok.find_last_not_of(" \t");
            if (a != std::string::npos) out.push_back(tok.substr(a, b - a + 1));
        }
    }
}

}  // namespace

ComponentTemplate read_component_template_cif(const std::string& path,
                                                bool with_dye_metadata) {
    ComponentTemplate tmpl;
    // The Python reader named unnamed templates after the file stem.
    auto slash = path.find_last_of('/');
    auto basename = (slash != std::string::npos) ? path.substr(slash + 1) : path;
    auto dot = basename.find_last_of('.');
    tmpl.name = (dot != std::string::npos) ? basename.substr(0, dot) : basename;

    FILE* fh = std::fopen(path.c_str(), "r");
    if (!fh) IMP_THROW("Cannot read " << path, IOException);
    const int fd = fileno(fh);

    TemplateCtx c;
    c.t = &tmpl;
    ihm_file* ifile = ihm_file_new_from_fd(fd);
    ihm_reader* reader = ihm_reader_new(ifile, false);  // text, not binary

#define TC_CAT(name, fn) ihm_category_new(reader, name, fn, NULL, NULL, &c, NULL)
    ihm_category* cat;
    cat = TC_CAT("_cgdye_template", tc_on_template);
    c.name = ihm_keyword_str_new(cat, "name");

    cat = TC_CAT("_cgdye_feature", tc_on_feature);
    c.f_id = ihm_keyword_str_new(cat, "feature_id");
    c.f_type = ihm_keyword_str_new(cat, "feature_type");
    c.f_rb = ihm_keyword_bool_new(cat, "rb");
    c.f_md = ihm_keyword_bool_new(cat, "md_fixed");
    c.f_color = ihm_keyword_str_new(cat, "region_color");

    cat = TC_CAT("_cgdye_feature_atom", tc_on_feature_atom);
    c.fa_id = ihm_keyword_str_new(cat, "feature_id");
    c.fa_atom = ihm_keyword_str_new(cat, "atom_name");
    c.fa_occ = ihm_keyword_int_new(cat, "occurrence");

    cat = TC_CAT("_cgdye_improper", tc_on_improper);
    c.i_center = ihm_keyword_str_new(cat, "center_atom");
    c.i_type = ihm_keyword_str_new(cat, "type");

    if (with_dye_metadata) {
        cat = TC_CAT("_cgdye_metadata", tc_on_metadata);
        c.m_key = ihm_keyword_str_new(cat, "key");
        c.m_value = ihm_keyword_str_new(cat, "value");
    }
#undef TC_CAT

    bool more = false;
    ihm_error* err = NULL;
    const bool ok = ihm_read_file(reader, &more, &err);
    if (!ok) {
        const std::string message = err && err->msg ? err->msg : "parse failed";
        if (err) ihm_error_free(err);
        ihm_reader_free(reader);
        std::fclose(fh);
        IMP_THROW("reading " << path << ": " << message, IOException);
    }
    ihm_reader_free(reader);
    std::fclose(fh);
    return tmpl;
}

void write_component_template_cif(const std::string& path,
                                    const ComponentTemplate& tmpl) {
    std::ofstream out(path);
    if (!out) IMP_THROW("cannot open " << path << " for writing", IOException);
    CifWriter w(out);
    std::string name = tmpl.name.empty() ? "cgdye_template" : tmpl.name;
    w.start_block(name);

    // _cgdye_template
    w.write_category("_cgdye_template", {{"name", name}});

    // _cgdye_feature
    if (!tmpl.features.empty()) {
        std::vector<std::vector<std::string>> rows;
        for (auto& [fid, spec] : tmpl.features) {
            rows.push_back({fid, spec.feature_type,
                            spec.rb ? "YES" : "NO",
                            spec.md_fixed ? "YES" : "NO",
                            spec.region_color.empty() ? "." : spec.region_color});
        }
        w.write_loop("_cgdye_feature",
                      {"feature_id", "feature_type", "rb", "md_fixed", "region_color"}, rows);
    }

    // _cgdye_feature_atom
    {
        std::vector<std::vector<std::string>> rows;
        for (auto& [fid, spec] : tmpl.features) {
            for (auto& a : spec.atoms) {
                rows.push_back({fid, a.name, std::to_string(a.occurrence)});
            }
        }
        if (!rows.empty())
            w.write_loop("_cgdye_feature_atom", {"feature_id", "atom_name", "occurrence"}, rows);
    }

    // _cgdye_improper
    {
        if (!tmpl.impropers.empty()) {
            std::vector<std::vector<std::string>> rows;
            for (auto& imp : tmpl.impropers) {
                rows.push_back({imp.center_atom, imp.type});
            }
            w.write_loop("_cgdye_improper", {"center_atom", "type"}, rows);
        }
    }
}

ComponentTemplate read_dye_template_cif(const std::string& path) {
    return read_component_template_cif(path, true);
}

void write_dye_template_cif(const std::string& path, const ComponentTemplate& tmpl) {
    std::ofstream out(path);
    if (!out) IMP_THROW("cannot open " << path << " for writing", IOException);
    CifWriter w(out);
    std::string name = tmpl.name.empty() ? "dye_template" : tmpl.name;
    w.start_block(name);

    // _cgdye_template
    w.write_category("_cgdye_template", {{"name", name}});

    // _cgdye_metadata
    {
        std::vector<std::vector<std::string>> rows;
        if (!tmpl.center_atom.empty()) rows.push_back({"center_atom", tmpl.center_atom});
        if (!tmpl.dipole_atom_1.empty()) rows.push_back({"dipole_atom_1", tmpl.dipole_atom_1});
        if (!tmpl.dipole_atom_2.empty()) rows.push_back({"dipole_atom_2", tmpl.dipole_atom_2});
        if (!tmpl.positive_atoms.empty()) {
            std::string joined;
            for (size_t i = 0; i < tmpl.positive_atoms.size(); i++) {
                if (i) joined += ", ";
                joined += tmpl.positive_atoms[i];
            }
            rows.push_back({"positive_atoms", joined});
        }
        if (!tmpl.negative_atoms.empty()) {
            std::string joined;
            for (size_t i = 0; i < tmpl.negative_atoms.size(); i++) {
                if (i) joined += ", ";
                joined += tmpl.negative_atoms[i];
            }
            rows.push_back({"negative_atoms", joined});
        }
        if (!rows.empty())
            w.write_loop("_cgdye_metadata", {"key", "value"}, rows);
    }

    // _cgdye_feature
    if (!tmpl.features.empty()) {
        std::vector<std::vector<std::string>> rows;
        for (auto& [fid, spec] : tmpl.features) {
            rows.push_back({fid, spec.feature_type,
                            spec.rb ? "YES" : "NO",
                            spec.md_fixed ? "YES" : "NO",
                            spec.region_color.empty() ? "." : spec.region_color});
        }
        w.write_loop("_cgdye_feature",
                      {"feature_id", "feature_type", "rb", "md_fixed", "region_color"}, rows);
    }

    // _cgdye_feature_atom
    {
        std::vector<std::vector<std::string>> rows;
        for (auto& [fid, spec] : tmpl.features) {
            for (auto& a : spec.atoms) {
                rows.push_back({fid, a.name, std::to_string(a.occurrence)});
            }
        }
        if (!rows.empty())
            w.write_loop("_cgdye_feature_atom", {"feature_id", "atom_name", "occurrence"}, rows);
    }

    // _cgdye_improper
    {
        if (!tmpl.impropers.empty()) {
            std::vector<std::vector<std::string>> rows;
            for (auto& imp : tmpl.impropers) {
                rows.push_back({imp.center_atom, imp.type});
            }
            w.write_loop("_cgdye_improper", {"center_atom", "type"}, rows);
        }
    }
}

std::map<std::string, std::string> region_features(const ComponentTemplate& tmpl) {
    std::map<std::string, std::string> out;
    for (auto& [fid, spec] : tmpl.features) {
        if (!spec.region_color.empty()) out[fid] = spec.region_color;
    }
    return out;
}

// --------------------------------------------------------------------------
// Rotamer library IO
// --------------------------------------------------------------------------

namespace {

std::string base_path(const std::string& p) {
    if (p.size() > 4 && (p.substr(p.size() - 4) == ".rmf" || p.substr(p.size() - 4) == ".npy"))
        return p.substr(0, p.size() - 4);
    return p;
}

// Minimal .npy reader for float64 arrays. The numpy .npy format is:
// magic \x93NUMPY, version, header length, header dict, then raw data.
std::vector<double> read_npy(const std::string& path, int& n0, int& n1, int& n2) {
    std::ifstream f(path, std::ios::binary);
    if (!f) IMP_THROW("cannot open " << path, IOException);
    char magic[6];
    f.read(magic, 6);
    if (std::string(magic, 6) != "\x93NUMPY") IMP_THROW("not a .npy file: " << path, IOException);
    char version[2];
    f.read(version, 2);
    int header_len = version[0] == 1 ? 0 : 0;
    if (version[0] == 1) { uint16_t h; f.read(reinterpret_cast<char*>(&h), 2); header_len = h; }
    else if (version[0] == 2) { uint32_t h; f.read(reinterpret_cast<char*>(&h), 4); header_len = h; }
    std::string header(header_len, '\0');
    f.read(&header[0], header_len);
    // Parse shape from header: "'shape': (N, M, 3),"
    // This is a minimal parser; numpy headers are Python dict literals.
    n0 = n1 = n2 = 0;
    auto shape_pos = header.find("'shape'");
    if (shape_pos == std::string::npos) shape_pos = header.find("\"shape\"");
    if (shape_pos != std::string::npos) {
        auto paren = header.find('(', shape_pos);
        if (paren != std::string::npos) {
            std::string s = header.substr(paren + 1);
            // parse comma-separated ints
            std::istringstream ss(s);
            char c;
            std::vector<int> dims;
            int d;
            while (ss >> d) {
                dims.push_back(d);
                ss >> c;
                if (c != ',') break;
            }
            if (dims.size() >= 1) n0 = dims[0];
            if (dims.size() >= 2) n1 = dims[1];
            if (dims.size() >= 3) n2 = dims[2];
        }
    }
    // Read data
    size_t count = 1;
    if (n0 > 0) count *= n0;
    if (n1 > 0) count *= n1;
    if (n2 > 0) count *= n2;
    std::vector<double> data(count);
    f.read(reinterpret_cast<char*>(data.data()), count * sizeof(double));
    return data;
}

void write_npy(const std::string& path, const std::vector<double>& data,
               int n0, int n1, int n2) {
    std::ofstream f(path, std::ios::binary);
    if (!f) IMP_THROW("cannot open " << path, IOException);
    // Build header
    std::string shape = "(";
    if (n0 > 0) shape += std::to_string(n0);
    if (n1 > 0) shape += ", " + std::to_string(n1);
    if (n2 > 0) shape += ", " + std::to_string(n2);
    if (n0 > 0 && n1 == 0) shape += ", ";  // numpy trailing comma for 1-D
    shape += ")";
    std::string dict = "{'descr': '<f8', 'fortran_order': False, 'shape': " + shape + ", }";
    // Pad to multiple of 64 bytes (numpy convention)
    int header_len = dict.size() + 1;  // +1 for newline
    int padded = ((header_len + 10 + 63) / 64) * 64;
    int pad = padded - 10 - dict.size() - 1;
    std::string padding(pad, ' ');
    std::string full_dict = dict + padding + "\n";

    f.write("\x93NUMPY", 6);
    char version[2] = {1, 0};
    f.write(version, 2);
    uint16_t hlen = full_dict.size();
    f.write(reinterpret_cast<char*>(&hlen), 2);
    f.write(full_dict.data(), full_dict.size());
    f.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(double));
}

}  // namespace

RotamerLibraryData read_rotamer_library(const std::string& path) {
    std::string base = base_path(path);
    RotamerLibraryData lib;
    int n0, n1, n2;
    lib.coords = read_npy(base + "_coords.npy", n0, n1, n2);
    lib.n_rotamers = n0;
    lib.n_atoms = n1;

    // weights
    std::string wpath = base + "_weights.txt";
    std::ifstream wf(wpath);
    if (wf) {
        std::string line;
        while (std::getline(wf, line)) {
            if (!line.empty()) lib.weight.push_back(std::stod(line));
        }
    } else {
        lib.weight.assign(lib.n_rotamers, 1.0 / lib.n_rotamers);
    }

    // atom names
    std::string apath = base + "_atoms.txt";
    std::ifstream af(apath);
    if (af) {
        std::string line;
        while (std::getline(af, line)) {
            if (!line.empty()) lib.atom_names.push_back(line);
        }
    } else {
        for (int i = 0; i < lib.n_atoms; i++)
            lib.atom_names.push_back("AT" + std::to_string(i));
    }

    lib.id.resize(lib.n_rotamers);
    for (int i = 0; i < lib.n_rotamers; i++) lib.id[i] = i + 1;
    return lib;
}

void write_rotamer_library(const std::string& path, const RotamerLibraryData& lib) {
    std::string base = base_path(path);
    write_npy(base + "_coords.npy", lib.coords, lib.n_rotamers, lib.n_atoms, 3);

    std::ofstream wf(base + "_weights.txt");
    for (auto& w : lib.weight) wf << w << "\n";

    std::ofstream af(base + "_atoms.txt");
    for (auto& name : lib.atom_names) af << name << "\n";
}

RotamerLibraryData normalize_weights(const RotamerLibraryData& lib) {
    RotamerLibraryData out = lib;
    double total = 0;
    for (double w : out.weight) total += w;
    if (total > 0) {
        for (double& w : out.weight) w /= total;
    }
    return out;
}

IMPBFF_END_NAMESPACE
