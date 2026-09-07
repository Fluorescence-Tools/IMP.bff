/**
 *  \file ComponentTemplate.cpp
 *  \brief The `_cgprobe_*` template, read through IMP's ihm parser and
 *         written through this module's one CIF writer.
 */

#include <IMP/bff/ComponentTemplate.h>

#include <IMP/bff/internal/CifReader.h>

#include <IMP/bff/Base.h>

#include "ihm_format.h"

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

using internal::flag;
using internal::has;
using internal::integer;
using internal::txt;

// --------------------------------------------------------------------------
// Component template CIF reader (through the ihm C reader)
// --------------------------------------------------------------------------

namespace {

// The template CIF reader: same ihm callback pattern ForceFieldCIF.cpp uses.

struct TemplateCtx {
    ComponentTemplate* t;
    //! Whether any `_cgprobe_*` category fired. An unnamed template takes the
    //! file's stem for a name, so a name is no evidence that anything parsed.
    bool saw_category;
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
    c->saw_category = true;
    const std::string n = txt(c->name);
    if (!n.empty()) c->t->name = n;
}

void tc_on_feature(ihm_reader*, int, void* d, ihm_error**) {
    TemplateCtx* c = (TemplateCtx*) d;
    c->saw_category = true;
    const std::string fid = txt(c->f_id);
    if (fid.empty()) return;
    ComponentTemplate::Feature& f = c->feature(fid);
    const std::string ft = txt(c->f_type);
    if (!ft.empty()) f.feature_type = ft;
    f.rb = flag(c->f_rb, false);
    f.md_fixed = flag(c->f_md, false);
    const std::string rc = txt(c->f_color);
    if (!rc.empty() && rc != ".") f.region_color = rc;
}

void tc_on_feature_atom(ihm_reader*, int, void* d, ihm_error**) {
    TemplateCtx* c = (TemplateCtx*) d;
    c->saw_category = true;
    const std::string fid = txt(c->fa_id);
    const std::string anm = txt(c->fa_atom);
    if (fid.empty() || anm.empty()) return;
    ComponentTemplate::FeatureAtom a;
    a.name = anm;
    a.occurrence = integer(c->fa_occ, 1);
    c->feature(fid).atoms.push_back(a);
}

void tc_on_improper(ihm_reader*, int, void* d, ihm_error**) {
    TemplateCtx* c = (TemplateCtx*) d;
    c->saw_category = true;
    const std::string ca = txt(c->i_center);
    const std::string it = txt(c->i_type);
    if (ca.empty() || it.empty()) return;
    ComponentTemplate::Improper imp;
    imp.center_atom = ca;
    imp.type = it;
    c->t->impropers.push_back(imp);
}

void tc_on_metadata(ihm_reader*, int, void* d, ihm_error**) {
    TemplateCtx* c = (TemplateCtx*) d;
    c->saw_category = true;
    const std::string k = txt(c->m_key);
    const std::string v = txt(c->m_value);
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
                                                bool with_probe_metadata) {
    ComponentTemplate tmpl;
    // An unnamed template takes its name from the file stem.
    auto slash = path.find_last_of('/');
    auto basename = (slash != std::string::npos) ? path.substr(slash + 1) : path;
    auto dot = basename.find_last_of('.');
    tmpl.name = (dot != std::string::npos) ? basename.substr(0, dot) : basename;

    FILE* fh = std::fopen(path.c_str(), "r");
    if (!fh) IMP_THROW("Cannot read " << path, IOException);
    const int fd = fileno(fh);

    TemplateCtx c;
    c.t = &tmpl;
    c.saw_category = false;
    ihm_file* ifile = ihm_file_new_from_fd(fd);
    ihm_reader* reader = ihm_reader_new(ifile, false);  // text, not binary

#define TC_CAT(name, fn) ihm_category_new(reader, name, fn, NULL, NULL, &c, NULL)
    ihm_category* cat;
    cat = TC_CAT("_cgprobe_template", tc_on_template);
    c.name = ihm_keyword_str_new(cat, "name");

    cat = TC_CAT("_cgprobe_feature", tc_on_feature);
    c.f_id = ihm_keyword_str_new(cat, "feature_id");
    c.f_type = ihm_keyword_str_new(cat, "feature_type");
    c.f_rb = ihm_keyword_bool_new(cat, "rb");
    c.f_md = ihm_keyword_bool_new(cat, "md_fixed");
    c.f_color = ihm_keyword_str_new(cat, "region_color");

    cat = TC_CAT("_cgprobe_feature_atom", tc_on_feature_atom);
    c.fa_id = ihm_keyword_str_new(cat, "feature_id");
    c.fa_atom = ihm_keyword_str_new(cat, "atom_name");
    c.fa_occ = ihm_keyword_int_new(cat, "occurrence");

    cat = TC_CAT("_cgprobe_improper", tc_on_improper);
    c.i_center = ihm_keyword_str_new(cat, "center_atom");
    c.i_type = ihm_keyword_str_new(cat, "type");

    if (with_probe_metadata) {
        cat = TC_CAT("_cgprobe_metadata", tc_on_metadata);
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
    if (!c.saw_category) {
        IMP_THROW("no `_cgprobe_*` categories in "
                          << path
                          << " -- a template written before 2026-08-27 spells "
                          << "them `_cgdye_*`; `sed -i '' 's/_cgdye_/_cgprobe_/g'`"
                          << " migrates it",
                  IOException);
    }
    return tmpl;
}


ComponentTemplate read_probe_template_cif(const std::string& path) {
    return read_component_template_cif(path, true);
}


std::map<std::string, std::string> region_features(const ComponentTemplate& tmpl) {
    std::map<std::string, std::string> out;
    for (auto& [fid, spec] : tmpl.features) {
        if (!spec.region_color.empty()) out[fid] = spec.region_color;
    }
    return out;
}

IMPBFF_END_NAMESPACE
