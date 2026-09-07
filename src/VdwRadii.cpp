/**
 *  \file src/VdwRadii.cpp
 *  \brief Olga's van der Waals table, vendored.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/bff/VdwRadii.h>

#include <IMP/atom/Atom.h>
#include <IMP/bff/Base.h>

#include <algorithm>
#include <cstdio>
#include <set>
#include <sstream>

IMPBFF_BEGIN_NAMESPACE

/* Olga's radii, converted from the upstream nanometres to Angstrom **once,
   here**, so that no call site can forget the factor. The conversion is not a
   guess about magnitudes: `Olga/src/AV/Position.cpp:118-124` scales the
   coordinate and the radius by the same 10.0f on one line when it builds the
   (x, y, z, vdW) array `calculateAV()` consumes, because pteros stores
   coordinates in nm and the AV kernel works in Angstrom. See VdwRadii.h.

   Six distinct values over 128 names; the upstream file has 131 lines, of
   which "H1", "H2" and "H3" appear twice with identical values (JSON keeps the
   last, so the duplicates are inert). File-local helpers in this translation
   unit are prefixed `vdwradii_`: src/*.cpp is compiled as one unit here, so an
   anonymous namespace isolates nothing. */
static const std::map<std::string, double> &vdwradii_olga_table() {
    static const std::map<std::string, double> table = {
    // hydrogen (0.100 nm). Olga rasterises explicit hydrogens; a structure
    // read with a heavy-atom selector simply has none of these, which is the
    // single largest difference between this set and IMP's united-atom one.
    {"H", 1.0000}, {"H1", 1.0000}, {"H1'", 1.0000}, {"H2", 1.0000},
    {"H2'", 1.0000}, {"H21", 1.0000}, {"H22", 1.0000}, {"H3", 1.0000},
    {"H3'", 1.0000}, {"H4'", 1.0000}, {"H41", 1.0000}, {"H42", 1.0000},
    {"H5", 1.0000}, {"H5'", 1.0000}, {"H5''", 1.0000}, {"H6", 1.0000},
    {"H61", 1.0000}, {"H62", 1.0000}, {"H8", 1.0000}, {"HA", 1.0000},
    {"HA2", 1.0000}, {"HA3", 1.0000}, {"HB", 1.0000}, {"HB1", 1.0000},
    {"HB2", 1.0000}, {"HB3", 1.0000}, {"HD1", 1.0000}, {"HD11", 1.0000},
    {"HD12", 1.0000}, {"HD13", 1.0000}, {"HD2", 1.0000}, {"HD21", 1.0000},
    {"HD22", 1.0000}, {"HD23", 1.0000}, {"HD3", 1.0000}, {"HE", 1.0000},
    {"HE1", 1.0000}, {"HE2", 1.0000}, {"HE21", 1.0000}, {"HE22", 1.0000},
    {"HE3", 1.0000}, {"HG", 1.0000}, {"HG1", 1.0000}, {"HG11", 1.0000},
    {"HG12", 1.0000}, {"HG13", 1.0000}, {"HG2", 1.0000}, {"HG21", 1.0000},
    {"HG22", 1.0000}, {"HG23", 1.0000}, {"HG3", 1.0000}, {"HH", 1.0000},
    {"HH11", 1.0000}, {"HH12", 1.0000}, {"HH2", 1.0000}, {"HH21", 1.0000},
    {"HH22", 1.0000}, {"HO2'", 1.0000}, {"HO3'", 1.0000}, {"HO5'", 1.0000},
    {"HZ", 1.0000}, {"HZ1", 1.0000}, {"HZ2", 1.0000}, {"HZ3", 1.0000},
    // oxygen (0.1490 nm)
    {"O", 1.4900}, {"O2", 1.4900}, {"O2'", 1.4900}, {"O3'", 1.4900},
    {"O4", 1.4900}, {"O4'", 1.4900}, {"O5'", 1.4900}, {"O6", 1.4900},
    {"OD1", 1.4900}, {"OD2", 1.4900}, {"OE1", 1.4900}, {"OE2", 1.4900},
    {"OG", 1.4900}, {"OG1", 1.4900}, {"OH", 1.4900}, {"OP1", 1.4900},
    {"OP2", 1.4900}, {"OXT", 1.4900},
    // nitrogen (0.1625 nm)
    {"N", 1.6250}, {"N1", 1.6250}, {"N2", 1.6250}, {"N3", 1.6250},
    {"N4", 1.6250}, {"N6", 1.6250}, {"N7", 1.6250}, {"N9", 1.6250},
    {"ND1", 1.6250}, {"ND2", 1.6250}, {"NE", 1.6250}, {"NE1", 1.6250},
    {"NE2", 1.6250}, {"NH1", 1.6250}, {"NH2", 1.6250}, {"NZ", 1.6250},
    // carbon (0.1700 nm) = Bondi
    {"C", 1.7000}, {"C1'", 1.7000}, {"C2", 1.7000}, {"C2'", 1.7000},
    {"C3'", 1.7000}, {"C4", 1.7000}, {"C4'", 1.7000}, {"C5", 1.7000},
    {"C5'", 1.7000}, {"C6", 1.7000}, {"C8", 1.7000}, {"CA", 1.7000},
    {"CB", 1.7000}, {"CD", 1.7000}, {"CD1", 1.7000}, {"CD2", 1.7000},
    {"CE", 1.7000}, {"CE1", 1.7000}, {"CE2", 1.7000}, {"CE3", 1.7000},
    {"CG", 1.7000}, {"CG1", 1.7000}, {"CG2", 1.7000}, {"CH2", 1.7000},
    {"CZ", 1.7000}, {"CZ2", 1.7000}, {"CZ3", 1.7000},
    // sulphur (0.1782 nm)
    {"SD", 1.7820}, {"SG", 1.7820},
    // phosphorus (0.1860 nm). Worth noticing: Olga's *built-in* fallback map
    // in loadvdWRadii() carries P = 0.1 nm (= 1.0 A, a hydrogen), and the
    // shipped JSON overrides it with 0.186 nm. The JSON is what a real run
    // loads, so 1.86 A is Olga's phosphorus and 1.0 A is a bug in the
    // fallback nobody reaches.
    {"P", 1.8600},
    };
    return table;
}

const std::map<std::string, double> &get_olga_vdw_radii() {
    return vdwradii_olga_table();
}

double get_olga_vdw_fallback_radius() {
    /* `vdWRMap.value(QString::fromStdString(system.atom(i).name), 0.15)`
       -- Olga/src/AV/Position.cpp:110, in nanometres, so 1.5 A. It is a
       property of Olga's *code*, not of its table, which is why it is a
       constant here and not a row in data/olga_vdw_radii.csv. */
    return 1.5;
}

double olga_vdw_radius(const std::string &atom_name) {
    const std::map<std::string, double> &t = vdwradii_olga_table();
    std::map<std::string, double>::const_iterator it = t.find(atom_name);
    if (it == t.end()) return get_olga_vdw_fallback_radius();
    return it->second;
}

unsigned int get_number_of_olga_vdw_radii() {
    return (unsigned int) vdwradii_olga_table().size();
}

std::vector<std::string> get_olga_vdw_atom_names() {
    std::vector<std::string> out;
    out.reserve(vdwradii_olga_table().size());
    for (const auto &kv : vdwradii_olga_table()) out.push_back(kv.first);
    return out;
}

double olga_vdw_particle_radius(IMP::Particle *p) {
    if (p == nullptr) return get_olga_vdw_fallback_radius();
    if (!IMP::atom::Atom::get_is_setup(p)) {
        return get_olga_vdw_fallback_radius();
    }
    return olga_vdw_radius(
            IMP::atom::Atom(p).get_atom_type().get_string());
}

std::vector<double> olga_vdw_radii(const IMP::ParticlesTemp &ps) {
    std::vector<double> out;
    out.reserve(ps.size());
    for (std::size_t i = 0; i < ps.size(); ++i) {
        out.push_back(olga_vdw_particle_radius(ps[i]));
    }
    return out;
}

std::vector<std::string> olga_vdw_unknown_atom_names(
        const IMP::ParticlesTemp &ps) {
    const std::map<std::string, double> &t = vdwradii_olga_table();
    std::set<std::string> missing;
    for (std::size_t i = 0; i < ps.size(); ++i) {
        if (!IMP::atom::Atom::get_is_setup(ps[i])) continue;
        const std::string n =
                IMP::atom::Atom(ps[i]).get_atom_type().get_string();
        if (t.find(n) == t.end()) missing.insert(n);
    }
    return std::vector<std::string>(missing.begin(), missing.end());
}

AVRadiiSource av_radii_source_from_string(const std::string &s) {
    if (s == "imp") return AV_RADII_IMP;
    if (s == "olga") return AV_RADII_OLGA;
    /* "model" was this value's spelling under schema 1.5, for one day. It is
       rejected rather than aliased so that a file carrying it is read as the
       1.5 file it is, loudly, instead of being silently reinterpreted. */
    IMP_THROW("radii source must be \"imp\" or \"olga\", not \"" << s << "\"",
              IMP::ValueException);
}

std::string av_radii_source_to_string(AVRadiiSource s) {
    return s == AV_RADII_OLGA ? "olga" : "imp";
}

std::string olga_vdw_radii_csv() {
    /* data/olga_vdw_radii.csv is *derived* from the table above, the way
       data/fps_json_schema.json is derived from FPSSchema.h's tables: the C++
       is the definition, the data file is the shipped, human-readable record
       of provenance, and a test regenerates it and fails on drift. */
    static const char *header = R"CSVHDR(# Olga's van der Waals radii, keyed by PDB atom name.
#
# Provenance
# ----------
#   upstream file : Olga/src/vdWRadii.json  (128 entries; 131 lines, of which
#                   "H1", "H2" and "H3" appear twice with identical values)
#   upstream repo : Fluorescence-Tools/Olga -- the FRET-restrained structural
#                   modelling program by Mykola Dimura et al., the program that
#                   produced the Zenodo 3376527 T4L screening table
#                   (Sanabria et al., Nat. Commun. 11, 1231 (2020)).
#   local checkout: ../ucfret/thirdparty/olga/src/vdWRadii.json
#   loader        : Olga/src/AV/Position.cpp:79-110, loadvdWRadii/pterosVDW.
#   licence       : Olga is GPL-3.0. Only the numeric table is vendored here,
#                   as data; no Olga code is copied into this module.
#
# Unit -- nanometre upstream, Angstrom here
# -----------------------------------------
# The upstream JSON is in NANOMETRES: pteros stores coordinates in nm, and
# Olga converts both coordinate and radius on the same line when it builds the
# (x, y, z, vdW) array the AV kernel consumes:
#
#     Olga/src/AV/Position.cpp:118-124, coordsVdW()
#         xyzw.emplace_back(frame.coord.at(i)[0] * 10.0f,
#                           frame.coord.at(i)[1] * 10.0f,
#                           frame.coord.at(i)[2] * 10.0f,
#                           pterosVDW(system, i) * 10.0f);
#
# so calculateAV() -- and every linker length, linker width and dye radius
# beside it -- works in ANGSTROM. IMP.bff works in Angstrom too, so the table
# is stored here already multiplied by 10: carbon 0.17 nm -> 1.70 A, which is
# the Bondi carbon radius and the same number FPS's own vdW.txt carries.
#
# Unknown atom names
# ------------------
# Olga does NOT fall back to an element lookup. pterosVDW is
#     vdWRMap.value(name, 0.15)
# (Position.cpp:110), i.e. a flat 0.15 nm = 1.50 A for any atom name not in
# this table -- smaller than every heavy atom in it. That fallback is Olga's
# code, not its table, so it is not a row here; see
# IMP::bff::get_olga_vdw_fallback_radius().
#
# This file is DERIVED from the C++ table in src/VdwRadii.cpp by
# IMP::bff::olga_vdw_radii_csv(); test/representation/test_olga_vdw_radii.py
# regenerates it and fails on drift. Edit the C++ table, not this file.
#
atom_name,radius_angstrom
)CSVHDR";
    std::ostringstream out;
    out << header;
    char buf[64];
    for (const auto &kv : vdwradii_olga_table()) {
        std::snprintf(buf, sizeof(buf), "%.4f", kv.second);
        out << kv.first << "," << buf << "\n";
    }
    return out.str();
}

IMPBFF_END_NAMESPACE
