/**
 * \file ProbeLibrary.cpp
 * \brief The dye as a species: its spectra, and the Förster radius they give.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/internal/Text.h>
#include <IMP/bff/ProbeLibrary.h>

#include <IMP/bff/Base.h>

// Vendored with IMP and exported by libimp_atom; ForceFieldCIF.cpp reads its
// eighteen `_ff_*` categories through the same parser, in the same text mode.
#include <IMP/bff/internal/Cif.h>

#include "ihm_format.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fcntl.h>
#include <limits>
#include <sys/stat.h>
#include <unistd.h>

IMPBFF_BEGIN_NAMESPACE

using internal::dbl;
using internal::has;
using internal::txt;

std::string probe_type_to_string(ProbeType type) {
    switch (type) {
        case PROBE_DYE: return "dye";
        case PROBE_FLUORESCENT_PROTEIN: return "fluorescent_protein";
        case PROBE_SPIN_LABEL: return "spin_label";
        default: return "unspecified";
    }
}

ProbeType probe_type_from_string(std::string name) {
    if (name == "dye") return PROBE_DYE;
    if (name == "fluorescent_protein") return PROBE_FLUORESCENT_PROTEIN;
    if (name == "spin_label") return PROBE_SPIN_LABEL;
    return PROBE_UNSPECIFIED;
}

bool probe_is_fluorescent(ProbeType type) {
    // PROBE_UNSPECIFIED counts as fluorescent. A library that does not state
    // a type is a dye library -- that is what every one written before the
    // column existed is -- so refusing it would break callers over a field
    // they never had. A spin label, by contrast, has to be declared.
    return type != PROBE_SPIN_LABEL;
}

const double DEFAULT_REFRACTIVE_INDEX = 1.4;
const char* const PROBE_LIBRARY_CIF = "probe_library.cif";

// IMP compiles the module as one translation unit (`bff_all.cpp`), so an
// anonymous namespace is shared with every other .cpp -- and `ForceFieldCIF.cpp`
// reads the same parser with helpers of the same obvious names. The named inner
// namespace is what keeps `has`, `txt` and `Ctx` from colliding.
namespace probe_cif {
using IMP::bff::internal::nan_value;


//! Numerical factor of the Förster expression with R0 in nm, the overlap
//! integral in M^-1 cm^-1 nm^4 and wavelengths in nm. This is the constant as
//! the literature states it; #forster_radius converts its result once.
const double R0_FACTOR = 0.02108;

//! Spectra are nanometres -- that is the unit a spectrometer, a dye
//! datasheet and the overlap integral are all in, and the wavelength grids
//! stay in it. Coordinates are Angstrom, so an R0 that meets a distance is
//! Angstrom. This is where the two meet, once.
const double NM_TO_ANGSTROM = 10.0;


//! One point of one dye's curves, before the rows are grouped and sorted.
struct SpectrumPoint {
    double wavelength, excitation, emission;
    bool operator<(const SpectrumPoint& o) const { return wavelength < o.wavelength; }
};

struct ProbeRow {
    std::string name, type, vendor;
    double extinction_coefficient, quantum_yield;
};

struct Ctx {
    std::vector<ProbeRow> probes;
    std::map<std::string, std::vector<SpectrumPoint> > curves;

    ihm_keyword *d_name, *d_type, *d_vendor, *d_eps, *d_qy;
    ihm_keyword *s_name, *s_wavelength, *s_excitation, *s_emission;
};

void on_probe(ihm_reader*, int, void* d, ihm_error**) {
    Ctx* c = (Ctx*) d;
    if (!has(c->d_name)) return;
    ProbeRow row;
    row.name = txt(c->d_name);
    // Both optional: a library written before the columns existed still reads,
    // and an absent type is PROBE_UNSPECIFIED rather than a parse failure.
    row.type = has(c->d_type) ? txt(c->d_type) : std::string();
    row.vendor = has(c->d_vendor) ? txt(c->d_vendor) : std::string();
    row.extinction_coefficient = dbl(c->d_eps);
    row.quantum_yield = dbl(c->d_qy);
    c->probes.push_back(row);
}

void on_spectrum(ihm_reader*, int, void* d, ihm_error**) {
    Ctx* c = (Ctx*) d;
    if (!has(c->s_name)) return;
    SpectrumPoint p;
    p.wavelength = dbl(c->s_wavelength);
    p.excitation = dbl(c->s_excitation);
    p.emission = dbl(c->s_emission);
    c->curves[txt(c->s_name)].push_back(p);
}

//! Trapezoidal integral of `y` over the non-uniform axis `x`.
double trapezoid(const std::vector<double>& y, const std::vector<double>& x) {
    double total = 0.0;
    for (std::size_t i = 1; i < x.size(); ++i) {
        total += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]);
    }
    return total;
}

//! What a cached parse is keyed on: the resolved path, its mtime and its size.
struct CacheKey {
    std::string path;
    long long mtime_ns;
    long long size;
    bool operator<(const CacheKey& o) const {
        if (path != o.path) return path < o.path;
        if (mtime_ns != o.mtime_ns) return mtime_ns < o.mtime_ns;
        return size < o.size;
    }
};

std::map<CacheKey, std::map<std::string, Probe> >& library_cache() {
    static std::map<CacheKey, std::map<std::string, Probe> > cache;
    return cache;
}

std::map<std::string, Probe> parse_probe_library(const std::string& path) {
    const int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) IMP_THROW("cannot open " << path, IOException);

    Ctx c;
    ihm_file* fh = ihm_file_new_from_fd(fd);
    ihm_reader* reader = ihm_reader_new(fh, false);          // text, not binary

    ihm_category* cat =
        ihm_category_new(reader, "_bff_probe", on_probe, NULL, NULL, &c, NULL);
    c.d_name = ihm_keyword_str_new(cat, "chromophore_name");
    c.d_type = ihm_keyword_str_new(cat, "probe_type");
    c.d_vendor = ihm_keyword_str_new(cat, "vendor");
    c.d_eps = ihm_keyword_float_new(cat, "extinction_coefficient");
    c.d_qy = ihm_keyword_float_new(cat, "quantum_yield");

    cat = ihm_category_new(reader, "_bff_probe_spectrum", on_spectrum, NULL, NULL,
                           &c, NULL);
    c.s_name = ihm_keyword_str_new(cat, "chromophore_name");
    c.s_wavelength = ihm_keyword_float_new(cat, "wavelength");
    c.s_excitation = ihm_keyword_float_new(cat, "excitation");
    c.s_emission = ihm_keyword_float_new(cat, "emission");

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

    std::map<std::string, Probe> out;
    for (std::size_t i = 0; i < c.probes.size(); ++i) {
        const ProbeRow& row = c.probes[i];
        Probe dye;
        dye.name = row.name;
        dye.probe_type = probe_type_from_string(row.type);
        dye.vendor = row.vendor;
        dye.extinction_coefficient = row.extinction_coefficient;
        dye.quantum_yield = row.quantum_yield;

        std::map<std::string, std::vector<SpectrumPoint> >::iterator it =
            c.curves.find(row.name);
        if (it != c.curves.end() && !it->second.empty()) {
            // The rows are not required to arrive in wavelength order, and the
            // trapezoidal integral below assumes they are.
            std::sort(it->second.begin(), it->second.end());
            for (std::size_t p = 0; p < it->second.size(); ++p) {
                dye.spectrum.wavelength.push_back(it->second[p].wavelength);
                dye.spectrum.excitation.push_back(it->second[p].excitation);
                dye.spectrum.emission.push_back(it->second[p].emission);
            }
        }
        out[row.name] = dye;
    }
    return out;
}

//! What a dye template carries beyond the component topology.
struct MetaCtx {
    std::string center_atom, dipole_atom_1, dipole_atom_2;
    std::vector<std::string> positive_atoms, negative_atoms;
    ihm_keyword *m_key, *m_value;
};

//! Split a comma-separated atom list, trimming each name.
std::vector<std::string> split_atoms(const std::string& value) {
    std::vector<std::string> out;
    std::size_t start = 0;
    while (start <= value.size()) {
        const std::size_t comma = value.find(',', start);
        const std::size_t end = comma == std::string::npos ? value.size() : comma;
        std::string item = value.substr(start, end - start);
        const std::size_t a = item.find_first_not_of(" \t\n\r");
        const std::size_t b = item.find_last_not_of(" \t\n\r");
        if (a != std::string::npos) out.push_back(item.substr(a, b - a + 1));
        if (comma == std::string::npos) break;
        start = comma + 1;
    }
    return out;
}

void on_metadata(ihm_reader*, int, void* d, ihm_error**) {
    MetaCtx* c = (MetaCtx*) d;
    const std::string key = txt(c->m_key);
    const std::string value = txt(c->m_value);
    if (key.empty() || value.empty()) return;
    if (key == "center_atom") c->center_atom = value;
    else if (key == "dipole_atom_1") c->dipole_atom_1 = value;
    else if (key == "dipole_atom_2") c->dipole_atom_2 = value;
    else if (key == "positive_atoms") {
        const std::vector<std::string> a = split_atoms(value);
        c->positive_atoms.insert(c->positive_atoms.end(), a.begin(), a.end());
    } else if (key == "negative_atoms") {
        const std::vector<std::string> a = split_atoms(value);
        c->negative_atoms.insert(c->negative_atoms.end(), a.begin(), a.end());
    }
}

//! Fold a dye template's `_cgprobe_metadata` into a species.
void apply_template_metadata(const std::string& path, Probe& dye) {
    const int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) IMP_THROW("cannot open " << path, IOException);

    MetaCtx c;
    ihm_file* fh = ihm_file_new_from_fd(fd);
    ihm_reader* reader = ihm_reader_new(fh, false);
    ihm_category* cat = ihm_category_new(reader, "_cgprobe_metadata", on_metadata,
                                         NULL, NULL, &c, NULL);
    c.m_key = ihm_keyword_str_new(cat, "key");
    c.m_value = ihm_keyword_str_new(cat, "value");

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

    dye.dipole_atoms.clear();
    if (!c.dipole_atom_1.empty() && !c.dipole_atom_2.empty()) {
        dye.dipole_atoms.push_back(c.dipole_atom_1);
        dye.dipole_atoms.push_back(c.dipole_atom_2);
    }
    dye.chromophore_center_atom = c.center_atom;
    dye.positive_atoms = c.positive_atoms;
    dye.negative_atoms = c.negative_atoms;
}

}  // namespace probe_cif

Probe::Probe()
    : probe_type(PROBE_UNSPECIFIED),
      extinction_coefficient(probe_cif::nan_value()),
      quantum_yield(probe_cif::nan_value()),
      lifetime(probe_cif::nan_value()),
      radius(probe_cif::nan_value()),
      hydrodynamic_radius(probe_cif::nan_value()) {}

void Probe::show(std::ostream& out) const {
    out << "Probe(name=" << name;
    if (quantum_yield == quantum_yield) out << ", QY=" << quantum_yield;
    if (extinction_coefficient == extinction_coefficient)
        out << ", eps=" << extinction_coefficient;
    out << ", spectrum=" << (has_spectrum() ? "yes" : "no") << ")";
}

std::map<std::string, Probe> read_probe_library(std::string path) {
    if (path.empty()) {
        path = IMP::bff::get_data_path("rotamer_library") + "/R0/" +
               std::string(PROBE_LIBRARY_CIF);
    }

    // A missing file is not cached: stat fails, and the parse below raises the
    // error the caller is entitled to see on every call rather than once.
    struct stat info;
    if (stat(path.c_str(), &info) == 0) {
        probe_cif::CacheKey key;
        key.path = path;
#ifdef __APPLE__
        key.mtime_ns = (long long) info.st_mtimespec.tv_sec * 1000000000LL +
                       info.st_mtimespec.tv_nsec;
#else
        key.mtime_ns = (long long) info.st_mtim.tv_sec * 1000000000LL +
                       info.st_mtim.tv_nsec;
#endif
        key.size = (long long) info.st_size;

        std::map<probe_cif::CacheKey, std::map<std::string, Probe> >& cache = probe_cif::library_cache();
        std::map<probe_cif::CacheKey, std::map<std::string, Probe> >::const_iterator hit =
            cache.find(key);
        if (hit != cache.end()) return hit->second;

        std::map<std::string, Probe> parsed = probe_cif::parse_probe_library(path);
        cache[key] = parsed;
        return parsed;
    }
    return probe_cif::parse_probe_library(path);
}

unsigned int probe_library_cache_size() {
    return static_cast<unsigned int>(probe_cif::library_cache().size());
}

double spectral_overlap(const Probe& donor, const Probe& acceptor) {
    const std::vector<double>& wavelengths = donor.spectrum.wavelength;
    const std::vector<double>& emission = donor.spectrum.emission;
    const std::vector<double>& excitation = acceptor.spectrum.excitation;
    if (excitation.size() != wavelengths.size()) {
        IMP_THROW("donor and acceptor spectra must share one wavelength grid",
                  ValueException);
    }
    const double emission_integral = probe_cif::trapezoid(emission, wavelengths);
    if (emission_integral == 0.0) return 0.0;

    std::vector<double> weighted(wavelengths.size());
    for (std::size_t i = 0; i < wavelengths.size(); ++i) {
        const double l2 = wavelengths[i] * wavelengths[i];
        weighted[i] = emission[i] * acceptor.extinction_coefficient *
                      excitation[i] * l2 * l2;
    }
    return probe_cif::trapezoid(weighted, wavelengths) / emission_integral;
}

double forster_radius(const Probe& donor, const Probe& acceptor, double k2,
                      double refractive_index) {
    // A spin label is not a failed dye; it is a probe of a kind that absorbs
    // and emits nothing, and "needs a spectrum" is a confusing thing to tell
    // someone who asked for the Forster radius of a nitroxide.
    if (!probe_is_fluorescent(donor.probe_type) ||
        !probe_is_fluorescent(acceptor.probe_type)) {
        IMP_THROW("a Forster radius needs two fluorescent probes; "
                      << donor.name << " and " << acceptor.name
                      << " are not both (a spin label has no spectra)",
                  ValueException);
    }
    if (!donor.has_spectrum() || !acceptor.has_spectrum()) {
        IMP_THROW("both probes need a spectrum to derive R0 ("
                      << donor.name << ": " << donor.has_spectrum() << ", "
                      << acceptor.name << ": " << acceptor.has_spectrum() << ")",
                  ValueException);
    }
    // NaN is the "not in the library" value, and `x != x` is the one test that
    // catches it without <cmath>'s isnan being ambiguous under -ffast-math.
    if (donor.quantum_yield != donor.quantum_yield) {
        IMP_THROW(donor.name << " has no quantum yield", ValueException);
    }
    if (acceptor.extinction_coefficient != acceptor.extinction_coefficient) {
        IMP_THROW(acceptor.name << " has no extinction coefficient",
                  ValueException);
    }
    const double n4 = refractive_index * refractive_index * refractive_index *
                      refractive_index;
    const double overlap = spectral_overlap(donor, acceptor);
    return probe_cif::NM_TO_ANGSTROM * probe_cif::R0_FACTOR *
           std::pow(k2 * donor.quantum_yield / n4 * overlap, 1.0 / 6.0);
}

std::vector<std::string> normalize_probe_name(std::string probe_name) {
    // Trim, then split on the last run of whitespace: the trailing token is the
    // number ("488", "Thio12") and everything before it is the type.
    const std::string space = " \t\n\r\f\v";
    const std::size_t first = probe_name.find_first_not_of(space);
    const std::size_t last = probe_name.find_last_not_of(space);
    const std::string name =
        first == std::string::npos ? std::string()
                                   : probe_name.substr(first, last - first + 1);

    std::string compact;
    for (std::size_t i = 0; i < name.size(); ++i) {
        if (space.find(name[i]) == std::string::npos) compact += name[i];
    }

    std::vector<std::string> out;
    const std::size_t split = name.find_last_of(space);
    if (split == std::string::npos) {
        out.push_back(compact);
        out.push_back(compact);
        out.push_back(compact);
        return out;
    }
    std::string type = name.substr(0, split);
    std::string number = name.substr(split + 1);
    // The number must be alphanumeric; anything else means the name does not
    // split.
    for (std::size_t i = 0; i < number.size(); ++i) {
        if (!std::isalnum(static_cast<unsigned char>(number[i]))) {
            out.push_back(compact);
            out.push_back(compact);
            out.push_back(compact);
            return out;
        }
    }
    const std::size_t type_end = type.find_last_not_of(space);
    type = type_end == std::string::npos ? std::string()
                                         : type.substr(0, type_end + 1);
    out.push_back(type);
    out.push_back(number);
    out.push_back(type + number);
    return out;
}

namespace {

//! Lowercase, alphanumeric only -- the form two spellings are compared in.
std::string dl_fold(const std::string& s) {
    std::string out;
    for (std::size_t i = 0; i < s.size(); ++i) {
        const unsigned char c = static_cast<unsigned char>(s[i]);
        if (std::isalnum(c)) out += static_cast<char>(std::tolower(c));
    }
    return out;
}

//! The vendor spellings a common abbreviation stands for.
std::vector<std::string> dl_expansions(const std::string& folded) {
    std::vector<std::string> out;
    out.push_back(folded);
    // "alexa488" is what everyone writes; "alexafluor488" is what the vendor
    // and therefore the library calls it.
    if (folded.compare(0, 5, "alexa") == 0 &&
        folded.compare(0, 10, "alexafluor") != 0) {
        out.push_back("alexafluor" + folded.substr(5));
    }
    // Cyanines are keyed by their supplier in the shipped library.
    if (folded.compare(0, 2, "cy") == 0 &&
        folded.compare(0, 11, "lumiprobecy") != 0) {
        out.push_back("lumiprobecy" + folded.substr(2));
    }
    return out;
}

}  // namespace

std::string probe_fold_name(const std::string& name) { return dl_fold(name); }

bool probe_name_matches(const std::string& folded, const std::string& key) {
    const std::vector<std::string> wanted = dl_expansions(folded);
    const std::string k = dl_fold(key);
    for (std::size_t i = 0; i < wanted.size(); ++i) {
        if (k == wanted[i]) return true;
    }
    return false;
}

std::string resolve_probe_name(std::string name, std::string library_cif) {
    const std::map<std::string, Probe> library = read_probe_library(library_cif);
    const std::string compact = normalize_probe_name(name)[2];

    // An exact key always wins, so nothing that worked before changes.
    if (library.count(compact)) return compact;

    const std::vector<std::string> wanted = dl_expansions(dl_fold(compact));
    std::vector<std::string> hits;
    for (std::map<std::string, Probe>::const_iterator it = library.begin();
         it != library.end(); ++it) {
        const std::string key = dl_fold(it->first);
        for (std::size_t w = 0; w < wanted.size(); ++w) {
            if (key == wanted[w]) { hits.push_back(it->first); break; }
        }
    }
    if (hits.size() == 1) return hits[0];
    if (hits.size() > 1) {
        std::string all;
        for (std::size_t i = 0; i < hits.size(); ++i) {
            all += (i ? ", " : "") + hits[i];
        }
        IMP_THROW("'" << name << "' matches more than one probe (" << all
                  << "); name it exactly rather than have one picked",
                  ValueException);
    }
    std::string known;
    std::map<std::string, Probe>::const_iterator it = library.begin();
    for (int i = 0; i < 3 && it != library.end(); ++i, ++it) {
        known += (i ? ", " : "") + it->first;
    }
    IMP_THROW("'" << name << "' is not in the probe library; known names "
                     "look like " << known, ValueException);
}

Probe find_probe(std::string name, std::string library_cif, std::string template_cif) {
    const std::map<std::string, Probe> library = read_probe_library(library_cif);
    std::map<std::string, Probe>::const_iterator hit =
            library.find(resolve_probe_name(name, library_cif));
    if (hit == library.end()) {
        IMP_THROW("'" << name << "' is not in the probe library", ValueException);
    }
    Probe dye = hit->second;
    if (!template_cif.empty()) probe_cif::apply_template_metadata(template_cif, dye);
    return dye;
}

std::vector<std::string> available_probes(std::string library_cif) {
    const std::map<std::string, Probe> library = read_probe_library(library_cif);
    std::vector<std::string> out;
    for (std::map<std::string, Probe>::const_iterator it = library.begin();
         it != library.end(); ++it) {
        out.push_back(it->first);   // std::map is already sorted by key
    }
    return out;
}

double forster_radius_from_spectra(std::string donor, std::string acceptor,
                                   double k2, std::string library_cif,
                                   double refractive_index) {
    const std::map<std::string, Probe> library = read_probe_library(library_cif);
    // Resolve the way find_probe does. This used to match on whitespace alone,
    // so the name-based route rejected `Alexa488` while the Probe-based route
    // accepted it -- one spelling working through one door and not the other
    // is worse than neither working.
    std::map<std::string, Probe>::const_iterator d = library.end();
    std::map<std::string, Probe>::const_iterator a = library.end();
    try {
        d = library.find(resolve_probe_name(donor, library_cif));
        a = library.find(resolve_probe_name(acceptor, library_cif));
    } catch (const IMP::ValueException&) {
        // Fall through to the shared message below, which names both dyes.
    }
    if (d == library.end() || a == library.end() ||
        !d->second.has_spectrum() || !a->second.has_spectrum()) {
        IMP_THROW("No spectra for '" << donor << "' / '" << acceptor << "'",
                  ValueException);
    }
    return forster_radius(d->second, a->second, k2, refractive_index);
}

std::map<std::string, std::string> probe_flrcif_items() {
    std::map<std::string, std::string> m;
    m["name"] = "_flr_probe_list.chromophore_name";
    m["reactive_probe_name"] = "_flr_probe_list.reactive_probe_name";
    m["probe_origin"] = "_flr_probe_list.probe_origin";
    m["probe_link_type"] = "_flr_probe_list.probe_link_type";
    m["chromophore_center_atom"] = "_flr_probe_descriptor.chromophore_center_atom";
    m["lifetime"] = "_flr_reference_measurement_lifetime.lifetime";
    // no flrCIF item exists for these bff-native categories
    return m;
}

std::map<std::string, std::string> forster_radius_flrcif_items() {
    std::map<std::string, std::string> m;
    m["forster_radius"] = "_flr_fret_forster_radius.forster_radius";
    m["k2"] = "_flr_fret_forster_radius.kappa_squared";
    m["refractive_index"] = "_flr_fret_forster_radius.index_of_refraction";
    m["donor"] = "_flr_fret_forster_radius.donor_probe_id";
    m["acceptor"] = "_flr_fret_forster_radius.acceptor_probe_id";
    // no dictionary in the stack has an item for spectral overlap
    return m;
}

IMPBFF_END_NAMESPACE
