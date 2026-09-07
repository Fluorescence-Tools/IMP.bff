

/**
 * \file ProbeLibrary.cpp
 * \brief The dye as a species, its .pto library, and its coarse-grained system.
 *
 * Sections in the order of IMP/bff/ProbeLibrary.h; each is marked with the
 * file it came from.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

// -------- from ProbeLibrary.cpp --------

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

Probe get_probe(std::string name, std::string library_cif, std::string template_cif) {
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
    // Resolve the way get_probe does. This used to match on whitespace alone,
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

// -------- from ProbeContainer.cpp --------
/**
 * (formerly ProbeContainer.cpp, now a section of this file)
 * \brief A dye library as one `.mmfdb.pto`.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/Pto.h>
#include <IMP/bff/bff_config.h>
#include <IMP/bff/internal/json.h>


#include <cstring>
#include <sstream>

IMPBFF_BEGIN_NAMESPACE

const char* const PROBE_PTO_PROBES = "probes.json";
const char* const PROBE_PTO_SPECTRA = "spectra.json";
const char* const PROBE_PTO_MODEL = "probe_model.json";

namespace {

const char* const DYE_README =
    "PTO.MFDB dye container\n"
    "======================\n"
    "\n"
    "An EBML document (RFC 8794), DocType \"pto\". Each payload is an\n"
    "AttachedFile with a FileName, a PtoKind, a PtoEncoding and its bytes.\n"
    "Nothing is compressed, encrypted, or stored outside the file.\n"
    "\n"
    "Objects:\n"
    "  README        this text\n"
    "  probes.json   one row per dye. Identity in flrCIF terms\n"
    "                (_flr_probe_list.chromophore_name and friends); scalars\n"
    "                as _mmfdb_optical_property triples of property_name,\n"
    "                property_value and unit. A property that was not\n"
    "                measured is ABSENT -- never zero, because a quantum\n"
    "                yield of 0 is a real claim about a dye.\n"
    "  spectra.json  one row per (dye, spectrum_type, wavelength) point,\n"
    "                _mmfdb_spectrum.spectrum_type being excitation or\n"
    "                emission. Wavelengths are nm; intensities are\n"
    "                normalised so the peak is 1, which is NOT unit area --\n"
    "                a spectral overlap integral has to renormalise.\n"
    "  probe_model.json  the container, artifact and operation tags.\n"
    "\n"
    "Every controlled value is a term from mmfdb_flr_ext.dic; the container\n"
    "tags record which version. A Forster radius computed from this file is\n"
    "traceable to the quantum yield and the two spectra it came from.\n";

//! A scalar, if the library actually carries it.
void dc_add_property(nlohmann::json& into, const std::string& name,
                     double value, const std::string& unit) {
    // NaN is how the library spells "not carried". Writing it as 0 would be a
    // measurement nobody made.
    if (value != value) return;
    // Declared but unchecked is how a vocabulary rots; check at the write.
    mfdb_check_term("_mmfdb_optical_property.property_name", name);
    if (!unit.empty()) mfdb_check_term("_mmfdb_optical_property.unit", unit);
    nlohmann::json p;
    p["_mmfdb_optical_property.property_name"] = name;
    p["_mmfdb_optical_property.property_value"] = value;
    if (!unit.empty()) p["_mmfdb_optical_property.unit"] = unit;
    into.push_back(p);
}

std::vector<MfdbColumn> dc_probe_columns() {
    std::vector<MfdbColumn> c;
    c.push_back(MfdbColumn("_flr_probe_list.chromophore_name", "",
                           "_flr_probe_list.chromophore_name",
                           "The probe, as the library keys it."));
    c.push_back(MfdbColumn("_mmfdb_probe.probe_type", "",
                           "_mmfdb_probe.probe_type",
                           "Dye, fluorescent protein or spin label. Absent "
                           "means unspecified, which is what a container "
                           "written before this column says."));
    c.push_back(MfdbColumn("_mmfdb_probe.vendor", "",
                           "_mmfdb_probe.vendor",
                           "Who supplies it. Not an identifier."));
    c.push_back(MfdbColumn("_mmfdb_probe.dipole_atoms", "",
                           "_mmfdb_probe.dipole_atoms",
                           "The two atoms whose separation is the transition "
                           "dipole. Absent means isotropic, which is a claim, "
                           "not missing data."));
    c.push_back(MfdbColumn("_mmfdb_probe.positive_atoms", "",
                           "_mmfdb_probe.positive_atoms",
                           "Atoms carrying a formal positive charge."));
    c.push_back(MfdbColumn("_mmfdb_probe.negative_atoms", "",
                           "_mmfdb_probe.negative_atoms",
                           "Atoms carrying a formal negative charge."));
    c.push_back(MfdbColumn("_flr_probe_list.probe_origin", "",
                           "_flr_probe_list.probe_origin",
                           "Intrinsic or extrinsic."));
    c.push_back(MfdbColumn("_flr_probe_list.probe_link_type", "",
                           "_flr_probe_list.probe_link_type",
                           "How it attaches."));
    c.push_back(MfdbColumn("_flr_probe_list.reactive_probe_name", "",
                           "_flr_probe_list.reactive_probe_name",
                           "The reagent form."));
    c.push_back(MfdbColumn("_flr_probe_descriptor.chromophore_center_atom", "",
                           "_flr_probe_descriptor.chromophore_center_atom",
                           "Which atom the chromophore is centred on."));
    c.push_back(MfdbColumn("optical_properties", "",
                           "_mmfdb_optical_property.property_name",
                           "The dye's scalars, as name/value/unit triples. "
                           "Absent means not measured."));
    return c;
}

std::vector<MfdbColumn> dc_spectrum_columns() {
    std::vector<MfdbColumn> c;
    c.push_back(MfdbColumn("_mmfdb_spectrum.probe_id", "",
                           "_mmfdb_spectrum.probe_id",
                           "The dye, by chromophore name."));
    c.push_back(MfdbColumn("_mmfdb_spectrum.spectrum_type", "",
                           "_mmfdb_spectrum.spectrum_type",
                           "excitation or emission."));
    c.push_back(MfdbColumn("_mmfdb_spectrum.wavelengths", "nanometres",
                           "_mmfdb_spectrum.wavelengths",
                           "The wavelength axis, as an array."));
    c.push_back(MfdbColumn("_mmfdb_spectrum.intensity_values", "dimensionless",
                           "_mmfdb_spectrum.intensity_values",
                           "Paired elementwise with the wavelengths. Normalised so the peak is 1; not unit area."));
    return c;
}

//! Find a dye in a container by exact key, then by the library's fold rules.
/*! Duplicating the fold would be worse than sharing it, so this asks
    resolve_probe_name for the canonical spelling of each key in *this*
    container and matches the query against that. */
std::map<std::string, Probe>::const_iterator dc_lookup(
        const std::map<std::string, Probe>& dyes, const std::string& name) {
    std::map<std::string, Probe>::const_iterator hit = dyes.find(name);
    if (hit != dyes.end()) return hit;

    // Build a one-entry library view so the shared resolver does the work.
    std::vector<std::string> keys;
    for (std::map<std::string, Probe>::const_iterator it = dyes.begin();
         it != dyes.end(); ++it) {
        keys.push_back(it->first);
    }
    const std::string folded = probe_fold_name(name);
    std::map<std::string, Probe>::const_iterator found = dyes.end();
    int matches = 0;
    for (std::map<std::string, Probe>::const_iterator it = dyes.begin();
         it != dyes.end(); ++it) {
        if (probe_name_matches(folded, it->first)) { found = it; ++matches; }
    }
    return matches == 1 ? found : dyes.end();
}

std::string dc_object_text(const std::string& path, const std::string& name) {
    PtoReader r(path);
    const int i = r.find(name);
    if (i < 0) return "";
    const std::vector<unsigned char> d = r.data(r.objects()[i]);
    return std::string(d.begin(), d.end());
}

nlohmann::json dc_tags(const std::vector<MfdbTag>& tags) {
    nlohmann::json j = nlohmann::json::object();
    for (std::size_t i = 0; i < tags.size(); ++i) j[tags[i].item] = tags[i].value;
    return j;
}

}  // namespace

// ---------------------------------------------------------------------------
// Writing
// ---------------------------------------------------------------------------

void probe_write_pto(const std::string& path,
                   const std::map<std::string, Probe>& dyes,
                   const std::string& source) {
    nlohmann::json probes = nlohmann::json::array();
    nlohmann::json spectra = nlohmann::json::array();

    for (std::map<std::string, Probe>::const_iterator it = dyes.begin();
         it != dyes.end(); ++it) {
        const Probe& d = it->second;

        nlohmann::json row;
        row["_flr_probe_list.chromophore_name"] = d.name.empty() ? it->first
                                                                 : d.name;
        // Written even when unspecified: "the library did not say" is a fact
        // about the source, and leaving the column out makes it unreadable
        // from a reader that would otherwise assume a default.
        const std::string type = probe_type_to_string(d.probe_type);
        mfdb_check_term("_mmfdb_probe.probe_type", type);
        row["_mmfdb_probe.probe_type"] = type;
        if (!d.vendor.empty()) row["_mmfdb_probe.vendor"] = d.vendor;
        if (!d.probe_origin.empty())
            row["_flr_probe_list.probe_origin"] = d.probe_origin;
        if (!d.probe_link_type.empty())
            row["_flr_probe_list.probe_link_type"] = d.probe_link_type;
        if (!d.reactive_probe_name.empty())
            row["_flr_probe_list.reactive_probe_name"] = d.reactive_probe_name;
        if (!d.chromophore_center_atom.empty())
            row["_flr_probe_descriptor.chromophore_center_atom"] =
                    d.chromophore_center_atom;
        // Three lists that are easy to drop on the floor. Losing
        // `dipole_atoms` is the one that does damage without looking like it:
        // an anisotropic probe comes back isotropic, kappa^2 quietly falls
        // back to 2/3, and the R0 that follows is wrong by a factor no reader
        // can see. The charges decide how a probe behaves beside a charged
        // surface. None of the three had an flrCIF item, which is why they
        // were left out; they have `_mmfdb_probe.*` ones now.
        if (!d.dipole_atoms.empty())
            row["_mmfdb_probe.dipole_atoms"] = d.dipole_atoms;
        if (!d.positive_atoms.empty())
            row["_mmfdb_probe.positive_atoms"] = d.positive_atoms;
        if (!d.negative_atoms.empty())
            row["_mmfdb_probe.negative_atoms"] = d.negative_atoms;

        nlohmann::json props = nlohmann::json::array();
        dc_add_property(props, "quantum_yield", d.quantum_yield,
                        "dimensionless");
        dc_add_property(props, "extinction_coefficient",
                        d.extinction_coefficient, "per_molar_per_centimetre");
        dc_add_property(props, "fluorescence_lifetime", d.lifetime,
                        "nanoseconds");
        dc_add_property(props, "steric_radius", d.radius, "angstroms");
        dc_add_property(props, "hydrodynamic_radius", d.hydrodynamic_radius,
                        "angstroms");
        row["optical_properties"] = props;
        probes.push_back(row);

        // One row per (probe, spectrum_type), carrying arrays -- which is what
        // `_mmfdb_spectrum.wavelengths` and `.intensity_values` are declared
        // as, what the `spectra` table stores, and about a tenth the size of
        // the same data written a point per row (9.4 MB against 0.9 MB for
        // this library, almost all of it repeated JSON keys).
        const Spectrum& s = d.spectrum;
        const std::string key = d.name.empty() ? it->first : d.name;
        for (int which = 0; which < 2; ++which) {
            const std::vector<double>& y = which == 0 ? s.excitation
                                                      : s.emission;
            if (y.empty()) continue;
            const char* const stype = which == 0 ? "excitation" : "emission";
            mfdb_check_term("_mmfdb_spectrum.spectrum_type", stype);

            nlohmann::json wl = nlohmann::json::array();
            nlohmann::json iv = nlohmann::json::array();
            for (std::size_t k = 0; k < s.wavelength.size() && k < y.size();
                 ++k) {
                wl.push_back(s.wavelength[k]);
                iv.push_back(y[k]);
            }
            nlohmann::json p;
            p["_mmfdb_spectrum.probe_id"] = key;
            p["_mmfdb_spectrum.spectrum_type"] = stype;
            p["_mmfdb_spectrum.wavelengths"] = wl;
            p["_mmfdb_spectrum.intensity_values"] = iv;
            p["_mmfdb_spectrum.wavelength_unit"] = "nm";
            p["_mmfdb_spectrum.intensity_unit"] = "normalized";
            spectra.push_back(p);
        }
    }

    nlohmann::json probe_doc, spectrum_doc;
    probe_doc["rows"] = probes;
    spectrum_doc["rows"] = spectra;
    const std::string probe_bytes = probe_doc.dump(1);
    const std::string spectrum_bytes = spectrum_doc.dump(1);

    nlohmann::json settings;
    settings["dyes"] = static_cast<int>(dyes.size());
    settings["spectrum_normalisation"] = "peak";
    settings["wavelength_unit"] = "nm";
    const std::string settings_json = settings.dump();

    nlohmann::json model;
    model["_container"] = dc_tags(mfdb_container_tags());
    model["_operation"] = dc_tags(mfdb_operation_tags(
            "import", "probe_library", settings_json, "IMP.bff",
            get_module_version()));
    model["settings"] = settings;

    nlohmann::json artifacts = nlohmann::json::object();
    artifacts[PROBE_PTO_PROBES] = dc_tags(mfdb_artifact_tags(
            PROBE_PTO_PROBES, "parameter_table", "json", "species",
            static_cast<long>(probes.size()), mfdb_checksum(probe_bytes),
            dc_probe_columns()));
    artifacts[PROBE_PTO_SPECTRA] = dc_tags(mfdb_artifact_tags(
            PROBE_PTO_SPECTRA, "parameter_table", "json", "spectrum",
            static_cast<long>(spectra.size()), mfdb_checksum(spectrum_bytes),
            dc_spectrum_columns()));
    model["_artifacts"] = artifacts;

    // A spectrum belongs to a probe, and the join is by name rather than by
    // row order -- the two tables count different things.
    nlohmann::json edges = nlohmann::json::array();
    edges.push_back(dc_tags(mfdb_edge_tags(
            PROBE_PTO_PROBES, PROBE_PTO_SPECTRA, "maps_rows_of",
            "_flr_probe_list.chromophore_name", "_mmfdb_spectrum.probe_id")));
    model["_edges"] = edges;

    if (!source.empty()) {
        MfdbAttribution a;
        a.source = source;
        // The bundled tables come from FRETpredict's compilation of vendor and
        // FPbase data; stating terms is required and guessing them is not
        // allowed, so this says what is actually known.
        a.terms = "as published by the upstream compilation; see source";
        model["_attribution"] = dc_tags(mfdb_attribution_tags(a));
    }
    const std::string model_bytes = model.dump(1);

    PtoWriter w(path);
    w.add("README", "readme", "text", DYE_README, std::strlen(DYE_README));
    w.add(PROBE_PTO_PROBES, "mfdb.probes", "json", probe_bytes.data(),
          probe_bytes.size());
    w.add(PROBE_PTO_SPECTRA, "mfdb.spectra", "json", spectrum_bytes.data(),
          spectrum_bytes.size());
    w.add(PROBE_PTO_MODEL, "mfdb.model", "json", model_bytes.data(),
          model_bytes.size());
    w.close();
}

// ---------------------------------------------------------------------------
// Reading
// ---------------------------------------------------------------------------

std::map<std::string, Probe> probe_read_pto(const std::string& path) {
    const std::string probe_text = dc_object_text(path, PROBE_PTO_PROBES);
    if (probe_text.empty()) {
        IMP_THROW("probe_read_pto: " << path << " carries no " << PROBE_PTO_PROBES,
                  IOException);
    }
    std::map<std::string, Probe> out;
    const nlohmann::json probes = nlohmann::json::parse(probe_text);
    const nlohmann::json& rows = probes.at("rows");
    for (nlohmann::json::const_iterator it = rows.begin(); it != rows.end();
         ++it) {
        Probe d;
        d.name = it->at("_flr_probe_list.chromophore_name").get<std::string>();
        if (it->count("_mmfdb_probe.probe_type"))
            d.probe_type = probe_type_from_string(
                (*it)["_mmfdb_probe.probe_type"].get<std::string>());
        if (it->count("_mmfdb_probe.vendor"))
            d.vendor = (*it)["_mmfdb_probe.vendor"].get<std::string>();
        if (it->count("_mmfdb_probe.dipole_atoms"))
            d.dipole_atoms = (*it)["_mmfdb_probe.dipole_atoms"]
                    .get<std::vector<std::string> >();
        if (it->count("_mmfdb_probe.positive_atoms"))
            d.positive_atoms = (*it)["_mmfdb_probe.positive_atoms"]
                    .get<std::vector<std::string> >();
        if (it->count("_mmfdb_probe.negative_atoms"))
            d.negative_atoms = (*it)["_mmfdb_probe.negative_atoms"]
                    .get<std::vector<std::string> >();
        if (it->count("_flr_probe_list.probe_origin"))
            d.probe_origin = (*it)["_flr_probe_list.probe_origin"].get<std::string>();
        if (it->count("_flr_probe_list.probe_link_type"))
            d.probe_link_type =
                    (*it)["_flr_probe_list.probe_link_type"].get<std::string>();
        if (it->count("_flr_probe_list.reactive_probe_name"))
            d.reactive_probe_name =
                    (*it)["_flr_probe_list.reactive_probe_name"].get<std::string>();
        if (it->count("_flr_probe_descriptor.chromophore_center_atom"))
            d.chromophore_center_atom =
                    (*it)["_flr_probe_descriptor.chromophore_center_atom"]
                            .get<std::string>();

        if (it->count("optical_properties")) {
            const nlohmann::json& props = (*it)["optical_properties"];
            for (nlohmann::json::const_iterator p = props.begin();
                 p != props.end(); ++p) {
                const std::string n =
                        p->at("_mmfdb_optical_property.property_name")
                                .get<std::string>();
                if (!p->count("_mmfdb_optical_property.property_value")) continue;
                const double v =
                        (*p)["_mmfdb_optical_property.property_value"].get<double>();
                if (n == "quantum_yield") d.quantum_yield = v;
                else if (n == "extinction_coefficient") d.extinction_coefficient = v;
                else if (n == "fluorescence_lifetime") d.lifetime = v;
                else if (n == "steric_radius") d.radius = v;
                else if (n == "hydrodynamic_radius") d.hydrodynamic_radius = v;
            }
        }
        out[d.name] = d;
    }

    const std::string spectrum_text = dc_object_text(path, PROBE_PTO_SPECTRA);
    if (!spectrum_text.empty()) {
        const nlohmann::json spectra = nlohmann::json::parse(spectrum_text);
        const nlohmann::json& srows = spectra.at("rows");
        // Two rows per dye -- excitation and emission -- sharing one axis.
        std::map<std::string, std::map<double, std::pair<double, double> > > acc;
        for (nlohmann::json::const_iterator it = srows.begin();
             it != srows.end(); ++it) {
            const std::string key =
                    it->at("_mmfdb_spectrum.probe_id").get<std::string>();
            const std::string type =
                    it->at("_mmfdb_spectrum.spectrum_type").get<std::string>();
            const nlohmann::json& wl = it->at("_mmfdb_spectrum.wavelengths");
            const nlohmann::json& iv =
                    it->at("_mmfdb_spectrum.intensity_values");
            // A container written before the array shape carries one scalar
            // per row here. Say so, rather than letting the JSON library
            // report `operator[] with a numeric argument` from three frames
            // down -- which is what it did, and it cost a while to place.
            if (!wl.is_array() || !iv.is_array()) {
                IMP_THROW("probe_read_pto: " << path << " stores a scalar where "
                          << "_mmfdb_spectrum.wavelengths must be an array. "
                          << "It predates the array shape; rewrite it with "
                          << "probe_library_to_pto.", IOException);
            }
            const std::size_t n = std::min(wl.size(), iv.size());
            for (std::size_t k = 0; k < n; ++k) {
                std::pair<double, double>& slot = acc[key][wl[k].get<double>()];
                if (type == "emission") slot.second = iv[k].get<double>();
                else slot.first = iv[k].get<double>();
            }
        }
        for (std::map<std::string,
                      std::map<double, std::pair<double, double> > >::const_iterator
                     a = acc.begin(); a != acc.end(); ++a) {
            std::map<std::string, Probe>::iterator d = out.find(a->first);
            if (d == out.end()) continue;
            for (std::map<double, std::pair<double, double> >::const_iterator
                         p = a->second.begin(); p != a->second.end(); ++p) {
                d->second.spectrum.wavelength.push_back(p->first);
                d->second.spectrum.excitation.push_back(p->second.first);
                d->second.spectrum.emission.push_back(p->second.second);
            }
        }
    }
    return out;
}

int probe_library_to_pto(const std::string& path,
                       const std::string& library_cif) {
    const std::map<std::string, Probe> dyes = read_probe_library(library_cif);
    probe_write_pto(path, dyes,
                  library_cif.empty() ? std::string("probe_library.cif")
                                      : library_cif);
    return static_cast<int>(dyes.size());
}

double probe_pto_forster_radius(const std::string& path,
                              const std::string& donor,
                              const std::string& acceptor, double kappa2,
                              double refractive_index) {
    const std::map<std::string, Probe> dyes = probe_read_pto(path);

    // Resolution has to work against *this* container rather than the shipped
    // library, so the same folding rules are applied to its keys: a caller
    // naming `Alexa488` must reach `AlexaFluor488` here exactly as it does
    // through get_probe, or the container is harder to use than the CIF it
    // replaces.
    std::map<std::string, Probe>::const_iterator d = dc_lookup(dyes, donor);
    std::map<std::string, Probe>::const_iterator a = dc_lookup(dyes, acceptor);
    if (d == dyes.end() || a == dyes.end()) {
        std::string names;
        int shown = 0;
        for (std::map<std::string, Probe>::const_iterator it = dyes.begin();
             it != dyes.end() && shown < 3; ++it, ++shown) {
            names += (shown ? ", " : "") + it->first;
        }
        IMP_THROW("probe_pto_forster_radius: '"
                  << (d == dyes.end() ? donor : acceptor)
                  << "' is not in " << path << "; names look like " << names,
                  ValueException);
    }
    // Angstrom, from `forster_radius`, which converts once for the whole
    // package (2026-08-25). This used to multiply by ten here, as five other
    // call sites did.
    return forster_radius(d->second, a->second, kappa2, refractive_index);
}

IMPBFF_END_NAMESPACE

// -------- from ProbeForceField.cpp --------
/**
 * (formerly ProbeForceField.cpp, now a section of this file)
 * \brief A coarse-grained dye/protein system: sites, bonds, and their terms.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/MolecularGraph.h>


#include <set>

IMPBFF_BEGIN_NAMESPACE

namespace {

double ffjson_num(const nlohmann::json& j, const char* k, double def) {
    if (j.contains(k) && j[k].is_number()) return j[k].get<double>();
    return def;
}

int ffjson_int(const nlohmann::json& j, const char* k, int def) {
    if (j.contains(k) && j[k].is_number_integer()) return j[k].get<int>();
    return def;
}

//! The first non-empty string of two keys.
std::string ffjson_pick(const nlohmann::json& j, const char* a, const char* b) {
    for (const char* k : {a, b}) {
        if (j.contains(k) && j[k].is_string() && !j[k].get<std::string>().empty())
            return j[k].get<std::string>();
    }
    return std::string();
}

//! A value as a string: strings verbatim, numbers by their JSON spelling,
//! null as empty.
std::string ffjson_val(const nlohmann::json& v) {
    if (v.is_string()) return v.get<std::string>();
    if (v.is_null()) return std::string();
    return v.dump();
}

//! A string field, whatever the JSON type: a numeric id (`"id": 0`) is the
//! string "0", not an error.
std::string ffjson_str(const nlohmann::json& j, const char* k) {
    return j.contains(k) ? ffjson_val(j[k]) : std::string();
}

}  // namespace

ProbeForceFieldSystem forcefield_system_from_json(const std::string& json) {
    nlohmann::json d = nlohmann::json::parse(json);
    if (!d.is_object()) {
        IMP_THROW("forcefield_system_from_json needs a JSON object",
                  ValueException);
    }
    ProbeForceFieldSystem out(ffjson_str(d, "name"));

    std::map<std::string, FFComponent> comps;
    if (d.contains("components") && d["components"].is_object()) {
        for (auto it = d["components"].begin(); it != d["components"].end();
             ++it) {
            FFComponent c;
            c.mol2_path = ffjson_pick(it.value(), "mol2_path", "mol2");
            c.pdb_path = ffjson_pick(it.value(), "pdb_path", "pdb");
            c.role = ffjson_pick(it.value(), "role", "");
            comps[it.key()] = c;
        }
    }
    out.set_components(comps);

    const nlohmann::json sites_in =
            d.contains("sites") && d["sites"].is_array() ? d["sites"]
                                                         : nlohmann::json::array();
    // A group member or bonded term may name a site by its number rather than
    // its id; resolve the ones the sites section declares.
    std::map<int, std::string> by_no;
    for (int i = 0; i < static_cast<int>(sites_in.size()); ++i) {
        const std::string id = ffjson_str(sites_in[i], "id");
        by_no[i] = id;
        if (sites_in[i].contains("site_no") &&
            sites_in[i]["site_no"].is_number_integer()) {
            by_no[sites_in[i]["site_no"].get<int>()] = id;
        }
    }
    auto token = [&](const nlohmann::json& x) -> std::string {
        if (x.is_string()) return x.get<std::string>();
        if (x.is_number_integer()) {
            auto f = by_no.find(x.get<int>());
            return f != by_no.end() ? f->second : x.dump();
        }
        return x.dump();
    };

    std::vector<FFSite> sites;
    for (int i = 0; i < static_cast<int>(sites_in.size()); ++i) {
        const nlohmann::json& row = sites_in[i];
        FFSite st;
        st.id = ffjson_str(row, "id");
        st.component = ffjson_str(row, "component");
        st.atom_name = ffjson_str(row, "atom_name");
        st.element = ffjson_str(row, "element");
        st.site_no = row.contains("site_no") &&
                             row["site_no"].is_number_integer()
                             ? row["site_no"].get<int>()
                             : i;
        st.site_serial = ffjson_int(row, "site_serial", 0);
        st.radius = ffjson_num(row, "radius", 1.7);
        st.mass = ffjson_num(row, "mass", 12.0);
        sites.push_back(st);
    }
    out.set_sites(sites);

    auto groups_from = [&](const char* key) {
        std::map<std::string, std::vector<std::string> > m;
        if (d.contains(key) && d[key].is_object()) {
            for (auto it = d[key].begin(); it != d[key].end(); ++it) {
                std::vector<std::string> v;
                if (it.value().is_array()) {
                    for (const auto& x : it.value()) v.push_back(token(x));
                }
                m[it.key()] = v;
            }
        }
        return m;
    };
    out.set_groups(groups_from("groups"));
    out.set_rb_groups(groups_from("rb_groups"));
    out.set_md_fixed_groups(groups_from("md_fixed_groups"));

    std::vector<std::string> fixed;
    if (d.contains("fixed_groups") && d["fixed_groups"].is_array()) {
        for (const auto& g : d["fixed_groups"]) {
            fixed.push_back(g.is_string() ? g.get<std::string>() : g.dump());
        }
    }
    out.set_fixed_groups(fixed);

    auto scalar_types = [&](const char* key) {
        std::map<std::string, double> m;
        if (d.contains(key) && d[key].is_object()) {
            for (auto it = d[key].begin(); it != d[key].end(); ++it) {
                m[it.key()] = ffjson_num(it.value(), "k", 0.0);
            }
        }
        return m;
    };
    out.set_bond_types(scalar_types("bond_types"));
    out.set_angle_types(scalar_types("angle_types"));

    auto torsion_types = [&](const char* key) {
        std::map<std::string, FFTorsionType> m;
        if (d.contains(key) && d[key].is_object()) {
            for (auto it = d[key].begin(); it != d[key].end(); ++it) {
                FFTorsionType t;
                t.periodicity = ffjson_int(it.value(), "periodicity", 1);
                t.phase = ffjson_num(it.value(), "phase_rad", 0.0);
                t.k = ffjson_num(it.value(), "k", 0.0);
                m[it.key()] = t;
            }
        }
        return m;
    };
    out.set_torsion_types(torsion_types("torsion_types"));
    out.set_improper_types(torsion_types("improper_types"));

    std::map<std::string, FFLJType> lj;
    if (d.contains("lj_types") && d["lj_types"].is_object()) {
        for (auto it = d["lj_types"].begin(); it != d["lj_types"].end(); ++it) {
            FFLJType t;
            t.element = ffjson_str(it.value(), "element");
            t.rmin_half = ffjson_num(it.value(), "rmin_half", 0.0);
            t.epsilon = ffjson_num(it.value(), "epsilon", 0.0);
            lj[it.key()] = t;
        }
    }
    out.set_lj_types(lj);

    if (d.contains("bonds") && d["bonds"].is_array()) {
        std::vector<FFBond> rows;
        for (const auto& row : d["bonds"]) {
            if (!row.is_array() || row.size() != 4) {
                IMP_THROW("a bond is [site_a, site_b, length, type_id]",
                          ValueException);
            }
            FFBond r;
            r.site_a = token(row[0]);
            r.site_b = token(row[1]);
            r.length = row[2].is_number() ? row[2].get<double>() : 0.0;
            r.type_id = ffjson_val(row[3]);
            rows.push_back(r);
        }
        out.set_bonds(rows);
    }
    if (d.contains("angles") && d["angles"].is_array()) {
        std::vector<FFAngle> rows;
        for (const auto& row : d["angles"]) {
            if (!row.is_array() || row.size() != 5) {
                IMP_THROW("an angle is [site_a, site_b, site_c, theta, type_id]",
                          ValueException);
            }
            FFAngle r;
            r.site_a = token(row[0]);
            r.site_b = token(row[1]);
            r.site_c = token(row[2]);
            r.theta = row[3].is_number() ? row[3].get<double>() : 0.0;
            r.type_id = ffjson_val(row[4]);
            rows.push_back(r);
        }
        out.set_angles(rows);
    }
    auto torsions_from = [&](const char* key) {
        std::vector<FFTorsion> rows;
        if (d.contains(key) && d[key].is_array()) {
            for (const auto& row : d[key]) {
                if (!row.is_array() || row.size() != 5) {
                    IMP_THROW("a torsion is [site_a, site_b, site_c, site_d,"
                              " type_id]",
                              ValueException);
                }
                FFTorsion r;
                r.site_a = token(row[0]);
                r.site_b = token(row[1]);
                r.site_c = token(row[2]);
                r.site_d = token(row[3]);
                r.type_id = ffjson_val(row[4]);
                rows.push_back(r);
            }
        }
        return rows;
    };
    out.set_dihedrals(torsions_from("dihedrals"));
    out.set_impropers(torsions_from("impropers"));

    if (d.contains("probes") && d["probes"].is_array()) {
        std::vector<FFProbe> probes;
        for (const auto& row : d["probes"]) {
            FFProbe pr;
            pr.id = ffjson_int(row, "id", 0);
            pr.name = ffjson_str(row, "name");
            pr.origin = ffjson_str(row, "origin");
            pr.link_type = ffjson_str(row, "link_type");
            if (pr.origin.empty()) pr.origin = "extrinsic";
            if (pr.link_type.empty()) pr.link_type = "covalent";
            probes.push_back(pr);
        }
        out.set_probes(probes);
    }

    FFNonbonded nb;
    const nlohmann::json empty = nlohmann::json::object();
    const nlohmann::json& nbs =
            d.contains("nonbonded") && d["nonbonded"].is_object()
                    ? d["nonbonded"] : empty;
    nb.enabled = !(nbs.contains("enabled") && nbs["enabled"].is_boolean() &&
                   !nbs["enabled"].get<bool>());
    nb.k = ffjson_num(nbs, "k", 5.0);
    nb.cutoff = ffjson_num(nbs, "cutoff_A", 6.0);
    out.set_nonbonded(nb);

    const nlohmann::json& sp =
            d.contains("sampling") && d["sampling"].is_object()
                    ? d["sampling"] : empty;
    FFSampling s;
    s.temperature_K = ffjson_num(sp, "temperature_K", 300.0);
    s.friction_ps = ffjson_num(sp, "friction_ps", 10.0);
    s.timestep_fs = ffjson_num(sp, "timestep_fs", 0.25);
    s.n_steps = ffjson_int(sp, "n_steps", 500000);
    s.write_every = ffjson_int(sp, "write_every", 1000);
    s.minimize_steps = ffjson_int(sp, "minimize_steps", 200);
    out.set_sampling(s);
    return out;
}

std::vector<int> ProbeForceFieldSystem::get_component_sites(
        const std::string& component) const {
    std::vector<int> out;
    for (std::size_t i = 0; i < sites_.size(); ++i) {
        if (sites_[i].component == component) {
            out.push_back(static_cast<int>(i));
        }
    }
    return out;
}

namespace {
std::vector<std::string> with_role(
        const std::map<std::string, FFComponent>& components,
        const std::string& role) {
    std::vector<std::string> out;
    for (std::map<std::string, FFComponent>::const_iterator it =
                 components.begin(); it != components.end(); ++it) {
        if (it->second.role == role) out.push_back(it->first);
    }
    return out;
}
}  // namespace

std::vector<std::string> ProbeForceFieldSystem::get_fixed_components() const {
    return with_role(components_, "fixed");
}

std::vector<std::string> ProbeForceFieldSystem::get_mobile_components() const {
    return with_role(components_, "mobile");
}


// --------------------------------------------------------------------------
// The molecular graph
// --------------------------------------------------------------------------

std::map<std::string, std::vector<std::string> >
ProbeForceFieldSystem::get_bonded_neighbors() const {
    std::map<std::string, std::set<std::string> > adj;
    for (size_t i = 0; i < bonds_.size(); ++i) {
        adj[bonds_[i].site_a].insert(bonds_[i].site_b);
        adj[bonds_[i].site_b].insert(bonds_[i].site_a);
    }
    std::map<std::string, std::vector<std::string> > out;
    for (std::map<std::string, std::set<std::string> >::const_iterator it = adj.begin();
         it != adj.end(); ++it) {
        out[it->first] = std::vector<std::string>(it->second.begin(), it->second.end());
    }
    return out;
}

namespace {
std::pair<std::string, std::string> ordered_pair(const std::string& a,
                                                 const std::string& b) {
    return a <= b ? std::make_pair(a, b) : std::make_pair(b, a);
}
}

std::set<std::pair<std::string, std::string> >
ProbeForceFieldSystem::get_exclusions(bool include_impropers) const {
    std::set<std::pair<std::string, std::string> > excl;
    for (size_t i = 0; i < bonds_.size(); ++i)
        excl.insert(ordered_pair(bonds_[i].site_a, bonds_[i].site_b));
    for (size_t i = 0; i < angles_.size(); ++i)
        excl.insert(ordered_pair(angles_[i].site_a, angles_[i].site_c));
    for (size_t i = 0; i < dihedrals_.size(); ++i)
        excl.insert(ordered_pair(dihedrals_[i].site_a, dihedrals_[i].site_d));
    if (include_impropers) {
        for (size_t i = 0; i < impropers_.size(); ++i) {
            const FFTorsion& t = impropers_[i];
            excl.insert(ordered_pair(t.site_a, t.site_c));
            excl.insert(ordered_pair(t.site_a, t.site_d));
            excl.insert(ordered_pair(t.site_b, t.site_d));
        }
    }
    return excl;
}

std::vector<std::vector<std::string> >
ProbeForceFieldSystem::find_rings(int max_len) const {
    // through `MolecularGraph`, so the ring walk exists once: site ids are
    // numbered in sorted order and mapped back, which is what makes the
    // canonical form come out ascending by site id rather than by number.
    std::vector<std::string> ids;
    std::map<std::string, int> index;
    const std::map<std::string, std::vector<std::string> > adj = get_bonded_neighbors();
    for (std::map<std::string, std::vector<std::string> >::const_iterator it = adj.begin();
         it != adj.end(); ++it) {
        index[it->first] = (int)ids.size();
        ids.push_back(it->first);
    }
    std::vector<std::pair<int, int> > edges;
    for (size_t i = 0; i < bonds_.size(); ++i) {
        std::map<std::string, int>::const_iterator a = index.find(bonds_[i].site_a);
        std::map<std::string, int>::const_iterator b = index.find(bonds_[i].site_b);
        if (a != index.end() && b != index.end())
            edges.push_back(std::make_pair(a->second, b->second));
    }
    const MolecularGraph g(edges);
    const std::vector<std::vector<int> > rings = g.get_rings(max_len);
    std::vector<std::vector<std::string> > out;
    out.reserve(rings.size());
    for (size_t i = 0; i < rings.size(); ++i) {
        std::vector<std::string> r;
        r.reserve(rings[i].size());
        for (size_t j = 0; j < rings[i].size(); ++j) r.push_back(ids[rings[i][j]]);
        out.push_back(r);
    }
    return out;
}

bool ProbeForceFieldSystem::is_within_bonds(const std::string& a,
                                          const std::string& b,
                                          int max_depth) const {
    if (a == b) return true;
    const std::map<std::string, std::vector<std::string> > adj = get_bonded_neighbors();
    std::set<std::string> seen;
    seen.insert(a);
    std::vector<std::string> frontier(1, a);
    for (int depth = 0; depth < max_depth && !frontier.empty(); ++depth) {
        std::vector<std::string> next;
        for (size_t i = 0; i < frontier.size(); ++i) {
            std::map<std::string, std::vector<std::string> >::const_iterator it =
                adj.find(frontier[i]);
            if (it == adj.end()) continue;
            for (size_t j = 0; j < it->second.size(); ++j) {
                const std::string& n = it->second[j];
                if (n == b) return true;
                if (seen.insert(n).second) next.push_back(n);
            }
        }
        frontier.swap(next);
    }
    return false;
}

std::string ProbeForceFieldSystem::get_inconsistency() const {
    if (components_.empty()) return "system requires components";
    if (sites_.empty()) return "system requires sites";

    std::set<std::string> ids;
    for (size_t i = 0; i < sites_.size(); ++i) {
        if (!ids.insert(sites_[i].id).second)
            return "duplicate site ids";
    }
    for (size_t i = 0; i < sites_.size(); ++i) {
        if (components_.find(sites_[i].component) == components_.end())
            return "site " + sites_[i].id + " unknown component " + sites_[i].component;
    }
    for (size_t i = 0; i < bonds_.size(); ++i) {
        if (!ids.count(bonds_[i].site_a) || !ids.count(bonds_[i].site_b))
            return "bond references unknown site";
    }
    for (size_t i = 0; i < angles_.size(); ++i) {
        if (!ids.count(angles_[i].site_a) || !ids.count(angles_[i].site_b) ||
            !ids.count(angles_[i].site_c))
            return "angle references unknown site";
    }
    const std::vector<FFTorsion>* blocks[2] = {&dihedrals_, &impropers_};
    const char* names[2] = {"dihedrals", "impropers"};
    for (int b = 0; b < 2; ++b) {
        for (size_t i = 0; i < blocks[b]->size(); ++i) {
            const FFTorsion& t = (*blocks[b])[i];
            if (!ids.count(t.site_a) || !ids.count(t.site_b) ||
                !ids.count(t.site_c) || !ids.count(t.site_d))
                return std::string(names[b]) + " references unknown site";
        }
    }
    return "";
}

IMPBFF_END_NAMESPACE
