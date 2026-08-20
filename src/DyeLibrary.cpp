/**
 * \file DyeLibrary.cpp
 * \brief The dye as a species: its spectra, and the Förster radius they give.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/DyeLibrary.h>

#include <IMP/exception.h>

// Vendored with IMP and exported by libimp_atom; ForceFieldCIF.cpp reads its
// eighteen `_ff_*` categories through the same parser, in the same text mode.
#include "ihm_format.h"

#include <algorithm>
#include <cmath>
#include <fcntl.h>
#include <limits>
#include <sys/stat.h>
#include <unistd.h>

IMPBFF_BEGIN_NAMESPACE

const double DEFAULT_REFRACTIVE_INDEX = 1.4;
const char* const DYE_LIBRARY_CIF = "dye_library.cif";

// IMP compiles the module as one translation unit (`bff_all.cpp`), so an
// anonymous namespace is shared with every other .cpp -- and `ForceFieldCIF.cpp`
// reads the same parser with helpers of the same obvious names. The named inner
// namespace is what keeps `has`, `txt` and `Ctx` from colliding.
namespace dye_cif {

//! Numerical factor of the Förster expression with R0 in nm, the overlap
//! integral in M^-1 cm^-1 nm^4 and wavelengths in nm.
const double R0_FACTOR = 0.02108;

double nan_value() { return std::numeric_limits<double>::quiet_NaN(); }

//! A keyword carries a value only if it is in the file and neither '.' nor '?'.
bool has(ihm_keyword* k) { return k && k->in_file && !k->omitted && !k->unknown; }

std::string txt(ihm_keyword* k) {
    return has(k) && k->data.str ? std::string(k->data.str) : std::string();
}
double dbl(ihm_keyword* k) { return has(k) ? k->data.fval : nan_value(); }

//! One point of one dye's curves, before the rows are grouped and sorted.
struct SpectrumPoint {
    double wavelength, excitation, emission;
    bool operator<(const SpectrumPoint& o) const { return wavelength < o.wavelength; }
};

struct DyeRow {
    std::string name;
    double extinction_coefficient, quantum_yield;
};

struct Ctx {
    std::vector<DyeRow> dyes;
    std::map<std::string, std::vector<SpectrumPoint> > curves;

    ihm_keyword *d_name, *d_eps, *d_qy;
    ihm_keyword *s_name, *s_wavelength, *s_excitation, *s_emission;
};

void on_dye(ihm_reader*, int, void* d, ihm_error**) {
    Ctx* c = (Ctx*) d;
    if (!has(c->d_name)) return;
    DyeRow row;
    row.name = txt(c->d_name);
    row.extinction_coefficient = dbl(c->d_eps);
    row.quantum_yield = dbl(c->d_qy);
    c->dyes.push_back(row);
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

std::map<CacheKey, std::map<std::string, Dye> >& library_cache() {
    static std::map<CacheKey, std::map<std::string, Dye> > cache;
    return cache;
}

std::map<std::string, Dye> parse_dye_library(const std::string& path) {
    const int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) IMP_THROW("cannot open " << path, IOException);

    Ctx c;
    ihm_file* fh = ihm_file_new_from_fd(fd);
    ihm_reader* reader = ihm_reader_new(fh, false);          // text, not binary

    ihm_category* cat =
        ihm_category_new(reader, "_bff_dye", on_dye, NULL, NULL, &c, NULL);
    c.d_name = ihm_keyword_str_new(cat, "chromophore_name");
    c.d_eps = ihm_keyword_float_new(cat, "extinction_coefficient");
    c.d_qy = ihm_keyword_float_new(cat, "quantum_yield");

    cat = ihm_category_new(reader, "_bff_dye_spectrum", on_spectrum, NULL, NULL,
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

    std::map<std::string, Dye> out;
    for (std::size_t i = 0; i < c.dyes.size(); ++i) {
        const DyeRow& row = c.dyes[i];
        Dye dye;
        dye.name = row.name;
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

//! Fold a dye template's `_cgdye_metadata` into a species.
void apply_template_metadata(const std::string& path, Dye& dye) {
    const int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) IMP_THROW("cannot open " << path, IOException);

    MetaCtx c;
    ihm_file* fh = ihm_file_new_from_fd(fd);
    ihm_reader* reader = ihm_reader_new(fh, false);
    ihm_category* cat = ihm_category_new(reader, "_cgdye_metadata", on_metadata,
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

}  // namespace dye_cif

Dye::Dye()
    : extinction_coefficient(dye_cif::nan_value()),
      quantum_yield(dye_cif::nan_value()),
      lifetime(dye_cif::nan_value()),
      radius(dye_cif::nan_value()),
      hydrodynamic_radius(dye_cif::nan_value()) {}

void Dye::show(std::ostream& out) const {
    out << "Dye(name=" << name;
    if (quantum_yield == quantum_yield) out << ", QY=" << quantum_yield;
    if (extinction_coefficient == extinction_coefficient)
        out << ", eps=" << extinction_coefficient;
    out << ", spectrum=" << (has_spectrum() ? "yes" : "no") << ")";
}

std::map<std::string, Dye> read_dye_library(std::string path) {
    if (path.empty()) {
        path = IMP::bff::get_data_path("rotamer_library") + "/R0/" +
               std::string(DYE_LIBRARY_CIF);
    }

    // A missing file is not cached: stat fails, and the parse below raises the
    // error the caller is entitled to see on every call rather than once.
    struct stat info;
    if (stat(path.c_str(), &info) == 0) {
        dye_cif::CacheKey key;
        key.path = path;
#ifdef __APPLE__
        key.mtime_ns = (long long) info.st_mtimespec.tv_sec * 1000000000LL +
                       info.st_mtimespec.tv_nsec;
#else
        key.mtime_ns = (long long) info.st_mtim.tv_sec * 1000000000LL +
                       info.st_mtim.tv_nsec;
#endif
        key.size = (long long) info.st_size;

        std::map<dye_cif::CacheKey, std::map<std::string, Dye> >& cache = dye_cif::library_cache();
        std::map<dye_cif::CacheKey, std::map<std::string, Dye> >::const_iterator hit =
            cache.find(key);
        if (hit != cache.end()) return hit->second;

        std::map<std::string, Dye> parsed = dye_cif::parse_dye_library(path);
        cache[key] = parsed;
        return parsed;
    }
    return dye_cif::parse_dye_library(path);
}

unsigned int dye_library_cache_size() {
    return static_cast<unsigned int>(dye_cif::library_cache().size());
}

double spectral_overlap(const Dye& donor, const Dye& acceptor) {
    const std::vector<double>& wavelengths = donor.spectrum.wavelength;
    const std::vector<double>& emission = donor.spectrum.emission;
    const std::vector<double>& excitation = acceptor.spectrum.excitation;
    if (excitation.size() != wavelengths.size()) {
        IMP_THROW("donor and acceptor spectra must share one wavelength grid",
                  ValueException);
    }
    const double emission_integral = dye_cif::trapezoid(emission, wavelengths);
    if (emission_integral == 0.0) return 0.0;

    std::vector<double> weighted(wavelengths.size());
    for (std::size_t i = 0; i < wavelengths.size(); ++i) {
        const double l2 = wavelengths[i] * wavelengths[i];
        weighted[i] = emission[i] * acceptor.extinction_coefficient *
                      excitation[i] * l2 * l2;
    }
    return dye_cif::trapezoid(weighted, wavelengths) / emission_integral;
}

double forster_radius(const Dye& donor, const Dye& acceptor, double k2,
                      double refractive_index) {
    if (!donor.has_spectrum() || !acceptor.has_spectrum()) {
        IMP_THROW("both dyes need a spectrum to derive R0 ("
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
    return dye_cif::R0_FACTOR *
           std::pow(k2 * donor.quantum_yield / n4 * overlap, 1.0 / 6.0);
}

std::vector<std::string> normalize_dye_name(std::string dye_name) {
    // Trim, then split on the last run of whitespace: the trailing token is the
    // number ("488", "Thio12") and everything before it is the type.
    const std::string space = " \t\n\r\f\v";
    const std::size_t first = dye_name.find_first_not_of(space);
    const std::size_t last = dye_name.find_last_not_of(space);
    const std::string name =
        first == std::string::npos ? std::string()
                                   : dye_name.substr(first, last - first + 1);

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
    // The number must be alphanumeric, as the pattern this replaces required;
    // anything else means the name does not split.
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

Dye find_dye(std::string name, std::string library_cif, std::string template_cif) {
    const std::string compact = normalize_dye_name(name)[2];
    const std::map<std::string, Dye> library = read_dye_library(library_cif);
    std::map<std::string, Dye>::const_iterator hit = library.find(compact);
    if (hit == library.end()) {
        std::string known;
        std::map<std::string, Dye>::const_iterator it = library.begin();
        for (int i = 0; i < 3 && it != library.end(); ++i, ++it) {
            known += (i ? ", " : "") + it->first;
        }
        IMP_THROW("'" << name << "' is not in the dye library; known names "
                         "look like " << known, ValueException);
    }
    Dye dye = hit->second;
    if (!template_cif.empty()) dye_cif::apply_template_metadata(template_cif, dye);
    return dye;
}

std::vector<std::string> available_dyes(std::string library_cif) {
    const std::map<std::string, Dye> library = read_dye_library(library_cif);
    std::vector<std::string> out;
    for (std::map<std::string, Dye>::const_iterator it = library.begin();
         it != library.end(); ++it) {
        out.push_back(it->first);   // std::map is already sorted by key
    }
    return out;
}

double forster_radius_from_spectra(std::string donor, std::string acceptor,
                                   double k2, std::string library_cif,
                                   double refractive_index) {
    const std::map<std::string, Dye> library = read_dye_library(library_cif);
    const std::string d_name = normalize_dye_name(donor)[2];
    const std::string a_name = normalize_dye_name(acceptor)[2];
    std::map<std::string, Dye>::const_iterator d = library.find(d_name);
    std::map<std::string, Dye>::const_iterator a = library.find(a_name);
    if (d == library.end() || a == library.end() ||
        !d->second.has_spectrum() || !a->second.has_spectrum()) {
        IMP_THROW("No spectra for '" << donor << "' / '" << acceptor << "'",
                  ValueException);
    }
    return forster_radius(d->second, a->second, k2, refractive_index);
}

IMPBFF_END_NAMESPACE
