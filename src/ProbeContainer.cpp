/**
 * \file ProbeContainer.cpp
 * \brief A dye library as one `.mmfdb.pto`.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/ProbeContainer.h>

#include <IMP/bff/Pto.h>
#include <IMP/bff/bff_config.h>
#include <IMP/bff/internal/json.h>

#include <IMP/bff/Base.h>

#include <cmath>
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
    // through find_probe, or the container is harder to use than the CIF it
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
