/**
 *  \file IMP/bff/ProbeLibrary.h
 *  \brief A probe as a species: what it is, independent of how it is modelled.
 *
 * A **probe** is anything attached to a molecule to report on it: an organic
 * dye, a genetically encoded fluorescent protein, a nitroxide spin label for
 * EPR. They differ in what properties they have — a spin label has no quantum
 * yield and no spectra — but not in what a library of them is for, so one type
 * carries them all and #IMP::bff::ProbeType says which kind a row is.
 *
 * This was `ProbeLibrary`, and modelling only dyes was already untrue of the
 * package: `ProbeAttachment.h` spoke of probes while holding a `Dye`, the
 * shipped rotamer stores are partitioned dyes / spinlabels / sidechains, and
 * the flrCIF vocabulary these fields follow calls all of them probes. The old
 * header was deleted rather than kept as a shim -- one name for the thing.
 *
 * A probe's spectra, quantum yield, extinction coefficient, transition dipole and
 * formal charges are what it *is*; an accessible volume, a rotamer library, a
 * coarse-grained conformer set and an MD trajectory are ways of *representing*
 * it. Only the former is here.
 *
 * The practical consequence is that **R0 is derived, not supplied**:
 * `forster_radius` takes two dyes, kappa^2 and the medium's refractive index,
 * so no call site has to carry `forster_radius=52.0`.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_PROBELIBRARY_H
#define IMPBFF_PROBELIBRARY_H

#include <IMP/bff/bff_config.h>

#include <IMP/bff/Base.h>

#include <map>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Excitation and emission on a shared wavelength grid.
/*!
    One grid for both curves because the overlap integral needs them
    point-for-point; the reader enforces it.
*/
struct IMPBFFEXPORT Spectrum {
    //! nm.
    std::vector<double> wavelength;
    //! Normalised to a maximum of 1.
    std::vector<double> excitation;
    //! Normalised to a maximum of 1.
    std::vector<double> emission;

    bool empty() const { return wavelength.empty(); }
    unsigned int size() const { return static_cast<unsigned int>(wavelength.size()); }

    IMP_SHOWABLE_INLINE(Spectrum, out << "Spectrum(n=" << wavelength.size() << ")");
};

//! A dye species, independent of how it is modelled.
/*!
    Optional numbers are **NaN** when the library does not carry them, which is
    the C++ spelling of the `None` this carried as a dataclass. `has_spectrum()`
    answers the one that is not a number.

    Field names follow flrCIF where the dictionary has an item. Three do not
    exist in any dictionary in the stack -- quantum yield, extinction
    coefficient and a spectrum -- checked across all ten `.dic` files in
    `../mmfdb/src/mmfdb/data`; the only near matches are
    `_em_detector.detective_quantum_efficiency` and the NMR spectral
    categories. Those three are bff-native, deliberately rather than by
    omission, and should be proposed for `mmfdb_flr_ext.dic`.

    Nothing here says how the dye's positions are enumerated. That is the
    representation's business.
*/
//! What kind of probe a library row describes.
/*!
    Not the same question as `_flr_probe_list.probe_origin`, which flrCIF uses
    for **how it got there** — a fluorescent protein is `intrinsic` because it
    is expressed as part of the molecule, a dye and a spin label are
    `extrinsic` because they are attached. This says **what it is**, which is
    what decides whether a given property can exist at all.
*/
enum ProbeType {
    //! An organic fluorophore, attached through a reactive group.
    PROBE_DYE = 0,
    //! Genetically encoded, and therefore part of the molecule.
    PROBE_FLUORESCENT_PROTEIN = 1,
    //! A nitroxide or similar for EPR. **Not fluorescent**: no quantum yield,
    //! no extinction coefficient, no spectra, and no Förster radius.
    PROBE_SPIN_LABEL = 2,
    //! The library did not say.
    PROBE_UNSPECIFIED = 3
};

//! Is this the kind of probe a Förster radius can be computed for?
/*! True for a dye or a fluorescent protein; false for a spin label, which
    absorbs and emits nothing. `PROBE_UNSPECIFIED` counts as fluorescent: a
    library that does not state a type is a dye library, and refusing it would
    break callers over a column they never had. A spin label has to say so. */
IMPBFFEXPORT bool probe_is_fluorescent(ProbeType type);

//! The controlled name of a probe type: `dye`, `fluorescent_protein`,
//! `spin_label`, or `unspecified`.
/*! These are the strings the `_bff_probe.probe_type` column holds and the
    strings a container round-trips, so they are the one spelling. */
IMPBFFEXPORT std::string probe_type_to_string(ProbeType type);

//! The inverse. Anything unrecognised -- including an empty string, which is
//! what a library written before the column existed gives -- is
//! #IMP::bff::PROBE_UNSPECIFIED rather than an error: a probe whose kind was
//! not stated is still a usable species.
IMPBFFEXPORT ProbeType probe_type_from_string(std::string name);

struct IMPBFFEXPORT Probe {
    //! What kind of probe this is; see #IMP::bff::ProbeType.
    ProbeType probe_type;
    //! The compact name, e.g. `AlexaFluor488`. flrCIF
    //! `_flr_probe_list.chromophore_name`.
    std::string name;
    //! Who sells it -- `ATTO`, `Lumiprobe`, `AlexaFluor`. Empty when the
    //! library does not say. Not an identifier: #name is.
    std::string vendor;
    //! Excitation/emission curves; empty when unknown.
    Spectrum spectrum;
    //! Molar extinction at the absorption maximum, M^-1 cm^-1. NaN if unknown.
    double extinction_coefficient;
    //! Fluorescence quantum yield of the free dye. NaN if unknown.
    double quantum_yield;
    //! Unquenched fluorescence lifetime, ns. NaN if unknown. flrCIF
    //! `_flr_reference_measurement_lifetime.lifetime`.
    double lifetime;
    //! Steric radius, Angstrom. NaN if unknown.
    double radius;
    //! Angstrom, for the rotational correlation time and the translational
    //! diffusion coefficient. NaN if unknown.
    double hydrodynamic_radius;
    //! The two atom names whose separation defines the transition dipole.
    //! Empty for an isotropic model; needed for kappa^2.
    std::vector<std::string> dipole_atoms;
    //! flrCIF `_flr_probe_descriptor.chromophore_center_atom`.
    std::string chromophore_center_atom;
    //! The labelling reagent, e.g. the maleimide form. flrCIF distinguishes the
    //! reactive probe from its chromophore; `name` is the chromophore.
    std::string reactive_probe_name;
    //! flrCIF `_flr_probe_list.probe_origin`.
    std::string probe_origin;
    //! flrCIF `_flr_probe_list.probe_link_type`.
    std::string probe_link_type;
    std::vector<std::string> positive_atoms;
    std::vector<std::string> negative_atoms;

    Probe();

    bool has_spectrum() const { return !spectrum.empty(); }

    IMP_SHOWABLE(Probe);
};

IMP_VALUES(Probe, Probes);
IMP_VALUES(Spectrum, Spectrums);

//! Refractive index of the medium when the caller does not say.
/*! 1.4 is the conventional value for a dye on a protein surface -- between
    water (1.33) and protein interior (~1.6). It was written as the literal
    `1.4**4` inside the calculation before PRD-113, which is why nothing could
    ask what R0 would be in a different solvent. */
IMPBFFEXPORT extern const double DEFAULT_REFRACTIVE_INDEX;

//! Name of the bundled library inside the rotamer-library data directory.
IMPBFFEXPORT extern const char* const PROBE_LIBRARY_CIF;

//! Every dye in a library CIF, keyed by chromophore name.
/*!
    Two bff-native categories: `_bff_probe` (`chromophore_name`, `probe_type`, `vendor`,
    `chromophore_number`, `extinction_coefficient`, `quantum_yield`) and
    `_bff_probe_spectrum` (`chromophore_name`, `wavelength`, `excitation`,
    `emission`), read through the `ihm` C reader IMP vendors. A dye with a table
    entry but no curves comes back without a spectrum -- it is still a usable
    species, it just cannot derive R0.

    **Cached on (resolved path, mtime, size).** The bundled library is a shipped
    read-only file and a spectrum lookup does not change it, but the parse is
    not cheap: profiling one FRETpredict comparison found it read 21 times for
    the same file, 2.2 s of a 13 s run. mtime and size are in the key so editing
    the file during a session is picked up; only an edit that changes neither
    would be missed, which is not a thing that happens to a data file.

    \param[in] path the CIF file; empty for the bundled one
    \return `{chromophore_name: Probe}`
*/
IMPBFFEXPORT std::map<std::string, Probe> read_probe_library(std::string path = "");

//! How many parsed libraries the cache is holding. For tests.
IMPBFFEXPORT unsigned int probe_library_cache_size();

//! The overlap integral J, in M^-1 cm^-1 nm^4.
/*!
    Donor emission against acceptor extinction, weighted by lambda^4 and
    normalised by the donor's emission integral. Trapezoidal, on the shared
    wavelength grid.
*/
IMPBFFEXPORT double spectral_overlap(const Probe& donor, const Probe& acceptor);

//! R0 in **Angstrom** for a pair of dyes, in a medium.
/*!
    The species-level entry point: everything it needs is a property of the two
    dyes and the solvent, so nothing has to be threaded in from a call site.

    **Angstrom, like every other length in this package.** The spectra it
    integrates are nanometres -- that is what a spectrometer, a datasheet and
    the overlap integral are in, and the wavelength grids stay in it -- but an
    R0 exists to be compared with a distance, and distances here come from
    coordinates. The conversion happens once, inside this function, rather than
    at each of the six call sites that used to multiply by ten. It read
    nanometres until 2026-08-25, and a caller shipped the ten-times-too-small
    value most of the way to a program before catching it.

    \param[in] donor needs a spectrum and a quantum yield
    \param[in] acceptor needs a spectrum and an extinction coefficient
    \param[in] k2 orientation factor; the isotropic 2/3 by default
    \param[in] refractive_index of the medium between the dyes
*/
IMPBFFEXPORT double forster_radius(
        const Probe& donor, const Probe& acceptor, double k2 = 2.0 / 3.0,
        double refractive_index = 1.4);

//! Assemble a dye from the bundled data.
/*!
    \param[in] name e.g. `"AlexaFluor 488"` or `"AlexaFluor488"`
    \param[in] library_cif optional path to a dye-library CIF; empty for the
        bundled one
    \param[in] template_cif optional dye template, for the molecular half --
        transition-dipole atoms, centre atom, formal charges. Reads the
        `_cgprobe_metadata` category.
    \throw ValueException if the dye is not in the table
*/
IMPBFFEXPORT Probe find_probe(std::string name, std::string library_cif = "",
                          std::string template_cif = "");

//! Chromophore names of every dye in a library, sorted.
IMPBFFEXPORT std::vector<std::string> available_probes(std::string library_cif = "");

//! R0 in **Angstrom** for a dye pair *named by string*.
/*!
    The name-based route, kept because the rotamer code addresses dyes by
    string. `forster_radius` is the same calculation over two `Probe` values and
    is what new code should use.

    Names resolve through #IMP::bff::resolve_probe_name, exactly as
    #IMP::bff::find_probe does, so `"Alexa488"` reaches `AlexaFluor488` here too.
    It used to match on whitespace alone, which meant a spelling worked through
    one door and not the other -- worse than neither working, because the
    failure looked like a missing dye rather than a missing alias.

    \param[in] donor,acceptor dye names, e.g. `"AlexaFluor 488"`
    \param[in] k2 orientation factor
    \param[in] library_cif optional path to a dye-library CIF
    \param[in] refractive_index of the medium between the dyes
*/
IMPBFFEXPORT double forster_radius_from_spectra(
        std::string donor, std::string acceptor, double k2,
        std::string library_cif = "", double refractive_index = 1.4);

//! Split a FRETpredict-style dye name into `(type, number, compact)`.
/*! `"AlexaFluor 488"` becomes `("AlexaFluor", "488", "AlexaFluor488")`. A name
    that does not split is returned three times with its spaces removed. */
//! The library name a dye is known by, from whatever a caller called it.
/*!
    The library keys are the vendor's full spellings -- `AlexaFluor488`,
    `LumiprobeCy3`, `ATTO647N` -- and nobody types those. An experimenter, a
    paper and the Labelizer backend all say `Alexa488`, `Cy3`, `Atto647N`, and
    an exact-match lookup rejects every one of them while holding the right
    row.

    Resolution is by three rules, applied in order and required to land on
    exactly one entry:

    1. ignore case and any non-alphanumeric character;
    2. expand a known vendor abbreviation (`Alexa` -> `AlexaFluor`,
       `Cy` -> `LumiprobeCy`);
    3. accept only a unique hit.

    Ambiguity is an error, never a guess: `Cy5` resolves because
    `LumiprobeCy5` matches exactly and `LumiprobeCy55` does not, but a query
    matching two entries raises rather than picking one.

    \param[in] name what the caller called it
    \param[in] library_cif the library, or empty for the shipped one
    \return the library key
    \throw ValueException when nothing matches, or more than one does
*/
//! A dye name folded to the form two spellings are compared in.
/*! Lowercase, alphanumeric only. Exposed so a caller holding its own dye table
    -- a container rather than the shipped library -- can resolve names by the
    same rules instead of reimplementing them slightly differently. */
IMPBFFEXPORT std::string probe_fold_name(const std::string& name);

//! Does \p folded name the dye whose library key is \p key?
/*! Applies the vendor-abbreviation expansions (`Alexa` -> `AlexaFluor`,
    `Cy` -> `LumiprobeCy`) as well as the fold. */
IMPBFFEXPORT bool probe_name_matches(const std::string& folded,
                                   const std::string& key);

IMPBFFEXPORT std::string resolve_probe_name(std::string name,
                                          std::string library_cif = "");

IMPBFFEXPORT std::vector<std::string> normalize_probe_name(std::string probe_name);

//! The flrCIF item for each `Probe` field, or empty where none exists.
/*!
    flrCIF's word for a dye is "probe", and for the fluorescent moiety
    "chromophore". Quantum yield, extinction coefficient and a spectrum have no
    dictionary item in the stack -- a bff-native category -- so those keys map
    to the empty string, not to a fabricated item.
*/
IMPBFFEXPORT std::map<std::string, std::string> probe_flrcif_items();

//! The flrCIF items for an R0 derivation's inputs and output.
/*! `index_of_refraction` and `kappa_squared` are not in the upstream IHM-FLR
    dictionary -- they are added by `mmfdb_flr_ext.dic`, which is what makes a
    stored R0 reproducible rather than a bare number. */
IMPBFFEXPORT std::map<std::string, std::string> forster_radius_flrcif_items();

IMPBFF_END_NAMESPACE

#endif //IMPBFF_PROBELIBRARY_H
