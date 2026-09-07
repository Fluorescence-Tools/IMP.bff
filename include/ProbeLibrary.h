#ifndef IMPBFF_PROBELIBRARY_H
#define IMPBFF_PROBELIBRARY_H

/**
 *  \file IMP/bff/ProbeLibrary.h
 *  \brief A probe as a species: what it is, its library file, and the
 *         coarse-grained system it is simulated as.
 *
 *  Three sections: the **species** (the former body of this file), the
 *  **container** (formerly `ProbeContainer.h`: a dye library as one
 *  `.mmfdb.pto`) and the **force field** (formerly `ProbeForceField.h`: a
 *  coarse-grained dye/protein system -- sites, bonds, and their terms).
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

// -------- from ProbeLibrary.h --------
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

// -------- from ProbeContainer.h --------
/**
 *  (formerly IMP/bff/ProbeContainer.h, now a section of this file)
 *  \brief A dye library as one `.mmfdb.pto`: probes, their optical properties
 *         and their spectra, in the vocabulary the database already uses.
 *
 * The bundled library is a CIF with two bff-native categories, `_bff_probe` and
 * `_bff_probe_spectrum`, invented because no dictionary in the stack declared an
 * item for a quantum yield, an extinction coefficient or a fluorescence
 * spectrum. Those items exist now — `_mmfdb_optical_property.*` and
 * `_mmfdb_spectrum.*`, dictionary 1.8 — and they describe the `probes`,
 * `optical_properties` and `spectra` tables this stack's databases already
 * have. So a dye library can be written in the same words a database row uses,
 * and this writes it.
 *
 * The point is not the file format. It is that a labelling run reads its dyes
 * from something whose every name resolves to a dictionary item, so
 * \f$R_0\f$ can be traced back to the quantum yield and the two spectra it was
 * derived from rather than taken on trust from a number in a config.
 *
 * Kinds are under `mfdb.` — the generic probe namespace, not `label.`, because
 * nothing here is specific to label-site scoring; `drot.` and `rot.bbdep.`
 * belong to the rotamer containers.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */



IMPBFF_BEGIN_NAMESPACE

//! The probe table: one row per dye, at `species` grain.
IMPBFFEXPORT extern const char* const PROBE_PTO_PROBES;
//! The spectra, one row per curve, at `spectrum` grain.
IMPBFFEXPORT extern const char* const PROBE_PTO_SPECTRA;
//! Container, artifact, operation and edge tags.
IMPBFFEXPORT extern const char* const PROBE_PTO_MODEL;

//! Write dyes as a PTO.MFDB container.
/*!
    Two tabular objects, because they count different things and the profile
    forbids joining by position:

    - `probes` — one row per dye (`species`), carrying the flrCIF identity
      (`_flr_probe_list.chromophore_name`, `probe_origin`, `probe_link_type`,
      `reactive_probe_name`, `_flr_probe_descriptor.chromophore_center_atom`)
      and its scalars as `_mmfdb_optical_property` triples of
      `property_name` / `property_value` / `unit`;
    - `spectra` — one row per (dye, `spectrum_type`) curve (`spectrum`),
      carrying the wavelength axis and the intensities as **arrays**, which is
      what `_mmfdb_spectrum.wavelengths` and `.intensity_values` are declared
      as and what the `spectra` table stores. A point-per-row serialisation of
      the same data was ten times larger, almost all of it repeated keys.

    A scalar that is NaN in the library is **absent** from the container, not
    written as zero: a quantum yield of 0 is a real and unusual claim about a
    dye, and the reference libraries carry several.

    \param[in] path the container to write
    \param[in] dyes the dyes, keyed as `_flr_probe_list.chromophore_name`
    \param[in] source what the data was made from, for the attribution tags
    \throw IOException when the container cannot be written
    \throw ValueException when a value is not a dictionary term
*/
IMPBFFEXPORT void probe_write_pto(const std::string& path,
                                const std::map<std::string, Probe>& dyes,
                                const std::string& source = "");

//! Read a dye container back.
/*! \param[in] path a container written by #probe_write_pto
    \return the dyes, keyed by chromophore name
    \throw IOException when the file is not a PTO document or holds no probes */
IMPBFFEXPORT std::map<std::string, Probe> probe_read_pto(const std::string& path);

//! Convert the bundled CIF library to a container.
/*!
    \param[in] path where to write
    \param[in] library_cif the source library, or empty for the shipped one
    \return how many dyes were written
*/
IMPBFFEXPORT int probe_library_to_pto(const std::string& path,
                                    const std::string& library_cif = "");

//! The Forster radius of a donor/acceptor pair named in a container, Angstrom.
/*!
    The convenience a labelling run wants: two names in, \f$R_0\f$ out, with the
    spectra and the quantum yield coming from a file whose vocabulary is
    checkable. Names resolve through #IMP::bff::resolve_probe_name, so
    `"Alexa488"` works as well as `"AlexaFluor488"`.

    **Angstrom**, like every length in this package --
    #IMP::bff::LlFretOptions::forster_radius (default 52),
    #IMP::bff::av_distance, #IMP::bff::fret_efficiency. This function used to
    be the one place nanometres met Angstrom, because
    #IMP::bff::forster_radius returned nm; since 2026-08-25 that conversion
    lives inside `forster_radius` itself, so this one simply forwards. The
    reason it mattered either way: a caller who forgets the factor of ten
    scores every pair against an \f$R_0\f$ ten times too small, which does not
    look like an error -- it looks like a protein where no pair is worth
    measuring.

    \param[in] path a dye container
    \param[in] donor,acceptor the two dyes
    \param[in] kappa2 the orientation factor; 2/3 is the isotropic average and
               is an assumption, not a measurement
    \param[in] refractive_index of the medium between them
    \return \f$R_0\f$ in **Angstrom**
    \throw ValueException when either dye is missing or lacks a spectrum
*/
IMPBFFEXPORT double probe_pto_forster_radius(const std::string& path,
                                           const std::string& donor,
                                           const std::string& acceptor,
                                           double kappa2 = 2.0 / 3.0,
                                           double refractive_index = 1.4);

IMPBFF_END_NAMESPACE

 /* IMPBFF_PROBECONTAINER_H */

// -------- from ProbeForceField.h --------
/**
 *  (formerly IMP/bff/ProbeForceField.h, now a section of this file)
 *  \brief A coarse-grained dye/protein system: sites, bonds, and their terms.
 *
 * The topology and parameters of one simulation system — what
 * #IMP::bff::cgprobe builds, writes to mmCIF, reads back, and hands to IMP to
 * become particles and restraints.
 *
 * Every term is a named field of a typed record: a bond carries `site_a`,
 * `site_b`, `length` and `type_id` by name, so what a slot means is stated
 * rather than conventional, and a consumer reads `system.bonds` without
 * having to guard against the key being absent.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */


#include <set>

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
    //! Key into ProbeForceFieldSystem::get_bond_types().
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
    (`site_id`, or `site_id_start`/`site_id_end`). Carrying whichever the file
    happened to use would hand a consumer `4` or `"CX4/S1"` for the same site
    with no way to tell which; one canonical form is the point of having a type
    at all, and resolving an index to its id is the reader's job.
*/
class IMPBFFEXPORT ProbeForceFieldSystem {
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
    ProbeForceFieldSystem(const std::string& name = "") : name_(name) {}

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

    //! Component names the sampler may **not** move.
    std::vector<std::string> get_fixed_components() const;
    //! Component names the sampler may move -- the dyes.
    std::vector<std::string> get_mobile_components() const;

    // ----------------------------------------------------------------------
    // The molecular graph
    //
    // Derived from the bonds this object already holds, so it belongs here
    // rather than in three Python copies over the same data. `cgprobe.sim`,
    // `cgprobe.topology` and `scoring` each had their own; they agreed on the
    // shipped system and had drifted in shape (a set of tuples, a set of
    // frozensets, a defaultdict of sets).
    // ----------------------------------------------------------------------

    //! Site ids adjacent to each site, keyed by site id.
    std::map<std::string, std::vector<std::string> > get_bonded_neighbors() const;

    //! The 1-2, 1-3 and 1-4 pairs, each as a sorted (a, b).
    /** Bonds give 1-2, angle ends 1-3, dihedral ends 1-4. Impropers contribute
        their a-c, a-d and b-d pairs when `include_impropers` is true. That
        flag is moot for a combined system today, because no builder of one
        fills `impropers`: okf/validation/impropers_are_dropped.md.

        A `std::set`, not a vector: its caller tests membership and takes
        subsets, and a deduplicated, sorted container is what the answer is. */
    std::set<std::pair<std::string, std::string> > get_exclusions(
            bool include_impropers = true) const;

    //! Simple cycles of at most `max_len` sites, each as sorted site ids.
    std::vector<std::vector<std::string> > find_rings(int max_len = 8) const;

    //! Whether `a` reaches `b` in at most `max_depth` bonds.
    bool is_within_bonds(const std::string& a, const std::string& b,
                         int max_depth) const;

    //! Why this system is not self-consistent, or an empty string.
    /** Components and sites present, site ids unique, every site's component
        declared, and every bonded term referring to sites that exist. The
        caller decides whether that is an error -- `cgprobe.sim` raises; a reader
        may prefer to report. */
    std::string get_inconsistency() const;

    IMP_SHOWABLE_INLINE(ProbeForceFieldSystem,
                        out << "ProbeForceFieldSystem(\"" << name_ << "\", "
                            << sites_.size() << " sites, " << bonds_.size()
                            << " bonds, " << components_.size() << " components)");
};
IMP_VALUES(ProbeForceFieldSystem, ProbeForceFieldSystems);

//! Build a typed system from a JSON object (the dict-literal convenience).
/*!
    The topology builder assembles a system as a plain JSON-shaped object --
    easier to author and to diff than eighteen typed vectors -- and this is
    the one conversion into #IMP::bff::ProbeForceFieldSystem. The accepted shape
    mirrors the reader: `name`, `components` (`{id: {mol2_path|mol2,
    pdb_path|pdb, role}}`), `sites` (a list of `{id, component, atom_name,
    element, site_no, site_serial, radius, mass}`), `groups`/`rb_groups`/
    `md_fixed_groups` (`{gid: [site-id...]}`, integers resolved to site ids by
    `site_no`), `fixed_groups`, `bond_types`/`angle_types` (`{tid: {k}}`),
    `torsion_types`/`improper_types` (`{tid: {periodicity, phase_rad, k}}`),
    `lj_types` (`{tid: {element, rmin_half, epsilon}}`), `bonds`
    (`[a, b, length, tid]`), `angles` (`[a, b, c, theta, tid]`), `dihedrals`/
    `impropers` (`[a, b, c, d, tid]`), `probes`, `nonbonded`
    (`{enabled, k, cutoff_A}`) and `sampling` (`{temperature_K, friction_ps,
    timestep_fs, n_steps, write_every, minimize_steps}`).
    Missing keys take the same defaults the reader does.
    \throw ValueException for a row of the wrong arity or a non-object root
*/
IMPBFFEXPORT ProbeForceFieldSystem forcefield_system_from_json(
        const std::string& json);

IMPBFF_END_NAMESPACE


#endif  // IMPBFF_PROBELIBRARY_H
