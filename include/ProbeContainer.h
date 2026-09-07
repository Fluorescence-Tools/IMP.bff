/**
 *  \file IMP/bff/ProbeContainer.h
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
#ifndef IMPBFF_PROBECONTAINER_H
#define IMPBFF_PROBECONTAINER_H

#include <IMP/bff/bff_config.h>

#include <IMP/bff/ProbeLibrary.h>

#include <map>
#include <string>
#include <vector>

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

#endif /* IMPBFF_PROBECONTAINER_H */
