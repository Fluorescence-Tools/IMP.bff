/**
 *  \file IMP/bff/VdwRadii.h
 *  \brief Olga's van der Waals radii, keyed by atom name, and the AV's choice
 *         of radii set.
 *
 * An accessible volume is a geometric object: its size is set entirely by the
 * linker parameters and by **how big the obstacles are**. Three radii sets are
 * in circulation in this ecosystem and they are not interchangeable:
 *
 *  - **FPS**: element-keyed Bondi, in Angstrom
 *    (`../chisurf/junk/fps/Fps/data/vdW.txt`: C 1.70, N 1.55, O 1.52,
 *    P 1.80, S 1.80);
 *  - **IMP**: united-atom radii assigned by `IMP::atom::read_pdb` from the
 *    CHARMM type, carrying *implicit hydrogens* -- a carbon is 1.85-2.275 A;
 *  - **Olga**: keyed by **atom name**, not element, from `vdWRadii.json`, with
 *    a flat fallback for a name the table does not carry.
 *
 * Mixing a constant calibrated against one set with another set is a defect,
 * not a preference. Measured: with FPS's `ClashTolerance` constants against
 * IMP's united-atom radii, HIV-RT's protein-DNA interface put 102.9 of a 134.8
 * parent score into the clash term and the bootstrap spread collapsed to
 * 0.000 +/- 0.000 A (PRD-121 G2). The volume must use **one** set.
 *
 * That set is **IMP's** -- the radius each particle carries (owner,
 * 2026-09-01). Not because it reproduces Olga's published numbers best; it
 * does not, and `okf/validation/fps_screening_ab.md` states the price. The
 * reason is that an accessible volume is one half of a docking score and the
 * other half is the excluded-volume term, `clash_container`, which reads
 * #IMP::core::XYZR -- the *particles'* radii. Inflate the volume's obstacles
 * by Olga's table while the clash term measures overlap with IMP's, and the
 * two halves of one score disagree about how big an atom is, silently, in a
 * quantity neither of them reports. Consistency with IMP docking is worth
 * more than agreement with Olga's table, which is still one call away
 * (#IMP::bff::AV_RADII_OLGA) and is what an Olga-era reproduction wants.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_VDWRADII_H
#define IMPBFF_VDWRADII_H

#include <IMP/bff/bff_config.h>


#include <map>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Which van der Waals radii an accessible volume inflates its obstacles by.
enum AVRadiiSource {
    //! Whatever radius each particle carries: IMP's own, which after
    //! `IMP::atom::read_pdb` is the CHARMM-derived **united-atom** set, and
    //! which is what `IMP::core::XYZR` -- and so every clash term in a
    //! docking score -- reads. **The default**; spelled `"imp"`.
    AV_RADII_IMP = 0,
    //! Olga's name-keyed table -- #IMP::bff::olga_vdw_radius. Spelled
    //! `"olga"`; what an Olga-era reproduction asks for.
    AV_RADII_OLGA = 1
};

//! `"imp"` / `"olga"` -> #IMP::bff::AVRadiiSource; throws on anything else.
/*! `"imp"` was spelled `"model"` under schema 1.5, for the one day that
    spelling existed. It is not accepted as an alias: two names for one
    quantity is the failure this module is trying not to have, and no file
    outside this repository was ever written with the old one. */
IMPBFFEXPORT AVRadiiSource av_radii_source_from_string(const std::string &s);

//! #IMP::bff::AVRadiiSource -> `"imp"` / `"olga"`.
IMPBFFEXPORT std::string av_radii_source_to_string(AVRadiiSource s);

//! Olga's van der Waals radius for a PDB atom name, in **Angstrom**.
/*!
    The table is `Olga/src/vdWRadii.json` (128 entries), vendored as
    `data/olga_vdw_radii.csv`; the C++ table in `src/VdwRadii.cpp` is the
    definition and the data file is regenerated from it by
    #IMP::bff::olga_vdw_radii_csv().

    **Unit.** Upstream the numbers are **nanometres** -- pteros, which Olga
    reads structures through, stores coordinates in nm. The conversion is not
    inferred from the magnitudes; it is on one line of Olga's own code, where
    the coordinate and the radius are scaled together before the AV kernel sees
    either (`Olga/src/AV/Position.cpp:118-124`, `coordsVdW()`):

    \code
        xyzw.emplace_back(frame.coord.at(i)[0] * 10.0f,
                          frame.coord.at(i)[1] * 10.0f,
                          frame.coord.at(i)[2] * 10.0f,
                          pterosVDW(system, i) * 10.0f);
    \endcode

    `calculateAV()` therefore works in Angstrom, as do the linker length,
    linker width and dye radii beside it -- and as does `IMP.bff`. So the
    values here are upstream's times ten: carbon 0.17 nm -> **1.70 A**, which
    is Bondi's carbon and the same number FPS's `vdW.txt` carries. Getting the
    factor wrong would be silent and catastrophic, which is why it is pinned
    (`test/representation/test_olga_vdw_radii.py`).

    **Unknown names.** Olga does *not* fall back to an element lookup. Its
    accessor is `vdWRMap.value(name, 0.15)` (`Position.cpp:110`) -- a flat
    0.15 nm = **1.50 A** for any name the table misses, which is smaller than
    every heavy atom in it and only marginally larger than its oxygen (1.49 A).
    That is reproduced here rather than improved on: the point of adopting a
    radii set is that the result is the set's, including where the set is
    coarse. #IMP::bff::olga_vdw_unknown_atom_names() reports which atoms of a
    structure take it.
 */
IMPBFFEXPORT double olga_vdw_radius(const std::string &atom_name);

//! The radius Olga gives an atom name it does not know: **1.50 A**.
IMPBFFEXPORT double get_olga_vdw_fallback_radius();

//! Number of entries in Olga's table (128).
IMPBFFEXPORT unsigned int get_number_of_olga_vdw_radii();

//! Every atom name in Olga's table, sorted.
IMPBFFEXPORT std::vector<std::string> get_olga_vdw_atom_names();

//! `data/olga_vdw_radii.csv`, regenerated from the table above.
/*! The data file is derived, not hand-edited; a test regenerates it and fails
    on drift, exactly as `data/fps_json_schema.json` is handled. */
IMPBFFEXPORT std::string olga_vdw_radii_csv();

//! Olga's radius for a particle, by its #IMP::atom::Atom name (Angstrom).
/*! Falls back to #IMP::bff::get_olga_vdw_fallback_radius() for a particle that
    is not an atom, or whose name the table does not carry. */

//! The atom names of `ps` that Olga's table does not carry, sorted, unique.
/*! Every one of them takes the 1.50 A fallback. Worth looking at before
    trusting a volume built over an unusual residue set: on the HIV-RT fixture
    the misses are `C7` (thymine's methyl, named `C5M` in the table's era) and
    `O1P`/`O2P` (the phosphate oxygens, `OP1`/`OP2` in the table) -- 94 atoms
    of 9023, and all but the thymine methyl within 0.01 A of the radius the
    table would have given them. */

#ifndef SWIG
//! The whole table, atom name -> radius in Angstrom.
IMPBFFEXPORT const std::map<std::string, double> &get_olga_vdw_radii();

#endif

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_VDWRADII_H
