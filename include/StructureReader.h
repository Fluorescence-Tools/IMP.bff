/**
 *  \file IMP/bff/StructureReader.h
 *  \brief IMP's PDB and mmCIF readers, behind a flat table.
 *
 * Reading a structure *well* means IMP's readers: they carry the selectors,
 * the alternate-location rules, the element table, the CHARMM topology that
 * decides how big an atom is, and the masses. The core has its own PDB parser
 * (#IMP::bff::read_pdb_records) and it is the right thing when geometry and a
 * van der Waals radius are all that is wanted -- ~20x faster, because it never
 * builds a hierarchy in order to walk back out of it. It is the wrong thing
 * when the radius has to be the one a docking score measures clashes against,
 * which is the particle's radius after `IMP::atom::read_pdb`
 * (see #IMP/bff/VdwRadii.h for why the two may not be mixed), and it does not
 * read mmCIF at all.
 *
 * This is the other road. It belongs to the connection layer -- it exists
 * only where IMP is linked -- but it names no IMP type, which is what lets it
 * be wrapped without IMP's own SWIG interfaces and offered by a package that
 * carries IMP as a private library (PRD-139), the same rule
 * #IMP::bff::DyeSimulation is built on.
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_STRUCTUREREADER_H
#define IMPBFF_STRUCTUREREADER_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/Base.h>
#include <IMP/bff/StructureTable.h>

#include <string>

IMPBFF_BEGIN_NAMESPACE

//! Read a PDB or mmCIF into a #IMP::bff::StructureTable through IMP's readers.
/*!
    The formats are taken from the extension: `.pdb` and `.ent` go to
    `IMP::atom::read_pdb`, `.cif` and `.mmcif` to `IMP::atom::read_mmcif`. A
    `.gz` suffix is **not** handled -- IMP's `TextInput` does not decompress,
    and silently reading a gzip stream as text yields an empty structure
    rather than an error, so it is refused by name instead.

    Alternate locations are dropped either way, which is what both of IMP's
    selectors here do: a multi-conformer file does not yield overlapping
    copies of a residue.

    \param[in] path the coordinate file
    \param[in] keep_water keep solvent. The default drops it, which is what
               the modelling code wants; a viewer showing the deposited model
               wants it.
    \param[in] only_standard_residues drop residues that are not standard
               amino acids or nucleotides -- ligands, sugars, modified
               residues. The list is #IMP::bff::is_standard_residue.
    \param[in] radius_no_interaction report the radius at which the 6-12
               Lennard-Jones potential is **zero**, which is
               \f$2^{-1/6} R_{min}\f$, rather than \f$R_{min}\f$ itself. The
               default matches what the accessible-volume code has always
               been given.
    \throw ValueException for an extension this does not read, including a
           compressed one
    \throw IOException when the file cannot be opened
*/
IMPBFFEXPORT StructureTable read_structure_table(
        const std::string& path, bool keep_water = false,
        bool only_standard_residues = true,
        bool radius_no_interaction = true);

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_STRUCTUREREADER_H
