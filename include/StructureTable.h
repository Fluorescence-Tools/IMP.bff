/**
 *  \file IMP/bff/StructureTable.h
 *  \brief A coordinate file read into flat columns -- without naming IMP.
 *
 * Reading a PDB or an mmCIF *well* means IMP's readers: they carry the
 * selectors, the alternate-location rules, the element table, the CHARMM
 * topology that decides how big an atom is, and the masses. The core has its
 * own PDB parser (#IMP::bff::read_pdb_records) and it is the right thing when
 * geometry and a van der Waals radius are all that is wanted -- it is ~20x
 * faster because it never builds a hierarchy in order to walk back out of it.
 * It is the wrong thing when the radius has to be the one a docking score
 * measures clashes against, which is the particle's radius after
 * `IMP::atom::read_pdb` (see #IMP/bff/VdwRadii.h for why the two may not be
 * mixed), and it does not read mmCIF at all.
 *
 * So this is the *other* road, as one flat table: the columns an application
 * puts in a record array, filled by IMP's readers, with no IMP type anywhere
 * in this header. That is what lets it be wrapped without IMP's own SWIG
 * interfaces and shipped in a package that carries IMP as a private library
 * (PRD-139) -- the same rule #IMP::bff::DyeSimulation is built on.
 *
 * The walk out of the hierarchy is done **here, in C++**. That is not an
 * implementation detail: doing it in Python costs a measured ~100 us per atom
 * in SWIG traffic (twelve decorator constructions and a `Vector3D` handed to
 * numpy, per atom), which is 1.3 s of the 1.4 s it takes to read a 9315-atom
 * structure. The parsing was never the expensive part.
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_STRUCTURETABLE_H
#define IMPBFF_STRUCTURETABLE_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/Base.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! A coordinate file as parallel columns, one entry per atom.
/*!
    Every column has the same length, #get_n_atoms, except \c xyz which has
    three entries per atom. An application builds its record array by
    assigning the columns; nothing here is IMP-shaped, so nothing has to be
    walked one atom at a time on the Python side.
*/
struct IMPBFFEXPORT StructureTable {
    //! Flat, three per atom, Angstrom.
    std::vector<double> xyz;
    //! IMP's radius for the atom, Angstrom -- see #radius_no_interaction.
    std::vector<double> radius;
    //! Atomic mass, Dalton.
    std::vector<double> mass;
    //! The PDB temperature factor.
    std::vector<double> bfactor;
    //! The atom's index within its structure, as IMP reports it.
    std::vector<int> atom_id;
    //! The residue sequence number.
    std::vector<int> res_id;
    std::vector<std::string> chain, res_name, atom_name, element;

    StructureTable() {}

    int get_n_atoms() const { return static_cast<int>(radius.size()); }
    void get_xyz(double** out_view, int* n_out_view) const;
    void get_radius(double** out_view, int* n_out_view) const;
    void get_mass(double** out_view, int* n_out_view) const;
    void get_bfactor(double** out_view, int* n_out_view) const;
    void get_atom_id(int** out_view_i, int* n_out_view_i) const;
    void get_res_id(int** out_view_i, int* n_out_view_i) const;

    IMP_SHOWABLE_INLINE(StructureTable,
                        out << "StructureTable(" << get_n_atoms()
                            << " atoms)");
};
IMP_VALUES(StructureTable, StructureTables);

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

//! Is this a standard amino-acid or nucleotide residue name?
/*! The test #read_structure_table applies when \c only_standard_residues is
    set, exposed because a caller that filters its own table wants the same
    answer rather than a second list that drifts from this one. */
IMPBFFEXPORT bool is_standard_residue(const std::string& res_name);

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_STRUCTURETABLE_H
