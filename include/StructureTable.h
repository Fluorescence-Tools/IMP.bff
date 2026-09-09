/**
 *  \file IMP/bff/StructureTable.h
 *  \brief A structure as flat columns: the record every reader and writer
 *         here agrees on.
 *
 * One structure, one table, whatever produced it. #IMP::bff::read_structure_table
 * (#IMP/bff/StructureReader.h) fills it from a PDB or an mmCIF through IMP's
 * readers; #IMP::bff::RmfStructureWriter writes it as a trajectory. It is a
 * plain record -- no IMP type, no RMF type -- which is why it can sit between
 * the two and be wrapped without either one's SWIG interfaces.
 *
 * The columns are the ones an application puts in a record array. Filling
 * them is done in C++ on purpose: walking a hierarchy from Python costs a
 * measured ~100 us per atom in SWIG traffic (twelve decorator constructions
 * and a `Vector3D` handed to numpy, per atom), 1.3 s of the 1.4 s it took to
 * read a 9315-atom structure.
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

    //! The coordinate column is the one that counts.
    /*! A table under construction has its columns assigned one at a time, so
        any other choice makes the count depend on assignment order. */
    int get_n_atoms() const { return static_cast<int>(xyz.size() / 3); }

    void get_xyz(double** out_view, int* n_out_view) const;
    void get_radius(double** out_view, int* n_out_view) const;
    void get_mass(double** out_view, int* n_out_view) const;
    void get_bfactor(double** out_view, int* n_out_view) const;
    void get_atom_id(int** out_view_i, int* n_out_view_i) const;
    void get_res_id(int** out_view_i, int* n_out_view_i) const;

    // A table is built as well as read: a caller with coordinates of its own
    // -- a trajectory frame, a generated chain, an application's record array
    // -- assembles one column by column and hands it to a writer.
    void set_xyz(const std::vector<double>& v) { xyz = v; }
    void set_radius(const std::vector<double>& v) { radius = v; }
    void set_mass(const std::vector<double>& v) { mass = v; }
    void set_bfactor(const std::vector<double>& v) { bfactor = v; }
    void set_atom_id(const std::vector<int>& v) { atom_id = v; }
    void set_res_id(const std::vector<int>& v) { res_id = v; }

    IMP_SHOWABLE_INLINE(StructureTable,
                        out << "StructureTable(" << get_n_atoms()
                            << " atoms)");
};
IMP_VALUES(StructureTable, StructureTables);

//! Is this a standard amino-acid or nucleotide residue name?
/*! The test #read_structure_table applies when \c only_standard_residues is
    set, exposed because a caller that filters its own table wants the same
    answer rather than a second list that drifts from this one. */
IMPBFFEXPORT bool is_standard_residue(const std::string& res_name);

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_STRUCTURETABLE_H
