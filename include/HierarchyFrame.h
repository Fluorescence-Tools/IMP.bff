/**
 *  \file IMP/bff/HierarchyFrame.h
 *  \brief Reading one frame out of an IMP hierarchy, in one call.
 *
 * Pulling coordinates and atom metadata out of a hierarchy from Python costs
 * about a dozen SWIG round trips **per atom** — a decorator construction, three
 * coordinate getters, a type lookup, two parent walks. On a few thousand atoms
 * that is tens of thousands of crossings to move a few kilobytes, and it
 * dominated every rotamer-library load.
 *
 * These do the walk on the C++ side and hand back arrays. Two calls rather than
 * one with out-parameters, because a *returned* `std::vector` becomes a Python
 * tuple that costs ~35-40 ns per element to build and walk back, while an out-parameter stays a wrapper
 * object numpy walks one element at a time.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_HIERARCHYFRAME_H
#define IMPBFF_HIERARCHYFRAME_H

#include <IMP/bff/bff_config.h>


#include <IMP/bff/Base.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! An atom's name as its structure spells it, without a `HET:` prefix.
/*! IMP names a heteroatom `HET: C3 `; every reader in this module wants `C3`.
    It was written out three times in three files before it was written down
    once here. */

//! Names of every XYZ leaf, in the same order as #hierarchy_atom_coordinates.
//! One frame of a structure: coordinates and the labels that identify atoms.
/*! What the rotamer scorer needs of a protein and no more, as six parallel
    arrays. A frame is one model of a multi-MODEL PDB, or one frame of a
    trajectory. */
struct IMPBFFEXPORT ProteinFrame {
    //! Flat, three per atom.
    std::vector<double> coords;
    std::vector<std::string> atom_names, atom_types, resnames, chain_ids;
    //! One per atom, -1 where the atom has no residue.
    std::vector<int> residue_indices;

    ProteinFrame() {}

    int get_n_atoms() const { return static_cast<int>(coords.size() / 3); }
    void get_coords(double** out_view, int* n_out_view) const;
    void set_coords(const std::vector<double>& v) { coords = v; }

    IMP_SHOWABLE_INLINE(ProteinFrame,
                        out << "ProteinFrame(" << get_n_atoms() << " atoms)");
};
IMP_VALUES(ProteinFrame, ProteinFrames);

//! Every XYZ leaf of a hierarchy, as a frame.
/*! The bridge for the formats this module does not read itself: a caller that
    has a hierarchy -- from RMF, from a CHARMM build, from its own sampling --
    gets a frame without going through a file. */

//! Every model of a PDB, as frames.
/*!
    A multi-MODEL PDB (a trajectory written as models) gives one frame per
    model and a plain PDB gives one. Waters are not read: they are not
    obstacles a dye is screened against in this model, and reading them costs
    the pair sum for nothing.

    \param[in] path a `.pdb` or `.ent` file
    \param[in] max_frames stop after this many; negative reads all
    \throw ValueException for any other extension. RMF is not read here --
           this module does not depend on `IMP.rmf` -- so an RMF caller loads
           its hierarchies and calls #protein_frame_from_hierarchy.
*/
IMPBFFEXPORT std::vector<ProteinFrame> load_protein_frames(
        const std::string& path, int max_frames = -1);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_HIERARCHYFRAME_H
