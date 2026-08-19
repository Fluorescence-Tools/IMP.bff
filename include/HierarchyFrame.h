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
 * tuple that costs ~66 ns per element to build and walk back, while an out-parameter stays a wrapper
 * object numpy walks one element at a time.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_HIERARCHYFRAME_H
#define IMPBFF_HIERARCHYFRAME_H

#include <IMP/bff/bff_config.h>

#include <IMP/atom/Hierarchy.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Coordinates and residue index of every XYZ leaf, in hierarchy order.
/*!
    \param[in] hierarchy the frame to read
    \return four values per atom: x, y, z, and the residue index as a double
            (-1 where the leaf has no residue parent)
*/
IMPBFFEXPORT std::vector<double> hierarchy_atom_coordinates(
        IMP::atom::Hierarchy hierarchy);

//! Names of every XYZ leaf, in the same order as #hierarchy_atom_coordinates.
/*!
    \param[in] hierarchy the frame to read
    \return four strings per atom: atom name, atom type, residue name, chain id.
            Empty strings where the leaf carries no such information.

    The atom name is the last whitespace-separated field of IMP's particle name,
    matching what the Python did: IMP names an atom `"Atom CB"`, and the caller
    wants `CB`.
*/
IMPBFFEXPORT std::vector<std::string> hierarchy_atom_metadata(
        IMP::atom::Hierarchy hierarchy);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_HIERARCHYFRAME_H
