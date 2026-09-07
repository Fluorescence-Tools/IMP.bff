/**
 * \file IMP/bff/Mol2IO.h
 * \brief Reading a MOL2 component into force-field sites and bonds.
 *
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_MOL2IO_H
#define IMPBFF_MOL2IO_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/ProbeLibrary.h>

#include <IMP/bff/Base.h>

#include <map>
#include <string>
#include <utility>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! One MOL2 atom, as the force-field builder needs it.
struct IMPBFFEXPORT Mol2Atom {
    int serial;
    std::string component;
    std::string atom_name;
    std::string resname;
    std::string element;
    double x, y, z;

    Mol2Atom() : serial(0), x(0), y(0), z(0) {}

    IMP_SHOWABLE_INLINE(Mol2Atom, out << "Mol2Atom(" << serial << ", "
                                      << atom_name << ")");
};

//! The atoms and bonds of a MOL2 file.
struct IMPBFFEXPORT Mol2Component {
    std::vector<Mol2Atom> atoms;
    std::vector<std::pair<int, int> > bonds;   //!< sorted serial pairs
};

//! Read `path` as one component named `component`.
/** The atom *name* comes from column 2 of `@<TRIPOS>ATOM` rather than from the
    TRIPOS type: IMP's `read_mol2` maps `C.3` to `C3` and loses `C12`/`N1`,
    which are the names a template's features and impropers refer to. The
    element is derived from that name. */
IMPBFFEXPORT Mol2Component read_mol2_component(const std::string& path,
                                               const std::string& component);

//! `{serial: atom name}` straight from the `@<TRIPOS>ATOM` records.
IMPBFFEXPORT std::map<int, std::string> read_mol2_atom_names(const std::string& path);

//! The element implied by a MOL2 atom name.
/** The **first letter** of the leading alphabetic run, uppercased -- so `C12`
    and `CL3` are both carbon. That is wrong for two-letter elements and is
    kept deliberately: every site's LJ type is keyed on it, so correcting it
    would move the force field.
    `C` when the name starts with no letter. */
IMPBFFEXPORT std::string element_from_atom_name(const std::string& atom_name);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_MOL2IO_H
