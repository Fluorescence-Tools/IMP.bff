/**
 *  \file IMP/bff/StructureIO.h
 *  \brief Structures in and out: PDB, MOL2, mmCIF, and the tables beside them.
 *
 * Nothing here is about fluorescence. It is serialisation, and it is a separate
 * header from `FPSIO.h` because a labelling file and a coordinate file are
 * different formats that change for different reasons — they shared one module
 * until PRD-113 stage 7, which is how "the fps reader" came to be the thing that
 * writes PDBs.
 *
 * The RMF door is **not** here. Writing RMF needs `IMP.rmf`, which is not one of
 * this module's `required_modules`, and widening the dependency graph for two
 * functions is the wrong trade; they are `%pythoncode` in
 * `pyext/IMP_bff.structureio.i` with a lazy `import RMF`, the same shape
 * `AVNetworkRestraintWrapper` has.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_STRUCTUREIO_H
#define IMPBFF_STRUCTUREIO_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/AVBuilder.h>

#include <IMP/Model.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/showable_macros.h>
#include <IMP/value_macros.h>

#include <map>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

// --------------------------------------------------------------------------
// MOL2: what a dye topology has to be in before a force field can be built
// --------------------------------------------------------------------------

//! A bond, by the two atom **serials** it joins, lower first.
struct IMPBFFEXPORT AtomBond {
    int a, b;
    AtomBond(int a = 0, int b = 0) : a(a < b ? a : b), b(a < b ? b : a) {}
    bool operator<(const AtomBond& o) const {
        return a != o.a ? a < o.a : b < o.b;
    }
    bool operator==(const AtomBond& o) const { return a == o.a && b == o.b; }
    IMP_SHOWABLE_INLINE(AtomBond, out << "AtomBond(" << a << ", " << b << ")");
};
IMP_VALUES(AtomBond, AtomBonds);

// The ATOM/HETATM records of a PDB are #IMP::bff::read_pdb_records, in
// `AVBuilder.h`. It is cached on (path, mtime, size), which the reader this
// replaced was not, and it now carries the serial, the residue name and the
// element that reader carried -- there were two records and two parsers for one
// file format until 2026-08-20.

//! The bonds a PDB's `CONECT` records declare. Empty when it has none.
IMPBFFEXPORT std::vector<AtomBond> parse_conect_bonds(const std::string& path);

//! Bonds from geometry, by a covalent-radius cutoff.
/*! A guess, and that is why it is not the default: `CONECT` records are what a
    PDB says, and this is what a distance suggests.
    \param[in] scale multiplies the sum of the two covalent radii */
IMPBFFEXPORT std::vector<AtomBond> infer_bonds(
        const std::vector<PDBAtomRecord>& atoms, double scale = 1.22);

//! Write a Tripos MOL2 file.
/*! The library could **parse** MOL2 but not write one, and MOL2 is the format a
    dye topology has to be in before a force field can be built for it. The
    writer lived in a CLI driver, which is why it was invisible: a capability
    reachable only by running a script is a capability nobody finds. */
IMPBFFEXPORT void write_mol2(const std::string& path,
                             const std::vector<PDBAtomRecord>& atoms,
                             const std::vector<AtomBond>& bonds,
                             const std::string& mol_name = "MOL");

// --------------------------------------------------------------------------
// PMI stat files: live docking progress
// --------------------------------------------------------------------------

//! Completed frames recorded in a PMI `stat.*.out`.
/*! One header line plus one line per frame. Zero for a missing or unreadable
    file — a progress bar asking how far a job has got wants a number. */
IMPBFFEXPORT int count_frames(const std::string& stat_path);

//! `(frames, scores)` from a PMI stat file or a `frame,score` CSV.
/*!
    The `ReplicaExchange` macro appends one `repr(dict)` per frame to
    `stat.0.out`, after a first header line whose values name the columns. Only
    two columns are wanted — the frame index and the total score — so the parsing
    is deliberately tolerant: an unreadable line is skipped, including the header
    itself, which carries an `environ(...)` call no literal parser accepts.

    A partial file from a running job is fine, and so is the plain
    `frame,score` CSV a minimisation writes.

    \param[out] out_view,n_out_view two values per point: frame, then score
*/
IMPBFFEXPORT void read_score_series(const std::string& stat_path,
                                    double** out_view, int* n_out_view);

// --------------------------------------------------------------------------
// coordinates
// --------------------------------------------------------------------------

//! Write `(N, 3)` coordinates as a single-model PDB.
/*!
    \param[in] coords flat, three per atom
    \param[in] transform empty, three values (a translation), nine (a rotation)
               or sixteen (a homogeneous 4x4, row-major)
    \throw ValueException for any other transform length
*/
IMPBFFEXPORT void write_pdb(const std::vector<double>& coords,
                            const std::string& path,
                            const std::string& chain = "A",
                            const std::string& res_name = "ALA",
                            const std::vector<double>& transform =
                                    std::vector<double>(),
                            int model_index = 0);

//! Apply a transform to `(N, 3)` coordinates.
/*! \param[in] transform three values (a translation), nine (a rotation) or
           sixteen (a homogeneous 4x4, row-major); empty is the identity
    \throw ValueException for any other length */
IMPBFFEXPORT void apply_transform(const std::vector<double>& coords,
                                  const std::vector<double>& transform,
                                  double** out_view, int* n_out_view);

//! Read a PDB into \p m and return its hierarchy.
/*! `NonWaterPDBSelector`, which is what every caller of this wanted. */
IMPBFFEXPORT IMP::atom::Hierarchy read_pdb_hierarchy(const std::string& path,
                                                     IMP::Model* m);

//! `(N, 3)` coordinates of a hierarchy's XYZ leaves, in hierarchy order.
IMPBFFEXPORT void structure_coordinates(IMP::atom::Hierarchy hierarchy,
                                        double** out_view, int* n_out_view);

//! Load a PDB and return its `(N, 3)` coordinates.
IMPBFFEXPORT void load_structure(const std::string& path, double** out_view,
                                 int* n_out_view);

//! RMSD between two coordinate sets, optionally after superposition.
/*!
    \param[in] coords_a,coords_b flat, three per atom, same length
    \param[in] selection_mask one per atom; only the selected atoms enter the
               RMSD **and** the superposition. Empty selects everything.
    \param[in] superpose Kabsch-align \p coords_a onto \p coords_b first. The
               rotation is fitted on the selection and applied to *all* of
               `coords_a`, which is what makes a selection a *reference* rather
               than a crop.
    \throw ValueException on a length mismatch
*/
IMPBFFEXPORT double compute_rmsd(const std::vector<double>& coords_a,
                                 const std::vector<double>& coords_b,
                                 const std::vector<int>& selection_mask =
                                         std::vector<int>(),
                                 bool superpose = false);

//! Convert a dye PDB to mmCIF with `_atom_site` records.
/*! \param[in] dye_id the data block name; empty takes the PDB's stem */
IMPBFFEXPORT void convert_pdb_to_cif(const std::string& pdb_path,
                                     const std::string& cif_path,
                                     const std::string& dye_id = "");

// --------------------------------------------------------------------------
// the tables beside a structure
// --------------------------------------------------------------------------

//! One row of a crosslink table.
struct IMPBFFEXPORT CrossLink {
    std::string protein_1, protein_2;
    int residue_1, residue_2;

    CrossLink() : residue_1(0), residue_2(0) {}

    IMP_SHOWABLE_INLINE(CrossLink, out << "CrossLink(" << protein_1 << ":"
                                       << residue_1 << " - " << protein_2 << ":"
                                       << residue_2 << ")");
};
IMP_VALUES(CrossLink, CrossLinks);

//! Read a crosslink table: a header line, then `protein,res,protein,res`.
/*! A row whose residue fields are not integers is skipped, not an error. */
IMPBFFEXPORT std::vector<CrossLink> read_xlink_table(const std::string& path);

//! An atom named the way a flexible-residue file names one.
struct IMPBFFEXPORT AtomReference {
    std::string chain_identifier, atom_name;
    int residue_seq_number;

    AtomReference() : residue_seq_number(0) {}
    AtomReference(const std::string& chain, int residue,
                  const std::string& atom = "")
        : chain_identifier(chain), atom_name(atom),
          residue_seq_number(residue) {}

    IMP_SHOWABLE_INLINE(AtomReference,
                        out << "AtomReference(" << chain_identifier << ":"
                            << residue_seq_number << ":" << atom_name << ")");
};
IMP_VALUES(AtomReference, AtomReferences);

//! Resolve the flexible residues a flexfit file names.
/*! \param[in] residues one entry per flexible residue; `atom_name` is ignored
    \return the first selected particle of each, in the order given */
IMPBFFEXPORT IMP::ParticleIndexes select_flexible_residues(
        IMP::atom::Hierarchy hierarchy,
        const std::vector<AtomReference>& residues);

//! Create the single bonds a flexfit file names between two named atoms.
/*! \param[in] atom_pairs two entries per bond, in order */
IMPBFFEXPORT IMP::ParticleIndexes create_named_bonds(
        IMP::atom::Hierarchy hierarchy,
        const std::vector<AtomReference>& atom_pairs);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_STRUCTUREIO_H
