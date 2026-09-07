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
 * RMF is **not** here: it is in #IMP::bff::RmfIO.h, beside the rotamer-library
 * pair that shares its file format. `rmf` is one of this module's
 * `required_modules`, so those readers and writers are C++ like these.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_STRUCTUREIO_H
#define IMPBFF_STRUCTUREIO_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/ZMatrix.h>

#include <IMP/atom/bond_decorators.h>
#include <IMP/bff/AVBuilder.h>

#include <IMP/Model.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/bff/Base.h>

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

    The rule is #perceive_bonds': bonded when the separation is under
    \f$r_i + r_j + \mathrm{tolerance}\f$ over \c covalent_radius. This used to
    be a second table (Cordero radii) under a multiplicative rule, which agreed
    with the first on every molecule in the tree but not on the cutoffs -- C-H
    at 1.31 A against 1.49 A -- so a stretched bond was a bond or not depending
    on which function was asked.

    \param[in] tolerance added to the sum of the two covalent radii */
IMPBFFEXPORT std::vector<AtomBond> infer_bonds(
        const std::vector<PDBAtomRecord>& atoms,
        double tolerance = 0.35);

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

// --------------------------------------------------------------------------
// point clouds and grids out, in the formats a viewer opens
// --------------------------------------------------------------------------

//! Write a weighted point cloud as an XYZ file.
/*! One line per point, `element x y z weight`. Every molecular viewer reads
    XYZ, and the weight rides in a fifth column that viewers ignore and numpy
    does not.

    \param[in] path where to write
    \param[in] points flat, four per point: x, y, z, weight
    \param[in] element the element symbol every point is given
    \param[in] comment the second line of the file
    \throw IOException when \p path cannot be opened
*/
IMPBFFEXPORT void write_points_xyz(const std::string& path,
                                   const std::vector<double>& points,
                                   const std::string& element = "He",
                                   const std::string& comment = "");

//! Write a weighted point cloud as a PQR file.
/*! PDB `ATOM` records with the weight in the **charge** column and \p radius
    in the radius column, which is what PQR is: a PDB whose last two columns
    are occupancy-shaped but mean charge and radius. PyMOL and VMD both read
    it, and colouring by charge then colours by weight.

    \param[in] path where to write
    \param[in] points flat, four per point
    \param[in] radius the radius written for every point, A
    \throw IOException when \p path cannot be opened
*/
IMPBFFEXPORT void write_points_pqr(const std::string& path,
                                   const std::vector<double>& points,
                                   double radius = 1.0);

//! Write a scalar grid as an OpenDX map.
/*! The volumetric format PyMOL (`load map.dx`), VMD and APBS all read. Values
    run in **x-fastest** order, which is the order OpenDX specifies and the
    opposite of the C order a numpy `(nx, ny, nz)` array iterates in -- the
    transposition is done here so a caller never has to know.

    \param[in] path where to write
    \param[in] density the grid, `nx * ny * nz`, in C order
    \param[in] nx,ny,nz the grid shape
    \param[in] origin three coordinates of the lowest corner, A
    \param[in] spacing the grid step, A
    \throw ValueException on a shape mismatch
    \throw IOException when \p path cannot be opened
*/
IMPBFFEXPORT void write_opendx(const std::string& path,
                               const std::vector<double>& density,
                               int nx, int ny, int nz,
                               const std::vector<double>& origin,
                               double spacing);

//! Convert a dye PDB to mmCIF with `_atom_site` records.
/*! \param[in] probe_id the data block name; empty takes the PDB's stem */
IMPBFFEXPORT void convert_pdb_to_cif(const std::string& pdb_path,
                                     const std::string& cif_path,
                                     const std::string& probe_id = "");

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

//! A PDB read into a model: the hierarchy, its leaves and their coordinates.
class IMPBFFEXPORT LoadedStructure {
    IMP::atom::Hierarchy hierarchy_;
    IMP::ParticlesTemp leaves_;
    std::vector<double> coords_;

 public:
    LoadedStructure() {}
    LoadedStructure(IMP::atom::Hierarchy h, const IMP::ParticlesTemp& leaves,
                    const std::vector<double>& coords)
        : hierarchy_(h), leaves_(leaves), coords_(coords) {}

    //! By value: IMP's value types may not be handed out by non-const ref.
    IMP::atom::Hierarchy get_hierarchy() const { return hierarchy_; }
    //! The leaf particles, in hierarchy order -- three coordinates each.
    IMP::ParticlesTemp get_leaves() const { return leaves_; }
    int get_number_of_leaves() const {
        return static_cast<int>(leaves_.size());
    }
    //! `(n_leaves, 3)` in Angstrom.
    void get_coords(double** output, int* n_output1, int* n_output2) const;
};

//! Read a PDB into \p model and return its hierarchy, leaves and coordinates.
/*! The model is the caller's because it owns the particles: drop it and the
    hierarchy decorates nothing. */
IMPBFFEXPORT LoadedStructure load_structure_with_particles(
        const std::string& path, IMP::Model* model);

//! What a FlexFit block names: flexible residues and the bonds to add.
class IMPBFFEXPORT FlexFitSelection {
    IMP::ParticlesTemp residues_;
    IMP::atom::Bonds bonds_;

 public:
    FlexFitSelection() {}
    FlexFitSelection(const IMP::ParticlesTemp& residues,
                     const IMP::atom::Bonds& bonds)
        : residues_(residues), bonds_(bonds) {}

    //! One particle per flexible residue, in the order the block lists them.
    IMP::ParticlesTemp get_residues() const { return residues_; }
    //! The bonds the block declares between named atoms, created in the model.
    IMP::atom::Bonds get_bonds() const { return bonds_; }
};

//! Resolve one FlexFit block against a hierarchy.
/*! \param[in] hier the structure the names resolve against
    \param[in] flexfit_json one FlexFit set as JSON -- the `"Flexible
               residues"` and `"Bonds"` arrays of an angle file
    \throw ValueException when the JSON does not parse or lacks those keys */
IMPBFFEXPORT FlexFitSelection read_angle_file(IMP::atom::Hierarchy hier,
                                              const std::string& flexfit_json);

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
