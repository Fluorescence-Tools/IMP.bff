/**
 *  \file IMP/bff/RmfIO.h
 *  \brief Structures, rotamer libraries and trajectories through RMF.
 *
 * `rmf` is one of this module's modules (see `dependencies.py`), so these sit
 * beside the readers they belong with and the library they read is the same
 * #IMP::bff::RotamerLibrary every other reader returns.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_RMFIO_H
#define IMPBFF_RMFIO_H

#include <IMP/bff/bff_config.h>

// RMF is optional (PRD-139), and opt-*out*: every build has it unless one
// says otherwise, so IMP's own tooling -- which knows nothing of this flag --
// keeps building what it always built. It is one file and four functions,
// called from `bin/imp_bff` and nowhere in the library, and it costs a wheel
// dearly: librmf pulls Boost.Iostreams, which as packaged (conda and
// Homebrew alike) pulls Boost.Regex and all of ICU -- some 40 MB against
// 3 MB of library. Without RMF the four declarations below are not made, so
// a caller finds out at compile time rather than at link time.
#ifndef IMPBFF_NO_RMF
#include <IMP/bff/HierarchyFrame.h>
#include <IMP/bff/RotamerLibrary.h>
#include <IMP/bff/StructureTable.h>

#include <memory>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Write `(N, 3)` coordinates to an RMF file as one frame of spheres.
/*!
    Each row becomes a particle with an #IMP::core::XYZR and an
    #IMP::atom::Mass, under one hierarchy root -- enough for a viewer to open
    the file and for #IMP::rmf to read it back.

    \param[in] coords,n_atoms,n_dim the positions, `(N, 3)`
    \param[in] path the file to write
    \param[in] model_name names the hierarchy root
    \param[in] metadata_json recorded as the file's description; a caller with
               a dict passes `json.dumps` of it. Empty writes no description.
    \param[in] radius of every sphere, Angstrom
    \throw ValueException when \p n_dim is not 3
*/
IMPBFFEXPORT void write_rmf(double* coords, int n_atoms, int n_dim,
                            const std::string& path,
                            std::string model_name = "structure",
                            std::string metadata_json = "",
                            double radius = 1.5);

//! Write a rotamer library to an RMF file, one frame per conformer.
/*!
    The weight of a conformer is a per-frame value in the `score` category;
    the jump counts, when the library has them, are one static `Ints` value in
    the `kinetic` category. Both are what #IMP::bff::read_rotamer_library_rmf
    reads back.

    \param[in] path the file to write; `.rmf3` is appended when missing
    \param[in] library the conformers, their weights and their atom names
    \throw ValueException when the library has no conformers, or its
           coordinates do not match `n_rotamers * n_atoms * 3`
*/
IMPBFFEXPORT void write_rotamer_library_rmf(const std::string& path,
                                            const RotamerLibrary& library);

//! Read a rotamer library from an RMF file.
/*!
    \param[in] path the file; `.rmf3` is appended when the path as given is
               not there
    \return the library. `transitions` is empty when the file carries none,
            which is what a container of static conformers looks like.
    \throw IOException when neither \p path nor `path + ".rmf3"` exists
*/
IMPBFFEXPORT RotamerLibrary read_rotamer_library_rmf(const std::string& path);

//! Every frame of an RMF trajectory, as #IMP::bff::ProteinFrame values.
/*! The PDB side of this is #IMP::bff::load_protein_frames and the frame
    itself is #IMP::bff::protein_frame_from_hierarchy; what is here is opening
    an RMF and stepping it.

    \param[in] path the trajectory
    \param[in] max_frames stop after this many; negative reads all
    \throw IOException when \p path cannot be opened
*/
IMPBFFEXPORT std::vector<ProteinFrame> protein_frames_from_rmf(
        const std::string& path, int max_frames = -1);

//! A structure written as an RMF trajectory: one hierarchy, many frames.
/*!
    The shape is PMI's, and deliberately so -- root, one node per chain, one
    per residue, one per atom -- because that is what the viewers and the
    analysis scripts around this ecosystem open. It is built through **RMF's
    own decorators**, not through `IMP::rmf`: the connection layer's IMP does
    not carry that module, and a writer that needed it could not run in a
    package that links IMP privately.

    Per-frame scalars go into RMF's `stat` category, which is where PMI's
    `IMP.pmi.output.Output` puts them, so a stat plot drawn from a PMI run
    reads one of these unchanged. They arrive as a JSON object because that
    is what carries a set of names the writer cannot know in advance; a
    number becomes a float or int key by its JSON type, anything else a
    string. A key is created once, on the frame that first mentions it.

    \see #IMP::bff::write_rmf for the single-frame sphere writer, which is a
         different thing: it has no residues, no chains and no names.
*/
class IMPBFFEXPORT RmfStructureWriter {
    struct Impl;
    std::shared_ptr<Impl> impl_;

 public:
    // No default constructor: a writer without a file is not a useful object,
    // and a second constructor would make this one an *overload*, which is
    // where SWIG switches keyword arguments off -- silently, so the Python
    // signature quietly becomes positional-only.

    //! Open \p path and build the hierarchy \p structure describes.
    /*!
        \param[in] path the file to write; `.rmf3` is appended when missing
        \param[in] structure the atoms -- their names, residues and chains fix
                   the tree, their radii and masses the particle decorators.
                   The coordinates it carries are *not* written: a frame is
                   #append.
        \param[in] root_name names the hierarchy root
        \param[in] metadata_json the file's description, and the stat keys to
                   create up front. Empty writes no description.
        \throw ValueException when the table has no atoms, or its columns
               disagree on how many there are
        \throw IOException when \p path cannot be opened
    */
    RmfStructureWriter(const std::string& path, const StructureTable& structure,
                       const std::string& root_name = "structure",
                       const std::string& metadata_json = "");

    //! Append one frame.
    /*!
        \param[in] coords flat, three per atom, in the table's atom order
        \param[in] frame_name empty takes the frame index
        \param[in] metadata_json per-frame scalars for the `stat` category
        \throw ValueException when \p coords is not three per atom
    */
    void append(const std::vector<double>& coords,
                const std::string& frame_name = "",
                const std::string& metadata_json = "");

    int get_n_atoms() const;
    int get_number_of_frames() const;

    //! Flush and release the file. Destruction does this too.
    void close();

    IMP_SHOWABLE_INLINE(RmfStructureWriter,
                        out << "RmfStructureWriter(" << get_n_atoms()
                            << " atoms, " << get_number_of_frames()
                            << " frames)");
};

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_NO_RMF

#endif //IMPBFF_RMFIO_H
