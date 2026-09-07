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
#include <IMP/bff/HierarchyFrame.h>
#include <IMP/bff/RotamerLibrary.h>

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

IMPBFF_END_NAMESPACE

#endif //IMPBFF_RMFIO_H
