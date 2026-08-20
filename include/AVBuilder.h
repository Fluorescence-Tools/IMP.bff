/**
 *  \file IMP/bff/AVBuilder.h
 *  \brief Building an accessible volume: the two front doors, and the PDB read.
 *
 * #IMP::bff::AV is the solver and #IMP::bff::AccessibleVolume is the result;
 * this is what stands between them. There are two front doors and they differ
 * only in what identifies the label:
 *
 * - #IMP::bff::compute_av takes **raw arrays** — atoms, radii and an attachment
 *   coordinate. No file, no fps position. It is what a restraint has.
 * - #IMP::bff::compute_av_from_structure takes a **PDB and an fps position
 *   definition**, applies the FPS strip and resolves the attachment atom by
 *   `(chain, residue, atom name)`.
 *
 * Both drive the same `AV`/`PathMap` core through #IMP::bff::resample_av, which
 * is the one place the decorator is set up and the map read. Sharing it is what
 * stops the two doors from each learning the same lessons: the array door had
 * found that IMP numbers voxels with *x* fastest, so a C-order reshape into
 * `(nx, ny, nz)` returns the volume **transposed**; the structure door had not,
 * and shipped a mirrored density for as long as it existed — 71 % of its own
 * point cloud landing on an occupied voxel instead of 100 %, and the mean donor
 * lifetime on T4L A132 out by 4.9 % (PRD-113 stage 3a). The structure door had
 * found that the AV must be decorated onto its **own** particle with the source
 * passed separately, or the resampled map sits at the coordinate origin.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_AVBUILDER_H
#define IMPBFF_AVBUILDER_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/AVModel.h>

#include <IMP/value_macros.h>
#include <IMP/showable_macros.h>

#include <IMP/Model.h>
#include <IMP/Particle.h>

#include <map>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! Radius around the attachment atom inside which obstacles are ignored.
/*!
    So the linker can leave the atom it is tied to. Below roughly 2 A the
    attachment atom's own neighbours block every starting voxel and the path
    search returns an **empty volume without raising**. 2.1 A is carbon's CHARMM
    radius, so it sits just clear of that cliff and matches the atom a linker is
    normally tied to.
*/
IMPBFFEXPORT extern const double DEFAULT_ALLOWED_SPHERE_RADIUS;

//! Van der Waals radius per atomic number, Angstrom; 1.70 for anything else.
IMPBFFEXPORT std::map<int, double> vdw_radii();

//! The element symbol of a PDB ATOM/HETATM record.
/*!
    Columns 77-78 carry it in a modern PDB and are used verbatim when present.
    Older files — the shipped FPS screening structures are 66-character records
    — leave them empty, so the element comes from the atom-name field, where it
    is *right-justified in columns 13-14*: a blank or numeric column 13 means a
    one-letter element, so `" CA "` is an alpha-carbon while `"CA  "` is
    calcium. A two-letter reading is additionally rejected when columns 15-16
    contain a digit, which is how a four-character hydrogen name such as
    `"HE21"` is written — helium would otherwise win.
*/
IMPBFFEXPORT std::string element_symbol_from_pdb_line(const std::string& line);

//! One ATOM/HETATM record: where it is, how big, and what names it.
struct IMPBFFEXPORT PDBAtomRecord {
    std::string chain;
    int resseq;
    std::string atom_name;
    double x, y, z;
    double vdw_radius;
    PDBAtomRecord() : resseq(0), x(0), y(0), z(0), vdw_radius(1.70) {}
    IMP_SHOWABLE_INLINE(PDBAtomRecord,
                        out << "PDBAtomRecord(" << chain << resseq << ":"
                            << atom_name << ")");
};
IMP_VALUES(PDBAtomRecord, PDBAtomRecords);

//! Every ATOM/HETATM record of a PDB, with a van der Waals radius each.
/*!
    Parsed directly rather than through `IMP::atom::read_pdb`, which warns per
    unsupported HETATM residue and there are thousands in a labelled structure.

    **Cached on (path, mtime, size)**: an fps.json with a dozen positions reads
    the same file a dozen times, and the strip below reads it again per site.

    \throw IOException when the file has no coordinates at all
*/
IMPBFFEXPORT std::vector<PDBAtomRecord> read_pdb_records(
        const std::string& pdb_path);

//! `(N, 4)` of x, y, z, vdW radius — the obstacle array the solver takes.
IMPBFFEXPORT void load_structure_with_vdw(const std::string& pdb_path,
                                          double** out_view, int* n_out_view);

//! The coordinates of an attachment atom, resolved by identity.
/*!
    A miss stays a miss: an unresolvable site returns an empty view rather than
    a positional guess, because a dye attached to an unrelated atom yields a
    plausible and entirely wrong accessible volume.

    \param[in] pdb_path the structure the site names
    \param[in] chain empty matches any chain
    \param[in] resseq,atom_name the site
    \param[out] out_view,n_out_view three coordinates, or zero-length
*/
IMPBFFEXPORT void find_attachment_point(const std::string& pdb_path,
                                        const std::string& chain, int resseq,
                                        const std::string& atom_name,
                                        double** out_view, int* n_out_view);

//! Decorate a fresh particle as an AV, resample it, and read the map.
/*!
    The AV is decorated onto its **own** particle with the source passed
    separately. Setting it up on the source particle leaves the resampled map at
    the coordinate origin — header origin (0, 0, 0), obstacles nowhere near the
    search region, and most of the grid reported accessible.

    \param[in] model holds the obstacle hierarchy and the attachment particle
    \param[in] source_particle the attachment atom, a real `XYZR` in that model
    \param[in] linker_length,linker_width,r1,r2,r3,disc_step the AV parameters
    \param[in] allowed_sphere_radius obstacles inside this radius of the
               attachment are ignored; without it the search starts inside the
               attachment atom's own neighbourhood and returns an empty volume,
               reporting nothing
    \param[in] contact_volume_thickness,contact_volume_trapped_fraction the ACV
               split, when one is wanted
    \param[in] search_stencil Dijkstra neighbour stencil. `0` leaves the
               decorator's own default, which is **74** — the LabelLib reference
               metric. `26` is the speed option: ~1.9x faster for ~20 % less
               volume, because a {1,2,3} stencil's isopath surface is cubic
               rather than spherical. `30` is the historical variant.
    \return the volume, carrying its cloud, its grid and the source coordinate
*/
IMPBFFEXPORT AccessibleVolume resample_av(
        IMP::Model* model, IMP::Particle* source_particle,
        double linker_length, double linker_width, double r1, double r2,
        double r3, double disc_step, double allowed_sphere_radius,
        double contact_volume_thickness = 0.0,
        double contact_volume_trapped_fraction = -1.0,
        int search_stencil = 0);

//! An accessible volume from raw atomic coordinates.
/*!
    The **array** front door: no PDB file and no fps position definition, which
    is what a restraint has.

    \param[in] atoms_xyzr,n_atoms,n_cols obstacles, `(N, 4)` of x, y, z, radius
    \param[in] source_xyz where the linker is tied, three coordinates
    \param[in] linker_length,linker_width the linker
    \param[in] r1,r2,r3 dye radii; `(3.5, 0, 0)` is the AV1 single-sphere model
    \param[in] grid_resolution voxel spacing, A
    \param[in] allowed_sphere_radius see #DEFAULT_ALLOWED_SPHERE_RADIUS
    \param[in] search_stencil see #IMP::bff::resample_av
*/
IMPBFFEXPORT AccessibleVolume compute_av(
        double* atoms_xyzr, int n_atoms, int n_cols,
        const std::vector<double>& source_xyz, double linker_length = 20.0,
        double linker_width = 0.5, double r1 = 3.5, double r2 = 0.0,
        double r3 = 0.0, double grid_resolution = 1.5,
        double allowed_sphere_radius = 2.1, int search_stencil = 0);

//! An accessible volume from a structure and an fps position definition.
/*!
    The **structure** front door. It applies the FPS strip before measuring
    anything: fps.json positions are calibrated for a structure whose attachment
    residue does not wall in its own dye. A declared `strip_mask` outside the
    fps dialect raises here — loudly, because the alternative is computing
    against obstacles the document said to remove.

    Source clearance scales with the linker: the path search inflates obstacles
    by half the linker width, so the free sphere has to clear that inflation
    plus a grid step of slack or the source tile is walled in and the volume
    comes back empty. A position that declares `allowed_sphere_radius` keeps its
    own value — honoured exactly, never escalated: an empty volume at the
    declared parameters is the answer, not a signal to retry at invented ones.

    \param[in] pdb_path the structure
    \param[in] chain,resseq,atom_name the attachment site
    \param[in] linker_length,linker_width,r1,r2,r3,disc_step the AV parameters
    \param[in] strip_mask an fps `strip_mask`; empty means the default strip
    \param[in] allowed_sphere_radius negative derives it from the linker width
    \param[in] contact_volume_thickness,contact_volume_trapped_fraction the ACV
               split, when one is wanted
    \throw ValueException when the attachment site is not in the structure
*/
IMPBFFEXPORT AccessibleVolume compute_av_from_structure(
        const std::string& pdb_path, const std::string& chain, int resseq,
        const std::string& atom_name, double linker_length, double linker_width,
        double r1, double r2, double r3, double disc_step = 1.5,
        const std::string& strip_mask = "",
        double allowed_sphere_radius = -1.0,
        double contact_volume_thickness = 0.0,
        double contact_volume_trapped_fraction = -1.0);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_AVBUILDER_H
