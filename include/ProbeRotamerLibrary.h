/**
 *  \file IMP/bff/ProbeRotamerLibrary.h
 *  \brief One rotamer library value, whatever file it came from.
 *
 * There were three of these: `RotamerLibraryData` (the numpy/text pair, whose
 * reader and writer are now at the bottom of this file), `ProbeRotamerLibrary` (a PDB plus a trajectory, `ProbeSampling.h`) and
 * `DrotLibrary` (a `.drot` container, `DrotReader.h`). All three held the same
 * thing -- an ensemble of conformers, a weight each, and the atom names the
 * coordinates are in the order of -- and differed only in which of the
 * *optional* columns their reader happened to fill. A caller that took one
 * could not be handed another, so each reader grew its own consumers, and
 * `id` (which was always `1..n`) existed in one of them alone.
 *
 * One value, three readers: #IMP::bff::read_probe_rotamer_library,
 * #IMP::bff::read_probe_rotamer_drot and #IMP::bff::load_rotamer_library_trajectory. A
 * column a format does not carry is left empty rather than invented, and
 * every consumer states which columns it needs.
 *
 * The `.drot` container -- the format this package writes its own libraries
 * in -- lives here too, reader and writer both. It was two files of its own
 * (`DrotReader.h`, `DrotWriter.h`) named after the extension rather than the
 * thing; but a `.drot` is not a kind of file this package happens to read,
 * it is how a #IMP::bff::ProbeRotamerLibrary is spelled on disk, and the value and
 * its serialisation are one subject.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_PROBEROTAMERLIBRARY_H
#define IMPBFF_PROBEROTAMERLIBRARY_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/MMFDBProfile.h>

#include <IMP/bff/IMPCompatibility.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! An ensemble of conformers with a weight each: a rotamer library.
/*! The conformers may be a dye's rotamers, a spin label's, a side chain's or
    the frames of a trajectory that was clustered into one -- the value says
    nothing about which, because nothing that reads it needs to know. */
struct IMPBFFEXPORT ProbeRotamerLibrary {
    //! Flat, `n_rotamers * n_atoms * 3`, in the library's own frame.
    std::vector<double> coords;
    //! One per conformer. Readers that can normalise do; `.drot` keeps the
    //! populations its encoder found, so check rather than assume.
    std::vector<double> weights;
    //! One per atom, the order `coords` is in.
    std::vector<std::string> atom_names;
    //! Residue name per atom; empty when the format carries none.
    /*! Atom names repeat between a dye residue and its linker, so a selector
        that must tell `C13` of `A48` from `C13` of `C1R` needs this. */
    std::vector<std::string> resnames;
    //! Element per atom; empty when the format carries none.
    std::vector<std::string> elements;
    //! The registry entry for this library, as JSON, or `{}` when it has none.
    /*! The dye's selectors -- which atoms are the chromophore centre, which
        two span the transition dipole, which carry charge -- are here and not
        in the coordinates, so a library read from an off-registry path has an
        empty object and its caller supplies them. */
    std::string metadata;
    //! What it was read from: a path, or a `container::library` locator.
    std::string path;
    //! Jump counts between conformers, flat `n_rotamers * n_rotamers`.
    /*! Empty unless the library came from a run that saw the conformers in
        *order* -- a Metropolis trajectory, a molecular-dynamics one -- because
        only then is there such a thing as a jump. A library read from a
        container of static conformers has none, and that is not a defect. */
    std::vector<int> transitions;
    int n_rotamers, n_atoms;

    ProbeRotamerLibrary() : n_rotamers(0), n_atoms(0) {}

    void get_coords(double** out_view, int* n_out_view) const;
    void get_weights(double** out_view, int* n_out_view) const;
    //! Replace the coordinates; `n_rotamers` and `n_atoms` say what shape
    //! they are, and this does not touch them.
    void set_coords(const std::vector<double>& v) { coords = v; }
    void set_weights(const std::vector<double>& v) { weights = v; }

    IMP_SHOWABLE_INLINE(ProbeRotamerLibrary,
                        out << "RotamerLibrary(" << n_rotamers << " x "
                            << n_atoms << " atoms)");
};
IMP_VALUES(ProbeRotamerLibrary, ProbeRotamerLibraries);

//! The same library with its weights summing to one.
/*! A copy, because the library is a value and there is no in-place spelling
    across the language boundary. Weights that sum to zero are left alone --
    there is no normalisation of nothing, and inventing a uniform one here
    would hide a library that failed to load. */
IMPBFFEXPORT ProbeRotamerLibrary normalize_probe_rotamer_weights(const ProbeRotamerLibrary& lib);

//! Weights summing to one, in place.
/*! A set that sums to nothing becomes **uniform**, not zero: a library whose
    weights were never filled in is one where every conformer counts the same,
    and leaving them at zero would silently drop it out of any average. */
IMPBFFEXPORT void normalize_weights_in_place(std::vector<double>& weights);


// ---------------------------------------------------------------------------
// The `.drot` container: reading
// ---------------------------------------------------------------------------
//! Read a `.drot` (v5/v7/v8) rotamer library and reconstruct every conformer.
/*!
    The whole library is rebuilt from internal coordinates: bond lengths and
    angles from the embedded template, per-conformer dihedrals and base-atom
    coordinates from the grids -- the same construction the encoder measured
    exact to ~0.014 A RMSD against the source trajectories over the shipped
    corpus (transition-dipole error ~0.1 deg median).

    \throw IOException when the file is not a .drot, a member is missing or
           inconsistent with the header, or the brotli stream is corrupt.
*/
//! Read a library, addressed by path or by locator.
/*!
    A **locator** is how one library inside a family container is named:
    `dyes.drot.pto::A48_C1R_cutoff10` -- the container, `::`, the library.
    A plain path with no `::` means the container holds one library and that
    is the one wanted. Everything that hands a rotamer library around in this
    package -- the resolver, the loader, the reference-file finder -- passes
    one of these two, so a caller never has to know which shape it has.
*/
IMPBFFEXPORT ProbeRotamerLibrary read_probe_rotamer_drot(const std::string& path);

//! Split a locator into its container and its library (empty when there is
//! no `::`).
IMPBFFEXPORT std::vector<std::string> split_probe_rotamer_drot_locator(
        const std::string& locator);

//! The libraries a container holds, in the order they were written.
/*!
    A single-library container answers with one empty name -- the library is
    the file. A family container (`dyes.drot.pto`) answers with the names its
    `drot.catalog` lists, which are the names `read_probe_rotamer_drot(path, library)`
    takes. Only the catalog is read, so asking is cheap whatever the file
    weighs.
*/
IMPBFFEXPORT std::vector<std::string> probe_rotamer_drot_catalog(const std::string& path);

//! Read one library out of a container that holds several.
/*!
    Reads that library's objects and no others: the framing is walked by
    seeking, so pulling one library out of a family container costs its own
    size, not the file's.

    \throw IOException when the container has no such library
*/
IMPBFFEXPORT ProbeRotamerLibrary read_probe_rotamer_drot(const std::string& path,
                                   const std::string& library);

//! How a library was made, as the #MfdbTag pairs `write_probe_rotamer_drot` was given.
/*!
    Empty when the container carries no `drot.provenance` object, which every
    library written before provenance existed does -- absence is "not stated",
    never "made the default way".

    \param[in] path the container
    \param[in] library which library in a family container; empty for a
               single-library container

    The tag worth reading first is `_mmfdb_operation.settings_hash`: two
    libraries of the same probe that disagree are distinguishable by it without
    a reader having to know what any particular setting means.
*/
IMPBFFEXPORT MfdbTags probe_rotamer_drot_provenance(const std::string& path,
                                      const std::string& library = "");

// ---------------------------------------------------------------------------
// The `.drot` container: writing
// ---------------------------------------------------------------------------
//! How a `.drot` stores its numbers.
struct IMPBFFEXPORT ProbeRotamerDrotEncoding {
    //! float32 grids (the default): exact to ~1e-6 A.
    bool lossless;
    //! Compact rung only: step for base coordinates and bond lengths, A.
    double grid_a;
    //! Compact rung only: step for bond angles and dihedrals, degrees.
    double grid_deg;
    //! brotli quality, 0-11 (11 is what the shipped libraries use). The
    //! window is the codec's, the largest standard one (24) at every
    //! quality -- what the shipped libraries are pinned against.
    int quality;

    ProbeRotamerDrotEncoding()
        : lossless(true), grid_a(0.001), grid_deg(0.01), quality(11) {}

    IMP_SHOWABLE_INLINE(ProbeRotamerDrotEncoding,
                        out << "DrotEncoding("
                            << (lossless ? "lossless f32" : "compact i16")
                            << ", brotli q" << quality << ")");
};
IMP_VALUES(ProbeRotamerDrotEncoding, ProbeRotamerDrotEncodings);

//! Write conformers and weights as a `.drot` library.
/*!
    The file carries everything a reader needs, so nothing has to sit beside
    it: the atom names, elements and (when given) residue names travel in the
    template, and the weights travel as written -- normalisation is the
    loader's business, not the format's.

    The Z-matrix is derived from the first conformer: bonds are perceived
    geometrically (\c perceive_bonds), the spanning tree is rooted at `CA`
    when the names have one and at the first atom otherwise, and the tree
    rides in the file (`rows.json`), so a reader never re-derives it.

    \param[in] path output file
    \param[in] rotamer_coords,n_rotamer_coords flat `(n_frames, n_atoms, 3)`
               in Angstrom; `n_atoms` is the length of \p atom_names
    \param[in] atom_names one per atom, the order of the coordinates
    \param[in] elements element symbols parallel to \p atom_names, used for
               bond perception
    \param[in] resnames residue names parallel to \p atom_names; may be empty,
               and then atom selectors that name a residue cannot be resolved
               from the file alone
    \param[in] rotamer_weights,n_rotamer_weights one per conformer; empty
               means all conformers weigh 1
    \param[in] encoding lossless float32 (default) or the compact int16 rung
    \param[in] provenance how this ensemble was made, as #MfdbTag pairs from
               `mfdb_operation_tags` and friends; written as a
               `drot.provenance` object and empty by default

    A library that does not say how it was made is a number without a unit.
    These libraries are the case in point: two of the same probe can differ in
    conformer count (727 against 711 for A48_C1R) purely because one was
    clustered with a periodic dihedral metric and the other was not, and
    nothing in the coordinates says which. `mfdb_operation_tags` carries the
    algorithm and a `settings_hash`, so a reader can tell two libraries apart
    without parsing anyone's settings keys.

    The object survives `write_probe_rotamer_drot_bundle`, which copies objects verbatim
    under a `<library>/` prefix, so a family container carries the provenance
    of each library it holds rather than one statement for the bundle.

    \throw ValueException when the shapes disagree, when the molecule is not
           connected, or when a compact-rung value does not fit its grid
    \throw IOException when the file cannot be written
*/
//! Bundle single-library containers into one family container.
/*!
    The rotamer libraries ship one file per *family* -- every dye library in
    `dyes.drot.pto`, every spin label in `spinlabels.drot.pto` -- rather than
    one file per library. That is a container operation and nothing more: each
    source's objects are copied across verbatim, payload bytes untouched,
    under names prefixed `<library>/`, and a `drot.catalog` object at the head
    lists what is inside. Nothing is re-encoded, so a bundle is bit-identical
    to its parts and building one costs a copy.

    It works because PTO addresses objects individually. The old objection to
    one archive -- that reading one library means inflating all of them -- was
    a property of the tar, not of bundling: `read_probe_rotamer_drot(path, library)` walks
    the framing and reads that library's objects alone.

    \param[in] sources the containers to bundle, each holding one library
    \param[in] names the library name to file each source under, parallel to
               \p sources; a name may not contain `/`
    \param[in] path the family container to write
    \throw ValueException when the two vectors disagree or a name is unusable
    \throw IOException when a source cannot be read or the file not written
*/
IMPBFFEXPORT void write_probe_rotamer_drot_bundle(const std::vector<std::string>& sources,
                                    const std::vector<std::string>& names,
                                    const std::string& path);

IMPBFFEXPORT void write_probe_rotamer_drot(const std::string& path,
                             double* rotamer_coords, int n_rotamer_coords,
                             const std::vector<std::string>& atom_names,
                             const std::vector<std::string>& elements,
                             const std::vector<std::string>& resnames,
                             double* rotamer_weights, int n_rotamer_weights,
                             const ProbeRotamerDrotEncoding& encoding = ProbeRotamerDrotEncoding());

//! `write_probe_rotamer_drot`, and say how the ensemble was made.
/*!
    A separate name rather than a defaulted parameter on `write_probe_rotamer_drot`, because
    the default does not survive SWIG. `write_probe_rotamer_drot` takes its
    weights through a `(double*, int)` typemap that consumes one Python
    argument, and SWIG's per-arity wrappers mis-bind around it: with the
    provenance argument defaulted, `write_probe_rotamer_drot(..., encoding)` aborted the
    interpreter rather than running. `%feature("compactdefaultargs")` did not
    help. Two names with no defaults generate two independent wrappers and no
    dispatch, which is worth more than the tidier signature.

    \param[in] provenance #MfdbTag pairs from `mfdb_operation_tags` and
               friends, written as a `drot.provenance` object

    Everything else is `write_probe_rotamer_drot`, including every throw.
*/
IMPBFFEXPORT void write_probe_rotamer_drot_with_provenance(
        const std::string& path,
        double* rotamer_coords, int n_rotamer_coords,
        const std::vector<std::string>& atom_names,
        const std::vector<std::string>& elements,
        const std::vector<std::string>& resnames,
        double* rotamer_weights, int n_rotamer_weights,
        const ProbeRotamerDrotEncoding& encoding, const MfdbTags& provenance);

// --------------------------------------------------------------------------


//! Read a rotamer library from numpy/text files.
IMPBFFEXPORT ProbeRotamerLibrary read_probe_rotamer_library(const std::string& path);

IMPBFF_END_NAMESPACE

#endif // IMPBFF_PROBEROTAMERLIBRARY_H
