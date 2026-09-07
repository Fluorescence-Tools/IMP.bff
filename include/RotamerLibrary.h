/**
 *  \file IMP/bff/RotamerLibrary.h
 *  \brief One rotamer library value, whatever file it came from.
 *
 * There were three of these: `RotamerLibraryData` (the numpy/text pair, whose
 * reader and writer are now at the bottom of this file), `RotamerLibrary` (a PDB plus a trajectory, `ProbeSampling.h`) and
 * `DrotLibrary` (a `.drot` container, `DrotReader.h`). All three held the same
 * thing -- an ensemble of conformers, a weight each, and the atom names the
 * coordinates are in the order of -- and differed only in which of the
 * *optional* columns their reader happened to fill. A caller that took one
 * could not be handed another, so each reader grew its own consumers, and
 * `id` (which was always `1..n`) existed in one of them alone.
 *
 * One value, three readers: #IMP::bff::read_rotamer_library,
 * #IMP::bff::read_drot and #IMP::bff::load_rotamer_library_trajectory. A
 * column a format does not carry is left empty rather than invented, and
 * every consumer states which columns it needs.
 *
 * The `.drot` container -- the format this package writes its own libraries
 * in -- lives here too, reader and writer both. It was two files of its own
 * (`DrotReader.h`, `DrotWriter.h`) named after the extension rather than the
 * thing; but a `.drot` is not a kind of file this package happens to read,
 * it is how a #IMP::bff::RotamerLibrary is spelled on disk, and the value and
 * its serialisation are one subject.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_ROTAMERLIBRARY_H
#define IMPBFF_ROTAMERLIBRARY_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/PtoProfile.h>

#include <IMP/bff/Base.h>

#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! An ensemble of conformers with a weight each: a rotamer library.
/*! The conformers may be a dye's rotamers, a spin label's, a side chain's or
    the frames of a trajectory that was clustered into one -- the value says
    nothing about which, because nothing that reads it needs to know. */
struct IMPBFFEXPORT RotamerLibrary {
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

    RotamerLibrary() : n_rotamers(0), n_atoms(0) {}

    void get_coords(double** out_view, int* n_out_view) const;
    void get_weights(double** out_view, int* n_out_view) const;
    //! Replace the coordinates; `n_rotamers` and `n_atoms` say what shape
    //! they are, and this does not touch them.
    void set_coords(const std::vector<double>& v) { coords = v; }
    void set_weights(const std::vector<double>& v) { weights = v; }

    IMP_SHOWABLE_INLINE(RotamerLibrary,
                        out << "RotamerLibrary(" << n_rotamers << " x "
                            << n_atoms << " atoms)");
};
IMP_VALUES(RotamerLibrary, RotamerLibraries);

//! The same library with its weights summing to one.
/*! A copy, because the library is a value and there is no in-place spelling
    across the language boundary. Weights that sum to zero are left alone --
    there is no normalisation of nothing, and inventing a uniform one here
    would hide a library that failed to load. */
IMPBFFEXPORT RotamerLibrary normalize_weights(const RotamerLibrary& lib);

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
IMPBFFEXPORT RotamerLibrary read_drot(const std::string& path);

//! Split a locator into its container and its library (empty when there is
//! no `::`).
IMPBFFEXPORT std::vector<std::string> split_drot_locator(
        const std::string& locator);

//! The libraries a container holds, in the order they were written.
/*!
    A single-library container answers with one empty name -- the library is
    the file. A family container (`dyes.drot.pto`) answers with the names its
    `drot.catalog` lists, which are the names `read_drot(path, library)`
    takes. Only the catalog is read, so asking is cheap whatever the file
    weighs.
*/
IMPBFFEXPORT std::vector<std::string> drot_catalog(const std::string& path);

//! Read one library out of a container that holds several.
/*!
    Reads that library's objects and no others: the framing is walked by
    seeking, so pulling one library out of a family container costs its own
    size, not the file's.

    \throw IOException when the container has no such library
*/
IMPBFFEXPORT RotamerLibrary read_drot(const std::string& path,
                                   const std::string& library);

//! How a library was made, as the #MfdbTag pairs `write_drot` was given.
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
IMPBFFEXPORT MfdbTags drot_provenance(const std::string& path,
                                      const std::string& library = "");

// ---------------------------------------------------------------------------
// The `.drot` container: writing
// ---------------------------------------------------------------------------
//! How a `.drot` stores its numbers.
struct IMPBFFEXPORT DrotEncoding {
    //! float32 grids (the default): exact to ~1e-6 A.
    bool lossless;
    //! Compact rung only: step for base coordinates and bond lengths, A.
    double grid_a;
    //! Compact rung only: step for bond angles and dihedrals, degrees.
    double grid_deg;
    //! brotli quality, 0-11 (11 is what the shipped libraries use).
    int quality;
    //! brotli window bits, 10-24.
    int window;

    DrotEncoding()
        : lossless(true), grid_a(0.001), grid_deg(0.01), quality(11),
          window(24) {}

    IMP_SHOWABLE_INLINE(DrotEncoding,
                        out << "DrotEncoding("
                            << (lossless ? "lossless f32" : "compact i16")
                            << ", brotli q" << quality << ")");
};
IMP_VALUES(DrotEncoding, DrotEncodings);

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

    The object survives `write_drot_bundle`, which copies objects verbatim
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
    a property of the tar, not of bundling: `read_drot(path, library)` walks
    the framing and reads that library's objects alone.

    \param[in] sources the containers to bundle, each holding one library
    \param[in] names the library name to file each source under, parallel to
               \p sources; a name may not contain `/`
    \param[in] path the family container to write
    \throw ValueException when the two vectors disagree or a name is unusable
    \throw IOException when a source cannot be read or the file not written
*/
IMPBFFEXPORT void write_drot_bundle(const std::vector<std::string>& sources,
                                    const std::vector<std::string>& names,
                                    const std::string& path);

IMPBFFEXPORT void write_drot(const std::string& path,
                             double* rotamer_coords, int n_rotamer_coords,
                             const std::vector<std::string>& atom_names,
                             const std::vector<std::string>& elements,
                             const std::vector<std::string>& resnames,
                             double* rotamer_weights, int n_rotamer_weights,
                             const DrotEncoding& encoding = DrotEncoding());

//! `write_drot`, and say how the ensemble was made.
/*!
    A separate name rather than a defaulted parameter on `write_drot`, because
    the default does not survive SWIG. `write_drot` takes its
    weights through a `(double*, int)` typemap that consumes one Python
    argument, and SWIG's per-arity wrappers mis-bind around it: with the
    provenance argument defaulted, `write_drot(..., encoding)` aborted the
    interpreter rather than running. `%feature("compactdefaultargs")` did not
    help. Two names with no defaults generate two independent wrappers and no
    dispatch, which is worth more than the tidier signature.

    \param[in] provenance #MfdbTag pairs from `mfdb_operation_tags` and
               friends, written as a `drot.provenance` object

    Everything else is `write_drot`, including every throw.
*/
IMPBFFEXPORT void write_drot_with_provenance(
        const std::string& path,
        double* rotamer_coords, int n_rotamer_coords,
        const std::vector<std::string>& atom_names,
        const std::vector<std::string>& elements,
        const std::vector<std::string>& resnames,
        double* rotamer_weights, int n_rotamer_weights,
        const DrotEncoding& encoding, const MfdbTags& provenance);

// --------------------------------------------------------------------------


//! Read a rotamer library from numpy/text files.
IMPBFFEXPORT RotamerLibrary read_rotamer_library(const std::string& path);

// --------------------------------------------------------------------------
// The backbone-dependent side-chain library
// --------------------------------------------------------------------------

/*! The second payload the container carries, and a different kind of thing
    from everything above. A probe or a spin label is an **ensemble**: a list
    of conformers with a weight each, which is what the `.drot` grid stores. A
    side-chain library is a **distribution over backbone conformation** -- a
    table indexed by amino acid and by the (phi, psi) bin the residue sits in,
    holding for each rotamer a probability, its chi angles and their standard
    deviations. Same envelope, different payload, which is what `PtoKind` is
    for: `rot.bbdep.records` sits beside `drot.grid` in one grammar.

    The table is the Dunbrack-2010 backbone-dependent library as FASPR ships
    it (`dun2010bbdep.bin`; Shapovalov & Dunbrack, Structure 2011;19:844-858,
    free for academic use -- the notice travels inside the container, in
    `bbdep.json`). **The conversion keeps FASPR's own 20-byte records byte for
    byte**: a record is `float32 probability` then eight `int16` in tenths of a
    degree (chi1..chi4, then their sigmas), and one amino acid's records are
    1296 (phi, psi) bins x its rotamer count, contiguous. Nothing is
    re-derived, so #IMP::bff::write_dunbrack_bin reproduces the input file
    exactly -- which is both the correctness proof and what lets the vendored
    FASPR engine, which seeks in that file by byte offset, run against a
    shipped container.

    \note Querying one is also #IMP::rotamer::RotamerLibrary's job, and that is
    what #IMP::bff::get_anchor_cb_position uses. The two read different files:
    IMP's reads the Dunbrack *text* library, this reads the binary FASPR ships
    and the container that carries it. Neither can read the other's, which is
    the whole reason both are here. */

//! The rotamers of one residue type at one backbone conformation.
struct IMPBFFEXPORT DunbrackRotamers {
    //! `n_rotamers * n_chi` chi angles in degrees, rotamer-major.
    std::vector<double> chi;
    //! The same shape: each chi's standard deviation, degrees.
    std::vector<double> sigma;
    //! One per rotamer, as the library stores it (they sum to ~1 per bin).
    std::vector<double> probability;
    int n_rotamers, n_chi;

    DunbrackRotamers() : n_rotamers(0), n_chi(0) {}

    void get_chi(double** out_view, int* n_out_view) const;
    void get_sigma(double** out_view, int* n_out_view) const;
    void get_probability(double** out_view, int* n_out_view) const;

    IMP_SHOWABLE_INLINE(DunbrackRotamers,
                        out << "DunbrackRotamers(" << n_rotamers << " x "
                            << n_chi << " chi)");
};
IMP_VALUES(DunbrackRotamers, DunbrackRotamersList);

//! Convert FASPR's `dun2010bbdep.bin` into a `.drot.pto` container.
/*!
    Writes one `rot.bbdep.records` object per amino acid -- that residue's
    records exactly as they lie in the binary -- plus a `rot.bbdep.header`
    describing the bin scheme, the record layout and the library's
    provenance, and a `drot.catalog` naming the residues. Each object is
    brotli'd on its own, which takes the 13.76 MB binary to about 3.3 MB.

    \param[in] bin_path FASPR's `dun2010bbdep.bin`
    \param[in] path the container to write
    \throw IOException when the binary is missing or is not the expected size
*/
IMPBFFEXPORT void write_dunbrack_library(const std::string& bin_path,
                                         const std::string& path);

//! Write the library back out as `dun2010bbdep.bin`.
/*!
    The inverse of `write_dunbrack_library`, byte for byte. Two things need
    it: proving the conversion lost nothing, and handing the vendored FASPR
    engine the file it insists on opening by name.
*/
IMPBFFEXPORT void write_dunbrack_bin(const std::string& path,
                                     const std::string& bin_path);

//! The rotamers of `residue` at backbone angles (`phi`, `psi`).
/*!
    Reads that residue's object and no other. The two cuts are FASPR's own
    defaults, applied in FASPR's order: rotamers come sorted by probability,
    the first one below \p probability_min ends the list, and the list also
    ends once the accumulated probability passes \p probability_accumulated.
    Pass 0 and 1 to take the bin whole.

    \param[in] path a container written by `write_dunbrack_library`
    \param[in] residue one-letter code; `A` and `G` have no rotamers
    \param[in] phi,psi degrees
    \throw ValueException for a residue this library does not carry
*/
IMPBFFEXPORT DunbrackRotamers read_dunbrack_rotamers(
        const std::string& path, char residue, double phi, double psi,
        double probability_min = 0.01,
        double probability_accumulated = 0.97);

// --------------------------------------------------------------------------
// Side-chain packing over that library
// --------------------------------------------------------------------------

//! Repack a protein's side chains to the global-minimum-energy conformation.
/*!
    Reads \p pdb_in (backbone N, CA, C, O required per residue; any side chains
    present are ignored), repacks every rotatable side chain, and writes the
    result to \p pdb_out with the residue numbering and ordering preserved.
    The search is deterministic: two calls on the same input and library give
    identical coordinates.

    \param[in] pdb_in the input backbone PDB
    \param[in] pdb_out where the repacked PDB is written
    \param[in] rotamer_library the Dunbrack-2010 library -- either FASPR's
               binary `dun2010bbdep.bin` or a `.drot.pto` container carrying
               it, which is unpacked to a scratch copy because the engine
               seeks in that file by byte offset
    \param[in] verbose echo the packing log to stdout

    \throw IOException when the input or the library cannot be read
    \throw std::runtime_error on a fatal packing error -- an incomplete
           backbone, an unknown library layout

    The algorithm is FASPR (Huang, Pearce & Zhang, *Bioinformatics*
    2020;36:3758-3765, MIT, github.com/tommyhuangthu/FASPR): CHARMM19 van der
    Waals, hydrogen bonds, disulfides and a rotamer prior, minimised by
    DEE-Goldstein and DEE-split elimination and then tree decomposition. It is
    vendored in `src/Faspr*.{h,cpp}`, which carry the original notice, and
    pinned against the reference executable by
    `test/faspr/test_faspr_port.py`.
*/
IMPBFFEXPORT void faspr_pack(const std::string& pdb_in,
                             const std::string& pdb_out,
                             const std::string& rotamer_library,
                             bool verbose = false);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_ROTAMERLIBRARY_H
