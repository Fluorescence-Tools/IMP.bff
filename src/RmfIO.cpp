/**
 * \file RmfIO.cpp
 * \brief Structures, rotamer libraries and trajectories through RMF.
 *
 * Written against **RMF's own API**, not `IMP.rmf`. RMF is a separate library
 * that IMP happens to vendor (`modules/rmf/dependency/RMF`), and everything
 * this file needs of it -- `ParticleFactory` for coordinates, radius and mass,
 * `ChainFactory` and `ResidueFactory` for the two levels above -- is RMF's.
 * Going through `IMP.rmf` bought convenience and cost the module a dependency
 * on `IMP.rmf`, and through it `isd`, `saxs`, `em` and `statistics`: five IMP
 * modules for one file.
 *
 * The on-disk layout is unchanged, because those RMF decorators *are* the
 * format -- `IMP.rmf` only ever mapped IMP's decorators onto them. Files this
 * writes are read by `IMP.rmf` and vice versa, and the round trip is a test
 * rather than a hope: `test/io/test_rmf_io.py` checks the shipped `.rmf3`
 * rotamer templates read back with the coordinates and weights they had.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/RmfIO.h>
#include <IMP/bff/internal/Text.h>

#include <RMF/FileConstHandle.h>
#include <RMF/FileHandle.h>
#include <RMF/Nullable.h>
#include <RMF/decorator/physics.h>
#include <RMF/decorator/sequence.h>
#include <RMF/infrastructure_macros.h>

#include <cmath>
#include <sstream>

IMPBFF_BEGIN_NAMESPACE

namespace {

//! The library's two extra channels, spelled once for reader and writer.
const char* const kWeightCategory = "score";
const char* const kWeightKey = "weight";
const char* const kKineticCategory = "kinetic";
const char* const kTransitionKey = "transition_matrix";

//! \p path, with `.rmf3` appended when it does not end in it.
std::string with_rmf3(const std::string& path) {
    return internal::ends_with(path, ".rmf3") ? path : path + ".rmf3";
}

//! Compose two rotations, parent then child, as quaternions (w, x, y, z).
RMF::Vector4 quat_mul(const RMF::Vector4& a, const RMF::Vector4& b) {
    return RMF::Vector4(
            a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
            a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
            a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
            a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0]);
}

//! Rotate \p v by the quaternion \p q.
RMF::Vector3 quat_rot(const RMF::Vector4& q, const RMF::Vector3& v) {
    const double w = q[0], x = q[1], y = q[2], z = q[3];
    const double m00 = 1 - 2 * (y * y + z * z), m01 = 2 * (x * y - z * w),
                 m02 = 2 * (x * z + y * w);
    const double m10 = 2 * (x * y + z * w), m11 = 1 - 2 * (x * x + z * z),
                 m12 = 2 * (y * z - x * w);
    const double m20 = 2 * (x * z - y * w), m21 = 2 * (y * z + x * w),
                 m22 = 1 - 2 * (x * x + y * y);
    return RMF::Vector3(
            static_cast<float>(m00 * v[0] + m01 * v[1] + m02 * v[2]),
            static_cast<float>(m10 * v[0] + m11 * v[1] + m12 * v[2]),
            static_cast<float>(m20 * v[0] + m21 * v[1] + m22 * v[2]));
}

//! A leaf, with the residue and chain it was found under.
struct LeafInfo {
    RMF::NodeConstHandle node;
    int residue_index;
    std::string resname;
    std::string chain_id;
    //! The reference-frame nodes above this leaf, outermost first.
    /*! The *nodes*, not their composition: a reference frame is frame data, so
        in a trajectory the transforms move while the tree does not. Composing
        once and reusing it put every atom of the T4L docking trajectory up to
        43 A from where IMP had it. The chain is found once; the composition is
        redone per frame from these nodes' current values. */
    std::vector<RMF::NodeConstHandle> frames;
};

//! Every representation leaf under \p n, depth first, with its ancestry.
/*! The order matters and is the whole reason this is written out rather than
    inlined: it has to be the order `IMP::atom::get_leaves` produced, because
    that is the order the atom names, coordinates and residue indices of every
    already-written file are in. Depth-first over children, taking nodes with
    none, is exactly what `get_leaves` did over the hierarchy `IMP.rmf` built
    from this same tree.

    Residue and chain are carried *down* rather than looked up afterwards:
    RMF stores them on the residue and chain nodes, a leaf carries neither, and
    `NodeConstHandle` has no parent link to walk back up. Threading them
    through the recursion costs nothing and is the only direction available. */
void collect_leaves(RMF::NodeConstHandle n,
                    RMF::decorator::ResidueFactory& rf,
                    RMF::decorator::ChainFactory& cf,
                    RMF::decorator::ReferenceFrameFactory& ff,
                    int residue_index, std::string resname,
                    std::string chain_id,
                    std::vector<RMF::NodeConstHandle> frames,
                    std::vector<LeafInfo>& out) {
    if (ff.get_is(n)) frames.push_back(n);
    if (rf.get_is(n)) {
        residue_index = rf.get(n).get_residue_index();
        resname = rf.get(n).get_residue_type();
    }
    if (cf.get_is(n)) chain_id = cf.get(n).get_chain_id();

    const RMF::NodeConstHandles children = n.get_children();
    if (children.empty()) {
        // Only representations: RMF's decorators usage-check the node type,
        // and a file may carry BOND or ORGANIZATIONAL leaves that are not
        // particles and were never in `get_leaves`'s answer either.
        if (n.get_type() == RMF::REPRESENTATION) {
            LeafInfo li;
            li.node = n;
            li.residue_index = residue_index;
            li.resname = resname;
            li.chain_id = chain_id;
            li.frames = frames;
            out.push_back(li);
        }
        return;
    }
    for (std::size_t i = 0; i < children.size(); ++i) {
        collect_leaves(children[i], rf, cf, ff, residue_index, resname,
                       chain_id, frames, out);
    }
}

//! A leaf's coordinates in the file's frame, not its own.
/*! Composed from the leaf's reference-frame ancestors **at the current frame**,
    outermost first: translate by what is accumulated so far, then fold in that
    frame's rotation. This is what `IMP.rmf` did when it turned these nodes into
    rigid bodies, and it is checked against IMP's own answer rather than
    reasoned about -- `test/io/test_rmf_io.py` compares all 159 atoms of the
    T4L docking trajectory across frames. */
RMF::Vector3 leaf_coordinates(const LeafInfo& li,
                              RMF::decorator::ParticleFactory& pf,
                              RMF::decorator::ReferenceFrameFactory& ff) {
    RMF::Vector4 rot(1, 0, 0, 0);
    RMF::Vector3 trans(0, 0, 0);
    for (std::size_t i = 0; i < li.frames.size(); ++i) {
        RMF::decorator::ReferenceFrameConst d = ff.get(li.frames[i]);
        const RMF::Vector3 rt = quat_rot(rot, d.get_translation());
        trans = RMF::Vector3(trans[0] + rt[0], trans[1] + rt[1],
                             trans[2] + rt[2]);
        rot = quat_mul(rot, d.get_rotation());
    }
    const RMF::Vector3 local =
            pf.get_is(li.node) ? pf.get(li.node).get_coordinates()
                               : RMF::Vector3(0, 0, 0);
    const RMF::Vector3 r = quat_rot(rot, local);
    return RMF::Vector3(trans[0] + r[0], trans[1] + r[1], trans[2] + r[2]);
}

//! The leaves of the first representation root, or of the file root.
/*! `IMP.rmf::create_hierarchies` takes the children of the file root that are
    representations; a file written by this module has exactly one. Falling
    back to the root itself keeps files that were written flat readable. */
std::vector<LeafInfo> library_leaves(RMF::FileConstHandle fh) {
    RMF::decorator::ResidueFactory rf(fh);
    RMF::decorator::ChainFactory cf(fh);
    RMF::decorator::ReferenceFrameFactory ff(fh);
    const std::vector<RMF::NodeConstHandle> no_frames;
    std::vector<LeafInfo> leaves;
    const RMF::NodeConstHandles roots = fh.get_root_node().get_children();
    for (std::size_t i = 0; i < roots.size(); ++i) {
        if (roots[i].get_type() == RMF::REPRESENTATION) {
            collect_leaves(roots[i], rf, cf, ff, -1, "", "", no_frames, leaves);
            if (!leaves.empty()) return leaves;
        }
    }
    collect_leaves(fh.get_root_node(), rf, cf, ff, -1, "", "", no_frames,
                   leaves);
    return leaves;
}

}  // namespace

void write_rmf(double* coords, int n_atoms, int n_dim, const std::string& path,
               std::string model_name, std::string metadata_json,
               double radius) {
    if (n_dim != 3) {
        IMP_THROW("write_rmf takes (N, 3) coordinates, not (N, " << n_dim
                                                                 << ")",
                  ValueException);
    }
    RMF::FileHandle fh = RMF::create_rmf_file(path);
    if (!metadata_json.empty()) fh.set_description(metadata_json);

    RMF::decorator::ParticleFactory pf(fh);
    RMF::NodeHandle root = fh.get_root_node().add_child(model_name,
                                                        RMF::REPRESENTATION);
    for (int i = 0; i < n_atoms; ++i) {
        std::ostringstream name;
        name << "p" << i;
        RMF::NodeHandle n = root.add_child(name.str(), RMF::REPRESENTATION);
        RMF::decorator::Particle d = pf.get(n);
        d.set_coordinates(RMF::Vector3(static_cast<float>(coords[3 * i]),
                                       static_cast<float>(coords[3 * i + 1]),
                                       static_cast<float>(coords[3 * i + 2])));
        d.set_radius(static_cast<float>(radius));
        // Mass is set because RMF's particle decorator carries it and readers
        // that expect a particle expect all three; the value is not meaningful
        // here, which is why it is a flat 1.0 as it always was.
        d.set_mass(1.0f);
    }
    fh.add_frame("frame_0", RMF::FRAME);
}

void write_rotamer_library_rmf(const std::string& path,
                               const RotamerLibrary& library) {
    const int n_rotamers = library.n_rotamers;
    const int n_atoms = library.n_atoms;
    if (n_rotamers <= 0 || n_atoms <= 0) {
        IMP_THROW("an empty rotamer library has no frames to write",
                  ValueException);
    }
    if (library.coords.size() !=
        static_cast<std::size_t>(n_rotamers) * n_atoms * 3) {
        IMP_THROW("the library's coordinates are "
                          << library.coords.size()
                          << " numbers, not n_rotamers * n_atoms * 3 = "
                          << n_rotamers * n_atoms * 3,
                  ValueException);
    }
    if (library.atom_names.size() != static_cast<std::size_t>(n_atoms)) {
        IMP_THROW("the library names " << library.atom_names.size()
                                       << " atoms and has " << n_atoms,
                  ValueException);
    }

    RMF::FileHandle fh = RMF::create_rmf_file(with_rmf3(path));
    RMF::decorator::ParticleFactory pf(fh);
    RMF::decorator::ChainFactory cf(fh);
    RMF::decorator::ResidueFactory rf(fh);

    // root -> chain "A" -> residue "DYE" -> one node per atom. The shape is
    // the one IMP.rmf produced from the equivalent IMP hierarchy, kept so that
    // every already-written library still reads.
    RMF::NodeHandle root =
            fh.get_root_node().add_child("rotamers", RMF::REPRESENTATION);
    RMF::NodeHandle chain = root.add_child("A", RMF::REPRESENTATION);
    cf.get(chain).set_chain_id("A");
    RMF::NodeHandle res = chain.add_child("DYE", RMF::REPRESENTATION);
    RMF::decorator::Residue rd = rf.get(res);
    rd.set_residue_index(1);
    rd.set_residue_type("DYE");

    std::vector<RMF::NodeHandle> atoms;
    atoms.reserve(n_atoms);
    for (int a = 0; a < n_atoms; ++a) {
        RMF::NodeHandle n = res.add_child(library.atom_names[a],
                                          RMF::REPRESENTATION);
        RMF::decorator::Particle d = pf.get(n);
        d.set_radius(1.0f);
        d.set_mass(1.0f);
        atoms.push_back(n);
    }

    RMF::FloatKey weight_key = fh.get_key<RMF::FloatTag>(
            fh.get_category(kWeightCategory), kWeightKey);
    RMF::IntsKey trans_key = fh.get_key<RMF::IntsTag>(
            fh.get_category(kKineticCategory), kTransitionKey);
    RMF::NodeHandle root_node = fh.get_root_node();
    if (!library.transitions.empty()) {
        root_node.set_static_value(
                trans_key, RMF::Ints(library.transitions.begin(),
                                     library.transitions.end()));
    }

    for (int r = 0; r < n_rotamers; ++r) {
        std::ostringstream name;
        name << (r + 1);
        fh.add_frame(name.str(), RMF::FRAME);
        for (int a = 0; a < n_atoms; ++a) {
            const std::size_t o = (static_cast<std::size_t>(r) * n_atoms + a) * 3;
            pf.get(atoms[a]).set_coordinates(
                    RMF::Vector3(static_cast<float>(library.coords[o]),
                                 static_cast<float>(library.coords[o + 1]),
                                 static_cast<float>(library.coords[o + 2])));
        }
        root_node.set_frame_value(
                weight_key,
                static_cast<float>(r < static_cast<int>(library.weights.size())
                                           ? library.weights[r]
                                           : 0.0));
    }
}

RotamerLibrary read_rotamer_library_rmf(const std::string& path) {
    std::string file = path;
    if (!internal::file_exists(file)) {
        if (internal::file_exists(with_rmf3(file))) {
            file = with_rmf3(file);
        } else {
            IMP_THROW("RMF library not found: " << path, IOException);
        }
    }

    RMF::FileConstHandle fh = RMF::open_rmf_file_read_only(file);
    RMF::decorator::ParticleFactory pf(fh);
    RMF::decorator::ReferenceFrameFactory ff(fh);
    // Position on a real frame before walking the tree. `get_is` is answered
    // from the *current* frame's values, and on the file RMF opens by default
    // three of this trajectory's twenty reference frames do not answer yes --
    // so a walk done first builds ancestor chains with three transforms
    // missing, and every frame but the one that happens to need none comes out
    // tens of angstroms wrong.
    if (fh.get_number_of_frames() > 0) fh.set_current_frame(RMF::FrameID(0));
    const std::vector<LeafInfo> atoms = library_leaves(fh);
    if (atoms.empty()) IMP_THROW("no hierarchy in " << file, ValueException);

    RMF::FloatKey weight_key = fh.get_key<RMF::FloatTag>(
            fh.get_category(kWeightCategory), kWeightKey);
    RMF::IntsKey trans_key = fh.get_key<RMF::IntsTag>(
            fh.get_category(kKineticCategory), kTransitionKey);
    RMF::NodeConstHandle root_node = fh.get_root_node();

    RotamerLibrary lib;
    lib.path = file;
    lib.n_atoms = static_cast<int>(atoms.size());
    lib.n_rotamers = static_cast<int>(fh.get_number_of_frames());
    for (std::size_t a = 0; a < atoms.size(); ++a) {
        lib.atom_names.push_back(atoms[a].node.get_name());
    }
    lib.coords.reserve(static_cast<std::size_t>(lib.n_rotamers) * lib.n_atoms * 3);
    for (int f = 0; f < lib.n_rotamers; ++f) {
        fh.set_current_frame(RMF::FrameID(f));
        RMF::Nullable<float> w = root_node.get_frame_value(weight_key);
        lib.weights.push_back(w.get_is_null() ? 0.0 : static_cast<double>(w.get()));
        for (std::size_t a = 0; a < atoms.size(); ++a) {
            const RMF::Vector3 c = leaf_coordinates(atoms[a], pf, ff);
            lib.coords.push_back(c[0]);
            lib.coords.push_back(c[1]);
            lib.coords.push_back(c[2]);
        }
    }

    RMF::Nullable<RMF::Ints> t = root_node.get_static_value(trans_key);
    if (!t.get_is_null()) {
        const RMF::Ints& flat = t.get();
        lib.transitions.assign(flat.begin(), flat.end());
    }
    return lib;
}

std::vector<ProteinFrame> protein_frames_from_rmf(const std::string& path,
                                                  int max_frames) {
    if (!internal::file_exists(path)) {
        IMP_THROW("RMF trajectory not found: " << path, IOException);
    }
    RMF::FileConstHandle fh = RMF::open_rmf_file_read_only(path);
    RMF::decorator::ParticleFactory pf(fh);
    RMF::decorator::ReferenceFrameFactory ff(fh);
    // Position on a real frame before walking the tree. `get_is` is answered
    // from the *current* frame's values, and on the file RMF opens by default
    // three of this trajectory's twenty reference frames do not answer yes --
    // so a walk done first builds ancestor chains with three transforms
    // missing, and every frame but the one that happens to need none comes out
    // tens of angstroms wrong.
    if (fh.get_number_of_frames() > 0) fh.set_current_frame(RMF::FrameID(0));
    const std::vector<LeafInfo> atoms = library_leaves(fh);
    if (atoms.empty()) IMP_THROW("no hierarchy in " << path, ValueException);

    // The labels do not change between frames, so they are read once. The
    // hierarchy walk that finds them is the expensive part of this function --
    // doing it per frame made a 1000-frame trajectory pay for it a thousand
    // times over for an answer that never moved.
    const std::size_t n_atoms = atoms.size();
    std::vector<std::string> atom_names(n_atoms), atom_types(n_atoms);
    std::vector<std::string> resnames(n_atoms), chain_ids(n_atoms);
    std::vector<int> residue_indices(n_atoms);
    for (std::size_t a = 0; a < n_atoms; ++a) {
        atom_names[a] = atoms[a].node.get_name();
        atom_types[a] = "";
        residue_indices[a] = atoms[a].residue_index;
        resnames[a] = atoms[a].resname;
        chain_ids[a] = atoms[a].chain_id;
    }

    int n = static_cast<int>(fh.get_number_of_frames());
    if (max_frames >= 0 && max_frames < n) n = max_frames;
    std::vector<ProteinFrame> frames;
    frames.reserve(n);
    for (int f = 0; f < n; ++f) {
        fh.set_current_frame(RMF::FrameID(f));
        ProteinFrame frame;
        frame.atom_names = atom_names;
        frame.atom_types = atom_types;
        frame.resnames = resnames;
        frame.chain_ids = chain_ids;
        frame.residue_indices = residue_indices;
        frame.coords.reserve(n_atoms * 3);
        for (std::size_t a = 0; a < n_atoms; ++a) {
            const RMF::Vector3 c = leaf_coordinates(atoms[a], pf, ff);
            frame.coords.push_back(c[0]);
            frame.coords.push_back(c[1]);
            frame.coords.push_back(c[2]);
        }
        frames.push_back(frame);
    }
    return frames;
}

IMPBFF_END_NAMESPACE
