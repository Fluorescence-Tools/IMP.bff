/**
 * \file RmfIO.cpp
 * \brief Structures, rotamer libraries and trajectories through RMF.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/RmfIO.h>
#include <IMP/bff/internal/Text.h>

#include <IMP/Model.h>
#include <IMP/Particle.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/atom/Mass.h>
#include <IMP/core/XYZR.h>
#include <IMP/rmf/atom_io.h>
#include <IMP/rmf/frames.h>

#include <RMF/FileConstHandle.h>
#include <RMF/FileHandle.h>
#include <RMF/Nullable.h>
#include <RMF/decorator/physics.h>
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

}  // namespace

void write_rmf(double* coords, int n_atoms, int n_dim, const std::string& path,
               std::string model_name, std::string metadata_json,
               double radius) {
    if (n_dim != 3) {
        IMP_THROW("write_rmf takes (N, 3) coordinates, not (N, " << n_dim
                                                                 << ")",
                  ValueException);
    }
    IMP_NEW(IMP::Model, m, ());
    IMP::atom::Hierarchy root = IMP::atom::Hierarchy::setup_particle(
            new IMP::Particle(m, model_name));
    for (int i = 0; i < n_atoms; ++i) {
        std::ostringstream name;
        name << "p" << i;
        IMP::Particle* p = new IMP::Particle(m, name.str());
        IMP::core::XYZR::setup_particle(
                p, IMP::algebra::Sphere3D(
                           IMP::algebra::Vector3D(coords[3 * i],
                                                  coords[3 * i + 1],
                                                  coords[3 * i + 2]),
                           radius));
        IMP::atom::Mass::setup_particle(p, 1.0);
        root.add_child(IMP::atom::Hierarchy::setup_particle(p));
    }

    RMF::FileHandle fh = RMF::create_rmf_file(path);
    if (!metadata_json.empty()) fh.set_description(metadata_json);
    IMP::rmf::add_hierarchy(fh, root);
    IMP::rmf::save_frame(fh, "frame_0");
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

    IMP_NEW(IMP::Model, m, ());
    IMP::atom::Hierarchy root =
            IMP::atom::Hierarchy::setup_particle(new IMP::Particle(m, "rotamers"));
    IMP::atom::Chain chain = IMP::atom::Chain::setup_particle(
            new IMP::Particle(m, "A"), "A");
    root.add_child(chain);
    IMP::atom::Residue res = IMP::atom::Residue::setup_particle(
            new IMP::Particle(m, "DYE"), IMP::atom::ResidueType("DYE"), 1);
    chain.add_child(res);

    IMP::Particles particles;
    for (int a = 0; a < n_atoms; ++a) {
        IMP::Particle* p = new IMP::Particle(m, library.atom_names[a]);
        // Mass rather than Atom, so the particle name survives into the RMF.
        IMP::atom::Mass::setup_particle(p, 1.0);
        IMP::core::XYZR::setup_particle(p).set_radius(1.0);
        res.add_child(IMP::atom::Hierarchy::setup_particle(p));
        particles.push_back(p);
    }

    RMF::FileHandle fh = RMF::create_rmf_file(with_rmf3(path));
    IMP::rmf::add_hierarchies(fh, IMP::atom::Hierarchies(1, root));
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
        for (int a = 0; a < n_atoms; ++a) {
            const std::size_t o = (static_cast<std::size_t>(r) * n_atoms + a) * 3;
            IMP::core::XYZ(particles[a]).set_coordinates(
                    IMP::algebra::Vector3D(library.coords[o],
                                           library.coords[o + 1],
                                           library.coords[o + 2]));
        }
        std::ostringstream name;
        name << (r + 1);
        IMP::rmf::save_frame(fh, name.str());
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

    IMP_NEW(IMP::Model, m, ());
    RMF::FileConstHandle fh = RMF::open_rmf_file_read_only(file);
    IMP::atom::Hierarchies roots = IMP::rmf::create_hierarchies(fh, m);
    if (roots.empty()) IMP_THROW("no hierarchy in " << file, ValueException);
    const IMP::atom::Hierarchies atoms = IMP::atom::get_leaves(roots[0]);

    RMF::FloatKey weight_key = fh.get_key<RMF::FloatTag>(
            fh.get_category(kWeightCategory), kWeightKey);
    RMF::IntsKey trans_key = fh.get_key<RMF::IntsTag>(
            fh.get_category(kKineticCategory), kTransitionKey);
    RMF::NodeConstHandle root_node = fh.get_root_node();

    RotamerLibrary lib;
    lib.path = file;
    lib.n_atoms = static_cast<int>(atoms.size());
    lib.n_rotamers = static_cast<int>(fh.get_number_of_frames());
    for (unsigned int a = 0; a < atoms.size(); ++a) {
        lib.atom_names.push_back(atoms[a]->get_name());
    }
    lib.coords.reserve(static_cast<std::size_t>(lib.n_rotamers) * lib.n_atoms * 3);
    for (int f = 0; f < lib.n_rotamers; ++f) {
        IMP::rmf::load_frame(fh, RMF::FrameID(f));
        RMF::Nullable<float> w = root_node.get_frame_value(weight_key);
        lib.weights.push_back(w.get_is_null() ? 0.0 : static_cast<double>(w.get()));
        for (unsigned int a = 0; a < atoms.size(); ++a) {
            const IMP::algebra::Vector3D c =
                    IMP::core::XYZ(atoms[a]).get_coordinates();
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
    IMP_NEW(IMP::Model, m, ());
    RMF::FileConstHandle fh = RMF::open_rmf_file_read_only(path);
    IMP::atom::Hierarchies hierarchies = IMP::rmf::create_hierarchies(fh, m);
    if (hierarchies.empty()) {
        IMP_THROW("no hierarchy in " << path, ValueException);
    }

    int n = static_cast<int>(fh.get_number_of_frames());
    if (max_frames >= 0 && max_frames < n) n = max_frames;
    std::vector<ProteinFrame> frames;
    frames.reserve(n);
    for (int f = 0; f < n; ++f) {
        IMP::rmf::load_frame(fh, RMF::FrameID(f));
        frames.push_back(protein_frame_from_hierarchy(hierarchies[0]));
    }
    return frames;
}

IMPBFF_END_NAMESPACE
