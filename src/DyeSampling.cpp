/**
 * \file DyeSampling.cpp
 * \brief Sampling a tethered dye: the walk, the photons, and the frame writers.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/DyeSampling.h>

#include <IMP/bff/BrownianWalk.h>
#include <IMP/bff/TrajectoryIO.h>
#include <IMP/bff/internal/OutputView.h>

#include <IMP/atom/Atom.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/atom/pdb.h>
#include <IMP/exception.h>

#include <cstdlib>
#include <fstream>
#include <random>
#include <sstream>
#include <sys/stat.h>

IMPBFF_BEGIN_NAMESPACE

// --------------------------------------------------------------------------
// the walk
// --------------------------------------------------------------------------

void DyeDiffusionTrajectory::get_xyz(double** out_view, int* n_out_view) const {
    internal::copy_to_view(xyz_, out_view, n_out_view);
}

void DyeDiffusionTrajectory::get_accepted(int** out_view_i,
                                          int* n_out_view_i) const {
    int* out = internal::new_int_view(accepted_.size(), out_view_i, n_out_view_i);
    if (out == NULL) return;
    for (std::size_t i = 0; i < accepted_.size(); ++i) out[i] = accepted_[i];
}

double DyeDiffusionTrajectory::get_acceptance_ratio() const {
    const int total = n_accepted_ + n_rejected_;
    return total > 0 ? static_cast<double>(n_accepted_) / total : 0.0;
}

DyeDiffusionTrajectory simulate_dye_diffusion(int* density, int n_density,
                                              double* mobility, int n_mobility,
                                              int ng, double dg, double t_max,
                                              double t_step,
                                              double diffusion_coefficient,
                                              int seed) {
    std::vector<int> counts;
    double* packed = NULL;
    int n_packed = 0;
    brownian_walk_in_volume(density, n_density, mobility, n_mobility, ng, dg,
                            t_max, t_step, diffusion_coefficient, seed, counts,
                            &packed, &n_packed);

    std::vector<double> xyz;
    std::vector<int> accepted;
    const int n_frames = n_packed / 4;
    xyz.reserve(static_cast<std::size_t>(n_frames) * 3);
    accepted.reserve(n_frames);
    for (int i = 0; i < n_frames; ++i) {
        xyz.push_back(packed[i * 4 + 0]);
        xyz.push_back(packed[i * 4 + 1]);
        xyz.push_back(packed[i * 4 + 2]);
        accepted.push_back(packed[i * 4 + 3] != 0.0 ? 1 : 0);
    }
    std::free(packed);

    const int n_accepted = counts.size() == 2 ? counts[0] : 0;
    const int n_rejected = counts.size() == 2 ? counts[1] : 0;

    if (n_frames == 0) {
        // No accessible starting voxel. An empty trajectory of the right length
        // rather than of zero length: every consumer indexes it against the
        // time axis it asked for.
        const int n_steps = static_cast<int>(t_max / t_step);
        xyz.assign(static_cast<std::size_t>(n_steps) * 3, 0.0);
        accepted.assign(n_steps, 0);
    }
    return DyeDiffusionTrajectory(xyz, accepted, n_accepted, n_rejected);
}

void equilibrium_occupancy(double* diffusion_map, int n_diffusion_map,
                           int* bounds, int n_bounds,
                           const std::string& flux_form, double** out_view,
                           int* n_out_view) {
    const bool smoluchowski = flux_form == "smoluchowski";
    if (!smoluchowski && flux_form != "ito") {
        IMP_THROW("flux_form must be 'smoluchowski' or 'ito', not '"
                          << flux_form << "'",
                  ValueException);
    }
    if (n_bounds != n_diffusion_map) {
        IMP_THROW("the bounds and the diffusion map must have the same shape: "
                          << n_bounds << " against " << n_diffusion_map,
                  ValueException);
    }

    double* out = internal::new_double_view(n_diffusion_map, out_view, n_out_view);
    if (out == NULL) return;

    double total = 0.0;
    for (int i = 0; i < n_diffusion_map; ++i) {
        out[i] = 0.0;
        if (bounds[i] <= 0) continue;
        if (smoluchowski) {
            out[i] = 1.0;
        } else if (diffusion_map[i] > 0.0) {
            out[i] = 1.0 / diffusion_map[i];
        }
        total += out[i];
    }
    if (total > 0.0) {
        for (int i = 0; i < n_diffusion_map; ++i) out[i] /= total;
    }
}

// --------------------------------------------------------------------------
// the rotamer library
// --------------------------------------------------------------------------

void RotamerLibrary::get_coords(double** out_view, int* n_out_view) const {
    internal::copy_to_view(coords, out_view, n_out_view);
}

void RotamerLibrary::get_weights(double** out_view, int* n_out_view) const {
    internal::copy_to_view(weights, out_view, n_out_view);
}

namespace {

bool file_exists(const std::string& path) {
    struct stat info;
    return stat(path.c_str(), &info) == 0;
}

}  // namespace

RotamerLibrary load_rotamer_library(const std::string& pdb_path,
                                    const std::string& trajectory_path,
                                    const std::string& weights_path,
                                    int max_frames) {
    RotamerLibrary out;

    IMP_NEW(IMP::Model, model, ());
    IMP::atom::Hierarchy hierarchy = IMP::atom::read_pdb(
            pdb_path, model, new IMP::atom::AllPDBSelector());
    const IMP::atom::Hierarchies leaves = IMP::atom::get_leaves(hierarchy);
    for (unsigned int i = 0; i < leaves.size(); ++i) {
        std::string name =
                IMP::atom::Atom(leaves[i]).get_atom_type().get_string();
        const std::size_t a = name.find_first_not_of(' ');
        const std::size_t b = name.find_last_not_of(' ');
        out.atom_names.push_back(a == std::string::npos
                                         ? std::string()
                                         : name.substr(a, b - a + 1));
    }
    out.n_atoms = static_cast<int>(out.atom_names.size());

    // BinaryCIF is what this package stores: the rotamer libraries were
    // re-encoded on 2026-08-19 and 44.78 MB of DCD and XTC became 17.99 MB.
    double* frames = NULL;
    int n_values = 0;
    read_bcif_trajectory(trajectory_path, out.n_atoms, "_rotamer_coord",
                         &frames, &n_values);
    int n_frames_read = out.n_atoms > 0 ? n_values / (out.n_atoms * 3) : 0;
    if (max_frames > 0 && n_frames_read > max_frames) n_frames_read = max_frames;
    if (n_frames_read == 0) {
        std::free(frames);
        IMP_THROW("No frames found in " << trajectory_path, IOException);
    }
    out.n_frames = n_frames_read;
    out.coords.assign(frames, frames + static_cast<std::size_t>(n_frames_read) *
                                               out.n_atoms * 3);
    std::free(frames);

    if (!weights_path.empty() && file_exists(weights_path)) {
        std::ifstream in(weights_path.c_str());
        std::string line;
        while (std::getline(in, line)) {
            std::istringstream item(line);
            double w = 0.0;
            if (item >> w) out.weights.push_back(w);
        }
        if (max_frames > 0 &&
            out.weights.size() > static_cast<std::size_t>(max_frames)) {
            out.weights.resize(max_frames);
        }
        if (out.weights.size() != static_cast<std::size_t>(out.n_frames)) {
            // Truncate to the shorter: a library whose weights and frames
            // disagree has extra frames with no weight, and dropping them is
            // the only reading that does not invent one.
            const int n = std::min<int>(static_cast<int>(out.weights.size()),
                                        out.n_frames);
            out.weights.resize(n);
            out.coords.resize(static_cast<std::size_t>(n) * out.n_atoms * 3);
            out.n_frames = n;
        }
    } else {
        out.weights.assign(out.n_frames, 1.0);
    }

    double total = 0.0;
    for (std::size_t i = 0; i < out.weights.size(); ++i) total += out.weights[i];
    if (!(total > 0.0)) {
        out.weights.assign(out.n_frames, 1.0);
        total = out.n_frames;
    }
    for (std::size_t i = 0; i < out.weights.size(); ++i) out.weights[i] /= total;
    return out;
}

std::vector<std::string> find_reference_rotamer_files(
        const std::string& lib_dir, const std::string& dye_name, int cutoff) {
    std::ostringstream pdb, traj, weights;
    pdb << lib_dir << "/" << dye_name << ".pdb";
    traj << lib_dir << "/" << dye_name << "_cutoff" << cutoff << ".bcif";
    weights << lib_dir << "/" << dye_name << "_cutoff" << cutoff
            << "_weights.txt";

    if (!file_exists(pdb.str()) || !file_exists(traj.str())) {
        IMP_THROW("Missing required reference files for " << dye_name,
                  IOException);
    }
    std::vector<std::string> out;
    out.push_back(pdb.str());
    out.push_back(traj.str());
    out.push_back(file_exists(weights.str()) ? weights.str() : std::string());
    return out;
}

int sample_rotamer_index(const std::vector<double>& weights, int seed) {
    if (weights.empty()) return 0;
    std::mt19937 engine(seed < 0 ? std::random_device()() 
                                 : static_cast<unsigned int>(seed));
    std::discrete_distribution<int> draw(weights.begin(), weights.end());
    return draw(engine);
}

IMPBFF_END_NAMESPACE
