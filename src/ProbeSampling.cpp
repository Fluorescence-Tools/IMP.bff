/**
 * \file ProbeSampling.cpp
 * \brief Probe sampling kernels and reference-library helpers.
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/ProbeSampling.h>
#include <IMP/bff/BrownianWalk.h>
#include <IMP/bff/TrajectoryIO.h>
#include <IMP/bff/internal/PdbFrames.h>
#include <IMP/bff/internal/GridShape.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/bff/internal/Text.h>
#include <IMP/algebra/Vector3D.h>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <random>
#include <sstream>
#include <sys/stat.h>
#include <thread>

IMPBFF_BEGIN_NAMESPACE

const int MAX_PARALLEL_TRAJECTORIES = 8;

int default_trajectory_count() {
    const unsigned n = std::thread::hardware_concurrency();
    const unsigned cap = static_cast<unsigned>(MAX_PARALLEL_TRAJECTORIES);
    const unsigned capped = n == 0 ? 1u : (n < cap ? n : cap);
    return static_cast<int>(capped);
}

int resolve_trajectory_count(int requested) {
    if (requested < 0) return default_trajectory_count();
    if (requested < 1) return 1;
    return requested < MAX_PARALLEL_TRAJECTORIES ? requested
                                                 : MAX_PARALLEL_TRAJECTORIES;
}

std::vector<int> trajectory_seeds(int random_seed, int n_trajectories) {
    const int n = n_trajectories > 0 ? n_trajectories : 1;
    std::vector<int> seeds;
    seeds.reserve(n);
    if (random_seed < 0) {
        // One walk gets -1 -- "draw freely" -- rather than a seed of its own,
        // so an unseeded single walk stays as unconstrained as it was.
        if (n <= 1) { seeds.push_back(-1); return seeds; }
        std::random_device source;
        for (int i = 0; i < n; ++i) {
            seeds.push_back(static_cast<int>(source() & 0x7fffffffu));
        }
        return seeds;
    }
    // A large prime stride, so trajectories from one base seed do not share the
    // low-order pattern `base + i` would give them.
    const long long max_seed = (1LL << 31) - 1;
    for (int i = 0; i < n; ++i) {
        seeds.push_back(static_cast<int>(
                (static_cast<long long>(random_seed) + i * 104729LL) % max_seed));
    }
    return seeds;
}

// --------------------------------------------------------------------------
// the walk
// --------------------------------------------------------------------------

void ProbeDiffusionTrajectory::get_xyz(double** out_view, int* n_out_view) const {
    internal::copy_to_view(xyz_, out_view, n_out_view);
}

void ProbeDiffusionTrajectory::get_accepted(int** out_view_i,
                                          int* n_out_view_i) const {
    int* out = internal::new_int_view(accepted_.size(), out_view_i, n_out_view_i);
    if (out == NULL) return;
    for (std::size_t i = 0; i < accepted_.size(); ++i) out[i] = accepted_[i];
}

double ProbeDiffusionTrajectory::get_acceptance_ratio() const {
    const int total = n_accepted_ + n_rejected_;
    return total > 0 ? static_cast<double>(n_accepted_) / total : 0.0;
}

namespace {
//! `ng` from a flat `ng^3` density; 0 when not a perfect cube.
int grid_side(const std::vector<int>& density) {
    const std::size_t n = density.size();
    const long side =
            static_cast<long>(std::lround(std::cbrt(static_cast<double>(n))));
    if (side <= 0 || static_cast<std::size_t>(side) * side * side != n) return 0;
    return static_cast<int>(side);
}
}  // namespace

ProbeDiffusionTrajectory simulate_probe_diffusion(
        const std::vector<int>& density, const std::vector<int>& slow_density,
        double dg, double t_max, double t_step, double D,
        const std::vector<double>& slow_fact, int random_seed) {
    const int ng = grid_side(density);
    if (ng == 0) {
        IMP_THROW("the density must be a flat cube of accessible voxels, "
                          "nonempty with an integer cube-root side",
                  ValueException);
    }
    const std::size_t n_vox = static_cast<std::size_t>(ng) * ng * ng;
    if (!slow_density.empty() && slow_density.size() != n_vox) {
        IMP_THROW("slow_density has "
                          << slow_density.size() << " entries against density "
                          << density.size(),
                  ValueException);
    }
    if (slow_fact.size() > 1 && slow_fact.size() != n_vox) {
        IMP_THROW("slow_fact has " << slow_fact.size()
                                   << " entries and must be the flat ng^3 field "
                                      "or one scalar",
                  ValueException);
    }

    std::vector<double> mobility;
    if (slow_fact.size() > 1) {
        mobility = slow_fact;
    } else {
        const double scalar = slow_fact.empty() ? 1.0 : slow_fact[0];
        if (scalar != 1.0 && !slow_density.empty()) {
            mobility.resize(n_vox);
            for (std::size_t i = 0; i < n_vox; ++i)
                mobility[i] = slow_density[i] != 0 ? scalar : 1.0;
        }
    }

    std::vector<int> counts;
    double* packed = NULL;
    int n_packed = 0;
    brownian_walk_in_volume(
            density.empty() ? NULL : const_cast<int*>(&density[0]),
            static_cast<int>(density.size()), mobility.data(),
            static_cast<int>(mobility.size()), ng, dg, t_max, t_step, D,
            random_seed, counts, &packed, &n_packed);

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
    return ProbeDiffusionTrajectory(xyz, accepted, n_accepted, n_rejected);
}

void equilibrium_occupancy(const std::vector<double>& diffusion_map,
                           const std::vector<int>& bounds,
                           const std::string& flux_form, double** out_view,
                           int* n_out_view) {
    const bool smoluchowski = flux_form == "smoluchowski";
    if (!smoluchowski && flux_form != "ito") {
        IMP_THROW("flux_form must be 'smoluchowski' or 'ito', not '"
                          << flux_form << "'",
                  ValueException);
    }
    if (diffusion_map.size() != bounds.size()) {
        IMP_THROW("the bounds and the diffusion map must have the same shape: "
                          << bounds.size() << " against " << diffusion_map.size(),
                  ValueException);
    }

    const int n = static_cast<int>(diffusion_map.size());
    double* out = internal::new_double_view(n, out_view, n_out_view);
    if (out == NULL) return;

    double total = 0.0;
    for (int i = 0; i < n; ++i) {
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
        for (int i = 0; i < n; ++i) out[i] /= total;
    }
}

// --------------------------------------------------------------------------
// the rotamer library
// --------------------------------------------------------------------------

namespace {

using IMP::bff::internal::ends_with;
using IMP::bff::internal::file_exists;

}  // namespace

RotamerLibrary load_rotamer_library_trajectory(const std::string& pdb_path,
                                    const std::string& trajectory_path,
                                    const std::string& weights_path,
                                    int max_frames) {
    RotamerLibrary out;

    // A `.drot` is the whole library in one file (PRD-118): it carries the
    // atom names and the weights, so the PDB and the weights file beside it
    // -- when there are any -- have nothing left to say.
    // A `.drot` path, or a locator into a family container -- `read_drot`
    // takes either, so this only has to recognise one.
    const std::string container = split_drot_locator(trajectory_path)[0];
    if (ends_with(container, ".drot") || ends_with(container, ".drot.pto")) {
        const RotamerLibrary drot = read_drot(trajectory_path);
        out.atom_names = drot.atom_names;
        out.n_atoms = drot.n_atoms;
        out.n_rotamers = max_frames > 0 && drot.n_rotamers > max_frames
                ? max_frames : drot.n_rotamers;
        out.coords.assign(drot.coords.begin(),
                          drot.coords.begin() + static_cast<std::size_t>(
                                  out.n_rotamers) * out.n_atoms * 3);
        out.weights.assign(drot.weights.begin(),
                           drot.weights.begin() + out.n_rotamers);
        double total = 0.0;
        for (std::size_t i = 0; i < out.weights.size(); ++i) {
            total += out.weights[i];
        }
        if (!(total > 0.0)) {
            out.weights.assign(out.n_rotamers, 1.0);
            total = out.n_rotamers;
        }
        for (std::size_t i = 0; i < out.weights.size(); ++i) {
            out.weights[i] /= total;
        }
        return out;
    }

    // every ATOM/HETATM of the first model (IMP's AllPDBSelector), the atom
    // type string trimmed -- the core reader, no Model
    const std::vector<ProteinFrame> pdb_frames =
            internal::read_pdb_frames(pdb_path, internal::PDB_ALL, true, 1);
    for (std::size_t i = 0; i < pdb_frames[0].atom_types.size(); ++i) {
        const std::string& name = pdb_frames[0].atom_types[i];
        const std::size_t a = name.find_first_not_of(' ');
        const std::size_t b = name.find_last_not_of(' ');
        out.atom_names.push_back(a == std::string::npos
                                         ? std::string()
                                         : name.substr(a, b - a + 1));
    }
    out.n_atoms = static_cast<int>(out.atom_names.size());

    // Whatever the frame store is: BinaryCIF for the two libraries still
    // shipped that way, DCD for a library a user brought.
    double* frames = NULL;
    int n_values = 0;
    read_trajectory(trajectory_path, out.n_atoms, -1, &frames, &n_values);
    int n_frames_read = out.n_atoms > 0 ? n_values / (out.n_atoms * 3) : 0;
    if (max_frames > 0 && n_frames_read > max_frames) n_frames_read = max_frames;
    if (n_frames_read == 0) {
        std::free(frames);
        IMP_THROW("No frames found in " << trajectory_path, IOException);
    }
    out.n_rotamers = n_frames_read;
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
        if (out.weights.size() != static_cast<std::size_t>(out.n_rotamers)) {
            // Truncate to the shorter: a library whose weights and frames
            // disagree has extra frames with no weight, and dropping them is
            // the only reading that does not invent one.
            const int n = std::min<int>(static_cast<int>(out.weights.size()),
                                        out.n_rotamers);
            out.weights.resize(n);
            out.coords.resize(static_cast<std::size_t>(n) * out.n_atoms * 3);
            out.n_rotamers = n;
        }
    } else {
        out.weights.assign(out.n_rotamers, 1.0);
    }

    double total = 0.0;
    for (std::size_t i = 0; i < out.weights.size(); ++i) total += out.weights[i];
    if (!(total > 0.0)) {
        out.weights.assign(out.n_rotamers, 1.0);
        total = out.n_rotamers;
    }
    for (std::size_t i = 0; i < out.weights.size(); ++i) out.weights[i] /= total;
    return out;
}

std::vector<std::string> get_reference_rotamer_files(
        const std::string& lib_dir, const std::string& probe_name, int cutoff) {
    std::ostringstream pdb, stem, weights;
    pdb << lib_dir << "/" << probe_name << ".pdb";
    stem << lib_dir << "/" << probe_name << "_cutoff" << cutoff;
    weights << stem.str() << "_weights.txt";

    // `.drot` is the shipped store, in the PTO container since v10 and
    // filed in a family container since the three families landed there;
    // `.bcif` remains for the libraries that have no `.drot` and for a set a
    // user points this at.
    const std::string library = probe_name + "_cutoff" + std::to_string(cutoff);
    std::string traj = stem.str() + ".drot.pto";
    if (!file_exists(traj)) traj = stem.str() + ".drot";
    if (!file_exists(traj)) {
        const std::string family = lib_dir + "/dyes.drot.pto";
        if (file_exists(family)) {
            const std::vector<std::string> listed = drot_catalog(family);
            for (std::size_t i = 0; i < listed.size(); ++i) {
                if (listed[i] == library) {
                    traj = family + "::" + library;
                    break;
                }
            }
        }
    }
    if (!file_exists(split_drot_locator(traj)[0])) {
        traj = stem.str() + ".bcif";
    }

    if (!file_exists(pdb.str()) || !file_exists(split_drot_locator(traj)[0])) {
        IMP_THROW("Missing required reference files for " << probe_name,
                  IOException);
    }
    std::vector<std::string> out;
    out.push_back(pdb.str());
    out.push_back(traj);
    out.push_back(file_exists(weights.str()) ? weights.str() : std::string());
    return out;
}

int sample_weighted_index(const std::vector<double>& weights, int seed) {
    if (weights.empty()) return 0;
    std::mt19937 engine(seed < 0 ? std::random_device()() 
                                 : static_cast<unsigned int>(seed));
    std::discrete_distribution<int> draw(weights.begin(), weights.end());
    return draw(engine);
}


IMPBFF_END_NAMESPACE
