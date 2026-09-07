/**
 * \file ProbeSampling.cpp
 * \brief Sampling a tethered dye: the walk as an object, then the kernels.
 *
 * Sections in the order of IMP/bff/ProbeSampling.h; each is marked with the
 * file it came from.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

// -------- from ProbeDiffusion.cpp --------
/**
 * (formerly ProbeDiffusion.cpp, now a section of this file)
 * \brief The dye's Brownian walk in its accessible volume, as an object.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/ProbeSampling.h>
#include <IMP/bff/BrownianWalk.h>
#include <IMP/bff/internal/GridShape.h>
#include <IMP/bff/internal/OutputView.h>
#include <IMP/bff/Base.h>

#include <cmath>
#include <cstdlib>
#include <random>
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

ProbeDiffusionSimulation::ProbeDiffusionSimulation(
        const std::vector<int>& density, double dg,
        const std::vector<double>& x0,
        const std::vector<int>& slow_density,
        const std::vector<double>& slow_factor_map,
        const std::vector<double>& quenching_rate_map)
        : density_(density), slow_density_(slow_density),
          slow_factor_map_(slow_factor_map),
          quenching_rate_map_(quenching_rate_map),
          x0_(x0.empty() ? std::vector<double>(3, 0.0) : x0),
          dg_(dg), t_step_(0.0), ng_(internal::cube_side(density.size())),
          n_accepted_(0), n_rejected_(0) {
    if (!density_.empty() && ng_ == 0) {
        IMP_THROW("an occupancy grid must be cubic; " << density_.size()
                          << " values are not a whole cube",
                  IMP::ValueException);
    }
    if (x0_.size() != 3) {
        IMP_THROW("the grid anchor is three coordinates", IMP::ValueException);
    }
}

int ProbeDiffusionSimulation::run(double D, double slow_fact, double t_step,
                                double t_max, int n_trajectories, int seed) {
    t_step_ = t_step;
    if (n_trajectories < 0) n_trajectories = default_trajectory_count();
    if (n_trajectories < 1) n_trajectories = 1;

    // Mobility is a *field*: one scaling per voxel. The scalar-plus-mask form
    // is one way to build it, and building it here is what collapsed two
    // near-identical kernels into one.
    std::vector<double> mobility;
    if (!slow_factor_map_.empty()) {
        mobility = slow_factor_map_;
    } else if (slow_fact != 1.0 && !slow_density_.empty()) {
        mobility.resize(density_.size(), 1.0);
        for (std::size_t i = 0; i < mobility.size() && i < slow_density_.size(); ++i) {
            if (slow_density_[i]) mobility[i] = slow_fact;
        }
    }

    // Distinct seeds per walk, derived from one. Reusing the base seed would
    // run n identical trajectories and concatenate them, which looks like n
    // times the sampling and is none of it.
    //
    // The stride is a large prime, so trajectories from one base do not share
    // the low-order pattern that `base + i` would give them. This is exactly
    // a fixed derivation, digit for digit, so a walk seeded here reproduces a
    // walk seeded anywhere else that follows it.
    std::vector<int> seeds(n_trajectories);
    if (seed < 0) {
        std::random_device rd;
        for (int i = 0; i < n_trajectories; ++i) {
            seeds[i] = static_cast<int>(rd() & 0x7fffffff);
        }
    } else {
        const long long max_seed = (1LL << 31) - 1;
        for (int i = 0; i < n_trajectories; ++i) {
            seeds[i] = static_cast<int>(
                    (static_cast<long long>(seed) + i * 104729LL) % max_seed);
        }
    }

    std::vector<std::vector<double>> parts(n_trajectories);
    std::vector<std::vector<int>> counts(n_trajectories);

    auto one = [&](int i) {
        double* buf = nullptr;
        int n = 0;
        counts[i].assign(2, 0);
        brownian_walk_in_volume(
                density_.empty() ? nullptr : const_cast<int*>(density_.data()),
                static_cast<int>(density_.size()),
                mobility.empty() ? nullptr : mobility.data(),
                static_cast<int>(mobility.size()),
                ng_, dg_, t_max, t_step, D, seeds[i], counts[i], &buf, &n);
        if (buf != nullptr) {
            parts[i].assign(buf, buf + n);
            std::free(buf);
        }
    };

    if (n_trajectories == 1) {
        one(0);
    } else {
        // std::thread, not OpenMP: this build leaves OpenMP_CXX_FLAGS empty, so
        // every `#pragma omp` in the module is inert and a trajectory fan-out
        // written with one would silently run on a single core.
        std::vector<std::thread> pool;
        pool.reserve(n_trajectories);
        for (int i = 0; i < n_trajectories; ++i) pool.emplace_back(one, i);
        for (auto& t : pool) t.join();
    }

    trajectory_.clear();
    n_accepted_ = n_rejected_ = 0;
    for (int i = 0; i < n_trajectories; ++i) {
        if (counts[i].size() < 2 || counts[i][0] <= 0) continue;   // no start found
        n_accepted_ += counts[i][0];
        n_rejected_ += counts[i][1];
        // The walk reports (x, y, z, accepted) per step; only the coordinates
        // are kept, and shifted into the structure's frame here so that every
        // consumer reads absolute positions.
        const std::size_t n_steps = parts[i].size() / 4;
        const std::size_t base = trajectory_.size();
        trajectory_.resize(base + n_steps * 3);
        for (std::size_t s = 0; s < n_steps; ++s) {
            trajectory_[base + 3 * s + 0] = parts[i][4 * s + 0] + x0_[0];
            trajectory_[base + 3 * s + 1] = parts[i][4 * s + 1] + x0_[1];
            trajectory_[base + 3 * s + 2] = parts[i][4 * s + 2] + x0_[2];
        }
    }
    return get_n_frames();
}

void ProbeDiffusionSimulation::get_trajectory(double** out_view,
                                            int* n_out_view) const {
    internal::copy_to_view(trajectory_, out_view, n_out_view);
}

void ProbeDiffusionSimulation::get_x0(double** out_view, int* n_out_view) const {
    internal::copy_to_view(x0_, out_view, n_out_view);
}

void ProbeDiffusionSimulation::get_mean_position(double** out_view,
                                               int* n_out_view) const {
    const int n = get_n_frames();
    if (n == 0) {
        IMP_THROW("Run the simulation first.", IMP::ValueException);
    }
    std::vector<double> m(3, 0.0);
    for (int i = 0; i < n; ++i) {
        m[0] += trajectory_[3 * i + 0];
        m[1] += trajectory_[3 * i + 1];
        m[2] += trajectory_[3 * i + 2];
    }
    for (int k = 0; k < 3; ++k) m[k] /= n;
    internal::copy_to_view(m, out_view, n_out_view);
}

void ProbeDiffusionSimulation::sample_grid(const std::vector<double>& grid, int ng,
                                         double** out_view,
                                         int* n_out_view) const {
    const int n = get_n_frames();
    if (n == 0) {
        IMP_THROW("Run the simulation first.", IMP::ValueException);
    }
    double* out = internal::new_double_view(static_cast<std::size_t>(n),
                                            out_view, n_out_view);
    if (out == nullptr) return;
    const int centre = grid_center_index(ng);
    for (int f = 0; f < n; ++f) {
        bool inside = true;
        std::size_t flat = 0;
        for (int axis = 0; axis < 3; ++axis) {
            // floor, not trunc: trunc maps [-1, 0) to 0, so a position up to one
            // voxel below the grid would be read as voxel 0.
            const double u = (trajectory_[3 * f + axis] - x0_[axis]) / dg_ + centre;
            const long idx = static_cast<long>(std::floor(u));
            if (idx < 0 || idx >= ng) { inside = false; break; }
            flat = flat * static_cast<std::size_t>(ng) + static_cast<std::size_t>(idx);
        }
        out[f] = (inside && flat < grid.size()) ? grid[flat] : 0.0;
    }
}

void ProbeDiffusionSimulation::get_density(double** out_view,
                                         int* n_out_view) const {
    std::vector<double> d(density_.begin(), density_.end());
    internal::copy_to_view(d, out_view, n_out_view);
}

void ProbeDiffusionSimulation::set_density(const std::vector<int>& density) {
    const int ng = internal::cube_side(density.size());
    if (!density.empty() && ng == 0) {
        IMP_THROW("an occupancy grid must be cubic; " << density.size()
                          << " values are not a whole cube",
                  IMP::ValueException);
    }
    density_ = density;
    ng_ = ng;
    trajectory_.clear();
    n_accepted_ = n_rejected_ = 0;
}

void ProbeDiffusionSimulation::get_slow_factor_map(double** out_view,
                                                 int* n_out_view) const {
    internal::copy_to_view(slow_factor_map_, out_view, n_out_view);
}

void ProbeDiffusionSimulation::get_slow_density(double** out_view,
                                              int* n_out_view) const {
    std::vector<double> d(slow_density_.begin(), slow_density_.end());
    internal::copy_to_view(d, out_view, n_out_view);
}

void ProbeDiffusionSimulation::get_quenching_rate_map(double** out_view,
                                                    int* n_out_view) const {
    internal::copy_to_view(quenching_rate_map_, out_view, n_out_view);
}

void ProbeDiffusionSimulation::get_k_quench(double** out_view,
                                          int* n_out_view) const {
    const int n = get_n_frames();
    if (n == 0) {
        IMP_THROW("Run the simulation first.", IMP::ValueException);
    }
    if (quenching_rate_map_.empty()) {
        internal::new_double_view(static_cast<std::size_t>(n), out_view, n_out_view);
        return;                                   // calloc already zeroed it
    }
    sample_grid(quenching_rate_map_, internal::cube_side(quenching_rate_map_.size()),
                out_view, n_out_view);
    // Rounded through `float`, and deliberately. The fused kernel holds its own
    // trace as `std::vector<float>` -- it halves the memory the photon race
    // walks, and that race makes tens of millions of random reads into it -- so
    // the two paths see bit-identical rates rather than rates that agree to
    // 2.6e-8 and could in principle disagree about a photon. The precision costs
    // nothing real: a PET rate constant is a transferable starting value known
    // to perhaps two significant figures, and float carries seven.
    //
    // Here rather than in a caller, so every consumer sees the same rates and
    // the equality holds for all of them.
    if (*out_view != nullptr) {
        for (int i = 0; i < *n_out_view; ++i) {
            (*out_view)[i] = static_cast<float>((*out_view)[i]);
        }
    }
}

double ProbeDiffusionSimulation::get_collision_fraction() const {
    const int n = get_n_frames();
    if (n == 0) return 0.0;
    double* k = nullptr;
    int nk = 0;
    get_k_quench(&k, &nk);
    if (k == nullptr) return 0.0;
    int hit = 0;
    for (int i = 0; i < nk; ++i) if (k[i] > 0.0) ++hit;
    std::free(k);
    return static_cast<double>(hit) / n;
}

IMPBFF_END_NAMESPACE

// -------- from ProbeSampling.cpp --------


#include <IMP/bff/RotamerLibrary.h>
#include <IMP/bff/TrajectoryIO.h>
#include <IMP/bff/internal/Text.h>

#include <IMP/atom/Atom.h>
#include <IMP/atom/Hierarchy.h>
#include <IMP/atom/pdb.h>
#include <IMP/algebra/Vector3D.h>
#include <IMP/core/XYZ.h>

#include <cctype>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

IMPBFF_BEGIN_NAMESPACE

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

std::vector<std::string> find_reference_rotamer_files(
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

void apply_coordinates(const IMP::atom::Hierarchy hierarchy,
                               const std::vector<double>& coords) {
    IMP::atom::Hierarchies leaves = IMP::atom::get_leaves(hierarchy);
    const std::size_t n_atoms = leaves.size();
    if (coords.size() != n_atoms * 3) {
        IMP_THROW("Atom count mismatch: coords=" << coords.size() / 3
                          << " hierarchy=" << n_atoms,
                  ValueException);
    }
    for (std::size_t i = 0; i < n_atoms; ++i) {
        IMP::core::XYZ atom(leaves[i]);
        atom.set_coordinates(IMP::algebra::Vector3D(
                coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]));
    }
}

IMPBFF_END_NAMESPACE
