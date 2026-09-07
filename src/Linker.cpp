/**
 * \file Linker.cpp
 * \brief A dye's linker: the geometry and the Metropolis sampler.
 *
 * Sections in the order of IMP/bff/Linker.h; each is marked with the file it
 * came from.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

// -------- from LinkerGeometry.cpp --------
/**
 * (formerly LinkerGeometry.cpp, now a section of this file)
 * \brief Applying a torsion/angle configuration to a linker.
 *
 * Copyright 2007-2022 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/Linker.h>

#include <cmath>
#include <stdexcept>

IMPBFF_BEGIN_NAMESPACE

namespace {

struct V3 {
    double x, y, z;
    V3() : x(0), y(0), z(0) {}
    V3(double a, double b, double c) : x(a), y(b), z(c) {}
    V3 operator-(const V3& o) const { return V3(x - o.x, y - o.y, z - o.z); }
    double norm() const { return std::sqrt(x * x + y * y + z * z); }
};

V3 cross(const V3& a, const V3& b) {
    return V3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

V3 at(const std::vector<double>& c, int i) {
    return V3(c[3 * i], c[3 * i + 1], c[3 * i + 2]);
}

//! Rodrigues rotation of `p` about the axis through `origin` along `axis`.
V3 rotate_about(const V3& p, const V3& origin, const V3& axis, double angle) {
    const double len = axis.norm();
    const V3 u(axis.x / len, axis.y / len, axis.z / len);
    const V3 d = p - origin;
    const double c = std::cos(angle), s = std::sin(angle);
    const V3 k = cross(u, d);
    const double dot = u.x * d.x + u.y * d.y + u.z * d.z;
    return V3(origin.x + d.x * c + k.x * s + u.x * dot * (1.0 - c),
              origin.y + d.y * c + k.y * s + u.y * dot * (1.0 - c),
              origin.z + d.z * c + k.z * s + u.z * dot * (1.0 - c));
}

}  // namespace

void LinkerGeometry::add_torsion(int fixed, int moving, const std::vector<int>& moves) {
    torsion_fixed_.push_back(fixed);
    torsion_moving_.push_back(moving);
    torsion_sets_.push_back(moves);
}

void LinkerGeometry::add_angle(int b, int c, int a, const std::vector<int>& moves) {
    angle_b_.push_back(b);
    angle_c_.push_back(c);
    angle_a_.push_back(a);
    angle_sets_.push_back(moves);
}

std::vector<double> LinkerGeometry::apply(const std::vector<double>& config) const {
    const size_t n_torsions = torsion_fixed_.size();
    const size_t n_angles = angle_b_.size();
    if (config.size() != n_torsions + n_angles) {
        throw std::invalid_argument(
            "config must hold one value per torsion followed by one per angle");
    }
    std::vector<double> c(base_);

    for (size_t t = 0; t < n_torsions; ++t) {
        const V3 cf = at(c, torsion_fixed_[t]);
        const V3 cm = at(c, torsion_moving_[t]);
        const V3 axis = cm - cf;
        // two atoms at the same point define no rotation
        if (axis.norm() < 1e-8) continue;
        const double angle = config[t];
        const std::vector<int>& moves = torsion_sets_[t];
        for (size_t i = 0; i < moves.size(); ++i) {
            const V3 p = rotate_about(at(c, moves[i]), cf, axis, angle);
            c[3 * moves[i]] = p.x; c[3 * moves[i] + 1] = p.y; c[3 * moves[i] + 2] = p.z;
        }
    }

    for (size_t a = 0; a < n_angles; ++a) {
        const V3 cb = at(c, angle_b_[a]);
        const V3 cc = at(c, angle_c_[a]);
        const V3 ca = at(c, angle_a_[a]);
        const V3 axis = cross(ca - cb, cc - cb);      // normal of the a-b-c plane
        if (axis.norm() < 1e-8) continue;
        const double angle = config[n_torsions + a];
        const std::vector<int>& moves = angle_sets_[a];
        for (size_t i = 0; i < moves.size(); ++i) {
            const V3 p = rotate_about(at(c, moves[i]), cb, axis, angle);
            c[3 * moves[i]] = p.x; c[3 * moves[i] + 1] = p.y; c[3 * moves[i] + 2] = p.z;
        }
    }
    return c;
}

IMPBFF_END_NAMESPACE

// -------- from LinkerSampling.cpp --------
/**
 * (formerly LinkerSampling.cpp, now a section of this file)
 * \brief Metropolis sampling of a dye's linker, and the library it makes.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/Clustering.h>
#include <IMP/bff/MolecularGraph.h>
#include <IMP/bff/Scoring.h>
#include <IMP/bff/TopologyBuild.h>
#include <IMP/bff/internal/OutputView.h>

#include <IMP/bff/Base.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>

#include <algorithm>
#include <cstdlib>
#include <map>
#include <set>

IMPBFF_BEGIN_NAMESPACE

void LinkerSamplingResult::get_coordinates(double** out_view,
                                           int* n_out_view) const {
    internal::copy_to_view(coordinates, out_view, n_out_view);
}

void LinkerSamplingResult::get_energies(double** out_view,
                                        int* n_out_view) const {
    internal::copy_to_view(energies, out_view, n_out_view);
}

double LinkerSamplingResult::get_acceptance() const {
    return n_steps > 0 ? static_cast<double>(n_accepted) /
                                 static_cast<double>(n_steps)
                       : 0.0;
}

namespace {

//! A hydrogen by name, which is what a MOL2 gives us to go on.
bool is_hydrogen(const std::string& atom_name) {
    return !atom_name.empty() && (atom_name[0] == 'H' || atom_name[0] == 'h');
}

}  // namespace

LinkerGeometry linker_geometry_from_mol2(const std::string& mol2_path,
                                         const std::string& anchor_atom) {
    const Mol2Component component = read_mol2_component(mol2_path, "dye");

    std::vector<int> serials;
    std::map<int, std::string> name_of, element_of;
    std::map<int, IMP::algebra::Vector3D> position_of;
    for (std::size_t i = 0; i < component.atoms.size(); ++i) {
        const Mol2Atom& a = component.atoms[i];
        serials.push_back(a.serial);
        name_of[a.serial] = a.atom_name;
        element_of[a.serial] = a.element;
        position_of[a.serial] = IMP::algebra::Vector3D(a.x, a.y, a.z);
    }
    std::sort(serials.begin(), serials.end());

    const std::vector<std::pair<int, int> >& bonds = component.bonds;
    const MolecularGraph graph(bonds);

    std::vector<std::string> elements;
    for (std::size_t i = 0; i < serials.size(); ++i) {
        elements.push_back(element_of[serials[i]]);
    }
    std::vector<int> ring_sizes;
    ring_sizes.push_back(5);
    ring_sizes.push_back(6);
    ring_sizes.push_back(7);
    const std::vector<int> ring_list =
            graph.get_ring_atoms(serials, elements, ring_sizes);
    const std::set<int> rings(ring_list.begin(), ring_list.end());

    // The anchor: the atom the whole library is expressed relative to.
    int anchor = serials.empty() ? 0 : serials[0];
    for (std::size_t i = 0; i < serials.size(); ++i) {
        if (name_of[serials[i]] == anchor_atom) {
            anchor = serials[i];
            break;
        }
    }

    std::map<int, int> row;
    LinkerGeometry geometry;
    std::vector<double> base;
    for (std::size_t i = 0; i < serials.size(); ++i) {
        row[serials[i]] = static_cast<int>(i);
        const IMP::algebra::Vector3D& v = position_of[serials[i]];
        base.push_back(v[0]);
        base.push_back(v[1]);
        base.push_back(v[2]);
    }
    geometry.set_coordinates(base);

    // -- rotatable bonds --------------------------------------------------
    std::set<std::string> anchor_names;
    anchor_names.insert("N");
    anchor_names.insert("CA");
    anchor_names.insert("C");
    std::vector<std::pair<int, int> > sorted_bonds = bonds;
    std::sort(sorted_bonds.begin(), sorted_bonds.end());
    for (std::size_t i = 0; i < sorted_bonds.size(); ++i) {
        const int a = sorted_bonds[i].first, b = sorted_bonds[i].second;
        const std::string na = name_of[a], nb = name_of[b];
        if (anchor_names.count(na) > 0 || anchor_names.count(nb) > 0) continue;
        if (rings.count(a) > 0 && rings.count(b) > 0) continue;
        if (is_hydrogen(na) || is_hydrogen(nb)) continue;
        const GraphRotor rotor = graph.get_bond_rotor(a, b, anchor);
        if (!rotor.is_rotatable || rotor.moving_nodes.empty()) continue;
        std::vector<int> moves;
        for (std::size_t k = 0; k < rotor.moving_nodes.size(); ++k) {
            std::map<int, int>::const_iterator r =
                    row.find(rotor.moving_nodes[k]);
            if (r != row.end()) moves.push_back(r->second);
        }
        if (moves.empty()) continue;
        geometry.add_torsion(row[rotor.fixed], row[rotor.moving], moves);
    }

    // -- rotatable angles -------------------------------------------------
    for (std::size_t i = 0; i < serials.size(); ++i) {
        const int b = serials[i];
        if (rings.count(b) > 0) continue;
        const std::vector<int> neighbours = graph.get_neighbors(b);
        for (std::size_t j = 0; j < neighbours.size(); ++j) {
            for (std::size_t k = j + 1; k < neighbours.size(); ++k) {
                const int a = neighbours[j], c = neighbours[k];
                if (is_hydrogen(name_of[a]) || is_hydrogen(name_of[c])) {
                    continue;
                }
                const GraphRotor rotor = graph.get_angle_rotor(a, b, c, anchor);
                if (!rotor.is_rotatable || rotor.moving_nodes.empty()) continue;
                std::vector<int> moves;
                for (std::size_t n = 0; n < rotor.moving_nodes.size(); ++n) {
                    std::map<int, int>::const_iterator r =
                            row.find(rotor.moving_nodes[n]);
                    if (r != row.end()) moves.push_back(r->second);
                }
                if (moves.empty()) continue;
                // The arm that stays is whichever of the two is not moving.
                const int fixed_arm = rotor.moving == c ? a : c;
                geometry.add_angle(row[rotor.fixed], row[rotor.moving],
                                   row[fixed_arm], moves);
            }
        }
    }
    return geometry;
}

LinkerSamplingResult sample_linker(const std::string& mol2_path, int n_steps,
                                   int write_every, double step_size_dih,
                                   double step_size_ang, double temperature,
                                   int seed, const std::string& anchor_atom) {
    const Mol2Component component = read_mol2_component(mol2_path, "dye");
    const LinkerGeometry geometry =
            linker_geometry_from_mol2(mol2_path, anchor_atom);
    const IntramolecularEnergy evaluator(internal_topology_system(component));

    const int n_atoms = static_cast<int>(component.atoms.size());
    const int n_torsions = static_cast<int>(geometry.get_number_of_torsions());
    const int n_angles = static_cast<int>(geometry.get_number_of_angles());

    LinkerSamplingResult out;
    out.n_atoms = n_atoms;
    out.n_steps = std::max(0, n_steps);

    std::vector<double> config(n_torsions + n_angles, 0.0);
    std::vector<double> coordinates = geometry.apply(config);
    double energy = evaluator.evaluate(coordinates, n_atoms);

    boost::mt19937 rng(static_cast<boost::uint32_t>(seed));
    boost::normal_distribution<double> unit_normal(0.0, 1.0);
    boost::uniform_real<double> unit(0.0, 1.0);
    const double beta = 1.0 / (kb_kcal() * temperature);
    const int stride = std::max(1, write_every);

    for (int step = 0; step < out.n_steps; ++step) {
        std::vector<double> proposal(config.size());
        for (std::size_t i = 0; i < config.size(); ++i) {
            const double width = static_cast<int>(i) < n_torsions
                                         ? step_size_dih
                                         : step_size_ang;
            proposal[i] = config[i] + width * unit_normal(rng);
        }
        const std::vector<double> trial = geometry.apply(proposal);
        const double trial_energy = evaluator.evaluate(trial, n_atoms);
        // Metropolis: downhill always, uphill with the Boltzmann probability.
        if (trial_energy <= energy ||
            unit(rng) < std::exp(-beta * (trial_energy - energy))) {
            config = proposal;
            coordinates = trial;
            energy = trial_energy;
            ++out.n_accepted;
        }
        if ((step + 1) % stride == 0) {
            out.coordinates.insert(out.coordinates.end(), coordinates.begin(),
                                   coordinates.end());
            out.energies.push_back(energy);
            ++out.n_frames;
        }
    }
    return out;
}

RotamerLibrary generate_linker_rotamers(
        const std::string& mol2_path, int n_steps, int write_every,
        double step_size_dih, double step_size_ang, double cluster_threshold,
        double temperature, int seed, const std::string& anchor_atom,
        const std::vector<double>& protein_coords,
        const std::vector<std::string>& protein_elements, double mean_field_k,
        int mean_field_n_iter) {
    const LinkerSamplingResult sampled =
            sample_linker(mol2_path, n_steps, write_every, step_size_dih,
                          step_size_ang, temperature, seed, anchor_atom);
    const Mol2Component component = read_mol2_component(mol2_path, "dye");

    RotamerLibrary library;
    library.path = mol2_path;
    library.n_atoms = sampled.n_atoms;
    std::vector<Mol2Atom> atoms = component.atoms;
    std::sort(atoms.begin(), atoms.end(),
              [](const Mol2Atom& a, const Mol2Atom& b) {
                  return a.serial < b.serial;
              });
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        library.atom_names.push_back(atoms[i].atom_name);
        library.elements.push_back(atoms[i].element);
    }
    if (sampled.n_frames == 0) return library;

    // The leaders are the library; every frame then joins the nearest one.
    internal::OwnedIntView centres;
    cluster_frames_leader(const_cast<double*>(&sampled.coordinates[0]),
                          sampled.n_frames, sampled.n_atoms, 3,
                          cluster_threshold, &centres.data, &centres.size);
    internal::OwnedIntView assignment;
    assign_frames_to_clusters(const_cast<double*>(&sampled.coordinates[0]),
                              sampled.n_frames, sampled.n_atoms, 3,
                              centres.data, centres.size, &assignment.data,
                              &assignment.size);
    const int n_clusters = centres.size;

    // A cluster is worth what its members are worth together: a broad shallow
    // basin holds more of the ensemble than a narrow deep one.
    const std::vector<double> frame_weights =
            boltzmann_weights(sampled.energies, temperature);
    std::vector<int> assignments(assignment.data,
                                 assignment.data + assignment.size);
    std::vector<double> weights =
            cluster_weights(assignments, frame_weights, n_clusters);

    for (int i = 0; i < n_clusters; ++i) {
        const std::size_t offset =
                static_cast<std::size_t>(centres.data[i]) * sampled.n_atoms * 3;
        library.coords.insert(
                library.coords.end(), sampled.coordinates.begin() + offset,
                sampled.coordinates.begin() + offset + sampled.n_atoms * 3);
    }
    library.n_rotamers = n_clusters;

    // With a protein given, the site decides which conformers survive.
    if (!protein_coords.empty() && n_clusters > 0) {
        std::vector<std::string> elements = protein_elements;
        if (elements.empty()) {
            elements.assign(protein_coords.size() / 3, "C");
        }
        weights = rotamer_mean_field_weights(
                library.coords, weights, protein_coords, library.elements,
                elements, mean_field_k, mean_field_n_iter);
    }
    library.weights = weights;

    // Jump counts along the trajectory: what the frames did in order, which
    // is the only thing a static library cannot say.
    library.transitions.assign(
            static_cast<std::size_t>(n_clusters) * n_clusters, 0);
    for (int i = 0; i + 1 < assignment.size; ++i) {
        const int from = assignment.data[i], to = assignment.data[i + 1];
        if (from < 0 || to < 0 || from >= n_clusters || to >= n_clusters) {
            continue;
        }
        ++library.transitions[static_cast<std::size_t>(from) * n_clusters + to];
    }
    return library;
}

IMPBFF_END_NAMESPACE
