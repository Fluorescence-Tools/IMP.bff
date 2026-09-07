/**
 * \file SpectrumNode.cpp
 * \brief The nodes that produce an interleaved lifetime spectrum.
 *
 * The arithmetic here is ChiSurf's `elte2` / `e1tn` spectrum algebra and its
 * `calculcate_spectrum`, transcribed rather than re-derived: the two have to
 * agree to the last bit or a polarised fit lands somewhere else, and the
 * ordering of the Cartesian product is observable in the spectrum a caller
 * reads back.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/SpectrumNode.h>
#include <IMP/bff/Distributions.h>
#include <IMP/bff/PolymerChain.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <stdexcept>

IMPBFF_BEGIN_NAMESPACE

namespace {

void add_scalar_port(Node* node, const std::string& key, double value,
                     Port** slot) {
  std::shared_ptr<Port> port(new Port(value));
  node->add_input_port(key, port);
  *slot = port.get();
}

//! Write a spectrum to the output port keyed by the node's own name.
/** Fit transport, so it does not sanitise: a NaN lifetime has to reach
    `ChiSquared` and make the misfit infinite rather than be floored to
    `tiny`, which in a fit reads as a *good* fit near zero. */
void publish(Node* node, const std::vector<double>& spectrum) {
  const std::shared_ptr<Port> out = node->get_output_port(node->get_name());
  if (!out) {
    throw std::domain_error(
        "'" + node->get_name() +
        "' writes its spectrum to the output port keyed by its own name, "
        "which this node does not have");
  }
  out->set_sanitize(false);
  out->set_value_vector(spectrum);
}

void check_interleaved(const std::string& who,
                       const std::vector<double>& spectrum) {
  if (spectrum.size() < 2 || spectrum.size() % 2 != 0) {
    std::ostringstream m;
    m << who << ": the spectrum has " << spectrum.size()
      << " entries; it is interleaved (a0, t0, a1, t1, ...) so the count is "
         "even and at least two";
    throw std::domain_error(m.str());
  }
}

}  // namespace

void interleaved_scale_amplitudes(const std::vector<double>& spectrum,
                                  double n, std::vector<double>* out) {
  const std::size_t pairs = spectrum.size() / 2;
  for (std::size_t i = 0; i < pairs; ++i) {
    out->push_back(spectrum[2 * i] * n);
    out->push_back(spectrum[2 * i + 1]);
  }
}

std::vector<double> interleaved_product(const std::vector<double>& first,
                                        const std::vector<double>& second) {
  const std::size_t n1 = first.size() / 2, n2 = second.size() / 2;
  std::vector<double> out;
  out.reserve(2 * n1 * n2);
  for (std::size_t i = 0; i < n1; ++i) {
    for (std::size_t j = 0; j < n2; ++j) {
      out.push_back(first[2 * i] * second[2 * j]);
      // 1 / (1/t1 + 1/t2). Written this way, and not as the algebraically
      // equal t1 t2 / (t1 + t2), because that spelling is what ChiSurf
      // evaluates and the two do not round identically.
      out.push_back(1.0 / (1.0 / first[2 * i + 1] + 1.0 / second[2 * j + 1]));
    }
  }
  return out;
}

// ------------------------------------------------------ LifetimeSpectrumNode

// Empty for the reason `TcspcDecay`'s is: `add_port` reaches for
// `shared_from_this()`, so a node cannot own ports until something owns it.
LifetimeSpectrumNode::LifetimeSpectrumNode(const std::string& name)
    : Node(name) {}

void LifetimeSpectrumNode::set_number_of_lifetimes(int n) {
  if (n <= 0) {
    throw std::domain_error(
        "LifetimeSpectrumNode::set_number_of_lifetimes: a spectrum has at "
        "least one species");
  }
  lifetime_ports_.clear();
  lifetime_ports_.reserve(static_cast<std::size_t>(2 * n));
  for (int i = 0; i < n; ++i) {
    std::ostringstream a, t;
    a << "a" << i;
    t << "t" << i;
    Port* amplitude = nullptr;
    Port* lifetime = nullptr;
    // Ports cannot be removed from a node, so a second call with a smaller
    // count leaves the surplus ports visible but unread -- same caveat, and
    // same remedy (build the node once), as `TcspcDecay`.
    if (get_input_port(a.str())) {
      amplitude = get_input_port(a.str()).get();
      lifetime = get_input_port(t.str()).get();
    } else {
      add_scalar_port(this, a.str(), 1.0, &amplitude);
      add_scalar_port(this, t.str(), 1.0, &lifetime);
    }
    lifetime_ports_.push_back(amplitude);
    lifetime_ports_.push_back(lifetime);
  }
  n_lifetimes_ = n;
  spectrum_.assign(static_cast<std::size_t>(2 * n), 0.0);
  set_valid(false);
}

void LifetimeSpectrumNode::evaluate() {
  if (n_lifetimes_ <= 0) {
    throw std::domain_error("LifetimeSpectrumNode '" + get_name() +
                            "' has no species; call "
                            "set_number_of_lifetimes() first");
  }
  spectrum_.resize(static_cast<std::size_t>(2 * n_lifetimes_));
  double sum = 0.0;
  for (int i = 0; i < n_lifetimes_; ++i) {
    const std::size_t k = static_cast<std::size_t>(2 * i);
    double amplitude = lifetime_ports_[k]->get_value();
    if (absolute_amplitudes_) amplitude = std::fabs(amplitude);
    spectrum_[k] = amplitude;
    sum += amplitude;
    // Absolute, for the reason `TcspcDecay` takes it: a lifetime walked
    // through zero is a growing exponential, not a decay.
    spectrum_[k + 1] = std::fabs(lifetime_ports_[k + 1]->get_value());
  }
  if (normalize_amplitudes_) {
    const double scale = std::fabs(sum);
    for (int i = 0; i < n_lifetimes_; ++i) {
      spectrum_[static_cast<std::size_t>(2 * i)] /= scale;
    }
  }
  publish(this, spectrum_);
  set_valid(true);
}

// ------------------------------------------------------- AnisotropySpectrum

AnisotropySpectrum::AnisotropySpectrum(const std::string& name) : Node(name) {}

void AnisotropySpectrum::set_number_of_rotations(int n) {
  if (n < 0) {
    throw std::domain_error(
        "AnisotropySpectrum::set_number_of_rotations: a negative count");
  }
  if (spectrum_port_ == nullptr) {
    std::shared_ptr<Port> incoming(new Port(std::vector<double>(2, 1.0)));
    incoming->set_sanitize(false);
    add_input_port(spectrum_port_key(), incoming);
    spectrum_port_ = incoming.get();

    add_scalar_port(this, "r0", 0.38, &r0_port_);
    add_scalar_port(this, "g", 1.0, &g_port_);
    add_scalar_port(this, "l1", 0.0, &l1_port_);
    add_scalar_port(this, "l2", 0.0, &l2_port_);
  }
  rotation_ports_.clear();
  rotation_ports_.reserve(static_cast<std::size_t>(2 * n));
  for (int i = 0; i < n; ++i) {
    std::ostringstream b, rho;
    b << "b" << i;
    rho << "rho" << i;
    Port* amplitude = nullptr;
    Port* time = nullptr;
    if (get_input_port(b.str())) {
      amplitude = get_input_port(b.str()).get();
      time = get_input_port(rho.str()).get();
    } else {
      add_scalar_port(this, b.str(), 1.0, &amplitude);
      add_scalar_port(this, rho.str(), 1.0, &time);
    }
    rotation_ports_.push_back(amplitude);
    rotation_ports_.push_back(time);
  }
  n_rotations_ = n;
  set_valid(false);
}

void AnisotropySpectrum::set_polarization(Polarization p) {
  polarization_ = p;
  set_valid(false);
}

void AnisotropySpectrum::set_polarization_name(const std::string& name) {
  std::string lowered;
  lowered.reserve(name.size());
  for (char c : name) {
    lowered.push_back(static_cast<char>(
        std::tolower(static_cast<unsigned char>(c))));
  }
  if (lowered == "vm") {
    set_polarization(VM);
  } else if (lowered == "vv") {
    set_polarization(VV);
  } else if (lowered == "vh") {
    set_polarization(VH);
  } else if (lowered == "vv/vh" || lowered == "vv_vh") {
    set_polarization(VV_VH);
  } else {
    // Not defaulted to VM: that would return the spectrum unchanged, which
    // is a model with no anisotropy at all rather than an error anyone sees.
    throw std::domain_error(
        "AnisotropySpectrum::set_polarization_name: '" + name +
        "' is not one of vm, vv, vh, vv/vh");
  }
}

std::string AnisotropySpectrum::get_polarization_name() const {
  switch (polarization_) {
    case VV:
      return "vv";
    case VH:
      return "vh";
    case VV_VH:
      return "vv/vh";
    default:
      return "vm";
  }
}

void AnisotropySpectrum::build_rotation_spectrum() {
  rotation_.resize(static_cast<std::size_t>(2 * n_rotations_));
  double sum = 0.0;
  for (int i = 0; i < n_rotations_; ++i) {
    const std::size_t k = static_cast<std::size_t>(2 * i);
    const double amplitude =
        std::fabs(rotation_ports_[k]->get_value());
    rotation_[k] = amplitude;
    sum += amplitude;
    rotation_[k + 1] = std::fabs(rotation_ports_[k + 1]->get_value());
  }
  // b_i <- |b_i| / sum|b| * r0, so `r0` alone sets r(0) and the b_i only
  // divide it up. ChiSurf normalises in the getter and writes the result
  // back to the parameters; the write-back is the application's business,
  // the normalisation is the model's.
  const double r0 = r0_port_->get_value();
  for (int i = 0; i < n_rotations_; ++i) {
    rotation_[static_cast<std::size_t>(2 * i)] *= r0 / sum;
  }
}

void AnisotropySpectrum::evaluate() {
  if (spectrum_port_ == nullptr) {
    throw std::domain_error("AnisotropySpectrum '" + get_name() +
                            "' has no ports; call set_number_of_rotations() "
                            "first, which is what builds them");
  }
  // Resolved by key rather than through the cached pointer: a caller who
  // hands in their *own* port for this key -- which is how the misfit node
  // downstream is wired, so it is the shape a builder copies -- replaces the
  // one built here, and the cached raw pointer would then be reading a port
  // nothing owns any more. The scalars keep their pointers; only this one is
  // ever replaced, and one map lookup per evaluation is not a cost a fit can
  // measure.
  const std::shared_ptr<Port> in = get_input_port(spectrum_port_key());
  if (!in) {
    throw std::domain_error("AnisotropySpectrum '" + get_name() +
                            "' has no '" + spectrum_port_key() + "' port");
  }
  const std::vector<double>& incoming = in->get_values_ref();
  check_interleaved("AnisotropySpectrum '" + get_name() + "'", incoming);

  if (polarization_ == VM || n_rotations_ <= 0) {
    // Magic angle has no anisotropy contribution, and building the rotation
    // spectrum first would be waste -- not cheap waste, since the caller's
    // normalisation writes back per component on every evaluation.
    spectrum_ = incoming;
    rotation_.clear();
    publish(this, spectrum_);
    set_valid(true);
    return;
  }

  build_rotation_spectrum();

  // The anisotropy decay times the fluorescence decay, as a spectrum. The
  // rotation spectrum is the *first* argument, so the product runs
  // rotation-major -- the order ChiSurf's `elte2(a, f)` produces, and the
  // order is visible to anyone who reads the spectrum back.
  const std::vector<double> product = interleaved_product(rotation_, incoming);

  const double g = g_port_->get_value();
  const double l1 = l1_port_->get_value();
  const double l2 = l2_port_->get_value();

  //  f_VV = f (1 + 2 r),  f_VH = f (1 - r) / G, as spectra: the unmodified
  //  spectrum followed by the product scaled by +2 (resp. -1), the whole of
  //  VH then divided by G.
  std::vector<double> vv(incoming);
  interleaved_scale_amplitudes(product, 2.0, &vv);

  std::vector<double> vh_unscaled(incoming);
  interleaved_scale_amplitudes(product, -1.0, &vh_unscaled);
  std::vector<double> vh;
  vh.reserve(vh_unscaled.size());
  // G divides. Multiplying here instead put the two forward models a factor
  // G^2 apart in VH, which is a bug this stack has already paid for once.
  interleaved_scale_amplitudes(vh_unscaled, 1.0 / g, &vh);

  // A mixed channel is the *union* of the scaled VV and VH components, not
  // their element-wise sum: adding would also add the time constants and so
  // double the decay times.
  spectrum_.clear();
  if (polarization_ == VV || polarization_ == VV_VH) {
    interleaved_scale_amplitudes(vv, 1.0 - l1, &spectrum_);
    interleaved_scale_amplitudes(vh, l1, &spectrum_);
  }
  if (polarization_ == VH || polarization_ == VV_VH) {
    interleaved_scale_amplitudes(vv, l2, &spectrum_);
    interleaved_scale_amplitudes(vh, 1.0 - l2, &spectrum_);
  }
  publish(this, spectrum_);
  set_valid(true);
}

std::string AnisotropySpectrum::describe() const {
  std::ostringstream out;
  out << "polarization   : " << get_polarization_name() << "\n"
      << "rotations      : " << n_rotations_ << "\n"
      << "species out    : " << spectrum_.size() / 2 << "\n";
  return out.str();
}

// -------------------------------------------------------- GaussianDistances

// ------------------------------------------------------- PolymerDistances

PolymerDistances::PolymerDistances(const std::string& name) : Node(name) {}

void PolymerDistances::set_mode(const std::string& mode) {
  if (!mode_.empty()) {
    throw std::domain_error("PolymerDistances '" + get_name() +
                            "': the mode is set once");
  }
  parameter_ports_.clear();
  Port* p = nullptr;
  if (mode == "worm_like_chain" || mode == "worm_like_chain_linker") {
    add_scalar_port(this, "chain_length", 100.0, &p);
    parameter_ports_.push_back(p);
    add_scalar_port(this, "persistence_length", 30.0, &p);
    parameter_ports_.push_back(p);
    // Both modes carry the linker width. Without the linker it is inert --
    // evaluate() never reads it -- which mirrors chisurf exactly: the
    // model's `w` is a fitting parameter whether or not the linker is on,
    // and with it off the numpy path fits an inert parameter too. A port
    // the optimiser can claim is what keeps a free-but-inert `w` from
    // refusing the whole graph.
    add_scalar_port(this, "sigma_linker", 6.0, &p);
    parameter_ports_.push_back(p);
  } else if (mode == "saw_nu") {
    add_scalar_port(this, "r_rms", 50.0, &p);
    parameter_ports_.push_back(p);
    add_scalar_port(this, "nu", 0.588, &p);
    parameter_ports_.push_back(p);
  } else if (mode == "ising_chain") {
    const char* keys[] = {"n_residues", "b_structured", "b_unstructured",
                          "coupling", "field"};
    const double defaults[] = {50.0, 4.0, 8.0, 1.5, 0.0};
    for (int i = 0; i < 5; ++i) {
      add_scalar_port(this, keys[i], defaults[i], &p);
      parameter_ports_.push_back(p);
    }
  } else {
    throw std::domain_error("PolymerDistances '" + get_name() +
                            "': unknown mode '" + mode + "'");
  }
  mode_ = mode;
  set_valid(false);
}

void PolymerDistances::set_axis(const std::vector<double>& axis) {
  if (axis.size() < 2) {
    throw std::domain_error("PolymerDistances '" + get_name() +
                            "': a distance axis needs at least two points");
  }
  axis_ = axis;
  set_valid(false);
}

void PolymerDistances::set_axis_array(double* in_axis, int n_axis) {
  set_axis(std::vector<double>(in_axis, in_axis + n_axis));
}

void PolymerDistances::evaluate() {
  if (mode_.empty()) {
    throw std::domain_error("PolymerDistances '" + get_name() +
                            "' has no mode; call set_mode() first");
  }
  if (axis_.empty()) {
    throw std::domain_error("PolymerDistances '" + get_name() +
                            "' has no distance axis");
  }
  // The kernels publish malloc'd views (their numpy contract); copy and
  // free. They are the same functions chisurf's rdf.py forwarders call, so
  // the graph path and the Python path evaluate one implementation.
  double* view = nullptr;
  int n_view = 0;
  if (mode_ == "worm_like_chain" || mode_ == "worm_like_chain_linker") {
    const double chain_length = parameter_ports_[0]->get_value();
    const double persistence_length = parameter_ports_[1]->get_value();
    if (!(chain_length > 0.0)) {
      throw std::domain_error("PolymerDistances '" + get_name() +
                              "': chain_length is not positive");
    }
    // The kernel takes the dimensionless kappa; the ports carry what the
    // model fits. Same derivation as chisurf's model property.
    const double kappa = persistence_length / chain_length;
    if (mode_ == "worm_like_chain") {
      // `distance = false`, always: chisurf's forwarder has never applied
      // the r^2 factor its signature advertises -- the flag was accepted
      // and dropped on the floor, every fit in the stack was made against
      // that behaviour, and preserving the answer is the port's contract
      // (rdf.py says the same over the same call; owner decision pending
      // in PRD-105 on the flag itself).
      worm_like_chain(axis_, kappa, chain_length, true, false, &view, &n_view);
    } else {
      worm_like_chain_linker(axis_, kappa, chain_length,
                             parameter_ports_[2]->get_value(), true,
                             &view, &n_view);
    }
  } else if (mode_ == "saw_nu") {
    saw_nu(axis_, parameter_ports_[0]->get_value(),
           parameter_ports_[1]->get_value(), 1.1615, &view, &n_view);
  } else {  // ising_chain
    const int n_residues = static_cast<int>(
        std::lround(parameter_ports_[0]->get_value()));
    ising_chain(axis_, n_residues, parameter_ports_[1]->get_value(),
                parameter_ports_[2]->get_value(),
                parameter_ports_[3]->get_value(),
                parameter_ports_[4]->get_value(), n_k_, &view, &n_view);
  }
  if (view == nullptr || n_view != static_cast<int>(axis_.size())) {
    std::free(view);
    throw std::domain_error("PolymerDistances '" + get_name() +
                            "': the kernel returned no distribution");
  }
  spectrum_.resize(2 * axis_.size());
  for (std::size_t j = 0; j < axis_.size(); ++j) {
    spectrum_[2 * j] = view[j];
    spectrum_[2 * j + 1] = axis_[j];
  }
  std::free(view);
  publish(this, spectrum_);
  set_valid(true);
}

// ------------------------------------------------------- GaussianDistances

GaussianDistances::GaussianDistances(const std::string& name) : Node(name) {}

void GaussianDistances::set_number_of_components(int n) {
  if (n <= 0) {
    throw std::domain_error(
        "GaussianDistances::set_number_of_components: a distribution has at "
        "least one component");
  }
  component_ports_.clear();
  component_ports_.reserve(static_cast<std::size_t>(4 * n));
  static const char* kNames[4] = {"mean", "sigma", "shape", "amplitude"};
  static const double kDefaults[4] = {50.0, 6.0, 0.0, 1.0};
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < 4; ++k) {
      std::ostringstream key;
      key << kNames[k] << i;
      Port* port = nullptr;
      if (get_input_port(key.str())) {
        port = get_input_port(key.str()).get();
      } else {
        add_scalar_port(this, key.str(), kDefaults[k], &port);
      }
      component_ports_.push_back(port);
    }
  }
  n_components_ = n;
  set_valid(false);
}

void GaussianDistances::set_axis(const std::vector<double>& axis) {
  if (axis.size() < 2) {
    throw std::domain_error(
        "GaussianDistances::set_axis: a distance axis needs at least two "
        "points");
  }
  axis_ = axis;
  density_.assign(axis_.size(), 0.0);
  set_valid(false);
}

void GaussianDistances::set_axis_array(double* in_axis, int n_axis) {
  set_axis(std::vector<double>(in_axis, in_axis + n_axis));
}

void GaussianDistances::evaluate() {
  if (axis_.empty()) {
    throw std::domain_error("GaussianDistances '" + get_name() +
                            "' has no distance axis");
  }
  const std::size_t n = axis_.size();
  std::vector<double> means, sigmas, shapes, amplitudes;
  means.reserve(n_components_);
  for (int c = 0; c < n_components_; ++c) {
    const std::size_t base = static_cast<std::size_t>(4 * c);
    // The mean and the amplitude are taken absolute, as ChiSurf's getters do:
    // a distance is positive, and a fit that walks one through zero would
    // otherwise put the distribution on the wrong side of the axis.
    means.push_back(std::fabs(component_ports_[base]->get_value()));
    sigmas.push_back(component_ports_[base + 1]->get_value());
    shapes.push_back(component_ports_[base + 2]->get_value());
    amplitudes.push_back(std::fabs(component_ports_[base + 3]->get_value()));
  }

  // One kernel, shared with the free function ChiSurf's own
  // `Gaussians.distribution` calls -- so the graph path and the Python path
  // cannot drift apart. `normalize_components` follows the branch, because the
  // reference is asymmetric between its two; see `Distributions.h`.
  //
  // ChiSurf normalises the *weights* in `Gaussians.amplitude` and the combined
  // density again in `combine_distributions`; the first is a division by a
  // constant that the second undoes, so normalising once here is the same
  // number.
  density_ = gaussian_distance_mixture_impl(
      axis_, means, sigmas, shapes, amplitudes,
      distance_between_gaussians_ ? GAUSSIAN_MIXTURE_DISTANCE_BETWEEN_GAUSSIANS
                                  : GAUSSIAN_MIXTURE_GENERALIZED_NORMAL,
      !distance_between_gaussians_, true);

  spectrum_.resize(2 * n);
  for (std::size_t j = 0; j < n; ++j) {
    spectrum_[2 * j] = density_[j];
    spectrum_[2 * j + 1] = axis_[j];
  }
  publish(this, spectrum_);
  set_valid(true);
}

// ------------------------------------------------------------ FretSpectrum

FretSpectrum::FretSpectrum(const std::string& name) : Node(name) {}

void FretSpectrum::build_ports() {
  if (donor_port_ != nullptr) return;
  std::shared_ptr<Port> donor(new Port(std::vector<double>{1.0, 4.0}));
  donor->set_sanitize(false);
  add_input_port(donor_port_key(), donor);
  donor_port_ = donor.get();

  std::shared_ptr<Port> distances(new Port(std::vector<double>{1.0, 52.0}));
  distances->set_sanitize(false);
  add_input_port(distance_port_key(), distances);
  distance_port_ = distances.get();

  add_scalar_port(this, "x_donly", 0.0, &x_donly_port_);
  add_scalar_port(this, "forster_radius", 52.0, &forster_radius_port_);
  add_scalar_port(this, "tau0", 4.0, &tau0_port_);
  add_scalar_port(this, "kappa2", 2.0 / 3.0, &kappa2_port_);
  set_valid(false);
}

void FretSpectrum::evaluate() {
  if (donor_port_ == nullptr) {
    throw std::domain_error("FretSpectrum '" + get_name() +
                            "' has no ports; call build_ports() first");
  }
  // By key, not through the cached pointer, for the reason
  // `AnisotropySpectrum` gives: these two are the ports a caller replaces.
  const std::shared_ptr<Port> donor_in = get_input_port(donor_port_key());
  const std::shared_ptr<Port> distance_in = get_input_port(distance_port_key());
  if (!donor_in || !distance_in) {
    throw std::domain_error("FretSpectrum '" + get_name() +
                            "' is missing one of its two spectrum ports");
  }
  const std::vector<double>& donor = donor_in->get_values_ref();
  const std::vector<double>& distances = distance_in->get_values_ref();
  check_interleaved("FretSpectrum '" + get_name() + "' (donor)", donor);
  check_interleaved("FretSpectrum '" + get_name() + "' (distances)",
                    distances);

  const double x_donly = std::fabs(x_donly_port_->get_value());
  const double forster_radius = forster_radius_port_->get_value();
  const double tau0 = tau0_port_->get_value();
  const double kappa2 = kappa2_port_->get_value();
  if (!(tau0 > 0.0)) {
    throw std::domain_error("FretSpectrum '" + get_name() +
                            "': tau0 is not positive");
  }

  // The donor as *rates*, which is the space the two spectra combine in.
  const std::size_t n_donor = donor.size() / 2;
  std::vector<double> donor_rates(2 * n_donor);
  for (std::size_t i = 0; i < n_donor; ++i) {
    donor_rates[2 * i] = donor[2 * i];
    donor_rates[2 * i + 1] = 1.0 / donor[2 * i + 1];
  }

  // Each distance as a transfer rate: k = 3/2 kappa^2 / tau0 * (R0/r)^6.
  const std::size_t n_distances = distances.size() / 2;
  std::vector<double> fret_rates(2 * n_distances);
  const double prefactor = 1.5 * kappa2 / tau0;
  for (std::size_t i = 0; i < n_distances; ++i) {
    const double ratio = forster_radius / distances[2 * i + 1];
    const double ratio3 = ratio * ratio * ratio;
    fret_rates[2 * i] = distances[2 * i];
    fret_rates[2 * i + 1] = prefactor * ratio3 * ratio3;
  }

  // Rates add: a quenched species decays at its own rate plus the transfer
  // rate, so the quenched spectrum is the Cartesian product with the
  // amplitudes multiplied. FRET-major, which is the order ChiSurf's
  // `ere2(fret, donor)` produces.
  std::vector<double> combined;
  combined.reserve(2 * (n_distances * n_donor + n_donor));
  for (std::size_t i = 0; i < n_distances; ++i) {
    for (std::size_t j = 0; j < n_donor; ++j) {
      combined.push_back(fret_rates[2 * i] * donor_rates[2 * j] *
                         (1.0 - x_donly));
      combined.push_back(fret_rates[2 * i + 1] + donor_rates[2 * j + 1]);
    }
  }
  // The donor-only fraction, appended unconditionally -- including at
  // `x_donly == 0`, where it carries zero amplitude. A node that sometimes
  // returns a shorter spectrum is a second code path for no gain, and the
  // length is observable.
  for (std::size_t j = 0; j < n_donor; ++j) {
    combined.push_back(donor_rates[2 * j] * x_donly);
    combined.push_back(donor_rates[2 * j + 1]);
  }

  // Back to lifetimes, which is what the instrument node reconvolves.
  spectrum_.resize(combined.size());
  for (std::size_t i = 0; i < combined.size() / 2; ++i) {
    spectrum_[2 * i] = combined[2 * i];
    spectrum_[2 * i + 1] = 1.0 / combined[2 * i + 1];
  }
  publish(this, spectrum_);
  set_valid(true);
}

std::string FretSpectrum::describe() const {
  std::ostringstream out;
  out << "species out    : " << spectrum_.size() / 2 << "\n"
      << "R0 / tau0      : " << forster_radius_port_->get_value() << " / "
      << tau0_port_->get_value() << "\n";
  return out.str();
}

IMPBFF_END_NAMESPACE
