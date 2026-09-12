#ifndef IMPBFF_INTERNAL_SPECTRUMNODEHELPERS_H
#define IMPBFF_INTERNAL_SPECTRUMNODEHELPERS_H

#include <IMP/bff/GraphNode.h>
#include <IMP/bff/GraphPort.h>
#include <sstream>
#include <stdexcept>

IMPBFF_BEGIN_NAMESPACE
namespace spectrum_node_detail {

inline void add_scalar_port(GraphNode* node, const std::string& key, double value,
                     GraphPort** slot) {
  std::shared_ptr<GraphPort> port(new GraphPort(value));
  node->add_input_port(key, port);
  *slot = port.get();
}

//! Write a spectrum to the output port keyed by the node's own name.
/** Fit transport, so it does not sanitise: a NaN lifetime has to reach
    `FitChiSquared` and make the misfit infinite rather than be floored to
    `tiny`, which in a fit reads as a *good* fit near zero. */
inline void publish(GraphNode* node, const std::vector<double>& spectrum) {
  const std::shared_ptr<GraphPort> out = node->get_output_port(node->get_name());
  if (!out) {
    throw std::domain_error(
        "'" + node->get_name() +
        "' writes its spectrum to the output port keyed by its own name, "
        "which this node does not have");
  }
  out->set_sanitize(false);
  out->set_value_vector(spectrum);
}

inline void check_interleaved(const std::string& who,
                       const std::vector<double>& spectrum) {
  if (spectrum.size() < 2 || spectrum.size() % 2 != 0) {
    std::ostringstream m;
    m << who << ": the spectrum has " << spectrum.size()
      << " entries; it is interleaved (a0, t0, a1, t1, ...) so the count is "
         "even and at least two";
    throw std::domain_error(m.str());
  }
}

}  // namespace spectrum_node_detail
IMPBFF_END_NAMESPACE

#endif
