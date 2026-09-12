/**
 * \file FitJointChiSquared.cpp
 * \brief One misfit over several datasets: the grouping, as a node.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#include <IMP/bff/FitJointChiSquared.h>

#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

IMPBFF_BEGIN_NAMESPACE

FitJointChiSquared::FitJointChiSquared(const std::string& name) : GraphNode(name) {}

int FitJointChiSquared::add_member(std::shared_ptr<GraphNode> member,
                                const std::string& residual_key) {
  if (!member) {
    throw std::domain_error(
        "FitJointChiSquared::add_member: the member is a null pointer");
  }
  const std::shared_ptr<GraphPort> source = member->get_output_port(residual_key);
  if (!source) {
    throw std::domain_error(
        "FitJointChiSquared::add_member: member '" + member->get_name() +
        "' has no output port '" + residual_key + "' carrying its residuals");
  }

  // The block is an ordinary linked input, so `GraphNode::update()` is what walks
  // the group: it evaluates each member whose node is out of date, copies the
  // residuals in, and only then evaluates this node. Nothing here has to know
  // how deep a member's own graph goes.
  std::ostringstream key;
  key << "block_" << members_.size();
  std::shared_ptr<GraphPort> block(new GraphPort(std::vector<double>(1, 0.0)));
  // A member's residuals are fit transport: a NaN must survive the copy.
  block->set_sanitize(false);
  source->set_sanitize(false);
  block->set_link(source);
  add_input_port(key.str(), block);

  members_.push_back(member);
  blocks_.push_back(block);
  return static_cast<int>(members_.size()) - 1;
}

std::shared_ptr<GraphNode> FitJointChiSquared::get_member(int index) const {
  if (index < 0 || static_cast<std::size_t>(index) >= members_.size()) {
    std::ostringstream m;
    m << "FitJointChiSquared::get_member: index " << index << " of "
      << members_.size() << " members";
    throw std::domain_error(m.str());
  }
  return members_[static_cast<std::size_t>(index)];
}

std::vector<std::string> FitJointChiSquared::get_member_names() const {
  std::vector<std::string> names;
  names.reserve(members_.size());
  for (const std::shared_ptr<GraphNode>& m : members_) names.push_back(m->get_name());
  return names;
}

std::vector<int> FitJointChiSquared::get_block_offsets() const {
  std::vector<int> offsets;
  offsets.reserve(block_sizes_.size());
  int at = 0;
  for (int n : block_sizes_) {
    offsets.push_back(at);
    at += n;
  }
  return offsets;
}

double FitJointChiSquared::get_chi2r(int n_free) const {
  const double dof =
      static_cast<double>(wres_.size()) - static_cast<double>(n_free) - 1.0;
  return chi2_ / dof;
}

void FitJointChiSquared::evaluate() {
  block_sizes_.clear();
  block_sizes_.reserve(blocks_.size());
  std::size_t total = 0;
  for (const std::shared_ptr<GraphPort>& block : blocks_) {
    const std::size_t n = block->get_values_ref().size();
    block_sizes_.push_back(static_cast<int>(n));
    total += n;
  }

  wres_.clear();
  wres_.reserve(total);
  for (const std::shared_ptr<GraphPort>& block : blocks_) {
    const std::vector<double>& v = block->get_values_ref();
    wres_.insert(wres_.end(), v.begin(), v.end());
  }

  chi2_ = 0.0;
  for (double r : wres_) chi2_ += r * r;
  // A NaN anywhere makes the whole group infinitely bad, which is
  // `FitChiSquared`'s convention and what makes a sampler reject rather than
  // propagate the NaN into the posterior.
  if (std::isnan(chi2_)) chi2_ = std::numeric_limits<double>::infinity();

  const std::shared_ptr<GraphPort> out = get_output_port(get_name());
  if (!out) {
    throw std::domain_error(
        "FitJointChiSquared '" + get_name() +
        "' writes chi-square to the output port keyed by its own name, "
        "which this node does not have");
  }
  out->set_value(chi2_);

  // The concatenated residuals, when the graph asked for them. Absent by
  // default, exactly as in `FitChiSquared`: a sampler wants the scalar and
  // would otherwise pay for a copy of every dataset's residuals per move.
  const std::shared_ptr<GraphPort> res = get_output_port(residuals_key_);
  if (res) {
    res->set_sanitize(false);
    res->set_value_vector(wres_);
  }
  set_valid(true);
}

std::string FitJointChiSquared::describe() const {
  std::ostringstream out;
  out << "members        : " << members_.size() << "\n";
  const std::vector<int> offsets = get_block_offsets();
  for (std::size_t i = 0; i < members_.size(); ++i) {
    out << "  [" << i << "] " << members_[i]->get_name();
    if (i < block_sizes_.size()) {
      out << " -- " << block_sizes_[i] << " residuals at " << offsets[i];
    }
    out << "\n";
  }
  out << "residuals      : " << wres_.size() << "\n"
      << "chi2           : " << chi2_ << "\n";
  return out.str();
}

IMPBFF_END_NAMESPACE
