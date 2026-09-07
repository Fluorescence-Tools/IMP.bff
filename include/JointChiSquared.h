/**
 * \file IMP/bff/JointChiSquared.h
 * \brief One misfit over several datasets: the grouping, as a node.
 *
 * A fit of several datasets at once is two things, and bff already had one
 * of them. **Sharing a parameter** between datasets is `Port::set_link` --
 * the member models' ports follow one master port, and the link graph is
 * enforced acyclic -- so the coupling that makes a joint fit *joint* has
 * been on this side of the boundary since the Port runtime landed. What was
 * missing is the other half: **one objective over all of them.**
 *
 * That is all this class is. Its residual is the members' residuals laid end
 * to end, and its chi-square is their sum, which is exactly what minimising
 * a joint fit means: one Levenberg-Marquardt step moves the shared
 * parameters using the curvature of every dataset at once, rather than each
 * dataset in turn hoping they agree.
 *
 * It is a `Node`, so the members are reached the ordinary way: each member's
 * residual output port is *linked* to one of this node's input ports, and
 * `Node::update()` therefore evaluates the whole tree -- every member's
 * model, every member's misfit, then this -- from one call, with nothing
 * crossing into Python. A `Minimizer` pointed at this node optimises the
 * group.
 *
 * **Not called a fit group**, and not shaped like one. ChiSurf's `FitGroup`
 * is a container of `Fit` objects that also owns a selection, a result
 * history, plots and a run policy; this is the arithmetic underneath that
 * and nothing else. The layering rule in `AGENTS.md` puts the objective here
 * and the bookkeeping in the application, so the name says what the object
 * *is* -- the joint chi-square -- rather than which application assembles
 * it.
 *
 * The members do not have to be `ChiSquared` nodes. Anything presenting a
 * residual vector on an output port qualifies, including a Python `Node`
 * director wrapping a model this library cannot represent, so a group may
 * mix representable and unrepresentable members and still take one step.
 *
 * \see ChiSquared, Minimizer, Expression, Port
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_JOINTCHISQUARED_H
#define IMPBFF_JOINTCHISQUARED_H

#include <IMP/bff/bff_config.h>

#include <memory>
#include <string>
#include <vector>

#include <IMP/bff/Node.h>
#include <IMP/bff/Port.h>

IMPBFF_BEGIN_NAMESPACE

//! The misfit of several datasets at once: their residuals, end to end.
class IMPBFFEXPORT JointChiSquared : public Node {
 public:
  explicit JointChiSquared(const std::string& name = "joint");

  //! Add a member, linking its residual port to a new input of this node.
  /**
      \param member the node producing one dataset's residuals
      \param residual_key the member's output port carrying them
      \return the index of the block this member occupies

      The member keeps its own data, fit range, mask and noise model -- a
      group whose members share none of those is the normal case, not an
      awkward one. Only the *parameters* are shared, and they are shared by
      linking ports, which happens outside this class.

      Members are evaluated in the order they were added, and the joint
      residual is in that order, so a caller can map a block back to the
      dataset that produced it.
   */
  int add_member(std::shared_ptr<Node> member,
                 const std::string& residual_key = "residuals");

  //! How many members the group holds.
  unsigned int get_number_of_members() const {
    return static_cast<unsigned int>(members_.size());
  }

  //! One member, by the index `add_member` returned.
  /** Reached one at a time because `std::vector<std::shared_ptr<Node> >` is
      not a template this module may name (see `IMP_bff.types.i`), and
      wrapping it anyway leaks -- SWIG finds no destructor for it. */
  std::shared_ptr<Node> get_member(int index) const;

  //! The members' names, in the order they were added.
  std::vector<std::string> get_member_names() const;

  //! How long each member's residual block was in the last evaluation.
  /** Empty until the node has been evaluated. A block length can change
      between evaluations -- a member is free to truncate against a short
      model curve, which `ChiSquared` does -- so this describes the last
      evaluation rather than a promise about the next. */
  std::vector<int> get_block_sizes() const { return block_sizes_; }

  //! Where each member's block starts in the joint residual.
  std::vector<int> get_block_offsets() const;

  //! The key of the output port carrying the concatenated residuals.
  void set_residuals_port_key(const std::string& key) { residuals_key_ = key; }
  const std::string& get_residuals_port_key() const { return residuals_key_; }

  //! The joint residuals from the last evaluation.
  const std::vector<double>& get_weighted_residuals() const { return wres_; }

  //! The sum of the members' chi-squares, from the last evaluation.
  double get_chi2() const { return chi2_; }

  //! `chi2 / (n_residuals - n_free - 1)`, over the whole group.
  double get_chi2r(int n_free) const;

  //! Number of residuals the last evaluation produced, over all members.
  unsigned int get_number_of_residuals() const {
    return static_cast<unsigned int>(wres_.size());
  }

  //! Read every member's residuals, concatenate, write the sum out.
  void evaluate() override;

  std::string describe() const;

 private:
  std::vector<std::shared_ptr<Node> > members_;
  //! This node's input ports, one per member, in member order.
  std::vector<std::shared_ptr<Port> > blocks_;
  std::vector<int> block_sizes_;
  std::vector<double> wres_;
  std::string residuals_key_ = "residuals";
  double chi2_ = 0.0;
};

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_JOINTCHISQUARED_H
