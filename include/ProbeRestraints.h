/**
 *  \file IMP/bff/ProbeRestraints.h
 *  \brief Chi-squared scoring of labelling data, with and without volumes.
 *
 * Two restraints over the same asymmetric chi-squared, differing only in what
 * they take as the *model* distance:
 *
 * - #IMP::bff::SimpleProbeNetworkRestraint scores a network of accessible-volume
 *   pair distances. Unlike #IMP::bff::ProbeNetworkRestraint it needs no fps.json —
 *   volumes and measurements are added programmatically.
 * - #IMP::bff::DirectProbeRestraint scores attachment-atom distances
 *   directly, with no volume at all. It is for MCMC inner loops, where the AV
 *   build is the whole cost.
 *
 * \authors Thomas-Otavio Peulen
 *  Copyright 2007-2026 IMP Inventors. All rights reserved.
 *
 */
#ifndef IMPBFF_PROBERESTRAINTS_H
#define IMPBFF_PROBERESTRAINTS_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/AVModel.h>

#include <IMP/bff/Base.h>

#include <map>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! An experimental distance measurement between a pair of accessible volumes.
struct IMPBFFEXPORT AVMeasurement {
    //! Names the volumes were registered under.
    std::string av1_name, av2_name;
    //! Experimental distance, Angstrom.
    double distance;
    //! Asymmetric error, Angstrom; `error_neg` applies when model < experiment.
    double error_neg, error_pos;
    //! \f$R_0\f$, used by `"RDAMeanE"`.
    double forster_radius;
    //! `"RDAMean"` (mean inter-point distance), `"RDAMeanE"` (FRET-averaged)
    //! or `"Rmp"` (between the clouds' mean positions).
    std::string distance_type;

    AVMeasurement(std::string av1_name = "", std::string av2_name = "",
                  double distance = 0.0, double error_neg = 3.0,
                  double error_pos = 5.0, double forster_radius = 52.0,
                  std::string distance_type = "RDAMean");

    IMP_SHOWABLE_INLINE(AVMeasurement,
                        out << "AVMeasurement(" << av1_name << "-" << av2_name
                            << ", " << distance << " A)");
};

//! A labelling site with an experimental distance to a partner.
/*!
    Paired distances are encoded by adding two sites, one per residue, and
    letting the restraint score every unique pair. Use #IMP::bff::AVMeasurement
    when a particular pair has to be left out.
*/
struct IMPBFFEXPORT ProbeSite {
    //! Row index into the restraint's coordinate array.
    int residue_seq_number;
    //! Attachment atom, for callers that resolve sites by name.
    std::string atom_name;
    //! Experimental distance to a partner, Angstrom.
    double distance;
    //! Asymmetric error, Angstrom.
    double error_neg, error_pos;
    //! \f$R_0\f$; unused by direct labelling, carried so a site can move
    //! between the two restraints unchanged.
    double forster_radius;

    ProbeSite(int residue_seq_number = 0, std::string atom_name = "CB",
                 double distance = 0.0, double error_neg = 3.0,
                 double error_pos = 5.0, double forster_radius = 52.0);

    IMP_SHOWABLE_INLINE(ProbeSite,
                        out << "ProbeSite(" << residue_seq_number << " "
                            << atom_name << ")");
};

IMP_VALUES(AVMeasurement, AVMeasurements);
IMP_VALUES(ProbeSite, ProbeSites);

//! Chi-squared over a network of accessible-volume pair distances.
class IMPBFFEXPORT SimpleProbeNetworkRestraint {
    std::string name_;
    double weight_;
    std::map<std::string, AccessibleVolume> avs_;
    std::vector<AVMeasurement> measurements_;

    double model_distance(const AVMeasurement& m) const;

public:
    SimpleProbeNetworkRestraint(std::string name = "SimpleProbeNetworkRestraint")
        : name_(name), weight_(1.0) {}

    //! Register a pre-computed accessible volume under a name.
    void add_av(std::string name, const AccessibleVolume& av) { avs_[name] = av; }

    void add_measurement(const AVMeasurement& measurement) {
        measurements_.push_back(measurement);
    }

    //! Weighted sum of the per-measurement chi-squared contributions.
    double evaluate() const;

    //! Model distance per measurement, keyed `"<av1>_<av2>"`.
    std::map<std::string, double> get_model_distances() const;

    std::string get_name() const { return name_; }
    double get_weight() const { return weight_; }
    void set_weight(double w) { weight_ = w; }
    unsigned int get_n_avs() const { return static_cast<unsigned int>(avs_.size()); }
    unsigned int get_n_measurements() const {
        return static_cast<unsigned int>(measurements_.size());
    }
    bool has_av(std::string name) const { return avs_.count(name) > 0; }

    IMP_SHOWABLE_INLINE(SimpleProbeNetworkRestraint,
                        out << "SimpleProbeNetworkRestraint(name=" << name_
                            << ", n_avs=" << avs_.size() << ", n_measurements="
                            << measurements_.size() << ")");
};

//! Chi-squared over attachment-atom distances, with no accessible volume.
/*!
    The model distance for a pair of sites is the Euclidean distance between the
    two rows of the coordinate array their `residue_seq_number`s index. Every
    unique pair of added sites contributes.
*/
class IMPBFFEXPORT DirectProbeRestraint {
    std::vector<double> xyz_;
    int n_atoms_;
    double weight_;
    std::vector<ProbeSite> sites_;

public:
    //! \param[in] xyz,n_atoms,n_dim coordinates, `(N, 3)`
    //! \param[in] weight overall restraint weight
    DirectProbeRestraint(double* xyz = 0, int n_atoms = 0, int n_dim = 0,
                            double weight = 1.0);

    void add_site(const ProbeSite& site) { sites_.push_back(site); }

    //! Weighted sum of chi-squared over every unique pair of sites.
    /*! Zero when there are fewer than two sites, or no coordinates. */
    double evaluate() const;

    double get_weight() const { return weight_; }
    void set_weight(double w) { weight_ = w; }
    unsigned int get_n_sites() const {
        return static_cast<unsigned int>(sites_.size());
    }

    IMP_SHOWABLE_INLINE(DirectProbeRestraint,
                        out << "DirectProbeRestraint(n_sites="
                            << sites_.size() << ")");
};

IMP_VALUES(SimpleProbeNetworkRestraint, SimpleProbeNetworkRestraints);
IMP_VALUES(DirectProbeRestraint, DirectProbeRestraints);

IMPBFF_END_NAMESPACE

#endif //IMPBFF_PROBERESTRAINTS_H
