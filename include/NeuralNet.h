/**
 *  \file IMP/bff/NeuralNet.h
 *  \brief Evaluating a dense network, on the CPU or through the compute door.
 *
 * The kernels are `internal/MlpCore.h`, vendored verbatim from tttrlib so a
 * network trained there evaluates identically here; this header is the face
 * that bff and its bindings use. What it adds is where the arithmetic runs:
 * `predict()` offers the whole forward pass to whatever
 * `IMP/bff/Compute.h` has loaded, and runs it on the CPU when nothing is
 * loaded, when the backend declines, or when the batch is too small to be
 * worth the trip.
 *
 * \authors Thomas-Otavio Peulen
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */

#ifndef IMPBFF_NEURALNET_H
#define IMPBFF_NEURALNET_H

#include <IMP/bff/bff_config.h>
#include <IMP/bff/Base.h>

#include <memory>
#include <string>
#include <vector>

IMPBFF_BEGIN_NAMESPACE

//! A dense multilayer perceptron, evaluated in batches.
/*!
    Constructed from a `tttrlib.neural_net` JSON document -- the same one
    `NeuralNet::train` writes there and the same one scikit-learn's
    `MLPRegressor` converts to -- so the weights need no bff-specific format.

    **The batch is the unit.** One `predict()` call crosses the binding once
    and, if an accelerator is loaded, crosses the plugin boundary once, however
    deep the network. Calling it per row would repeat the mistake that once
    made the quenching suite forty times slower.
*/
class IMPBFFEXPORT NeuralNet {
public:
    //! Build from a `tttrlib.neural_net` JSON document.
    /*! \throws IMP::ValueException if the document is not one, or if the
        layers do not chain. */
    explicit NeuralNet(const std::string& json);
    ~NeuralNet();

    //! Inputs the network takes.
    int get_n_inputs() const;
    //! Outputs it produces.
    int get_n_outputs() const;
    //! Layers it has.
    int get_n_layers() const;

    //! Where the last predict() ran: `cpu`, or the backend's name.
    /*! A loaded accelerator may still decline a batch as too small, so this
        is a fact about the last call rather than about the machine. */
    std::string get_last_backend() const;

    //! Evaluate a batch.
    /*!
        \param[in] x the batch, row-major, `n_rows * get_n_inputs()` long
        \param[in] n_rows rows in the batch
        \param[out] out_view,n_out_view `n_rows * get_n_outputs()`, a managed
                    numpy view; in C++ free it with `std::free`
    */
    void predict(const std::vector<double>& x, int n_rows,
                 double** out_view, int* n_out_view) const;

private:
    struct Impl;
    std::shared_ptr<Impl> impl_;
};

IMPBFF_END_NAMESPACE

#endif  // IMPBFF_NEURALNET_H
