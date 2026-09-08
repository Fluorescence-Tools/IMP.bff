/**
 * \file NeuralNet.cpp
 * \brief Evaluating a dense network, on the CPU or through the compute door.
 *
 * Copyright 2007-2026 IMP Inventors. All rights reserved.
 */
#include <IMP/bff/NeuralNet.h>

#include <IMP/bff/Compute.h>
#include <IMP/bff/internal/OutputView.h>

#define TTTRLIB_MLPCORE_NAMESPACE IMP::bff::internal
#include <IMP/bff/internal/MlpCore.h>
#include <IMP/bff/internal/json.h>

#include <cstdlib>

IMPBFF_BEGIN_NAMESPACE

namespace mc = IMP::bff::internal::mlpcore;

struct NeuralNet::Impl {
    internal::MlpModel model;
    // Flattened once at construction: the door takes plain arrays, and a
    // network is evaluated far more often than it is loaded.
    std::vector<int> n_in, n_out, activation;
    std::vector<double> weights, biases;
    mutable std::string last_backend;
};

NeuralNet::NeuralNet(const std::string& json) : impl_(new Impl) {
    nlohmann::json j = nlohmann::json::parse(json, nullptr, false);
    if (j.is_discarded()) IMP_THROW("NeuralNet: the model is not JSON", IMP::ValueException);
    try {
        impl_->model = mc::model_from_json(j);
        impl_->model.validate();
    } catch (const std::exception& e) {
        IMP_THROW(e.what(), IMP::ValueException);
    }
    for (std::size_t i = 0; i < impl_->model.layers.size(); ++i) {
        const internal::DenseLayer& l = impl_->model.layers[i];
        impl_->n_in.push_back(l.n_in);
        impl_->n_out.push_back(l.n_out);
        impl_->activation.push_back(static_cast<int>(l.activation));
        impl_->weights.insert(impl_->weights.end(), l.weight.begin(), l.weight.end());
        impl_->biases.insert(impl_->biases.end(), l.bias.begin(), l.bias.end());
    }
    impl_->last_backend = "cpu";
}

NeuralNet::~NeuralNet() {}

int NeuralNet::get_n_inputs() const { return impl_->model.n_inputs(); }
int NeuralNet::get_n_outputs() const { return impl_->model.n_outputs(); }
int NeuralNet::get_n_layers() const { return static_cast<int>(impl_->model.layers.size()); }
std::string NeuralNet::get_last_backend() const { return impl_->last_backend; }

void NeuralNet::predict(const std::vector<double>& x, int n_rows,
                        double** out_view, int* n_out_view) const {
    const int n_in = get_n_inputs(), n_out = get_n_outputs();
    if (n_rows < 0) IMP_THROW("NeuralNet::predict: n_rows must be >= 0", IMP::ValueException);
    if (x.size() != static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_in))
        IMP_THROW("NeuralNet::predict: x must be n_rows * n_inputs long",
                  IMP::ValueException);
    std::vector<double> y;
    impl_->last_backend = "cpu";

    // The accelerator is offered the whole forward pass and may decline it --
    // no device, or a batch small enough that the round trip costs more than
    // the arithmetic. Nothing is written when it declines.
    const ImpBffComputeBackend* backend = get_compute_backend();
    bool done = false;
    if (n_rows > 0 && backend && backend->abi >= 2 && backend->mlp_forward) {
        std::vector<double> xs(x);
        mc::detail::scale_in(xs, n_rows, n_in, impl_->model.x_scaler);
        y.assign(static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_out), 0.0);
        const int rc = backend->mlp_forward(
                get_n_layers(), impl_->n_in.data(), impl_->n_out.data(),
                impl_->activation.data(), impl_->weights.data(),
                impl_->biases.data(), xs.data(), n_rows, y.data());
        if (rc == 0) {
            mc::detail::unscale_out(y, n_rows, n_out, impl_->model.y_scaler);
            impl_->last_backend = get_compute_backend_name();
            done = true;
        } else {
            y.clear();
        }
    }
    if (!done) {
        std::vector<double> d1, d2;
        mc::model_predict<mc::PortableGemm>(impl_->model, x.data(), n_rows,
                                                  0, nullptr, y, d1, d2);
    }
    internal::copy_to_view(y, out_view, n_out_view);
}

IMPBFF_END_NAMESPACE
