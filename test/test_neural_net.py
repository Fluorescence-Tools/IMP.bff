"""A dense network, evaluated in batches, on whatever backend is loaded.

`NeuralNet` wraps the vendored `internal/MlpCore.h` kernels and offers the
whole forward pass to the compute door. The tests that matter are therefore
about *agreement*: the accelerator may take a batch or decline it, and the
answer must be the same either way to the precision it works in.

`get_last_backend()` reports which one ran, so an accelerated result cannot be
mistaken for a CPU one, and the GPU tests below skip where nothing is loaded.
"""

import json

import numpy as np
import pytest

import IMP.bff


def _model(widths, activations, seed=0):
    rng = np.random.default_rng(seed)
    layers = []
    for k in range(len(widths) - 1):
        n_in, n_out = widths[k], widths[k + 1]
        layers.append({
            "n_in": n_in, "n_out": n_out,
            # row-major n_out x n_in, as DenseLayer holds it
            "weight": rng.normal(0.0, 1.0 / np.sqrt(n_in), n_in * n_out).tolist(),
            "bias": rng.normal(0.0, 0.1, n_out).tolist(),
            "activation": activations[k],
        })
    return json.dumps({"format": "tttrlib.neural_net", "layers": layers})


def _numpy_forward(spec, X):
    """The same network in numpy, so the C++ is checked against something that
    shares no code with it."""
    m = json.loads(spec)
    A = np.asarray(X, dtype=float)
    for layer in m["layers"]:
        W = np.asarray(layer["weight"]).reshape(layer["n_out"], layer["n_in"])
        b = np.asarray(layer["bias"])
        Z = A @ W.T + b
        act = layer["activation"]
        if act == "identity": A = Z
        elif act == "relu": A = np.maximum(Z, 0.0)
        elif act == "tanh": A = np.tanh(Z)
        elif act == "logistic": A = 1.0 / (1.0 + np.exp(-Z))
        elif act == "softplus": A = np.logaddexp(0.0, Z)
        elif act == "silu": A = Z / (1.0 + np.exp(-Z))
        elif act == "sin": A = np.sin(Z)
        else: raise AssertionError(act)
    return A


def test_it_reads_the_shape_off_the_document():
    net = IMP.bff.NeuralNet(_model([3, 8, 8, 2], ["tanh", "relu", "identity"]))
    assert net.get_n_inputs() == 3
    assert net.get_n_outputs() == 2
    assert net.get_n_layers() == 3


def test_it_refuses_a_document_that_is_not_one():
    with pytest.raises(ValueError):
        IMP.bff.NeuralNet("{not json")
    with pytest.raises(ValueError):
        IMP.bff.NeuralNet(json.dumps({"format": "tttrlib.neural_net", "layers": []}))
    # layers that do not chain
    bad = json.loads(_model([3, 8, 2], ["relu", "identity"]))
    bad["layers"][1]["n_in"] = 7
    with pytest.raises(ValueError):
        IMP.bff.NeuralNet(json.dumps(bad))


@pytest.mark.parametrize("activation", ["identity", "relu", "tanh", "logistic",
                                        "softplus", "silu", "sin"])
def test_every_activation_matches_numpy(activation):
    """The GPU writes each of these out again in WGSL, so each is checked."""
    spec = _model([3, 16, 2], [activation, "identity"])
    net = IMP.bff.NeuralNet(spec)
    X = np.random.default_rng(4).normal(0.0, 1.5, (40, 3))
    y = np.asarray(net.predict(X.ravel().tolist(), 40)).reshape(40, 2)
    np.testing.assert_allclose(y, _numpy_forward(spec, X), rtol=1e-9, atol=1e-9)


def test_a_batch_of_no_rows_is_empty_not_an_error():
    net = IMP.bff.NeuralNet(_model([3, 8, 2], ["relu", "identity"]))
    assert len(np.asarray(net.predict([], 0))) == 0


def test_the_batch_length_is_checked():
    net = IMP.bff.NeuralNet(_model([3, 8, 2], ["relu", "identity"]))
    with pytest.raises(ValueError):
        net.predict([1.0, 2.0], 1)          # 2 values for 3 inputs


@pytest.mark.skipif(IMP.bff.get_compute_backend_name() == "cpu",
                    reason="no GPU backend loaded")
def test_the_accelerator_and_the_cpu_agree():
    """Large enough that the backend takes it: the two must agree to the
    single precision the GPU works in."""
    spec = _model([4, 64, 64, 2], ["tanh", "relu", "identity"])
    net = IMP.bff.NeuralNet(spec)
    n_rows = 30000
    X = np.random.default_rng(7).normal(0.0, 1.0, (n_rows, 4))
    xl = X.ravel().tolist()
    y_gpu = np.asarray(net.predict(xl, n_rows))
    assert net.get_last_backend() != "cpu", "the batch should have been taken"
    IMP.bff.reset_compute_backend()
    try:
        y_cpu = np.asarray(net.predict(xl, n_rows))
        assert net.get_last_backend() == "cpu"
    finally:
        IMP.bff.enable_gpu(quiet=True)
    assert np.abs(y_gpu - y_cpu).max() / np.abs(y_cpu).max() < 1e-5


@pytest.mark.skipif(IMP.bff.get_compute_backend_name() == "cpu",
                    reason="no GPU backend loaded")
def test_a_batch_too_small_to_be_worth_the_trip_stays_on_the_cpu():
    net = IMP.bff.NeuralNet(_model([3, 8, 2], ["relu", "identity"]))
    X = np.random.default_rng(8).normal(0.0, 1.0, (10, 3))
    net.predict(X.ravel().tolist(), 10)
    assert net.get_last_backend() == "cpu"
