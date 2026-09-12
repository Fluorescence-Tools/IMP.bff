"""The reconvolved basis, so a caller can differentiate the model.

`TcspcDecay` produces the finished curve, which is the sum over species of
`amplitude * (lifetime (*) response)`. A caller that needs a Jacobian wants
the terms rather than the sum, because **the model is linear in the
amplitudes given the basis** -- so the derivative with respect to all of them
is that matrix, and a Jacobian is one matrix product instead of one finite
difference per parameter. The consumer that asked for this
(`../ucfret`) has 113 parameters and a Fisher-scoring loop; without the
basis it would evaluate the curve 113 times per iteration to learn what one
matrix already says.

The contract these pin: row-major bins x species, columns matching the input
spectrum **one for one**, with the amplitude compaction disabled while the
basis is requested. A compacted basis plus an index map is right the day it
is written and drifts afterwards, and an amplitude indexed against the wrong
column fails silently.
"""

import numpy as np
import pytest

import IMP.bff


def _decay(amplitudes, lifetimes, n=64, emit=True):
    d = IMP.bff.TCSPCDecay("decay")
    d.set_number_of_lifetimes(len(lifetimes))
    # the curve goes to the port keyed by the node's own name
    d.add_output_port("decay", IMP.bff.GraphPort([0.0]))
    d.set_emit_basis(emit)
    # a narrow response, one bin wide, so the convolution is easy to reason about
    irf = np.zeros(n)
    irf[2] = 1.0
    d.set_response(irf.tolist())
    d.set_data(np.zeros(n).tolist(), np.ones(n).tolist())
    for i, (a, t) in enumerate(zip(amplitudes, lifetimes)):
        d.get_port("a%d" % i).set_value(float(a))
        d.get_port("t%d" % i).set_value(float(t))
    d.get_port("scatter").set_value(0.0)
    d.get_port("background").set_value(0.0)
    d.set_autoscale(False)
    d.update()
    return d


def _basis(d, n_species):
    v = np.asarray(d.get_output_port(IMP.bff.TCSPCDecay.basis_port_key())
                   .get_value_view(), dtype=float)
    return v.reshape(-1, n_species)


def test_the_basis_has_a_column_per_species():
    lifetimes = [0.5, 2.0, 4.0]
    d = _decay([0.2, 0.3, 0.5], lifetimes)
    basis = _basis(d, len(lifetimes))
    curve = np.asarray(d.get_output_port("decay").get_value_view(), dtype=float)
    assert basis.shape == (len(curve), len(lifetimes))


def test_the_curve_is_the_basis_times_the_amplitudes():
    """The property the whole thing rests on. If this holds, a Jacobian with
    respect to the amplitudes is the basis itself."""
    amps = [0.2, 0.3, 0.5]
    lifetimes = [0.5, 2.0, 4.0]
    d = _decay(amps, lifetimes)
    basis = _basis(d, len(lifetimes))
    curve = np.asarray(d.get_output_port("decay").get_value_view(), dtype=float)
    np.testing.assert_allclose(basis @ np.asarray(amps), curve, rtol=1e-12, atol=1e-12)


def test_the_columns_survive_an_amplitude_of_zero():
    """Compaction would drop a zero-amplitude species and shift every column
    after it. The caller's amplitude vector has a fixed layout, so the columns
    must not move."""
    lifetimes = [0.5, 2.0, 4.0]
    d = _decay([0.5, 0.0, 0.5], lifetimes)
    basis = _basis(d, len(lifetimes))
    assert basis.shape[1] == 3, "the zero-amplitude species keeps its column"
    assert np.abs(basis[:, 1]).max() > 0.0, "and the column is its response, not zeros"


def test_a_column_is_that_species_alone():
    """Column i must be the response of lifetime i at unit amplitude,
    independent of what the amplitudes are."""
    lifetimes = [0.5, 2.0, 4.0]
    a = _basis(_decay([0.2, 0.3, 0.5], lifetimes), 3)
    b = _basis(_decay([9.0, -1.0, 0.25], lifetimes), 3)
    np.testing.assert_allclose(a, b, rtol=1e-12, atol=1e-12)
    one = _basis(_decay([1.0], [2.0]), 1)
    np.testing.assert_allclose(a[:, 1], one[:, 0], rtol=1e-12, atol=1e-12)


def test_the_curve_is_unchanged_by_asking_for_the_basis():
    """Turning it on must not move the answer -- only disable an optimisation
    that was bit-exact at the default threshold anyway."""
    amps, lifetimes = [0.2, 0.3, 0.5], [0.5, 2.0, 4.0]
    with_basis = np.asarray(_decay(amps, lifetimes, emit=True)
                            .get_output_port("decay").get_value_view())
    without = np.asarray(_decay(amps, lifetimes, emit=False)
                         .get_output_port("decay").get_value_view())
    np.testing.assert_allclose(with_basis, without, rtol=1e-12, atol=1e-12)


def test_the_port_is_absent_until_it_is_asked_for():
    d = IMP.bff.TCSPCDecay("decay")
    d.set_number_of_lifetimes(2)
    d.add_output_port("decay", IMP.bff.GraphPort([0.0]))
    assert not d.get_emit_basis()
    assert d.get_output_port(IMP.bff.TCSPCDecay.basis_port_key()) is None
    d.set_emit_basis(True)
    assert d.get_emit_basis()
    assert d.get_output_port(IMP.bff.TCSPCDecay.basis_port_key()) is not None
