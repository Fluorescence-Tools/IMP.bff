"""The chinet-shaped keyword constructor on bff.Port.

chisurf's Parameter builds its ports through chinet's Python constructor
surface -- ``Port(value=..., name=..., lb=..., ub=..., is_bounded=...)`` --
so phase 3 of removing chinet needs the same call to land on
``IMP.bff.Port`` unchanged. SWIG cannot generate keyword arguments for an
overloaded constructor, so the keyword form is a shim over the positional
overloads; these tests pin the shim to chinet's accepted inputs and
semantics: python and numpy scalars, lists, tuples and arrays, the flag
kwargs, the prior dict, and the ``link`` property being assignable.
"""

import numpy as np
import pytest

import IMP.bff as bff


def test_kwargs_ctor_float_scalar():
    p = bff.Port(value=3.5, name="amp")
    assert p.value == 3.5
    assert p.name == "amp"
    assert p.fixed is False
    assert p.get_value_type() == 1  # a float argument promotes, as chinet


def test_kwargs_ctor_int_scalar_stays_int():
    p = bff.Port(value=5)
    assert p.value == 5
    assert p.get_value_type() == 0


def test_kwargs_ctor_numpy_scalars():
    p = bff.Port(value=np.float64(2.5))
    assert p.value == 2.5
    assert p.get_value_type() == 1

    q = bff.Port(value=np.int64(7))
    assert q.value == 7
    assert q.get_value_type() == 0


def test_kwargs_ctor_list_tuple_and_arrays():
    for v in ([1.5, 2.5], (1.5, 2.5), np.array([1.5, 2.5]),
              np.array([1.5, 2.5], dtype=np.float64)):
        p = bff.Port(value=v)
        assert p.get_is_vector() is True
        assert np.allclose(p.value, [1.5, 2.5])


def test_kwargs_ctor_chisurf_shape():
    """Exactly the call chisurf's Parameter makes."""
    p = bff.Port(value=np.atleast_1d(3.5).astype(np.float64),
                 name="gamma", lb=0.0, ub=10.0, is_bounded=True)
    assert p.get_is_vector() is True
    assert p.get_value_type() == 1  # chinet's atleast_1d port is code 1
    assert p.name == "gamma"
    assert p.bounded is True
    assert p.bounds == (0.0, 10.0)
    assert np.allclose(p.value, [3.5])


def test_kwargs_ctor_bounds_clip_initial_value():
    p = bff.Port(value=np.array([50.0]), lb=0.0, ub=10.0, is_bounded=True)
    assert np.allclose(p.value, [10.0])


def test_kwargs_ctor_flags():
    p = bff.Port(value=1.0, fixed=True, is_output=True, is_reactive=True)
    assert p.fixed is True
    assert p.is_output is True
    assert p.is_reactive is True
    p.value = 5.0  # fixed: a no-op
    assert p.value == 1.0


def test_kwargs_ctor_prior_dict():
    p = bff.Port(value=1.0, prior={"kind": "normal", "mu": 1.0,
                                   "sigma": 0.5})
    assert p.prior == {"kind": "normal", "mu": 1.0, "sigma": 0.5}
    p.prior = None
    assert p.prior is None


def test_kwargs_ctor_int_array_is_int_vector():
    p = bff.Port(value=np.array([1, 2, 3], dtype=np.int64))
    assert p.get_value_type() == 2
    assert np.allclose(p.value, [1, 2, 3])


def test_kwargs_ctor_value_type_override():
    p = bff.Port(value=[1, 2, 3], value_type=3)
    assert p.get_value_type() == 3


def test_kwargs_ctor_unknown_kwarg_raises():
    with pytest.raises(TypeError):
        bff.Port(value=1.0, nonsense=True)


def test_kwargs_ctor_positional_still_works():
    p = bff.Port(23.0, True)
    assert p.value == 23.0
    assert p.fixed is True
    q = bff.Port()
    assert q.value == 0


def test_link_property_is_assignable():
    master = bff.Port(value=7.0)
    follower = bff.Port(value=0.0)
    follower.link = master
    assert follower.is_linked()
    assert follower.value == 7.0
    follower.link = None
    assert not follower.is_linked()


def test_get_json_read_json_round_trip():
    p = bff.Port(value=np.array([2.0]), name="x", lb=-1.0, ub=5.0,
                 is_bounded=True, fixed=True,
                 prior={"kind": "lognormal", "mu": 0.0, "sigma": 0.4})
    doc = p.get_json()
    q = bff.Port(value=0.0)
    assert q.read_json(doc) is True
    assert q.name == "x"
    assert q.fixed is True
    assert q.bounded is True
    assert q.bounds == (-1.0, 5.0)
    assert np.allclose(q.value, [2.0])
    assert q.prior == {"kind": "lognormal", "mu": 0.0, "sigma": 0.4}
    assert q.get_value_type() == p.get_value_type()


def test_get_json_read_json_int_port():
    p = bff.Port(value=3)
    q = bff.Port(value=0.0)
    q.read_json(p.get_json())
    assert q.value == 3
    assert q.get_value_type() == 0  # the document's code, not the data's


def test_read_json_restores_value_of_fixed_port():
    """chinet's set_document wrote the data directly, fixed or not."""
    p = bff.Port(value=9.0, fixed=True)
    q = bff.Port(value=0.0, fixed=True)
    q.read_json(p.get_json())
    assert q.value == 9.0
    assert q.fixed is True
