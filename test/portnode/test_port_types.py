"""Port value types: int, float, bool, and each of them as a vector.

The type code existed before this file (0 int, 1 float, 2 int vector, 3 float
vector) but it was a bare number with no bool in it, and a port answered every
read as a double -- `Port(value=True).value` was `1.0`, which is only half a
type. This pins the six types, the coercion rules, and the two traps found
while adding them.

The codes are **not** free to change: `Session` writes them into the chinet
document and `test/session/chinet_fixture.jsonl` pins 0-3. Bool took fresh
codes rather than renumbering, which is what `test_the_chinet_codes_are_fixed`
is here to keep true.
"""

import numpy as np
import pytest

import IMP.bff as bff
from IMP.bff import (
    PORT_BOOL, PORT_BOOL_VECTOR, PORT_FLOAT, PORT_FLOAT_VECTOR, PORT_INT,
    PORT_INT_VECTOR, Port, port_value_type_element, port_value_type_is_vector,
    port_value_type_name, port_value_type_of,
)


def test_the_chinet_codes_are_fixed():
    """0-3 are chinet's and are written into saved documents verbatim."""
    assert (PORT_INT, PORT_FLOAT, PORT_INT_VECTOR, PORT_FLOAT_VECTOR) == (0, 1, 2, 3)
    assert (PORT_BOOL, PORT_BOOL_VECTOR) == (4, 5)


@pytest.mark.parametrize("code,element,is_vector,name", [
    (PORT_INT, PORT_INT, False, "int"),
    (PORT_FLOAT, PORT_FLOAT, False, "float"),
    (PORT_INT_VECTOR, PORT_INT, True, "int[]"),
    (PORT_FLOAT_VECTOR, PORT_FLOAT, True, "float[]"),
    (PORT_BOOL, PORT_BOOL, False, "bool"),
    (PORT_BOOL_VECTOR, PORT_BOOL, True, "bool[]"),
])
def test_the_code_helpers_decompose_every_type(code, element, is_vector, name):
    assert port_value_type_element(code) == element
    assert port_value_type_is_vector(code) is is_vector
    assert port_value_type_name(code) == name
    assert port_value_type_of(element, is_vector) == code


# --- scalars read back in their own type ----------------------------------


def test_a_scalar_reads_back_as_its_own_type():
    assert Port(value=12).value == 12
    assert isinstance(Port(value=12).value, int)
    assert isinstance(Port(value=1.5).value, float)
    assert Port(value=True).value is True
    assert Port(value=False).value is False


def test_a_bool_port_is_built_from_a_python_bool_not_an_int():
    """`bool` is a subclass of `int`, so the dispatch has to test it first.

    Before this it did not, and every flag became an integer port on the way
    in -- the reason `_set_value` checks `isinstance(v, bool)` before it looks
    for `__index__`.
    """
    assert Port(value=True).get_value_type() == PORT_BOOL
    assert Port(value=1).get_value_type() == PORT_INT
    assert Port(value=np.bool_(True)).get_value_type() == PORT_BOOL


# --- bool is declared, not inferred ---------------------------------------


def test_a_bool_port_stays_bool_when_written_a_float():
    """Bool was the first declared type; every type is one now.

    This test predates that: bool was carved out as an exception because
    promoting a flag on assignment obviously destroys what the port means. The
    same argument applies to int, and since 2026-09-02 it is the general rule
    rather than a special case for bool.
    """
    p = Port(value=True)
    p.value = 3.7
    assert p.value is True
    assert p.get_value_type() == PORT_BOOL
    p.value = 0
    assert p.value is False
    assert p.get_value_type() == PORT_BOOL


def test_an_int_port_accepts_a_float_and_coerces_it():
    """A write of another kind is accepted, not refused -- and not promoted.

    chinet inferred the dtype from every write, so 1.5 into an integer port
    made it a float port. That is what a declared type is for: whatever held a
    reference to this port had already decided what it was, and a write is not
    entitled to change that answer underneath it.
    """
    p = Port(value=12)
    p.value = 1.5
    assert p.get_value_type() == PORT_INT
    assert p.value == 1


def test_only_set_value_type_takes_a_port_out_of_bool():
    p = Port(value=True)
    p.value = 5
    assert p.value is True
    p.set_value_type(PORT_INT)
    assert p.get_value_type() == PORT_INT
    assert p.value == 1          # the stored flag, now read as an integer


def test_set_value_type_converts_the_stored_data():
    p = Port(value=[2.7, -1.2, 0.0])
    p.set_value_type(PORT_INT_VECTOR)
    assert list(p.value) == [2, -1, 0]          # truncation, as astype(int64)
    p.set_value_type(PORT_BOOL_VECTOR)
    assert list(p.value) == [True, True, False]  # != 0


# --- vectors ---------------------------------------------------------------


def test_a_float_vector_is_not_truncated_by_the_int_path():
    """The regression that removing the `vector<int>` overload fixed.

    A Python list of floats converts to `std::vector<int>` as happily as to
    `std::vector<double>`, so an overload let SWIG's dispatch pick the element
    type -- and it picked wrong: `Port(value=[1.5, 2.5])` came back `[1, 2]`.
    """
    for value in ([1.5, 2.5], (1.5, 2.5), np.array([1.5, 2.5])):
        p = Port(value=value)
        assert p.get_is_vector() is True
        np.testing.assert_allclose(p.value, [1.5, 2.5])


def test_an_integer_vector_keeps_its_type():
    """Closes the divergence the Port.h file comment used to record.

    chinet kept int64 for a list of ints; every vector written through the
    double entry point was float input, so the port float-typed. The named
    integer entry point is what makes this possible.
    """
    p = Port(value=[1, 2, 3])
    assert p.get_value_type() == PORT_INT_VECTOR
    assert p.value.dtype == np.int64
    assert list(p.value) == [1, 2, 3]


def test_set_value_vector_int_is_named_not_overloaded():
    p = Port(value=0)                       # an int port, not a float one
    p.set_value_vector_int([4, 5, 6])
    assert p.get_value_type() == PORT_INT_VECTOR
    assert list(p.value) == [4, 5, 6]


def test_an_integer_write_into_a_float_port_stays_float():
    """The entry point says what the *input* is, not what the port becomes.

    `set_value_vector_int` on a float port leaves it float and stores the
    values as doubles. Only `set_value_type` changes a declared type.
    """
    p = Port(value=0.0)                     # float port
    p.set_value_vector_int([4, 5, 6])
    assert p.get_value_type() == PORT_FLOAT_VECTOR
    p.set_value_type(PORT_INT_VECTOR)
    assert p.get_value_type() == PORT_INT_VECTOR


def test_a_bool_vector_coerces_every_element():
    p = Port(value=True)
    p.value = [1, 0, 2, -3]
    assert p.get_value_type() == PORT_BOOL_VECTOR
    assert p.value.dtype == bool
    assert list(p.value) == [True, False, True, True]


def test_a_vector_port_is_a_port_input_and_an_output():
    """Vectors have to work in both directions, through a link."""
    source = Port(value=[1.0, 2.0, 3.0], is_output=True)
    sink = Port(value=0.0)
    sink.set_link(source)
    np.testing.assert_allclose(sink.value, [1.0, 2.0, 3.0])
    source.value = [4.0, 5.0]
    np.testing.assert_allclose(sink.value, [4.0, 5.0])


# --- an integer port holds an integer -------------------------------------
#
# The type code alone was only half a type: storage was `double` for every
# port, so an integer port was a double wearing a label. Two things went wrong
# and both are pinned here. `Port(value=2**53+1)` came back **typed float**,
# because a Python int wider than 32 bits does not convert to C++ `int` and
# SWIG fell through to the `double` overload; and even a correctly typed int
# port rounded, because the value had passed through a double to get in.

BIG = [2 ** 53 + 1, 2 ** 60 + 7, -(2 ** 61) - 3, 9223372036854775807]


@pytest.mark.parametrize("value", BIG)
def test_a_scalar_integer_survives_past_two_to_the_53(value):
    p = Port(value=value)
    assert p.get_value_type() == PORT_INT, "a wide int must not land on double"
    assert p.value == value


@pytest.mark.parametrize("value", BIG)
def test_writing_a_wide_integer_after_construction_is_exact(value):
    p = Port(value=0)
    p.value = value
    assert p.get_value_type() == PORT_INT
    assert p.value == value


def test_an_integer_vector_survives_past_two_to_the_53():
    p = Port(value=BIG)
    assert p.get_value_type() == PORT_INT_VECTOR
    assert [int(x) for x in p.value] == BIG
    p.value = [2 ** 62, 5]
    assert [int(x) for x in p.value] == [2 ** 62, 5]


def test_a_wide_integer_survives_a_link():
    """Links copy through the port API, so exactness has to survive that too."""
    source = Port(value=2 ** 60 + 7, is_output=True)
    sink = Port(value=0)
    sink.set_link(source)
    assert sink.value == 2 ** 60 + 7


def test_a_wide_integer_survives_the_document():
    p = Port(value=2 ** 60 + 7)
    restored = Port()
    restored.read_json(p.get_json())
    assert restored.get_value_type() == PORT_INT
    assert restored.value == 2 ** 60 + 7


def test_a_wide_integer_port_stays_exact_across_a_float_write():
    """The payoff of an immutable type, stated as exactness.

    Under the old inferred promotion, one stray float write turned an integer
    port into a float port and every later value was rounded. Now the port is
    still an integer port afterwards, so the next wide value is still exact.
    """
    p = Port(value=2 ** 60 + 7)
    p.value = 1.5
    assert p.get_value_type() == PORT_INT
    assert p.value == 1
    p.value = 2 ** 62 + 5
    assert p.value == 2 ** 62 + 5


def test_a_numpy_bool_array_makes_a_bool_port_not_an_int_one():
    p = Port(value=np.array([True, False, True]))
    assert p.get_value_type() == PORT_BOOL_VECTOR
    assert p.value.dtype == bool
    assert list(p.value) == [True, False, True]


def test_a_wide_integer_survives_a_saved_session(tmp_path):
    """The C++ writer and reader, not just Port.get_json.

    `Session` emitted through `get_value()` -- a double -- and read back
    through `get<double>()`, so a saved integer port was rounded twice on a
    path that never touches Python. Both ends take the exact route now.
    """
    node = bff.Node("n")
    node.add_input_port("count", Port(value=2 ** 60 + 7))
    node.add_input_port("flag", Port(value=True))
    node.add_input_port("ids", Port(value=[2 ** 53 + 1, 3]))
    session = bff.Session()
    session.add_node("n", node)
    path = str(tmp_path / "s.jsonl")
    session.save(path)

    restored = bff.Session.load(path).get_node("n")
    assert restored.get_input_port("count").value == 2 ** 60 + 7
    assert restored.get_input_port("flag").value is True
    assert [int(x) for x in restored.get_input_port("ids").value] == [2 ** 53 + 1, 3]


# --- the declared type is immutable ---------------------------------------


@pytest.mark.parametrize("declared,written,expected_code,expected", [
    (PORT_INT, 1.5, PORT_INT, 1),
    (PORT_INT, True, PORT_INT, 1),
    (PORT_FLOAT, 3, PORT_FLOAT, 3.0),
    (PORT_FLOAT, True, PORT_FLOAT, 1.0),
    (PORT_BOOL, 3.7, PORT_BOOL, True),
    (PORT_BOOL, 0, PORT_BOOL, False),
])
def test_a_write_of_any_kind_is_accepted_and_never_retypes(
        declared, written, expected_code, expected):
    """Accept every kind, coerce to the declared one, change type for none.

    The whole matrix, because the rule is only worth anything if it has no
    exceptions -- bool used to be the one type that behaved this way and the
    other two inferred themselves from the data.
    """
    p = Port()
    p.set_value_type(declared)
    p.value = written
    assert p.get_value_type() == expected_code
    assert p.value == expected
    assert type(p.value) is type(expected)


@pytest.mark.parametrize("declared,vector_code", [
    (PORT_INT, PORT_INT_VECTOR),
    (PORT_FLOAT, PORT_FLOAT_VECTOR),
    (PORT_BOOL, PORT_BOOL_VECTOR),
])
def test_shape_follows_the_data_but_the_element_type_does_not(
        declared, vector_code):
    """Vector-ness is shape, not type: an array write is still allowed."""
    p = Port()
    p.set_value_type(declared)
    p.value = [1, 0, 2]
    assert p.get_value_type() == vector_code
    assert p.get_is_vector() is True


def test_a_bare_port_defaults_to_float():
    """chinet defaulted to the int code because the first write retyped it.

    With a declared type that default would truncate the first float anybody
    stored, so a bare port is a float port -- which is what nearly every
    caller wants (fitting parameters, spectra, residuals).
    """
    p = Port()
    assert p.get_value_type() == PORT_FLOAT
    p.value = [1.5, 2.5]
    assert list(p.value) == [1.5, 2.5]


def test_set_value_type_is_the_only_way_to_change_a_declared_type():
    """It is the declaration path, which is why loading a document uses it."""
    p = Port(value=5)
    p.value = 9.9
    assert p.get_value_type() == PORT_INT and p.value == 9
    p.set_value_type(PORT_FLOAT)
    p.value = 9.9
    assert p.value == 9.9


# --- one store, and a derived view of it ----------------------------------
#
# Storage is a single slot buffer plus the element type: a float slot holds
# the double, an int or bool slot holds the int64 bit pattern. The double view
# a non-float port hands out is built on demand and thrown away on every
# write. These tests are about that view being derived rather than a second
# store -- the previous design kept two arrays and had to answer "which is
# right" on every path, and answered wrong three times.


def test_an_integer_port_reads_as_its_value_not_its_bit_pattern():
    """The slot holds an int64 bit pattern; reading it as a double must not.

    Reinterpreting the slot rather than converting it would come back as a
    denormal near 1e-308, not as 5.
    """
    p = Port(value=5)
    assert p.get_value() == 5.0
    assert list(np.asarray(Port(value=[5, 6]).get_value_view())) == [5.0, 6.0]


def test_the_double_view_of_an_integer_port_follows_writes():
    """A cache that survived a write would be the old bug in a new place."""
    p = Port(value=5)
    assert p.get_value() == 5.0
    p.value = 9
    assert p.get_value() == 9.0
    p.value = [11, 12]
    assert list(np.asarray(p.get_value_view())) == [11.0, 12.0]


def test_the_double_view_follows_a_type_change():
    p = Port(value=[3, 4])
    assert list(np.asarray(p.get_value_view())) == [3.0, 4.0]
    p.set_value_type(PORT_BOOL_VECTOR)
    assert list(np.asarray(p.get_value_view())) == [1.0, 1.0]


def test_the_double_view_of_a_linked_integer_port_is_the_source():
    source = Port(value=7, is_output=True)
    sink = Port(value=0.0)
    sink.set_link(source)
    assert sink.get_value() == 7.0
    source.value = 8
    assert sink.get_value() == 8.0


def test_a_wide_integer_survives_being_pushed_to_a_follower():
    """propagate_to_followers pushes typed, not through the double view.

    Pushing doubles would round the value on the way into the follower's own
    storage -- the same mistake the old mirror made, one layer further out.
    """
    source = Port(value=0, is_output=True)
    follower = Port(value=0)
    follower.set_link(source)
    source.value = 2 ** 60 + 7
    assert follower.value == 2 ** 60 + 7


# --- the document round-trip ----------------------------------------------


@pytest.mark.parametrize("value", [True, False])
def test_a_bool_scalar_round_trips_through_json(value):
    """`true` is not a number, so the loader has to accept booleans itself.

    The value is restored *before* the type code, so a loader that skipped a
    JSON boolean would read every saved flag back as `False`.
    """
    p = Port(value=value)
    document = p.get_json()
    assert ('"value": %s' % str(value).lower()) in document.replace("'", '"')
    restored = Port()
    restored.read_json(document)
    assert restored.value is value
    assert restored.get_value_type() == PORT_BOOL


def test_a_bool_vector_round_trips_through_json():
    p = Port(value=True)
    p.value = [1, 0, 1]
    restored = Port()
    restored.read_json(p.get_json())
    assert restored.get_value_type() == PORT_BOOL_VECTOR
    assert list(restored.value) == [True, False, True]


@pytest.mark.parametrize("value,code", [
    (7, PORT_INT), (1.25, PORT_FLOAT), ([1, 2], PORT_INT_VECTOR),
    # A float vector built by the *constructor* reports the scalar float code,
    # not PORT_FLOAT_VECTOR. That is chinet's quirk and it is load-bearing --
    # `test_kwargs_ctor_chisurf_shape` pins it as the shape chisurf's
    # Parameter builds -- so it is asserted rather than corrected. Only the
    # setter path writes code 3.
    ([1.5, 2.5], PORT_FLOAT),
])
def test_the_other_types_still_round_trip(value, code):
    p = Port(value=value)
    assert p.get_value_type() == code
    assert p.get_is_vector() is isinstance(value, list)
    restored = Port()
    restored.read_json(p.get_json())
    assert restored.get_value_type() == code
    np.testing.assert_allclose(np.atleast_1d(restored.value),
                               np.atleast_1d(np.asarray(value, dtype=float)))


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q", "-p", "no:cacheprovider"]))
