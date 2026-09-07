"""Deterministic A/B parity: bff's C++ runtime against frozen Python chinet.

Every test drives both runtimes through the same operations and compares the
full observable state after each one. A failure here is a bug in bff.
"""

import unittest

from ab_driver import CYCLE_REJECTED, Pair


class ABScenarioTests(unittest.TestCase):
    def setUp(self):
        self.pair = Pair()

    def test_construction_parity(self):
        """Positional, keyword and mixed construction agree."""
        self.pair.both("make_port", "a", 1.0)
        self.pair.both("make_port", "b", 3, name="counter")
        self.pair.both("make_port", "v", [1.0, 2.0, 3.0], name="spectrum")
        self.pair.compare("construction")

    def test_scalar_writes(self):
        self.pair.both("make_port", "a", 1.0)
        for value in (0.0, -2.5, 1e6, 3, 0.1):
            self.pair.both("@a.set_value", value)
            self.pair.compare(f"set_value({value})")

    def test_shape_changes(self):
        """Scalar -> vector -> scalar, the shape-following quirk included."""
        self.pair.both("make_port", "a", 1.0)
        self.pair.both("@a.set_value", [1.0, 2.0, 3.0])
        self.pair.compare("scalar to vector")
        self.pair.both("@a.set_value", [4.0, 5.0])
        self.pair.compare("vector resize")
        self.pair.both("@a.set_value", 7.0)
        self.pair.compare("vector to scalar")

    def test_bounds_and_clipping(self):
        self.pair.both("make_port", "a", 0.5)
        self.pair.both("@a.set_bounds", 0.0, 1.0)
        self.pair.both("@a.set_is_bounded", True)
        self.pair.compare("bounded")
        for value in (3.5, -2.0, 0.25):
            self.pair.both("@a.set_value", value)
            self.pair.compare(f"clip({value})")

    def test_unbounded_reports_none_on_both(self):
        """The documented nan-vs-None divergence, mapped at the boundary."""
        self.pair.both("make_port", "a", 0.5)
        self.pair.both("@a.set_bounds", 0.0, 1.0)
        self.pair.both("@a.set_is_bounded", False)
        self.pair.compare("bounds set but not enforced")
        self.assertEqual(self.pair.ref.port("a").state({})["bounds"], (None, None))

    def test_flags(self):
        self.pair.both("make_port", "a", 1.0)
        for flag in (True, False, True):
            self.pair.both("@a.set_fixed", flag)
            self.pair.compare(f"fixed={flag}")
            self.pair.both("@a.set_reactive", flag)
            self.pair.compare(f"reactive={flag}")

    def test_link_follow_and_unlink(self):
        self.pair.both("make_port", "up", 1.0)
        self.pair.both("make_port", "down", 9.0)
        self.pair.both("@down.link", "@up")
        self.pair.compare("linked")
        self.pair.both("@up.set_value", 4.25)
        self.pair.compare("upstream write propagates")
        self.pair.both("@down.unlink")
        self.pair.compare("unlinked")

    def test_link_follows_vector_shape(self):
        """chinet's quirk: a scalar follower reads its upstream's vector."""
        self.pair.both("make_port", "up", 1.0)
        self.pair.both("make_port", "down", 9.0)
        self.pair.both("@down.link", "@up")
        self.pair.both("@up.set_value", [1.0, 2.0, 3.0])
        self.pair.compare("follower reads vector")

    def test_direct_cycle_rejected_by_both(self):
        self.pair.both("make_port", "a", 1.0)
        self.pair.both("make_port", "b", 2.0)
        self.pair.both("@b.link", "@a")
        results = self.pair.both("@a.link", "@b")
        self.assertEqual(results, [CYCLE_REJECTED, CYCLE_REJECTED])
        self.pair.compare("cycle refused")

    def test_long_cycle_rejected_by_both(self):
        labels = [f"p{i}" for i in range(6)]
        for label in labels:
            self.pair.both("make_port", label, 1.0)
        for downstream, upstream in zip(labels[1:], labels):
            self.pair.both(f"@{downstream}.link", f"@{upstream}")
        self.pair.compare("chain built")
        results = self.pair.both(f"@{labels[0]}.link", f"@{labels[-1]}")
        self.assertEqual(results, [CYCLE_REJECTED, CYCLE_REJECTED])
        self.pair.compare("long cycle refused")

    def test_relink_to_a_different_upstream(self):
        for label, value in (("a", 1.0), ("b", 2.0), ("c", 3.0)):
            self.pair.both("make_port", label, value)
        self.pair.both("@c.link", "@a")
        self.pair.compare("linked to a")
        self.pair.both("@c.link", "@b")
        self.pair.compare("relinked to b")

    def test_linked_chain_propagates_through_depth(self):
        labels = [f"n{i}" for i in range(8)]
        for label in labels:
            self.pair.both("make_port", label, 0.0)
        for downstream, upstream in zip(labels[1:], labels):
            self.pair.both(f"@{downstream}.link", f"@{upstream}")
        self.pair.both(f"@{labels[0]}.set_value", 5.5)
        self.pair.compare("depth-8 propagation")

    def test_bounded_follower_keeps_its_own_bounds(self):
        self.pair.both("make_port", "up", 10.0)
        self.pair.both("make_port", "down", 0.5)
        self.pair.both("@down.set_bounds", 0.0, 1.0)
        self.pair.both("@down.set_is_bounded", True)
        self.pair.both("@down.link", "@up")
        self.pair.compare("bounded follower on out-of-range upstream")

    def test_nodes_with_ports(self):
        self.pair.both("make_node", "fit")
        self.pair.both("make_port", "x", 3.0)
        self.pair.both("make_port", "out", 0.0)
        self.pair.both("@fit.add_input_port", "x", "@x")
        self.pair.both("@fit.add_output_port", "out", "@out")
        self.pair.compare("node wired")


if __name__ == "__main__":
    unittest.main()
