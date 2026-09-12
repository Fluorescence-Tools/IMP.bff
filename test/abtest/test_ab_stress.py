"""Lifetime and robustness stress: the failure mode that killed the original C++.

chinet was C++ once and crashed constantly. Ownership across the language
boundary is the usual culprit, so these tests hold the runtime under garbage
collection pressure, dense link churn and deep chains, and check that bff both
survives and still agrees with the frozen Python reference.

The known ownership rule is asserted rather than worked around: a node holds
its ports strongly but its upstream *nodes* only weakly, so a graph's
intermediate nodes must be kept alive by the caller.
"""

import gc
import random
import unittest

from ab_driver import CYCLE_REJECTED, Pair
from IMP import bff

N_PORTS = 400
CHAIN_DEPTH = 200


class ABStressTests(unittest.TestCase):
    def test_dense_link_churn_under_gc(self):
        """Constant cycle pressure with collection between every operation."""
        rng = random.Random(20260831)
        pair = Pair()
        labels = [f"s{i}" for i in range(12)]
        for label in labels:
            pair.both("make_port", label, rng.uniform(-5.0, 5.0))
        for index in range(300):
            source, target = rng.sample(labels, 2)
            results = pair.both(f"@{source}.link", f"@{target}")
            self.assertEqual(results[0], results[1],
                             f"op {index}: link outcome differs")
            if rng.random() < 0.25:
                pair.both(f"@{rng.choice(labels)}.unlink")
            if rng.random() < 0.2:
                pair.both(f"@{rng.choice(labels)}.set_value",
                          rng.uniform(-5.0, 5.0))
            gc.collect()
            pair.compare(f"gc op {index}")

    def test_many_ports_survive_collection(self):
        """Hundreds of live ports, collected repeatedly, still readable."""
        ports = [bff.GraphPort(float(i), name=f"bulk{i}") for i in range(N_PORTS)]
        for _ in range(5):
            gc.collect()
        for index, port in enumerate(ports):
            self.assertEqual(port.value, float(index))

    def test_deep_chain_propagates_after_collection(self):
        """A depth-200 link chain still follows its head after gc."""
        pair = Pair()
        labels = [f"d{i}" for i in range(CHAIN_DEPTH)]
        for label in labels:
            pair.both("make_port", label, 0.0)
        for downstream, upstream in zip(labels[1:], labels):
            pair.both(f"@{downstream}.link", f"@{upstream}")
        gc.collect()
        pair.both(f"@{labels[0]}.set_value", 42.5)
        gc.collect()
        pair.compare("deep chain after gc")

    def test_dropping_a_python_reference_keeps_the_link_alive(self):
        """A follower holds its upstream port strongly, so gc cannot orphan it.

        This is the ownership guarantee bff must provide: chisurf links
        parameters and then lets the local handle go out of scope.
        """
        upstream = bff.GraphPort(3.0, name="upstream")
        follower = bff.GraphPort(0.0, name="follower")
        follower.link = upstream
        upstream.value = 7.25
        del upstream
        for _ in range(3):
            gc.collect()
        self.assertEqual(follower.value, 7.25)
        self.assertEqual(follower.get_link().name, "upstream")

    def test_unlink_storm(self):
        """Repeated link/unlink of the same pair leaves no residue."""
        pair = Pair()
        pair.both("make_port", "a", 1.0)
        pair.both("make_port", "b", 2.0)
        for index in range(200):
            pair.both("@b.link", "@a")
            pair.both("@b.unlink")
            if index % 25 == 0:
                gc.collect()
                pair.compare(f"unlink storm {index}")
        pair.compare("unlink storm end")

    def test_cycle_pressure_never_corrupts_topology(self):
        """Rejected cycles must leave the graph exactly as it was."""
        pair = Pair()
        labels = [f"c{i}" for i in range(8)]
        for label in labels:
            pair.both("make_port", label, 1.0)
        for downstream, upstream in zip(labels[1:], labels):
            pair.both(f"@{downstream}.link", f"@{upstream}")
        before = pair.bff.state()
        for _ in range(50):
            results = pair.both(f"@{labels[0]}.link", f"@{labels[-1]}")
            self.assertEqual(results, [CYCLE_REJECTED, CYCLE_REJECTED])
            gc.collect()
        self.assertEqual(pair.bff.state(), before)
        pair.compare("after 50 refused cycles")

    def test_bulk_graph_with_nodes_survives_gc(self):
        """Nodes keep their ports alive across collection."""
        nodes = []
        for index in range(100):
            node = bff.GraphNode()
            node.name = f"n{index}"
            node.add_input_port("x", bff.GraphPort(float(index), name=f"x{index}"))
            node.add_output_port(f"n{index}", bff.GraphPort(0.0, name=f"o{index}"))
            nodes.append(node)
        for _ in range(3):
            gc.collect()
        for index, node in enumerate(nodes):
            self.assertEqual(node.get_input_port("x").value, float(index))


if __name__ == "__main__":
    unittest.main()
