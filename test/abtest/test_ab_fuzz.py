"""Randomised A/B fuzzing: identical op sequences, both runtimes, compared each step.

Deterministic scenarios only cover what someone thought to write down. These
sequences are seeded random, so a failure names the exact op index that first
diverged and reproduces from the seed.
"""

import random
import unittest

from ab_driver import Pair

N_SEQUENCES = 150
N_OPS = 200
N_PORTS = 6


def _random_value(rng):
    kind = rng.random()
    if kind < 0.45:
        return rng.uniform(-100.0, 100.0)
    if kind < 0.6:
        return rng.randint(-50, 50)
    length = rng.randint(1, 5)
    return [rng.uniform(-10.0, 10.0) for _ in range(length)]


class ABFuzzTests(unittest.TestCase):
    def _run_sequence(self, seed):
        rng = random.Random(seed)
        pair = Pair()
        labels = [f"p{i}" for i in range(N_PORTS)]
        for label in labels:
            pair.both("make_port", label, rng.uniform(-1.0, 1.0))
        pair.compare(f"seed={seed} construction")

        for index in range(N_OPS):
            label = rng.choice(labels)
            choice = rng.random()
            if choice < 0.30:
                op = f"set_value({label})"
                pair.both(f"@{label}.set_value", _random_value(rng))
            elif choice < 0.42:
                op = f"set_bounds({label})"
                low = rng.uniform(-20.0, 0.0)
                pair.both(f"@{label}.set_bounds", low, low + rng.uniform(0.1, 40.0))
            elif choice < 0.52:
                flag = rng.random() < 0.5
                op = f"is_bounded({label},{flag})"
                pair.both(f"@{label}.set_is_bounded", flag)
            elif choice < 0.62:
                flag = rng.random() < 0.5
                op = f"fixed({label},{flag})"
                pair.both(f"@{label}.set_fixed", flag)
            elif choice < 0.70:
                flag = rng.random() < 0.5
                op = f"reactive({label},{flag})"
                pair.both(f"@{label}.set_reactive", flag)
            elif choice < 0.90:
                other = rng.choice(labels)
                op = f"link({label}->{other})"
                if other != label:
                    results = pair.both(f"@{label}.link", f"@{other}")
                    self.assertEqual(
                        results[0], results[1],
                        f"seed={seed} op {index} {op}: link outcome differs "
                        f"(ref={results[0]!r} bff={results[1]!r})")
            else:
                op = f"unlink({label})"
                pair.both(f"@{label}.unlink")
            pair.compare(f"seed={seed} op {index}: {op}")

    def test_fuzzed_operation_sequences(self):
        for seed in range(N_SEQUENCES):
            with self.subTest(seed=seed):
                self._run_sequence(seed)

    def test_self_link_refused_identically(self):
        """A port linking to itself is the degenerate cycle."""
        for seed in range(5):
            pair = Pair()
            pair.both("make_port", "a", float(seed))
            results = pair.both("@a.link", "@a")
            self.assertEqual(results[0], results[1])
            pair.compare(f"self-link seed={seed}")


if __name__ == "__main__":
    unittest.main()
