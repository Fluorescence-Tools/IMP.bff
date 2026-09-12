"""A node as a consumer writes one: imports, a module constant, an array."""
import numpy as np

import IMP.bff as bff

SCALE = 3.0
OFFSET = np.array([0.5])


class Scaled(bff.GraphNode):
    def evaluate(self):
        v = np.asarray(self.get_input_port("x").get_value_view(), dtype=float)
        self.get_output_port("y").set_value_vector((v * SCALE + OFFSET).tolist())
