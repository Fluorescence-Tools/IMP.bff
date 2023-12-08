import utils
import os
import unittest
import numpy as np

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)

import IMP.bff
from constants import *


class Tests(unittest.TestCase):

    def test_session_init(self):
        n1 = IMP.bff.CnNode(
            ports={
                'portA': IMP.bff.CnPort(1),
                'portB': IMP.bff.CnPort(2),
                'portC': IMP.bff.CnPort(3)
            },
        )
        n2 = IMP.bff.CnNode(
            ports={
                'inA': IMP.bff.CnPort(5),
                'inB': IMP.bff.CnPort(7),
                'outA': IMP.bff.CnPort(11, False, True),
                'outB': IMP.bff.CnPort(13, False, True)
            }
        )

        s1 = IMP.bff.Session()
        s1.add_node("nodeA", n1)
        s1.add_node("nodeB", n2)

        nodes = {
            "nodeA": n1,
            "nodeB": n2
        }
        s2 = IMP.bff.Session(nodes)

        na1 = s1.nodes['nodeA']
        na2 = s2.nodes['nodeA']
        na1.name = "new name"

        # the nodes are references
        self.assertEqual(na1.name, na2.name)

    def test_read_session_template(self):
        template_file = "examples/session_template.json"

        with open(template_file, 'r') as fp:
            json_string = fp.read()

            s = IMP.bff.Session()
            s.read_session_template(json_string)

            nodeA = s.nodes[list(s.nodes.keys())[0]]
            nodeB = s.nodes[list(s.nodes.keys())[1]]

            # test reading of nodes
            self.assertListEqual(
                list(s.nodes.keys()),
                ['nodeA', 'nodeB']
            )

            # test reading of ports
            self.assertListEqual(
                nodeA.get_input_ports().keys(),
                ['inA']
            )

            # test reading of port values
            self.assertListEqual(
                list(np.hstack([d.value for d in nodeA.get_input_ports().values()])),
                [1., 2., 3., 4., 5.]
            )

            # test evaluate
            nodeA.evaluate()

            nodeA_inA = nodeA.get_port("inA")
            nodeB_inA = nodeB.get_port("inA")

            # test links
            self.assertListEqual(
                list(nodeA_inA.value),
                list(nodeB_inA.value)
            )

            v = list(np.hstack([d.value for d in nodeA.get_output_ports().values()]))
            self.assertListEqual(v, [1., 2., 3., 4., 5.])

    def test_create_node_template(self):
        template_file = "examples/node_template.json"
        json_string = ""
        with open(template_file, 'r') as fp:
            json_string = fp.read()
            s = IMP.bff.Session()
            s.read_session_template(json_string)
        nodeA = s.create_node(json_string, "node_name")

        v = list(nodeA.get_output_port('outA').value)
        self.assertListEqual(v, [1., 2., 4., 8., 16.])


if __name__ == '__main__':
    unittest.main()