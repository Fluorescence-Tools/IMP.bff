import utils
import os
import unittest
import sys
import json

import numba as nb
import numpy as np
import IMP
import IMP.bff

from constants import *

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
utils.set_search_paths(TOPDIR)


class CallbackNodePassOn(IMP.bff.CnNodeCallback):

    def __init__(self, *args, **kwargs):
        super(CallbackNodePassOn, self).__init__(*args, **kwargs)

    def run(self, inputs, outputs):
        outputs["outA"].value = inputs["inA"].value


class NodeCallbackMultiply(IMP.bff.CnNodeCallback):

    def __init__(self, *args, **kwargs):
        super(NodeCallbackMultiply, self).__init__(*args, **kwargs)

    def run(self, inputs, outputs):
        print("run:NodeCallbackMultiply")
        mul = 1.0
        for key in inputs:
            mul *= inputs[key].value
            print(mul)
        outputs["portC"].value = mul


def node_cb(func, ):

    class NodeCB(IMP.bff.CnNodeCallback):

        def __init__(self, *args, **kwargs):
            super(NodeCB, self).__init__(*args, **kwargs)

        def run(self, inputs, outputs):
            mul = 1.0
            func(inputs[key].value)
            for in_key in inputs:
                mul *= inputs[key].value
                print(mul)
            outputs["portC"].value = mul

    def wrapper():
        print("Something is happening before the function is called.")
        func()
        print("Something is happening after the function is called.")
    return wrapper


class Tests(unittest.TestCase):

    def test_node_init(self):
        port_d = {
            "value": 0.0,
            "fixed": False,
            "is_output": True
        }
        node_d = {
            'name': 'NodeName',
            'ports': {
                'inA': IMP.bff.CnPort(7.0),
                'inB': IMP.bff.CnPort(13.0),
                'outC': IMP.bff.CnPort(**port_d)
            }
        }
        node_with_ports = IMP.bff.CnNode(**node_d)
        self.assertEqual(
            node_with_ports.get_input_ports().keys(),
            ['inA', 'inB']
        )
        values = np.hstack(
            [v.value for n, v in node_with_ports.ports]
        )
        np.testing.assert_allclose(values, np.array([7.0, 13.0, 0.0]))

    def test_node_ports(self):
        node = IMP.bff.CnNode()
        portA = IMP.bff.CnPort(17, name="portA")
        node.add_input_port(portA)
        portB = IMP.bff.CnPort(23, name="portB")
        node.add_output_port(portB)

        self.assertEqual(portA, node.get_input_port(portA.name))
        self.assertEqual(portB, node.get_output_port(portB.name))

    def test_node_array_callback(self):
        inA = IMP.bff.CnPort([2., 3., 4.], name="inA")
        inB = IMP.bff.CnPort([2., 3., 4.], name="inB")
        outA = IMP.bff.CnPort([0., 0., 0.], name="outA")

        node = IMP.bff.CnNode()
        node.add_input_port(inA)
        node.add_input_port(inB)
        node.add_output_port(outA)
        multiply = lambda inA, inB: inA * inB
        node.callback_function = multiply

        node.evaluate()
        # np.testing.assert_allclose(inA.value * inB.value, outA.value)

    def test_callback_2(self):
        node = IMP.bff.CnNode(
            ports={
                "inA": IMP.bff.CnPort(
                    value=[2., 3., 4.],
                    fixed=False,
                    is_output=False
                ),
                "inB": IMP.bff.CnPort(
                    value=[2., 3., 4.],
                    fixed=False,
                    is_output=False
                ),
                "outA": IMP.bff.CnPort(
                    value=0,
                    fixed=False,
                    is_output=True
                )
            }
        )
        multiply = lambda inA, inB: inA * inB
        node.set_python_callback_function(multiply)
        node.callback_function = multiply
        self.assertEqual(node.is_valid, False)

        node.evaluate()
        self.assertEqual(node.is_valid, True)

        outA = node.outputs["outA"]
        inA = node.inputs["inA"]
        inB = node.inputs["inB"]

        np.testing.assert_allclose(inA.value * inB.value, outA.value)

    def test_node_python_callback_1(self):
        """Test chinet Node class python callbacks"""
        node = IMP.bff.CnNode()
        portIn1 = IMP.bff.CnPort(55, name="portA")
        node.add_input_port(portIn1)
        portIn2 = IMP.bff.CnPort(2, name="portB")
        node.add_input_port(portIn2)
        portOut1 = IMP.bff.CnPort(name="portC")
        node.add_output_port(portOut1)

        cb = NodeCallbackMultiply()
        cb.thisown = 0
        node.set_callback(cb)
        node.evaluate()

        self.assertEqual(portOut1.value, portIn1.value * portIn2.value)

    def test_node_python_callback_2(self):
        """Test chinet Node class python callbacks"""
        node = IMP.bff.CnNode()
        portIn1 = IMP.bff.CnPort(55, name="portA")
        node.add_input_port(portIn1)
        portIn2 = IMP.bff.CnPort(2, name="portB")
        node.add_input_port(portIn2)
        portOut1 = IMP.bff.CnPort(name="portC")
        node.add_output_port(portOut1)
        m = NodeCallbackMultiply()
        node.set_callback(m)
        # node.set_callback(NodeCallbackMultiply().__disown__())
        node.evaluate()
        self.assertEqual(
            portOut1.value,
            portIn1.value * portIn2.value
        )

    @unittest.skipUnless(CONNECTS, "Cloud not connect to DB")
    def test_node_write_to_db(self):
        node = IMP.bff.CnNode(
            ports={
                'portA': IMP.bff.CnPort(55),
                'portB': IMP.bff.CnPort(2),
                'portC': IMP.bff.CnPort()
            }
        )
        m = NodeCallbackMultiply()
        node.set_callback(m)
        self.assertEqual(node.connect_to_db(**DB_DICT), True)
        self.assertEqual(node.write_to_db(), True)

    @unittest.skipUnless(CONNECTS, "Cloud not connect to DB")
    def test_node_restore_from_db(self):
        # Make new node that will be written to the DB
        node = IMP.bff.CnNode(
            ports={
                'portA': IMP.bff.CnPort(13.0),
                'portB': IMP.bff.CnPort(2.0),
                'portC': IMP.bff.CnPort(1.0, is_output=True)
            }
        )
        m = NodeCallbackMultiply()
        node.set_callback(m)
        node.connect_to_db(**DB_DICT)
        node.write_to_db()

        # Restore the Node
        node_restore = IMP.bff.CnNode()
        node_restore.connect_to_db(**DB_DICT)
        node_restore.read_from_db(node.oid)

        # compare the dictionaries of the nodes
        dict_restore = json.loads(node_restore.get_json())
        dict_original = json.loads(node.get_json())
        self.assertEqual(dict_restore, dict_original)

    def test_node_valid(self):
        """
        In this test the node node_1 has one inputs (inA) and one output:

                 (inA)-(node_1)-(outA)

        In this example all ports are "non-reactive" meaning, when the input
        of a node changes, the node is set to invalid. A node is set to valid
        when it is evaluated. When a node is initialized it is invalid.
        """
        out_node_1 = IMP.bff.CnPort(
            1.0,
            fixed=False,
            is_reactive=True,
            is_output=True
        )
        in_node_1 = IMP.bff.CnPort(3.0)
        node_1 = IMP.bff.CnNode(
            ports={
                'inA': in_node_1,
                'outA': out_node_1
            }
        )
        cb = CallbackNodePassOn()
        node_1.set_callback(cb)

        self.assertEqual(out_node_1.value, 1.0)
        self.assertEqual(in_node_1.value, 3.0)
        self.assertEqual(node_1.is_valid, False)
        node_1.evaluate()
        self.assertEqual(node_1.is_valid, True)

    def test_node_valid_reactive_port(self):
        """
        In this test the node node_1 has one inputs (inA) and one output:

                 (inA)-(node_1)-(outA)

        The input inA is reactive, meaning, when the value of inA changes
        the node is evaluated.
        """

        out_node_1 = IMP.bff.CnPort(
            value=1.0,
            fixed=False,
            is_output=True
        )
        in_node_1 = IMP.bff.CnPort(
            value=3.0,
            fixed=False,
            is_output=False,
            is_reactive=True
        )
        node_1 = IMP.bff.CnNode(
            ports={
                'inA': in_node_1,
                'outA': out_node_1
            }
        )
        cb = CallbackNodePassOn()
        node_1.set_callback(cb)

        self.assertEqual(node_1.is_valid, False)

        # A reactive port calls Node::evaluate when its value changes
        in_node_1.value = 12
        self.assertEqual(node_1.is_valid, True)
        self.assertEqual(out_node_1.value, 12)

    def test_node_valid_connected_nodes(self):
        """
        In this test the two nodes node_1 and node_2 are connected.

        node_1 has one inputs (inA) and one output:

                 (inA)-(node_1)-(outA)

        node_2 has one input (inA) and one output. The input of node_2 is
        connected to the output of node_1:

                (inA->(Node1, outA))-(node_2)-(outA)

        In this example all ports are "non-reactive" meaning, when the input
        of a node changes, the node is set to invalid. A node is set to valid
        when it is evaluated. When a node is initialized it is invalid.
        """
        in_node_1 = IMP.bff.CnPort(
            value=3.0,
            fixed=False,
            is_output=False,
            is_reactive=False
        )
        out_node_1 = IMP.bff.CnPort(
            value=1.0,
            fixed=False,
            is_output=True
        )
        node_1 = IMP.bff.CnNode(
            ports={
                'inA': in_node_1,
                'outA': out_node_1
            }
        )
        f = lambda inA: inA
        node_1.callback_function = f

        self.assertEqual(in_node_1.value, 3.0)
        self.assertEqual(out_node_1.value, 1.0)
        self.assertEqual(node_1.is_valid, False)

        node_1.evaluate()
        self.assertEqual(node_1.is_valid, True)
        self.assertEqual(in_node_1.value, out_node_1.value)

        in_node_2 = IMP.bff.CnPort(
            value=13.0,
            fixed=False,
            is_output=False,
            is_reactive=False
        )
        out_node_2 = IMP.bff.CnPort(
            value=1.0,
            fixed=False,
            is_output=True
        )
        node_2 = IMP.bff.CnNode(
            ports={
                'inA': in_node_2,
                'outA': out_node_2
            }
        )
        node_2.callback_function = f
        in_node_2.link = out_node_1

        self.assertEqual(node_2.is_valid, False)

        node_2.evaluate()

        self.assertEqual(node_2.is_valid, True)
        self.assertEqual(out_node_2.value, 3.0)

        in_node_1.value = 13
        self.assertEqual(node_1.is_valid, False)
        self.assertEqual(node_2.is_valid, False)

        node_1.evaluate()
        self.assertEqual(out_node_1.value, 13.0)

        node_2.evaluate()
        self.assertEqual(out_node_2.value, 13.0)

    def test_node_valid_connected_nodes_reactive_ports(self):
        """
        In this test the two nodes node_1 and node_2 are connected.

        node_1 has one inputs (inA) and one output:

                 (inA)-(node_1)-(outA)

        node_2 has one input (inA) and one output. The input of node_2 is
        connected to the output of node_1:

                (inA->(node_1, outA))-(node_2)-(outA)

        In this example all ports are "non-reactive" meaning, when the input
        of a node changes, the node is set to invalid. A node is set to valid
        when it is evaluated. When a node is initialized it is invalid.
        """
        in_node_1 = IMP.bff.CnPort(
            value=3.0,
            fixed=False,
            is_output=False,
            is_reactive=True
        )
        out_node_1 = IMP.bff.CnPort(
            value=1.0,
            fixed=False,
            is_output=True
        )
        node_1 = IMP.bff.CnNode(
            ports={
                'inA': in_node_1,
                'outA': out_node_1
            }
        )
        f = lambda inA: inA
        node_1.callback_function = f

        in_node_2 = IMP.bff.CnPort(
            value=1.0,
            fixed=False,
            is_output=False,
            is_reactive=True
        )
        out_node_2 = IMP.bff.CnPort(
            value=1.0,
            fixed=False,
            is_output=True
        )
        node_2 = IMP.bff.CnNode(
            ports={
                'inA': in_node_2,
                'outA': out_node_2
            }
        )
        node_2.callback_function = f
        in_node_2.link = out_node_1
        in_node_1.value = 13

        self.assertEqual(node_1.is_valid, True)
        self.assertEqual(node_2.is_valid, True)
        self.assertEqual(out_node_2.value, 13.0)

    def test_node_init_with_callback_function(self):

        def h(x, y):
            # type: (np.ndarray, np.ndarray)
            return x * y

        node = IMP.bff.CnNode(
            callback_function=h,
            name="NodeName"
        )
        v = np.arange(10, dtype=np.double)
        node.inputs['x'].value = v
        node.inputs['y'].value = v
        node.evaluate()
        np.testing.assert_allclose(node.outputs['out_00'].value, v * v)

    def call_back_setter(self):
        """Tests the setter Node.function_callback that takes Python functions

        The Ports and the Callback function of a Node can be initialized
        using a normal Python function. The parameters names of the function
        are taken as names of the input Ports. The output Port names are
        either numbered from out_00 to out_xx if the function returns a
        tuples or a single object (out_00). Or the names correspond to the
        keys of a returned dictionary.

        :return:
        """
        # one input to one output
        import IMP.bff
        node = IMP.bff.CnNode()
        f = lambda x: 2.*x
        node.callback_function = f
        x = node.inputs['x']
        out = node.outputs['out_00']
        x.reactive = True
        x.value = 11.0
        self.assertEqual(
            out.value,
            f(x.value)
        )

        # one input to many outputs
        node = IMP.bff.CnNode()

        def g(x):
            return x // 4.0, x % 4.0

        node.callback_function = g
        x = node.inputs['x']
        out_0 = node.outputs['out_00']
        out_1 = node.outputs['out_01']
        x.value = 11.0
        self.assertEqual(out_0, x // 4.0)
        self.assertEqual(out_1, x % 4.0)

        # use of numba decorated function as a Node callback
        node = IMP.bff.CnNode()

        @nb.jit
        def h(x, y):
            # type: (np.ndarray, np.ndarray) -> np.ndarray
            return x * y

        node.callback_function = h
        x = node.inputs['x']
        y = node.inputs['y']
        z = node.outputs['out_00']
        y.value = 2.0
        x.reactive = True
        x.value = np.arange(100, dtype=np.double)
        self.assertEqual(z.value, x.value * y.value)


if __name__ == '__main__':
    unittest.main()
