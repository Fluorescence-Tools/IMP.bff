import IMP.bff
from typing import Dict


class Session(IMP.bff.MongoObject):

    nodes: Dict[str, IMP.bff.Node] = None

    def __init__(
            self,
            nodes: Dict[str, IMP.bff.Node] = None
    ):
        super(Session, self).__init__()
        self.nodes = nodes

    def add_node(self, name: str, node: IMP.bff.Node):
        pass

    def get_nodes(self) -> Dict[str, IMP.bff.Node]:
        pass

    def write_to_db(self):
        pass

    def read_from_db(self, oid_string) -> bool:
        rv = super(Session, self).read_from_db(oid_string)
        rv &= self.create_and_connect_objects_from_oid_array(
            self.document, "nodes", self.nodes)
        return rv

    def create_port(self, port_template: dict, port_key: str):
        pass

    def create_node(self, node_template: dict, node_key: str):
        pass

    def read_session_template(self, json_string: str):
        pass

    def get_session_template(self) -> str:
        pass

    def link_nodes(
        self,
        node_name: str,
        port_name: str,
        target_node_name: str,
        target_port_name: str
    ) -> bool:
        pass

