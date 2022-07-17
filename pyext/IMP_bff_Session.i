%shared_ptr(IMP::bff::Session)

%template(MapStringNode) std::map<std::string, std::shared_ptr<IMP::bff::Node>>;
%ignore create_port(json port_template, std::string &port_key);
%ignore create_node(json node_template, std::string &node_key);

%include "IMP/bff/Session.h"


%extend IMP::bff::Session{

    %pythoncode %{

        @property
        def nodes(self):
        # type: () -> dict
            return dict(self.get_nodes())

        def __str__(self):
            return self.get_json(indent=2)

   %}
}
