/* Director needed for wrapping of callback  functions*/
%feature("director") CnNodeCallback;

%shared_ptr(IMP::bff::CnPort)
%shared_ptr(IMP::bff::CnNode)
%shared_ptr(IMP::bff::CnNodeCallback)

%template(PairStringCnPortPtr) std::pair<std::string, std::shared_ptr<IMP::bff::CnPort>>;
%template(VecPairStringCnPortPtr) std::vector<std::pair<std::string, std::shared_ptr<IMP::bff::CnPort>>>;
%template(MapStringCnPortPtr) std::map<std::string, std::shared_ptr<IMP::bff::CnPort>>;
%template(ListCnNodePtr) std::list<std::shared_ptr<IMP::bff::CnNode>>;
%template(MapStringDouble) std::map<std::string, double>;

%attribute(IMP::bff::CnNode, bool, is_valid, is_valid);

%include "IMP/bff/CnNodeCallback.h"
%include "IMP/bff/CnNode.h"


%extend IMP::bff::CnNode {

    std::string __repr__(){
        std::ostringstream os;
        os << "CnNode(";
        for (auto &v : $self->get_input_ports()) os << (v.first)<< ",";
        os << "->";
        for (auto &v : $self->get_output_ports()) os << (v.first)<< ",";
        os << ")";
        return os.str();
    }

    %pythoncode %{

        @property
        def inputs(self):
            return dict(self.get_input_ports())

        @property
        def outputs(self):
            return dict(self.get_output_ports())

        @property
        def ports(self):
            return list(self.get_ports())
        
        @ports.setter
        def ports(self, v):
            for key in v:
                if isinstance(v[key], IMP.bff.CnPort):
                    p = v[key]
                elif isinstance(v[key], dict):
                    p = IMP.bff.CnPort(**v[key])
                else:
                    p = IMP.bff.CnPort(value=v[key])
                p.name = key
                self.add_port(p, p.is_output, True)

        def set_python_callback_function(self, cb, reactive_inputs = True, reactive_outputs = True):
            # type: (Callable, bool, bool) -> None
            import sys
            if sys.version_info[0] > 2:
                import inspect
                from inspect import signature as sig
                empty_name = inspect.Signature.empty
                empty_type = inspect._empty
                _use_py2 = True
            else:
                import funcsigs
                from funcsigs import signature as sig
                empty_name = funcsigs._empty
                empty_type = funcsigs._empty
                _use_py2 = False

            class CnNodeCallbackPython(IMP.bff.CnNodeCallback):
        
                def __init__(self, cb_function, *args, **kwargs):
                    super(CnNodeCallbackPython, self).__init__(*args, **kwargs)
                    self._cb = cb_function
        
                def run(self, inputs, outputs):
                    input_dict = dict(
                        [(inputs[i].name, inputs[i].value) for i in inputs]
                    )
                    output = self._cb(**input_dict)
                    if isinstance(output, dict):
                        for key, ov in zip(outputs, output):
                            outputs[key].value = output[key]
                    elif isinstance(output, tuple):
                        for key, ov in zip(outputs, output):
                            outputs[key].value = ov
                    else:
                        next(iter(outputs.values())).value = output
        
            if isinstance(cb, str):
                code_obj = compile(cb, '<string>', 'exec')
                self.set_string('source', cb)
                for o in code_obj.co_consts:
                    if isinstance(o, types.CodeType):
                        cb = types.FunctionType(o, globals())
                        break
            else:
                self.set_string('source', inspect.getsource(cb))

            _sig = sig(cb)
            input_ports = list()
            for parameter_name in _sig.parameters:
                o = _sig.parameters[parameter_name]
                if o.default is empty_name:
                    if o.annotation is empty_type:
                        print("WARNING no type specified for %s IMP.bff.CnPort." % parameter_name)
                        value = np.array([1.0], dtype=np.float64)
                    elif o.annotation == 'int':
                        value = 1
                    elif 'ndarray' == o.annotation:
                        value = np.array([1.0], dtype=np.float64)
                    elif o.annotation == 'float':
                        value = 1.0
                    else:
                        print("WARNING no type specified for %s IMP.bff.CnPort." % parameter_name)
                        value = np.array([1.0], dtype=np.float64)
                else:
                    value = o.default
                names = [str(n) for n, _ in self.ports]
                if str(o.name) not in names:
                    p = IMP.bff.CnPort(name=o.name, value=value, is_output=False, is_reactive=reactive_inputs)
                    self.add_input_port(p)
                    input_ports.append(p)
                else:
                    input_ports.append(self.get_input_ports()[o.name])

            input_values = dict([(p.name, p.value) for p in input_ports])
            returned_value = cb(**input_values)

            output_values = list()
            output_names = list()
            if isinstance(returned_value, dict):
                for ok in returned_value:
                    output_values.append(returned_value[ok])
                    output_names.append(ok)
            elif isinstance(returned_value, tuple):
                output_values = returned_value
                for i in range(len(returned_value)):
                    output_names.append("out_" + str(i).zfill(2))
            else:
                output_names = 'out_00',
                output_values = returned_value,
            for on, ov in zip(output_names, output_values):
                if on not in self.get_output_ports().keys():
                    po = IMP.bff.CnPort(
                            name=on,
                            value=ov,
                            is_output=True,
                            is_reactive=reactive_outputs
                    )
                    self.add_output_port(po)
            cb_instance = CnNodeCallbackPython(cb_function=cb)
            cb_instance.__disown__()
            self.set_callback(cb_instance)
        
        callback_function = property(None, set_python_callback_function)
        
        def __getattr__(self, name):
            in_input = name in self.inputs.keys()
            in_output = name in self.outputs.keys()
            if in_input and in_output:
                raise KeyError("Ambiguous access. %s an input and output parameter.")
            elif not in_output and not in_input:
                raise AttributeError
            else:
                if in_input:
                    return self.inputs[name].value
                else:
                    return self.outputs[name].value
        
        def __setattr__(self, name, value):
            try:
                in_input = name in self.inputs.keys()
                in_output = name in self.outputs.keys()
                if in_input and in_output:
                    raise KeyError("Ambiguous access. %s an input and output parameter.")
                elif not in_output and not in_input:
                    raise AttributeError
                else:
                    if in_input:
                        self.inputs[name].value = value
                    else:
                        self.outputs[name].value = value
            except:
                IMP.bff.CnMongoObject.__setattr__(self, name, value)
        
        def __call__(self):
            return self.evaluate()
        
        def __init__(self,
                     obj=None,
                     name = "", #type: str
                     ports = None, #type: dict
                     callback_function=None,
                     reactive_inputs=True,
                     reactive_outputs=True,
                     *args, **kwargs
        ):
            this = _IMP_bff.new_CnNode(*args, **kwargs)
            try:
                self.this.append(this)
            except:
                self.this = this
            self.register_instance(None)
            if isinstance(ports, dict):
                self.ports = ports
            if isinstance(obj, str):
                name = obj
            elif callable(obj):
                callback_function = obj
            if callable(callback_function) or isinstance(callback_function, str):
                self.set_python_callback_function(
                    cb=callback_function,
                    reactive_inputs=reactive_inputs,
                    reactive_outputs=reactive_outputs
                )
            if len(name) > 0:
                self.name = name
            else:
                if callable(callback_function):
                    self.name = callback_function.__name__

   %}
}

