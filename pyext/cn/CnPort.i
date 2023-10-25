%shared_ptr(IMP::bff::CnPort);
%shared_ptr(IMP::bff::CnNode);

%attribute(IMP::bff::CnPort, bool, fixed, is_fixed, set_fixed);
%attribute(IMP::bff::CnPort, bool, is_output, is_output, set_port_type);
%attribute(IMP::bff::CnPort, bool, reactive, is_reactive, set_reactive);
%attribute(IMP::bff::CnPort, bool, is_linked, is_linked);
%attribute(IMP::bff::CnPort, bool, bounded, is_bounded, set_bounded);
%attribute_np(IMP::bff::CnPort, std::vector<char>, bytes, get_bytes, set_bytes);
%attribute_py(IMP::bff::CnPort, IMP::bff::CnPort, link, get_link, set_link);
%attribute_py(IMP::bff::CnPort, IMP::bff::CnNode*, node, get_node, set_node);

%include "IMP/bff/CnPort.h"

%extend IMP::bff::CnPort{

    public:
        %template(set_value_d) set_value<double>;
        %template(get_value_d) get_value<double>;
        %template(get_value_i) get_value<long>;
        %template(set_value_i) set_value<long>;
        %template(get_value_c) get_value<unsigned char>;
        %template(set_value_c) set_value<unsigned char>;

    %pythoncode %{

        @property
        def value(self):
            if self.get_value_type() == IMP.bff.BFF_PORT_VECTOR_INT:
                v = self.get_value_i()
            elif self.get_value_type() == IMP.bff.BFF_PORT_VECTOR_FLOAT:
                v = self.get_value_d()
            elif self.get_value_type() == IMP.bff.BFF_PORT_INT:
                v = self.get_value_i()[0]
            elif self.get_value_type() == IMP.bff.BFF_PORT_FLOAT:
                v = self.get_value_d()[0]
            else:
                v = None
            return v

        @value.setter
        def value(self, v):
            if not isinstance(v, np.ndarray):
                v = np.atleast_1d(v)
            if v.dtype.kind == 'i':
                self.set_value_i(v)
            else:
                self.set_value_d(v)

        @property
        def bounds(self):
            if self.bounded:
                return self.get_bounds()
            else:
                return None, None

        @bounds.setter
        def bounds(self, v):
            self.set_bounds(np.array(v, dtype=np.float64))

        def __init__(self, value=[], fixed=False, *args, **kwargs):
            this = _IMP_bff.new_CnPort(*args, **kwargs)
            try:
                self.this.append(this)
            except:
                self.this = this
            if not isinstance(value, np.ndarray):
                value = np.atleast_1d(value)
            self.value = value
            self.fixed = fixed

        def __str__(self):
            return self.get_json(indent=2)

   %}
}
