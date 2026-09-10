/* [bff] The shared swig macros, defined once, sorted first. */

/* The standard-library typemaps the topic files rely on. In the standalone
   build the order happens to work out; in the module build the entry
   %includes the topics alphabetically, and avbuilder.i -- the alphabetically
   first consumer of std::map -- reached its %template before types.i
   delivered the map typemaps. SWIG's library files are idempotent. */
%include <std_map.i>
%include <std_attribute.i>

%define %attribute_np(Class, Type, Name, GetMethod, SetMethod...)
    %extend Class {
    #if #SetMethod != ""
        %pythoncode
        {
            Name = property(
                lambda x: np.array(x.GetMethod()),
                SetMethod
            )
        }
    #else
        %pythoncode
        {
            Name = property(
                    lambda x: np.array(x.GetMethod())
            )
        }
    #endif
    }
%enddef

// The same, for a getter that returns a flat view of a `(n, cols)` table:
// the reshape is the one thing a numpy caller cannot be left to guess, and
// doing it here means it is written once rather than in every wrapper (the
// C++ side keeps the flat view, which is what its own callers want).
// The same for a `(rows, -1, cols)` cube, where `rows` is another attribute
// of the object: a per-conformer atom set is `(n_rotamers, n_atoms, 3)` and
// neither trailing dimension can be inferred from the flat length alone.
// The same where the column count is another attribute of the object rather
// than a literal: a pair matrix is `(n1, n2)` and both come from the value.
// `reshape(-1, cols)` and not `reshape(rows, cols)`, so an array a run left
// empty -- `k_fret` without a lifetime -- still answers with its own shape.
%define %attribute_np2v(Class, Type, Name, GetMethod, ColsAttribute)
    %extend Class {
        %pythoncode
        {
            Name = property(
                    lambda x: np.asarray(x.GetMethod()).reshape(
                            -1, x.ColsAttribute)
            )
        }
    }
%enddef

%define %attribute_np3(Class, Type, Name, GetMethod, RowsAttribute, Cols,
                       SetMethod...)
    %extend Class {
    #if #SetMethod != ""
        %pythoncode
        {
            Name = property(
                    lambda x: np.asarray(x.GetMethod()).reshape(
                            x.RowsAttribute, -1, Cols),
                    SetMethod
            )
        }
    #else
        %pythoncode
        {
            Name = property(
                    lambda x: np.asarray(x.GetMethod()).reshape(
                            x.RowsAttribute, -1, Cols)
            )
        }
    #endif
    }
%enddef

%define %attribute_np2(Class, Type, Name, GetMethod, Cols, SetMethod...)
    %extend Class {
    #if #SetMethod != ""
        %pythoncode
        {
            Name = property(
                    lambda x: np.asarray(x.GetMethod()).reshape(-1, Cols),
                    SetMethod
            )
        }
    #else
        %pythoncode
        {
            Name = property(
                    lambda x: np.asarray(x.GetMethod()).reshape(-1, Cols)
            )
        }
    #endif
    }
%enddef

%define %attribute_py(Class, Type, Name, GetMethod, SetMethod...)
%extend Class {
#if #SetMethod != ""
        %pythoncode
        {
            Name = property(GetMethod, SetMethod)
        }
#else
        %pythoncode
        {
            Name = property(GetMethod)
        }
#endif
}
%enddef


%define %class_callable(Class, Method)
    %extend Class {
        %pythoncode
        {
            def __call__(self, *args, **kwargs):
                self.Method(*args, **kwargs)
        }
    }
%enddef

