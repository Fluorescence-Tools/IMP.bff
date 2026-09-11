/* [bff] The shared swig macros, defined once, sorted first. */

/* Keyword arguments, as the standalone entry grants them globally: the
   generated module entry does not, and 432 test failures on the first
   full lane run were one thing -- every call passing a keyword
   (`linker_length=`, `points=`, `density=`, ...) rejected. This file
   sorts ahead of every bff topic, so the feature covers them all; the
   kernel and IMP's own modules declared above are unaffected (their
   entry block grants or withholds kwargs on its own). */
%feature("kwargs", 1);


%{
// numpy.i owns PyArray_API in this translation unit only when this is
// defined; without it NO_IMPORT_ARRAY makes the symbol an undefined extern
// and the module build's link fails (the standalone entry defines it in its
// own preamble; the generated module entry does not -- this file sorts first
// among the module's bff topics, so the define lands before numpy.i's
// include. The old bff module solved it the same way, in its swig.i-in).
#ifndef SWIG_FILE_WITH_INIT
#define SWIG_FILE_WITH_INIT
#endif
%}

%init %{
// Owning the pointer is only half of it: someone must initialize it. The
// standalone entry does this in its own %init; the generated module entry
// runs only IMP kernel's numpy init, which initializes kernel's
// PyArray_API, not this TU's -- every numpy-array typemap then dereferenced
// NULL and segfaulted (found on cn1: every Port.get_value_view() call).
    import_array();
%}

/* The standard-library typemaps the topic files rely on. In the standalone
   build the order happens to work out; in the module build the entry
   %includes the topics alphabetically, and avbuilder.i -- the alphabetically
   first consumer of std::map -- reached its %template before types.i
   delivered the map typemaps. SWIG's library files are idempotent. */
%include <std_map.i>
// SWIG 4.4+ ships the attribute support as python/attribute.i (it wraps
// typemaps/attribute.swg); the older std_attribute.i spelling is gone.
%include <attribute.i>

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

