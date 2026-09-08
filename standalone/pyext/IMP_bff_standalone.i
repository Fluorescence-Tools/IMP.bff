/*
 * The standalone build's SWIG entry: IMP.bff without IMP (PRD-137 step 6c).
 *
 * IMP's module tooling writes a module's entry around swig.i-in -- the kernel
 * interface, every dependency's, the exception translation, numpy. This file
 * is that entry for the core alone: the preamble, the IMP_SWIG_* equivalents
 * (IMP_bff_standalone.macros.i), and then core.i verbatim -- the same file
 * the IMP build wraps before layer.i. The generated shadow module is
 * installed as IMP/bff/__init__.py next to _IMP_bff.
 */
%module(directors="1", moduleimport="from . import _IMP_bff") IMP_bff
%feature("autodoc", 1);
%feature("kwargs", 1);

%{
#define SWIG_FILE_WITH_INIT
#include <IMP/bff/bff_config.h>
#include <IMP/bff/Base.h>
#include "bff_core_headers.h"
%}

/* The SWIG side of the module's vocabulary: the export/namespace macros and
   Base.h's standalone branch (IMP_VALUES, IMP_SHOWABLE_INLINE, ...). IMP's
   kernel interface provides these to a module's SWIG; here the headers do. */
/* numpy.i (pulled in by core.i) needs the C API initialised once per module;
   IMP's kernel does this in IMP_kernel.import_numpy.i. */
%{
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
%}
%init %{
    import_array();
%}

/* Stream insertion is a C++ convenience (IMP_SHOWABLE_INLINE writes it as a
   friend); Python has show()/__str__, never __lshift__. */
%ignore operator<<;
%ignore IMP::operator<<;
%ignore IMP::bff::operator<<;

%include <std_string.i>
%implicitconv;  /* a dict is a ProbeMap, a list an FFComponentSpecList -- as under IMP's tooling */

/* The IMP_SWIG_* equivalents and, first of all, the exception translation:
   %exception covers only what is wrapped after it, and bff_config.h
   (get_data_path, which throws) comes next. */
%include "IMP_bff_standalone.macros.i"

%include <IMP/bff/bff_config.h>
/* Base.h hides its C++ exception classes from SWIG; the module makes Python
   ones (IMP_bff_standalone.macros.i), as IMP's kernel does. */
%ignore IMP::Pointer;
%ignore IMP::Object::ref;
%ignore IMP::Object::unref;
%ignore IMP::Object::release_ref;
%include <IMP/Object.h>
%include <IMP/bff/Base.h>

%pythoncode %{
# The standalone core: no IMP.atom, no IMP.Model. What is here is what
# fits fluorescence data, samples labels and reads structures on its own.
IMPBFF_STANDALONE = True
%}

%include "IMP_bff.core.i"

%pythoncode %{
def get_module_version():
    return _IMP_bff.get_module_version()

def get_module_name():
    return "IMP::bff"

def get_build():
    """Which IMP.bff this is: "core" (no IMP, this build) or "imp" (the IMP
    module, which adds the connection layer)."""
    return "core"

def _warn_if_this_replaced_the_imp_module():
    # The pip core and the conda imp.bff package both own site-packages/IMP/bff/.
    # When an IMP kernel is importable next to this core, the pip install has
    # overwritten the IMP module's files: the connection layer's names are gone
    # until the conda package is reinstalled.
    import importlib.util
    try:
        found = importlib.util.find_spec("IMP.atom") is not None
    except Exception:
        found = False
    if found:
        import warnings
        warnings.warn(
            "IMP.bff is the IMP-free core (package 'bff'), but IMP itself is "
            "installed here: the connection layer (get_av_from_structure, the "
            "restraints, ...) is unavailable. `conda install imp.bff` restores "
            "the IMP module build.", RuntimeWarning, stacklevel=2)


_warn_if_this_replaced_the_imp_module()
%}
