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
#ifdef IMPBFF_WITH_IMP
// IMP's VectorD calls boost::distance in a constructor template, and its own
// translation units get the declaration from headers this one does not
// include. Ahead of everything, so the order the wrapper includes IMP in
// does not matter.
#include <boost/range/distance.hpp>
#endif
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

/* RMF is optional in this build (PRD-139); RmfIO.h hides its declarations
   when it is absent, and CMake passes IMPBFF_NO_RMF to SWIG as well, so the
   wrapper and the library agree on what exists. */

%include <std_string.i>
%implicitconv;  /* a dict is a ProbeMap, a list an FFComponentSpecList -- as under IMP's tooling */

/* The IMP_SWIG_* equivalents and, first of all, the exception translation:
   %exception covers only what is wrapped after it, and bff_config.h
   (get_data_path, which throws) comes next. */
%include "IMP_bff_standalone.macros.i"

%include <IMP/bff/bff_config.h>

#ifndef IMPBFF_WITH_IMP
/* Base.h hides its C++ exception classes from SWIG; the module makes Python
   ones (IMP_bff_standalone.macros.i), as IMP's kernel does. */
%ignore IMP::Pointer;
%ignore IMP::Object::ref;
%ignore IMP::Object::unref;
%ignore IMP::Object::release_ref;
%include <IMP/Object.h>
%include <IMP/bff/Base.h>
#else
/* Where IMP is linked, Base.h forwards to IMP's own headers, and those need
   IMP's SWIG interfaces to parse -- which is exactly what this module does
   not want (they make the extension import _IMP_kernel). SWIG does not
   compile anything, though: it only has to *parse* the declarations, so it
   is enough to tell it what these macros expand to. The compiler still sees
   IMP's real definitions through the %{ %} block above.
   `IMP::Object` is then an unknown base class, which SWIG says out loud and
   then ignores; the three classes that derive from it are directors, and
   their reference counting comes from IMP_SWIG_OBJECT rather than the base. */
#define IMP_SHOWABLE_INLINE(Name, how_to_show) void show(std::ostream &out=std::cout) const
#define IMP_SHOWABLE(Name) void show(std::ostream &out=std::cout) const
/* IMP's own spellings, so that what SWIG writes into the wrapper is the
   type the compiler will see. IMP_SWIG_VALUE gives IMP::Vector the same
   typemaps as std::vector. */
#define IMP_VALUES(Name, PluralName) typedef IMP::Vector<Name> PluralName;
#define IMP_OBJECTS(Name, PluralName) typedef IMP::Vector<IMP::Pointer<Name> > PluralName; typedef IMP::Vector<Name*> PluralName##Temp;
#define IMP_OBJECT_METHODS(Name) virtual ~Name();
#define IMP_OBJECT_SERIALIZE_DECL(Name)
#define IMP_DECORATORS(Name, PluralName, Parent)
#define IMP_USAGE_CHECK(expr, message)
#define IMP_WARN(message)
#define IMP_THROW(message, exception_name)
#endif

%pythoncode %{
# The standalone core: no IMP.atom, no IMP.Model. What is here is what
# fits fluorescence data, samples labels and reads structures on its own.
IMPBFF_STANDALONE = True
%}

#ifdef IMPBFF_WITH_IMP
%pythoncode %{
_IMP_BFF_WITH_IMP = True
%}
#else
%pythoncode %{
_IMP_BFF_WITH_IMP = False
%}
#endif

%include "IMP_bff.core.i"

/* The dye roads by file path, where this build links IMP (PRD-139). Their
   signatures name no IMP type, so they need none of IMP's own interfaces --
   see IMP_bff.dyedynamics.i. */
#ifdef IMPBFF_WITH_IMP
%include "IMP_bff.dyedynamics.i"
/* IMP's PDB and mmCIF readers filling the core's structure table -- the
   radius a docking score measures clashes against, and the one format the
   core cannot parse at all. Names no IMP type, so it needs none of IMP's
   interfaces either. */
%include "IMP_bff.structurereader.i"
#endif

%pythoncode %{
def get_module_version():
    return _IMP_bff.get_module_version()

def get_module_name():
    return "IMP::bff"

def get_build():
    """Which IMP.bff this is.

    "core"      -- no IMP at all: the fitting stack, volumes, rotamer dyes.
    "core+imp"  -- the same, plus the roads that run on an atomistic model
                   (attach_dye_to_pdb, run_dye_langevin). IMP is linked as a
                   private library; there is still no IMP in Python.
    "imp"       -- the IMP module build, where IMP's own Python is there too.
    """
    return "core+imp" if _IMP_BFF_WITH_IMP else "core"

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
