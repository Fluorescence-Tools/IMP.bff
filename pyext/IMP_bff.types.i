%include "stl.i";
%include "typemaps.i";

%include "std_shared_ptr.i";
%include "std_map.i";
%include "std_set.i";
%include "std_vector.i";
%include "std_pair.i";
%include "std_string.i";
%include "std_list.i";
%include attribute.i

// Member variables of class type come back **by value**.
//
// Without this SWIG hands out a pointer *into* the owning object, and a
// caller who does not keep that object alive reads freed memory:
//
//     IMP.bff.read_forcefield_cif(path).components.values()
//
// -- the system is a temporary, it dies as the expression unwinds, and the
// map proxy that outlives it returns a path of zero bytes on a good day and
// segfaults on a bad one. A copy costs a copy; this is metadata, and the
// arrays that are worth not copying are numpy views (`%attribute_np`) which
// are unaffected.
%naturalvar;

// A `const T&` return is a **borrowed** view of the owner's insides, and SWIG
// hands it to Python as a proxy that does not keep the owner alive:
//
//     IMP.bff.read_forcefield_cif(path).components.values()
//
// -- the system dies as the expression unwinds and the map proxy is left
// pointing at freed memory, which returns a string of zero bytes on a good
// day and segfaults on a bad one. The accessors cannot return by value: they
// are called inside C++ loops (`Scoring.cpp` walks `get_bonds()` per
// iteration), and a copy per iteration is quadratic.
//
// So the copy is made **at the boundary**, where the cost is one copy per
// Python access, and the lifetime becomes Python's.
// The local typedef is not decoration: `SWIG_NewPointerObj` is a C macro, and
// a template argument list carrying a comma would be read as two arguments of
// it. `%arg()` does the same for the commas in the invocations below.
%define %owned_container_out(Type...)
%typemap(out) const Type& {
    typedef Type bff_owned_result_t;
    $result = SWIG_NewPointerObj(new bff_owned_result_t(*$1),
                                 $descriptor(Type*), SWIG_POINTER_OWN);
}
%enddef

// Not `std::vector<std::string>`, `<double>` or `<int>`: those already come
// back converted -- a Python list, a numpy array -- which is a copy, so they
// were never borrowed. Claiming them here replaces that conversion with a raw
// proxy, and a caller comparing a result to a tuple stops recognising it.
%owned_container_out(%arg(std::map<std::string, double>))
%owned_container_out(%arg(std::map<std::string, std::string>))
%owned_container_out(%arg(std::map<std::string, std::vector<std::string> >))

// ...and the ones whose element is a value of this module's own. The list is
// every `const container&` a `ProbeForceFieldSystem` accessor returns; a type
// left off it is not wrong, only still borrowed.
%owned_container_out(%arg(std::map<std::string, IMP::bff::FFComponent>))
%owned_container_out(%arg(std::map<std::string, IMP::bff::FFLJType>))
%owned_container_out(%arg(std::map<std::string, IMP::bff::FFTorsionType>))
%owned_container_out(std::vector<IMP::bff::FFSite>)
%owned_container_out(std::vector<IMP::bff::FFBond>)
%owned_container_out(std::vector<IMP::bff::FFAngle>)
%owned_container_out(std::vector<IMP::bff::FFTorsion>)
%owned_container_out(std::vector<IMP::bff::FFProbe>)

// Pairs
%template(PairFloatFloat) std::pair<float,float>;

// Vectors.
//
// `std::vector` of `std::string`, `double`, `float` and `int` are **not**
// instantiated here, and cannot be: SWIG wraps a type once across a module and
// the modules it imports, and this one imports `rmf` (RmfIO.h), which brings
// `isd`, which brings `saxs` -- whose `%template(DistBase) std::vector<double>`
// is declared first. A `%template(VectorDouble)` here is silently skipped, so
// naming one would promise a class that never appears. `IMP.saxs.DistBase` is
// `std::vector<double>`, `RMF.Strings` is `std::vector<std::string>`.
//
// Nothing in this module's surface needs the *class*: a list or a numpy array
// converts either way, and the two functions that used to make a caller build
// one -- `diffusion_propagate` and `wobbling_kappa2_distribution` -- return
// numpy views and a #IMP::bff::Kappa2Distribution instead.
%template(VectorLong) std::vector<long>;
// One entry per dye, for the multi-dye mean field: a list of numpy arrays on
// the Python side, a vector of flat vectors here.
%template(VectorVectorDouble) std::vector<std::vector<double> >;
%template(VectorVectorString) std::vector<std::vector<std::string> >;

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



/*---------------------------------------------------------------------------*
 * A contiguous ndarray into a `const std::vector<double>&` -- in one memcpy.
 *
 * SWIG's default conversion for a vector parameter walks the sequence and calls
 * SWIG_AsVal_double on every item. Handed a numpy array that is the *slowest*
 * of the three ways in, because each element access mints a fresh Python float:
 *
 *     ndarray -> const std::vector<double>&     32-37 ns/element
 *     list    -> const std::vector<double>&      8-13 ns/element
 *     ndarray -> (double*, int) via IN_ARRAY1        free
 *
 * Measured 2026-08-19 on this box, on `quenched_decay` and
 * `lifetime_spectrum_decay` at 10k-400k elements. Every adapter in
 * A caller that hands a kernel `np.ascontiguousarray(x).ravel()` pays the
 * worst of the three. `diffusion_propagate` takes four ng^3 grids: at
 * ng = 41 that is 9.5 ms of boundary crossing before any work -- the cost of
 * 53 solver steps, and `equilibrium_occupancy` chunks that call up to 200
 * times.
 *
 * This typemap takes the fast path when the argument is a 1-D contiguous
 * float64 array and falls back to SWIG's own converter otherwise, so lists,
 * tuples and VectorDouble keep working exactly as before. It still *copies* --
 * the kernels take a vector and some of them keep it -- but a memcpy is about
 * 0.5 ns/element, so the conversion stops being visible.
 *
 * It is a copy on purpose. Aliasing numpy's buffer into a std::vector is not
 * expressible, and a kernel that outlives the call would hold a dangling
 * pointer. Kernels on the hot path should still take (double*, int) through
 * IN_ARRAY1 and read the caller's memory directly; this typemap is what makes
 * the rest of them cheap without touching their signatures.
 *---------------------------------------------------------------------------*/
%typemap(in, fragment="NumPy_Macros")
        const std::vector<double>& (std::vector<double> imp_bff_tmp,
                                    std::vector<double>* imp_bff_ptr = 0,
                                    PyArrayObject* imp_bff_arr = 0) {
    if (is_numpy_array($input) && array_type($input) == NPY_DOUBLE &&
        array_is_contiguous($input)) {
        const double* imp_bff_data = (const double*) array_data($input);
        const npy_intp imp_bff_n = PyArray_SIZE((PyArrayObject*)$input);
        // Grid kernels take a flat `ng^3` buffer; a caller's cube is C-
        // contiguous, so copying it row-major is exactly a ravel.
        // Any-dimensional, contiguous, and flat out.
        imp_bff_tmp.assign(imp_bff_data, imp_bff_data + imp_bff_n);
        $1 = &imp_bff_tmp;
    } else if ($input == Py_None) {
        // An absent optional grid: `None` means "use the default", and the
        // empty vector is what the kernels read that as.
        $1 = &imp_bff_tmp;
    } else if (SWIG_IsOK(SWIG_ConvertPtr($input, (void**) &imp_bff_ptr,
                                         $descriptor(std::vector<double>*), 0))
               && imp_bff_ptr) {
        // An already-wrapped VectorDouble: pass it through untouched. Taking
        // this case explicitly rather than through swig::asptr matters --
        // overriding the typemap keeps std_vector.i's traits specialisation
        // from ever being emitted, and the generic asptr only knows sequences.
        $1 = imp_bff_ptr;
    } else {
        // Any other array-like -- an int8/uint8 mask, a float32 or an
        // F-contiguous layout -- is converted to a flat float64 copy.
        // swig::asptr does not iterate numpy buffers of the wrong flavour,
        // so it cannot serve here.
        imp_bff_arr = (PyArrayObject*) PyArray_FROMANY(
                $input, NPY_DOUBLE, 0, 0,
                NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_FORCECAST |
                        NPY_ARRAY_ENSUREARRAY);
        if (imp_bff_arr == NULL) {
            SWIG_exception_fail(
                SWIG_ValueError,
                "in method '$defname', argument $argnum of type '$1_type'");
        }
        const npy_intp imp_bff_n = PyArray_SIZE(imp_bff_arr);
        const double* imp_bff_src = (const double*) PyArray_DATA(imp_bff_arr);
        imp_bff_tmp.assign(imp_bff_src, imp_bff_src + imp_bff_n);
        Py_DECREF(imp_bff_arr);
        $1 = &imp_bff_tmp;
    }
}

// The `in` typemap above declares its own temporaries; std_vector.i's default
// `freearg` frees ones it no longer declares, so it has to go with it.
%typemap(freearg) const std::vector<double>& {}

/*---------------------------------------------------------------------------*
 * The same, one level up: a *sequence of* arrays into
 * `const std::vector<std::vector<double> >&`.
 *
 * The multi-dye mean field takes one flat coordinate block per dye, so a
 * caller writes `[rot_1, rot_2]` with each entry an `(n_conf, n_atoms, 3)`
 * array. SWIG's own nested conversion walks the outer sequence and then asks
 * its traits to turn each element into a `std::vector<double>` -- and those
 * traits only know sequences, so a 3-D ndarray arrives as a sequence of 2-D
 * ndarrays and the conversion fails with a type error naming the whole nested
 * type. Flattening at the call site would work; doing it here means the C++
 * signature says what it means and the caller still passes the arrays it
 * has.
 *---------------------------------------------------------------------------*/
%typemap(in, fragment="NumPy_Macros")
        const std::vector<std::vector<double> >&
                (std::vector<std::vector<double> > imp_bff_outer) {
    PyObject* imp_bff_seq = PySequence_Fast(
            $input, "expected a sequence of arrays");
    if (imp_bff_seq == NULL) {
        SWIG_exception_fail(
            SWIG_TypeError,
            "in method '$defname', argument $argnum of type '$1_type'");
    }
    const Py_ssize_t imp_bff_len = PySequence_Fast_GET_SIZE(imp_bff_seq);
    imp_bff_outer.resize((std::size_t) imp_bff_len);
    for (Py_ssize_t imp_bff_i = 0; imp_bff_i < imp_bff_len; ++imp_bff_i) {
        PyObject* imp_bff_item =
                PySequence_Fast_GET_ITEM(imp_bff_seq, imp_bff_i);
        PyArrayObject* imp_bff_arr = (PyArrayObject*) PyArray_FROMANY(
                imp_bff_item, NPY_DOUBLE, 0, 0,
                NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_FORCECAST |
                        NPY_ARRAY_ENSUREARRAY);
        if (imp_bff_arr == NULL) {
            Py_DECREF(imp_bff_seq);
            SWIG_exception_fail(
                SWIG_ValueError,
                "in method '$defname', argument $argnum of type '$1_type'");
        }
        const npy_intp imp_bff_n = PyArray_SIZE(imp_bff_arr);
        const double* imp_bff_src = (const double*) PyArray_DATA(imp_bff_arr);
        imp_bff_outer[(std::size_t) imp_bff_i].assign(imp_bff_src,
                                                      imp_bff_src + imp_bff_n);
        Py_DECREF(imp_bff_arr);
    }
    Py_DECREF(imp_bff_seq);
    $1 = &imp_bff_outer;
}

%typemap(freearg) const std::vector<std::vector<double> >& {}

// ... and back out again as a list of numpy arrays, which is what a caller
// that passed arrays in expects to get. Without this SWIG returns a tuple of
// tuples of floats and every consumer starts with `np.asarray`.
%typemap(out, fragment="NumPy_Macros") std::vector<std::vector<double> > {
    $result = PyList_New((Py_ssize_t) $1.size());
    if ($result == NULL) SWIG_fail;
    for (std::size_t imp_bff_i = 0; imp_bff_i < $1.size(); ++imp_bff_i) {
        npy_intp imp_bff_dim = (npy_intp) $1[imp_bff_i].size();
        PyObject* imp_bff_arr = PyArray_SimpleNew(1, &imp_bff_dim, NPY_DOUBLE);
        if (imp_bff_arr == NULL) {
            Py_DECREF($result);
            SWIG_fail;
        }
        if (imp_bff_dim > 0) {
            std::memcpy(PyArray_DATA((PyArrayObject*) imp_bff_arr),
                        &$1[imp_bff_i][0],
                        (std::size_t) imp_bff_dim * sizeof(double));
        }
        PyList_SET_ITEM($result, (Py_ssize_t) imp_bff_i, imp_bff_arr);
    }
}

// The int sibling: an occupancy/field grid `std::vector<int>`. Anything
// array-like is converted through `PyArray_FromAny` to a flat C-contiguous
// int32 copy, so a uint8 mask, a float grid, or an F-contiguous float32
// layout all land here as
// one flat int vector (zeros where a value does not cast). `None` is the empty
// vector; an already-wrapped `std::vector<int>` passes through untouched.
%typemap(in, fragment="NumPy_Macros")
        const std::vector<int>& (std::vector<int> imp_bff_tmp,
                                 std::vector<int>* imp_bff_ptr = 0) {
    if ($input == Py_None) {
        $1 = &imp_bff_tmp;
    } else if (SWIG_IsOK(SWIG_ConvertPtr($input, (void**) &imp_bff_ptr,
                                         $descriptor(std::vector<int>*), 0))
               && imp_bff_ptr) {
        $1 = imp_bff_ptr;
    } else {
        PyArrayObject* imp_bff_arr_tmp = (PyArrayObject*) PyArray_FROMANY(
                $input, NPY_INT, 0, 0,
                NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_FORCECAST |
                        NPY_ARRAY_ENSUREARRAY);
        if (imp_bff_arr_tmp == NULL) {
            SWIG_exception_fail(
                SWIG_ValueError,
                "in method '$symname', argument $argnum of type '$1_type'");
        }
        const npy_intp imp_bff_n = PyArray_SIZE(imp_bff_arr_tmp);
        const int* imp_bff_src = (const int*) PyArray_DATA(imp_bff_arr_tmp);
        imp_bff_tmp.resize(static_cast<std::size_t>(imp_bff_n));
        for (npy_intp imp_bff_i = 0; imp_bff_i < imp_bff_n; ++imp_bff_i) {
            imp_bff_tmp[static_cast<std::size_t>(imp_bff_i)] =
                    imp_bff_src[imp_bff_i];
        }
        Py_DECREF(imp_bff_arr_tmp);
        $1 = &imp_bff_tmp;
    }
}

%typemap(freearg) const std::vector<int>& {}

%typemap(typecheck, precedence=SWIG_TYPECHECK_VECTOR, fragment="NumPy_Macros")
        const std::vector<int>& {
    void* imp_bff_vp = 0;
    $1 = (is_numpy_array($input) && array_is_contiguous($input)) ? 1
       : ($input == Py_None) ? 1
       : SWIG_IsOK(SWIG_ConvertPtr($input, &imp_bff_vp,
                                   $descriptor(std::vector<int>*), 0)) ? 1
       : (PySequence_Check($input) ? 1 : 0);
}

%typemap(typecheck, precedence=SWIG_TYPECHECK_VECTOR, fragment="NumPy_Macros")
        const std::vector<double>& {
    void* imp_bff_vp = 0;
    $1 = (is_numpy_array($input) && array_type($input) == NPY_DOUBLE &&
          array_is_contiguous($input)) ? 1
       : ($input == Py_None) ? 1
       : SWIG_IsOK(SWIG_ConvertPtr($input, &imp_bff_vp,
                                   $descriptor(std::vector<double>*), 0)) ? 1
       : (PySequence_Check($input) ? 1 : 0);
}

/*---------------------*/
// Generic numpy arrays
/*---------------------*/

// Inplace arrays
/*---------------------*/

// Float/Double
%apply(double* INPLACE_ARRAY1, int DIM1) {(double* inplace_output, int n_output)}
%apply (double* INPLACE_ARRAY1, int DIM1) {(double *input, int n_input)}
%apply (long* INPLACE_ARRAY1, int DIM1) {(long *input, int n_input)}
%apply (unsigned char* INPLACE_ARRAY1, int DIM1) {(unsigned char *input, int n_input)}

// Input array
/*---------------------*/

// Float/Double
%apply(double* IN_ARRAY1, int DIM1) {(double *input, int n_input)}
%apply(double* IN_ARRAY2, int DIM1, DIM2) {(double *input, int n_input1, int n_input2)}

// Integers
%apply(char* IN_ARRAY1, int DIM1) {(char *input, int n_input)}
%apply(short* IN_ARRAY1, int DIM1) {(short* input, int n_input)}
%apply(unsigned short* IN_ARRAY1, int DIM1) {(unsigned short* input, int n_input)}
%apply(int* IN_ARRAY1, int DIM1) {(int* input, int n_input)}

// Zero-copy inputs for the kernels whose arrays are grids. Converting a numpy
// array into a std::vector costs about 34 ns per element -- on a 101^3
// occupancy grid that is 34 ms of pure marshalling per call, which for a short
// walk is the entire wall clock. These typemaps hand the kernel numpy's own
// buffer instead. The parameter names below must match the C++ exactly.
%apply(int* IN_ARRAY1, int DIM1) {(int* occupancy, int n_occupancy)}
%apply(double* IN_ARRAY1, int DIM1) {(double* mobility, int n_mobility)}
%apply(double* IN_ARRAY1, int DIM1) {(double* rate_map, int n_rate_map)}
%apply(double* IN_ARRAY1, int DIM1) {(double* coords, int n_coords)}
%apply(double* IN_ARRAY1, int DIM1) {(double* rotamer_coords, int n_rotamer_coords)}
%apply(double* IN_ARRAY1, int DIM1) {(double* rotamer_weights, int n_rotamer_weights)}
%apply(double* IN_ARRAY1, int DIM1) {(double* protein_coords, int n_protein_coords)}
%apply(double* IN_ARRAY1, int DIM1) {(double* rmin_ij, int n_rmin_ij)}
%apply(double* IN_ARRAY1, int DIM1) {(double* eps_ij, int n_eps_ij)}
%apply(double* IN_ARRAY1, int DIM1) {(double* points1, int n_points1)}
%apply(double* IN_ARRAY1, int DIM1) {(double* points2, int n_points2)}
%apply(double* IN_ARRAY1, int DIM1) {(double* mu1, int n_mu1)}
%apply(double* IN_ARRAY1, int DIM1) {(double* mu2, int n_mu2)}
%apply(double* IN_ARRAY1, int DIM1) {(double* r, int n_r)}
%apply(double* IN_ARRAY1, int DIM1) {(double* kappa2, int n_kappa2)}
%apply(double* IN_ARRAY1, int DIM1) {(double* coords_a, int n_coords_a)}
%apply(double* IN_ARRAY2, int DIM1, int DIM2) {(double* z_values, int n_frames, int n_pair)}
%apply(double* IN_ARRAY1, int DIM1) {(double* values, int n_values)}
%apply(double* IN_ARRAY1, int DIM1) {(double* weights, int n_weights)}
%apply(double* IN_ARRAY1, int DIM1) {(double* coords_b, int n_coords_b)}
%apply(int* IN_ARRAY1, int DIM1) {(int* cluster_centers, int n_cluster_centers)}
// (n_frames, n_atoms, 3): the conformer stack a rotamer library is clustered
// from. 3-D so the kernel knows the frame stride without being told twice.
%apply(double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(double* cluster_coords, int n_frames, int n_atoms, int n_dim)}
%apply(long long* IN_ARRAY1, int DIM1) {(long long *input, int n_input)}
%apply(unsigned long long* IN_ARRAY1, int DIM1) {(unsigned long long *input, int n_input)}

// `std::map<std::string, std::string>` -- an AV's `params`, a distance-type
// table, anything keyed and valued by name. Instantiated here rather than in
// the first `.i` that happens to need it, because SWIG resolves a template at
// the point of use and a second consumer earlier in the include order gets a
// bare `SwigPyObject` instead of a mapping.
%template(MapStringString) std::map<std::string, std::string>;

// Output arrays views
/*---------------------*/
// floating points
%apply(double** ARGOUTVIEW_ARRAY1, int* DIM1) {(double** output_view, int* n_output), (double **output, int *n_output)}
%apply (long** ARGOUTVIEW_ARRAY1, int* DIM1) {(long **output, int *n_output)}
%apply(float** ARGOUTVIEW_ARRAY1, int* DIM1, int* DIM2, int* DIM3, int* DIM4) {(float **output, int *dim1, int *dim2, int *dim3, int *dim4)}
%apply (unsigned char** ARGOUTVIEW_ARRAY1, int* DIM1) {(unsigned char **output, int *n_output)}

// Generic output memory managed arrays
/*---------------------*/

// float and double
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** output, int* n_output)}

// `out_view` rather than `output`, and deliberately. `(double **output, int
// *n_output)` is claimed TWICE above -- once by ARGOUTVIEW (numpy does *not*
// own the buffer) and once by ARGOUTVIEWM (numpy frees it). Which one binds
// depends on the order of these lines, and getting the non-managed one means a
// silent leak of the whole array on every call. A name only the managed
// typemap claims cannot be resolved the wrong way by reordering.
%apply(double** ARGOUTVIEWM_ARRAY1, int* DIM1) {(double** out_view, int* n_out_view)}
%apply(int** ARGOUTVIEWM_ARRAY1, int* DIM1) {(int** out_view_i, int* n_out_view_i)}
%apply(double** ARGOUTVIEWM_ARRAY2, int* DIM1, int* DIM2) {(double** output, int* n_output1, int* n_output2)}
%apply (double** ARGOUTVIEWM_ARRAY3, int* DIM1, int* DIM2, int* DIM3) {(double** output, int* dim1, int* dim2, int* dim3)}
%apply (float** ARGOUTVIEWM_ARRAY3, int* DIM1, int* DIM2, int* DIM3) {(float **output, int *nx, int *ny, int *nz)}
%apply (float** ARGOUTVIEWM_ARRAY4, int* DIM1, int* DIM2, int* DIM3, int* DIM4) {(float **output, int *dim1, int *dim2, int *dim3, int *dim4)}

// integers
%apply(long long** ARGOUTVIEWM_ARRAY1, int* DIM1) {(long long **output, int *n_output)}
%apply(unsigned long long** ARGOUTVIEWM_ARRAY1, int* DIM1) {(unsigned long long** output, int* n_output)}
%apply(int** ARGOUTVIEWM_ARRAY1, int* DIM1) {(int** output, int* n_output)}
%apply(unsigned int** ARGOUTVIEWM_ARRAY1, int* DIM1) {(unsigned int** output, int* n_output)}
%apply(short** ARGOUTVIEWM_ARRAY1, int* DIM1) {(short** output, int* n_output)}
%apply(unsigned short** ARGOUTVIEWM_ARRAY1, int* DIM1) {(unsigned short** output, int* n_output)}
%apply(char** ARGOUTVIEWM_ARRAY1, int* DIM1) {(char** output, int* n_output)}
%apply(signed char** ARGOUTVIEW_ARRAY1, int* DIM1) {(signed char** output, int* n_output)}
%apply (unsigned int** ARGOUTVIEWM_ARRAY2, int* DIM1, int* DIM2) {(unsigned int** output, int* dim1, int* dim2)}
%apply (unsigned int** ARGOUTVIEWM_ARRAY3, int* DIM1, int* DIM2, int* DIM3) {(unsigned int** output, int* dim1, int* dim2, int* dim3)}
%apply (unsigned char** ARGOUTVIEWM_ARRAY4, int* DIM1, int* DIM2, int* DIM3, int* DIM4) {(unsigned char** output, int* dim1, int* dim2, int* dim3, int* dim4)}

