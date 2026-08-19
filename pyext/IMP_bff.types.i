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

// Pairs
%template(PairFloatFloat) std::pair<float,float>;

// Vectors
%template(VectorString) std::vector<std::string>;
%template(VectorDouble) std::vector<double>;
%template(VectorFloat) std::vector<float>;
%template(VectorInt) std::vector<int>;
%template(VectorLong) std::vector<long>;

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
 * `pyext/src` hands these kernels `np.ascontiguousarray(x).ravel()`, so the
 * package was paying the worst of the three everywhere. `diffusion_propagate`
 * takes four ng^3 grids: at ng = 41 one call spent 9.5 ms crossing the boundary
 * before doing any work -- the cost of 53 solver steps, and
 * `equilibrium_occupancy` chunks that call up to 200 times.
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
                                    int imp_bff_res = 0) {
    if (is_numpy_array($input) && array_type($input) == NPY_DOUBLE &&
        array_numdims($input) == 1 && array_is_contiguous($input)) {
        const double* imp_bff_data = (const double*) array_data($input);
        imp_bff_tmp.assign(imp_bff_data,
                           imp_bff_data + array_size($input, 0));
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
        imp_bff_res = swig::asptr($input, &imp_bff_ptr);
        if (!SWIG_IsOK(imp_bff_res) || !imp_bff_ptr) {
            SWIG_exception_fail(
                SWIG_ArgError(imp_bff_res),
                "in method '$symname', argument $argnum of type '$1_type'");
        }
        imp_bff_tmp = *imp_bff_ptr;
        if (SWIG_IsNewObj(imp_bff_res)) delete imp_bff_ptr;
        $1 = &imp_bff_tmp;
    }
}

// The `in` typemap above declares its own temporaries; std_vector.i's default
// `freearg` frees ones it no longer declares, so it has to go with it.
%typemap(freearg) const std::vector<double>& {}

%typemap(typecheck, precedence=SWIG_TYPECHECK_VECTOR, fragment="NumPy_Macros")
        const std::vector<double>& {
    void* imp_bff_vp = 0;
    $1 = (is_numpy_array($input) && array_type($input) == NPY_DOUBLE &&
          array_numdims($input) == 1) ? 1
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
%apply(double* IN_ARRAY1, int DIM1) {(double* coords_b, int n_coords_b)}
%apply(int* IN_ARRAY1, int DIM1) {(int* cluster_centers, int n_cluster_centers)}
// (n_frames, n_atoms, 3): the conformer stack a rotamer library is clustered
// from. 3-D so the kernel knows the frame stride without being told twice.
%apply(double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(double* cluster_coords, int n_frames, int n_atoms, int n_dim)}
%apply(long long* IN_ARRAY1, int DIM1) {(long long *input, int n_input)}
%apply(unsigned long long* IN_ARRAY1, int DIM1) {(unsigned long long *input, int n_input)}

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

