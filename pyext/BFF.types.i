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

