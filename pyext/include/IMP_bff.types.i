%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"

%init %{
    import_array();
%}

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

// Numpy arrays
%apply (double* INPLACE_ARRAY1, int DIM1) {(double *input, int n_input)}
%apply (double** ARGOUTVIEW_ARRAY1, int* DIM1) {(double **output, int *n_output)}
%apply (long* INPLACE_ARRAY1, int DIM1) {(long *input, int n_input)}
%apply (long** ARGOUTVIEW_ARRAY1, int* DIM1) {(long **output, int *n_output)}
%apply (unsigned char* INPLACE_ARRAY1, int DIM1) {(unsigned char *input, int n_input)}
%apply (unsigned char** ARGOUTVIEW_ARRAY1, int* DIM1) {(unsigned char **output, int *n_output)}
%apply (float** ARGOUTVIEWM_ARRAY3, int* DIM1, int* DIM2, int* DIM3) {(float **output, int *nx, int *ny, int *nz)}
