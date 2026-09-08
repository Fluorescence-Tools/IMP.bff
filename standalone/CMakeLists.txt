# The standalone (IMP-free) build of IMP.bff's core: libimp_bff without the
# connection layer. PRD-137.
#
#   cmake -S standalone -B build-standalone -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
#   cmake --build build-standalone
#
# What it compiles: every source in src/ except the layer (src/imp/, reached
# by src/ImpLayer.cpp), plus src/standalone/ and the vendored python-ihm
# parser. What it needs: Eigen3, cereal, Boost headers (random), RMF (the
# library and its headers), OpenMP when available. What it does not need:
# IMP. The IMP-module build is the top-level CMakeLists.txt; this file never
# touches it.
cmake_minimum_required(VERSION 3.18)
project(imp_bff_core VERSION 0.13.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_VISIBILITY_PRESET default)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(IMPBFF_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/..")
option(IMPBFF_UNITY "Compile the core as one translation unit, as IMP's build does" ON)
set(IMPBFF_DATA_DIR "${CMAKE_INSTALL_PREFIX}/share/IMP/bff/data" CACHE PATH "Where get_data_path() looks at run time")
set(IMPBFF_EXAMPLE_DIR "${CMAKE_INSTALL_PREFIX}/share/IMP/bff/examples" CACHE PATH "Where get_example_path() looks at run time")

# ---- dependencies
find_package(Eigen3 3.3 REQUIRED CONFIG)
find_package(Boost 1.70 REQUIRED)            # headers: random
find_package(cereal CONFIG QUIET)
if(NOT cereal_FOUND)
  find_path(CEREAL_INCLUDE_DIR cereal/cereal.hpp REQUIRED)
endif()
find_path(RMF_INCLUDE_DIR RMF/FileHandle.h REQUIRED)
find_library(RMF_LIBRARY NAMES RMF REQUIRED)
find_package(OpenMP QUIET)

# ---- the headers, laid out as IMP/bff/... like IMP's build tree
set(IMPBFF_GEN_INCLUDE "${CMAKE_BINARY_DIR}/include")
file(GLOB _public_headers "${IMPBFF_ROOT}/include/*.h")
file(GLOB _internal_headers "${IMPBFF_ROOT}/include/internal/*.h")
foreach(h ${_public_headers})
  get_filename_component(n "${h}" NAME)
  configure_file("${h}" "${IMPBFF_GEN_INCLUDE}/IMP/bff/${n}" COPYONLY)
endforeach()
foreach(h ${_internal_headers})
  get_filename_component(n "${h}" NAME)
  configure_file("${h}" "${IMPBFF_GEN_INCLUDE}/IMP/bff/internal/${n}" COPYONLY)
endforeach()

# ---- the sources: the core, not the layer
file(GLOB _core_sources "${IMPBFF_ROOT}/src/*.cpp")
list(FILTER _core_sources EXCLUDE REGEX "/ImpLayer\\.cpp$")
file(GLOB _standalone_sources "${IMPBFF_ROOT}/src/standalone/*.cpp")
set(_ihm_sources "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/ihm/ihm_format.c")
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/ihm/cmp.c")
  list(APPEND _ihm_sources "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/ihm/cmp.c")
endif()
if(IMPBFF_UNITY)
  set(_unity "${CMAKE_BINARY_DIR}/bff_core_all.cpp")
  set(_content "// generated: the core as one translation unit\n")
  list(SORT _core_sources)
  foreach(s ${_core_sources} ${_standalone_sources})
    string(APPEND _content "#include <${s}>\n")
  endforeach()
  file(WRITE "${_unity}" "${_content}")
  set(_cxx_sources "${_unity}")
else()
  set(_cxx_sources ${_core_sources} ${_standalone_sources})
endif()

add_library(imp_bff SHARED ${_cxx_sources} ${_ihm_sources})
target_compile_definitions(imp_bff PRIVATE IMPBFF_COMPILATION IMPBFF_EXPORTS
  IMPBFF_STANDALONE_VERSION="${PROJECT_VERSION}"
  IMPBFF_STANDALONE_DATA_DIR="${IMPBFF_DATA_DIR}"
  IMPBFF_STANDALONE_EXAMPLE_DIR="${IMPBFF_EXAMPLE_DIR}")
target_compile_definitions(imp_bff PUBLIC IMPBFF_STANDALONE=1)
target_include_directories(imp_bff
  PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"   # the IMP/ shims, first
         "$<BUILD_INTERFACE:${IMPBFF_GEN_INCLUDE}>"
         "$<INSTALL_INTERFACE:include>"
  PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/ihm" "${RMF_INCLUDE_DIR}")
target_link_libraries(imp_bff PUBLIC Eigen3::Eigen Boost::headers PRIVATE "${RMF_LIBRARY}")
if(cereal_FOUND)
  target_link_libraries(imp_bff PUBLIC cereal::cereal)
else()
  target_include_directories(imp_bff PUBLIC "${CEREAL_INCLUDE_DIR}")
endif()
if(OpenMP_CXX_FOUND)
  target_link_libraries(imp_bff PRIVATE OpenMP::OpenMP_CXX)
endif()
set_target_properties(imp_bff PROPERTIES OUTPUT_NAME imp_bff VERSION ${PROJECT_VERSION} SOVERSION 0)

# ---- install: the library, the headers (shims first, then the module's), the data
include(GNUInstallDirs)
install(TARGETS imp_bff EXPORT imp_bffTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/" DESTINATION include)
install(DIRECTORY "${IMPBFF_GEN_INCLUDE}/" DESTINATION include)
install(DIRECTORY "${IMPBFF_ROOT}/data/" DESTINATION share/IMP/bff/data)
install(DIRECTORY "${IMPBFF_ROOT}/examples/" DESTINATION share/IMP/bff/examples)

# ---- the Python module: IMP/bff/{__init__.py,_IMP_bff.so}, importable without IMP
option(IMPBFF_PYTHON "Build the Python module (SWIG, Python3, numpy)" ON)
if(IMPBFF_PYTHON)
  find_package(SWIG 4.0 REQUIRED)
  include(UseSWIG)
  find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module NumPy)
  # every core header, for the wrapper's %{ %} block
  set(_agg "${CMAKE_BINARY_DIR}/bff_core_headers.h")
  file(READ "${IMPBFF_ROOT}/src/imp/Headers.cmake" _hc)
  string(REGEX MATCH "\"([^\"]*)\"" _m "${_hc}")
  string(REPLACE ";" " " _layer_list "${CMAKE_MATCH_1}")
  separate_arguments(_layer_list)
  set(_agg_content "// generated: every public core header\n")
  foreach(h ${_public_headers})
    get_filename_component(n "${h}" NAME)
    list(FIND _layer_list "${n}" _idx)
    if(_idx EQUAL -1 AND NOT n STREQUAL "Base.h")
      string(APPEND _agg_content "#include <IMP/bff/${n}>\n")
    endif()
  endforeach()
  file(WRITE "${_agg}" "${_agg_content}")
  set(_pydir "${CMAKE_BINARY_DIR}/python/IMP/bff")
  file(MAKE_DIRECTORY "${_pydir}")
  set(_entry "${CMAKE_CURRENT_SOURCE_DIR}/pyext/IMP_bff_standalone.i")
  set_property(SOURCE "${_entry}" PROPERTY CPLUSPLUS ON)
  set_property(SOURCE "${_entry}" PROPERTY SWIG_MODULE_NAME IMP_bff)
  set_property(SOURCE "${_entry}" PROPERTY INCLUDE_DIRECTORIES
    "${IMPBFF_ROOT}/pyext" "${CMAKE_CURRENT_SOURCE_DIR}/pyext" "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${IMPBFF_GEN_INCLUDE}" "${CMAKE_BINARY_DIR}" "${Python3_NumPy_INCLUDE_DIRS}")
  set_property(SOURCE "${_entry}" PROPERTY COMPILE_OPTIONS -castmode -naturalvar -fvirtual -Wextra -Dfinal= -DPySwigIterator=IMP_BFF_PySwigIterator -DSwigPyIterator=IMP_BFF_SwigPyIterator)
  swig_add_library(_IMP_bff LANGUAGE python TYPE MODULE
    OUTPUT_DIR "${_pydir}" OUTFILE_DIR "${CMAKE_BINARY_DIR}/swig" SOURCES "${_entry}")
  target_include_directories(_IMP_bff PRIVATE "${CMAKE_BINARY_DIR}" "${Python3_NumPy_INCLUDE_DIRS}" "${RMF_INCLUDE_DIR}")
  target_link_libraries(_IMP_bff PRIVATE imp_bff Python3::Module)
  set_target_properties(_IMP_bff PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${_pydir}" PREFIX "" SUFFIX "${Python3_SOABI_SUFFIX}")
  if(APPLE)
    set_target_properties(_IMP_bff PROPERTIES SUFFIX ".so" INSTALL_RPATH "@loader_path/../../lib;@loader_path" BUILD_RPATH "${CMAKE_BINARY_DIR}")
  else()
    set_target_properties(_IMP_bff PROPERTIES SUFFIX ".so" INSTALL_RPATH "$ORIGIN/../../lib;$ORIGIN" BUILD_RPATH "${CMAKE_BINARY_DIR}")
  endif()
  # the shadow module is the package
  add_custom_command(TARGET _IMP_bff POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy "${_pydir}/IMP_bff.py" "${_pydir}/__init__.py"
    COMMENT "IMP/bff/__init__.py from the SWIG shadow module")
  # In the build tree only: an IMP/__init__.py so that this IMP/ wins over an
  # installed IMP on PYTHONPATH (a regular package beats a namespace one). It
  # is not installed -- an installed core is a namespace portion of IMP/.
  file(WRITE "${CMAKE_BINARY_DIR}/python/IMP/__init__.py" "# build-tree marker: the standalone IMP.bff (PRD-137)\n")
  set(IMPBFF_PYTHON_INSTALL_DIR "${Python3_SITEARCH}" CACHE PATH "Where IMP/bff/ is installed")
  install(TARGETS _IMP_bff LIBRARY DESTINATION "${IMPBFF_PYTHON_INSTALL_DIR}/IMP/bff")
  install(FILES "${_pydir}/__init__.py" DESTINATION "${IMPBFF_PYTHON_INSTALL_DIR}/IMP/bff")
endif()
