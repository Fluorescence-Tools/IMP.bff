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
# The version: the wheel build hands it in from the git tag
# (scikit-build-core's SKBUILD_PROJECT_VERSION); otherwise the release the
# conda recipe names. project() wants digits and dots only.
set(_imp_bff_version "0.15.0")
if(DEFINED SKBUILD_PROJECT_VERSION)
  string(REGEX MATCH "^[0-9]+(\\.[0-9]+)*" _imp_bff_version "${SKBUILD_PROJECT_VERSION}")
endif()
project(imp_bff_core VERSION ${_imp_bff_version} LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)

# MSVC's <cmath> hides M_PI and friends behind _USE_MATH_DEFINES, and its
# <windows.h> leaks min/max into every translation unit that includes it.
if(MSVC)
  add_compile_definitions(_USE_MATH_DEFINES NOMINMAX)
  # The unity build (bff_core_all.cpp) exceeds MSVC's object-file section
  # count without this (C1128).
  add_compile_options(/bigobj)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# Default, not hidden, and that is a decision rather than an oversight:
# hidden visibility stops exporting the typeinfo of the exception classes,
# and the `catch (const IMP::ValueException&)` that translates them lives in
# the extension module while the `throw` is in this library -- across that
# boundary the match then fails and every refusal arrives as a bare
# RuntimeError. Measured: it saved 0.23 MB of a 5.5 MB wheel and broke four
# tests. Marking each thrown type visible would buy it back, if the size ever
# matters more than the risk.
set(CMAKE_CXX_VISIBILITY_PRESET default)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(IMPBFF_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/..")
option(IMPBFF_UNITY "Compile the core as one translation unit, as IMP's build does" ON)
# The wheel layout (PRD-137 step 6d): everything under site-packages/IMP/bff/
# -- the extension, the shadow module, libimp_bff beside them (RPATH
# @loader_path / $ORIGIN), and the small data set in IMP/bff/data. The rest
# of data/ is fetched on demand (IMP.bff.fetch_data). scikit-build-core sets
# SKBUILD, so a `pip wheel .` needs no flag; a conda or plain install keeps
# the usual lib/ + include/ + share/IMP/bff layout.
option(IMPBFF_WHEEL "Lay the install out as a Python wheel" "${SKBUILD}")
if(IMPBFF_WHEEL)
  # no compiled-in directory: the wheel's __init__ sets IMP_BFF_DATA
  set(IMPBFF_DATA_DIR "")
  set(IMPBFF_EXAMPLE_DIR "")
endif()
set(IMPBFF_DATA_DIR "${CMAKE_INSTALL_PREFIX}/share/IMP/bff/data" CACHE PATH "Where get_data_path() looks at run time")
set(IMPBFF_EXAMPLE_DIR "${CMAKE_INSTALL_PREFIX}/share/IMP/bff/examples" CACHE PATH "Where get_example_path() looks at run time")

# ---- the connection layer, when a pre-built IMP is there (PRD-139)
#
# Two worlds, one source tree. Without IMP the shims in standalone/include
# stand in for IMP's vocabulary and only the core is built. With
# IMPBFF_WITH_IMP the shims are off the path, IMP's own headers are used, and
# src/imp/ is compiled in -- so the same sources produce the connection layer
# too, including the dye Langevin road. The two cannot be mixed in one
# library: both define IMP::Object, IMP::Exception and IMP::algebra::VectorD,
# so this is a switch, not an addition.
#
# The modules this needs, measured (PRD-139): kernel, algebra, display,
# score_functor, core, container, atom -- 10 MB -- plus em and rotamer for
# two single functions (create_density_map, the Dunbrack side-chain
# placement), 1.4 MB more. Nothing else of IMP is named anywhere in the tree.
option(IMPBFF_WITH_IMP "Build the IMP connection layer against a pre-built IMP" OFF)
set(IMPBFF_IMP_MODULES "kernel;algebra;display;score_functor;core;container;atom;em;statistics;rotamer"
    CACHE STRING "IMP modules to link when IMPBFF_WITH_IMP is on")

# ---- dependencies
# Eigen: no version asked. A "3.3" request makes Eigen 5's config refuse the
# match -- its PACKAGE_VERSION_COMPATIBLE holds only for a same-major
# request -- and conda-forge has moved to 5.0.1, so the version floor here
# became a build breaker. What this library uses of Eigen is stable across
# the 3/5 boundary; a real API break will surface as a compile error, not
# as a silent wrong answer.
find_package(Eigen3 REQUIRED CONFIG)
find_package(Boost 1.66 REQUIRED)            # headers: random (1.66 is what manylinux_2_28 ships)
find_package(cereal CONFIG QUIET)
if(NOT cereal_FOUND)
  find_path(CEREAL_INCLUDE_DIR cereal/cereal.hpp REQUIRED)
endif()
# RMF: one source file and four functions (write_rmf, the two rotamer-library
# doors, protein_frames_from_rmf), called from bin/imp_bff and nowhere in the
# library. Measured (PRD-139): librmf is 3.2 MB, but it links Boost.Iostreams,
# which in some builds pulls Boost.Regex and all of ICU -- 35 MB a wheel would
# carry for a file format most callers never open. On by default; a wheel that
# wants to stay small turns it off.
option(IMPBFF_WITH_RMF "Read and write RMF files (librmf)" ON)
if(IMPBFF_WITH_RMF)
  find_path(RMF_INCLUDE_DIR RMF/FileHandle.h REQUIRED)
  find_library(RMF_LIBRARY NAMES RMF REQUIRED)
endif()
# OpenMP, and the AppleClang hole it usually falls into.
#
# `find_package(OpenMP)` reports success on this platform and hands back an
# *empty* `OpenMP_CXX_FLAGS`: AppleClang does not accept `-fopenmp` directly,
# so CMake fills in the C flags and leaves C++ blank. A `#pragma` the compiler
# does not recognise is a comment, and nothing warns -- which is how every
# speed figure ever recorded for this library came to be single-threaded with
# 19 `#pragma omp` in the sources (okf/validation/openmp_is_not_enabled.md).
#
# So: take what CMake found when it is usable, and otherwise look for libomp
# ourselves and pass `-Xclang -fopenmp`, which is what AppleClang wants.
# `IMP::bff::built_with_openmp()` and `parallel_threads()` report the outcome
# at run time, so a build that silently lost it can be caught from Python.
find_package(OpenMP QUIET)
set(IMPBFF_OPENMP_FLAGS "")
set(IMPBFF_OPENMP_LIBS "")
set(IMPBFF_OPENMP_INCLUDE "")
if(OpenMP_CXX_FOUND AND NOT "${OpenMP_CXX_FLAGS}" STREQUAL "")
  # a string, and `-Xclang -fopenmp` is two arguments: passed whole it
  # reaches the compiler as one unknown option
  separate_arguments(IMPBFF_OPENMP_FLAGS NATIVE_COMMAND "${OpenMP_CXX_FLAGS}")
  set(IMPBFF_OPENMP_LIBS OpenMP::OpenMP_CXX)
  message(STATUS "IMP.bff: OpenMP through CMake (${OpenMP_CXX_FLAGS})")
else()
  find_path(IMPBFF_OMP_INCLUDE_DIR omp.h)
  find_library(IMPBFF_OMP_LIBRARY NAMES omp gomp iomp5)
  if(IMPBFF_OMP_INCLUDE_DIR AND IMPBFF_OMP_LIBRARY)
    set(IMPBFF_OPENMP_FLAGS "-Xclang;-fopenmp")
    set(IMPBFF_OPENMP_LIBS "${IMPBFF_OMP_LIBRARY}")
    set(IMPBFF_OPENMP_INCLUDE "${IMPBFF_OMP_INCLUDE_DIR}")
    message(STATUS "IMP.bff: OpenMP by hand (-Xclang -fopenmp, ${IMPBFF_OMP_LIBRARY})")
  else()
    message(STATUS "IMP.bff: no OpenMP; the kernels run on one thread")
  endif()
endif()

if(IMPBFF_WITH_IMP)
  # Found by hand rather than through IMPConfig.cmake: a wheel builds these
  # libraries itself (utility/wheel_deps.sh) and has no IMP installation,
  # only a prefix -- the same way RMF is found above.
  find_path(IMP_INCLUDE_DIR IMP/Object.h REQUIRED)
  # IMP reads its own data at run time (IMP::atom's top.lib and par.lib for
  # the CHARMM topology, the element table) through IMP_DATA. A wheel has to
  # carry it; an installation beside a real IMP does not.
  find_path(IMP_DATA_DIR IMP/atom/top.lib PATHS "${IMP_INCLUDE_DIR}/../share" PATH_SUFFIXES share)
  if(NOT IMP_DATA_DIR)
    find_path(IMP_DATA_DIR atom/top.lib PATHS "${IMP_INCLUDE_DIR}/../share/IMP" NO_DEFAULT_PATH)
  endif()
  set(IMPBFF_IMP_LIBRARIES "")
  foreach(m ${IMPBFF_IMP_MODULES})
    find_library(IMP_${m}_LIBRARY NAMES imp_${m} REQUIRED)
    list(APPEND IMPBFF_IMP_LIBRARIES "${IMP_${m}_LIBRARY}")
  endforeach()
  message(STATUS "IMP.bff: with the connection layer, against ${IMP_INCLUDE_DIR}")
endif()

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

# ---- the sources: the core, and the layer when IMP is there
file(GLOB _core_sources "${IMPBFF_ROOT}/src/*.cpp")
# ImpLayer.cpp exists only for IMP's own unity build, which compiles no
# subdirectory; here the layer's sources are named directly.
list(FILTER _core_sources EXCLUDE REGEX "/ImpLayer\\.cpp$")
file(GLOB _standalone_sources "${IMPBFF_ROOT}/src/standalone/*.cpp")
if(IMPBFF_WITH_IMP)
  file(GLOB _layer_sources "${IMPBFF_ROOT}/src/imp/*.cpp")
  list(APPEND _core_sources ${_layer_sources})
  # This build has no IMP.bff module inside IMP, so the module's own config
  # header and its data paths are still this project's -- only the shims of
  # IMP's vocabulary are not. The header goes into the generated tree, which
  # is ahead of IMP's on the include path.
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/include/IMP/bff/bff_config.h"
                 "${IMPBFF_GEN_INCLUDE}/IMP/bff/bff_config.h" COPYONLY)
endif()
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

# The compute backend is optional at run time, so it must also be optional at
# build time: OFF both skips the plugin and compiles the loader out of
# Compute.cpp, leaving no dlopen anywhere. IMP's module tooling has no
# per-module options, so that build says the same thing with a compiler flag:
# `cmake .. -DCMAKE_CXX_FLAGS=-DIMPBFF_WITH_GPU=0`.
option(IMPBFF_WITH_GPU "Build and allow the run-time-loaded WGSL compute backend" ON)

add_library(imp_bff SHARED ${_cxx_sources} ${_ihm_sources})
if(NOT IMPBFF_WITH_GPU)
  target_compile_definitions(imp_bff PUBLIC IMPBFF_WITH_GPU=0)
endif()
if(NOT IMPBFF_WITH_RMF)
  target_compile_definitions(imp_bff PUBLIC IMPBFF_NO_RMF=1)
endif()
target_compile_definitions(imp_bff PRIVATE IMPBFF_COMPILATION IMPBFF_EXPORTS
  IMPBFF_STANDALONE_VERSION="${PROJECT_VERSION}"
  IMPBFF_STANDALONE_DATA_DIR="${IMPBFF_DATA_DIR}"
  IMPBFF_STANDALONE_EXAMPLE_DIR="${IMPBFF_EXAMPLE_DIR}")
if(IMPBFF_WITH_IMP)
  # Base.h takes its IMP branch, and IMP's headers are on the path instead of
  # the shims
  target_compile_definitions(imp_bff PUBLIC IMPBFF_WITH_IMP=1)
  target_include_directories(imp_bff
    PUBLIC "$<BUILD_INTERFACE:${IMPBFF_GEN_INCLUDE}>" "${IMP_INCLUDE_DIR}"
           "$<INSTALL_INTERFACE:include>"
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/ihm" ${RMF_INCLUDE_DIR})
  target_link_libraries(imp_bff PUBLIC Eigen3::Eigen Boost::headers
                        PRIVATE ${RMF_LIBRARY} ${IMPBFF_IMP_LIBRARIES})
else()
  target_compile_definitions(imp_bff PUBLIC IMPBFF_STANDALONE=1)
  target_include_directories(imp_bff
    PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"   # the IMP/ shims, first
           "$<BUILD_INTERFACE:${IMPBFF_GEN_INCLUDE}>"
           "$<INSTALL_INTERFACE:include>"
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/ihm" ${RMF_INCLUDE_DIR})
  target_link_libraries(imp_bff PUBLIC Eigen3::Eigen Boost::headers PRIVATE ${RMF_LIBRARY})
endif()
if(cereal_FOUND)
  target_link_libraries(imp_bff PUBLIC cereal::cereal)
else()
  target_include_directories(imp_bff PUBLIC "${CEREAL_INCLUDE_DIR}")
endif()
if(IMPBFF_OPENMP_FLAGS)
  target_compile_options(imp_bff PRIVATE ${IMPBFF_OPENMP_FLAGS})
  target_link_libraries(imp_bff PRIVATE ${IMPBFF_OPENMP_LIBS})
  if(IMPBFF_OPENMP_INCLUDE)
    target_include_directories(imp_bff PRIVATE "${IMPBFF_OPENMP_INCLUDE}")
  endif()
endif()
if(IMPBFF_WHEEL)
  # one file, no version symlinks: a wheel cannot carry symlinks
  set_target_properties(imp_bff PROPERTIES OUTPUT_NAME imp_bff)
else()
  set_target_properties(imp_bff PROPERTIES OUTPUT_NAME imp_bff VERSION ${PROJECT_VERSION} SOVERSION 0)
endif()

# ---- install, in three components: python (the module and what it loads),
# data (what get_data_path() reads), dev (headers, the CMake export).
include(GNUInstallDirs)
if(IMPBFF_WHEEL)
  set(IMPBFF_PY_PKG_DIR "IMP/bff")
  # (COMPONENT binds to the artifact group it follows, so each group names it)
  install(TARGETS imp_bff
    LIBRARY DESTINATION "${IMPBFF_PY_PKG_DIR}" COMPONENT python
    RUNTIME DESTINATION "${IMPBFF_PY_PKG_DIR}" COMPONENT python)
  # the small set travels with the wheel; rotamer_library/ and cgprobe/ (62 MB)
  # are listed in data/registry.json and fetched on first use
  install(DIRECTORY "${IMPBFF_ROOT}/data/" DESTINATION "${IMPBFF_PY_PKG_DIR}/data" COMPONENT data
    PATTERN "rotamer_library" EXCLUDE PATTERN "cgprobe" EXCLUDE)
  if(IMPBFF_WITH_IMP AND IMP_DATA_DIR)
    # score_functor's statistical potentials are 18 of IMP's 21 MB and the dye
    # roads never open them; atom's topology is what they read.
    install(DIRECTORY "${IMP_DATA_DIR}/" DESTINATION "${IMPBFF_PY_PKG_DIR}/imp_data"
            COMPONENT data PATTERN "score_functor" EXCLUDE)
  endif()
else()
  install(TARGETS imp_bff EXPORT imp_bffTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT python
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT dev
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT python)
  install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/" DESTINATION include COMPONENT dev)
  install(DIRECTORY "${IMPBFF_GEN_INCLUDE}/" DESTINATION include COMPONENT dev)
  install(DIRECTORY "${IMPBFF_ROOT}/data/" DESTINATION share/IMP/bff/data COMPONENT data)
  install(DIRECTORY "${IMPBFF_ROOT}/examples/" DESTINATION share/IMP/bff/examples COMPONENT data)
endif()

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
  # Layer headers whose declarations name no IMP type. Where the layer is
  # built they are wrapped like core headers -- that is the whole point of
  # them: a door shaped as arrays and file paths can be offered by a package
  # that carries IMP as a private C++ library and has no IMP in Python
  # (PRD-139). A header added here must keep that property, or SWIG needs
  # IMP's own .i files and the extension imports _IMP_kernel at load.
  set(_layer_array_doors "DyeDynamics.h" "StructureReader.h")
  set(_agg_content "// generated: every public core header\n")
  foreach(h ${_public_headers})
    get_filename_component(n "${h}" NAME)
    list(FIND _layer_list "${n}" _idx)
    if(IMPBFF_WITH_IMP)
      list(FIND _layer_array_doors "${n}" _door_idx)
      if(NOT _door_idx EQUAL -1)
        set(_idx -1)
      endif()
    endif()
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
    "${IMPBFF_ROOT}/pyext" "${IMPBFF_ROOT}/pyext/include" "${CMAKE_CURRENT_SOURCE_DIR}/pyext" "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${IMPBFF_GEN_INCLUDE}" "${CMAKE_BINARY_DIR}" "${Python3_NumPy_INCLUDE_DIRS}")
  set(_swig_defs "")
  if(NOT IMPBFF_WITH_RMF)
    list(APPEND _swig_defs IMPBFF_NO_RMF)
  endif()
  if(IMPBFF_WITH_IMP)
    list(APPEND _swig_defs IMPBFF_WITH_IMP)
  endif()
  if(_swig_defs)
    set_property(SOURCE "${_entry}" PROPERTY COMPILE_DEFINITIONS ${_swig_defs})
  endif()
  set_property(SOURCE "${_entry}" PROPERTY COMPILE_OPTIONS -castmode -naturalvar -fvirtual -Wextra -Dfinal= -DPySwigIterator=IMP_BFF_PySwigIterator -DSwigPyIterator=IMP_BFF_SwigPyIterator)
  swig_add_library(_IMP_bff LANGUAGE python TYPE MODULE
    OUTPUT_DIR "${_pydir}" OUTFILE_DIR "${CMAKE_BINARY_DIR}/swig" SOURCES "${_entry}")
  target_include_directories(_IMP_bff PRIVATE "${CMAKE_BINARY_DIR}" "${Python3_NumPy_INCLUDE_DIRS}" "${RMF_INCLUDE_DIR}")
  target_link_libraries(_IMP_bff PRIVATE imp_bff Python3::Module)
  set_target_properties(_IMP_bff PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${_pydir}" PREFIX "" SUFFIX "${Python3_SOABI_SUFFIX}")
  if(WIN32)
    # Windows Python imports .pyd, never .so (the wheel landed with
    # _IMP_bff.so beside __init__.py and import failed as a circular-import
    # report -- the tell is the extension naming in the install log).
    set_target_properties(_IMP_bff PROPERTIES SUFFIX ".pyd")
  elseif(APPLE)
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
  if(IMPBFF_WHEEL)
    set(_py_install_dir "${IMPBFF_PY_PKG_DIR}")
  else()
    set(IMPBFF_PYTHON_INSTALL_DIR "${Python3_SITEARCH}" CACHE PATH "Where IMP/bff/ is installed")
    set(_py_install_dir "${IMPBFF_PYTHON_INSTALL_DIR}/IMP/bff")
  endif()
  # The GPU backend sits beside the extension, because that is where
  # `IMP.bff.enable_gpu` looks for it. It is a plugin loaded at run time and
  # links no wgpu, so building it costs nothing and requires nothing.
  #
  # It has an off switch because it is *optional at run time* and therefore
  # must not be able to break the build that does not need it. As of
  # 2026-09-08 it has only ever been compiled and run on macOS/Metal; the
  # first Linux and Windows builds will be CI's.
  if(IMPBFF_WITH_GPU)
    set(IMPBFF_GPU_OUTPUT_DIR "${_pydir}")
    set(IMPBFF_GPU_INSTALL_DIR "${_py_install_dir}")
    add_subdirectory("${IMPBFF_ROOT}/gpu" "${CMAKE_BINARY_DIR}/gpu")
  endif()

  install(TARGETS _IMP_bff LIBRARY DESTINATION "${_py_install_dir}" COMPONENT python)
  install(FILES "${_pydir}/__init__.py" DESTINATION "${_py_install_dir}" COMPONENT python)
  if(IMPBFF_WHEEL AND APPLE)
    # install_name_tool rewrites the binaries' load commands at install and
    # leaves them with a stale ad-hoc signature, which arm64 macOS refuses
    # to load (the process is killed). Sign again, in place, after that.
    install(CODE "
      foreach(f \"\${CMAKE_INSTALL_PREFIX}/${_py_install_dir}/_IMP_bff.so\" \"\${CMAKE_INSTALL_PREFIX}/${_py_install_dir}/libimp_bff.dylib\")
        if(EXISTS \"\${f}\")
          execute_process(COMMAND /usr/bin/codesign --force --sign - \"\${f}\" RESULT_VARIABLE _rc)
          if(NOT _rc EQUAL 0)
            message(WARNING \"codesign failed for \${f}\")
          endif()
        endif()
      endforeach()" COMPONENT python)
  endif()
endif()
