set (CMAKE_CXX_STANDARD 14)
# Are we running cmake from this directory (out of tree build) ?
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  cmake_minimum_required(VERSION 2.8.12)
  project(imp_module)

  # Default Python3_FIND_STRATEGY to LOCATION
  if(POLICY CMP0094)
    cmake_policy(SET CMP0094 NEW)
  endif(POLICY CMP0094)
  cmake_policy(SET CMP0058 NEW)

  # Use new-style variables/escapes
  if(POLICY CMP0053)
    cmake_policy(SET CMP0053 NEW)
  endif(POLICY CMP0053)

  # Quiet warning from cmake 3.1 or later about quoted variables in if()
  if(POLICY CMP0054)
    cmake_policy(SET CMP0054 NEW)
  endif(POLICY CMP0054)

  # Honor link flags in try_compile() tests
  if(POLICY CMP0056)
    cmake_policy(SET CMP0056 NEW)
  endif(POLICY CMP0056)

  SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/tools)
  SET(SWIG_EXECUTABLE swig CACHE STRING "Swig program")

  # IMP
  ###########################
  FIND_PACKAGE(IMP REQUIRED)
  MESSAGE ("-- IMP_USE_FILE \"${IMP_USE_FILE}\"")
  MESSAGE ("-- IMP_INCLUDE_DIR \"${IMP_INCLUDE_DIR}\"")
  include(${IMP_USE_FILE})

  enable_testing()

  set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

  # Don't set MACOSX_RPATH by default
  set(CMAKE_MACOSX_RPATH 0)

  include_directories(SYSTEM ${IMP_INCLUDE_DIR})
  include_directories("${CMAKE_BINARY_DIR}/include")

  imp_build_module(${CMAKE_SOURCE_DIR} bff)

else()
  include(ModuleBuild.cmake)
endif()

include(Setup.cmake)
