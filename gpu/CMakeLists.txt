# The GPU backend: a plugin, not a dependency.
#
# It links **nothing**. `webgpu.h` is vendored beside this file and every wgpu
# entry point is resolved with dlopen/GetProcAddress at load, from a library
# the Python side finds (wgpu-py ships one for every platform bff targets).
# So this builds on a machine with no wgpu present, the IMP module build never
# sees it, and `dependencies.py` and the conda recipes do not change -- which
# is the constraint PRD-110 defends.
#
# The name matters: `IMP.bff.enable_gpu` looks for libimp_bff_wgpu.dylib,
# libimp_bff_wgpu.so or imp_bff_wgpu.dll beside the extension module.

add_library(imp_bff_wgpu SHARED "${CMAKE_CURRENT_SOURCE_DIR}/imp_bff_wgpu.c")
target_include_directories(imp_bff_wgpu PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}" "${IMPBFF_ROOT}/include")
set_target_properties(imp_bff_wgpu PROPERTIES
  C_STANDARD 11 C_VISIBILITY_PRESET hidden POSITION_INDEPENDENT_CODE ON)

if(NOT WIN32)
  target_link_libraries(imp_bff_wgpu PRIVATE ${CMAKE_DL_LIBS})
endif()
if(APPLE)
  # The plugin is dlopened by libimp_bff; nothing links it, so it needs no
  # install name beyond its own.
  set_target_properties(imp_bff_wgpu PROPERTIES MACOSX_RPATH ON)
endif()

if(IMPBFF_GPU_OUTPUT_DIR)
  set_target_properties(imp_bff_wgpu PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${IMPBFF_GPU_OUTPUT_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${IMPBFF_GPU_OUTPUT_DIR}")
endif()
if(IMPBFF_GPU_INSTALL_DIR)
  install(TARGETS imp_bff_wgpu
          LIBRARY DESTINATION "${IMPBFF_GPU_INSTALL_DIR}" COMPONENT python
          RUNTIME DESTINATION "${IMPBFF_GPU_INSTALL_DIR}" COMPONENT python)
endif()
