
# libmongoc
###########################
FIND_PACKAGE (mongoc-1.0 1.7 REQUIRED)
LINK_LIBRARIES (mongo::mongoc_shared)

# SIMD
###########################
IF(NOT APPLE)
    IF(CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "x86_64")
        SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/dependency)
        INCLUDE(simd)
        IF(${AVX_FOUND})
            MESSAGE("BUILD WITH SIMD")
            add_compile_definitions(WITH_AVX)
            IF (MSVC)
                # https://devblogs.microsoft.com/cppblog/simd-extension-to-c-openmp-in-visual-studio/
                # /Oi is for intrinsics
                SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX -openmp:experimental /Oi")
            ELSE(MSVC)
                STRING(APPEND CMAKE_CXX_FLAGS " -Ofast -mavx -mfma -ffast-math -funroll-loops")
            ENDIF (MSVC)
        ELSE(${AVX_FOUND})
            UNSET(WITH_AVX)
        ENDIF(${AVX_FOUND})
    ENDIF()
ENDIF()
