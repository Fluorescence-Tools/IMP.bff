# SIMD
###########################
IF(NOT APPLE)
    IF(CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "x86_64")
        SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/dependency)
        INCLUDE(simd)
        IF(${AVX_FOUND})
            MESSAGE("BUILD WITH SIMD")
            ADD_DEFINITIONS(-DWITH_AVX)
            IF (MSVC)
                # https://devblogs.microsoft.com/cppblog/simd-extension-to-c-openmp-in-visual-studio/
                # /Oi is for intrinsics
                SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX /Oi")
            ELSE(MSVC)
                SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx -mfma")
            ENDIF (MSVC)
        ELSE(${AVX_FOUND})
            UNSET(WITH_AVX)
        ENDIF(${AVX_FOUND})
    ENDIF()
ENDIF()
