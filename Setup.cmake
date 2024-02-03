# SIMD
###########################

option(WITH_AVX "Enable AVX support" OFF)

if (WITH_AVX AND NOT APPLE AND CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "x86_64")
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/dependency)
    include(simd)

    if (${AVX_FOUND})
        message("BUILD WITH SIMD")
        add_definitions(-DWITH_AVX)

        if (MSVC)
            # https://devblogs.microsoft.com/cppblog/simd-extension-to-c-openmp-in-visual-studio/
            # /Oi is for intrinsics
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX /Oi")
        else (MSVC)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx -mfma")
        endif (MSVC)
    else (${AVX_FOUND})
        unset(WITH_AVX)
    endif (${AVX_FOUND})
else ()
    # AVX only on x86_64 and not on Apple
    unset(WITH_AVX)
endif()