# libmongoc
###########################
FIND_PACKAGE (mongoc-1.0 1.7 REQUIRED)
LINK_LIBRARIES (mongo::mongoc_shared)

# SIMD
###########################
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/dependency)
INCLUDE(simd)
IF(${AVX_FOUND} AND NOT APPLE)
    MESSAGE("BUILD WITH SIMD")
    IF (MSVC)
        STRING(APPEND CMAKE_CXX_FLAGS " /arch:AVX")
    ELSE(MSVC)
        STRING(APPEND CMAKE_CXX_FLAGS " -Ofast -mavx -mfma -ffast-math -funroll-loops")
    ENDIF (MSVC)
ENDIF(${AVX_FOUND} AND NOT APPLE)
