FIND_PACKAGE(PythonInterp)
FIND_PACKAGE(PythonLibs)

EXECUTE_PROCESS(
        COMMAND ${PYTHON_EXECUTABLE} -c "from sysconfig import get_paths as gp; print(gp()['include'])"
        OUTPUT_VARIABLE PYTHON_INCLUDE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
)
EXECUTE_PROCESS(
        COMMAND ${PYTHON_EXECUTABLE} -c "from __future__ import print_function; import numpy; print(numpy.get_include())"
        OUTPUT_VARIABLE Python_NumPy_PATH
        OUTPUT_STRIP_TRAILING_WHITESPACE
)


INCLUDE_DIRECTORIES(BEFORE ${Python_NumPy_PATH} ${PYTHON_INCLUDE_DIR})
LINK_LIBRARIES(${PYTHON_LIBRARY})

MESSAGE(STATUS "PYTHON_EXECUTABLE='${PYTHON_EXECUTABLE}'")
MESSAGE(STATUS "Python_NumPy_PATH='${Python_NumPy_PATH}'")
MESSAGE(STATUS "PYTHON_INCLUDE_DIR='${PYTHON_INCLUDE_DIR}'")
MESSAGE(STATUS "PYTHON_LIBRARIES='${PYTHON_LIBRARIES}'")

IF(BUILD_PYTHON_DOCS)
    # Doxygen: To create the python interface documentation
    FIND_PACKAGE(Doxygen REQUIRED)
    MESSAGE(STATUS "Building documentation for Python interface")
    MESSAGE(STATUS "Doxygen executable " ${DOXYGEN_EXECUTABLE})
    #This will be the main output of our command
    execute_process(COMMAND ${DOXYGEN_EXECUTABLE} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc OUTPUT_QUIET)
    execute_process(
            COMMAND python ${CMAKE_SOURCE_DIR}/tools/doxy2swig.py ${CMAKE_SOURCE_DIR}/doc/_build/xml/index.xml ${CMAKE_SOURCE_DIR}/ext/documentation.i
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
ENDIF()
