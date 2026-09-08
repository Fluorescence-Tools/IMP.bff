@echo off
:: The core (package `bff`) from standalone\CMakeLists.txt on Windows: the
:: same build as build_core.sh, with conda's Library layout.

:: The large data (data\registry.json) is not in git; the package ships it.
python "%SRC_DIR%\utility\data_registry.py" --fetch "%SRC_DIR%\data" --quiet
if errorlevel 1 exit 1

mkdir build-core
cd build-core

cmake ..\standalone -G Ninja ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DCMAKE_PREFIX_PATH="%PREFIX:\=/%;%PREFIX:\=/%/Library" ^
  -DCMAKE_INSTALL_PREFIX="%LIBRARY_PREFIX%" ^
  -DCMAKE_INSTALL_LIBDIR=bin ^
  -DPython3_EXECUTABLE="%PYTHON%" ^
  -DIMPBFF_PYTHON_INSTALL_DIR="%SP_DIR%"
if errorlevel 1 exit 1

:: See bld.bat: conda-forge's Windows hosts occasionally run the compiler out
:: of heap space, and a retry (finally single-threaded) gets through.
ninja
if errorlevel 1 ninja
if errorlevel 1 ninja -j 1
if errorlevel 1 exit 1

ninja install
if errorlevel 1 exit 1
