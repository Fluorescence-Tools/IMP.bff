echo off

:: echo "Patch IMPConfig.cmake"
set IMP_CMAKE="%LIBRARY_LIB%\cmake\IMP\IMPConfig.cmake"
powershell -Command "(gc %IMP_CMAKE%) -replace '.dll', '.lib' | Out-File -encoding ASCII %IMP_CMAKE%"

echo "Build documentation.i"
cd doc
mkdir _build\html\stable\api
doxygen
python "%RECIPE_DIR%\doxy2swig.py" _build/xml/index.xml ../pyext/documentation.i
cd ..

echo "Build the compute backend (optional; absence costs nothing)"

:: A run-time-loaded plugin: one C file linking nothing, found by path at run
:: time, so IMP's module tooling never has to build it. See build.sh for why.
:: Failure is deliberately not fatal -- the CPU path is the default.
set IMPBFF_PY=%PREFIX%\Lib\site-packages\IMP\bff
if "%IMPBFF_WITH_GPU%"=="0" (
  echo "IMP.bff: IMPBFF_WITH_GPU=0; no compute backend built"
) else if exist "%IMPBFF_PY%" (
  cl /nologo /LD /O2 /I "%SRC_DIR%\gpu" "%SRC_DIR%\gpu\imp_bff_wgpu.c" ^
     /Fo"%TEMP%\imp_bff_wgpu.obj" /Fe:"%IMPBFF_PY%\imp_bff_wgpu.dll"
  if errorlevel 1 echo "IMP.bff: no compute backend built; the kernels stay on the CPU"
)

echo "Build app wrapper"

:: build app wrapper
copy "%RECIPE_DIR%\app_wrapper.c" .
cl app_wrapper.c shell32.lib
if errorlevel 1 exit 1

:: The module's swig run %includes RMF.i (through IMP's rmf fragments); the
:: fragments live in RMF's own swig share, not in IMP's, so put them where
:: the swig flags already point. See build.sh for the longer story.
if exist "%PREFIX%\share\RMF\swig" xcopy /E /I /Y "%PREFIX%\share\RMF\swig" "%PREFIX%\share\IMP\swig\" >nul

mkdir build
cd build

echo on
cmake -G Ninja .. ^
      -DCMAKE_BUILD_TYPE=Release ^
      -DCMAKE_PREFIX_PATH="%PREFIX:\=/%;%PREFIX:\=/%\Library" ^
      -DCMAKE_INSTALL_PREFIX="%LIBRARY_PREFIX%" ^
	  -DCMAKE_INSTALL_LIBDIR=bin ^
      -DCMAKE_INSTALL_PYTHONDIR="%SP_DIR%" ^
      -DCMAKE_CXX_FLAGS="/DBOOST_ALL_DYN_LINK /EHsc /DWIN32 /DMSMPI_NO_DEPRECATE_20 /bigobj /DBOOST_ZLIB_BINARY=kernel32"

if errorlevel 1 exit 1

:: Occasionally builds fail on Windows on conda-forge's build hosts
:: due to the compiler running out of heap space. If this happens, try
:: the build again; if it still fails, restrict to one core.
ninja install -k 0
if errorlevel 1 ninja install -k 0
if errorlevel 1 ninja install -k 0 -j 1
if errorlevel 1 exit 1

:: Add wrappers to path for each Python command line tool
:: (all files without an extension)
cd %PREFIX%\Library\bin
for /f %%f in ('dir /b *.') do copy "%SRC_DIR%\app_wrapper.exe" "%PREFIX%\Library\bin\%%f.exe"
if errorlevel 1 exit 1

:: echo "Copy examples"
mkdir -p %PREFIX%\Library\share\doc\IMP\examples\bff
xcopy /e /k /h /i %SRC_DIR%\examples %PREFIX%\Library\share\doc\IMP\examples\bff
:: if errorlevel 1 exit 1
