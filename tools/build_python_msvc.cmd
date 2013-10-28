cd igraph-*-msvc
:: used so that the python does not try to setup again
set DISTUTILS_USE_SDK=1 
:: use the right environment per architecture
:: VS90COMNTOOLS must probably changed for python 3.3, as another version is 
:: used to compile python
if %PROCESSOR_ARCHITECTURE%==x86 (
  call "%VS90COMNTOOLS%..\..\VC\bin\vcvars32.bat"
) else (
  call "%VS90COMNTOOLS%..\..\VC\bin\vcvars64.bat"
)
vcbuild /upgrade
vcbuild igraph.vcproj "release|x64"
cd ..
cd interfaces\python
python setup.py bdist_wininst
::exit