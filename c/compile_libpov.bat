@echo off
cl random-sweep.c model.c /W4 /O2 /c /std:c11
link random-sweep.obj model.obj /DLL /out:libpov.dll gslcblas.lib gsl.lib
copy libpov.dll "../python"
del *.obj libpov.exp libpov.dll libpov.lib std ut