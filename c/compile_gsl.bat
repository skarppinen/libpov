cmake .. -G"Visual Studio 17 2022" -A x64 -DGSL_INSTALL_MULTI_CONFIG=ON -DBUILD_SHARED_LIBS=ON -DMSVC_RUNTIME_DYNAMIC=ON -DNO_AMPL_BINDINGS=1
cmake --build . --config Release
cmake --install .