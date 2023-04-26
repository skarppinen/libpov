
# Instructions for usage under Windows 10/11:

To use the library in this repository under Windows 10/11, one needs:

1. The MSVC C compiler (comes with Visual Studio build tools)
2. The GNU scientific library (GSL) (for random number generation)

## Installing dependencies

1. Install Visual Studio Build Tools 
    - This step will install the MSVC C compiler and necessary build tools.
    - During installation, on the left, check "Desktop development with C++" and uncheck everything else. 
    Then, on the right, under "Installation details", only select: 
    - "MSVC v143",
    - "C++ CMake tools for Windows", and
    - "Windows 11 SDK" (also if running Windows 10)

2. Download the archive containing the source code of the GSL library from .., then:
    - Extract the archive somewhere. 
    - Launch "Developer Command Prompt for VS 2022" (installed in the previous step) as administrator and navigate to the extracted folder.
    - In the extracted folder:
      1. Create a new folder named `build` using `mkdir build`.  
      2. Move to the `build` folder, and execute the following three commands (each takes some time):
      3. Run `cmake .. -DGSL_INSTALL_MULTI_CONFIG=ON -DBUILD_SHARED_LIBS=ON -DMSVC_RUNTIME_DYNAMIC=ON -DNO_AMPL_BINDINGS=1`. This creates files from which the library can be compiled.
      4. Run `cmake --build . --config Release` to compile the source code.
      5. Run `cmake --install .`. This installs the compiled library.

## Compiling libpov 

Once the previous steps are done, navigate to the "c" folder in this project and execute:
...
The above command should not produce any errors or warnings and results in the file "libpov.dll" (the library) being created.

## Running the example code in Python 3

To test that everything works, run the script "python/test_pov_computation.py" in Python 3.
