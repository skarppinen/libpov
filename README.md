# libpov 

This repository is a work in progress related to the project "Value of information in harvest planning of Nordic forests".
The repository contains C code that can be compiled to a C library that exposes a function for computing the posterior value of a certain type of forest inventory decision. Example code for calling the function in Python 3 is also included.

## Installing the library under 64bit Windows 10/11

Start by cloning or downloading this repository, and placing it somewhere on the file system.
Next, the following dependencies need to be installed. 

1. The MSVC C/C++ compiler and necessary build tools (comes with Visual Studio build tools)

2. The GNU scientific library (GSL) (for random number generation)

See below for instructions.

### Installing dependencies and the GSL library

1. Install [Build Tools for Visual Studio 2022](https://visualstudio.microsoft.com/downloads/) (scroll down and select "Build Tools for Visual Studio 2022" under "Tools for Visual Studio")
    - During installation, on the left, check "Desktop development with C++" and uncheck everything else. 
    Then, on the right, under "Installation details", only select: 
        1. "MSVC v143",
        2. "C++ CMake tools for Windows"
        3. "Windows 11 SDK" (also in the case that you are running Windows 10)

2. Download the archive containing the source code of the GSL library from [this link](https://github.com/ampl/gsl/archive/refs/heads/master.zip), then:
    - Extract the archive anywhere on the filesystem. 
    - Launch "x64 Native Tools Command Prompt for VS 2022" (installed with the build tools) *as administrator* (this is important) and navigate to the extracted folder using the `cd` command.
    - In the extracted folder:
        1. Create a new folder named "build" using `mkdir build`.  
        2. Move to the "build" folder using `cd build`. 
        3. From the build folder (this is important), execute the script "c/compile_gsl.bat" in this project (on your filesystem) by typing the full path to the file on the command line. 
        This will generate necessary files for compiling and installing GSL using CMAKE (installed with the build tools). 
        After this, the script will install GSL. Executing `compile_gsl.bat" will take a few minutes. 

### Compiling libpov.dll 

Once the previous steps are done, "libpov.dll" (the library) may be compiled. To do this, follow these steps:

1. Find out the folder where GSL was installed. By default, this should be "C:\Program Files\GSL". 
The folder where GSL was installed will be referred to as `GSL_FOLDER` below.

2. Create the environment variables `INCLUDE` and `LIB` and set them to point to the paths of GSL header files and library.
If the environment variables `INCLUDE` and `LIB` already exist, append the respective paths to these variables. 
The GSL headers (`INCLUDE`) are at "`GSL_FOLDER`\include" and the library (`LIB`) is at "GSL_FOLDER\lib\Release". 
The correct path for `INCLUDE` will contain a directory "gsl" full of .h files, and the correct path for `LIB` will contain the file "gsl.lib". 

2. Create another environment variable named `LIBPOV_DLLS` that contains two paths (separated by ';'):

    1. The path "`GSL_FOLDER`\bin\Release". This path should contain "gsl.dll" and "gslcblas.dll". 

    2. The path to the folder "python" in this project (on your filesystem). This path should point to the folder containing "libpov.dll" (compiled shortly, see below). 

2. Close and restart the Native Developer Command prompt for the changes to the environment variables to take effect.

3. Navigate to the "c" folder in this project (on your filesystem) and execute: `compile_libpov.bat` (just by typing the name). This compiles libpov and moves the library "libpov.dll" to the "python" folder in this project (on your filesystem).

### Running example code in Python 3 to see that everything works

Finally, to test that the library has been successfully installed, run the following commands in the "python" folder of this project:

1. Create a Python virtual environment with `python -m venv venv`. Here, `python` corresponds to a Python version of 3.10+, although lower versions are likely to work as well. This step creates a reproducible Python environment.

2. Activate the virtual environment with `venv\Scripts\activate`.

3. Install the required Python packages with `pip install -r requirements.txt`. This will install a few Python packages that are used by the Python files in this project.

4. Place "test_dataset_3.xlsx" to the "python" folder. (see also source code for "test_compute_pov.py") 

5. With the virtual environment activated, run the example script `python test_compute_pov.py`.

If everything is successful, the example script prints out the value of PoV for a particular configuration.
Inspecting the source code of the test script should make it clear how the PoV function should be called from arbitrary code.
Finally, note that the second path in the environment variable `LIBPOV_DLLS` above needs to point to the folder that contains the DLL "libpov.dll". This ensures that calls to the PoV function will work. It is therefore possible to place "libpov.dll" (and any files that call PoV) at arbitrary locations in the file system, provided that the environment variable is updated accordingly. 
