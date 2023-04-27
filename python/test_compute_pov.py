# This example script shows how the PoV function can be called from Python code.
# The following two imports are necessary.
from pov_data_structures import StandData, PoVStorage
from read_stand_data import get_stand_data

import numpy as np # Used for constructing inventory decision vector.

# Load stand data from Excel file.  
# `get_stand_data` is a quick and dirty function to load data from one of the spreadsheets.
# The function builds a `StandData` object that contains all relevant data from the spreadsheet that
# are needed in the PoV computations. 
# The first argument selects the XLSX file to read from, and the second argument 
# should be set to the index of the data (3 is the latest, 2 is one before that).
# Index 1 also exists but is not implemented currently.
# The index is used in 
data_index = 3
path_to_data = f"test-data-{data_index}.xlsx"
stda = get_stand_data(path_to_data, data_index)

# Build an object containing all relevant objects and data needed to call the PoV function.
# Internally, a handle to the shared library containing the relevant C code is constructed.
# This requires that relevant shared libraries, most notably, libpov.so or libpov.dll can be found by the OS.
# An error is thrown otherwise.
# `stda` is the stand data object constructed above,
# `nsamples` sets the number of Monte Carlo samples used in computing PoV (below), and
# `seed` sets a global random seed that is used to draw a random seed for each of the
# Monte Carlo samples drawn in the PoV computation.
nsamples = 30
povst = PoVStorage(stda = stda, nsamples = nsamples, seed = 1)

# The call to the C function `PoV` is implemented as a class method of the `PoVStorage` object.
# The argument `ninits` specifies the number of random initialisations in the random sweeping algorithm
# used internally.
# The argument `xI` is the inventory decision vector.
ninits = 250
inventory_method_for_all = 3 # Just as an example setting same inventory method for all stands here.
xI = (inventory_method_for_all - 1) * np.ones(stda.nS, dtype = "int32") # -1 is because of 0-based indexing!
pov_value = povst.PoV(xI = xI, ninits = ninits) # Computes PoV (number of Monte Carlo samples is set upon construction of `povst`).
 
print(f"Value of PoV is {pov_value} (using {nsamples} Monte Carlo samples and {ninits} random initialisations)")


