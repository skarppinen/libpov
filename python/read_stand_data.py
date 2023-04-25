from pov_data_structures import StandData
import pandas as pd
import numpy as np

def get_stand_data(path, num):
    """
    Function reads testing data number 3 from the XLSX file in path `path`. 
    The output is a `StandData` object, containing the relevant stand data.
    """
    assert num in (2, 3), "`num` should be in `(2, 3)`"
    if num == 3:
        sheet_names = ("Volume means and SDs", "Demand")
        skiprows = 3
        demand_range = range(1, 4)
        prior_mean_range = range(4, 7)
        prior_sd_range = range(7, 10)
        meas_sd_range = range(10, 19)
    else:
        sheet_names = ("Expected volume and variance", "Demand")
        skiprows = 4
        demand_range = range(1, 4)
        prior_mean_range = range(17, 20)
        prior_sd_range = range(21, 24)
        meas_sd_range = range(25, 34)

    volume_data = pd.read_excel(path, sheet_name = sheet_names[0], skiprows = skiprows, header = None)
    demands = pd.read_excel(path, sheet_name = sheet_names[1], header = None, skiprows = 2)  
   
    demands = np.array(demands[[*demand_range]])
    prior_means = np.array(volume_data[[*prior_mean_range]]) 
    prior_vars = np.array(volume_data[[*prior_sd_range]]) ** 2
    nS, nA = prior_means.shape
    nI = 3
    meas_sds = np.array(volume_data[[*meas_sd_range]]).reshape(nS, nA, nI)
    meas_sds = meas_sds.transpose((2, 0, 1)) # note: dimension order is inventory, stand, assortment 
    meas_vars = meas_sds ** 2
    stda = StandData(muprior = prior_means, varprior = prior_vars,
            varmeas = meas_vars, demands = demands)
    return stda
