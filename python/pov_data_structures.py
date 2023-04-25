import ctypes
from ctypes import c_uint, c_double, c_ulong, POINTER, byref
import numpy as np
import random

def as_flat_double_arr(x: np.ndarray):
    """
    Helper function that converts a numpy array to a column major 1D double array
    that can be used with the ctypes module.
    """
    return (c_double * (x.size))(*x.ravel("F"))

class StandData:
    """
    A class representing information regarding stands.
    `muprior[s, a]` contains the prior mean volume of stand `s` and assortment `a`.
    `varprior[s, a]` contains the prior volume variance of stand `s` and assortment `a`.
    `varmeas[s, a, i]` contains the measurement variance of stand `s` and assortment `a` using \
            inventory method `i`.
    `demands[t, a]` contains the demand of assortment `a` at time `t`. 
    """
    def __init__(self, muprior: np.ndarray, varprior: np.ndarray, 
            varmeas: np.ndarray, demands: np.ndarray):
        assert muprior.shape == varprior.shape, "invalid sizes, muprior and varprior"
        (nS, nA) = muprior.shape
        assert demands.shape[1] == nA, "invalid size, demands.shape[1]"
        assert varmeas.shape[1] == nS, "invalid size, varmeas.shape[1]"
        assert varmeas.shape[2] == nA, "invalid size, varmeas.shape[2]"
        self.nS = nS
        self.nA = nA
        self.nT = demands.shape[0]
        self.nI = varmeas.shape[0]
        self.muprior = muprior
        self.varprior = varprior
        self.varmeas = varmeas
        self.demands = demands
       
class StandDataC(ctypes.Structure):
    """
    C version of the `StandData` struct. This is used to pass data to the C implementation.
    """
    _fields_ = [("nS", c_uint),
                ("nT", c_uint),
                ("nA", c_uint),
                ("nI", c_uint),
                ("muprior", POINTER(c_double)),
                ("sigma2prior", POINTER(c_double)),
                ("sigma2meas", POINTER(c_double)),
                ("demands", POINTER(c_double))
                ]

    def __init__(self, stda: StandData):
        self.nS = c_uint(stda.nS)
        self.nT = c_uint(stda.nT)
        self.nA = c_uint(stda.nA)
        self.nI = c_uint(stda.nI)
        self.muprior = as_flat_double_arr(stda.muprior)
        self.sigma2prior = as_flat_double_arr(stda.varprior)
        self.sigma2meas = as_flat_double_arr(stda.varmeas)
        self.demands = as_flat_double_arr(stda.demands)


class VolumePosteriorC(ctypes.Structure):
    """
    Class used to pass a `VolumePosterior` object to the C code. This class mirrors
    the C struct `VolumePosterior`.
    """
    _fields_ = [("nS", c_uint),
                ("nA", c_uint),
                ("muplus", POINTER(c_double)),
                ("sigma2plus", POINTER(c_double)),
                ("y", POINTER(c_double))
                ]

    def __init__(self, nS: int, nA: int):
        self.nS = c_uint(nS)
        self.nA = c_uint(nA)
        nelems = nS * nA
        self.muplus = (c_double * nelems)(*np.zeros(nelems));
        self.sigma2plus = (c_double * nelems)(*np.zeros(nelems));
        self.y = (c_double * nelems)(*np.zeros(nelems));

class RandomSweepStorageC(ctypes.Structure):
    """
    Class used to pass a `RandomSweepStorage` object to the C code. This class mirrors
    the C struct `RandomSweepStorage` and is used in the random sweep algorithm.
    """
    _fields_ = [("nS", c_uint),
                ("nT", c_uint),
                ("X", POINTER(c_double)),
                ("Xopt", POINTER(c_double)),
                ("C", POINTER(c_double)),
                ("Qb", POINTER(c_double)),
                ("QbX", POINTER(c_double)),
                ("Qbx", POINTER(c_double)),
                ("order", POINTER(c_uint)),
                ("status_lu", POINTER(c_uint)),
                ("increments", POINTER(c_double)),
                ("r", c_double)
                ]
    def __init__(self, stda: StandData):
        self.nS = c_uint(stda.nS)
        self.nT = c_uint(stda.nT)
        nelems = stda.nS * stda.nT
        self.X = (c_double * nelems)(*np.zeros(nelems))
        self.Xopt = (c_double * nelems)(*np.zeros(nelems))
        self.C = (c_double * nelems)(*np.zeros(nelems))
        nelems = stda.nS * stda.nS
        self.Qb = (c_double * nelems)(*np.zeros(nelems))
        nelems = stda.nS * stda.nT
        self.QbX = (c_double * nelems)(*np.zeros(nelems)) 
        nelems = stda.nS
        self.Qbx = (c_double * nelems)(*np.zeros(nelems)) 
        self.order = (c_uint * nelems)(*np.zeros(nelems, dtype = "int32"))
        self.status_lu = (c_uint * nelems)(*np.zeros(nelems, dtype = "int32"))
        nelems = stda.nT + 1
        self.increments = (c_double * nelems)(*np.zeros(nelems))
        self.r = c_double(np.sum(stda.demands ** 2))

def rand_uint(rng):
    bits = ctypes.sizeof(c_ulong(0)) * 8 
    return rng.randint(0, 2 ** bits - 1) 

class SeedVectorC(ctypes.Structure):
    _fields_ = [("arr", POINTER(c_ulong)),
                ("n", c_uint)]

    def __init__(self, n: int, seed: int):
        assert n > 0, "`n` should be positive"
        rng = random.Random(seed)
        seeds = [rand_uint(rng) for _ in range(n)] 
        self.arr = (c_ulong * n)(*seeds)
        self.n = n
         
class PoVStorage:
    """
    A class wrapping all memory used in computing the C function `PoV`.
    The class method `PoV` of this class performs the call to the underlying C code.
    """
    def _setup_PoV(self):
        """
        Method setups up the argument types and return value of the function `PoV` as needed by
        the `ctypes` module. This function need not be called by the user.
        """
        self.lib.PoV.argtypes = [POINTER(c_uint), POINTER(RandomSweepStorageC),
                                 POINTER(VolumePosteriorC), POINTER(StandDataC),
                                 POINTER(SeedVectorC), c_uint, c_uint]
        self.lib.PoV.restype = c_double

    def __init__(self, libname: str, stda: StandData, nsamples: int, seed: int): 
        """
        Initialise a `PoVStorage` object. 

        The input arguments are as follows:
        `libname`: The name of the shared library containing the C function `PoV`. 
        Note that the shared library must be discoverable by the operating system or an error will be thrown. 
        `stda`: A `StandData` object that contains the relevant stand data. 
        See `help(pov_data_structures.StandData)` for more details.
        `nsamples`: The number of Monte Carlo samples used each subsequent call to the method `PoV`.
        `seed`: A random seed used to initialise `nsamples` random seeds used to draw each Monte Carlo sample internally. 
        """
        self.stdaC = StandDataC(stda)
        self.vpC = VolumePosteriorC(nS = stda.nS, nA = stda.nA)
        self.rsC = RandomSweepStorageC(stda)
        self.svC = SeedVectorC(n = nsamples, seed = seed)
        self.lib = ctypes.CDLL(libname)
        self._setup_PoV()
        
    def PoV(self, xI, ninits: int, maxsweeps:int = 9999):
        """
        Main method. Calls the C code for PoV with the following arguments:

        `xI`: the inventory decision, that is, a vector of length `nS` (the number of stands), 
        whose `i`th entry gives the inventory method applied to stand `i` (in 0-based indexing!). 
        Each element in this vector should be within the closed interval `[0, self.stdaC.nI - 1]`, where `self.stdaC.nI` corresponds to the number of possible inventory methods. The input `xI` is not checked when the function is called. 
        `ninits`: the number of random initialisations in the random sweep algorithm.
        `maxsweeps`: the number of maximum sweeps to use in the random sweep algorithm. The default,
        9999 basically runs the method until convergence. The default is recommended.
        """

        assert ninits > 0, "`ninits` should be > 0" 
        assert maxsweeps > 0, "`maxsweeps` should be > 0" 
        xIC = (c_uint * self.stdaC.nS)(*xI)
        return self.lib.PoV(xIC, byref(self.rsC), byref(self.vpC), byref(self.stdaC),
                            byref(self.svC), c_uint(maxsweeps), c_uint(ninits))
