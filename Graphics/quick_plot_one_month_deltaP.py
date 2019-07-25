from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import xarray as xr
import pandas as pd
import os

import sys
sys.path.insert(0, '/scratch/mp586/Code/PYCODES')
from plotting_routines_kav7 import *
import stats as st

nc = Dataset('/scratch/mp586/Isca_DATA/ISCA_HPC/full_continents_newbucket_fixedSSTs_zonally_symmetric/run0001/atmos_monthly.nc')
nc2 = Dataset('/scratch/mp586/Isca_DATA/ISCA_HPC/full_continents_newbucket_fixedSSTs_zonally_symmetric_commit7bb4387/run0001/atmos_monthly.nc')

p1 = nc.variables['precipitation'][:]
p2 = nc2.variables['precipitation'][:]
p1 = xr.DataArray(p1)
p2 = xr.DataArray(p2)
p1m = p1.mean('dim_0')
p2m = p2.mean('dim_0')
(p2m - p1m).plot()
plt.show()
