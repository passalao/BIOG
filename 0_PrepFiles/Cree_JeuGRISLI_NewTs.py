import NC_Resources as ncr
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib as mpl

#Import GRISLI data
GRISLI = netCDF4.Dataset('../../SourceData/WorkingFiles/GRISLIMappedonSMOS.nc')
Acc = GRISLI.variables['Acc']
T = GRISLI.variables['T']
H = GRISLI.variables['H']

#Import Ts data
GRISLI = netCDF4.Dataset('../../SourceData/WorkingFiles/TbSMOSandTsMODIS.nc')
Ts = GRISLI.variables['TsCrocus']

