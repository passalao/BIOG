#!/usr/bin/python
# -*- coding: cp1252 -*-
#
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4
import time
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys
sys.path.insert(0, "/home/passalao/Documents/SMOS-FluxGeo/BIOG")
import BIOG

Tz=np.array([-23.15,-23.15])
Thick=200.

Tb=BIOG.fun.GetTb_DMRTML(Tz, Thick, 2, 36.5e9, 50, 16)
print("Tb :", Tb)
