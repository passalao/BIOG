#!/usr/bin/python
# -*- coding: cp1252 -*-
#
import NC_Resources as ncr
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib as mpl

#nc_f = '../../SourceData/GRISLI/RGL-0194_an2040.nc'
nc_f = '../../SourceData/GRISLI/TB40S004_1.nc'
nc_fid = Dataset(nc_f, 'r')

nc_fT = '../../SourceData/GRISLI/TB40S004_1_T.nc'
nc_fidT = Dataset(nc_fT, 'r')

#Print out the properties of NetCDF file:
nc_attrs, nc_dims, nc_vars = ncr.ncdump(nc_fid)

#Extract data from NetCDF file
H = nc_fid.variables['H'][:]  # extract/copy the data
B = nc_fid.variables['B'][:]  # extract/copy the data
Uxbar = nc_fid.variables['UXBAR'][:]  # extract/copy the data
Uybar = nc_fid.variables['UYBAR'][:]  # extract/copy the data
Ux = nc_fid.variables['UX'][:]  # extract/copy the data
Uy = nc_fid.variables['UY'][:]  # extract/copy the data
T = nc_fidT.variables['T'][:]  # extract/copy the data
Bmelt = nc_fid.variables['BMELT'][:]  # extract/copy the data
#PhiG = nc_fid.variables['PhiG'][:]  # extract/copy the data

#Get the surface temperature
Ts=T[0,:,:]

#Get the Geothermal flux from Shapiro dataset
nc_gf = '../../SourceData/GRISLI/ghf_Ant_Lebrocq_40km.nc'
nc_gfid = Dataset(nc_gf, 'r')
PhiG = nc_gfid.variables['ghf'][:]
print(np.shape(PhiG))
lx=np.shape(PhiG)[1]
ly=np.shape(PhiG)[2]

#Get the accumulation
nc_acc = '../../SourceData/GRISLI/RACMO2.3p2_ANT27_smb_ltm_1979_2016_ant40.grd'
nc_accid = Dataset(nc_acc, 'r')
Acc = nc_accid.variables['smb'][:]


#Compute Z, elevation above bedrock, and Zeta, reduced height
Z=np.zeros((21,ly,lx))
Zeta=np.zeros((21,ly,lx))

for i in np.arange(0,20):
    Z[i,:,:]=0.05*(20-i)*H[:,:]
    Zeta[i,:,:]=0.05*(20-i)

#Compute hor. divergence, and reshape
GradUx=np.zeros((21,ly,lx))
GradUy=np.zeros((21,ly,lx))
LogUh=np.zeros((21,ly,lx))

j=0  #j for Y coordinates, i for X coordinates
for line in Uxbar[:]:
    i=0
    for node in line:
        k=0
        for u in Ux[:,j,i]:
            #print(k)
            LogUh[k,j,i] = math.log((Ux[k,j,i] ** 2 + Uy[k,j,i] ** 2) ** 0.5 + 1e-6)
            k=k+1

        if i!=0 and i!=lx-1 and j!=0 and j!=ly-1:
            k=0
            for z in GradUx[:,0,0]:
                GradUx[k,j,i]=(Ux[k,j,i+1] - Ux[k,j,i-1])/50000 # cell width of 15km
                GradUy[k,j,i]=(Uy[k,j+1,i] - Uy[k,j-1,i])/50000
                k=k+1
        i=i+1
    j=j+1

HDiv=GradUx+GradUy

'''#Plot the data
fig, ax = plt.subplots()
cmap = mpl.cm.Reds#coolwarm
norm = mpl.colors.Normalize(vmin=0.0, vmax=0.25)
myplot=plt.pcolormesh(Acc, cmap=cmap, norm=norm)
cbar = fig.colorbar(myplot, ticks=[0.0,0.1,0.25])
cbar.ax.set_yticklabels(['0.0', '0.1', '0.25'])  # vertically oriented colorbar
plt.axis('equal')
plt.show()'''

# Export of the enriched GRISLI dataset for KERAS
w_nc_fid = Dataset('../../SourceData/WorkingFiles/TB40S004_1_KERAS.nc', 'w', format='NETCDF4')
w_nc_fid.description = "GRISLI data enriched with accumulation, geothermal flux and horizontal divergence "
for dim in nc_dims:
    w_nc_fid.createDimension(dim, nc_fid.dimensions[dim].size)
nbZelts = nc_fid.dimensions['z'].size

for var in nc_vars:
    w_nc_fid.createVariable(var,nc_fid.variables[var].dtype, \
                            nc_fid.variables[var].dimensions)
    w_nc_fid.variables[var][:] = nc_fid.variables[var][:]

w_nc_fid.createVariable('Acc','float64' ,nc_accid.variables['smb'].dimensions)
w_nc_fid.variables['Acc'][:] = Acc
w_nc_fid.createVariable('PhiG','float64' , nc_gfid.variables['ghf'].dimensions)
w_nc_fid.variables['PhiG'][:] = PhiG
w_nc_fid.createVariable('HDiv','float64' ,nc_fid.variables['TDEP'].dimensions)
w_nc_fid.variables['HDiv'][:] = HDiv
w_nc_fid.createVariable('Zeta','float64' ,nc_fid.variables['TDEP'].dimensions)
w_nc_fid.variables['Zeta'][:] = Zeta
w_nc_fid.createVariable('LogUh','float64' ,nc_fid.variables['TDEP'].dimensions)
w_nc_fid.variables['LogUh'][:] = LogUh
w_nc_fid.close()  # close the new file