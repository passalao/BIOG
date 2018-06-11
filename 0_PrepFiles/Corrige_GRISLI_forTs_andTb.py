from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import NC_Resources as ncr
import netCDF4, math
import numpy as np

# Import Crocus data
ObsCrocus = netCDF4.Dataset('../../SourceData/WorkingFiles/TbSMOSandTsCrocus.nc')
Ts = ObsCrocus.variables['TsCrocus'][::-1,:]
TsCrocus=Ts-273.15
TsCrocus[TsCrocus==-273.15]=0.0

# Import temperature data
GRISLI = netCDF4.Dataset('../../SourceData/WorkingFiles/TB40S123_1_MappedonSMOS.nc')
#Print out the properties of NetCDF file:
nc_attrs, nc_dims, nc_vars = ncr.ncdump(GRISLI)
H = np.array(GRISLI.variables['H'])
S = GRISLI.variables['S']
Tz_gr = GRISLI.variables['T']
Zeta = GRISLI.variables['Zeta']
TsRACMO=Tz_gr[:,:,0]
DeltaTs=TsCrocus-TsRACMO

Tz_corr=np.zeros((201,225,21))
Tz_corr2=np.zeros((201,225,21))

for k in np.arange(0,21):
    Tz_corr[:,:,k]=Tz_gr[:,:,k]+DeltaTs*Zeta[0,0,k]
Ts=Tz_corr[:,:,0]

#Correction so that Tb=Tm
#1 Compute Tm
rho=917
g=9.81
Tm=-(0.074*rho*g*H/1e6+0.024)
DeltaTm=Tm-Tz_corr[:,:,20]

print(Tm[120,120])

#2 Spread the correction
for k in np.arange(0,21):
    Tz_corr2[:,:,k]=Tz_corr[:,:,k]+DeltaTm*(1-Zeta[0,0,k])
Error=Tz_corr2[:,:,20]-Tm
print(np.size(Error[Error<0]))

# Export of the enriched GRISLI dataset for KERAS
w_nc_fid = Dataset('../../SourceData/WorkingFiles/TB40S123_1_Corrected4TsandTb.nc', 'w', format='NETCDF4')
w_nc_fid.description = "GRISLI data enriched with accumulation, geothermal flux and horizontal divergence "
for dim in nc_dims:
    w_nc_fid.createDimension(dim, GRISLI.dimensions[dim].size)
nbZelts = GRISLI.dimensions['z'].size

for var in nc_vars:
    w_nc_fid.createVariable(var,GRISLI.variables[var].dtype, \
                            GRISLI.variables[var].dimensions)
    if var!="T":
        w_nc_fid.variables[var][:] = GRISLI.variables[var][:]
    else:
        w_nc_fid.variables[var][:] = np.array(Tz_corr2[:,:,:])
w_nc_fid.close()  # close the new file