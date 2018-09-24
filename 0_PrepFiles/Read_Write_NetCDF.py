
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import netCDF4
import NC_Resources as ncr

#######################################################################################
# Import NetCDF data
Model = netCDF4.Dataset('../../SourceData/GRISLI/Avec_FoxMaule/T3D-AN40C006-k0.nc')
ny_Mod = Model.dimensions['y'].size
nx_Mod = Model.dimensions['x'].size
Tz=Model.variables['T']
Ts=Tz[0,:,:] #Seulement les noeuds de surface

#Pratique pour extraire les informations du NetCDF (NC_Resources.py)
#J'utilise en général nc_moddims pour ensuite écrire un fichier netCDF au même format
nc_modattrs, nc_moddims, nc_modvars = ncr.ncdump(Model)

#######################################################################################
#Create new NetCDF file
nc_out = Dataset('../../SourceData/GRISLI/Avec_FoxMaule/Ts.nc', 'w', format='NETCDF4')
nc_out.description = "Output file in NetCDF format"

#Create NetCDF dimensions
for dim in nc_moddims:
    dim_name=dim
    if dim=="rows":
        dim_name="x"
    if dim=="cols":
        dim_name="y"
    nc_out.createDimension(dim_name, Model.dimensions[dim_name].size)

#Create variables and fill with data
nc_out.createVariable('Ts', 'float64', ('x','y'))
nc_out.variables['Ts'][:] = Ts
nc_out.close()