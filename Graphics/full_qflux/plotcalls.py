from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import xarray as xr
import pandas as pd
import os

import sys
sys.path.insert(0, '/scratch/mp586/Code/PYCODES')
import plotting_routines

testdir='full_qflux'
runmin=97 # January
runmax=241 # run 240=Dec (there is no run 241)


landfile=Dataset('/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/input/land.nc',mode='r')
landmask=landfile.variables['land_mask'][:]
lats=landfile.variables['lat'][:]
lons=landfile.variables['lon'][:]

[tsurf,tsurf_avg,tsurf_seasonal_avg,tsurf_month_avg,time]=plotting_routines.seasonal_surface_variable(testdir,runmin,runmax,'t_surf','K')
[convection_rain,convection_rain_avg,convection_rain_seasonal_avg,convection_rain_month_avg,time]=plotting_routines.seasonal_surface_variable(testdir,runmin,runmax,'convection_rain','kg/m^2s') # latent heat flux at surface (UP)
[condensation_rain,condensation_rain_avg,condensation_rain_seasonal_avg,condensation_rain_month_avg,time]=plotting_routines.seasonal_surface_variable(testdir,runmin,runmax,'condensation_rain','kg/m^2s') # latent heat flux at surface (UP)
precipitation_avg=convection_rain_avg+condensation_rain_avg
[net_lhe,net_lhe_avg,net_lhe_seasonal_avg,net_lhe_month_avg,time]=plotting_routines.seasonal_surface_variable(testdir,runmin,runmax,'flux_lhe','W/m^2') # latent heat flux at surface (UP)


PE_avg=precipitation_avg*86400-net_lhe_avg/28. # 28.=conversion from W/m^2 to mm/day using E=H/(rho*L), rho=1000kg/m3, L=2.5*10^6J/kg
# see www.ce.utexas.edu/prof/maidment/CE374KSpr12/.../Latent%20heat%20flux.pptx @30DegC

plotting_routines.globavg_tsurf_timeseries(testdir,runmin,runmax)


land_temp_global=tsurf_avg.where(landmask==1.).mean()
ocean_temp_global=tsurf_avg.where(landmask==0.).mean()
print('Average temperature over land (global) = '+str(land_temp_global))
print('Average temperature over ocean (global) = '+str(ocean_temp_global))

# for specified lats
landmaskxr=xr.DataArray(landmask,coords=[lats,lons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

minlat=input('Enter minimum latitude ')
maxlat=input('Enter maximum latitude ')

land_temp=tsurf_avg.sel(lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==1.) # .where does not recognize xarray as argument!
land_temp.plot()
plt.show()
ocean_temp=tsurf_avg.sel(lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==0.) # .where does not recognize xarray as argument!
ocean_temp.plot()
plt.show()
print('Average temperature over land between '+str(minlat)+' and '+str(maxlat)+' = '+str(land_temp.mean()))
print('Average temperature over ocean between '+str(minlat)+' and '+str(maxlat)+' = '+str(ocean_temp.mean()))



land_precip_global=precipitation_avg.where(landmask==1.).mean()
ocean_precip_global=precipitation_avg.where(landmask==0.).mean()
print('Average precipitation over land (global) = '+str(land_precip_global*86400)+' mm/day')
print('Average precipitation over ocean (global) = '+str(ocean_precip_global*86400)+' mm/day')


land_precip=precipitation_avg.sel(lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==1.) # .where does not recognize xarray as argument!
ocean_precip=precipitation_avg.sel(lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==0.) # .where does not recognize xarray as argument!
print('Average precipitation over land between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(land_precip.mean()*86400)+' mm/day')
print('Average precipitation over ocean between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(ocean_precip.mean()*86400) +' mm/day')
