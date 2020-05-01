from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import xarray as xr
import pandas as pd
import os
from matplotlib.patches import Rectangle

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import metpy as metpy

import sys
sys.path.insert(0, '/scratch/mp586/Code/PYCODES')
from plotting_routines_kav7_tabreplaced_ongv3_py3 import *
import stats as st

GFDL_BASE = os.environ['GFDL_BASE']
sys.path.insert(0, os.path.join(GFDL_BASE,'src/extra/python/scripts'))
import cell_area as ca

ctl_model = 'isca'
if (ctl_model == 'Isca') or (ctl_model == 'isca'): 
    control_model = 'Isca_DATA'
elif (ctl_model == 'gfdl') or (ctl_model == 'GFDL'):
    control_model = 'GFDL_DATA'

HPC = 'yes'
control_dir = 'withtv/aquaplanet_frierson_insolation_0qflux_mld20_commitd15c267'
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    control_dir= control_model + '/ISCA_HPC/' + control_dir
else: 
    control_dir= control_model + '/' + control_dir

#print control_dir
ctl_runmin=25
ctl_runmax=97

model = 'isca'
if (model == 'Isca') or (model == 'isca'): 
    model_data = 'Isca_DATA'
    output_dir1 = 'Isca'
elif (model == 'gfdl') or (model == 'GFDL'):
    model_data = 'GFDL_DATA'
    output_dir1 = ''

HPC = 'yes'
testdir_in1= 'withtv/square_South_America_frierson_insolation_lepref1_0qflux_samealbedo_to_01land_samehcp_landocean_commitd15c267' # this needs to be the dark patch simulation
runmin=24
runmax=96
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    exp1_name = 'ISCA_HPC_'+testdir_in1
    testdir = model_data + '/ISCA_HPC/' + testdir_in1
    testdir_in1 = '/ISCA_HPC/' + testdir_in1
else: 
    exp1_name = testdir_in1
    testdir = model_data + '/' + testdir_in1



outdir1 = 'Isca'+testdir_in1
data_pert1 = xr.open_mfdataset('/scratch/mp586/'+testdir+'/*/atmos_monthly_interp.nc')


model = 'isca'
if (model == 'Isca') or (model == 'isca'): 
    model_data = 'Isca_DATA'
    output_dir1 = 'Isca'
elif (model == 'gfdl') or (model == 'GFDL'):
    model_data = 'GFDL_DATA'
    output_dir1 = ''

HPC = 'yes'
testdir_in1= 'withtv/aquaplanet_frierson_insolation_0qflux_mld20_commitd15c267' # this needs to be the dark patch simulation
runmin=24
runmax=96
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    exp1_name = 'ISCA_HPC_'+testdir_in1
    testdir = model_data + '/ISCA_HPC/' + testdir_in1
    testdir_in1 = '/ISCA_HPC/' + testdir_in1
else: 
    exp1_name = testdir_in1
    testdir = model_data + '/' + testdir_in1

outdir2 = 'Isca'+testdir_in1
data_pert2 = xr.open_mfdataset('/scratch/mp586/'+testdir+'/*/atmos_monthly_interp.nc')

outdir = outdir1

area_array, dx, dy = ca.cell_area_all(t_res=42,base_dir='/scratch/mp586/Isca/') # added _all because then dx and dy are also returned 
area_array = xr.DataArray(area_array) # returned in units of m bzw m^2, because radius in cell_area.py is given in metres

area_array_3D = np.expand_dims(area_array, axis=0)
area_array_3D = np.repeat(area_array_3D, 40, axis = 0) # to make area_array 3D (pressure, lat, lon)



data_pert1 = data_pert1.metpy.parse_cf()
data_pert2 = data_pert2.metpy.parse_cf()
data_pert = data_pert1


data_pert_avg1 = data_pert1.mean('time')
data_pert_avg2= data_pert2.mean('time')

data_pert_avg = data_pert_avg1

pfull = data_pert.pfull

data = xr.open_dataset('/scratch/mp586/'+testdir+'/run0121/atmos_monthly.nc')
data = data.metpy.parse_cf()
dp = data.phalf.diff('phalf')*100
dp = xr.DataArray(dp[::-1], coords = [pfull.values], dims = ['pfull'])




dxmp, dymp = metpy.calc.grid_deltas_from_dataarray(data_pert.ucomp)



#conversion following https://www.ncl.ucar.edu/Document/Functions/Contributed/omega_to_w.shtml


cp_dry = 287.04/(2./7.) # units = J/kg/K, https://github.com/ExeClim/Isca/blob/77a3d49c5e3131dc6312d32b3698feac2cc8d156/postprocessing/plevel_interpolation/src/shared/constants/constants.F90
g = 9.81
L = 2.500e6 #J/kg, https://github.com/ExeClim/Isca/blob/77a3d49c5e3131dc6312d32b3698feac2cc8d156/postprocessing/plevel_interpolation/src/shared/constants/constants.F90

vcomp_temp_avg = (data_pert_avg.vcomp_temp)
sensible_v = np.sum((cp_dry/g)*vcomp_temp_avg*dp, axis = 0)
ucomp_temp_avg = (data_pert_avg.ucomp_temp)
sensible_u = np.sum((cp_dry/g)*ucomp_temp_avg*dp, axis = 0)
div_sens = metpy.calc.divergence(sensible_u,sensible_v,dxmp[0,0,:,:],dymp[0,0,:,:])

#same result 
#div_sens = metpy.calc.first_derivative(sensible_u, 'lon') + metpy.calc.first_derivative(sensible_v, 'lat')


# same as div_sens. Integral of div = div of integral because int is over dp, div is in x/y?
# div_sens_divin = np.empty_like(ucomp_temp_avg)
# for i in range(len(pfull)):
#     div_sens_divin[i] = metpy.calc.divergence(ucomp_temp_avg[i,:,:], vcomp_temp_avg[i,:,:], dxmp[0,0,:,:], dymp[0,0,:,:])
# div_sens_divin = xr.DataArray(div_sens_divin, coords = [pfull,data_pert.lat,data_pert.lon], dims = ['pfull','lat','lon'])
# div_sens_divint = np.sum((cp_dry/g)*div_sens_divin*dp, axis = 0)

vcomp_height_avg = (data_pert_avg.vcomp_height)
height_v = np.sum(vcomp_height_avg*dp, axis = 0)
ucomp_height_avg = (data_pert_avg.ucomp_height)
height_u = np.sum(ucomp_height_avg*dp, axis = 0)
div_height = metpy.calc.divergence(height_u, height_v, dxmp[0,0,:,:], dymp[0,0,:,:])


vcomp_latent_avg = (data_pert_avg.sphum_v)
latent_v = np.sum((L/g)*vcomp_latent_avg*dp, axis = 0)
ucomp_latent_avg = (data_pert_avg.sphum_u)
latent_u = np.sum((L/g)*ucomp_latent_avg*dp, axis = 0)
div_latent = metpy.calc.divergence(latent_u, latent_v, dxmp[0,0,:,:], dymp[0,0,:,:])

div_tot = div_latent + div_sens + div_height

div_tot = xr.DataArray(div_tot, coords = [data_pert.lat,data_pert.lon], dims = ['lat','lon'])
div_sens = xr.DataArray(div_sens, coords = [data_pert.lat,data_pert.lon], dims = ['lat','lon'])
div_latent = xr.DataArray(div_latent, coords = [data_pert.lat,data_pert.lon], dims = ['lat','lon'])
div_height = xr.DataArray(div_height, coords = [data_pert.lat,data_pert.lon], dims = ['lat','lon'])


vcomp_avg = (data_pert_avg.vcomp)
move_v = np.sum((1./g)*vcomp_avg*dp, axis = 0)
ucomp_avg = (data_pert_avg.ucomp)
move_u = np.sum((1./g)*ucomp_avg*dp, axis = 0)
div_move = metpy.calc.divergence(move_u, move_v, dxmp[0,0,:,:], dymp[0,0,:,:])



land = 'square_South_America'
landfile=Dataset(os.path.join(GFDL_BASE,'input/'+land+'/land.nc'),mode='r')

landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]
# for specified lats
landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it




sigma = 5.67*10**(-8)
net_lw_avg = sigma*(data_pert_avg.t_surf**4) - data_pert_avg.flux_lw # positive upward


SEB = - data_pert_avg.flux_sw + net_lw_avg + data_pert_avg.flux_lhe + data_pert_avg.flux_t

TOA = data_pert_avg.toa_sw - data_pert_avg.olr


N = TOA + SEB


SEB_areaavg = area_weighted_avg(SEB,area_array,landmaskxr,option='all_sfcs') #should be zero?
TOA_areaavg = area_weighted_avg(TOA,area_array,landmaskxr,option='all_sfcs') #should be zero?
area_weighted_avg(div_tot,area_array,landmaskxr,option='all_sfcs') #should be zero?

any_configuration_plot(outdir,121,481,-90.,90.,div_tot,area_array,'W/m2','div_TOT','rainnorm',landmaskxr,minval = -1000., maxval = 1000., month_annotate=0.)

diffo = np.asarray(div_tot) - np.asarray(N)
diffo = xr.DataArray(diffo, coords = [data_pert.lat,data_pert.lon], dims = ['lat','lon'])
any_configuration_plot(outdir,121,481,-90.,90.,diffo,area_array,'W/m2','div_TOT_minus_N','rainnorm',landmaskxr,minval = -1000., maxval = 1000., month_annotate=0.)


any_configuration_plot(outdir,121,481,-90.,90.,div_sens,area_array,'W/m2','div_sens','rainnorm',landmaskxr,minval = -1000., maxval = 1000., month_annotate=0.)
any_configuration_plot(outdir,121,481,-90.,90.,div_latent,area_array,'W/m2','div_latent','rainnorm',landmaskxr,minval = -1000., maxval = 1000., month_annotate=0.)
any_configuration_plot(outdir,121,481,-90.,90.,div_height,area_array,'W/m2','div_height','rainnorm',landmaskxr,minval = -1000., maxval = 1000., month_annotate=0.)

any_configuration_plot(outdir,121,481,-90.,90.,N,area_array,'W/m2','N','rainnorm',landmaskxr,minval = -100., maxval = 100., month_annotate=0.)

