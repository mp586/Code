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

GFDL_BASE = os.environ['GFDL_BASE']
sys.path.insert(0, os.path.join(GFDL_BASE,'src/extra/python/scripts'))
import cell_area as ca

ctl_model = 'isca'
if (ctl_model == 'Isca') or (ctl_model == 'isca'): 
    control_model = 'Isca_DATA'
elif (ctl_model == 'gfdl') or (ctl_model == 'GFDL'):
    control_model = 'GFDL_DATA'

HPC = 'yes'
control_dir = 'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387'
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    control_dir= control_model + '/ISCA_HPC/' + control_dir
else: 
    control_dir= control_model + '/' + control_dir

#print control_dir
ctl_runmin=121
ctl_runmax=481

model = 'isca'
if (model == 'Isca') or (model == 'isca'): 
    model_data = 'Isca_DATA'
    output_dir1 = 'Isca'
elif (model == 'gfdl') or (model == 'GFDL'):
    model_data = 'GFDL_DATA'
    output_dir1 = ''

HPC = 'yes'
testdir_in1= 'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387'
runmin=120
runmax=480
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    exp1_name = 'ISCA_HPC_'+testdir_in1
    testdir = model_data + '/ISCA_HPC/' + testdir_in1
    testdir_in1 = '/ISCA_HPC/' + testdir_in1
else: 
    exp1_name = testdir_in1
    testdir = model_data + '/' + testdir_in1

land = 'square_South_America'
landfile=Dataset(os.path.join(GFDL_BASE,'input/'+land+'/land.nc'),mode='r')

landmaskAM=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]
# for specified lats
landmaskAMxr=xr.DataArray(landmaskAM,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it


area_array, dx, dy = ca.cell_area_all(t_res=42,base_dir='/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/') # added _all because then dx and dy are also returned 
area_array = xr.DataArray(area_array) # returned in units of m bzw m^2, because radius in cell_area.py is given in metres

area_array_3D = np.expand_dims(area_array, axis=0)
area_array_3D = np.repeat(area_array_3D, 40, axis = 0) # to make area_array 3D (pressure, lat, lon)


[net_lheAM,net_lheAM_avg,net_lheAM_seasonal_avg,net_lheAM_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[precipitationAM,precipitationAM_avg,precipitationAM_seasonal_avg,precipitationAM_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'precipitation','mm/d', factor=86400)

[net_lheAM_ctl,net_lheAM_avg_ctl,net_lheAM_seasonal_avg_ctl,net_lheAM_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[precipitationAM_ctl,precipitationAM_avg_ctl,precipitationAM_seasonal_avg_ctl,precipitationAM_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'precipitation','mm/d', factor=86400)


testdir_in1= 'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387'
dire = testdir_in1
runmin=120
runmax=480
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    exp1_name = 'ISCA_HPC/'+testdir_in1
    testdir = model_data + '/ISCA_HPC/' + testdir_in1
    testdir_in1 = '/ISCA_HPC/' + testdir_in1
else: 
    exp1_name = testdir_in1
    testdir = model_data + '/' + testdir_in1

[net_lheAMveg,net_lheAMveg_avg,net_lheAMveg_seasonal_avg,net_lheAMveg_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[precipitationAMveg,precipitationAMveg_avg,precipitationAMveg_seasonal_avg,precipitationAMveg_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'precipitation','mm/d', factor=86400)




ctl_model = 'isca'
if (ctl_model == 'Isca') or (ctl_model == 'isca'): 
    control_model = 'Isca_DATA'
elif (ctl_model == 'gfdl') or (ctl_model == 'GFDL'):
    control_model = 'GFDL_DATA'

HPC = 'yes'
control_dir = 'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387'
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    control_dir= control_model + '/ISCA_HPC/' + control_dir
else: 
    control_dir= control_model + '/' + control_dir

#print control_dir
ctl_runmin=121
ctl_runmax=481

model = 'isca'
if (model == 'Isca') or (model == 'isca'): 
    model_data = 'Isca_DATA'
    output_dir1 = 'Isca'
elif (model == 'gfdl') or (model == 'GFDL'):
    model_data = 'GFDL_DATA'
    output_dir1 = ''

HPC = 'yes'
testdir_in1= 'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387'
runmin=120
runmax=480
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    exp1_name = 'ISCA_HPC_'+testdir_in1
    testdir = model_data + '/ISCA_HPC/' + testdir_in1
    testdir_in1 = '/ISCA_HPC/' + testdir_in1
else: 
    exp1_name = testdir_in1
    testdir = model_data + '/' + testdir_in1

land = 'square_Africa'
landfile=Dataset(os.path.join(GFDL_BASE,'input/'+land+'/land.nc'),mode='r')

landmaskAF=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]
# for specified lats
landmaskxrAF=xr.DataArray(landmaskAF,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it


area_array, dx, dy = ca.cell_area_all(t_res=42,base_dir='/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/') # added _all because then dx and dy are also returned 
area_array = xr.DataArray(area_array) # returned in units of m bzw m^2, because radius in cell_area.py is given in metres

area_array_3D = np.expand_dims(area_array, axis=0)
area_array_3D = np.repeat(area_array_3D, 40, axis = 0) # to make area_array 3D (pressure, lat, lon)


[net_lheAF,net_lheAF_avg,net_lheAF_seasonal_avg,net_lheAF_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[precipitationAF,precipitationAF_avg,precipitationAF_seasonal_avg,precipitationAF_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'precipitation','mm/d', factor=86400)

[net_lheAF_ctl,net_lheAF_avg_ctl,net_lheAF_seasonal_avg_ctl,net_lheAF_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[precipitationAF_ctl,precipitationAF_avg_ctl,precipitationAF_seasonal_avg_ctl,precipitationAF_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'precipitation','mm/d', factor=86400)


testdir_in1= 'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387'
runmin=120
runmax=480
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    exp1_name = 'ISCA_HPC/'+testdir_in1
    testdir = model_data + '/ISCA_HPC/' + testdir_in1
    testdir_in1 = '/ISCA_HPC/' + testdir_in1
else: 
    exp1_name = testdir_in1
    testdir = model_data + '/' + testdir_in1

[net_lheAFveg,net_lheAFveg_avg,net_lheAFveg_seasonal_avg,net_lheAFveg_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[precipitationAFveg,precipitationAFveg_avg,precipitationAFveg_seasonal_avg,precipitationAFveg_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'precipitation','mm/d', factor=86400)


fig = plt.figure()

plt.plot(precipitationAM_avg_ctl.where(landmaskAM == 1.).sel(lat = slice(-10., 10.)), net_lheAM_avg_ctl.where(landmaskAM == 1.).sel(lat = slice(-10., 10.)), 'c.', markersize = 4., label = 'AM SB ctl')
plt.plot(precipitationAM_avg.where(landmaskAM == 1.).sel(lat = slice(-10., 10.)), net_lheAM_avg.where(landmaskAM == 1.).sel(lat = slice(-10., 10.)), 'c+', markersize = 4., label = 'AM SB pert')
plt.plot(precipitationAMveg_avg.where(landmaskAM == 1.).sel(lat = slice(-10., 10.)), net_lheAMveg_avg.where(landmaskAM == 1.).sel(lat = slice(-10., 10.)), 'c1', markersize = 4., label = 'AM VP pert')

plt.plot(precipitationAF_avg_ctl.where(landmaskAF == 1.).sel(lat = slice(-10., 10.)), net_lheAF_avg_ctl.where(landmaskAF == 1.).sel(lat = slice(-10., 10.)), 'm.', markersize = 4., label = 'AF SB ctl')
plt.plot(precipitationAF_avg.where(landmaskAF == 1.).sel(lat = slice(-10., 10.)), net_lheAF_avg.where(landmaskAF == 1.).sel(lat = slice(-10., 10.)), 'm+', markersize = 4., label = 'AF SB pert')
plt.plot(precipitationAFveg_avg.where(landmaskAF == 1.).sel(lat = slice(-10., 10.)), net_lheAFveg_avg.where(landmaskAF == 1.).sel(lat = slice(-10., 10.)), 'm1', markersize = 4. , label = 'AF VP pert')

plt.xlabel('P (mm/day)')
plt.ylabel('E (mm/day)')

fig.legend(bbox_to_anchor=(0.85,0.4))

fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/'+dire+'/budyko_curve.pdf', bbox_inches = 'tight', dpi=400)
fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/'+dire+'/budyko_curve.png', bbox_inches = 'tight', dpi=400)
