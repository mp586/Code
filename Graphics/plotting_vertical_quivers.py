
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

ctl_model = input('Enter model 1 name as string ')
if (ctl_model == 'Isca') or (ctl_model == 'isca'): 
    control_model = 'Isca_DATA'
elif (ctl_model == 'gfdl') or (ctl_model == 'GFDL'):
    control_model = 'GFDL_DATA'

HPC = input('Data in ISCA_HPC ? Yes or No? ')
control_dir = input('Enter control directory from exp 1 as string ')
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    control_dir= control_model + '/ISCA_HPC/' + control_dir
else: 
    control_dir= control_model + '/' + control_dir

#print control_dir
ctl_runmin=input('Enter runmin number for exp 1 ')  # Should be a January month for seasonal variables to be correct
ctl_runmax=input('Enter runmax number for comparison 1 ')
ctl_timeseries_max = input('Enter end of ctl timeseries month for exp 1 ')

model = input('Enter model exp 1 ')
if (model == 'Isca') or (model == 'isca'): 
    model_data = 'Isca_DATA'
    output_dir1 = 'Isca'
elif (model == 'gfdl') or (model == 'GFDL'):
    model_data = 'GFDL_DATA'
    output_dir1 = ''

HPC = input('Data in ISCA_HPC ? Yes or No? ')
testdir_in1= input('Enter perturbed data directory name as string for exp 1 ')
runmin=input('Enter runmin number for exp 1 ')  # Should be a January month for seasonal variables to be correct
runmax=input('Enter runmax number for exp 1 ')
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    exp1_name = 'ISCA_HPC_'+testdir_in1
    testdir = model_data + '/ISCA_HPC/' + testdir_in1
    testdir_in1 = '/ISCA_HPC/' + testdir_in1
else: 
    exp1_name = testdir_in1
    testdir = model_data + '/' + testdir_in1

land = input('Which landmask? ')
landfile=Dataset(os.path.join(GFDL_BASE,'input/'+land+'/land.nc'),mode='r')

landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]
# for specified lats
landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

level = input('Which Level? ')


area_array, dx, dy = ca.cell_area_all(t_res=42,base_dir='/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/') # added _all because then dx and dy are also returned 
area_array = xr.DataArray(area_array) # returned in units of m bzw m^2, because radius in cell_area.py is given in metres

area_array_3D = np.expand_dims(area_array, axis=0)
area_array_3D = np.repeat(area_array_3D, 40, axis = 0) # to make area_array 3D (pressure, lat, lon)


[omega1,omega1_avg,omega1_seasonal_avg,omega1_month_avg,omega1_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'omega','Pa/s')
[omega1_ctl,omega1_avg_ctl,omega1_seasonal_avg_ctl,omega1_month_avg_ctl,omega1_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'omega','Pa/s')
[rh1_ctl,rh1_avg_ctl,rh1_seasonal_avg_ctl,rh1_month_avg_ctl,rh1_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'rh','%')
[rh1,rh1_avg,rh1_seasonal_avg,rh1_month_avg,rh1_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'rh','%')
[sphum1_ctl,sphum1_avg_ctl,sphum1_seasonal_avg_ctl,sphum1_month_avg_ctl,sphum1_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'sphum','kg/kg')
[sphum1,sphum1_avg,sphum1_seasonal_avg,sphum1_month_avg,sphum1_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'sphum','kg/kg')


[ucomp1,ucomp1_avg,ucomp1_seasonal_avg,ucomp1_month_avg,ucomp1_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'ucomp','m/s')
[ucomp1_ctl,ucomp1_avg_ctl,ucomp1_seasonal_avg_ctl,ucomp1_month_avg_ctl,ucomp1_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'ucomp','m/s')



################ read in data from exp 2 ###############################

ctl_model = input('Enter model 1 name as string ')
if (ctl_model == 'Isca') or (ctl_model == 'isca'): 
    control_model = 'Isca_DATA'
elif (ctl_model == 'gfdl') or (ctl_model == 'GFDL'):
    control_model = 'GFDL_DATA'

HPC = input('Data in ISCA_HPC ? Yes or No? ')
control_dir = input('Enter control directory from exp 1 as string ')
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    control_dir= control_model + '/ISCA_HPC/' + control_dir
else: 
    control_dir= control_model + '/' + control_dir

#print control_dir
ctl_runmin=input('Enter runmin number for exp 1 ')  # Should be a January month for seasonal variables to be correct
ctl_runmax=input('Enter runmax number for comparison 1 ')
ctl_timeseries_max = input('Enter end of ctl timeseries month for exp 1 ')

model = input('Enter model exp 1 ')
if (model == 'Isca') or (model == 'isca'): 
    model_data = 'Isca_DATA'
    output_dir1 = 'Isca'
elif (model == 'gfdl') or (model == 'GFDL'):
    model_data = 'GFDL_DATA'
    output_dir1 = ''

HPC = input('Data in ISCA_HPC ? Yes or No? ')
testdir_in1= input('Enter perturbed data directory name as string for exp 1 ')
runmin=input('Enter runmin number for exp 1 ')  # Should be a January month for seasonal variables to be correct
runmax=input('Enter runmax number for exp 1 ')
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    exp1_name = 'ISCA_HPC_'+testdir_in1
    testdir = model_data + '/ISCA_HPC/' + testdir_in1
    testdir_in1 = '/ISCA_HPC/' + testdir_in1
else: 
    exp1_name = testdir_in1
    testdir = model_data + '/' + testdir_in1

land = input('Which landmask? ')
landfile=Dataset(os.path.join(GFDL_BASE,'input/'+land+'/land.nc'),mode='r')

landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]
# for specified lats
landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

level = input('Which Level? ')


[omega2,omega2_avg,omega2_seasonal_avg,omega2_month_avg,omega2_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'omega','Pa/s')
[omega2_ctl,omega2_avg_ctl,omega2_seasonal_avg_ctl,omega2_month_avg_ctl,omega2_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'omega','Pa/s')
[rh2_ctl,rh2_avg_ctl,rh2_seasonal_avg_ctl,rh2_month_avg_ctl,rh2_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'rh','%')
[rh2,rh2_avg,rh2_seasonal_avg,rh2_month_avg,rh2_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'rh','%')
[sphum2_ctl,sphum2_avg_ctl,sphum2_seasonal_avg_ctl,sphum2_month_avg_ctl,sphum2_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'sphum','kg/kg')
[sphum2,sphum2_avg,sphum2_seasonal_avg,sphum2_month_avg,sphum2_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'sphum','kg/kg')

[ucomp2,ucomp2_avg,ucomp2_seasonal_avg,ucomp2_month_avg,ucomp2_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'ucomp','m/s')
[ucomp2_ctl,ucomp2_avg_ctl,ucomp2_seasonal_avg_ctl,ucomp2_month_avg_ctl,ucomp2_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'ucomp','m/s')

ucomp2_ctl_zonavg = np.expand_dims(ucomp2_avg_ctl.mean(dim = 'lon'), axis = 2)
ucomp2_ctl_zonavg = np.repeat(ucomp2_ctl_zonavg, 128, axis = 2)
np.shape(ucomp2_ctl_zonavg)

omega2_ctl_zonavg = np.expand_dims(omega2_avg_ctl.mean(dim = 'lon'), axis = 2)
omega2_ctl_zonavg = np.repeat(omega2_ctl_zonavg, 128, axis = 2)

vert_horiz_winds(outdir,runmin,runmax,'u w*80',ucomp2_avg_ctl - ucomp2_ctl_zonavg,(omega2_avg_ctl - omega2_ctl_zonavg)*300.,rh2_avg_ctl.sel(lat = slice(-10.,10.)).mean(dim = 'lat'),20.,90.,veclen=5,units_numerator = '', units_denom = '',save = False)


# don't need to invert vertical axis for RH because I invert yaxis in plot, but that doesnt work for the quivers so need to invert those!
vert_horiz_winds(outdir,runmin,runmax,'u w*80',(((ucomp2_avg - ucomp2_avg_ctl) - (ucomp1_avg - ucomp1_avg_ctl))/ucomp2_avg.max())[::-1,:,:],(((omega2_avg - omega2_avg_ctl) - (omega1_avg - omega1_avg_ctl))/omega2_avg.max())[::-1,:,:],((rh2_avg - rh2_avg_ctl) - (rh1_avg - rh1_avg_ctl)).sel(lat = slice(-10.,10.)).mean(dim = 'lat'),-10.,10.,veclen=1,units_numerator = 'Pa m', units_denom = 's s',save = False)


#this needs to produce vectors that are at a 45 degree angle - angle in quiver has to be set to 'uv' (which is the default) or not set at all duh. Wasn't working because I had set angle to 'xy' so that it would invert the yaxis. Instead can just feed the info in reversed omega coords! 
vert_horiz_winds(outdir,runmin,runmax,'u w*80',(ucomp2_avg),ucomp2_avg,((rh2_avg).sel(lat = slice(-10.,10.))).mean(dim = 'lat'),0,80.,veclen=10,units_numerator = 'Pa m', units_denom = 's s',save = False)
