# hovmoeller plot
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

ctl_model = input('Enter model name as string ')
if (ctl_model == 'Isca') or (ctl_model == 'isca'): 
    control_model = 'Isca_DATA'
elif (ctl_model == 'gfdl') or (ctl_model == 'GFDL'):
    control_model = 'GFDL_DATA'


control_dir= control_model + '/' + input('Enter control directory name as string ')
#print control_dir
ctl_runmin=input('Enter runmin number ')  # Should be a January month for seasonal variables to be correct
ctl_runmax=input('Enter runmax number for comparison ')
ctl_timeseries_max = input('Enter end of ctl timeseries month ')

model = input('Enter model ')
if (model == 'Isca') or (model == 'isca'): 
    model_data = 'Isca_DATA'
    output_dir = 'Isca'
elif (model == 'gfdl') or (model == 'GFDL'):
    model_data = 'GFDL_DATA'
    output_dir = ''

testdir= input('Enter data directory name as string ')
runmin=input('Enter runmin number ')  # Should be a January month for seasonal variables to be correct
runmax=input('Enter runmax number ')


outdir = output_dir + '/' + testdir
testdir = model_data + '/' + testdir

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

# Read in variables MONTHLY
[precipitation,precipitation_avg,precipitation_seasonal_avg,precipitation_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'precipitation','mm/d', factor=86400)
#[CIWV,CIWV_avg,CIWV_seasonal_avg,CIWV_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'sphum','kg/kg',level='all')

[precipitation_ctl,precipitation_avg_ctl,precipitation_seasonal_avg_ctl,precipitation_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'precipitation','mm/d', factor=86400)
#[CIWV_ctl,CIWV_avg_ctl,CIWV_seasonal_avg_ctl,CIWV_month_avg

precipitation_month_avg_ctl, lons_cyclic = addcyclic(precipitation_month_avg_ctl, precipitation_month_avg_ctl.lon)
precipitation_month_avg_ctl,lons_shift = shiftgrid(180., precipitation_month_avg_ctl, lons_cyclic, start = False, cyclic = np.max(lons_cyclic)) 


fig, ax = plt.subplots(1,1) 
ax.contourf(precipitation_month_avg.month, lons_shift, precipitation_month_avg_ctl[:,40,:].T, cmap = 'Blues') 
fig.savefig('/scratch/mp586/Code/Graphics/hovmoeller_shifted.png', bbox_inches='tight') 


# #################################################

# Read in variables DAILY
#[precipitation,precipitation_avg,time]=seasonal_surface_variable_daily(testdir,model,runmin,runmax,'precipitation','mm/d', factor=86400)
#[CIWV,CIWV_avg,time]=seasonal_surface_variable_daily(testdir,model,runmin,runmax,'sphum','kg/kg',level='all')
days = np.linspace(1,2880,2880)
[precipitation_ctl,precipitation_avg_ctl,time]=seasonal_surface_variable_daily(control_dir,ctl_model,ctl_runmin,ctl_runmax,'precipitation','mm/d', factor=86400)
[CIWV_ctl,CIWV_avg_ctl,time]=seasonal_surface_variable_daily(control_dir,ctl_model,ctl_runmin,ctl_runmax,'sphum','kg/kg',level='all')
#precip_daily_avg = precipitation_ctl.groupby('time.day').mean('time')
#precip_ctl, lons_cyclic = addcyclic(precip_daily_avg, precipitation_ctl.lon)
precip_ctl, lons_cyclic = addcyclic(precipitation_ctl, precipitation_ctl.lon)
precip_ctl,lons_shift = shiftgrid(180., precip_ctl, lons_cyclic, start = False, cyclic = np.max(lons_cyclic)) 
fig, ax = plt.subplots(1,1) 
ax.contourf(lons_shift, days, precip_ctl[:,40,:], cmap = 'Blues')
fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/squareland_newbucket_fixedSSTs_from_realworld_zonallysymm_vegetation_vegpref0_plus_uniform_warming_and_2xCO2_spinup_121/hovmoeller_'+str(runmin)+'-'+str(runmax)+'precip_daily.png', format = 'png', bbox_inches='tight')

ciwv_ctl, lons_cyclic = addcyclic(CIWV_ctl, precipitation_ctl.lon)
ciwv_ctl,lons_shift = shiftgrid(180., ciwv_ctl, lons_cyclic, start = False, cyclic = np.max(lons_cyclic)) 
fig, ax = plt.subplots(1,1) 
ax.contourf(lons_shift, days[::5], ciwv_ctl[::5,40,:], cmap = 'Blues')
fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/squareland_newbucket_fixedSSTs_from_realworld_zonallysymm_vegetation_vegpref0_plus_uniform_warming_and_2xCO2_spinup_121/hovmoeller_'+str(runmin)+'-'+str(runmax)+'_CIWV_daily_lastyear.png', format = 'png', bbox_inches='tight')


### tried some sort of daily averaging -- average over the same 'day' in each month, although also realised then that this calendar assumes 31 days/month but we have 30 days per month
# precip_daily_avg = precipitation_ctl.groupby('time.day').mean('time')
# precip_ctl, lons_cyclic = addcyclic(precip_daily_avg, precipitation_ctl.lon)
# precip_ctl,lons_shift = shiftgrid(180., precip_ctl, lons_cyclic, start = False, cyclic = np.max(lons_cyclic)) 
# fig, ax = plt.subplots(1,1) 
# ax.contourf(precip_daily_avg.day, lons_shift, precip_ctl[:,40,:].T, cmap = 'Blues')
# fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/squareland_newbucket_fixedSSTs_from_realworld_zonallysymm_vegetation_vegpref0_plus_uniform_warming_and_2xCO2_spinup_121/hovmoeller_'+str(runmin)+'-'+str(runmax)+'precip_dailyavg.png', format = 'png', bbox_inches='tight')

# ciwv_daily_avg = CIWV_ctl.groupby('time.day').mean('time')
# ciwv_ctl, lons_cyclic = addcyclic(ciwv_daily_avg, precipitation_ctl.lon)
# ciwv_ctl,lons_shift = shiftgrid(180., ciwv_ctl, lons_cyclic, start = False, cyclic = np.max(lons_cyclic)) 
# fig, ax = plt.subplots(1,1) 
# ax.contourf(ciwv_daily_avg.day,lons_shift, ciwv_ctl[:,40,:].T, cmap = 'Blues')
# fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/squareland_newbucket_fixedSSTs_from_realworld_zonallysymm_vegetation_vegpref0_plus_uniform_warming_and_2xCO2_spinup_121/hovmoeller_'+str(runmin)+'-'+str(runmax)+'_CIWV_dailyavg.png', format = 'png', bbox_inches='tight')
