from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import xarray as xr
import pandas as pd
import os
from matplotlib.patches import Rectangle

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
dire = testdir_in1
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

landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]
# for specified lats
landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it


area_array, dx, dy = ca.cell_area_all(t_res=42,base_dir='/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/') # added _all because then dx and dy are also returned 
area_array = xr.DataArray(area_array) # returned in units of m bzw m^2, because radius in cell_area.py is given in metres

area_array_3D = np.expand_dims(area_array, axis=0)
area_array_3D = np.repeat(area_array_3D, 40, axis = 0) # to make area_array 3D (pressure, lat, lon)

data = xr.open_mfdataset('/scratch/mp586/'+control_dir+'/*/atmos_monthly_interp.nc')
data = data.mean('time')

omega1_avg_ctl = data.omega
rh1_avg_ctl = data.rh
sphum1_avg_ctl = data.sphum
ucomp1_avg_ctl = data.ucomp
temp1_avg_ctl = data.temp

data = xr.open_mfdataset('/scratch/mp586/'+testdir+'/*/atmos_monthly_interp.nc')
data = data.mean('time')

omega1_avg = data.omega
rh1_avg = data.rh
sphum1_avg = data.sphum
ucomp1_avg = data.ucomp
temp1_avg = data.temp


# [omega1,omega1_avg,omega1_seasonal_avg,omega1_month_avg,omega1_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'omega','Pa/s')
# [omega1_ctl,omega1_avg_ctl,omega1_seasonal_avg_ctl,omega1_month_avg_ctl,omega1_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'omega','Pa/s')
# [rh1_ctl,rh1_avg_ctl,rh1_seasonal_avg_ctl,rh1_month_avg_ctl,rh1_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'rh','%')
# [rh1,rh1_avg,rh1_seasonal_avg,rh1_month_avg,rh1_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'rh','%')
# [sphum1_ctl,sphum1_avg_ctl,sphum1_seasonal_avg_ctl,sphum1_month_avg_ctl,sphum1_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'sphum','kg/kg')
# [sphum1,sphum1_avg,sphum1_seasonal_avg,sphum1_month_avg,sphum1_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'sphum','kg/kg')


# [ucomp1,ucomp1_avg,ucomp1_seasonal_avg,ucomp1_month_avg,ucomp1_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'ucomp','m/s')
# [ucomp1_ctl,ucomp1_avg_ctl,ucomp1_seasonal_avg_ctl,ucomp1_month_avg_ctl,ucomp1_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'ucomp','m/s')

# [temp1,temp1_avg,temp1_seasonal_avg,temp1_month_avg,temp1_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'temp','K')
# [temp1_ctl,temp1_avg_ctl,temp1_seasonal_avg_ctl,temp1_month_avg_ctl,temp1_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'temp','K')


################ read in data from exp 2 ###############################

ctl_model = 'isca'
if (ctl_model == 'Isca') or (ctl_model == 'isca'): 
    control_model = 'Isca_DATA'
elif (ctl_model == 'gfdl') or (ctl_model == 'GFDL'):
    control_model = 'GFDL_DATA'

HPC = 'yes'
control_dir = 'two_continents_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387'
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
testdir_in1= 'two_continents_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387'
runmin=120
runmax=480
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    exp1_name = 'ISCA_HPC/'+testdir_in1
    testdir = model_data + '/ISCA_HPC/' + testdir_in1
    testdir_in1 = '/ISCA_HPC/' + testdir_in1
else: 
    exp1_name = testdir_in1
    testdir = model_data + '/' + testdir_in1

land = 'two_AM'
landfile=Dataset(os.path.join(GFDL_BASE,'input/'+land+'/land.nc'),mode='r')

landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]
# for specified lats
landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it


data = xr.open_mfdataset('/scratch/mp586/'+control_dir+'/*/atmos_monthly_interp.nc')
data = data.mean('time')

omega2_avg_ctl = data.omega
rh2_avg_ctl = data.rh
sphum2_avg_ctl = data.sphum
ucomp2_avg_ctl = data.ucomp
temp2_avg_ctl = data.temp

data = xr.open_mfdataset('/scratch/mp586/'+testdir+'/*/atmos_monthly_interp.nc')
data = data.mean('time')

omega2_avg = data.omega
rh2_avg = data.rh
sphum2_avg = data.sphum
ucomp2_avg = data.ucomp
temp2_avg = data.temp


# [omega2,omega2_avg,omega2_seasonal_avg,omega2_month_avg,omega2_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'omega','Pa/s')
# [omega2_ctl,omega2_avg_ctl,omega2_seasonal_avg_ctl,omega2_month_avg_ctl,omega2_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'omega','Pa/s')
# [rh2_ctl,rh2_avg_ctl,rh2_seasonal_avg_ctl,rh2_month_avg_ctl,rh2_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'rh','%')
# [rh2,rh2_avg,rh2_seasonal_avg,rh2_month_avg,rh2_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'rh','%')
# [sphum2_ctl,sphum2_avg_ctl,sphum2_seasonal_avg_ctl,sphum2_month_avg_ctl,sphum2_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'sphum','kg/kg')
# [sphum2,sphum2_avg,sphum2_seasonal_avg,sphum2_month_avg,sphum2_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'sphum','kg/kg')

# [ucomp2,ucomp2_avg,ucomp2_seasonal_avg,ucomp2_month_avg,ucomp2_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'ucomp','m/s')
# [ucomp2_ctl,ucomp2_avg_ctl,ucomp2_seasonal_avg_ctl,ucomp2_month_avg_ctl,ucomp2_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'ucomp','m/s')

# [temp2,temp2_avg,temp2_seasonal_avg,temp2_month_avg,temp2_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'temp','K')
# [temp2_ctl,temp2_avg_ctl,temp2_seasonal_avg_ctl,temp2_month_avg_ctl,temp2_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'temp','K')


################ read in data from exp 3 ###############################

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
ctl_runmin=121  # Should be a January month for seasonal variables to be correct
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
    exp1_name = 'ISCA_HPC/'+testdir_in1
    testdir = model_data + '/ISCA_HPC/' + testdir_in1
    testdir_in1 = '/ISCA_HPC/' + testdir_in1
else: 
    exp1_name = testdir_in1
    testdir = model_data + '/' + testdir_in1

land = 'square_Africa'
landfile=Dataset(os.path.join(GFDL_BASE,'input/'+land+'/land.nc'),mode='r')

landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]
# for specified lats
landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

data = xr.open_mfdataset('/scratch/mp586/'+control_dir+'/*/atmos_monthly_interp.nc')
data = data.mean('time')


omega3_avg_ctl = data.omega
rh3_avg_ctl = data.rh
sphum3_avg_ctl = data.sphum
ucomp3_avg_ctl = data.ucomp
temp3_avg_ctl = data.temp

data = xr.open_mfdataset('/scratch/mp586/'+testdir+'/*/atmos_monthly_interp.nc')
data = data.mean('time')

omega3_avg = data.omega
rh3_avg = data.rh
sphum3_avg = data.sphum
ucomp3_avg = data.ucomp
temp3_avg = data.temp

# [omega3,omega3_avg,omega3_seasonal_avg,omega3_month_avg,omega3_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'omega','Pa/s')
# [omega3_ctl,omega3_avg_ctl,omega3_seasonal_avg_ctl,omega3_month_avg_ctl,omega3_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'omega','Pa/s')
# [rh3_ctl,rh3_avg_ctl,rh3_seasonal_avg_ctl,rh3_month_avg_ctl,rh3_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'rh','%')
# [rh3,rh3_avg,rh3_seasonal_avg,rh3_month_avg,rh3_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'rh','%')
# [sphum3_ctl,sphum3_avg_ctl,sphum3_seasonal_avg_ctl,sphum3_month_avg_ctl,sphum3_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'sphum','kg/kg')
# [sphum3,sphum3_avg,sphum3_seasonal_avg,sphum3_month_avg,sphum3_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'sphum','kg/kg')

# [ucomp3,ucomp3_avg,ucomp3_seasonal_avg,ucomp3_month_avg,ucomp3_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'ucomp','m/s')
# [ucomp3_ctl,ucomp3_avg_ctl,ucomp3_seasonal_avg_ctl,ucomp3_month_avg_ctl,ucomp3_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'ucomp','m/s')

# [temp3,temp3_avg,temp3_seasonal_avg,temp3_month_avg,temp3_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'temp','K')
# [temp3_ctl,temp3_avg_ctl,temp3_seasonal_avg_ctl,temp3_month_avg_ctl,temp3_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'temp','K')



outdir = 'Isca' + '/' + exp1_name

# vert_horiz_winds(outdir,runmin,runmax,'delta wind',ucomp2_avg_ctl - ucomp2_ctl_zonavg,(omega2_avg_ctl - omega2_ctl_zonavg)*300.,rh2_avg_ctl.sel(lat = slice(-10.,10.)).mean(dim = 'lat'),20.,90.,veclen=5,units_numerator = '', units_denom = '',save = False)


# # don't need to invert vertical axis for RH because I invert yaxis in plot, but that doesnt work for the quivers so need to invert those!
# vert_horiz_winds(outdir,runmin,runmax.'weighted change',(((ucomp2_avg - ucomp2_avg_ctl) - (ucomp1_avg - ucomp1_avg_ctl))/ucomp2_avg.max())[::-1,:,:],(((omega2_avg - omega2_avg_ctl) - (omega1_avg - omega1_avg_ctl))/omega2_avg.max())[::-1,:,:],((rh2_avg - rh2_avg_ctl) - (rh1_avg - rh1_avg_ctl)).sel(lat = slice(-10.,10.)).mean(dim = 'lat'),-10.,10.,veclen=1,units_numerator = 'Pa m', units_denom = 's s',save = False)


# #this needs to produce vectors that are at a 45 degree angle - angle in quiver has to be set to 'uv' (which is the default) or not set at all duh. Wasn't working because I had set angle to 'xy' so that it would invert the yaxis. Instead can just feed the info in reversed omega coords! 
# vert_horiz_winds(outdir,runmin,runmax,'u u',(ucomp2_avg),ucomp2_avg,((rh2_avg).sel(lat = slice(-10.,10.))).mean(dim = 'lat'),0,80.,veclen=10,units_numerator = 'Pa m', units_denom = 's s',save = False)


g = 9.81
Rspec = 287.058
pfull = data.pfull * 100 # convert from hPa to Pa
pfull = np.expand_dims(pfull, axis = 1)
pfull = np.expand_dims(pfull, axis = 2)
pfull = np.repeat(pfull, 64, axis=1)
pfull = np.repeat(pfull, 128, axis = 2)
pres_lev = data.pfull
pfull = xr.DataArray(pfull, coords = [pres_lev, data.lat, data.lon], dims = ['pfull','lat','lon'])
wcomp1_avg_ctl = - (omega1_avg_ctl * temp1_avg_ctl * Rspec)/(pfull * g)
wcomp2_avg_ctl = - (omega2_avg_ctl * temp2_avg_ctl * Rspec)/(pfull * g)
wcomp1_avg = - (omega1_avg * temp1_avg * Rspec)/(pfull * g)
wcomp2_avg = - (omega2_avg * temp2_avg * Rspec)/(pfull * g)
wcomp3_avg_ctl = - (omega3_avg_ctl * temp3_avg_ctl * Rspec)/(pfull * g)
wcomp3_avg = - (omega3_avg * temp3_avg * Rspec)/(pfull * g)

#conversion following https://www.ncl.ucar.edu/Document/Functions/Contributed/omega_to_w.shtml


#vert_horiz_winds(outdir,runmin,runmax,'u w*80',(((ucomp2_avg - ucomp2_avg_ctl) - (ucomp1_avg - ucomp1_avg_ctl)))[::-1,:,:],(((wcomp2_avg - wcomp2_avg_ctl) - (wcomp1_avg - wcomp1_avg_ctl))*80.)[::-1,:,:],((rh2_avg - rh2_avg_ctl) - (rh1_avg - rh1_avg_ctl)).sel(lat = slice(-10.,10.)).mean(dim = 'lat'),-10.,10.,veclen=5,units_numerator = 'm', units_denom = 's',save = False)


# Panel plot with 4 cases: America only, Africa only, Two continents, and Two continents minus America - climate change for all 


# quiver_plots_4cases(runmin, runmax, 'ucorr_interplevs_quivers_4cases_deltarh_fct_latweights_', 'BrBG', dire, landmaskxr, '$\Delta$ r (%)', ucomp1_avg, ucomp1_avg_ctl, wcomp1_avg, wcomp1_avg_ctl, ucomp2_avg, ucomp2_avg_ctl, wcomp2_avg, wcomp2_avg_ctl, ucomp3_avg, ucomp3_avg_ctl, wcomp3_avg, wcomp3_avg_ctl, (rh1_avg - rh1_avg_ctl), (rh2_avg - rh2_avg_ctl), (rh3_avg - rh3_avg_ctl), minval=-10., maxval=10., vertmult = 8000, minlat=-10., maxlat=10.)

quiver_plots_4cases(runmin, runmax, 'ucorr_interplevs_quivers_4cases_deltasphum_gkg_fct_latweights_', 'BrBG', dire, landmaskxr, '$\Delta$ q (g/kg)', ucomp1_avg, ucomp1_avg_ctl, wcomp1_avg, wcomp1_avg_ctl, ucomp2_avg, ucomp2_avg_ctl, wcomp2_avg, wcomp2_avg_ctl, ucomp3_avg, ucomp3_avg_ctl, wcomp3_avg, wcomp3_avg_ctl, (sphum1_avg - sphum1_avg_ctl)*1000, (sphum2_avg - sphum2_avg_ctl)*1000, (sphum3_avg - sphum3_avg_ctl)*1000, minval=-3., maxval=3., vertmult = 8000, minlat=-10., maxlat=10.)

# quiver_plots_4cases(runmin, runmax, 'ucorr_interplevs_quivers_4cases_deltatemp_fct_latweights_', 'RdBu_r', dire, landmaskxr, '$\Delta$ T (K)', ucomp1_avg, ucomp1_avg_ctl, wcomp1_avg, wcomp1_avg_ctl, ucomp2_avg, ucomp2_avg_ctl, wcomp2_avg, wcomp2_avg_ctl, ucomp3_avg, ucomp3_avg_ctl, wcomp3_avg, wcomp3_avg_ctl, (temp1_avg - temp1_avg_ctl), (temp2_avg - temp2_avg_ctl), (temp3_avg - temp3_avg_ctl), minval=-10., maxval=10., vertmult = 8000, minlat=-10., maxlat=10.)
