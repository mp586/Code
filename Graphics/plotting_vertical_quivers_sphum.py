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

HPC = 'no'
control_dir = 'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm'
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

HPC = 'no'
testdir_in1= 'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361'
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


[omega1,omega1_avg,omega1_seasonal_avg,omega1_month_avg,omega1_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'omega','Pa/s')
[omega1_ctl,omega1_avg_ctl,omega1_seasonal_avg_ctl,omega1_month_avg_ctl,omega1_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'omega','Pa/s')
# [rh1_ctl,rh1_avg_ctl,rh1_seasonal_avg_ctl,rh1_month_avg_ctl,rh1_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'rh','%')
# [rh1,rh1_avg,rh1_seasonal_avg,rh1_month_avg,rh1_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'rh','%')
[sphum1_ctl,sphum1_avg_ctl,sphum1_seasonal_avg_ctl,sphum1_month_avg_ctl,sphum1_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'sphum','kg/kg')
[sphum1,sphum1_avg,sphum1_seasonal_avg,sphum1_month_avg,sphum1_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'sphum','kg/kg')


[ucomp1,ucomp1_avg,ucomp1_seasonal_avg,ucomp1_month_avg,ucomp1_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'ucomp','m/s')
[ucomp1_ctl,ucomp1_avg_ctl,ucomp1_seasonal_avg_ctl,ucomp1_month_avg_ctl,ucomp1_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'ucomp','m/s')

[temp1,temp1_avg,temp1_seasonal_avg,temp1_month_avg,temp1_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'temp','K')
[temp1_ctl,temp1_avg_ctl,temp1_seasonal_avg_ctl,temp1_month_avg_ctl,temp1_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'temp','K')


################ read in data from exp 2 ###############################

ctl_model = 'isca'
if (ctl_model == 'Isca') or (ctl_model == 'isca'): 
    control_model = 'Isca_DATA'
elif (ctl_model == 'gfdl') or (ctl_model == 'GFDL'):
    control_model = 'GFDL_DATA'

HPC = 'yes'
control_dir = 'two_continents_newbucket_fixedSSTs_from_realworld_zonallysymm'
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
testdir_in1= 'two_continents_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361'
runmin=120
runmax=480
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    exp1_name = 'ISCA_HPC/'+testdir_in1
    testdir = model_data + '/ISCA_HPC/' + testdir_in1
    testdir_in1 = '/ISCA_HPC/' + testdir_in1
else: 
    exp1_name = testdir_in1
    testdir = model_data + '/' + testdir_in1

land = 'two_continents'
landfile=Dataset(os.path.join(GFDL_BASE,'input/'+land+'/land.nc'),mode='r')

landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]
# for specified lats
landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it



[omega2,omega2_avg,omega2_seasonal_avg,omega2_month_avg,omega2_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'omega','Pa/s')
[omega2_ctl,omega2_avg_ctl,omega2_seasonal_avg_ctl,omega2_month_avg_ctl,omega2_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'omega','Pa/s')
# [rh2_ctl,rh2_avg_ctl,rh2_seasonal_avg_ctl,rh2_month_avg_ctl,rh2_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'rh','%')
# [rh2,rh2_avg,rh2_seasonal_avg,rh2_month_avg,rh2_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'rh','%')
[sphum2_ctl,sphum2_avg_ctl,sphum2_seasonal_avg_ctl,sphum2_month_avg_ctl,sphum2_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'sphum','kg/kg')
[sphum2,sphum2_avg,sphum2_seasonal_avg,sphum2_month_avg,sphum2_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'sphum','kg/kg')

[ucomp2,ucomp2_avg,ucomp2_seasonal_avg,ucomp2_month_avg,ucomp2_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'ucomp','m/s')
[ucomp2_ctl,ucomp2_avg_ctl,ucomp2_seasonal_avg_ctl,ucomp2_month_avg_ctl,ucomp2_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'ucomp','m/s')

[temp2,temp2_avg,temp2_seasonal_avg,temp2_month_avg,temp2_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'temp','K')
[temp2_ctl,temp2_avg_ctl,temp2_seasonal_avg_ctl,temp2_month_avg_ctl,temp2_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'temp','K')


################ read in data from exp 3 ###############################

ctl_model = 'isca'
if (ctl_model == 'Isca') or (ctl_model == 'isca'): 
    control_model = 'Isca_DATA'
elif (ctl_model == 'gfdl') or (ctl_model == 'GFDL'):
    control_model = 'GFDL_DATA'

HPC = 'yes'
control_dir = 'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_rerun'
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
testdir_in1= 'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_rerun_plus_uniform_warming_and_2xCO2_spinup_361'
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



[omega3,omega3_avg,omega3_seasonal_avg,omega3_month_avg,omega3_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'omega','Pa/s')
[omega3_ctl,omega3_avg_ctl,omega3_seasonal_avg_ctl,omega3_month_avg_ctl,omega3_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'omega','Pa/s')
# [rh3_ctl,rh3_avg_ctl,rh3_seasonal_avg_ctl,rh3_month_avg_ctl,rh3_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'rh','%')
# [rh3,rh3_avg,rh3_seasonal_avg,rh3_month_avg,rh3_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'rh','%')
[sphum3_ctl,sphum3_avg_ctl,sphum3_seasonal_avg_ctl,sphum3_month_avg_ctl,sphum3_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'sphum','kg/kg')
[sphum3,sphum3_avg,sphum3_seasonal_avg,sphum3_month_avg,sphum3_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'sphum','kg/kg')

[ucomp3,ucomp3_avg,ucomp3_seasonal_avg,ucomp3_month_avg,ucomp3_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'ucomp','m/s')
[ucomp3_ctl,ucomp3_avg_ctl,ucomp3_seasonal_avg_ctl,ucomp3_month_avg_ctl,ucomp3_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'ucomp','m/s')

[temp3,temp3_avg,temp3_seasonal_avg,temp3_month_avg,temp3_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'temp','K')
[temp3_ctl,temp3_avg_ctl,temp3_seasonal_avg_ctl,temp3_month_avg_ctl,temp3_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'temp','K')



outdir = 'Isca' + '/' + exp1_name

# vert_horiz_winds(outdir,runmin,runmax,'delta wind',ucomp2_avg_ctl - ucomp2_ctl_zonavg,(omega2_avg_ctl - omega2_ctl_zonavg)*300.,rh2_avg_ctl.sel(lat = slice(-10.,10.)).mean(dim = 'lat'),20.,90.,veclen=5,units_numerator = '', units_denom = '',save = False)


# # don't need to invert vertical axis for RH because I invert yaxis in plot, but that doesnt work for the quivers so need to invert those!
# vert_horiz_winds(outdir,runmin,runmax.'weighted change',(((ucomp2_avg - ucomp2_avg_ctl) - (ucomp1_avg - ucomp1_avg_ctl))/ucomp2_avg.max())[::-1,:,:],(((omega2_avg - omega2_avg_ctl) - (omega1_avg - omega1_avg_ctl))/omega2_avg.max())[::-1,:,:],((rh2_avg - rh2_avg_ctl) - (rh1_avg - rh1_avg_ctl)).sel(lat = slice(-10.,10.)).mean(dim = 'lat'),-10.,10.,veclen=1,units_numerator = 'Pa m', units_denom = 's s',save = False)


# #this needs to produce vectors that are at a 45 degree angle - angle in quiver has to be set to 'uv' (which is the default) or not set at all duh. Wasn't working because I had set angle to 'xy' so that it would invert the yaxis. Instead can just feed the info in reversed omega coords! 
# vert_horiz_winds(outdir,runmin,runmax,'u u',(ucomp2_avg),ucomp2_avg,((rh2_avg).sel(lat = slice(-10.,10.))).mean(dim = 'lat'),0,80.,veclen=10,units_numerator = 'Pa m', units_denom = 's s',save = False)


g = 9.81
Rspec = 287.058
pfull = temp1.pres_lev * 100 # convert from hPa to Pa
pfull = np.expand_dims(pfull, axis = 1)
pfull = np.expand_dims(pfull, axis = 2)
pfull = np.repeat(pfull, 64, axis=1)
pfull = np.repeat(pfull, 128, axis = 2)
pres_lev = temp1.pres_lev
pfull = xr.DataArray(pfull, coords = [pres_lev, temp1.lat, temp1.lon], dims = ['pres_lev','lat','lon'])
wcomp1_avg_ctl = - (omega1_avg_ctl * temp1_avg_ctl * Rspec)/(pfull * g)
wcomp2_avg_ctl = - (omega2_avg_ctl * temp2_avg_ctl * Rspec)/(pfull * g)
wcomp1_avg = - (omega1_avg * temp1_avg * Rspec)/(pfull * g)
wcomp2_avg = - (omega2_avg * temp2_avg * Rspec)/(pfull * g)
wcomp3_avg_ctl = - (omega3_avg_ctl * temp3_avg_ctl * Rspec)/(pfull * g)
wcomp3_avg = - (omega3_avg * temp3_avg * Rspec)/(pfull * g)

#conversion following https://www.ncl.ucar.edu/Document/Functions/Contributed/omega_to_w.shtml


#vert_horiz_winds(outdir,runmin,runmax,'u w*80',(((ucomp2_avg - ucomp2_avg_ctl) - (ucomp1_avg - ucomp1_avg_ctl)))[::-1,:,:],(((wcomp2_avg - wcomp2_avg_ctl) - (wcomp1_avg - wcomp1_avg_ctl))*80.)[::-1,:,:],((rh2_avg - rh2_avg_ctl) - (rh1_avg - rh1_avg_ctl)).sel(lat = slice(-10.,10.)).mean(dim = 'lat'),-10.,10.,veclen=5,units_numerator = 'm', units_denom = 's',save = False)


# Panel plot with 4 cases: America only, Africa only, Two continents, and Two continents minus America - climate change for all 



small = 18 #largefonts 14 # smallfonts 10 # medfonts = 14
med = 20 #largefonts 18 # smallfonts 14 # medfonts = 16
lge = 22 #largefonts 22 # smallfonts 18 # medfonts = 20


v = np.linspace(-.003,.003,21) # , endpoint=True)
minlat = -10.
maxlat = 10.
veclen = 10.
units_numerator = 'm'
units_denom = 's'
shiftby = 180. # = 180. --> 0 degrees in the middle, = 105. --> idealized continents overlap realistic ones 
fig, axes = plt.subplots(2, 2, sharey = True, figsize = (25,25))


# panel 1: Only South America 
uwind = (ucomp1_avg - ucomp1_avg_ctl)
wwind = ((wcomp1_avg - wcomp1_avg_ctl)*8000.)[::-1,:,:]
array = (sphum1_avg - sphum1_avg_ctl)
lons = uwind.lon 
lats = uwind.lat
pres = wwind.pres_lev
presar = array.pres_lev # need separate z coords for winds because read them in in reverse z order to flip axis
uwind, lons_cyclic = addcyclic(uwind, lons)
wwind, lons_cyclic = addcyclic(wwind, lons)

uwind = np.asarray(uwind)
wwind = np.asarray(wwind)
uwind,lons_shift = shiftgrid(np.max(lons_cyclic)-shiftby,uwind,lons_cyclic,start=False,
			   cyclic=np.max(lons_cyclic))
wwind,lons_shift = shiftgrid(np.max(lons_cyclic)-shiftby,wwind,lons_cyclic,start=False,
			   cyclic=np.max(lons_cyclic))  

array, lons_cyclic = addcyclic(array, lons)
array = np.asarray(array)
array, lons_shift = shiftgrid(np.max(lons_cyclic)-shiftby,array,lons_cyclic,
				 start=False,cyclic=np.max(lons_cyclic))

array = xr.DataArray(array,coords=[presar,lats,lons_shift],dims=['pres','lat','lon'])
uwind = xr.DataArray(uwind,coords=[pres,lats,lons_shift],dims=['pres','lat','lon'])
wwind = xr.DataArray(wwind,coords=[pres,lats,lons_shift],dims=['pres','lat','lon'])
# landmask,landlons = shiftgrid(np.max(landlons)-80.,landmask,landlons,start=False,cyclic=np.max(landlons))
# landmask, landlons = addcyclic(landmask, landlons)


Xar, Zar = np.meshgrid(lons_shift, presar)
X, Z = np.meshgrid(lons_shift, pres)


wwind_tropmean = wwind.sel(lat=slice(minlat,maxlat)).mean(dim='lat')
uwind_tropmean = uwind.sel(lat=slice(minlat,maxlat)).mean(dim='lat')
array_tropmean = array.sel(lat=slice(minlat,maxlat)).mean(dim='lat')

cset1 = axes[0,0].contourf(Xar, Zar, array_tropmean, v, cmap='BrBG', extend = 'both')
axes[0,0].set_ylabel('Pressure (hPa)', fontsize = med)

Q = axes[0,0].quiver(X[::2,::2], Z[::2,::2], uwind_tropmean[::2,::2], wwind_tropmean[::2,::2], scale = 50, scale_units = 'inches') # if angle isn't set to 'xy', can't invert yaxis on quivers, but angle 'xy' doesn't plot quivers of (u,u) in 45 degree angle! angle 'uv' which is the default does and that's what I want
# in order for the vectors to all be the same length on all panels and quiverkey to apply to all of them, set scale and scale_units 
axes[0,0].set_title('(a) America only (AM)', fontsize = med)
# Africa Only 

uwind = (ucomp3_avg - ucomp3_avg_ctl)
wwind = ((wcomp3_avg - wcomp3_avg_ctl)*8000.)[::-1,:,:]
array = (sphum3_avg - sphum3_avg_ctl)
lons = uwind.lon
lats = uwind.lat
pres = wwind.pres_lev
presar = array.pres_lev # need separate z coords for winds because read them in in reverse z order to flip axis
uwind, lons_cyclic = addcyclic(uwind, lons)
wwind, lons_cyclic = addcyclic(wwind, lons)

uwind = np.asarray(uwind)
wwind = np.asarray(wwind)
uwind,lons_shift = shiftgrid(np.max(lons_cyclic)-shiftby,uwind,lons_cyclic,start=False,
			   cyclic=np.max(lons_cyclic))
wwind,lons_shift = shiftgrid(np.max(lons_cyclic)-shiftby,wwind,lons_cyclic,start=False,
			   cyclic=np.max(lons_cyclic))  

array, lons_cyclic = addcyclic(array, lons)
array = np.asarray(array)
array, lons_shift = shiftgrid(np.max(lons_cyclic)-shiftby,array,lons_cyclic,
				 start=False,cyclic=np.max(lons_cyclic))

array = xr.DataArray(array,coords=[presar,lats,lons_shift],dims=['pres','lat','lon'])
uwind = xr.DataArray(uwind,coords=[pres,lats,lons_shift],dims=['pres','lat','lon'])
wwind = xr.DataArray(wwind,coords=[pres,lats,lons_shift],dims=['pres','lat','lon'])
# landmask,landlons = shiftgrid(np.max(landlons)-80.,landmask,landlons,start=False,cyclic=np.max(landlons))
# landmask, landlons = addcyclic(landmask, landlons)


Xar, Zar = np.meshgrid(lons_shift, presar)
X, Z = np.meshgrid(lons_shift, pres)


wwind_tropmean = wwind.sel(lat=slice(minlat,maxlat)).mean(dim='lat')
uwind_tropmean = uwind.sel(lat=slice(minlat,maxlat)).mean(dim='lat')
array_tropmean = array.sel(lat=slice(minlat,maxlat)).mean(dim='lat')

cset1 = axes[0,1].contourf(Xar, Zar, array_tropmean, v, cmap='BrBG', extend = 'both')


Q = axes[0,1].quiver(X[::2,::2], Z[::2,::2], uwind_tropmean[::2,::2], wwind_tropmean[::2,::2], scale = 50, scale_units = 'inches') # if angle isn't set to 'xy', can't invert yaxis on quivers, but angle 'xy' doesn't plot quivers of (u,u) in 45 degree angle! angle 'uv' which is the default does and that's what I want
#qk = axes[0,1].quiverkey(Q, 0.9, 0.9, veclen, str(veclen)+r'$\frac{'+units_numerator+'}{'+units_denom+'}$', labelpos='E', coordinates='figure')
axes[0,1].set_title('(b) Africa only (AF)', fontsize = med)

# Two continents 

uwind = (ucomp2_avg - ucomp2_avg_ctl)
wwind = ((wcomp2_avg - wcomp2_avg_ctl)*8000.)[::-1,:,:]
array = (sphum2_avg - sphum2_avg_ctl)
lons = uwind.lon
lats = uwind.lat
pres = wwind.pres_lev
presar = array.pres_lev # need separate z coords for winds because read them in in reverse z order to flip axis
uwind, lons_cyclic = addcyclic(uwind, lons)
wwind, lons_cyclic = addcyclic(wwind, lons)

uwind = np.asarray(uwind)
wwind = np.asarray(wwind)
uwind,lons_shift = shiftgrid(np.max(lons_cyclic)-shiftby,uwind,lons_cyclic,start=False,
			   cyclic=np.max(lons_cyclic))
wwind,lons_shift = shiftgrid(np.max(lons_cyclic)-shiftby,wwind,lons_cyclic,start=False,
			   cyclic=np.max(lons_cyclic))  

array, lons_cyclic = addcyclic(array, lons)
array = np.asarray(array)
array, lons_shift = shiftgrid(np.max(lons_cyclic)-shiftby,array,lons_cyclic,
				 start=False,cyclic=np.max(lons_cyclic))

array = xr.DataArray(array,coords=[presar,lats,lons_shift],dims=['pres','lat','lon'])
uwind = xr.DataArray(uwind,coords=[pres,lats,lons_shift],dims=['pres','lat','lon'])
wwind = xr.DataArray(wwind,coords=[pres,lats,lons_shift],dims=['pres','lat','lon'])
# landmask,landlons = shiftgrid(np.max(landlons)-80.,landmask,landlons,start=False,cyclic=np.max(landlons))
# landmask, landlons = addcyclic(landmask, landlons)


Xar, Zar = np.meshgrid(lons_shift, presar)
X, Z = np.meshgrid(lons_shift, pres)


wwind_tropmean = wwind.sel(lat=slice(minlat,maxlat)).mean(dim='lat')
uwind_tropmean = uwind.sel(lat=slice(minlat,maxlat)).mean(dim='lat')
array_tropmean = array.sel(lat=slice(minlat,maxlat)).mean(dim='lat')

cset1 = axes[1,0].contourf(Xar, Zar, array_tropmean, v, cmap='BrBG', extend = 'both')
axes[1,0].set_xlabel('Longitude E', fontsize = med)
axes[1,0].set_ylabel('Pressure (hPa)', fontsize = med)

Q = axes[1,0].quiver(X[::2,::2], Z[::2,::2], uwind_tropmean[::2,::2], wwind_tropmean[::2,::2], scale = 50, scale_units = 'inches') # if angle isn't set to 'xy', can't invert yaxis on quivers, but angle 'xy' doesn't plot quivers of (u,u) in 45 degree angle! angle 'uv' which is the default does and that's what I want
#qk = axes[1,0].quiverkey(Q, 0.9, 0.9, veclen, str(veclen)+r'$\frac{'+units_numerator+'}{'+units_denom+'}$', labelpos='E', coordinates='figure')
axes[1,0].set_title('(c) Two Continents (2C)', fontsize = med)

# Two continents - America 

uwind = ((ucomp2_avg - ucomp2_avg_ctl) - (ucomp1_avg - ucomp1_avg_ctl))[::-1,:,:]
wwind = (((wcomp2_avg - wcomp2_avg_ctl) - (wcomp1_avg - wcomp1_avg_ctl))*8000.)[::-1,:,:]
array = (sphum2_avg - sphum2_avg_ctl) - (sphum1_avg - sphum1_avg_ctl)
lons = uwind.lon
lats = uwind.lat
pres = wwind.pres_lev
presar = array.pres_lev # need separate z coords for winds because read them in in reverse z order to flip axis
uwind, lons_cyclic = addcyclic(uwind, lons)
wwind, lons_cyclic = addcyclic(wwind, lons)

uwind = np.asarray(uwind)
wwind = np.asarray(wwind)
uwind,lons_shift = shiftgrid(np.max(lons_cyclic)-shiftby,uwind,lons_cyclic,start=False,
			   cyclic=np.max(lons_cyclic))
wwind,lons_shift = shiftgrid(np.max(lons_cyclic)-shiftby,wwind,lons_cyclic,start=False,
			   cyclic=np.max(lons_cyclic))  

array, lons_cyclic = addcyclic(array, lons)
array = np.asarray(array)
array, lons_shift = shiftgrid(np.max(lons_cyclic)-shiftby,array,lons_cyclic,
				 start=False,cyclic=np.max(lons_cyclic))

array = xr.DataArray(array,coords=[presar,lats,lons_shift],dims=['pres','lat','lon'])
uwind = xr.DataArray(uwind,coords=[pres,lats,lons_shift],dims=['pres','lat','lon'])
wwind = xr.DataArray(wwind,coords=[pres,lats,lons_shift],dims=['pres','lat','lon'])
# landmask,landlons = shiftgrid(np.max(landlons)-80.,landmask,landlons,start=False,cyclic=np.max(landlons))
# landmask, landlons = addcyclic(landmask, landlons)


Xar, Zar = np.meshgrid(lons_shift, presar)
X, Z = np.meshgrid(lons_shift, pres)


wwind_tropmean = wwind.sel(lat=slice(minlat,maxlat)).mean(dim='lat')
uwind_tropmean = uwind.sel(lat=slice(minlat,maxlat)).mean(dim='lat')
array_tropmean = array.sel(lat=slice(minlat,maxlat)).mean(dim='lat')

cset1 = axes[1,1].contourf(Xar, Zar, array_tropmean, v, cmap='BrBG', extend = 'both')
axes[1,1].set_xlabel('Longitude E', fontsize = med)

Q = axes[1,1].quiver(X[::2,::2], Z[::2,::2], uwind_tropmean[::2,::2], wwind_tropmean[::2,::2], scale = 50, scale_units = 'inches') # if angle isn't set to 'xy', can't invert yaxis on quivers, but angle 'xy' doesn't plot quivers of (u,u) in 45 degree angle! angle 'uv' which is the default does and that's what I want
qk = axes[1,1].quiverkey(Q, 0.83, 0.87, veclen, str(veclen)+r'$\frac{'+units_numerator+'}{'+units_denom+'}$', labelpos='E', coordinates='figure', fontproperties={'size': med})
axes[1,1].set_title('(d) 2C - AM', fontsize = med)

fig.gca().invert_yaxis()
axes[0,0].tick_params(labelsize = small)
axes[1,0].tick_params(labelsize = small)
axes[0,1].tick_params(labelsize = small)
axes[1,1].tick_params(labelsize = small)

cbar = fig.colorbar(cset1,ax=axes)
cbar.ax.tick_params(labelsize=small)
cbar.set_label('$\Delta q$ (kg/kg)', size = med)





##### realistic continents (RC) #######

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
control_dir = 'full_continents_newbucket_fixedSSTs_zonally_symmetric'
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
testdir_in1= 'full_continents_newbucket_fixedSSTs_zonally_symmetric_plus_2pt52K_and_2xCO2_spinup_361_withtopo'
runmin=120
runmax=480
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    exp1_name = 'ISCA_HPC_'+testdir_in1
    testdir = model_data + '/ISCA_HPC/' + testdir_in1
    testdir_in1 = '/ISCA_HPC/' + testdir_in1
else: 
    exp1_name = testdir_in1
    testdir = model_data + '/' + testdir_in1

land = 'all_continents'
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


[omega1,omega1_avg,omega1_seasonal_avg,omega1_month_avg,omega1_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'omega','Pa/s')
[omega1_ctl,omega1_avg_ctl,omega1_seasonal_avg_ctl,omega1_month_avg_ctl,omega1_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'omega','Pa/s')
[rh1_ctl,rh1_avg_ctl,rh1_seasonal_avg_ctl,rh1_month_avg_ctl,rh1_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'rh','%')
[rh1,rh1_avg,rh1_seasonal_avg,rh1_month_avg,rh1_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'rh','%')
[sphum1_ctl,sphum1_avg_ctl,sphum1_seasonal_avg_ctl,sphum1_month_avg_ctl,sphum1_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'sphum','kg/kg')
[sphum1,sphum1_avg,sphum1_seasonal_avg,sphum1_month_avg,sphum1_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'sphum','kg/kg')


[ucomp1,ucomp1_avg,ucomp1_seasonal_avg,ucomp1_month_avg,ucomp1_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'ucomp','m/s')
[ucomp1_ctl,ucomp1_avg_ctl,ucomp1_seasonal_avg_ctl,ucomp1_month_avg_ctl,ucomp1_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'ucomp','m/s')

[temp1,temp1_avg,temp1_seasonal_avg,temp1_month_avg,temp1_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'temp','K')
[temp1_ctl,temp1_avg_ctl,temp1_seasonal_avg_ctl,temp1_month_avg_ctl,temp1_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'temp','K')


g = 9.81
Rspec = 287.058
pfull = temp1.pres_lev * 100 # convert from hPa to Pa
pfull = np.expand_dims(pfull, axis = 1)
pfull = np.expand_dims(pfull, axis = 2)
pfull = np.repeat(pfull, 64, axis=1)
pfull = np.repeat(pfull, 128, axis = 2)
pres_lev = temp1.pres_lev
pfull = xr.DataArray(pfull, coords = [pres_lev, temp1.lat, temp1.lon], dims = ['pres_lev','lat','lon'])
wcomp1_avg_ctl = - (omega1_avg_ctl * temp1_avg_ctl * Rspec)/(pfull * g)
wcomp1_avg = - (omega1_avg * temp1_avg * Rspec)/(pfull * g)

outdir = 'Isca' + '/' + exp1_name



small = 18 #largefonts 14 # smallfonts 10 # medfonts = 14
med = 20 #largefonts 18 # smallfonts 14 # medfonts = 16
lge = 22 #largefonts 22 # smallfonts 18 # medfonts = 20


v = np.linspace(-10.,10.,21) # , endpoint=True)
minlat = -10.
maxlat = 10.
veclen = 10.
units_numerator = 'm'
units_denom = 's'

uwind = (ucomp1_avg - ucomp1_avg_ctl)[::-1,:,:]
wwind = ((wcomp1_avg - wcomp1_avg_ctl)*8000.)[::-1,:,:]
array = (rh1_avg - rh1_avg_ctl)
lons = uwind.lon
lats = uwind.lat
pres = wwind.pres_lev
presar = array.pres_lev # need separate z coords for winds because read them in in reverse z order to flip axis
uwind, lons_cyclic = addcyclic(uwind, lons)
wwind, lons_cyclic = addcyclic(wwind, lons)

uwind = np.asarray(uwind)
wwind = np.asarray(wwind)
uwind,lons_shift = shiftgrid(np.max(lons_cyclic)-180.,uwind,lons_cyclic,start=False,
			   cyclic=np.max(lons_cyclic))
wwind,lons_shift = shiftgrid(np.max(lons_cyclic)-180.,wwind,lons_cyclic,start=False,
			   cyclic=np.max(lons_cyclic))  

array, lons_cyclic = addcyclic(array, lons)
array = np.asarray(array)
array, lons_shift = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,
				 start=False,cyclic=np.max(lons_cyclic))

array = xr.DataArray(array,coords=[presar,lats,lons_shift],dims=['pres','lat','lon'])
uwind = xr.DataArray(uwind,coords=[pres,lats,lons_shift],dims=['pres','lat','lon'])
wwind = xr.DataArray(wwind,coords=[pres,lats,lons_shift],dims=['pres','lat','lon'])
# landmask,landlons = shiftgrid(np.max(landlons)-80.,landmask,landlons,start=False,cyclic=np.max(landlons))
# landmask, landlons = addcyclic(landmask, landlons)


Xar, Zar = np.meshgrid(lons_shift, presar)
X, Z = np.meshgrid(lons_shift, pres)


wwind_tropmean = wwind.sel(lat=slice(minlat,maxlat)).mean(dim='lat')
uwind_tropmean = uwind.sel(lat=slice(minlat,maxlat)).mean(dim='lat')
array_tropmean = array.sel(lat=slice(minlat,maxlat)).mean(dim='lat')

cset1 = plt.contourf(Xar, Zar, array_tropmean, v, cmap='BrBG', extend = 'both')
plt.xlabel('Longitude E', fontsize = med)
plt.ylabel('Pressure (hPa)', fontsize = med)

Q = plt.quiver(X[::2,::2], Z[::2,::2], uwind_tropmean[::2,::2], wwind_tropmean[::2,::2], units='width') # if angle isn't set to 'xy', can't invert yaxis on quivers, but angle 'xy' doesn't plot quivers of (u,u) in 45 degree angle! angle 'uv' which is the default does and that's what I want
qk = plt.quiverkey(Q, 0.83, 0.87, veclen, str(veclen)+r'$\frac{'+units_numerator+'}{'+units_denom+'}$', labelpos='E', coordinates='figure', fontproperties={'size': med})
plt.title('RC', fontsize = med)

plt.gca().invert_yaxis()
plt.tick_params(labelsize = small)

cbar = plt.colorbar(cset1,orientation = 'vertical')
cbar.ax.tick_params(labelsize=small)
cbar.set_label('$\Delta RH$ (%)', size = med)