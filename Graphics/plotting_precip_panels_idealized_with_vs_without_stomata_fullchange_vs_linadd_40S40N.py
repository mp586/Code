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
    exp1_name = 'ISCA_HPC_'+testdir_in1
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


area_array, dx, dy = ca.cell_area_all(t_res=42,base_dir='/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/') # added _all because then dx and dy are also returned 
area_array = xr.DataArray(area_array) # returned in units of m bzw m^2, because radius in cell_area.py is given in metres

area_array_3D = np.expand_dims(area_array, axis=0)
area_array_3D = np.repeat(area_array_3D, 40, axis = 0) # to make area_array 3D (pressure, lat, lon)


[net_lhe1,net_lhe1_avg,net_lhe1_seasonal_avg,net_lhe1_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[precipitation1,precipitation1_avg,precipitation1_seasonal_avg,precipitation1_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'precipitation','mm/d', factor=86400)

[net_lhe1_ctl,net_lhe1_avg_ctl,net_lhe1_seasonal_avg_ctl,net_lhe1_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[precipitation1_ctl,precipitation1_avg_ctl,precipitation1_seasonal_avg_ctl,precipitation1_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'precipitation','mm/d', factor=86400)



################ read in data from exp 2 ###############################

model = 'isca'
if (model == 'Isca') or (model == 'isca'): 
    model_data = 'Isca_DATA'
    output_dir1 = 'Isca'
elif (model == 'gfdl') or (model == 'GFDL'):
    model_data = 'GFDL_DATA'
    output_dir1 = ''

HPC = 'yes'
testdir_in1= 'two_continents_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387'
dire = testdir_in1
runmin=120
runmax=360
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    exp1_name = 'ISCA_HPC/'+testdir_in1
    testdir = model_data + '/ISCA_HPC/' + testdir_in1
    testdir_in1 = '/ISCA_HPC/' + testdir_in1
else: 
    exp1_name = testdir_in1
    testdir = model_data + '/' + testdir_in1

[net_lhe2,net_lhe2_avg,net_lhe2_seasonal_avg,net_lhe2_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[precipitation2,precipitation2_avg,precipitation2_seasonal_avg,precipitation2_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'precipitation','mm/d', factor=86400)


model = 'isca'
if (model == 'Isca') or (model == 'isca'): 
    model_data = 'Isca_DATA'
    output_dir1 = 'Isca'
elif (model == 'gfdl') or (model == 'GFDL'):
    model_data = 'GFDL_DATA'
    output_dir1 = ''

HPC = 'yes'
testdir_in1= 'two_continents_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_commit7bb4387'
dire = testdir_in1
runmin=120
runmax=360
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    exp1_name = 'ISCA_HPC/'+testdir_in1
    testdir = model_data + '/ISCA_HPC/' + testdir_in1
    testdir_in1 = '/ISCA_HPC/' + testdir_in1
else: 
    exp1_name = testdir_in1
    testdir = model_data + '/' + testdir_in1

[net_lhe3,net_lhe3_avg,net_lhe3_seasonal_avg,net_lhe3_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[precipitation3,precipitation3_avg,precipitation3_seasonal_avg,precipitation3_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'precipitation','mm/d', factor=86400)

############ plotting ##############

small = 14 #largefonts 14 # smallfonts 10 # medfonts = 14
med = 20 #largefonts 18 # smallfonts 14 # medfonts = 16
lge = 24 #largefonts 22 # smallfonts 18 # medfonts = 20

v = np.linspace(-2.,2.,21)
nmb_contours = [-2.,1.,2.]

# RC 

array = precipitation2_avg  - precipitation1_avg_ctl
#array = precipitation2_avg - net_lhe2_avg - (precipitation2_avg_ctl - net_lhe2_avg_ctl)
ctl_array = precipitation1_avg_ctl - net_lhe1_avg_ctl

lats=array.lat
lons=array.lon

landlats = np.asarray(landmaskxr.lat)
landlons = np.asarray(landmaskxr.lon)

landmask = np.asarray(landmaskxr)


fig, axes = plt.subplots(3,1, figsize = (25,15))

axes[0].set_title('(a) Full Change', size = med)
#fig = plt.figure()

m = Basemap(projection='cyl',resolution='c', ax = axes[0],llcrnrlat=-40, urcrnrlat=40,llcrnrlon=-180, urcrnrlon=180)
array = xr.DataArray(array,coords=[lats,lons],dims=['lat','lon'])

array = np.asarray(array)
array, lons_cyclic = addcyclic(array, lons)
array,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))

array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])

ctl_array = np.asarray(ctl_array)
ctl_array, lons_cyclic = addcyclic(ctl_array, lons)
ctl_array,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,ctl_array,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))
ctl_array = xr.DataArray(ctl_array,coords=[lats,lons_cyclic],dims=['lat','lon'])

lons = lons_cyclic
m.drawparallels(np.arange(-40.,40.,20.),labels=[1,0,0,0], fontsize=small)
m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1], fontsize=small)

lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

cs = m.contourf(xi,yi,array, v, cmap='BrBG', extend = 'both')

cont = m.contour(xi,yi,ctl_array,nmb_contours, colors = 'k', linewidth=2) # if nmb_contours is not an int, it can be interpreted as an array specifying the contour levels


# Read landmask

# Add rectangles
#    landmask,landlons = shiftgrid(np.max(landlons)-100.,landmask,landlons,start=True,cyclic=np.max(landlons)) # this works when the array shift is commented....
landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))

landmask, lons_cyclic = addcyclic(landmask, landlons)

if np.any(landmask != 0.):
    m.contour(xi,yi,landmask, 1)

# RC07

# array = precipitation1_avg - net_lhe1_avg - ( precipitation1_avg_ctl - net_lhe1_avg_ctl)
array = precipitation1_avg - precipitation1_avg_ctl + precipitation3_avg - precipitation1_avg_ctl

lats=array.lat
lons=array.lon

landlats = np.asarray(landmaskxr.lat)
landlons = np.asarray(landmaskxr.lon)

landmask = np.asarray(landmaskxr)

axes[1].set_title('(b) Linear addition', size = med)
#fig = plt.figure()
m = Basemap(projection='cyl',resolution='c', ax = axes[1],llcrnrlat=-40, urcrnrlat=40,llcrnrlon=-180, urcrnrlon=180)
array = xr.DataArray(array,coords=[lats,lons],dims=['lat','lon'])

array = np.asarray(array)
array, lons_cyclic = addcyclic(array, lons)
array,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))

array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])

lons = lons_cyclic

m.drawparallels(np.arange(-40.,40.,20.),labels=[1,0,0,0], fontsize=small)
m.drawmeridians(np.arange(-180.,180.,60),labels=[0,0,0,1], fontsize=small)

lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

cs = m.contourf(xi,yi,array, v, cmap='BrBG', extend = 'both')



# Read landmask

# Add rectangles
#    landmask,landlons = shiftgrid(np.max(landlons)-100.,landmask,landlons,start=True,cyclic=np.max(landlons)) # this works when the array shift is commented....
landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))

landmask, lons_cyclic = addcyclic(landmask, landlons)

if np.any(landmask != 0.):
    m.contour(xi,yi,landmask, 1)

# array =  (precipitation2_avg - net_lhe2_avg - (precipitation2_avg_ctl - net_lhe2_avg_ctl)) - (precipitation1_avg - net_lhe1_avg - ( precipitation1_avg_ctl - net_lhe1_avg_ctl))
array =  (precipitation2_avg - precipitation1_avg_ctl) - (precipitation1_avg - precipitation1_avg_ctl + precipitation3_avg - precipitation1_avg_ctl)

lats=array.lat
lons=array.lon

landlats = np.asarray(landmaskxr.lat)
landlons = np.asarray(landmaskxr.lon)

landmask = np.asarray(landmaskxr)

axes[2].set_title('(c) Full minus linear addition', size = med)
#fig = plt.figure()

m = Basemap(projection='cyl',resolution='c', ax = axes[2],llcrnrlat=-40, urcrnrlat=40,llcrnrlon=-180, urcrnrlon=180)
array = xr.DataArray(array,coords=[lats,lons],dims=['lat','lon'])

array = np.asarray(array)
array, lons_cyclic = addcyclic(array, lons)
array,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))

array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])

lons = lons_cyclic
m.drawparallels(np.arange(-40.,40.,20.),labels=[1,0,0,0], fontsize=small)
m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1], fontsize=small)

lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

cs = m.contourf(xi,yi,array, v, cmap='BrBG', extend = 'both')

# Read landmask

# Add rectangles
#    landmask,landlons = shiftgrid(np.max(landlons)-100.,landmask,landlons,start=True,cyclic=np.max(landlons)) # this works when the array shift is commented....
landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))

landmask, lons_cyclic = addcyclic(landmask, landlons)

if np.any(landmask != 0.):
    m.contour(xi,yi,landmask, 1)



# Add Colorbar
cbar = fig.colorbar(cs, orientation = 'vertical', ax = axes, shrink = 0.5) # usually on right 
cbar.set_label('mm/d', size=med)
cbar.ax.tick_params(labelsize=med)

fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/'+dire+'/Pavg_minus_ctl_fullchange_vs_linear_addition_40S40N_120-480.png', bbox_inches='tight', dpi=100)
fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/'+dire+'/Pavg_minus_ctl_fullchange_vs_linear_addition_40S40N_120-480.svg', bbox_inches='tight', dpi=100)

# fig.savefig('/scratch/mp586/Code/Graphics/Isca/full_continents_newbucket_fixedSSTs_zonally_symmetric_vegetation_vegpref02_plus_2pt52K_and_2xCO2_spinup_361_witholr/Pavg_minus_ctl_bucket_vs_VP02_P-Econts_40S40N_ctl_25-121_pert_24-120.png', bbox_inches='tight', dpi=100)
# fig.savefig('/scratch/mp586/Code/Graphics/Isca/full_continents_newbucket_fixedSSTs_zonally_symmetric_vegetation_vegpref02_plus_2pt52K_and_2xCO2_spinup_361_witholr/Pavg_minus_ctl_bucket_vs_VP02_P-Econts_40S40N_ctl_25-121_pert_24-120.svg', bbox_inches='tight', dpi=100)
