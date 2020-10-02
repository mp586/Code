from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import xarray as xr
import pandas as pd
import os
import scipy

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
control_dir = 'full_continents_newbucket_fixedSSTs_zonally_symmetric_commit7bb4387'
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
testdir_in1= 'full_continents_newbucket_fixedSSTs_zonally_symmetric_plus_2pt52K_and_2xCO2_spinup_361_commit7bb4387'
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


[net_lhe1,net_lhe1_avg,net_lhe1_seasonal_avg,net_lhe1_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[precipitation1,precipitation1_avg,precipitation1_seasonal_avg,precipitation1_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'precipitation','mm/d', factor=86400)

[net_lhe1_ctl,net_lhe1_avg_ctl,net_lhe1_seasonal_avg_ctl,net_lhe1_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[precipitation1_ctl,precipitation1_avg_ctl,precipitation1_seasonal_avg_ctl,precipitation1_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'precipitation','mm/d', factor=86400)

pctl_annual_mean = precipitation1_ctl.groupby('time.year').mean(dim='time')
pavg_annual_mean = precipitation1.groupby('time.year').mean(dim='time')
[ttest_pctl_pavg,pval_pctl_pavg] = scipy.stats.ttest_ind(pctl_annual_mean,pavg_annual_mean,equal_var=False)
# pval says something about the significance, calculated with Welch's t-test 


small = 14 #largefonts 14 # smallfonts 10 # medfonts = 14
med = 22 #largefonts 18 # smallfonts 14 # medfonts = 16
lge = 26 #largefonts 22 # smallfonts 18 # medfonts = 20

v = np.linspace(-1.,1.,21)

# South America Only 

array = precipitation1_avg - precipitation1_avg_ctl
ctl_array = precipitation1_avg_ctl - net_lhe1_avg_ctl


lats=ctl_array.lat
lons=ctl_array.lon

landlats = np.asarray(landmaskxr.lat)
landlons = np.asarray(landmaskxr.lon)

landmask = np.asarray(landmaskxr)


fig, axes = plt.subplots(1, 1, figsize = (25,10))

#fig = plt.figure()

m = Basemap(projection='kav7',lon_0=0.,resolution='c')
array = xr.DataArray(array,coords=[lats,lons],dims=['lat','lon'])
pval_pctl_pavg = xr.DataArray(pval_pctl_pavg,coords=[lats,lons],dims=['lat','lon'])


array = np.asarray(array)
array, lons_cyclic = addcyclic(array, lons)
array,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))
array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])



pval_pctl_pavg = np.asarray(pval_pctl_pavg)
pval_pctl_pavg, lons_cyclic = addcyclic(pval_pctl_pavg, lons)
pval_pctl_pavg,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,pval_pctl_pavg,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))
pval_pctl_pavg = xr.DataArray(pval_pctl_pavg,coords=[lats,lons_cyclic],dims=['lat','lon'])

lons = lons_cyclic
m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0], fontsize=small)
m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1], fontsize=small)

lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

cs = m.contourf(xi,yi,array, v, cmap='BrBG', extend = 'both')
stip = m.contourf(xi,yi,pval_pctl_pavg.where(pval_pctl_pavg < 0.05), hatches = '.', alpha = 0)
# cont = m.contour(xi,yi,ctl_array,nmb_contours, colors = 'k', linewidth=2) # if nmb_contours is not an int, it can be interpreted as an array specifying the contour levels


# Read landmask

# Add rectangles
#    landmask,landlons = shiftgrid(np.max(landlons)-100.,landmask,landlons,start=True,cyclic=np.max(landlons)) # this works when the array shift is commented....
landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))

landmask, lons_cyclic = addcyclic(landmask, landlons)

if np.any(landmask != 0.):
    m.contour(xi,yi,landmask, 1)
cbar = fig.colorbar(cs, orientation = 'vertical', shrink = 0.8) # usually on right 
cbar.set_label('mm/d', size=med)
cbar.ax.tick_params(labelsize=med)

land = 'all_continents'
landfile=Dataset(os.path.join(GFDL_BASE,'input/'+land+'/land.nc'),mode='r')

landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]
# for specified lats
landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/'+dire+'/P_avg_minus_ctl_120-480_stippling.png', bbox_inches='tight', dpi=100)

fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/'+dire+'/P_avg_minus_ctl_120-480_stippling.svg', bbox_inches='tight', dpi=100)
fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/'+dire+'/P_avg_minus_ctl_120-480_stippling.pdf', bbox_inches='tight', dpi=400)





##### statistical tests 

# check that the datapoints are independent i.e. no auto-correlation:
rs = np.empty_like(pctl_annual_mean[0])

for i in range(len(precipitation1.lat)):
    for j in range(len(precipitation1.lon)):
        r, p = scipy.stats.pearsonr(pctl_annual_mean[1:30,i,j],pctl_annual_mean[0:29,i,j])
        rs[i,j] = r

plt.plot(np.sort(rs.flatten()))

rs_pearson = rs

rs = np.empty_like(pctl_annual_mean[0])

for i in range(len(precipitation1.lat)):
    for j in range(len(precipitation1.lon)):
        r, p = scipy.stats.pearsonr(pavg_annual_mean[1:30,i,j],pavg_annual_mean[0:29,i,j]) 
        rs[i,j] = r

plt.plot(np.sort(rs.flatten()))

rs_pearson_pert = rs


small = 14 #largefonts 14 # smallfonts 10 # medfonts = 14
med = 22 #largefonts 18 # smallfonts 14 # medfonts = 16
lge = 26 #largefonts 22 # smallfonts 18 # medfonts = 20

v = np.linspace(-0.3,0.3,7)
nmb_contours = [.05,.1]

# South America Only 

array = rs_pearson
ctl_array = precipitation1_avg_ctl - net_lhe1_avg_ctl

lats=ctl_array.lat
lons=ctl_array.lon

landlats = np.asarray(landmaskxr.lat)
landlons = np.asarray(landmaskxr.lon)

landmask = np.asarray(landmaskxr)


fig, axes = plt.subplots(1, 1, figsize = (28,10))

axes.set_title('Autocorrelation', size = lge)
#fig = plt.figure()

m = Basemap(projection='kav7',lon_0=0.,resolution='c')
array = xr.DataArray(array,coords=[lats,lons],dims=['lat','lon'])

array = np.asarray(array)
array, lons_cyclic = addcyclic(array, lons)
array,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))

array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])

# array2 = xr.DataArray(rs_pearson_pert,coords=[lats,lons],dims=['lat','lon'])

# array2 = np.asarray(array2)
# array2, lons_cyclic = addcyclic(array2, lons)
# array2,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))

# array2 = xr.DataArray(array2,coords=[lats,lons_cyclic],dims=['lat','lon'])

lons = lons_cyclic
m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0], fontsize=med)
m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1], fontsize=med)

lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

cs = m.contourf(xi,yi,array.where(np.abs(array) < 0.3), v, cmap='RdBu')

# cont = m.contour(xi,yi,ctl_array,nmb_contours, colors = 'k', linewidth=2) # if nmb_contours is not an int, it can be interpreted as an array specifying the contour levels


# Read landmask

# Add rectangles
#    landmask,landlons = shiftgrid(np.max(landlons)-100.,landmask,landlons,start=True,cyclic=np.max(landlons)) # this works when the array shift is commented....
landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))

landmask, lons_cyclic = addcyclic(landmask, landlons)

if np.any(landmask != 0.):
    m.contour(xi,yi,landmask, 1)
cbar = fig.colorbar(cs, orientation = 'vertical', shrink = 0.8) # usually on right 
cbar.set_label('Correlation coefficient', size=med)
cbar.ax.tick_params(labelsize=med)

land = 'all_continents'
landfile=Dataset(os.path.join(GFDL_BASE,'input/'+land+'/land.nc'),mode='r')

landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]
# for specified lats
landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

fig.savefig('/scratch/mp586/Code/Graphics/Pearson_r_strict_ctl.pdf', bbox_inches='tight')



# test whether the assumption of normal distribution applies 

rs = np.empty_like(pctl_annual_mean[0])

for i in range(len(precipitation1.lat)):
    for j in range(len(precipitation1.lon)):
        a, rs[i,j] = scipy.stats.shapiro(pctl_annual_mean[:,i,j])

rs_shapiro = rs

rs = np.empty_like(pctl_annual_mean[0])

for i in range(len(precipitation1.lat)):
    for j in range(len(precipitation1.lon)):
        a, rs[i,j] = scipy.stats.shapiro(pavg_annual_mean[:,i,j])

rs_shapiro_pert = rs

array = rs_shapiro
ctl_array = precipitation1_avg_ctl - net_lhe1_avg_ctl

lats=ctl_array.lat
lons=ctl_array.lon

landlats = np.asarray(landmaskxr.lat)
landlons = np.asarray(landmaskxr.lon)

landmask = np.asarray(landmaskxr)


fig, axes = plt.subplots(1, 1, figsize = (28,10))

axes.set_title('Shapiro Wilk test', size = lge)
#fig = plt.figure()

m = Basemap(projection='kav7',lon_0=0.,resolution='c')
array = xr.DataArray(array,coords=[lats,lons],dims=['lat','lon'])

array = np.asarray(array)
array, lons_cyclic = addcyclic(array, lons)
array,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))

array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])


array2 = rs_shapiro_pert
array2 = xr.DataArray(array2,coords=[lats,lons],dims=['lat','lon'])

array2 = np.asarray(array2)
array2, lons_cyclic = addcyclic(array2, lons)
array2,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,array2,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))

array2 = xr.DataArray(array2,coords=[lats,lons_cyclic],dims=['lat','lon'])


lons = lons_cyclic
m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0], fontsize=med)
m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1], fontsize=med)

lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

cs = m.contourf(xi,yi,array.where(array>0.05).where(array2>0.05), np.linspace(0.,1.,21), cmap='Blues')

# cont = m.contour(xi,yi,ctl_array,nmb_contours, colors = 'k', linewidth=2) # if nmb_contours is not an int, it can be interpreted as an array specifying the contour levels


# Read landmask

# Add rectangles
#    landmask,landlons = shiftgrid(np.max(landlons)-100.,landmask,landlons,start=True,cyclic=np.max(landlons)) # this works when the array shift is commented....
landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))

landmask, lons_cyclic = addcyclic(landmask, landlons)

if np.any(landmask != 0.):
    m.contour(xi,yi,landmask, 1)
cbar = fig.colorbar(cs, orientation = 'vertical', shrink = 0.8) # usually on right 
cbar.ax.tick_params(labelsize=med)
cbar.set_label('p-value', fontsize = lge)
land = 'all_continents'
landfile=Dataset(os.path.join(GFDL_BASE,'input/'+land+'/land.nc'),mode='r')

landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]
# for specified lats
landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it
fig.savefig('/scratch/mp586/Code/Graphics/Shapiro_p_ctlpert.pdf', bbox_inches='tight')



small = 14 #largefonts 14 # smallfonts 10 # medfonts = 14
med = 22 #largefonts 18 # smallfonts 14 # medfonts = 16
lge = 26 #largefonts 22 # smallfonts 18 # medfonts = 20

v = np.linspace(-1.,1.,21)

# South America Only 

array = precipitation1_avg - precipitation1_avg_ctl
ctl_array = precipitation1_avg_ctl - net_lhe1_avg_ctl


lats=ctl_array.lat
lons=ctl_array.lon

landlats = np.asarray(landmaskxr.lat)
landlons = np.asarray(landmaskxr.lon)

landmask = np.asarray(landmaskxr)


fig, axes = plt.subplots(1, 1, figsize = (28,10))

#fig = plt.figure()

m = Basemap(projection='kav7',lon_0=0.,resolution='c')
array = xr.DataArray(array,coords=[lats,lons],dims=['lat','lon'])
# pval_pctl_pavg = xr.DataArray(pval_pctl_pavg,coords=[lats,lons],dims=['lat','lon'])
rs_shapiro = xr.DataArray(rs_shapiro,coords=[lats,lons],dims=['lat','lon'])
rs_pearson = xr.DataArray(rs_pearson,coords=[lats,lons],dims=['lat','lon'])


array = np.asarray(array)
array, lons_cyclic = addcyclic(array, lons)
array,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))
array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])



pval_pctl_pavg = np.asarray(pval_pctl_pavg)
pval_pctl_pavg, lons_cyclic = addcyclic(pval_pctl_pavg, lons)
pval_pctl_pavg,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,pval_pctl_pavg,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))
pval_pctl_pavg = xr.DataArray(pval_pctl_pavg,coords=[lats,lons_cyclic],dims=['lat','lon'])

rs_pearson = np.asarray(rs_pearson)
rs_pearson, lons_cyclic = addcyclic(rs_pearson, lons)
rs_pearson,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,rs_pearson,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))
rs_pearson = xr.DataArray(rs_pearson,coords=[lats,lons_cyclic],dims=['lat','lon'])

rs_shapiro = np.asarray(rs_shapiro)
rs_shapiro, lons_cyclic = addcyclic(rs_shapiro, lons)
rs_shapiro,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,rs_shapiro,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))
rs_shapiro = xr.DataArray(rs_shapiro,coords=[lats,lons_cyclic],dims=['lat','lon'])

# rs_pearson_pert = np.asarray(rs_pearson_pert)
# rs_pearson_pert, lons_cyclic = addcyclic(rs_pearson_pert, lons)
# rs_pearson_pert,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,rs_pearson_pert,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))
# rs_pearson_pert = xr.DataArray(rs_pearson_pert,coords=[lats,lons_cyclic],dims=['lat','lon'])

# rs_shapiro_pert = np.asarray(rs_shapiro_pert)
# rs_shapiro_pert, lons_cyclic = addcyclic(rs_shapiro_pert, lons)
# rs_shapiro_pert,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,rs_shapiro_pert,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))
# rs_shapiro_pert = xr.DataArray(rs_shapiro_pert,coords=[lats,lons_cyclic],dims=['lat','lon'])


lons = lons_cyclic
m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0], fontsize=med)
m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1], fontsize=med)

lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

cs = m.contourf(xi,yi,array.where(np.abs(rs_pearson) < 0.3).where(rs_shapiro > 0.05), v, cmap='BrBG', extend = 'both')
stip = m.contourf(xi,yi,pval_pctl_pavg.where(pval_pctl_pavg < 0.05), hatches = '.', alpha = 0)
# cont = m.contour(xi,yi,ctl_array,nmb_contours, colors = 'k', linewidth=2) # if nmb_contours is not an int, it can be interpreted as an array specifying the contour levels


# Read landmask

# Add rectangles
#    landmask,landlons = shiftgrid(np.max(landlons)-100.,landmask,landlons,start=True,cyclic=np.max(landlons)) # this works when the array shift is commented....
landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))

landmask, lons_cyclic = addcyclic(landmask, landlons)

if np.any(landmask != 0.):
    m.contour(xi,yi,landmask, 1)
cbar = fig.colorbar(cs, orientation = 'vertical', shrink = 0.8) # usually on right 
cbar.set_label('mm/d', size=med)
cbar.ax.tick_params(labelsize=med)

land = 'all_continents'
landfile=Dataset(os.path.join(GFDL_BASE,'input/'+land+'/land.nc'),mode='r')

landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]
# for specified lats
landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/'+dire+'/P_avg_minus_ctl_120-480_stippling_shapiro_strictpearson.png', bbox_inches='tight', dpi=100)

fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/'+dire+'/P_avg_minus_ctl_120-480_stippling_shapiro_strictpearson.svg', bbox_inches='tight', dpi=100)
fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/'+dire+'/P_avg_minus_ctl_120-480_stippling_shapiro_strictpearson.pdf', bbox_inches='tight', dpi=400)



