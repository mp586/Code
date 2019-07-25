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
control_dir = 'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_vegetation_vegpref1'
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
testdir_in1= 'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_vegetation_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361'
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

ctl_model = 'isca'
if (ctl_model == 'Isca') or (ctl_model == 'isca'): 
    control_model = 'Isca_DATA'
elif (ctl_model == 'gfdl') or (ctl_model == 'GFDL'):
    control_model = 'GFDL_DATA'

HPC = 'yes'
control_dir = 'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_vegetation_vegpref1'
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
testdir_in1= 'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_vegetation_vegpref0_plus_uniform_warming_and_2xCO2_spinup_121'
runmin=24
runmax=120
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


[net_lhe2,net_lhe2_avg,net_lhe2_seasonal_avg,net_lhe2_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[precipitation2,precipitation2_avg,precipitation2_seasonal_avg,precipitation2_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'precipitation','mm/d', factor=86400)

[net_lhe2_ctl,net_lhe2_avg_ctl,net_lhe2_seasonal_avg_ctl,net_lhe2_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[precipitation2_ctl,precipitation2_avg_ctl,precipitation2_seasonal_avg_ctl,precipitation2_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'precipitation','mm/d', factor=86400)


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




outdir = 'Isca' + '/' + exp1_name



[net_lhe3,net_lhe3_avg,net_lhe3_seasonal_avg,net_lhe3_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[precipitation3,precipitation3_avg,precipitation3_seasonal_avg,precipitation3_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'precipitation','mm/d', factor=86400)

[net_lhe3_ctl,net_lhe3_avg_ctl,net_lhe3_seasonal_avg_ctl,net_lhe3_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[precipitation3_ctl,precipitation3_avg_ctl,precipitation3_seasonal_avg_ctl,precipitation3_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'precipitation','mm/d', factor=86400)



############ plotting ##############

small = 18 #largefonts 14 # smallfonts 10 # medfonts = 14
med = 20 #largefonts 18 # smallfonts 14 # medfonts = 16
lge = 22 #largefonts 22 # smallfonts 18 # medfonts = 20

v = np.linspace(-2.,2.,21)
nmb_contours = [-2.,1.,2.]

# South America Only 

array = precipitation1_avg - precipitation1_avg_ctl
ctl_array = precipitation1_avg_ctl - net_lhe1_avg_ctl

lats=array.lat
lons=array.lon

landlats = np.asarray(landmaskxr.lat)
landlons = np.asarray(landmaskxr.lon)

landmask = np.asarray(landmaskxr)

fig, axes = plt.subplots(2,2, figsize = (25,12))
fig.subplots_adjust(hspace = 0.2, wspace = 0.05)


axes[0,0].set_title('(a) AFVP05', size = med)
#fig = plt.figure()

m = Basemap(projection='kav7',lon_0=0.,resolution='c', ax = axes[0,0])
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
m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0], fontsize=small)
m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,0], fontsize=small)

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

# Africa Only 

array = precipitation2_avg - precipitation2_avg_ctl
ctl_array = precipitation2_avg_ctl - net_lhe2_avg_ctl

lats=array.lat
lons=array.lon

landlats = np.asarray(landmaskxr.lat)
landlons = np.asarray(landmaskxr.lon)

landmask = np.asarray(landmaskxr)

axes[0,1].set_title('(b) AFVP0', size = med)
#fig = plt.figure()

m = Basemap(projection='kav7',lon_0=0.,resolution='c', ax = axes[0,1])
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

m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0], fontsize=small)
m.drawmeridians(np.arange(-180.,180.,60),labels=[0,0,0,0], fontsize=small)

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


# Two continents 

array = precipitation3_avg - precipitation3_avg_ctl
ctl_array = precipitation3_avg_ctl - net_lhe3_avg_ctl

lats=array.lat
lons=array.lon

landlats = np.asarray(landmaskxr.lat)
landlons = np.asarray(landmaskxr.lon)

landmask = np.asarray(landmaskxr)


axes[1,0].set_title('(c) AFVP1', size = med)
#fig = plt.figure()

m = Basemap(projection='kav7',lon_0=0.,resolution='c', ax = axes[1,0])
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

m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0], fontsize=small)
m.drawmeridians(np.arange(-180.,180.,60),labels=[0,0,0,0], fontsize=small)

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


# Two continents - America only 

array = (precipitation3_avg - precipitation3_avg_ctl) - (precipitation1_avg - precipitation1_avg_ctl)
ctl_array = (precipitation3_avg_ctl - net_lhe3_avg_ctl) - (precipitation1_avg - net_lhe2_avg_ctl)

lats=array.lat
lons=array.lon

landlats = np.asarray(landmaskxr.lat)
landlons = np.asarray(landmaskxr.lon)

landmask1 = np.asarray(landmaskxr)


axes[1,1].set_title('(d) AFVP1 - AFVP05', size = med)
#fig = plt.figure()

m = Basemap(projection='kav7',lon_0=0.,resolution='c', ax = axes[1,1])
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

m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0], fontsize=small)
m.drawmeridians(np.arange(-180.,180.,60),labels=[0,0,0,0], fontsize=small)

lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

cs = m.contourf(xi,yi,array, v, cmap='BrBG', extend = 'both')

cont = m.contour(xi,yi,ctl_array,nmb_contours, colors = 'k', linewidth=2) # if nmb_contours is not an int, it can be interpreted as an array specifying the contour levels


# Read landmask

# Add rectangles
#    landmask,landlons = shiftgrid(np.max(landlons)-100.,landmask,landlons,start=True,cyclic=np.max(landlons)) # this works when the array shift is commented....
landmask1,landlons1 = shiftgrid(np.max(landlons)-180.,landmask1,landlons,start=False,cyclic=np.max(landlons))

landmask1, lons_cyclic1 = addcyclic(landmask1, landlons1)

m.contour(xi,yi,landmask1, 1)



# Add Colorbar
cbar = fig.colorbar(cs, orientation = 'vertical', ax = axes, shrink = 0.65) 
cbar.set_label('mm/d', size=med)
cbar.ax.tick_params(labelsize=med)

fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_vegetation_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361/P_avg_minus_ctl_VP0_24-120_VP1_VP05_120-480.png')

landmask = np.asarray(landmaskxr)


dP_SB =  (precipitation3_avg - precipitation3_avg_ctl).where(landmask == 1).sel(lat = slice(-30.,30.))
dP_VP05 = (precipitation1_avg - precipitation1_avg_ctl).where(landmask == 1).sel(lat = slice(-30.,30.))
dP_VP0 = (precipitation2_avg - precipitation2_avg_ctl).where(landmask == 1).sel(lat = slice(-30.,30.))

dP_SB = np.asarray(dP_SB).flatten()
dP_VP05 = np.asarray(dP_VP05).flatten()
dP_VP0 = np.asarray(dP_VP0).flatten()


fig, ax = plt.subplots(1, 2, sharex = True, sharey = True, figsize = (25,12) )

ax[0].plot(dP_SB, dP_VP05, 'k.')
ax[0].set_xlim(-5.,5.)
ax[0].set_ylim(-5.,5.)
ax[0].set_ylabel('$\Delta P$ (VP05) in mm/d')
ax[0].set_xlabel('$\Delta P$ (VP1) in mm/d')

mask = ~np.isnan(dP_SB)

[k,dy,r,p,stderr] = linreg(dP_SB[mask],dP_VP05[mask]) # aa = 8.4, dq = -32
x1 = np.linspace(np.min(dP_SB[mask]),np.max(dP_SB[mask]),500)
y = k*x1 + dy
ax[0].plot(x1,y,'k-')
ax[0].annotate('r = '+str("%.2f" % r)+', p = '+str("%.5f" % p), xy=(0.05,0.05), xycoords='axes fraction')


ax[1].plot(dP_SB, dP_VP0, 'k.')
ax[1].set_ylabel('$\Delta P$ (VP0) in mm/d')
ax[1].set_xlabel('$\Delta P$ (VP1) in mm/d')

mask = ~np.isnan(dP_SB)

[k,dy,r,p,stderr] = linreg(dP_SB[mask],dP_VP0[mask]) # aa = 8.4, dq = -32
x1 = np.linspace(np.min(dP_SB[mask]),np.max(dP_SB[mask]),500)
y = k*x1 + dy
ax[1].plot(x1,y,'k-')
ax[1].annotate('r = '+str("%.2f" % r)+', p = '+str("%.5f" % p), xy=(0.05,0.05), xycoords='axes fraction')


fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_vegetation_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361/dP_VP1_vs_dP_VP05_120-480_and_VP0_24-120.png')

