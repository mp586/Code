from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from netCDF4 import Dataset
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
ctl_model = 'isca'
if (ctl_model == 'Isca') or (ctl_model == 'isca'): 
    control_model = 'Isca_DATA'
elif (ctl_model == 'gfdl') or (ctl_model == 'GFDL'):
    control_model = 'GFDL_DATA'

HPC = 'yes'
#control_dir = 'full_continents_newbucket_fixedSSTs_zonally_symmetric_commit7bb4387'
control_dir = 'full_continents_newbucket_fixedSSTs_from_realworld_zonallysymm_commitfe93b9d'

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

#testdir_in1= 'full_continents_newbucket_fixedSSTs_zonally_symmetric_corrected_vegpref05_plus_2pt52K_and_2xCO2_spinup_361_commit7bb4387'
testdir_in1= 'full_continents_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commitfe93b9d'

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

[net_lhe2,net_lhe2_avg,net_lhe2_seasonal_avg,net_lhe2_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[precipitation2,precipitation2_avg,precipitation2_seasonal_avg,precipitation2_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'precipitation','mm/d', factor=86400)

[net_lhe2_ctl,net_lhe2_avg_ctl,net_lhe2_seasonal_avg_ctl,net_lhe2_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[precipitation2_ctl,precipitation2_avg_ctl,precipitation2_seasonal_avg_ctl,precipitation2_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'precipitation','mm/d', factor=86400)

small = 22
med = 24
lge = 28

v = np.linspace(-2.,2.,41)
nmb_contours = [0.,4.,8.]

# RC 



land = 'all_continents'
landfile=Dataset(os.path.join(GFDL_BASE,'input/'+land+'/land.nc'),mode='r')

landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]
# for specified lats
landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

area_array, dx, dy = ca.cell_area_all(t_res=42,base_dir='/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/') # added _all because then dx and dy are also returned 
area_array = xr.DataArray(area_array) # returned in units of m bzw m^2, because radius in cell_area.py is given in metres
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
testdir_in1= 'two_continents_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387'
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


[net_lhe4,net_lhe4_avg,net_lhe4_seasonal_avg,net_lhe4_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[precipitation4,precipitation4_avg,precipitation4_seasonal_avg,precipitation4_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'precipitation','mm/d', factor=86400)

[net_lhe4_ctl,net_lhe4_avg_ctl,net_lhe4_seasonal_avg_ctl,net_lhe4_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[precipitation4_ctl,precipitation4_avg_ctl,precipitation4_seasonal_avg_ctl,precipitation4_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'precipitation','mm/d', factor=86400)




array = precipitation2_avg  - precipitation2_avg_ctl

lats=array.lat
lons=array.lon

arrayid = precipitation4_avg - precipitation4_avg_ctl

land = 'two_continents'
landfile=Dataset(os.path.join(GFDL_BASE,'input/'+land+'/land.nc'),mode='r')

landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]
# for specified lats
landmaskxrid=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

fig, axes = plt.subplots(1,2, figsize = (9,5))
m = Basemap(projection='cyl',resolution='c', ax = axes[0],llcrnrlat=-40, urcrnrlat=40,llcrnrlon=-100, urcrnrlon=-20)

array = xr.DataArray(array,coords=[lats,lons],dims=['lat','lon'])

array = np.asarray(array)
array, lons_cyclic = addcyclic(array, lons)
array,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))

array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])
lons = lons_cyclic
m.drawparallels(np.arange(-30.,40.,10.),labels=[1,0,0,0], fontsize=16)
m.drawmeridians(np.arange(-80.,-20.,20.),labels=[0,0,0,1], fontsize=16)

lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

cs = m.contourf(xi,yi,array, v, cmap='BrBG', extend = 'both')

landlats = np.asarray(landmaskxr.lat)
landlons = np.asarray(landmaskxr.lon)

landmask = np.asarray(landmaskxr)
landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))

landmask, lons_cyclic = addcyclic(landmask, landlons)

m.contour(xi,yi,landmask, 1, linewidths = 2, colors='k')

m = Basemap(projection='cyl',resolution='c', ax = axes[1],llcrnrlat=-40, urcrnrlat=40,llcrnrlon=-20, urcrnrlon=60)

lats=arrayid.lat
lons=arrayid.lon

array = xr.DataArray(arrayid,coords=[lats,lons],dims=['lat','lon'])

array = np.asarray(array)
array, lons_cyclic = addcyclic(array, lons)
array,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))

array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])
lons = lons_cyclic
m.drawparallels(np.arange(-30.,40.,10.),labels=[])
m.drawmeridians(np.arange(0.,60.,20.),labels=[0,0,0,1], fontsize=16)

lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

cs = m.contourf(xi,yi,array, v, cmap='BrBG', extend = 'both')

landlats = np.asarray(landmaskxrid.lat)
landlons = np.asarray(landmaskxrid.lon)

landmask = np.asarray(landmaskxrid)
landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))

landmask, lons_cyclic = addcyclic(landmask, landlons)

m.contour(xi,yi,landmask, 1, linewidths = 2, colors='k')

cbar = fig.colorbar(cs, orientation = 'vertical', ax = axes, shrink=0.72) # usually on right 
cbar.set_label('mm/d', size=18)
cbar.ax.tick_params(labelsize=18)


axes[0].set_title('(a) RealCont $\Delta P_{50\%cond}$',fontsize=18)
axes[1].set_title('(b) 2Cont $\Delta P_{50\%cond}$',fontsize=18)

# plt.subplots_adjust(wspace=0.02)
# plt.tight_layout()

plt.savefig('/scratch/mp586/Code/Graphics/AM_SA_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commitfe93b9d.svg', bbox_inches='tight')

plt.close()




