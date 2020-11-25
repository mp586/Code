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

landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]
# for specified lats
landmaskxrSA=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it


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
testdir_in1= 'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387'
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

land = 'square_Africa'
landfile=Dataset(os.path.join(GFDL_BASE,'input/'+land+'/land.nc'),mode='r')

landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]
# for specified lats
landmaskxrAF=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it


[net_lhe3,net_lhe3_avg,net_lhe3_seasonal_avg,net_lhe3_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[precipitation3,precipitation3_avg,precipitation3_seasonal_avg,precipitation3_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'precipitation','mm/d', factor=86400)

[net_lhe3_ctl,net_lhe3_avg_ctl,net_lhe3_seasonal_avg_ctl,net_lhe3_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[precipitation3_ctl,precipitation3_avg_ctl,precipitation3_seasonal_avg_ctl,precipitation3_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'precipitation','mm/d', factor=86400)


################ read in data from exp 3 ###############################

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



############ plotting ##############

small = 14 #largefonts 14 # smallfonts 10 # medfonts = 14
med = 20 #largefonts 18 # smallfonts 14 # medfonts = 16
lge = 24 #largefonts 22 # smallfonts 18 # medfonts = 20

v = np.linspace(-2.,2.,21)
nmb_contours = [-2.,1.,2.]

# RC 

array = precipitation3_avg - precipitation3_avg_ctl - (precipitation1_avg  - precipitation1_avg_ctl)
#array = precipitation2_avg - net_lhe2_avg - (precipitation2_avg_ctl - net_lhe2_avg_ctl)

lats=array.lat
lons=array.lon

fig, axes = plt.subplots(3,1, figsize = (10,8))

axes[1].set_title('(b) 2C - AM (bucket)', size = med)
#fig = plt.figure()

m = Basemap(projection='cyl',resolution='c', ax = axes[1],llcrnrlat=-40, urcrnrlat=40,llcrnrlon=-180, urcrnrlon=180)
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


landlats = np.asarray(landmaskxrAF.lat)
landlons = np.asarray(landmaskxrAF.lon)

landmask1 = np.asarray(landmaskxrSA)
landmask2 = np.asarray(landmaskxrAF)

landmask1,landlons1 = shiftgrid(np.max(landlons)-180.,landmask1,landlons,start=False,cyclic=np.max(landlons))

landmask1, lons_cyclic1 = addcyclic(landmask1, landlons1)

landmask2,landlons2 = shiftgrid(np.max(landlons)-180.,landmask2,landlons,start=False,cyclic=np.max(landlons))

landmask2, lons_cyclic2 = addcyclic(landmask2, landlons2)


m.contour(xi,yi,landmask1, 1)
m.contour(xi,yi,landmask2, 1, linestyles = 'dotted')

# RC07

# array = precipitation1_avg - net_lhe1_avg - ( precipitation1_avg_ctl - net_lhe1_avg_ctl)
array = precipitation4_avg - precipitation4_avg_ctl - (precipitation2_avg - precipitation2_avg_ctl)

lats=array.lat
lons=array.lon

axes[0].set_title('(a) 2C - AM (50%cond)', size = med)
#fig = plt.figure()
m = Basemap(projection='cyl',resolution='c', ax = axes[0],llcrnrlat=-40, urcrnrlat=40,llcrnrlon=-180, urcrnrlon=180)
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




m.contour(xi,yi,landmask1, 1)
m.contour(xi,yi,landmask2, 1, linestyles = 'dotted')


# array =  (precipitation2_avg - net_lhe2_avg - (precipitation2_avg_ctl - net_lhe2_avg_ctl)) - (precipitation1_avg - net_lhe1_avg - ( precipitation1_avg_ctl - net_lhe1_avg_ctl))
array =  ((precipitation4_avg - precipitation4_avg_ctl) - (precipitation2_avg - precipitation2_avg_ctl)) - (precipitation3_avg - precipitation3_avg_ctl - (precipitation1_avg  - precipitation1_avg_ctl))

lats=array.lat
lons=array.lon

axes[2].set_title('(c) 2C - AM (50%cond - bucket)', size = med)
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

m.contour(xi,yi,landmask1, 1)
m.contour(xi,yi,landmask2, 1, linestyles = 'dotted')

# Add Colorbar
cbar = fig.colorbar(cs, orientation = 'vertical', ax = axes, shrink = 0.5) # usually on right 
cbar.set_label('mm/d', size=med)
cbar.ax.tick_params(labelsize=med)

fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/'+dire+'/Pavg_minus_ctl_bucket_vs_VP05_2C-AM_40S40N_120-480_paper.png', bbox_inches='tight', dpi=100)
fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/'+dire+'/Pavg_minus_ctl_bucket_vs_VP05_2C-AM_40S40N_120-480_paper.pdf', bbox_inches='tight', dpi=400)
fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/'+dire+'/Pavg_minus_ctl_bucket_vs_VP05_2C-AM_40S40N_120-480_paper.eps', bbox_inches='tight', dpi=600)


# fig.savefig('/scratch/mp586/Code/Graphics/Isca/full_continents_newbucket_fixedSSTs_zonally_symmetric_vegetation_vegpref05_plus_2pt52K_and_2xCO2_spinup_361_witholr/Pavg_minus_ctl_bucket_vs_VP05_P-Econts_40S40N_ctl_25-121_pert_24-120.png', bbox_inches='tight', dpi=100)
# fig.savefig('/scratch/mp586/Code/Graphics/Isca/full_continents_newbucket_fixedSSTs_zonally_symmetric_vegetation_vegpref05_plus_2pt52K_and_2xCO2_spinup_361_witholr/Pavg_minus_ctl_bucket_vs_VP05_P-Econts_40S40N_ctl_25-121_pert_24-120.svg', bbox_inches='tight', dpi=100)

plt.close('all')

landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]
# for specified lats
landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it


dP_SB =  (precipitation1_avg - precipitation1_avg_ctl).where(landmask == 1).sel(lat = slice(-30.,30.))
dP_VP05 = (precipitation2_avg - precipitation2_avg_ctl).where(landmask == 1).sel(lat = slice(-30.,30.))

dP_SB = np.asarray(dP_SB).flatten()
dP_VP05 = np.asarray(dP_VP05).flatten()


fig, ax = plt.subplots()
ax.plot(dP_SB, dP_VP05, 'k.')
ax.set_xlim(-6.,6.)
ax.set_ylim(-6.,6.)
ax.set_ylabel('$\Delta P$ (50%cond) in mm/d')
ax.set_xlabel('$\Delta P$ (bucket) in mm/d')

mask = ~np.isnan(dP_SB)

[k,dy,r,p,stderr] = linreg(dP_SB[mask],dP_VP05[mask]) # aa = 8.4, dq = -32
x1 = np.linspace(np.min(dP_SB[mask]),np.max(dP_SB[mask]),500)
y = k*x1 + dy
ax.plot(x1,y,'k-')
ax.annotate('r = '+str("%.2f" % r)+', k = '+str("%.2f" % k)+', dy = '+str("%.2f" % dy), xy=(0.05,0.05), xycoords='axes fraction')

fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/'+dire+'/dP_simple_bucket_versus_dP_VP05.png', bbox_inches='tight', dpi=100)
