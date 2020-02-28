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

# nc = Dataset('/scratch/mp586/Isca/input/sst_clim_amip_zonalsymm.nc')
nc = Dataset('/scratch/mp586/Isca/input/sst.nc')

ssts = xr.DataArray(nc.variables['sst_clim_amip_zonalsymm'][:]) - 273.15
lats = nc.variables['lat'][:]
lons = nc.variables['lon'][:]

latr = np.deg2rad(lats)
cos_1d = np.cos(latr)
cos_2d = np.expand_dims(cos_1d, axis = 1) # for lons
cos_2d = np.repeat(cos_2d, len(lons), axis = 1)	
cos_2d = xr.DataArray(cos_2d, coords=[lats,lons], dims = ['lat','lon'])

sstsm = ssts.mean(dim='dim_0')
ma =  np.ma.MaskedArray(sstsm, mask=np.isnan(sstsm))
w_avg = np.average(ma, axis=1,weights=cos_2d)


small = 18
med = 22
lge = 26


fig = plt.figure(figsize = (25,10))

ax1 = plt.subplot2grid((5,8), (0,1), colspan = 5, rowspan = 3)


m = Basemap(projection='kav7',lon_0=0.,llcrnrlon=-180.,llcrnrlat=-30.,urcrnrlon=180.,urcrnrlat=30.,resolution='c')
#    m = Basemap(projection='cyl',llcrnrlon=-180.,llcrnrlat=-30.,urcrnrlon=180.,urcrnrlat=30.,resolution='c') works, but not with kav7 projection
array = xr.DataArray(sstsm,coords=[lats,lons],dims=['lat','lon'])

zonavg_thin = w_avg

array = np.asarray(array) #- This line fixes the problem!
#the following line caused DataType error all of a sudden... Doesnt' accept xarray as input array for shiftgrid anymore.

array, lons_cyclic = addcyclic(array, lons)
array,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))

array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])

# shiftgrid lons0 (first entry) is the longitude in lons which should correspond with the center longitude of the map. start = False --> lonsin is the end latitude, not the beginning.
# this doesn't work for some reason
#array, lons = shiftgrid(np.max(lons)-100.,array,lons,start=True,cyclic=np.max(lons))

lons = lons_cyclic

m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0], fontsize=small)
m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1], fontsize=small)

lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

v = np.linspace(-30.,30.,21) # , endpoint=True)

cs = m.contourf(xi,yi,array, v, cmap='RdBu_r', extend = 'both')
#cs = m.contourf(xi,yi,array, norm=MidpointNormalize(midpoint=0.), cmap=plt.cm.RdBu_r,vmin = -10., vmax = 40.)

# Add Colorbar
cbar = m.colorbar(cs, location='right', pad="10%") # usually on right 
cbar.set_label('$^{\circ}$C', size=med)
cbar.ax.tick_params(labelsize=med) 
cbar.set_clim(-30., 30.)


ax2 = plt.subplot2grid((5,8), (0,6), rowspan = 3)
plt.plot(zonavg_thin,lats, 'k-')
plt.ylabel('Latitude', size=med)
plt.xlabel('$^{\circ}$C', size=med)
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position('right')
ax2.tick_params(axis='both', which='major', labelsize=small)
ax2.invert_xaxis()

plt.savefig('/scratch/mp586/Code/Graphics/sst_clim_amip_zonalsymm_plot.png', format = 'png', bbox_inches='tight')
plt.savefig('/scratch/mp586/Code/Graphics/sst_clim_amip_zonalsymm_plot.pdf', format = 'pdf', bbox_inches='tight')