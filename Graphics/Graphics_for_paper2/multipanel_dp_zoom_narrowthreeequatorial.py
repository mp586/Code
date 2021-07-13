from netCDF4 import Dataset
import numpy as np
import matplotlib as mpl 
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import xarray as xr
import pandas as pd
import os
from scipy.odr import *

import sys
sys.path.insert(0, '/scratch/mp586/Code/PYCODES')
from plotting_routines_kav7 import *
import stats as st

GFDL_BASE = os.environ['GFDL_BASE']
sys.path.insert(0, os.path.join(GFDL_BASE,'src/extra/python/scripts'))
import cell_area as ca





#mpl.rcParams["lines.linewidth"] = 0.5 # setting linewidth for landmask contour plot doesn't work otherwise 

variable = 'precipitation' # 'bucket_depth' # 'precipitation' # 'flux_lhe' # precipitation # t_surf
colormap = 'BrBG' # 'BrBG' # 'BrBG' # 'RdBu_r'
minval = -2. # -0.05 # -2. # -2. # -10.
maxval = 2. # 0.05 #  2. # 2. # 10. 
units ='mm/d' # 'm' # 'mm/d' # 'mm/d' # 'K'
factor = 86400.# 1. # 86400. # 1./28. # 86400 # 1.
low_lim = -4.# -0.1 # -4. # 0.
up_lim = 4.  # 0.1 # 4. # 5. 


control_sb_dirs = [
'narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387',
'narrow_three_equatorial_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387'
]

control_vp0_dirs = [
'narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref0_commit7bb4387',
'x']

control_vp02_dirs = [
'x',
'x'
]

control_vp05_dirs = [
'narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_commit7bb4387',
'narrow_three_equatorial_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_commit7bb4387'
]


vp0_dirs = [
'narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref0_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'x']

vp05_dirs = [
'narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'narrow_three_equatorial_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387'
]


simple_bucket_dirs = [
'narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'narrow_three_equatorial_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387'
]


ctl_dict = {'vp05_ctl': control_vp05_dirs,'sb_ctl': control_sb_dirs}

ctl_list = ['vp05_ctl','sb_ctl']

[precipitation_ctl,precipitation_avg_ctl,x,x,x]=seasonal_surface_variable('Isca_DATA/ISCA_HPC/narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387','isca',1,10,variable,units, factor=factor)


lats = precipitation_avg_ctl.lat
lons = precipitation_avg_ctl.lon


landmasks = ['narrow_three','narrow_three_equatorial']
landfile=Dataset(os.path.join(GFDL_BASE,'input/'+landmasks[0]+'/land.nc'),mode='r')
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]

landmask_array = np.zeros((len(landmasks), len(landlats), len(landlons)))

for i in range(len(landmasks)):
    landfile=Dataset(os.path.join(GFDL_BASE,'input/'+landmasks[i]+'/land.nc'),mode='r')
    landmask_array[i,:,:]=landfile.variables['land_mask'][:]




pert_dict = {'vp05': vp05_dirs,
'sb': simple_bucket_dirs}

pert_list = ['vp05','sb']

precip_change_matrix = np.zeros((len(vp05_dirs),len(pert_dict),len(lats),len(lons)))
precip_ctl_matrix = np.zeros((len(vp05_dirs),len(pert_dict),len(lats),len(lons)))
precip_pert_matrix = np.zeros((len(vp05_dirs),len(pert_dict),len(lats),len(lons)))

ucomp_change_matrix = np.zeros((len(vp05_dirs),len(pert_dict),len(lats),len(lons)))
vcomp_change_matrix = np.zeros((len(vp05_dirs),len(pert_dict),len(lats),len(lons)))


for i in range(len(control_sb_dirs)):
    [precipitation_ctl,precipitation_avg_ctl,x,x,x]=seasonal_surface_variable('Isca_DATA/ISCA_HPC/'+control_sb_dirs[i],'isca',121,481,variable,units, factor=factor)
    [x,ucomp_avg_ctl,x,x,x]=seasonal_surface_variable_interp('Isca_DATA/ISCA_HPC/'+control_sb_dirs[i],'isca',121,481,'ucomp','m/s', factor=1., level=2)
    [x,vcomp_avg_ctl,x,x,x]=seasonal_surface_variable_interp('Isca_DATA/ISCA_HPC/'+control_sb_dirs[i],'isca',121,481,'vcomp','m/s', factor=1., level=2)
    
    for j in range(len(pert_list)):
        testdir = pert_dict[pert_list[j]][i]
        if testdir != 'x':
            testdir = 'Isca_DATA/ISCA_HPC/'+pert_dict[pert_list[j]][i]
            [precipitation,precipitation_avg,x,x,x]=seasonal_surface_variable(testdir,'isca',120,480,variable,units, factor=factor)
            precip_change_matrix[i,j,:,:] = precipitation_avg - precipitation_avg_ctl
            precip_ctl_matrix[i,j,:,:] = precipitation_avg_ctl
            precip_pert_matrix[i,j,:,:] = precipitation_avg
            [x,ucomp_avg,x,x,x]=seasonal_surface_variable(testdir,'isca',120,480,'ucomp','m/s', factor=1.,level=2)
            ucomp_change_matrix[i,j,:,:] = ucomp_avg - ucomp_avg_ctl
            [x,vcomp_avg,x,x,x]=seasonal_surface_variable(testdir,'isca',120,480,'vcomp','m/s', factor=1.,level=2)
            vcomp_change_matrix[i,j,:,:] = vcomp_avg - vcomp_avg_ctl

        
small = 22 #largefonts 14 # smallfonts 10 # medfonts = 14
med = 24 #largefonts 18 # smallfonts 14 # medfonts = 16
lge = 26 #largefonts 22 # smallfonts 18 # medfonts = 20


names = ['$\Delta P_{50\%cond}$', '$\Delta P_{100\%cond}$']
conts = ['6$^{\circ}$ lon','6$^{\circ}$ lon, 5$^{\circ}$S-5$^{\circ}$N']

fig = plt.figure(figsize = (15,8))

m = Basemap(projection='cyl',resolution='c', llcrnrlat=-40, urcrnrlat=40,llcrnrlon=-30, urcrnrlon=170)

v = np.linspace(minval,maxval,41) # , endpoint=True)


for i in range(len(control_sb_dirs)):
    for j in range(len(pert_list)):
        testdir = pert_dict[pert_list[j]][i]
        ax = plt.subplot2grid((len(vp0_dirs),len(pert_dict)), (i,j))
        if i == 0:
            ax.set_title(names[j], size = med)
        if j == 0:
            ax.set_ylabel(conts[i], size = med)
        if testdir == 'x':
            ax.xaxis.set_visible(False)
            # make spines (the box) invisible
            plt.setp(ax.spines.values(), visible=False)
            # remove ticks and labels for the left axis
            ax.tick_params(left=False, labelleft=False)
            #remove background patch (only needed for non-white background)
            # ax.patch.set_visible(False)
        else:
            array = xr.DataArray(precip_change_matrix[i,j,:,:],coords=[lats,lons],dims=['lat','lon'])
            array = np.asarray(array)
            array, lons_cyclic = addcyclic(array, lons)
            array,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))

            uarray = xr.DataArray(ucomp_change_matrix[i,j,:,:],coords=[lats,lons],dims=['lat','lon'])
            uarray = np.asarray(uarray)
            uarray, lons_cyclic = addcyclic(uarray, lons)
            uarray,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,uarray,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))

            varray = xr.DataArray(vcomp_change_matrix[i,j,:,:],coords=[lats,lons],dims=['lat','lon'])
            varray = np.asarray(varray)
            varray, lons_cyclic = addcyclic(varray, lons)
            varray,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,varray,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))


            lon, lat = np.meshgrid(lons_cyclic, lats)
            xi, yi = m(lon, lat)

            array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])

            cs = m.contourf(xi,yi,array, v, cmap=colormap, extend = 'both')

            Q = ax.quiver(xi[::4,::4], yi[::4,::4], uarray[::4,::4], varray[::4,::4], scale=10, units='inches')


            landmask,landlons_shift = shiftgrid(np.max(landlons)-180.,landmask_array[i,:,:],landlons,start=False,cyclic=np.max(landlons))
            landmask,lons_cyclic = addcyclic(landmask, landlons_shift)
            m.contour(xi,yi,landmask, 1, colors = 'k', linewidths = 1.5)
            m.drawparallels(np.arange(-10.,20.,10.),labels=[], fontsize=small)


# plt.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0.02, hspace=0.02)
# cb_ax = plt.axes([0.83, 0.3, 0.01, 0.4])
# # add an axes, lower left corner in [0.83, 0.3] measured in figure coordinate with axes width 0.01 and height 0.4
# cbar = plt.colorbar(cs, cax=cb_ax)

plt.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.9, wspace=0.02, hspace=0.02)
cb_ax = plt.axes([0.3, 0.05, 0.4, 0.03])
cbar = plt.colorbar(cs, cax=cb_ax, orientation = 'horizontal')
cbar.set_label(units, size = med)
cbar.ax.tick_params(labelsize= 18)

qk = plt.quiverkey(Q, 0.8, 0.05, 5., '5 '+r'$\frac{m}{s}$', coordinates='figure', labelpos = 'E', fontproperties={'size': small})


plt.savefig('/scratch/mp586/Code/Graphics/multipanel_zoombucket_narrowthreeequatorial_avg_minus_ctl_120-480_lowcbar_grid.png', bbox_inches = 'tight', format = 'png', dpi = 400)
plt.savefig('/scratch/mp586/Code/Graphics/multipanel_zoombucket_narrowthreeequatorial_avg_minus_ctl_120-480_lowcbar_grid.pdf', bbox_inches = 'tight', format = 'pdf', dpi = 400)
plt.savefig('/scratch/mp586/Code/Graphics/multipanel_zoombucket_narrowthreeequatorial_avg_minus_ctl_120-480_lowcbar_grid.eps', bbox_inches = 'tight', format = 'eps', dpi = 600)


plt.close()



