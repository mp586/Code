# from netCDF4 import Dataset
# import numpy as np
# from matplotlib import pyplot as plt
# from mpl_toolkits.basemap import Basemap, cm
# import xarray as xr
# import pandas as pd
# import os

# import sys
# sys.path.insert(0, '/scratch/mp586/Code/PYCODES')
# from plotting_routines_kav7 import *
# import stats as st

# GFDL_BASE = os.environ['GFDL_BASE']
# sys.path.insert(0, os.path.join(GFDL_BASE,'src/extra/python/scripts'))
# import cell_area as ca

# ctl_model = input('Enter model 1 name as string ')
# if (ctl_model == 'Isca') or (ctl_model == 'isca'): 
#     control_model = 'Isca_DATA'
# elif (ctl_model == 'gfdl') or (ctl_model == 'GFDL'):
#     control_model = 'GFDL_DATA'

# HPC = input('Data in ISCA_HPC ? Yes or No? ')
# control_dir = input('Enter control directory from exp 1 as string ')
# if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
#     control_dir= control_model + '/ISCA_HPC/' + control_dir
# else: 
#     control_dir= control_model + '/' + control_dir

# #print control_dir
# ctl_runmin=input('Enter runmin number for exp 1 ')  # Should be a January month for seasonal variables to be correct
# ctl_runmax=input('Enter runmax number for comparison 1 ')
# ctl_timeseries_max = input('Enter end of ctl timeseries month for exp 1 ')

# model = input('Enter model exp 1 ')
# if (model == 'Isca') or (model == 'isca'): 
#     model_data = 'Isca_DATA'
#     output_dir1 = 'Isca'
# elif (model == 'gfdl') or (model == 'GFDL'):
#     model_data = 'GFDL_DATA'
#     output_dir1 = ''

# HPC = input('Data in ISCA_HPC ? Yes or No? ')
# testdir_in1= input('Enter perturbed data directory name as string for exp 1 ')
# runmin=input('Enter runmin number for exp 1 ')  # Should be a January month for seasonal variables to be correct
# runmax=input('Enter runmax number for exp 1 ')
# if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
#     exp1_name = 'ISCA_HPC_'+testdir_in1
#     testdir = model_data + '/ISCA_HPC/' + testdir_in1
#     testdir_in1 = '/ISCA_HPC/' + testdir_in1
# else: 
#     exp1_name = testdir_in1
#     testdir = model_data + '/' + testdir_in1

# land = input('Which landmask? ')
# landfile=Dataset(os.path.join(GFDL_BASE,'input/'+land+'/land.nc'),mode='r')

# landmask=landfile.variables['land_mask'][:]
# landlats=landfile.variables['lat'][:]
# landlons=landfile.variables['lon'][:]
# # for specified lats
# landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

# level = input('Which Level? ')


# area_array, dx, dy = ca.cell_area_all(t_res=42,base_dir='/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/') # added _all because then dx and dy are also returned 
# area_array = xr.DataArray(area_array) # returned in units of m bzw m^2, because radius in cell_area.py is given in metres

# area_array_3D = np.expand_dims(area_array, axis=0)
# area_array_3D = np.repeat(area_array_3D, 40, axis = 0) # to make area_array 3D (pressure, lat, lon)

# [ucomp1,ucomp1_avg,ucomp1_seasonal_avg,ucomp1_month_avg,ucomp1_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'ucomp','m/s')
# [ucomp1_ctl,ucomp1_avg_ctl,ucomp1_seasonal_avg_ctl,ucomp1_month_avg_ctl,ucomp1_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'ucomp','m/s')
# [vcomp1,vcomp1_avg,vcomp1_seasonal_avg,vcomp1_month_avg,vcomp1_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'vcomp','m/s')
# [vcomp1_ctl,vcomp1_avg_ctl,vcomp1_seasonal_avg_ctl,vcomp1_month_avg_ctl,vcomp1_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'vcomp','m/s')
# [gph1,gph1_avg,gph1_seasonal_avg,gph1_month_avg,gph1_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'height','m')
# [gph1_ctl,gph1_avg_ctl,gph1_seasonal_avg_ctl,gph1_month_avg_ctl,gph1_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'height','m')


# ################ read in data from exp 2 ###############################

# ctl_model = input('Enter model 1 name as string ')
# if (ctl_model == 'Isca') or (ctl_model == 'isca'): 
#     control_model = 'Isca_DATA'
# elif (ctl_model == 'gfdl') or (ctl_model == 'GFDL'):
#     control_model = 'GFDL_DATA'

# HPC = input('Data in ISCA_HPC ? Yes or No? ')
# control_dir = input('Enter control directory from exp 1 as string ')
# if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
#     control_dir= control_model + '/ISCA_HPC/' + control_dir
# else: 
#     control_dir= control_model + '/' + control_dir

# #print control_dir
# ctl_runmin=input('Enter runmin number for exp 1 ')  # Should be a January month for seasonal variables to be correct
# ctl_runmax=input('Enter runmax number for comparison 1 ')
# ctl_timeseries_max = input('Enter end of ctl timeseries month for exp 1 ')

# model = input('Enter model exp 1 ')
# if (model == 'Isca') or (model == 'isca'): 
#     model_data = 'Isca_DATA'
#     output_dir1 = 'Isca'
# elif (model == 'gfdl') or (model == 'GFDL'):
#     model_data = 'GFDL_DATA'
#     output_dir1 = ''

# HPC = input('Data in ISCA_HPC ? Yes or No? ')
# testdir_in1= input('Enter perturbed data directory name as string for exp 1 ')
# runmin=input('Enter runmin number for exp 1 ')  # Should be a January month for seasonal variables to be correct
# runmax=input('Enter runmax number for exp 1 ')
# if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
#     exp1_name = 'ISCA_HPC/'+testdir_in1
#     testdir = model_data + '/ISCA_HPC/' + testdir_in1
#     testdir_in1 = '/ISCA_HPC/' + testdir_in1
# else: 
#     exp1_name = testdir_in1
#     testdir = model_data + '/' + testdir_in1

# land = input('Which landmask? ')
# landfile=Dataset(os.path.join(GFDL_BASE,'input/'+land+'/land.nc'),mode='r')

# landmask=landfile.variables['land_mask'][:]
# landlats=landfile.variables['lat'][:]
# landlons=landfile.variables['lon'][:]
# # for specified lats
# landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

# level = input('Which Level? ')



# [ucomp2,ucomp2_avg,ucomp2_seasonal_avg,ucomp2_month_avg,ucomp2_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'ucomp','m/s')
# [ucomp2_ctl,ucomp2_avg_ctl,ucomp2_seasonal_avg_ctl,ucomp2_month_avg_ctl,ucomp2_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'ucomp','m/s')
# [vcomp2,vcomp2_avg,vcomp2_seasonal_avg,vcomp2_month_avg,vcomp2_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'vcomp','m/s')
# [vcomp2_ctl,vcomp2_avg_ctl,vcomp2_seasonal_avg_ctl,vcomp2_month_avg_ctl,vcomp2_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'vcomp','m/s')
# [gph2,gph2_avg,gph2_seasonal_avg,gph2_month_avg,gph2_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'height','m')
# [gph2_ctl,gph2_avg_ctl,gph2_seasonal_avg_ctl,gph2_month_avg_ctl,gph2_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'height','m')



# # 200hpa
# winds_one_level(outdir,runmin,runmax,'200hPa_winds_gph_diff_twoCs-America_climatechange',((ucomp2_avg - ucomp2_avg_ctl) -(ucomp1_avg - ucomp1_avg_ctl))[15,:,:],
# 	((vcomp2_avg - vcomp2_avg_ctl) - (vcomp1_avg - vcomp1_avg_ctl))[15,:,:],((gph2_avg - gph2_avg_ctl) - (gph1_avg - gph1_avg_ctl))[15,:,:],'slp','m',-30.,30.,
# 	landmaskxr,veclen=5,level=15,units_numerator = 'm', units_denom = 's', nmb_contours = [-15.,15.,30.])

# # 870 hpa 
# winds_one_level(outdir,runmin,runmax,'870hPa_winds_gph_diff_twoCs-America_climatechange',((ucomp2_avg - ucomp2_avg_ctl) -(ucomp1_avg - ucomp1_avg_ctl))[37,:,:],
# 	((vcomp2_avg - vcomp2_avg_ctl) - (vcomp1_avg - vcomp1_avg_ctl))[37,:,:],((gph2_avg - gph2_avg_ctl) - (gph1_avg - gph1_avg_ctl))[37,:,:],'slp','m',-5.,5.,
# 	landmaskxr,veclen=5,level=37,units_numerator = 'm', units_denom = 's', nmb_contours = [-2.,2.,5.])




########################### Plotting Geopotential Height for just one experiment, but as panel plot for 200 and 850 hPa ###############################################################

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

HPC = 'no'
control_dir = 'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm'
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
testdir_in1= 'two_continents_newbucket_fixedSSTs_from_realworld_zonallysymm'
runmin=120
runmax=480
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    exp1_name = 'ISCA_HPC/'+testdir_in1
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

[ucomp1,ucomp1_avg,ucomp1_seasonal_avg,ucomp1_month_avg,ucomp1_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'ucomp','m/s')
[ucomp1_ctl,ucomp1_avg_ctl,ucomp1_seasonal_avg_ctl,ucomp1_month_avg_ctl,ucomp1_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'ucomp','m/s')
[vcomp1,vcomp1_avg,vcomp1_seasonal_avg,vcomp1_month_avg,vcomp1_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'vcomp','m/s')
[vcomp1_ctl,vcomp1_avg_ctl,vcomp1_seasonal_avg_ctl,vcomp1_month_avg_ctl,vcomp1_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'vcomp','m/s')
[gph1,gph1_avg,gph1_seasonal_avg,gph1_month_avg,gph1_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'height','m')
[gph1_ctl,gph1_avg_ctl,gph1_seasonal_avg_ctl,gph1_month_avg_ctl,gph1_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'height','m')

## plotting ## 

small = 18 #largefonts 14 # smallfonts 10 # medfonts = 14
med = 20 #largefonts 18 # smallfonts 14 # medfonts = 16
lge = 22 #largefonts 22 # smallfonts 18 # medfonts = 20

landlats = np.asarray(landmaskxr.lat)
landlons = np.asarray(landmaskxr.lon)
landmask = np.asarray(landmaskxr)

landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))
landmask, landlons = addcyclic(landmask, landlons)


land2 = 'square_Africa'
landfile2=Dataset(os.path.join(GFDL_BASE,'input/'+land+'/land.nc'),mode='r')

landmask2=landfile2.variables['land_mask'][:]
landlats2=landfile2.variables['lat'][:]
landlons2=landfile2.variables['lon'][:]
# for specified lats
landmaskxr2=xr.DataArray(landmask2,coords=[landlats2,landlons2],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

landlats2 = np.asarray(landmaskxr2.lat)
landlons2 = np.asarray(landmaskxr2.lon)
landmask2 = np.asarray(landmaskxr2)

landmask2,landlons2 = shiftgrid(np.max(landlons2)-180.,landmask2,landlons2,start=False,cyclic=np.max(landlons2))
landmask2, landlons2 = addcyclic(landmask2, landlons2)

veclen=10
units_numerator = 'm'
units_denom = 's'

fig, ax = plt.subplots(1, 2, figsize = (25,10))

# if fig, ax = plt.subplots(0,0, figsize = (25,10)) then ax is a numpy array --> can't do quiver plot on it
level = [37, 15]
minmax = [-3.,3.,-3.,3.]
nmb_contours = [5,20]
plt_title = ['870 hPa', '200hPa']
j = 0 
for i in range(2):

    m = Basemap(projection='kav7',lon_0=0.,resolution='c', ax = ax[i])

    uwind = (ucomp1_avg - ucomp1_avg_ctl)[level[i],:,:]
    vwind = (vcomp1_avg - vcomp1_avg_ctl)[level[i],:,:]
    array = ((gph1_avg - gph1_avg_ctl)/gph1_avg_ctl)[level[i],:,:]*100. # change in gph in percent
    lons = uwind.lon
    lats = uwind.lat
    uwind, lons_cyclic = addcyclic(uwind, lons)
    vwind, lons_cyclic = addcyclic(vwind, lons)
    uwind = np.asarray(uwind)
    vwind = np.asarray(vwind)
    uwind,lons_shift = shiftgrid(np.max(lons_cyclic)-180.,uwind,lons_cyclic,start=False,
                   cyclic=np.max(lons_cyclic))
    vwind,lons_shift = shiftgrid(np.max(lons_cyclic)-180.,vwind,lons_cyclic,start=False,
                   cyclic=np.max(lons_cyclic))  

    m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0], fontsize=med)

    array, lons_cyclic = addcyclic(array, lons)
    array = np.asarray(array)
    array, lons_shift = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,
                     start=False,cyclic=np.max(lons_cyclic))
    array = xr.DataArray(array,coords=[lats,lons_shift],dims=['lat','lon'])



    lon, lat = np.meshgrid(lons_shift, lats)
    xi, yi = m(lon, lat)
    m.contour(xi,yi,landmask, 1)
    m.contour(xi,yi,landmask2, 1, linestyles = 'dotted')

    v = np.linspace(minmax[j], minmax[j+1],21) # , endpoint=True)

    cs = ax[i].contourf(xi,yi,array, v, cmap=plt.cm.coolwarm, extend = 'both')

#    if nmb_contours[i] != 0:  # add contours 
#        cont = ax[i].contour(xi,yi,array,nmb_contours[i], colors = 'k', linewidth=5)

    # Add Colorbar
    cbar = fig.colorbar(cs, ax = ax[i], orientation = 'horizontal')
    cbar.set_label('% (m)', size=med)
    cbar.ax.tick_params(labelsize=med) 
    cbar.set_clim(minmax[j], minmax[j+1])

    Q = ax[i].quiver(xi[::3,::3], yi[::3,::3], uwind[::3,::3], vwind[::3,::3], scale = 30, scale_units = 'inches')
    j += 2

    ax[i].set_title(plt_title[i], fontsize = med)

qk = plt.quiverkey(Q, 0.83, 0.83, veclen, str(veclen)+r'$\frac{'+units_numerator+'}{'+units_denom+'}$', labelpos='E', coordinates='figure', fontproperties={'size': med})

fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/two_continents_newbucket_fixedSSTs_from_realworld_zonallysymm/gph_2C-AM_control.png')



##########################################################################################################################################################


########################### Plotting Geopotential Height comparing two experiments, as panel plot for 200 and 850 hPa ###############################################################

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

HPC = 'no'
control_dir = 'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm'
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    control_dir= control_model + '/ISCA_HPC/' + control_dir
else: 
    control_dir= control_model + '/' + control_dir

#print control_dir
ctl_runmin=121 # Should be a January month for seasonal variables to be correct
ctl_runmax=481
ctl_timeseries_max = 361

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

level = 39


area_array, dx, dy = ca.cell_area_all(t_res=42,base_dir='/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/') # added _all because then dx and dy are also returned 
area_array = xr.DataArray(area_array) # returned in units of m bzw m^2, because radius in cell_area.py is given in metres

area_array_3D = np.expand_dims(area_array, axis=0)
area_array_3D = np.repeat(area_array_3D, 40, axis = 0) # to make area_array 3D (pressure, lat, lon)

[ucomp1,ucomp1_avg,ucomp1_seasonal_avg,ucomp1_month_avg,ucomp1_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'ucomp','m/s')
[ucomp1_ctl,ucomp1_avg_ctl,ucomp1_seasonal_avg_ctl,ucomp1_month_avg_ctl,ucomp1_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'ucomp','m/s')
[vcomp1,vcomp1_avg,vcomp1_seasonal_avg,vcomp1_month_avg,vcomp1_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'vcomp','m/s')
[vcomp1_ctl,vcomp1_avg_ctl,vcomp1_seasonal_avg_ctl,vcomp1_month_avg_ctl,vcomp1_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'vcomp','m/s')
[gph1,gph1_avg,gph1_seasonal_avg,gph1_month_avg,gph1_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'height','m')
[gph1_ctl,gph1_avg_ctl,gph1_seasonal_avg_ctl,gph1_month_avg_ctl,gph1_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'height','m')


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
ctl_timeseries_max = 361

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

outdir = testdir_in1


land = 'two_continents'
landfile=Dataset(os.path.join(GFDL_BASE,'input/'+land+'/land.nc'),mode='r')

landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]
# for specified lats
landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

level = 39



[ucomp2,ucomp2_avg,ucomp2_seasonal_avg,ucomp2_month_avg,ucomp2_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'ucomp','m/s')
[ucomp2_ctl,ucomp2_avg_ctl,ucomp2_seasonal_avg_ctl,ucomp2_month_avg_ctl,ucomp2_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'ucomp','m/s')
[vcomp2,vcomp2_avg,vcomp2_seasonal_avg,vcomp2_month_avg,vcomp2_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'vcomp','m/s')
[vcomp2_ctl,vcomp2_avg_ctl,vcomp2_seasonal_avg_ctl,vcomp2_month_avg_ctl,vcomp2_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'vcomp','m/s')
[gph2,gph2_avg,gph2_seasonal_avg,gph2_month_avg,gph2_annual_avg,time]=seasonal_4D_variable(testdir,model,runmin,runmax,'height','m')
[gph2_ctl,gph2_avg_ctl,gph2_seasonal_avg_ctl,gph2_month_avg_ctl,gph2_annual_avg_ctl,time]=seasonal_4D_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'height','m')


## plotting ## 

small = 18 #largefonts 14 # smallfonts 10 # medfonts = 14
med = 20 #largefonts 18 # smallfonts 14 # medfonts = 16
lge = 22 #largefonts 22 # smallfonts 18 # medfonts = 20

landlats = np.asarray(landmaskxr.lat)
landlons = np.asarray(landmaskxr.lon)
landmask = np.asarray(landmaskxr)

landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))
landmask, landlons = addcyclic(landmask, landlons)

land1 = 'square_South_America'
landfile1=Dataset(os.path.join(GFDL_BASE,'input/'+land1+'/land.nc'),mode='r')

landmask1=landfile1.variables['land_mask'][:]
landlats1=landfile1.variables['lat'][:]
landlons1=landfile1.variables['lon'][:]
# for specified lats
landmaskxr1=xr.DataArray(landmask1,coords=[landlats1,landlons1],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

landlats1 = np.asarray(landmaskxr1.lat)
landlons1 = np.asarray(landmaskxr1.lon)
landmask1 = np.asarray(landmaskxr1)

landmask1,landlons1 = shiftgrid(np.max(landlons1)-180.,landmask1,landlons1,start=False,cyclic=np.max(landlons1))
landmask1, landlons1 = addcyclic(landmask1, landlons1)

land2 = 'square_Africa'
landfile2=Dataset(os.path.join(GFDL_BASE,'input/'+land2+'/land.nc'),mode='r')

landmask2=landfile2.variables['land_mask'][:]
landlats2=landfile2.variables['lat'][:]
landlons2=landfile2.variables['lon'][:]
# for specified lats
landmaskxr2=xr.DataArray(landmask2,coords=[landlats2,landlons2],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

landlats2 = np.asarray(landmaskxr2.lat)
landlons2 = np.asarray(landmaskxr2.lon)
landmask2 = np.asarray(landmaskxr2)

landmask2,landlons2 = shiftgrid(np.max(landlons2)-180.,landmask2,landlons2,start=False,cyclic=np.max(landlons2))
landmask2, landlons2 = addcyclic(landmask2, landlons2)

veclen=5
units_numerator = 'm'
units_denom = 's'

fig, ax = plt.subplots(1, 2, figsize = (25,10))

# if fig, ax = plt.subplots(0,0, figsize = (25,10)) then ax is a numpy array --> can't do quiver plot on it
level = [37, 15]
minmax = [-5.,5.,-30.,30.]
nmb_contours = 10
plt_title = ['a) 870 hPa', 'b) 200hPa']
j = 0 
for i in range(2):

    m = Basemap(projection='kav7',lon_0=0.,resolution='c', ax = ax[i])

    uwind = ((ucomp2_avg - ucomp2_avg_ctl) - (ucomp1_avg - ucomp1_avg_ctl))[level[i],:,:]
    vwind = ((vcomp2_avg - vcomp2_avg_ctl) - (vcomp1_avg - vcomp1_avg_ctl))[level[i],:,:]
    array = ((((gph2_avg - gph2_avg_ctl) - (gph1_avg - gph1_avg_ctl))))[level[i],:,:] #*100. # change in gph in percent /(gph2_avg_ctl - gph1_avg_ctl)
    lons = uwind.lon
    lats = uwind.lat
    uwind, lons_cyclic = addcyclic(uwind, lons)
    vwind, lons_cyclic = addcyclic(vwind, lons)
    uwind = np.asarray(uwind)
    vwind = np.asarray(vwind)
    uwind,lons_shift = shiftgrid(np.max(lons_cyclic)-180.,uwind,lons_cyclic,start=False,
                   cyclic=np.max(lons_cyclic))
    vwind,lons_shift = shiftgrid(np.max(lons_cyclic)-180.,vwind,lons_cyclic,start=False,
                   cyclic=np.max(lons_cyclic))  

    m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0], fontsize=med)

    array, lons_cyclic = addcyclic(array, lons)
    array = np.asarray(array)
    array, lons_shift = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,
                     start=False,cyclic=np.max(lons_cyclic))
    array = xr.DataArray(array,coords=[lats,lons_shift],dims=['lat','lon'])



    lon, lat = np.meshgrid(lons_shift, lats)
    xi, yi = m(lon, lat)
    m.contour(xi,yi,landmask1, 1, linewidths = 2., colors = 'w')
    m.contour(xi,yi,landmask2, 1, linestyles = 'dashed', linewidths = 2., colors = 'w')

    v = np.linspace(minmax[j], minmax[j+1],21) # , endpoint=True)

    cs = ax[i].contourf(xi,yi,array, v, cmap=plt.cm.coolwarm, extend = 'both')

    if nmb_contours != 0:  # add contours 
        cont = ax[i].contour(xi,yi,array,nmb_contours, colors = 'dimgray', linewidths=1.)

    # Add Colorbar
    cbar = fig.colorbar(cs, ax = ax[i], orientation = 'horizontal')
    cbar.set_label('$\Delta gph$ (m)', size=med)
    cbar.ax.tick_params(labelsize=med) 
    cbar.set_clim(minmax[j], minmax[j+1])

    Q = ax[i].quiver(xi[::3,::3], yi[::3,::3], uwind[::3,::3], vwind[::3,::3], scale = 10, scale_units = 'inches')
    j += 2

    ax[i].set_title(plt_title[i], fontsize = med)

qk = plt.quiverkey(Q, 0.87, 0.83, veclen, str(veclen)+r'$\frac{'+units_numerator+'}{'+units_denom+'}$', labelpos='E', coordinates='figure', fontproperties={'size': med})

fig.savefig('/scratch/mp586/Code/Graphics/Isca'+outdir+'/gph_change_epx2_minus_exp1.png', bbox_inches = 'tight')
fig.savefig('/scratch/mp586/Code/Graphics/Isca'+outdir+'/gph_change_epx2_minus_exp1.svg', bbox_inches = 'tight')
fig.savefig('/scratch/mp586/Code/Graphics/Isca'+outdir+'/gph_change_epx2_minus_exp1.pdf', bbox_inches = 'tight')

