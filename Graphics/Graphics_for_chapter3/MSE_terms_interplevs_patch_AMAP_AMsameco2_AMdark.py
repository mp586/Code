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
control_dir = 'withtv/aquaplanet_frierson_insolation_0qflux_mld20_commitd15c267'
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    control_dir= control_model + '/ISCA_HPC/' + control_dir
else: 
    control_dir= control_model + '/' + control_dir

#print control_dir
ctl_runmin=25
ctl_runmax=97

model = 'isca'
if (model == 'Isca') or (model == 'isca'): 
    model_data = 'Isca_DATA'
    output_dir1 = 'Isca'
elif (model == 'gfdl') or (model == 'GFDL'):
    model_data = 'GFDL_DATA'
    output_dir1 = ''

HPC = 'yes'
testdir_in1= 'withtv/square_South_America_frierson_insolation_lepref1_0qflux_samealbedo_to_01land_samehcp_landocean_commitd15c267' # this needs to be the dark patch simulation
runmin=24
runmax=96
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    exp1_name = 'ISCA_HPC_'+testdir_in1
    testdir = model_data + '/ISCA_HPC/' + testdir_in1
    testdir_in1 = '/ISCA_HPC/' + testdir_in1
else: 
    exp1_name = testdir_in1
    testdir = model_data + '/' + testdir_in1



area_array, dx, dy = ca.cell_area_all(t_res=42,base_dir='/scratch/mp586/Isca/') # added _all because then dx and dy are also returned 
area_array = xr.DataArray(area_array) # returned in units of m bzw m^2, because radius in cell_area.py is given in metres

area_array_3D = np.expand_dims(area_array, axis=0)
area_array_3D = np.repeat(area_array_3D, 40, axis = 0) # to make area_array 3D (pressure, lat, lon)


[sphum_u1_ctl,sphum_u1_avg_ctl,sphum_u1_seasonal_avg_ctl,sphum_u1_month_avg_ctl,sphum_u1_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'sphum_u','kgm/kgs')
[sphum_u1,sphum_u1_avg,sphum_u1_seasonal_avg,sphum_u1_month_avg,sphum_u1_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'sphum_u','kgm/kgs')

[sphum_v1_ctl,sphum_v1_avg_ctl,sphum_v1_seasonal_avg_ctl,sphum_v1_month_avg_ctl,sphum_v1_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'sphum_v','kgm/kgs')
[sphum_v1,sphum_v1_avg,sphum_v1_seasonal_avg,sphum_v1_month_avg,sphum_v1_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'sphum_v','kgm/kgs')

[ucomp_temp1,ucomp_temp1_avg,ucomp_temp1_seasonal_avg,ucomp_temp1_month_avg,ucomp_temp1_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'ucomp_temp','Km/s')
[ucomp_temp1_ctl,ucomp_temp1_avg_ctl,ucomp_temp1_seasonal_avg_ctl,ucomp_temp1_month_avg_ctl,ucomp_temp1_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'ucomp_temp','Km/s')

[vcomp_temp1,vcomp_temp1_avg,vcomp_temp1_seasonal_avg,vcomp_temp1_month_avg,vcomp_temp1_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'vcomp_temp','Km/s')
[vcomp_temp1_ctl,vcomp_temp1_avg_ctl,vcomp_temp1_seasonal_avg_ctl,vcomp_temp1_month_avg_ctl,vcomp_temp1_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'vcomp_temp','Km/s')

[height1,height1_avg,height1_seasonal_avg,height1_month_avg,height1_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'height','m')
[height1_ctl,height1_avg_ctl,height1_seasonal_avg_ctl,height1_month_avg_ctl,height1_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'height','m')

[ucomp1,ucomp1_avg,ucomp1_seasonal_avg,ucomp1_month_avg,ucomp1_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'ucomp','m/s')
[ucomp1_ctl,ucomp1_avg_ctl,ucomp1_seasonal_avg_ctl,ucomp1_month_avg_ctl,ucomp1_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'ucomp','m/s')

[vcomp1,vcomp1_avg,vcomp1_seasonal_avg,vcomp1_month_avg,vcomp1_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'vcomp','m/s')
[vcomp1_ctl,vcomp1_avg_ctl,vcomp1_seasonal_avg_ctl,vcomp1_month_avg_ctl,vcomp1_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'vcomp','m/s')




[tsurf1,tsurf1_avg,tsurf1_seasonal_avg,tsurf1_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'t_surf','K')
[lhe_flux1,lhe_flux1_avg,lhe_flux1_seasonal_avg,lhe_flux1_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lhe','W/m2',factor = 1.) # latent heat flux at surface (UP)
[net_lhe1,net_lhe1_avg,net_lhe1_seasonal_avg,net_lhe1_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[net_sw1,net_sw1_avg,net_sw1_seasonal_avg,net_sw1_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_sw','W/m^2',factor = 1.) # Pretty sure this is actually pos. downward! 
[lw_down1,lw_down1_avg,lw_down1_seasonal_avg,lw_down1_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lw','W/m^2',factor = 1.) # 
[net_t1,net_t1_avg,net_t1_seasonal_avg,net_t1_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_t','W/m^2',factor = 1.) # 
[toa_sw1,toa_sw1_avg,toa_sw1_seasonal_avg,toa_sw1_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'toa_sw','W/m^2',factor = 1.) # positive DOWN
[olr1,olr1_avg,olr1_seasonal_avg,olr1_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'olr','W/m^2',factor = 1.) # positive DOWN

[tsurf1_ctl,tsurf1_avg_ctl,tsurf1_seasonal_avg_ctl,tsurf1_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'t_surf','K')
[lhe_flux1_ctl,lhe_flux1_avg_ctl,lhe_flux1_seasonal_avg_ctl,lhe_flux1_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_lhe','W/m2',factor = 1.) # latent heat flux at surface (UP)
[net_lhe1_ctl,net_lhe1_avg_ctl,net_lhe1_seasonal_avg_ctl,net_lhe1_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[toa_sw1_ctl,toa_sw1_avg_ctl,toa_sw1_seasonal_avg_ctl,toa_sw1_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'toa_sw','W/m^2',factor = 1.) # positive DOWN
[net_sw1_ctl,net_sw1_avg_ctl,net_sw1_seasonal_avg_ctl,net_sw1_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_sw','W/m^2',factor = 1.) # 
[lw_down1_ctl,lw_down1_avg_ctl,lw_down1_seasonal_avg_ctl,lw_down1_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_lw','W/m^2',factor = 1.) # 
[net_t1_ctl,net_t1_avg_ctl,net_t_seasonal1_avg_ctl,net_t1_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_t','W/m^2',factor = 1.) # 


################ read in data from exp 2 ###############################

ctl_model = 'isca'
if (ctl_model == 'Isca') or (ctl_model == 'isca'): 
    control_model = 'Isca_DATA'
elif (ctl_model == 'gfdl') or (ctl_model == 'GFDL'):
    control_model = 'GFDL_DATA'

HPC = 'yes'
control_dir = 'square_South_America_frierson_insolation_newbucket_0qflux_samealbedo_samehcp_landocean_commitd15c267'
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    control_dir= control_model + '/ISCA_HPC/' + control_dir
else: 
    control_dir= control_model + '/' + control_dir

#print control_dir
ctl_runmin=121
ctl_runmax=241

model = 'isca'
if (model == 'Isca') or (model == 'isca'): 
    model_data = 'Isca_DATA'
    output_dir1 = 'Isca'
elif (model == 'gfdl') or (model == 'GFDL'):
    model_data = 'GFDL_DATA'
    output_dir1 = ''

HPC = 'yes'
testdir_in1= 'square_South_America_frierson_insolation_newbucket_0qflux_samealbedo_samehcp_landocean_plus_2xCO2_spinup_361_commitd15c267'
runmin=120
runmax=240
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    exp1_name = 'ISCA_HPC/'+testdir_in1
    testdir = model_data + '/ISCA_HPC/' + testdir_in1
    testdir_in1 = '/ISCA_HPC/' + testdir_in1
else: 
    exp1_name = testdir_in1
    testdir = model_data + '/' + testdir_in1

[sphum2_ctl,sphum2_avg_ctl,sphum2_seasonal_avg_ctl,sphum2_month_avg_ctl,sphum2_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'sphum','kg/kg')
[sphum2,sphum2_avg,sphum2_seasonal_avg,sphum2_month_avg,sphum2_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'sphum','kg/kg')

[temp2,temp2_avg,temp2_seasonal_avg,temp2_month_avg,temp2_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'temp','K')
[temp2_ctl,temp2_avg_ctl,temp2_seasonal_avg_ctl,temp2_month_avg_ctl,temp2_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'temp','K')

[height2,height2_avg,height2_seasonal_avg,height2_month_avg,height2_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'height','m')
[height2_ctl,height2_avg_ctl,height2_seasonal_avg_ctl,height2_month_avg_ctl,height2_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'height','m')

################ read in data from exp 3 ###############################

# ctl_model = 'isca'
# if (ctl_model == 'Isca') or (ctl_model == 'isca'): 
#     control_model = 'Isca_DATA'
# elif (ctl_model == 'gfdl') or (ctl_model == 'GFDL'):
#     control_model = 'GFDL_DATA'

# HPC = 'yes'
# control_dir = ''
# if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
#     control_dir= control_model + '/ISCA_HPC/' + control_dir
# else: 
#     control_dir= control_model + '/' + control_dir

# #print control_dir
# ctl_runmin=121  # Should be a January month for seasonal variables to be correct
# ctl_runmax=241

# model = 'isca'
# if (model == 'Isca') or (model == 'isca'): 
#     model_data = 'Isca_DATA'
#     output_dir1 = 'Isca'
# elif (model == 'gfdl') or (model == 'GFDL'):
#     model_data = 'GFDL_DATA'
#     output_dir1 = ''

# HPC = 'yes'
# testdir_in1= ''
# runmin=120
# runmax=240
# if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
#     exp1_name = 'ISCA_HPC/'+testdir_in1
#     testdir = model_data + '/ISCA_HPC/' + testdir_in1
#     testdir_in1 = '/ISCA_HPC/' + testdir_in1
# else: 
#     exp1_name = testdir_in1
#     testdir = model_data + '/' + testdir_in1



[sphum3,sphum3_avg,sphum3_seasonal_avg,sphum3_month_avg,sphum3_annual_avg,time]=[sphum2_ctl,sphum2_avg_ctl,sphum2_seasonal_avg_ctl,sphum2_month_avg_ctl,sphum2_annual_avg_ctl,time]
[temp3,temp3_avg,temp3_seasonal_avg,temp3_month_avg,temp3_annual_avg,time]=[temp2_ctl,temp2_avg_ctl,temp2_seasonal_avg_ctl,temp2_month_avg_ctl,temp2_annual_avg_ctl,time]
[height3,height3_avg,height3_seasonal_avg,height3_month_avg,height3_annual_avg,time]=[height2_ctl,height2_avg_ctl,height2_seasonal_avg_ctl,height2_month_avg_ctl,height2_annual_avg_ctl,time]

[sphum3_ctl,sphum3_avg_ctl,sphum3_seasonal_avg_ctl,sphum3_month_avg_ctl,sphum3_annual_avg_ctl,time]=[sphum1_ctl,sphum1_avg_ctl,sphum1_seasonal_avg_ctl,sphum1_month_avg_ctl,sphum1_annual_avg_ctl,time]
[temp3_ctl,temp3_avg_ctl,temp3_seasonal_avg_ctl,temp3_month_avg_ctl,temp3_annual_avg_ctl,time]=[temp1_ctl,temp1_avg_ctl,temp1_seasonal_avg_ctl,temp1_month_avg_ctl,temp1_annual_avg_ctl,time]
[height3_ctl,height3_avg_ctl,height3_seasonal_avg_ctl,height3_month_avg_ctl,height3_annual_avg_ctl,time]=[height1_ctl,height1_avg_ctl,height1_seasonal_avg_ctl,height1_month_avg_ctl,height1_annual_avg_ctl,time]

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



################ read in data from exp 4 ###############################

# ctl_model = 'isca'
# if (ctl_model == 'Isca') or (ctl_model == 'isca'): 
#     control_model = 'Isca_DATA'
# elif (ctl_model == 'gfdl') or (ctl_model == 'GFDL'):
#     control_model = 'GFDL_DATA'

# HPC = 'yes'
# control_dir = 'square_South_America_frierson_insolation_newbucket_0qflux_samehcp_landocean_commitd15c267'
# if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
#     control_dir= control_model + '/ISCA_HPC/' + control_dir
# else: 
#     control_dir= control_model + '/' + control_dir

# #print control_dir
# ctl_runmin=121
# ctl_runmax=241

model = 'isca'
if (model == 'Isca') or (model == 'isca'): 
    model_data = 'Isca_DATA'
    output_dir1 = 'Isca'
elif (model == 'gfdl') or (model == 'GFDL'):
    model_data = 'GFDL_DATA'
    output_dir1 = ''

HPC = 'yes'
testdir_in1= 'square_South_America_frierson_insolation_newbucket_0qflux_samealbedo_to_01land_samehcp_landocean_commitd15c267'
dire = testdir_in1
runmin=120
runmax=240
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    exp1_name = 'ISCA_HPC/'+testdir_in1
    testdir = model_data + '/ISCA_HPC/' + testdir_in1
    testdir_in1 = '/ISCA_HPC/' + testdir_in1
else: 
    exp1_name = testdir_in1
    testdir = model_data + '/' + testdir_in1


[sphum4,sphum4_avg,sphum4_seasonal_avg,sphum4_month_avg,sphum4_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'sphum','kg/kg')
[temp4,temp4_avg,temp4_seasonal_avg,temp4_month_avg,temp4_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'temp','K')
[height4,height4_avg,height4_seasonal_avg,height4_month_avg,height4_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'height','K')

[sphum4_ctl,sphum4_avg_ctl,sphum4_seasonal_avg_ctl,sphum4_month_avg_ctl,sphum4_annual_avg_ctl,time] = [sphum2_ctl,sphum2_avg_ctl,sphum2_seasonal_avg_ctl,sphum2_month_avg_ctl,sphum2_annual_avg_ctl,time]
[temp4_ctl,temp4_avg_ctl,temp4_seasonal_avg_ctl,temp4_month_avg_ctl,temp4_annual_avg_ctl,time] = [temp2_ctl,temp2_avg_ctl,temp2_seasonal_avg_ctl,temp2_month_avg_ctl,temp2_annual_avg_ctl,time]
[height4_ctl,height4_avg_ctl,height4_seasonal_avg_ctl,height4_month_avg_ctl,height4_annual_avg_ctl,time] = [height2_ctl,height2_avg_ctl,height2_seasonal_avg_ctl,height2_month_avg_ctl,height2_annual_avg_ctl,time]


outdir = 'Isca' + '/' + exp1_name

data = xr.open_dataset('/scratch/mp586/'+testdir+'/run0001/atmos_monthly.nc')
dp = xr.DataArray(data.phalf.diff('phalf').values*100, coords=[('pfull', data.pfull)])

pfull = ucomp_temp1.pres_lev * 100 # convert from hPa to Pa
pfull = np.expand_dims(pfull, axis = 1)
pfull = np.expand_dims(pfull, axis = 2)
pfull = np.repeat(pfull, 64, axis=1)
pfull = np.repeat(pfull, 128, axis = 2)
pres_lev = ucomp_temp1.pres_lev
pfull = xr.DataArray(pfull, coords = [pres_lev, ucomp_temp1.lat, ucomp_temp1.lon], dims = ['pres_lev','lat','lon'])

dp = np.expand_dims(dp, axis = 1)
dp = np.expand_dims(dp, axis = 2)
dp = np.repeat(dp, 64, axis=1)
dp = np.repeat(dp, 128, axis = 2)
dp = xr.DataArray(dp, coords = [pres_lev, ucomp_temp1.lat, ucomp_temp1.lon], dims = ['pres_lev','lat','lon'])





#conversion following https://www.ncl.ucar.edu/Document/Functions/Contributed/omega_to_w.shtml


cp_dry = 287.04/(2./7.) # units = J/kg/K, https://github.com/ExeClim/Isca/blob/77a3d49c5e3131dc6312d32b3698feac2cc8d156/postprocessing/plevel_interpolation/src/shared/constants/constants.F90
g = 9.81
L = 2.500e6 #J/kg, https://github.com/ExeClim/Isca/blob/77a3d49c5e3131dc6312d32b3698feac2cc8d156/postprocessing/plevel_interpolation/src/shared/constants/constants.F90


sensible_u1_avg = np.sum(cp_dry*ucomp_temp1_avg*dp/g, axis = 0)

sensible_v1_avg = np.sum(cp_dry*vcomp_temp1_avg*dp/g, axis = 0)

latent_u1_avg = np.sum(L*sphum_u1_avg*dp/g, axis = 0)

latent_v1_avg = np.sum(L*sphum_v1_avg*dp/g, axis = 0)

potential_u1_avg = np.sum(height1_avg*ucomp1_avg*dp, axis=0)

potential_v1_avg = np.sum(height1_avg*vcomp1_avg*dp, axis=0)






land = 'square_South_America'
landfile=Dataset(os.path.join(GFDL_BASE,'input/'+land+'/land.nc'),mode='r')

landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]
# for specified lats
landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it



sigma = 5.67*10**(-8)
net_lw1_avg = sigma*(tsurf1_avg**4) - lw_down1_avg # positive upward

#positive upward: energy input to the atmosphere is pos.
SEB1_areaavg = - area_weighted_avg(net_sw1_avg,area_array,landmaskxr,option='all_sfcs') + area_weighted_avg(net_lw1_avg,area_array,landmaskxr,option='all_sfcs') + area_weighted_avg(lhe_flux1_avg,area_array,landmaskxr,option='all_sfcs') + area_weighted_avg(net_t1_avg,area_array,landmaskxr,option='all_sfcs') # optional, if there is a qflux + area_weighted_avg(flux_oceanq_avg,area_array,landmaskxr,option='ocean')
print('SEB1 (perturbed) in W/m2 = '+str(SEB1_areaavg))

SEB1_avg = - net_sw1_avg + net_lw1_avg + lhe_flux1_avg + net_t1_avg ##perturbed run 
TOA1_avg = toa_sw1_avg - olr1_avg # toa sw pos down!, olr pos up (? yep)  # TOA pos down

N1_avg = SEB1_avg + TOA1_avg #total energy input to the atmosphere 

globavg_N1_avg = area_weighted_avg(N1_avg,area_array,landmaskxr,option='all_sfcs')

net_lw1_avg_ctl = sigma*(tsurf1_avg_ctl**4) - lw_down1_avg_ctl # positive upward

#positive upward: energy input to the atmosphere is pos.
SEB1_areaavg_ctl = area_weighted_avg( - net_sw1_avg_ctl,area_array,landmaskxr,option='all_sfcs') + area_weighted_avg(net_lw1_avg_ctl,area_array,landmaskxr,option='all_sfcs') + area_weighted_avg(lhe_flux1_avg_ctl,area_array,landmaskxr,option='all_sfcs') + area_weighted_avg(net_t1_avg_ctl,area_array,landmaskxr,option='all_sfcs') # optional, if there is a qflux + area_weighted_avg(flux_oceanq_avg,area_array,landmaskxr,option='ocean')
print('SEB (control) in W/m2 = '+str(SEB_areaavg_ctl))

SEB_avg_ctl = - net_sw_avg_ctl + net_lw_avg_ctl + lhe_flux_avg_ctl + net_t_avg_ctl ##perturbed run 

TOA_avg_ctl = toa_sw_avg_ctl - olr_avg_ctl # toa sw pos down!, olr pos up (?) 
TOA_avg = toa_sw_avg - olr_avg # toa sw pos down!, olr pos up (?) 

deltaTOA_areawav = area_weighted_avg((TOA_avg - TOA_avg_ctl),area_array,landmaskxr,option='all_sfcs')



N_avg_ctl = SEB_avg_ctl + TOA_avg_ctl #total energy input to the atmosphere 

DeltaN = N_avg - N_avg_ctl

DeltaN_zonalmean = DeltaN.mean('lon')
DeltaN_zonalmean = np.expand_dims(DeltaN_zonalmean,axis = 1)
DeltaN_zonalmean = np.repeat(DeltaN_zonalmean, 128, axis = 1)
DeltaN_zonalmean = xr.DataArray(DeltaN_zonalmean, coords = [N_avg.lat, N_avg.lon], dims = ['lat','lon'])




