from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import xarray as xr
import pandas as pd
import os
from matplotlib.patches import Rectangle

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import metpy as metpy

import sys
sys.path.insert(0, '/scratch/mp586/Code/PYCODES')
from plotting_routines_kav7_tabreplaced_ongv3_py3 import *
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


data_pert = xr.open_mfdataset('/scratch/mp586/'+testdir+'/run003*/atmos_monthly_interp.nc')

print(data_pert)
data_pert = data_pert.metpy.parse_cf()





# ################ read in data from exp 2 ###############################

# ctl_model = 'isca'
# if (ctl_model == 'Isca') or (ctl_model == 'isca'): 
#     control_model = 'Isca_DATA'
# elif (ctl_model == 'gfdl') or (ctl_model == 'GFDL'):
#     control_model = 'GFDL_DATA'

# HPC = 'yes'
# control_dir = 'square_South_America_frierson_insolation_newbucket_0qflux_samealbedo_samehcp_landocean_commitd15c267'
# if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
#     control_dir= control_model + '/ISCA_HPC/' + control_dir
# else: 
#     control_dir= control_model + '/' + control_dir

# #print control_dir
# ctl_runmin=121
# ctl_runmax=241

# model = 'isca'
# if (model == 'Isca') or (model == 'isca'): 
#     model_data = 'Isca_DATA'
#     output_dir1 = 'Isca'
# elif (model == 'gfdl') or (model == 'GFDL'):
#     model_data = 'GFDL_DATA'
#     output_dir1 = ''

# HPC = 'yes'
# testdir_in1= 'square_South_America_frierson_insolation_newbucket_0qflux_samealbedo_samehcp_landocean_plus_2xCO2_spinup_361_commitd15c267'
# runmin=120
# runmax=240
# if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
#     exp1_name = 'ISCA_HPC/'+testdir_in1
#     testdir = model_data + '/ISCA_HPC/' + testdir_in1
#     testdir_in1 = '/ISCA_HPC/' + testdir_in1
# else: 
#     exp1_name = testdir_in1
#     testdir = model_data + '/' + testdir_in1

# [sphum2_ctl,sphum2_avg_ctl,sphum2_seasonal_avg_ctl,sphum2_month_avg_ctl,sphum2_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'sphum','kg/kg')
# [sphum2,sphum2_avg,sphum2_seasonal_avg,sphum2_month_avg,sphum2_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'sphum','kg/kg')

# [temp2,temp2_avg,temp2_seasonal_avg,temp2_month_avg,temp2_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'temp','K')
# [temp2_ctl,temp2_avg_ctl,temp2_seasonal_avg_ctl,temp2_month_avg_ctl,temp2_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'temp','K')

# [height2,height2_avg,height2_seasonal_avg,height2_month_avg,height2_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'height','m')
# [height2_ctl,height2_avg_ctl,height2_seasonal_avg_ctl,height2_month_avg_ctl,height2_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'height','m')

# ################ read in data from exp 3 ###############################

# # ctl_model = 'isca'
# # if (ctl_model == 'Isca') or (ctl_model == 'isca'): 
# #     control_model = 'Isca_DATA'
# # elif (ctl_model == 'gfdl') or (ctl_model == 'GFDL'):
# #     control_model = 'GFDL_DATA'

# # HPC = 'yes'
# # control_dir = ''
# # if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
# #     control_dir= control_model + '/ISCA_HPC/' + control_dir
# # else: 
# #     control_dir= control_model + '/' + control_dir

# # #print control_dir
# # ctl_runmin=121  # Should be a January month for seasonal variables to be correct
# # ctl_runmax=241

# # model = 'isca'
# # if (model == 'Isca') or (model == 'isca'): 
# #     model_data = 'Isca_DATA'
# #     output_dir1 = 'Isca'
# # elif (model == 'gfdl') or (model == 'GFDL'):
# #     model_data = 'GFDL_DATA'
# #     output_dir1 = ''

# # HPC = 'yes'
# # testdir_in1= ''
# # runmin=120
# # runmax=240
# # if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
# #     exp1_name = 'ISCA_HPC/'+testdir_in1
# #     testdir = model_data + '/ISCA_HPC/' + testdir_in1
# #     testdir_in1 = '/ISCA_HPC/' + testdir_in1
# # else: 
# #     exp1_name = testdir_in1
# #     testdir = model_data + '/' + testdir_in1



# [sphum3,sphum3_avg,sphum3_seasonal_avg,sphum3_month_avg,sphum3_annual_avg,time]=[sphum2_ctl,sphum2_avg_ctl,sphum2_seasonal_avg_ctl,sphum2_month_avg_ctl,sphum2_annual_avg_ctl,time]
# [temp3,temp3_avg,temp3_seasonal_avg,temp3_month_avg,temp3_annual_avg,time]=[temp2_ctl,temp2_avg_ctl,temp2_seasonal_avg_ctl,temp2_month_avg_ctl,temp2_annual_avg_ctl,time]
# [height3,height3_avg,height3_seasonal_avg,height3_month_avg,height3_annual_avg,time]=[height2_ctl,height2_avg_ctl,height2_seasonal_avg_ctl,height2_month_avg_ctl,height2_annual_avg_ctl,time]

# [sphum3_ctl,sphum3_avg_ctl,sphum3_seasonal_avg_ctl,sphum3_month_avg_ctl,sphum3_annual_avg_ctl,time]=[sphum1_ctl,sphum1_avg_ctl,sphum1_seasonal_avg_ctl,sphum1_month_avg_ctl,sphum1_annual_avg_ctl,time]
# [temp3_ctl,temp3_avg_ctl,temp3_seasonal_avg_ctl,temp3_month_avg_ctl,temp3_annual_avg_ctl,time]=[temp1_ctl,temp1_avg_ctl,temp1_seasonal_avg_ctl,temp1_month_avg_ctl,temp1_annual_avg_ctl,time]
# [height3_ctl,height3_avg_ctl,height3_seasonal_avg_ctl,height3_month_avg_ctl,height3_annual_avg_ctl,time]=[height1_ctl,height1_avg_ctl,height1_seasonal_avg_ctl,height1_month_avg_ctl,height1_annual_avg_ctl,time]

# # [omega3,omega3_avg,omega3_seasonal_avg,omega3_month_avg,omega3_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'omega','Pa/s')
# # [omega3_ctl,omega3_avg_ctl,omega3_seasonal_avg_ctl,omega3_month_avg_ctl,omega3_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'omega','Pa/s')
# # [rh3_ctl,rh3_avg_ctl,rh3_seasonal_avg_ctl,rh3_month_avg_ctl,rh3_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'rh','%')
# # [rh3,rh3_avg,rh3_seasonal_avg,rh3_month_avg,rh3_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'rh','%')
# # [sphum3_ctl,sphum3_avg_ctl,sphum3_seasonal_avg_ctl,sphum3_month_avg_ctl,sphum3_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'sphum','kg/kg')
# # [sphum3,sphum3_avg,sphum3_seasonal_avg,sphum3_month_avg,sphum3_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'sphum','kg/kg')

# # [ucomp3,ucomp3_avg,ucomp3_seasonal_avg,ucomp3_month_avg,ucomp3_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'ucomp','m/s')
# # [ucomp3_ctl,ucomp3_avg_ctl,ucomp3_seasonal_avg_ctl,ucomp3_month_avg_ctl,ucomp3_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'ucomp','m/s')

# # [temp3,temp3_avg,temp3_seasonal_avg,temp3_month_avg,temp3_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'temp','K')
# # [temp3_ctl,temp3_avg_ctl,temp3_seasonal_avg_ctl,temp3_month_avg_ctl,temp3_annual_avg_ctl,time]=seasonal_4D_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'temp','K')



# ################ read in data from exp 4 ###############################

# # ctl_model = 'isca'
# # if (ctl_model == 'Isca') or (ctl_model == 'isca'): 
# #     control_model = 'Isca_DATA'
# # elif (ctl_model == 'gfdl') or (ctl_model == 'GFDL'):
# #     control_model = 'GFDL_DATA'

# # HPC = 'yes'
# # control_dir = 'square_South_America_frierson_insolation_newbucket_0qflux_samehcp_landocean_commitd15c267'
# # if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
# #     control_dir= control_model + '/ISCA_HPC/' + control_dir
# # else: 
# #     control_dir= control_model + '/' + control_dir

# # #print control_dir
# # ctl_runmin=121
# # ctl_runmax=241

# model = 'isca'
# if (model == 'Isca') or (model == 'isca'): 
#     model_data = 'Isca_DATA'
#     output_dir1 = 'Isca'
# elif (model == 'gfdl') or (model == 'GFDL'):
#     model_data = 'GFDL_DATA'
#     output_dir1 = ''

# HPC = 'yes'
# testdir_in1= 'square_South_America_frierson_insolation_newbucket_0qflux_samealbedo_to_01land_samehcp_landocean_commitd15c267'
# dire = testdir_in1
# runmin=120
# runmax=240
# if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
#     exp1_name = 'ISCA_HPC/'+testdir_in1
#     testdir = model_data + '/ISCA_HPC/' + testdir_in1
#     testdir_in1 = '/ISCA_HPC/' + testdir_in1
# else: 
#     exp1_name = testdir_in1
#     testdir = model_data + '/' + testdir_in1


# [sphum4,sphum4_avg,sphum4_seasonal_avg,sphum4_month_avg,sphum4_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'sphum','kg/kg')
# [temp4,temp4_avg,temp4_seasonal_avg,temp4_month_avg,temp4_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'temp','K')
# [height4,height4_avg,height4_seasonal_avg,height4_month_avg,height4_annual_avg,time]=seasonal_4D_variable_interp(testdir,model,runmin,runmax,'height','K')

# [sphum4_ctl,sphum4_avg_ctl,sphum4_seasonal_avg_ctl,sphum4_month_avg_ctl,sphum4_annual_avg_ctl,time] = [sphum2_ctl,sphum2_avg_ctl,sphum2_seasonal_avg_ctl,sphum2_month_avg_ctl,sphum2_annual_avg_ctl,time]
# [temp4_ctl,temp4_avg_ctl,temp4_seasonal_avg_ctl,temp4_month_avg_ctl,temp4_annual_avg_ctl,time] = [temp2_ctl,temp2_avg_ctl,temp2_seasonal_avg_ctl,temp2_month_avg_ctl,temp2_annual_avg_ctl,time]
# [height4_ctl,height4_avg_ctl,height4_seasonal_avg_ctl,height4_month_avg_ctl,height4_annual_avg_ctl,time] = [height2_ctl,height2_avg_ctl,height2_seasonal_avg_ctl,height2_month_avg_ctl,height2_annual_avg_ctl,time]

pfull = data_pert.pfull

data = xr.open_dataset('/scratch/mp586/'+testdir+'/run0001/atmos_monthly.nc')
data = data.metpy.parse_cf()
dp = data.phalf.diff('phalf')*100
dp = xr.DataArray(dp[::-1], coords = [pfull.values], dims = ['pfull'])




dxmp, dymp = metpy.calc.grid_deltas_from_dataarray(data_pert.ucomp)



#conversion following https://www.ncl.ucar.edu/Document/Functions/Contributed/omega_to_w.shtml


cp_dry = 287.04/(2./7.) # units = J/kg/K, https://github.com/ExeClim/Isca/blob/77a3d49c5e3131dc6312d32b3698feac2cc8d156/postprocessing/plevel_interpolation/src/shared/constants/constants.F90
g = 9.81
L = 2.500e6 #J/kg, https://github.com/ExeClim/Isca/blob/77a3d49c5e3131dc6312d32b3698feac2cc8d156/postprocessing/plevel_interpolation/src/shared/constants/constants.F90

vcomp_temp_avg = (data_pert.vcomp_temp).mean('time')
sensible_v = np.sum((cp_dry/g)*vcomp_temp_avg*dp, axis = 0)

ucomp_temp_avg = (data_pert.ucomp_temp).mean('time')
sensible_u = np.sum((cp_dry/g)*ucomp_temp_avg*dp, axis = 0)

div_sens = metpy.calc.divergence(sensible_u, sensible_v, dxmp[0,0,:,:], dymp[0,0,:,:])






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




