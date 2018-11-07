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




ctl_model = input('Enter model 1 name as string ')
if (ctl_model == 'Isca') or (ctl_model == 'isca'): 
    control_model = 'Isca_DATA'
elif (ctl_model == 'gfdl') or (ctl_model == 'GFDL'):
    control_model = 'GFDL_DATA'

HPC = input('Data in ISCA_HPC ? Yes or No? ')
control_dir = input('Enter control directory from exp 1 as string ')
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    control_dir= control_model + '/ISCA_HPC/' + control_dir
else: 
    control_dir= control_model + '/' + control_dir

#print control_dir
ctl_runmin=input('Enter runmin number for exp 1 ')  # Should be a January month for seasonal variables to be correct
ctl_runmax=input('Enter runmax number for comparison 1 ')
ctl_timeseries_max = input('Enter end of ctl timeseries month for exp 1 ')

model = input('Enter model exp 1 ')
if (model == 'Isca') or (model == 'isca'): 
    model_data = 'Isca_DATA'
    output_dir1 = 'Isca'
elif (model == 'gfdl') or (model == 'GFDL'):
    model_data = 'GFDL_DATA'
    output_dir1 = ''

HPC = input('Data in ISCA_HPC ? Yes or No? ')
testdir_in1= input('Enter perturbed data directory name as string for exp 1 ')
runmin=input('Enter runmin number for exp 1 ')  # Should be a January month for seasonal variables to be correct
runmax=input('Enter runmax number for exp 1 ')
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    exp1_name = 'ISCA_HPC_'+testdir_in1
    testdir = model_data + '/ISCA_HPC/' + testdir_in1
    testdir_in1 = '/ISCA_HPC/' + testdir_in1
else: 
    exp1_name = testdir_in1
    testdir = model_data + '/' + testdir_in1

land = input('Which landmask? ')
landfile=Dataset(os.path.join(GFDL_BASE,'input/'+land+'/land.nc'),mode='r')

landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]
# for specified lats
landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

level = input('Which Level? ')


area_array, dx, dy = ca.cell_area_all(t_res=42,base_dir='/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/') # added _all because then dx and dy are also returned 
area_array = xr.DataArray(area_array) # returned in units of m bzw m^2, because radius in cell_area.py is given in metres

area_array_3D = np.expand_dims(area_array, axis=0)
area_array_3D = np.repeat(area_array_3D, 40, axis = 0) # to make area_array 3D (pressure, lat, lon)

total_sfc_area = np.sum(area_array)
#print ('total sfc area (*10^14) = '+str(np.sum(area_array/(10**14)))) # -- test: correct, equals sfc area of earth (5.1*10**14 m^2)
land_sfc_area = np.sum(area_array.where(landmask==1.))
#print ('land sfc area (*10^14) = '+str(land_sfc_area/(10**14)))
ocean_sfc_area = np.sum(area_array.where(landmask!=1.))
#print ('ocean sfc area (*10^14) = '+str(ocean_sfc_area/(10**14)))

# Read in variables 

[tsurf,tsurf_avg,tsurf_seasonal_avg,tsurf_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'t_surf','K')
[net_lhe,net_lhe_avg,net_lhe_seasonal_avg,net_lhe_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[precipitation,precipitation_avg,precipitation_seasonal_avg,precipitation_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'precipitation','mm/d', factor=86400)
[rh,rh_avg,rh_seasonal_avg,rh_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'rh','%',level=level)

[tsurf_ctl,tsurf_avg_ctl,tsurf_seasonal_avg_ctl,tsurf_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'t_surf','K')
[net_lhe_ctl,net_lhe_avg_ctl,net_lhe_seasonal_avg_ctl,net_lhe_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[precipitation_ctl,precipitation_avg_ctl,precipitation_seasonal_avg_ctl,precipitation_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'precipitation','mm/d', factor=86400)
[rh_ctl,rh_avg_ctl,rh_seasonal_avg_ctl,rh_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'rh','%',level=level)

[sphum,sphum_avg,sphum_seasonal_avg,sphum_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'sphum','kg/kg',level=level)
[sphum_ctl,sphum_avg_ctl,sphum_seasonal_avg_ctl,sphum_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'sphum','kg/kg',level=level)

# [div,div_avg,div_seasonal_avg,div_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'div','1/sec',level=level)
# [div_ctl,div_avg_ctl,div_seasonal_avg_ctl,div_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'div','1/sec',level=level)

PE_1 = precipitation - net_lhe
PE_ctl_1 = precipitation_ctl - net_lhe_ctl
DPE1 = PE_1 - PE_ctl_1 
[DPE1,DPE1_avg,DPE1_seasonal_avg,DPE1_month_avg,time] = make_var_seasonal(DPE1)


DP1 = precipitation - precipitation_ctl 
[DP1,DP1_avg,DP1_seasonal_avg,DP1_month_avg,time] = make_var_seasonal(DP1)

DT1 = tsurf - tsurf_ctl
[DT1,DT1_avg,DT1_seasonal_avg,DT1_month_avg,time] = make_var_seasonal(DT1)

DE1 = net_lhe - net_lhe_ctl
[DE1,DE1_avg,DE1_seasonal_avg,DE1_month_avg,time] = make_var_seasonal(DE1)

DRH1 = rh - rh_ctl
[DRH1,DRH1_avg,DRH1_seasonal_avg,DRH1_month_avg,time] = make_var_seasonal(DRH1)

DSH1 = sphum - sphum_ctl
[DSH1,DSH1_avg,DSH1_seasonal_avg,DSH1_month_avg,time] = make_var_seasonal(DSH1)



## ######### Enter stuff for experiment 2 ################

ctl_model = input('Enter model 2 name as string ')
if (ctl_model == 'Isca') or (ctl_model == 'isca'): 
    control_model = 'Isca_DATA'
elif (ctl_model == 'gfdl') or (ctl_model == 'GFDL'):
    control_model = 'GFDL_DATA'

HPC = input('Data in ISCA_HPC ? Yes or No? ')
control_dir = input('Enter control directory from exp 2 as string ')
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    control_dir= control_model + '/ISCA_HPC/' + control_dir
else: 
    control_dir= control_model + '/' + control_dir

#print control_dir
ctl_runmin=input('Enter runmin number for exp 2 ')  # Should be a January month for seasonal variables to be correct
ctl_runmax=input('Enter runmax number for comparison 2 ')
ctl_timeseries_max = input('Enter end of ctl timeseries month for exp 2 ')

model = input('Enter model exp 2 ')
if (model == 'Isca') or (model == 'isca'): 
    model_data = 'Isca_DATA'
    output_dir = 'Isca'
elif (model == 'gfdl') or (model == 'GFDL'):
    model_data = 'GFDL_DATA'
    output_dir = ''


HPC = input('Data in ISCA_HPC ? Yes or No? ')
testdir= input('Enter perturbed data directory name as string for exp 2 ')
runmin=input('Enter runmin number for exp 2 ')  # Should be a January month for seasonal variables to be correct
runmax=input('Enter runmax number for exp 2 ')
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    exp2_name = 'ISCA_HPC_'+testdir
    testdir = model_data + '/ISCA_HPC/' + testdir
else: 
    exp2_name = testdir
    testdir = model_data + '/' + testdir


if not os.path.exists('/scratch/mp586/Code/Graphics/'+output_dir1+'/'+testdir_in1+'/comparison_to_'+exp2_name):
    os.mkdir('/scratch/mp586/Code/Graphics/'+output_dir1+'/'+testdir_in1+'/comparison_to_'+exp2_name)
outdir_saving = output_dir1+'/'+testdir_in1+'/comparison_to_'+exp2_name
print('Will save figures to '+outdir_saving)


# Read in variables 

[tsurf,tsurf_avg,tsurf_seasonal_avg,tsurf_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'t_surf','K')
[net_lhe,net_lhe_avg,net_lhe_seasonal_avg,net_lhe_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[precipitation,precipitation_avg,precipitation_seasonal_avg,precipitation_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'precipitation','mm/d', factor=86400)
[rh,rh_avg,rh_seasonal_avg,rh_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'rh','%',level=level)

[tsurf_ctl,tsurf_avg_ctl,tsurf_seasonal_avg_ctl,tsurf_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'t_surf','K')
[net_lhe_ctl,net_lhe_avg_ctl,net_lhe_seasonal_avg_ctl,net_lhe_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[precipitation_ctl,precipitation_avg_ctl,precipitation_seasonal_avg_ctl,precipitation_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'precipitation','mm/d', factor=86400)
[rh_ctl,rh_avg_ctl,rh_seasonal_avg_ctl,rh_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'rh','%',level=level)

[sphum,sphum_avg,sphum_seasonal_avg,sphum_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'sphum','kg/kg',level=level)
[sphum_ctl,sphum_avg_ctl,sphum_seasonal_avg_ctl,sphum_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'sphum','kg/kg',level=level)

# [div,div_avg,div_seasonal_avg,div_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'div','1/sec',level=level)
# [div_ctl,div_avg_ctl,div_seasonal_avg_ctl,div_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'div','1/sec',level=level)

PE_2 = precipitation - net_lhe
PE_ctl_2 = precipitation_ctl - net_lhe_ctl
DPE2 = PE_2 - PE_ctl_2 
[DPE2,DPE2_avg,DPE2_seasonal_avg,DPE2_month_avg,time] = make_var_seasonal(DPE2)


DP2 = precipitation - precipitation_ctl 
[DP2,DP2_avg,DP2_seasonal_avg,DP2_month_avg,time] = make_var_seasonal(DP2)

DT2 = tsurf - tsurf_ctl
[DT2,DT2_avg,DT2_seasonal_avg,DT2_month_avg,time] = make_var_seasonal(DT2)


DE2 = net_lhe - net_lhe_ctl
[DE2,DE2_avg,DE2_seasonal_avg,DE2_month_avg,time] = make_var_seasonal(DE2)

DRH2 = rh - rh_ctl
[DRH2,DRH2_avg,DRH2_seasonal_avg,DRH2_month_avg,time] = make_var_seasonal(DRH2)

DSH2 = sphum - sphum_ctl
[DSH2,DSH2_avg,DSH2_seasonal_avg,DSH2_month_avg,time] = make_var_seasonal(DSH2)




any_configuration_plot(outdir_saving,runmin,runmax,-90.,90.,(DT2_avg - DT1_avg),area_array,'K','$\Delta$_$T_S$_difference','tempdiff',landmaskxr, minval = -5., maxval = 5., save_title = 'Delta_TS_')
any_configuration_plot(outdir_saving,runmin,runmax,-90.,90.,(DP2_avg - DP1_avg),area_array,'mm/d','$\Delta$_P_difference','rainnorm',landmaskxr, minval = -2., maxval = 2., save_title = 'Delta_P_')
any_configuration_plot(outdir_saving,runmin,runmax,-90.,90.,(DE2_avg - DE1_avg),area_array,'mm/d','$\Delta$_E_difference','rainnorm',landmaskxr, minval = -1., maxval = 1., save_title = 'Delta_E_')
any_configuration_plot(outdir_saving,runmin,runmax,-90.,90.,(DRH2_avg - DRH1_avg),area_array,'%','$\Delta$_rh_lev'+str(level)+'difference','rainnorm',landmaskxr, minval = -4., maxval = 4., save_title = 'Delta_rh_lev'+str(level))
any_configuration_plot(outdir_saving,runmin,runmax,-90.,90.,(DSH2_avg - DSH1_avg),area_array,'kg/kg','$\Delta$_sphum_lev'+str(level)+'_difference','rainnorm',landmaskxr, minval = -0.001, maxval = 0.001, save_title = 'Delta_sphum_lev'+str(level))
