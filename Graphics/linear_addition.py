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


# Model which both runs will be referenced to 
ref_model = input('Enter reference model name as string ')
if (ref_model == 'Isca') or (ref_model == 'isca'): 
    reference_model = 'Isca_DATA'
elif (ref_model == 'gfdl') or (ref_model == 'GFDL'):
    reference_model = 'GFDL_DATA'


reference_dir= reference_model + '/' + input('Enter reference directory name as string ')
#print reference_dir
ref_runmin=input('Enter runmin number ')  # Should be a January month for seasonal variables to be correct
ref_runmax=input('Enter runmax number for comparison ')




fst_model = input('Enter first model as string ')
if (fst_model == 'Isca') or (fst_model == 'isca'): 
    first_model = 'Isca_DATA'
    output_dir = 'Isca'
elif (fst_model == 'gfdl') or (fst_model == 'GFDL'):
    first_model = 'GFDL_DATA'
    output_dir = ''

exp_name = input('Enter first directory name as string ')

first_dir= first_model + '/' + exp_name
fst_outdir = output_dir + '/' + exp_name

#print first_dir
fst_runmin=input('Enter runmin number ')  # Should be a January month for seasonal variables to be correct
fst_runmax=input('Enter runmax number for comparison ')

scnd_model = input('Enter second model as string ')
if (scnd_model == 'Isca') or (scnd_model == 'isca'): 
    second_model = 'Isca_DATA'
    output_dir = 'Isca'

elif (scnd_model == 'gfdl') or (scnd_model == 'GFDL'):
    second_model = 'GFDL_DATA'
    output_dir = ''


exp_name = input('Enter second directory name as string ')

second_dir= second_model + '/' + exp_name
scnd_outdir = output_dir + '/' + exp_name

scnd_runmin=input('Enter runmin number ')  # Should be a January month for seasonal variables to be correct
scnd_runmax=input('Enter runmax number for comparison ')


thrd_model = input('Enter third model as string ')
if (thrd_model == 'Isca') or (thrd_model == 'isca'): 
    third_model = 'Isca_DATA'
    output_dir = 'Isca'

elif (thrd_model == 'gfdl') or (thrd_model == 'GFDL'):
    third_model = 'GFDL_DATA'
    output_dir = ''


exp_name = input('Enter third directory name as string ')

third_dir= third_model + '/' + exp_name
thrd_outdir = output_dir + '/' + exp_name

thrd_runmin=input('Enter runmin number ')  # Should be a January month for seasonal variables to be correct
thrd_runmax=input('Enter runmax number for comparison ')



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

# [slp_fst,slp_avg_fst,slp_seasonal_avg_fst,slp_month_avg_fst,time]=seasonal_surface_variable(first_dir,fst_model,fst_runmin,fst_runmax,'slp','hPa',factor=10**(-2))
# [tsurf_fst,tsurf_avg_fst,tsurf_seasonal_avg_fst,tsurf_month_avg_fst,time]=seasonal_surface_variable(first_dir,fst_model,fst_runmin,fst_runmax,'t_surf','K')
# [lhe_flux_fst,lhe_flux_avg_fst,lhe_flux_seasonal_avg_fst,lhe_flux_month_avg_fst,time]=seasonal_surface_variable(first_dir,fst_model,fst_runmin,fst_runmax,'flux_lhe','W/m2',factor = 1.) # latent heat flux at surface (UP)
# [net_lhe_fst,net_lhe_avg_fst,net_lhe_seasonal_avg_fst,net_lhe_month_avg_fst,time]=seasonal_surface_variable(first_dir,fst_model,fst_runmin,fst_runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)

# [precipitation_fst,precipitation_avg_fst,precipitation_seasonal_avg_fst,precipitation_month_avg_fst,time]=seasonal_surface_variable(first_dir,fst_model,fst_runmin,fst_runmax,'precipitation','mm/d', factor=86400)
# #[bucket_depth_fst,bucket_depth_avg_fst,bucket_depth_seasonal_avg_fst,bucket_depth_month_avg_fst,time]=seasonal_surface_variable(first_dir,fst_model,fst_runmin,fst_runmax,'bucket_depth','m')
# [rh_fst,rh_avg_fst,rh_seasonal_avg_fst,rh_month_avg_fst,time]=seasonal_surface_variable(first_dir,fst_model,fst_runmin,fst_runmax,'rh','%',level=level)
# # [flux_oceanq_fst,flux_oceanq_avg_fst,flux_oceanq_seasonal_avg_fst,flux_oceanq_month_avg_fst,time]=seasonal_surface_variable(first_dir,fst_model,fst_runmin,fst_runmax,'flux_oceanq','W/m^2')
# [net_sw_fst,net_sw_avg_fst,net_sw_seasonal_avg_fst,net_sw_month_avg_fst,time]=seasonal_surface_variable(first_dir,fst_model,fst_runmin,fst_runmax,'flux_sw','W/m^2',factor = 1.) # 
# [net_lw_fst,net_lw_avg_fst,net_lw_seasonal_avg_fst,net_lw_month_avg_fst,time]=seasonal_surface_variable(first_dir,fst_model,fst_runmin,fst_runmax,'flux_lw','W/m^2',factor = 1.) # 
# [net_t_fst,net_t_avg_fst,net_t_seasonal_avg_fst,net_t_month_avg_fst,time]=seasonal_surface_variable(first_dir,fst_model,fst_runmin,fst_runmax,'flux_t','W/m^2',factor = 1.) # 
# #[CIWV_fst,CIWV_avg_fst,CIWV_seasonal_avg_fst,CIWV_month_avg_fst,time]=seasonal_surface_variable(first_dir,fst_model,fst_runmin,fst_runmax,'sphum','kg/kg',level='all')

# [sphum_fst,sphum_avg_fst,sphum_seasonal_avg_fst,sphum_month_avg_fst,time]=seasonal_surface_variable(first_dir,fst_model,fst_runmin,fst_runmax,'sphum','kg/kg',level=level)
# [div_fst,div_avg_fst,div_seasonal_avg_fst,div_month_avg_fst,time]=seasonal_surface_variable(first_dir,fst_model,fst_runmin,fst_runmax,'div','1/sec',level=level)

# PE_fst = precipitation_fst - net_lhe_fst
# [PE_fst,PE_avg_fst,PE_seasonal_avg_fst,PE_month_avg_fst,time] = make_var_seasonal(PE_fst)


# 870 hPa 
[gph_fst,gph_avg_fst,gph_seasonal_avg_fst,gph_month_avg_fst,time]=seasonal_surface_variable(first_dir,fst_model,fst_runmin,fst_runmax,'height','m', level = 37) 
[ucomp_fst,ucomp_avg_fst,ucomp_seasonal_avg_fst,ucomp_month_avg_fst,time]=seasonal_surface_variable(first_dir,fst_model,fst_runmin,fst_runmax,'ucomp','m/s', level = 37)
[vcomp_fst,vcomp_avg_fst,vcomp_seasonal_avg_fst,vcomp_month_avg_fst,time]=seasonal_surface_variable(first_dir,fst_model,fst_runmin,fst_runmax,'vcomp','m/s', level = 37)


[gph_scnd,gph_avg_scnd,gph_seasonal_avg_scnd,gph_month_avg_scnd,time]=seasonal_surface_variable(second_dir,scnd_model,scnd_runmin,scnd_runmax,'height','m', level = 37) 
[ucomp_scnd,ucomp_avg_scnd,ucomp_seasonal_avg_scnd,ucomp_month_avg_scnd,time]=seasonal_surface_variable(second_dir,scnd_model,scnd_runmin,scnd_runmax,'ucomp','m/s', level = 37)
[vcomp_scnd,vcomp_avg_scnd,vcomp_seasonal_avg_scnd,vcomp_month_avg_scnd,time]=seasonal_surface_variable(second_dir,scnd_model,scnd_runmin,scnd_runmax,'vcomp','m/s', level = 37)


[gph_ref,gph_avg_ref,gph_seasonal_avg_ref,gph_month_avg_ref,time]=seasonal_surface_variable(reference_dir,ref_model,ref_runmin,ref_runmax,'height','m', level = 37) 
[ucomp_ref,ucomp_avg_ref,ucomp_seasonal_avg_ref,ucomp_month_avg_ref,time]=seasonal_surface_variable(reference_dir,ref_model,ref_runmin,ref_runmax,'ucomp','m/s', level = 37)
[vcomp_ref,vcomp_avg_ref,vcomp_seasonal_avg_ref,vcomp_month_avg_ref,time]=seasonal_surface_variable(reference_dir,ref_model,ref_runmin,ref_runmax,'vcomp','m/s', level = 37)

[gph_thrd,gph_avg_thrd,gph_seasonal_avg_thrd,gph_month_avg_thrd,time]=seasonal_surface_variable(third_dir,thrd_model,thrd_runmin,thrd_runmax,'height','m', level = 37) 
[ucomp_thrd,ucomp_avg_thrd,ucomp_seasonal_avg_thrd,ucomp_month_avg_thrd,time]=seasonal_surface_variable(third_dir,thrd_model,thrd_runmin,thrd_runmax,'ucomp','m/s', level = 37)
[vcomp_thrd,vcomp_avg_thrd,vcomp_seasonal_avg_thrd,vcomp_month_avg_thrd,time]=seasonal_surface_variable(third_dir,thrd_model,thrd_runmin,thrd_runmax,'vcomp','m/s', level = 37)

winds_one_level(fst_outdir,fst_runmin,fst_runmax,'870hPa_gph_winds_square_Africa_minus_AP',ucomp_avg_fst - ucomp_avg_ref,vcomp_avg_fst - vcomp_avg_ref,gph_avg_fst - gph_avg_ref,
    'slp','m',15., 30.,landmaskxr,veclen=10,level=37,units_numerator = 'm', units_denom = 's')

winds_one_level(scnd_outdir,scnd_runmin,scnd_runmax,'870hPa_gph_winds_two_continents_minus_AP',(ucomp_avg_thrd - ucomp_avg_scnd - ucomp_avg_fst + ucomp_avg_ref), 
    (vcomp_avg_thrd - vcomp_avg_scnd - vcomp_avg_fst + vcomp_avg_ref),(gph_avg_thrd - gph_avg_scnd - gph_avg_fst + gph_avg_ref),'slp','m',-30., 30., landmaskxr,
    veclen=5,level=37,units_numerator = 'm', units_denom = 's')

