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

sys.path.insert(0, '/scratch/mp586/Code/Graphics/Graphics_for_chapter3/MSE')
import gradients as gr, model_constants as mc

ctl_model = input('Enter model name as string ')
if (ctl_model == 'Isca') or (ctl_model == 'isca'): 
    control_model = 'Isca_DATA'
elif (ctl_model == 'gfdl') or (ctl_model == 'GFDL'):
    control_model = 'GFDL_DATA'


control_dir= control_model + '/' + input('Enter control directory name as string ')
#print control_dir
ctl_runmin=input('Enter runmin number ')  # Should be a January month for seasonal variables to be correct
ctl_runmax=input('Enter runmax number for comparison ')
ctl_timeseries_max = input('Enter end of ctl timeseries month ')

model = input('Enter model ')
if (model == 'Isca') or (model == 'isca'): 
    model_data = 'Isca_DATA'
    output_dir = 'Isca'
elif (model == 'gfdl') or (model == 'GFDL'):
    model_data = 'GFDL_DATA'
    output_dir = ''

testdir= input('Enter data directory name as string ')
runmin=input('Enter runmin number ')  # Should be a January month for seasonal variables to be correct
runmax=input('Enter runmax number ')


outdir = output_dir + '/' + testdir
testdir = model_data + '/' + testdir

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


area_array_3Dtime = np.expand_dims(area_array, axis=0)
area_array_3Dtime = np.repeat(area_array_3Dtime, 12, axis = 0) # to make area_array 3D (pressure, lat, lon)

total_sfc_area = np.sum(area_array)
#print ('total sfc area (*10^14) = '+str(np.sum(area_array/(10**14)))) # -- test: correct, equals sfc area of earth (5.1*10**14 m^2)
land_sfc_area = np.sum(area_array.where(landmask==1.))
#print ('land sfc area (*10^14) = '+str(land_sfc_area/(10**14)))
ocean_sfc_area = np.sum(area_array.where(landmask!=1.))
#print ('ocean sfc area (*10^14) = '+str(ocean_sfc_area/(10**14)))

# for plotting a spin up run ('control') timeseries followed by the timeseries from the perturbed experiment
globavg_var_timeseries_total_and_land_perturbed(testdir,outdir,model,area_array,'t_surf',1,runmax,1.,landmaskxr,control_dir,ctl_model,1,ctl_timeseries_max)
#globavg_var_timeseries_total_and_land_perturbed(testdir,outdir,model,area_array,'bucket_depth',1,runmax,1.,landmask,control_dir,ctl_model,1,ctl_timeseries_max,select='land')

# mass stream function adapted from J Penn
# [msf,msf_avg,msf_seasonal_avg,msf_month_avg] = mass_streamfunction(testdir,model,runmin,runmax) # 
# plot_streamfunction_seasonal(msf_seasonal_avg, outdir, runmin, runmax, ctl_pert = 'pert')

[msf_ctl,msf_avg_ctl,msf_seasonal_avg_ctl,msf_month_avg_ctl] = mass_streamfunction_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax) # 
plot_streamfunction_seasonal(msf_seasonal_avg_ctl, outdir, ctl_runmin, ctl_runmax, minval=-30., maxval=30., steps=21, ctl_pert = 'ctl')
plot_streamfunction_seasonal(msf_seasonal_avg_ctl.sel(lat=slice(-50.,50.)), outdir, ctl_runmin, ctl_runmax, minval=-30., maxval=30., steps=21, ctl_pert = 'ctl_tropics', figx=15)

# [CIWV_ctl,CIWV_avg_ctl,CIWV_seasonal_avg_ctl,CIWV_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'sphum','kg/kg',level='all')
# any_configuration_plot(outdir,runmin,runmax,-90.,90.,(CIWV_avg_ctl),area_array,'kg/kg','CIWV_ctl','fromwhite',landmaskxr,nmb_contours=0, minval = 0, maxval = 0.011)


# Read in variables 

[slp,slp_avg,slp_seasonal_avg,slp_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'slp','hPa',factor=10**(-2))
[tsurf,tsurf_avg,tsurf_seasonal_avg,tsurf_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'t_surf','K')
[lhe_flux,lhe_flux_avg,lhe_flux_seasonal_avg,lhe_flux_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lhe','W/m2',factor = 1.) # latent heat flux at surface (UP)
[net_lhe,net_lhe_avg,net_lhe_seasonal_avg,net_lhe_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)

[precipitation,precipitation_avg,precipitation_seasonal_avg,precipitation_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'precipitation','mm/d', factor=86400)
#[bucket_depth,bucket_depth_avg,bucket_depth_seasonal_avg,bucket_depth_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'bucket_depth','m')
[rh,rh_avg,rh_seasonal_avg,rh_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'rh','%',level=level)
[net_sw,net_sw_avg,net_sw_seasonal_avg,net_sw_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_sw','W/m^2',factor = 1.) # Pretty sure this is actually pos. downward! 
[lw_down,lw_down_avg,lw_down_seasonal_avg,lw_down_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lw','W/m^2',factor = 1.) # 
#[CIWV,CIWV_avg,CIWV_seasonal_avg,CIWV_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'sphum','kg/kg',level='all')
[net_t,net_t_avg,net_t_seasonal_avg,net_t_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_t','W/m^2',factor = 1.) # 
[toa_sw,toa_sw_avg,toa_sw_seasonal_avg,toa_sw_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'toa_sw','W/m^2',factor = 1.) # positive DOWN



[gph,gph_avg,gph_seasonal_avg,gph_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'height','m', level = 27) # geopotential height at 500 hPa
[gph_ctl,gph_avg_ctl,gph_seasonal_avg_ctl,gph_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'height','m', level = 27) # geopotential height at 500 hPa


[slp_ctl,slp_avg_ctl,slp_seasonal_avg_ctl,slp_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'slp','hPa',factor=10**(-2))
[tsurf_ctl,tsurf_avg_ctl,tsurf_seasonal_avg_ctl,tsurf_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'t_surf','K')
[lhe_flux_ctl,lhe_flux_avg_ctl,lhe_flux_seasonal_avg_ctl,lhe_flux_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_lhe','W/m2',factor = 1.) # latent heat flux at surface (UP)
[net_lhe_ctl,net_lhe_avg_ctl,net_lhe_seasonal_avg_ctl,net_lhe_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[toa_sw_ctl,toa_sw_avg_ctl,toa_sw_seasonal_avg_ctl,toa_sw_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'toa_sw','W/m^2',factor = 1.) # positive DOWN

[precipitation_ctl,precipitation_avg_ctl,precipitation_seasonal_avg_ctl,precipitation_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'precipitation','mm/d', factor=86400)
[bucket_depth_ctl,bucket_depth_avg_ctl,bucket_depth_seasonal_avg_ctl,bucket_depth_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'bucket_depth','m')
[rh_ctl,rh_avg_ctl,rh_seasonal_avg_ctl,rh_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'rh','%',level=level)
# [flux_oceanq_ctl,flux_oceanq_avg_ctl,flux_oceanq_seasonal_avg_ctl,flux_oceanq_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_oceanq','W/m^2')
[net_sw_ctl,net_sw_avg_ctl,net_sw_seasonal_avg_ctl,net_sw_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_sw','W/m^2',factor = 1.) # 
[lw_down_ctl,lw_down_avg_ctl,lw_down_seasonal_avg_ctl,lw_down_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_lw','W/m^2',factor = 1.) # 
[net_t_ctl,net_t_avg_ctl,net_t_seasonal_avg_ctl,net_t_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_t','W/m^2',factor = 1.) # 


[sphum,sphum_avg,sphum_seasonal_avg,sphum_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'sphum','kg/kg',level=level)
[sphum_ctl,sphum_avg_ctl,sphum_seasonal_avg_ctl,sphum_month_avg_ctl,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'sphum','kg/kg',level=level)


#############################

# Enter control directory name as string 'ISCA_HPC/square_South_America_newbucket_0qflux_samealbedo_samehcp_commitfe93b9d'
# Enter runmin number 121
# Enter runmax number for comparison 481
# Enter end of ctl timeseries month 361
# Enter data directory name as string 'ISCA_HPC/square_South_America_newbucket_0qflux_samealbedo_samehcp_plus_2xCO2_spinup_361_commitfe93b9d'

[temp_ctl,temp_avg_ctl,temp_seasonal_avg_ctl,temp_month_avg_ctl,time]=seasonal_surface_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'temp','K', level = 2)
[sphum_ctl,sphum_avg_ctl,sphum_seasonal_avg_ctl,sphum_month_avg_ctl,time]=seasonal_surface_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'sphum','kg/kg', level = 2)
[height_ctl,height_avg_ctl,height_seasonal_avg_ctl,height_month_avg_ctl,time]=seasonal_surface_variable_interp(control_dir,ctl_model,ctl_runmin,ctl_runmax,'height','m', level = 2)



mse_subcloud = mc.cp_air * temp_month_avg_ctl + mc.L * sphum_month_avg_ctl+ mc.grav * height_month_avg_ctl
mse_subcloud_zm_land = xr.DataArray(area_weighted_avg_4D(mse_subcloud,area_array_3Dtime,landmaskxr,'all_sfcs',minlat=-90.,maxlat=90, minlon = 20., maxlon = 40., axis = 2), coords = [mse_subcloud.month,mse_subcloud.lat], dims = ['month','lat'])
precip_zm_land = xr.DataArray(area_weighted_avg_4D(precipitation_month_avg_ctl,area_array_3Dtime,landmaskxr,'all_sfcs',minlat=-90.,maxlat=90, minlon = 20., maxlon = 40., axis = 2), coords = [mse_subcloud.month,mse_subcloud.lat], dims = ['month','lat'])
mse_subcloud_zm_ocean = xr.DataArray(area_weighted_avg_4D(mse_subcloud,area_array_3Dtime,landmaskxr,'ocean',minlat=-90.,maxlat=90, minlon = 0., maxlon = 360., axis = 2), coords = [mse_subcloud.month,mse_subcloud.lat], dims = ['month','lat'])
precip_zm_ocean = xr.DataArray(area_weighted_avg_4D(precipitation_month_avg_ctl,area_array_3Dtime,landmaskxr,'ocean',minlat=-90.,maxlat=90, minlon = 0., maxlon = 360., axis = 2), coords = [mse_subcloud.month,mse_subcloud.lat], dims = ['month','lat'])

maxlat_land = temp_ctl.lat[np.argmax(mse_subcloud_zm_land, axis = 1)]
maxlat_ocean = temp_ctl.lat[np.argmax(mse_subcloud_zm_ocean, axis = 1)]

T, Y = np.meshgrid(mse_subcloud.month, mse_subcloud.lat.sel(lat=slice(-30.,30.)))

v = np.linspace(0.,12.,13) # using different colourbar min/max/steps because it's seasonal P and not annual mean! 
fig, axes = plt.subplots(1,2, sharey = True, sharex = True, figsize = (20,10))
pp = axes[0].contourf(T, Y, (precip_zm_land.transpose()).sel(lat=slice(-30.,30.)),v,cmap = 'Blues', extend = 'max')
con = axes[0].contour(T, Y, (mse_subcloud_zm_land.transpose()).sel(lat=slice(-30.,30.))/(10**5), [3.0, 3.1, 3.2, 3.25, 3.28, 3.29, 3.3, 3.31], cmap = 'Greys_r')
axes[1].contourf(T, Y, (precip_zm_ocean.transpose()).sel(lat=slice(-30.,30.)),v,cmap = 'Blues', extend = 'max')
con2 = axes[1].contour(T, Y, (mse_subcloud_zm_ocean.transpose()).sel(lat=slice(-30.,30.))/(10**5), [3.0, 3.1, 3.2, 3.25, 3.28, 3.29, 3.3, 3.31], cmap = 'Greys_r')

# 3.25, 3.26, 3.27, 3.28, 3.29, 3.3, 3.31, 3.32, 3.33
axes[0].clabel(con, [3.1, 3.2, 3.25, 3.28, 3.29, 3.3, 3.31], fmt = '%1.2f', fontsize = 18)
axes[1].clabel(con2, [3.1, 3.2, 3.25, 3.28, 3.29, 3.3], fmt = '%1.2f', fontsize = 18)

axes[0].plot(mse_subcloud.month, maxlat_land, 'wx', markersize=10)
axes[1].plot(mse_subcloud.month, maxlat_ocean, 'wx', markersize=10)

cbar = fig.colorbar(pp,ax=axes) # ax = axes tells it to take space away from all the subplots. could adjust location by setting ax to axes[0,0] for example. 
cbar.set_label('Precipitation (mm/d)', fontsize = 22)
cbar.ax.tick_params(labelsize=22) 
# axes[0].set_xlabel('Month', fontsize = 22)
axes[0].set_xticks(np.linspace(1,12,12))
axes[0].set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'], fontsize = 18)
# axes[1].set_xlabel('Month', fontsize = 22)
axes[1].set_xticks(np.linspace(1,12,12))
axes[1].set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'], fontsize = 18)
axes[0].set_ylabel('Latitude ($^{\circ}$N)', fontsize = 22)
axes[0].set_title('(a) Land East', fontsize = 22)
axes[1].set_title('(b) Ocean', fontsize = 22)
plt.rcParams['ytick.labelsize']=22

# axes[0].plot(mse_subcloud.month, (np.ones(12))*30, 'k--')
# axes[0].plot(mse_subcloud.month, (np.ones(12))*-30, 'k--')
fig.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/P_ctl_mse_subc_zonmean_40s40n_landEast_120-480.png', dpi=100, bbox_inches = 'tight')
fig.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/P_ctl_mse_subc_zonmean_40s40n_landEast_120-480.pdf', dpi=400, bbox_inches = 'tight')

# mse_subcloud_seasons = mc.cp_air * temp_seasonal_avg_ctl + mc.L * sphum_seasonal_avg_ctl+ mc.grav * height_seasonal_avg_ctl
# any_configuration_plot_seasonal(outdir,runmin,runmax,-90.,90.,precipitation_seasonal_avg_ctl,area_array,'mm/day','P_ctl_mseb_conts','fromwhite',landmaskxr,nmb_contours = [3.0, 3.1, 3.2, 3.25, 3.28, 3.29, 3.3, 3.31], minval = 0., maxval = 12., steps = 13, array2=mse_subcloud_seasons/(10**5))
# any_configuration_plot_allmonths(outdir,runmin,runmax,-90.,90.,precipitation_month_avg_ctl,area_array,'mm/day','P_ctl','fromwhite',landmaskxr, minval = 0., maxval = 12., steps = 13)

# any_configuration_plot_seasonal(outdir,runmin,runmax,-90.,90.,mse_subcloud_seasons/(10**5),area_array,'x10^5 J/kg','mseb','bucket',landmaskxr, minval = 3.0, maxval = 3.4,steps = 31)

sens_subcloud_zm_land = xr.DataArray(area_weighted_avg_4D(mc.cp_air*temp_month_avg_ctl/(10**5),area_array_3Dtime,landmaskxr,'all_sfcs',minlat=-90.,maxlat=90, minlon = 20., maxlon = 40., axis = 2), coords = [mse_subcloud.month,mse_subcloud.lat], dims = ['month','lat'])
sens_subcloud_zm_ocean = xr.DataArray(area_weighted_avg_4D(mc.cp_air*temp_month_avg_ctl/(10**5),area_array_3Dtime,landmaskxr,'ocean',minlat=-90.,maxlat=90, minlon = 0., maxlon = 360., axis = 2), coords = [mse_subcloud.month,mse_subcloud.lat], dims = ['month','lat'])
latent_subcloud_zm_land = xr.DataArray(area_weighted_avg_4D(mc.L*sphum_month_avg_ctl/(10**4),area_array_3Dtime,landmaskxr,'all_sfcs',minlat=-90.,maxlat=90, minlon = 20., maxlon = 40., axis = 2), coords = [mse_subcloud.month,mse_subcloud.lat], dims = ['month','lat'])
latent_subcloud_zm_ocean = xr.DataArray(area_weighted_avg_4D(mc.L*sphum_month_avg_ctl/(10**4),area_array_3Dtime,landmaskxr,'ocean',minlat=-90.,maxlat=90, minlon = 0., maxlon = 360., axis = 2), coords = [mse_subcloud.month,mse_subcloud.lat], dims = ['month','lat'])
height_subcloud_zm_land = xr.DataArray(area_weighted_avg_4D(mc.grav*height_month_avg_ctl/(10**4),area_array_3Dtime,landmaskxr,'all_sfcs',minlat=-90.,maxlat=90, minlon = 20., maxlon = 40., axis = 2), coords = [mse_subcloud.month,mse_subcloud.lat], dims = ['month','lat'])
height_subcloud_zm_ocean = xr.DataArray(area_weighted_avg_4D(mc.grav*height_month_avg_ctl/(10**4),area_array_3Dtime,landmaskxr,'ocean',minlat=-90.,maxlat=90, minlon = 0., maxlon = 360., axis = 2), coords = [mse_subcloud.month,mse_subcloud.lat], dims = ['month','lat'])
# bucket_depth_zm_land = xr.DataArray(area_weighted_avg_4D(bucket_depth_month_avg_ctl/(0.15*0.75),area_array_3Dtime,landmaskxr,'all_sfcs',minlat=-90.,maxlat=90, minlon = 20., maxlon = 40., axis = 2), coords = [mse_subcloud.month,mse_subcloud.lat], dims = ['month','lat'])

fig, axes = plt.subplots(2,2, sharey = True, figsize = (20,15))
plt.rcParams['ytick.labelsize']=22
v1 = np.arange(2.78,3.00,0.02)
v2 = np.arange(0.60,3.3,0.3)
v3 = np.arange(1.16,1.27,0.01)
v4 = np.arange(3.00,3.4,0.04)
sns = axes[0,0].contourf(T, Y, (sens_subcloud_zm_land.transpose()).sel(lat=slice(-30.,30.)),v1, cmap = 'Reds', extend = 'both')
lat = axes[0,1].contourf(T, Y, (latent_subcloud_zm_land.transpose()).sel(lat=slice(-30.,30.)),v2, cmap = 'Purples', extend = 'both')
hgt = axes[1,0].contourf(T, Y, (height_subcloud_zm_land.transpose()).sel(lat=slice(-30.,30.)),v3, cmap = 'Greens', extend = 'both')
mseall = axes[1,1].contourf(T, Y, ((mse_subcloud_zm_land*10**-5).transpose()).sel(lat=slice(-30.,30.)),v4, cmap = 'Blues', extend = 'both')
cbar = fig.colorbar(sns,ax=axes[0,0],orientation = 'horizontal', format = '%.2f') # ax = axes tells it to take space away from all the subplots. could adjust location by setting ax to axes[0,0] for example. 
cbar.set_label('x10$^5$ (J/kg)', fontsize = 22)
cbar.ax.tick_params(labelsize=22) 
cbar = fig.colorbar(lat,ax=axes[0,1],orientation = 'horizontal', format = '%.2f') # ax = axes tells it to take space away from all the subplots. could adjust location by setting ax to axes[0,0] for example. 
cbar.set_label('x10$^4$ (J/kg)', fontsize = 22)
cbar.ax.tick_params(labelsize=22) 
cbar = fig.colorbar(hgt,ax=axes[1,0],orientation = 'horizontal', format = '%.2f') # ax = axes tells it to take space away from all the subplots. could adjust location by setting ax to axes[0,0] for example. 
cbar.set_label('x10$^4$ (J/kg)', fontsize = 22)
cbar.ax.tick_params(labelsize=22) 
cbar = fig.colorbar(mseall,ax=axes[1,1],orientation = 'horizontal', format = '%.2f') # ax = axes tells it to take space away from all the subplots. could adjust location by setting ax to axes[0,0] for example. 
cbar.set_label('x10$^5$ (J/kg)', fontsize = 22)
cbar.ax.tick_params(labelsize=22) 
axes[0,0].set_xticks(np.linspace(1,12,12))
axes[0,1].set_xticks(np.linspace(1,12,12))
axes[1,0].set_xticks(np.linspace(1,12,12))
axes[1,1].set_xticks(np.linspace(1,12,12))
axes[0,0].set_title('(a) Sensible Heat', fontsize = 28)
axes[0,1].set_title('(b) Latent Heat', fontsize = 28)
axes[1,0].set_title('(c) Potential Energy', fontsize = 28)
axes[1,1].set_title('(d) MSE$_b$', fontsize = 28)
axes[1,0].set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'], fontsize = 22)
axes[1,1].set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'], fontsize = 22)
axes[0,0].set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'], fontsize = 22)
axes[0,1].set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'], fontsize = 22)
fig.text(0.05, 0.5, 'Latitude ($^{\circ}$N)', va='center', rotation='vertical', fontsize = 28)
# fig.text(0.45, 0.05, 'Month', rotation='horizontal', fontsize = 28)

fig.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/MSEb_terms_40s40n_landEast_cbartest_120-360.png', dpi=100, bbox_inches = 'tight')
fig.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/MSEb_terms_40s40n_landEast_cbartest_120-360.pdf', bbox_inches = 'tight')


fig, axes = plt.subplots(2,2, sharey = True, figsize = (20,15))
plt.rcParams['ytick.labelsize']=22
v1 = np.arange(2.78,3.00,0.02)
v2 = np.arange(0.60,3.3,0.3)
v3 = np.arange(1.16,1.27,0.01)
v4 = np.arange(3.00,3.4,0.04)
sns = axes[0,0].contourf(T, Y, (sens_subcloud_zm_ocean.transpose()).sel(lat=slice(-30.,30.)),v1,cmap = 'Reds', extend = 'both')
lat = axes[0,1].contourf(T, Y, (latent_subcloud_zm_ocean.transpose()).sel(lat=slice(-30.,30.)),v2,cmap = 'Purples', extend = 'both')
hgt = axes[1,0].contourf(T, Y, (height_subcloud_zm_ocean.transpose()).sel(lat=slice(-30.,30.)),v3,cmap = 'Greens', extend = 'both')
mseall = axes[1,1].contourf(T, Y, ((mse_subcloud_zm_ocean*10**-5).transpose()).sel(lat=slice(-30.,30.)),v4,cmap = 'Blues', extend = 'both')
cbar = fig.colorbar(sns,ax=axes[0,0],orientation = 'horizontal', format = '%.2f') # ax = axes tells it to take space away from all the subplots. could adjust location by setting ax to axes[0,0] for example. 
cbar.set_label('x10$^5$ (J/kg)', fontsize = 22)
cbar.ax.tick_params(labelsize=22) 
cbar = fig.colorbar(lat,ax=axes[0,1],orientation = 'horizontal', format = '%.2f') # ax = axes tells it to take space away from all the subplots. could adjust location by setting ax to axes[0,0] for example. 
cbar.set_label('x10$^4$ (J/kg)', fontsize = 22)
cbar.ax.tick_params(labelsize=22) 
cbar = fig.colorbar(hgt,ax=axes[1,0],orientation = 'horizontal', format = '%.2f') # ax = axes tells it to take space away from all the subplots. could adjust location by setting ax to axes[0,0] for example. 
cbar.set_label('x10$^4$ (J/kg)', fontsize = 22)
cbar.ax.tick_params(labelsize=22) 
cbar = fig.colorbar(mseall,ax=axes[1,1],orientation = 'horizontal', format = '%.2f') # ax = axes tells it to take space away from all the subplots. could adjust location by setting ax to axes[0,0] for example. 
cbar.set_label('x10$^5$ (J/kg)', fontsize = 22)
cbar.ax.tick_params(labelsize=22) 
axes[0,0].set_xticks(np.linspace(1,12,12))
axes[0,1].set_xticks(np.linspace(1,12,12))
axes[1,0].set_xticks(np.linspace(1,12,12))
axes[1,1].set_xticks(np.linspace(1,12,12))
axes[0,0].set_title('(a) Sensible Heat', fontsize = 28)
axes[0,1].set_title('(b) Latent Heat', fontsize = 28)
axes[1,0].set_title('(c) Potential Energy', fontsize = 28)
axes[1,1].set_title('(d) MSE$_b$', fontsize = 28)
axes[1,0].set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'], fontsize = 22)
axes[1,1].set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'], fontsize = 22)
axes[0,0].set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'], fontsize = 22)
axes[0,1].set_xticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'], fontsize = 22)
fig.text(0.05, 0.5, 'Latitude ($^{\circ}$N)', va='center', rotation='vertical', fontsize = 28)
# fig.text(0.45, 0.05, 'Month', rotation='horizontal', fontsize = 28)

fig.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/MSEb_terms_40s40n_ocean_cbartest_120-360.png', dpi=100, bbox_inches = 'tight')
fig.savefig('/scratch/mp586/Code/Graphics/'+outdir+'/MSEb_terms_40s40n_ocean_cbartest_120-360.pdf', bbox_inches = 'tight')
