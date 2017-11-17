from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import xarray as xr
import pandas as pd
import os

import sys
sys.path.insert(0, '/scratch/mp586/Code/PYCODES')
import plotting_routines
import plotting_routines_kav7
import stats as st

control_dir= input('Enter control directory name as string ')
ctl_runmin=input('Enter runmin number ')  # Should be a January month for seasonal variables to be correct
ctl_runmax=input('Enter runmax number ')

testdir= input('Enter data directory name as string ')
runmin=input('Enter runmin number ')  # Should be a January month for seasonal variables to be correct
runmax=input('Enter runmax number ')


# landfile=Dataset('/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/input/two_continents/land_two_continents.nc',mode='r')
landfile=Dataset('/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/input/squareland/land_square.nc',mode='r')
# landfile=Dataset('/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/input/sqland_plus_antarctica/land_sqland_plus_antarctica.nc',mode='r')
# landfile=Dataset('/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/input/sqland_plus_antarctica/land_sqland_plus_antarctica_to35S.nc',mode='r')
# landfile=Dataset('/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/input/aquaplanet/land_aquaplanet.nc',mode='r')
# landfile=Dataset('/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/input/square_South_America/land_square_South_America.nc')

landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]


# for specified lats
landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

#plotting_routines_kav7.globavg_var_timeseries(testdir,'t_surf',109,122)

# for plotting a spin up run ('control') timeseries followed by the timeseries from the perturbed experiment
plotting_routines_kav7.globavg_var_timeseries_total_and_land_perturbed(testdir,'t_surf',1,runmax,1.,landmask,control_dir,1,ctl_runmax)
plotting_routines_kav7.globavg_var_timeseries_total_and_land_perturbed(testdir,'bucket_depth',1,runmax,1.,landmask,control_dir,1,ctl_runmax)


# plotting_routines_kav7.globavg_var_timeseries(testdir,'co2',1,runmax)


[tsurf,tsurf_avg,tsurf_seasonal_avg,tsurf_month_avg,time]=plotting_routines_kav7.seasonal_surface_variable(testdir,runmin,runmax,'t_surf','K')
[net_lhe,net_lhe_avg,net_lhe_seasonal_avg,net_lhe_month_avg,time]=plotting_routines_kav7.seasonal_surface_variable(testdir,runmin,runmax,'flux_lhe','W/m^2') # latent heat flux at surface (UP)
[precipitation,precipitation_avg,precipitation_seasonal_avg,precipitation_month_avg,time]=plotting_routines_kav7.seasonal_surface_variable(testdir,runmin,runmax,'precipitation','kg/m2s')
[convection_rain,convection_rain_avg,convection_rain_seasonal_avg,convection_rain_month_avg,time]=plotting_routines_kav7.seasonal_surface_variable(testdir,runmin,runmax,'convection_rain','kg/m2s')
[condensation_rain,condensation_rain_avg,condensation_rain_seasonal_avg,condensation_rain_month_avg,time]=plotting_routines_kav7.seasonal_surface_variable(testdir,runmin,runmax,'condensation_rain','kg/m2s')
#[bucket_depth,bucket_depth_avg,bucket_depth_seasonal_avg,bucket_depth_month_avg,time]=plotting_routines_kav7.seasonal_surface_variable(testdir,runmin,runmax,'bucket_depth','m')

[tsurf_ctl,tsurf_avg_ctl,tsurf_seasonal_avg_ctl,tsurf_month_avg_ctl,time_ctl]=plotting_routines_kav7.seasonal_surface_variable(control_dir,ctl_runmin,ctl_runmax,'t_surf','K')
[net_lhe_ctl,net_lhe_avg_ctl,net_lhe_seasonal_avg_ctl,net_lhe_month_avg_ctl,time_ctl]=plotting_routines_kav7.seasonal_surface_variable(control_dir,ctl_runmin,ctl_runmax,'flux_lhe','W/m^2') # latent heat flux at surface (UP)
[precipitation_ctl,precipitation_avg_ctl,precipitation_seasonal_avg_ctl,precipitation_month_avg_ctl,time_ctl]=plotting_routines_kav7.seasonal_surface_variable(control_dir,ctl_runmin,ctl_runmax,'precipitation','kg/m2s')




# [ucomp,ucomp_avg,ucomp_seasonal_avg,ucomp_month_avg,time]=plotting_routines_kav7.seasonal_4D_variable(testdir,runmin,runmax,'ucomp','m/s')
# [vcomp,vcomp_avg,vcomp_seasonal_avg,vcomp_month_avg,time]=plotting_routines_kav7.seasonal_4D_variable(testdir,runmin,runmax,'vcomp','m/s')
# [omega,omega_avg,omega_seasonal_avg,omega_month_avg,time]=plotting_routines_kav7.seasonal_4D_variable(testdir,runmin,runmax,'omega','Pa/s')


# [ucomp_ctl,ucomp_avg_ctl,ucomp_seasonal_avg_ctl,ucomp_month_avg_ctl,time]=plotting_routines_kav7.seasonal_4D_variable(control_dir,ctl_runmin,ctl_runmax,'ucomp','m/s')
# [vcomp_ctl,vcomp_avg_ctl,vcomp_seasonal_avg_ctl,vcomp_month_avg_ctl,time]=plotting_routines_kav7.seasonal_4D_variable(control_dir,ctl_runmin,ctl_runmax,'vcomp','m/s')
# [omega_ctl,omega_avg_ctl,omega_seasonal_avg_ctl,omega_month_avg_ctl,time]=plotting_routines_kav7.seasonal_4D_variable(control_dir,ctl_runmin,ctl_runmax,'omega','Pa/s')



# maxval_precip = np.absolute((precipitation_month_avg*86400).max())

# for i in range(1,13):

#     month_plot = plotting_routines_kav7.winds_at_heightlevel(ucomp_month_avg.sel(month=i),vcomp_month_avg.sel(month=i),39,precipitation_month_avg.sel(month=i)*86400,'fromwhite','mm/day',0,maxval_precip)
#     month_plot.savefig('/scratch/mp586/Code/Graphics/'+testdir+'/anim_plot_precip'+str(i)+'.png',bbox_inches='tight')
# os.system('convert -delay 100 /scratch/mp586/Code/Graphics/'+testdir+'/anim_plot_precip*.png /scratch/mp586/Code/Graphics/'+testdir+'/precip_wind_monthly_clim.gif')

# maxval_omega_surf = np.absolute((omega_month_avg[:,39,:,:]).max())
# minval_omega_surf = np.absolute((omega_month_avg[:,39,:,:]).min())


# for i in range(1,13):

#     month_plot = plotting_routines_kav7.winds_at_heightlevel(ucomp_month_avg.sel(month=i),vcomp_month_avg.sel(month=i),39,(omega_month_avg[:,39,:,:]).sel(month=i),'rainnorm','mm/day',minval_omega_surf,maxval_omega_surf)
#     month_plot.savefig('/scratch/mp586/Code/Graphics/'+testdir+'/anim_plot_omega'+str(i)+'.png',bbox_inches='tight')
# os.system('convert -delay 100 /scratch/mp586/Code/Graphics/'+testdir+'/anim_plot_omega*.png /scratch/mp586/Code/Graphics/'+testdir+'/omega_wind_monthly_clim.gif')


# JJA = 'JJA'
# DJF = 'DJF'
# MAM = 'MAM'
# SON = 'SON'

# summer_plot = plotting_routines_kav7.winds_at_heightlevel(ucomp_seasonal_avg.sel(season=JJA),vcomp_seasonal_avg.sel(season=JJA),39,precipitation_seasonal_avg.sel(season=JJA)*86400,'fromwhite','mm/day')
# summer_plot.savefig('anim_plot1.png',bbox_inches='tight')
# fall_plot = plotting_routines_kav7.winds_at_heightlevel(ucomp_seasonal_avg.sel(season=SON),vcomp_seasonal_avg.sel(season=SON),39,precipitation_seasonal_avg.sel(season=SON)*86400,'fromwhite','mm/day')
# fall_plot.savefig('anim_plot2.png',bbox_inches='tight')
# winter_plot = plotting_routines_kav7.winds_at_heightlevel(ucomp_seasonal_avg.sel(season=DJF),vcomp_seasonal_avg.sel(season=DJF),39,precipitation_seasonal_avg.sel(season=DJF)*86400,'fromwhite','mm/day')
# winter_plot.savefig('anim_plot3.png',bbox_inches='tight')
# spring_plot = plotting_routines_kav7.winds_at_heightlevel(ucomp_seasonal_avg.sel(season=MAM),vcomp_seasonal_avg.sel(season=MAM),39,precipitation_seasonal_avg.sel(season=MAM)*86400,'fromwhite','mm/day')
# spring_plot.savefig('anim_plot4.png',bbox_inches='tight')

# os.system('convert -delay 50 anim_plot*.png animation_seasons.gif')



#plotting_routines_kav7.animated_map(testdir,bucket_depth.where(landmask==1.),'m','bucket depth','bucket_depth','fromwhite',0,140) # need runmin = 1!


PE_avg=precipitation_avg*86400-net_lhe_avg/28. # 28.=conversion from W/m^# 2 to mm/day using E=H/(rho*L), rho=1000kg/m3, L=2.5*10^6J/kg
# # # see www.ce.utexas.edu/prof/maidment/CE374KSpr12/.../Latent%20heat%20flux.pptx @30DegC


# # # plotting_routines_kav7.any_configuration_plot(-90.,90.,convection_rain_avg*86400,'mm/day','convection_rain avg','fromwhite')
# # # plotting_routines_kav7.any_configuration_plot(-90.,90.,condensation_rain_avg*86400,'mm/day','condensation_rain avg','fromwhite')
# #plotting_routines_kav7.any_configuration_plot(-90.,90.,precipitation_avg.where(landmask==1.)*86400,'mm/day','P avg','fromwhite')
# # # #plotting_routines_kav7.any_configuration_plot(-90.,90.,bucket_depth_avg.where(landmask==1.),'m','bucket_depth','fromwhite')
# # # plotting_routines_kav7.any_configuration_plot(-90.,90.,net_lhe_avg.where(landmask==1.)/28.,'mm/day','E avg','fromwhite')
# plotting_routines_kav7.any_configuration_plot(-100.,100.,PE_avg,'mm/day','P-E avg','rainnorm',landmask,landlats,landlons)
# # # #plotting_routines_kav7.any_configuration_plot_minuszonavg(-90.,90.,PE_avg,'mm/day','P-E avg minus zonavg','rainnorm','P-E avg')
# plotting_routines_kav7.any_configuration_plot(-90.,90.,tsurf_avg,'K','$T_S$ avg','temp',landmask,landlats,landlons) # degrees C symbol : ...,u"\u00b0"+'C',...
plotting_routines_kav7.any_configuration_plot(-90.,90.,(tsurf_avg-tsurf_avg_ctl),'K','$T_S$ avg minus ctrl','tempdiff',landmask,landlats,landlons,contourson = True)
# # #plotting_routines_kav7.any_configuration_plot_minuszonavg(-90.,90.,tsurf_avg,'K','tsurf avg minus zonavg','temp','T avg')
plotting_routines_kav7.any_configuration_plot(-90.,90.,precipitation_avg*86400,'mm/day','P avg','fromwhite',landmask,landlats,landlons,contourson = True)
plotting_routines_kav7.any_configuration_plot(-90.,90.,(precipitation_avg - precipitation_avg_ctl)*86400,'mm/day','P avg minus ctrl','rainnorm',landmask,landlats,landlons,contourson = False)
plotting_routines_kav7.any_configuration_plot(-90.,90.,(net_lhe_avg - net_lhe_avg_ctl)/28.,'mm/day','E avg minus ctrl','rainnorm',landmask,landlats,landlons,contourson = False)
# # plotting_routines_kav7.any_configuration_plot(-90.,90.,(omega_avg[39,:,:] - omega_avg_ctl[39,:,:]),'Pa/s','Omega avg minus ctrl','rainnorm',landmask,landlats,landlons,contourson = False)
# # # #plotting_routines_kav7.any_configuration_plot_minuszonavg(-90.,90.,precipitation_avg*86400,'mm/day','P avg minus zonavg','rainnorm','P avg')
# # plotting_routines_kav7.any_configuration_plot(-90.,90.,net_lhe_avg/28.,'mm/day','E avg','fromwhite')

land_temp_global=tsurf_avg.where(landmask==1.).mean()
ocean_temp_global=tsurf_avg.where(landmask==0.).mean()
print('Average temperature over land (global) = '+str(land_temp_global))
print('Average temperature over ocean (global) = '+str(ocean_temp_global))


minlat=-30. #input('Enter minimum latitude ')
maxlat=30. #input('Enter maximum latitude ')

# difference between exp and control
land_temp=(tsurf_avg-tsurf_avg_ctl).sel(lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==1.) # .where does not recognize xarray as argument!
ocean_temp=(tsurf_avg-tsurf_avg_ctl).sel(lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==0.) # .where does not recognize xarray as argument!

print('Average temperature diff over land between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(land_temp.mean()))
print('Average temperature diff over ocean between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(ocean_temp.mean()))

# difference between exp and control
land_precip_global=(precipitation_avg-precipitation_avg_ctl).where(landmask==1.).mean() 
ocean_precip_global=(precipitation_avg-precipitation_avg_ctl).where(landmask==0.).mean()
print('Average precipitation diff over land (global) = '+str(land_precip_global*86400)+' mm/day')
print('Average precipitation diff over ocean (global) = '+str(ocean_precip_global*86400)+' mm/day')

# difference between exp and control
land_precip=(precipitation_avg-precipitation_avg_ctl).sel(lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==1.) # .where does not recognize xarray as argument!
ocean_precip=(precipitation_avg-precipitation_avg_ctl).sel(lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==0.) # .where does not recognize xarray as argument!
print('Average precipitation diff over land between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(land_precip.mean()*86400)+'+/-'+str(land_precip.std()*86400)+'mm/day') # spatial standard deviation
print('Average precipitation diff over ocean between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(ocean_precip.mean()*86400)+'+/-'+str(ocean_precip.std()*86400)+'mm/day')
print('Min precipitation diff over land between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(land_precip.min()*86400)+'mm/day') # spatial standard deviation
print('Max precipitation diff over land between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(land_precip.max()*86400)+'mm/day') # spatial standard deviation
print('Min precipitation diff over ocean between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(ocean_precip.min()*86400)+'mm/day') # spatial standard deviation
print('Max precipitation diff over ocean between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(ocean_precip.max()*86400)+'mm/day')

# difference between exp and control
land_lhe=(net_lhe_avg-net_lhe_avg_ctl).sel(lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==1.) # .where does not recognize xarray as argument!
ocean_lhe=(net_lhe_avg-net_lhe_avg_ctl).sel(lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==0.) # .where does not recognize xarray as argument!
print('Average net_lhe diff over land between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(land_lhe.mean()/28.)+'+/-'+str(land_lhe.std()/28.)+'mm/day') # spatial standard deviation
print('Average net_lhe diff over ocean between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(ocean_lhe.mean()/28.)+'+/-'+str(ocean_lhe.std()/28.)+'mm/day')
print('Min net_lhe diff over land between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(land_lhe.min()/28.)+'mm/day') # spatial standard deviation
print('Max net_lhe diff over land between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(land_lhe.max()/28.)+'mm/day') # spatial standard deviation
print('Min net_lhe diff over ocean between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(ocean_lhe.min()/28.)+'mm/day') # spatial standard deviation
print('Max net_lhe diff over ocean between '+str(minlat)+'N and '+str(maxlat)+'N = '+str(ocean_lhe.max()/28.)+'mm/day')




JJA = 'JJA'
DJF = 'DJF'
MAM = 'MAM'
SON = 'SON'

land_precip_JJA=precipitation_seasonal_avg.sel(season=JJA,lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==1.) # .where does not recognize xarray as argument!
ocean_precip_JJA=precipitation_seasonal_avg.sel(season=JJA,lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==0.) # .where does not recognize xarray as argument!
land_temp_JJA=tsurf_seasonal_avg.sel(season=JJA,lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==1.) # .where does not recognize xarray as argument!
ocean_temp_JJA=tsurf_seasonal_avg.sel(season=JJA,lat=slice(minlat,maxlat)).where(np.asarray(landmaskxr.sel(lat=slice(minlat,maxlat)))==0.) # .where does not recognize xarray as argument!

plotting_routines_kav7.several_vars_zonalavg2(tsurf_seasonal_avg.sel(season=JJA),'T_S (K)','Red',precipitation_seasonal_avg.sel(season=JJA)*86400,'P (mm/day)','Blue',net_lhe_seasonal_avg.sel(season=JJA)/28.,'E (mm/day)','g','JJA T_S, P and E') # 28.=conversion from W/m^2 to mm/day using E=H/(rho*L), rho=1000kg/m3, L=2.5*10^6J/kg
plotting_routines_kav7.several_vars_zonalavg2(tsurf_seasonal_avg.sel(season=DJF),'T_S (K)','Red',precipitation_seasonal_avg.sel(season=DJF)*86400,'P (mm/day)','Blue',net_lhe_seasonal_avg.sel(season=DJF)/28.,'E (mm/day)','g','DJF T_S, P and E')
plotting_routines_kav7.several_vars_zonalavg2(tsurf_avg,'T_S (K)','Red',precipitation_avg*86400,'P (mm/day)','Blue',net_lhe_avg/28.,'E (mm/day)','g','avg T_S, P and E')


plotting_routines_kav7.any_configuration_plot(-90.,90.,precipitation_month_avg.sel(month=7)*86400,'mm/day','P_July (mm/day)','rainnorm',landmask,landlats,landlons)
plotting_routines_kav7.any_configuration_plot(-90.,90.,precipitation_month_avg.sel(month=1)*86400,'mm/day','P_January (mm/day)','rainnorm',landmask,landlats,landlons)
plotting_routines_kav7.any_configuration_plot(-90.,90.,tsurf_month_avg.sel(month=7),'K','tsurf_July (K)','temp',landmask,landlats,landlons)
plotting_routines_kav7.any_configuration_plot(-90.,90.,tsurf_month_avg.sel(month=1),'K','tsurf_January (K)','temp',landmask,landlats,landlons)

plotting_routines_kav7.any_configuration_plot(-90.,90.,precipitation_month_avg.sel(month=7)*86400-net_lhe_month_avg.sel(month=7)/28.,'mm/day','P-E_July (mm/day)','rainnorm',landmask,landlats,landlons)
plotting_routines_kav7.any_configuration_plot(-90.,90.,precipitation_month_avg.sel(month=1)*86400-net_lhe_month_avg.sel(month=1)/28.,'mm/day','P-E_January (mm/day)','rainnorm',landmask,landlats,landlons)


print('JJA tsurf_avg (global) = '+str(tsurf_seasonal_avg.sel(season=JJA).mean()))
print('DJF tsurf_avg (global) = '+str(tsurf_seasonal_avg.sel(season=DJF).mean()))