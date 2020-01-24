from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import xarray as xr
import pandas as pd

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
control_dir = 'full_continents_newbucket_fixedSSTs_from_realworld_zonallysymm_commitfe93b9d'
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
testdir_in1= 'full_continents_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commitfe93b9d'
runmin=120
runmax=480
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    exp1_name = 'ISCA_HPC_'+testdir_in1
    testdir = model_data + '/ISCA_HPC/' + testdir_in1
    testdir_in1 = '/ISCA_HPC/' + testdir_in1
else: 
    exp1_name = testdir_in1
    testdir = model_data + '/' + testdir_in1

land = 'all_continents'
landfile=Dataset(os.path.join(GFDL_BASE,'input/'+land+'/land.nc'),mode='r')

landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]
# for specified lats
landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))

landmask, lons_cyclic = addcyclic(landmask, landlons)

landmaskxr = xr.DataArray(landmask,coords=[landlats,lons_cyclic],dims=['lat','lon'])



area_array, dx, dy = ca.cell_area_all(t_res=42,base_dir='/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/') # added _all because then dx and dy are also returned 
area_array = xr.DataArray(area_array) # returned in units of m bzw m^2, because radius in cell_area.py is given in metres
area_array = np.asarray(area_array)
area_array, lons_cyclic = addcyclic(area_array, landlons)
area_array,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,area_array,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))

area_array = xr.DataArray(area_array,coords=[landlats,lons_cyclic],dims=['lat','lon'])



[net_lhe,net_lhe_avg,net_lhe_seasonal_avg,net_lhe_month_avg,time]=seasonal_surface_variable_shift(testdir,model,runmin,runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[precipitation,precipitation_avg,precipitation_seasonal_avg,precipitation_month_avg,time]=seasonal_surface_variable_shift(testdir,model,runmin,runmax,'precipitation','mm/d', factor=86400)
[sw_net,sw_net_avg,sw_net_seasonal_avg,sw_net_month_avg,time]=seasonal_surface_variable_shift(testdir,model,runmin,runmax,'flux_sw','W/m^2',factor = 1.) # basically assume that all energy from Rnet goes into evaporation
[tsurf,tsurf_avg,tsurf_seasonal_avg,tsurf_month_avg,time]=seasonal_surface_variable_shift(testdir,model,runmin,runmax,'t_surf','K')
[lw_down,lw_down_avg,lw_down_seasonal_avg,lw_down_month_avg,time]=seasonal_surface_variable_shift(testdir,model,runmin,runmax,'flux_lw','W/m^2',factor = 1.) # 
[ep,ep_avg,ep_seasonal_avg,ep_month_avg,time]=seasonal_surface_variable_shift(testdir,model,runmin,runmax,'potential_evap','mm/day',factor = 86400) # latent heat flux at surface (UP)

sigma = 5.67*10**(-8)

net_lw_avg = sigma*(tsurf_avg**4) - lw_down_avg

# Ep = Rnet/L, to convert from flux_lhe to mm/day basically convert W/mz to mm/day... factor = 1./28. 
# 1/28.=conversion from W/m^2 to mm/day using E=H/(rho*L), rho=1000kg/m3, L=2.5*10^6J/kg, see www.ce.utexas.edu/prof/maidment/CE374KSpr12/.../Latent%20heat%20flux.pptx @30DegC


[net_lhe_ctl,net_lhe_avg_ctl,net_lhe_seasonal_avg_ctl,net_lhe_month_avg_ctl,time]=seasonal_surface_variable_shift(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[precipitation_ctl,precipitation_avg_ctl,precipitation_seasonal_avg_ctl,precipitation_month_avg_ctl,time]=seasonal_surface_variable_shift(control_dir,ctl_model,ctl_runmin,ctl_runmax,'precipitation','mm/d', factor=86400)
[sw_net_ctl,sw_net_avg_ctl,sw_net_seasonal_avg_ctl,sw_net_month_avg_ctl,time]=seasonal_surface_variable_shift(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_sw','W/m^2',factor = 1.) # latent heat flux at surface (UP)
[tsurf_ctl,tsurf_avg_ctl,tsurf_seasonal_avg_ctl,tsurf_month_avg_ctl,time]=seasonal_surface_variable_shift(control_dir,ctl_model,ctl_runmin,ctl_runmax,'t_surf','K')
[lw_down_ctl,lw_down_avg_ctl,lw_down_seasonal_avg_ctl,lw_down_month_avg_ctl,time]=seasonal_surface_variable_shift(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_lw','W/m^2',factor = 1.) # 
[ep_ctl,ep_avg_ctl,ep_seasonal_avg_ctl,ep_month_avg_ctl,time]=seasonal_surface_variable_shift(control_dir,ctl_model,ctl_runmin,ctl_runmax,'potential_evap','mm/day',factor = 86400) # latent heat flux at surface (UP)
# [net_t_ctl,net_t_avg_ctl,net_t_seasonal_avg_ctl,net_t_month_avg_ctl,time]=seasonal_surface_variable_shift(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_t','W/m^2',factor = 1.) # 

net_lw_avg_ctl = sigma*(tsurf_avg_ctl**4) - lw_down_avg_ctl
net_lw_ctl = sigma*(tsurf_ctl**4) - lw_down_ctl

# Ep_estim_ctl = (sw_net_avg_ctl - net_lw_avg_ctl - net_t_avg_ctl)/28. 
Ep_estim_avg_ctl = (sw_net_avg_ctl - net_lw_avg_ctl)/28. 
Ep_estim_ctl = (sw_net_ctl - net_lw_ctl)/28. 


testdir_in1= 'full_continents_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commitfe93b9d'
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

[net_lheveg,net_lheveg_avg,net_lheveg_seasonal_avg,net_lheveg_month_avg,time]=seasonal_surface_variable_shift(testdir,model,runmin,runmax,'flux_lhe','mm/d',factor = 1./28.) # latent heat flux at surface (UP)
[precipitationveg,precipitationveg_avg,precipitationveg_seasonal_avg,precipitationveg_month_avg,time]=seasonal_surface_variable_shift(testdir,model,runmin,runmax,'precipitation','mm/d', factor=86400)
# [sw_netveg,sw_netveg_avg,sw_netveg_seasonal_avg,sw_netveg_month_avg,time]=seasonal_surface_variable_shift(testdir,model,runmin,runmax,'flux_sw','mm/day',factor = 1.) # latent heat flux at surface (UP)
# [tsurfveg,tsurfveg_avg,tsurfveg_seasonal_avg,tsurfveg_month_avg,time]=seasonal_surface_variable_shift(testdir,model,runmin,runmax,'t_surf','K')
# [lw_downveg,lw_downveg_avg,lw_downveg_seasonal_avg,lw_downveg_month_avg,time]=seasonal_surface_variable_shift(testdir,model,runmin,runmax,'flux_lw','W/m^2',factor = 1.) # 
[epveg,epveg_avg,epveg_seasonal_avg,epveg_month_avg,time]=seasonal_surface_variable_shift(testdir,model,runmin,runmax,'potential_evap','mm/day',factor = 86400) # latent heat flux at surface (UP)


# finding omega-parameter

# america 

ep_ctl_yrs = ep_ctl.groupby('time.year').mean('time')
Ep_estim_ctl_yrs = Ep_estim_ctl.groupby('time.year').mean('time')
E_ctl_yrs = net_lhe_ctl.groupby('time.year').mean('time')
P_ctl_yrs = precipitation_ctl.groupby('time.year').mean('time')

rge = np.linspace(1.,5.,99)

obj = np.empty((len(rge)))

j = 0 

for omega in rge: # parameter sweep
    obj_years = np.empty(((np.shape(ep_ctl_yrs))[0]))
    for i in range((np.shape(ep_ctl_yrs))[0]):
        ep_ctl_yrs_awave = area_weighted_avg(ep_ctl_yrs[i,:,:], area_array, landmaskxr, 'land', minlat=-12.,maxlat=8., minlon = -70., maxlon = -45.)
        E_ctl_yrs_awave = area_weighted_avg(E_ctl_yrs[i,:,:], area_array, landmaskxr, 'land', minlat=-12.,maxlat=8., minlon = -70., maxlon = -45.)
        P_ctl_yrs_awave = area_weighted_avg(P_ctl_yrs[i,:,:], area_array, landmaskxr, 'land', minlat=-12.,maxlat=8., minlon = -70., maxlon = -45.)
        obj_years[i] = (E_ctl_yrs_awave/P_ctl_yrs_awave - (1 + ep_ctl_yrs_awave/P_ctl_yrs_awave - (1 + (ep_ctl_yrs_awave/P_ctl_yrs_awave)**omega)**(1./omega))**2)
    obj[j] = np.sum(obj_years)
    j += 1

omega_opt_AM = rge[np.argmin(np.abs(obj))]


obj = np.empty((len(rge)))

j = 0 

rge = np.linspace(1.,10.,99)

for omega in rge: # parameter sweep
    obj_years = np.empty(((np.shape(ep_ctl_yrs))[0]))
    for i in range((np.shape(ep_ctl_yrs))[0]):
        Ep_estim_ctl_yrs_awave = area_weighted_avg(Ep_estim_ctl_yrs[i,:,:], area_array, landmaskxr, 'land', minlat=-12.,maxlat=8., minlon = -70., maxlon = -45.)
        E_ctl_yrs_awave = area_weighted_avg(E_ctl_yrs[i,:,:], area_array, landmaskxr, 'land', minlat=-12.,maxlat=8., minlon = -70., maxlon = -45.)
        P_ctl_yrs_awave = area_weighted_avg(P_ctl_yrs[i,:,:], area_array, landmaskxr, 'land', minlat=-12.,maxlat=8., minlon = -70., maxlon = -45.)
        obj_years[i] = (E_ctl_yrs_awave/P_ctl_yrs_awave - (1 + Ep_estim_ctl_yrs_awave/P_ctl_yrs_awave - (1 + (Ep_estim_ctl_yrs_awave/P_ctl_yrs_awave)**omega)**(1./omega))**2)
    obj[j] = np.sum(obj_years)
    j += 1

omega_opt_AM_estim = rge[np.argmin(np.abs(obj))]


rge = np.linspace(1.,5.,99)

obj = np.empty((len(rge)))

j = 0 


for omega in rge: # parameter sweep
    obj_years = np.empty(((np.shape(ep_ctl_yrs))[0]))
    for i in range((np.shape(ep_ctl_yrs))[0]):
        ep_ctl_yrs_awave = area_weighted_avg(ep_ctl_yrs[i,:,:], area_array, landmaskxr, 'land',minlat=-18.,maxlat=18., minlon = -17., maxlon = 60.)
        E_ctl_yrs_awave = area_weighted_avg(E_ctl_yrs[i,:,:], area_array, landmaskxr, 'land', minlat=-18.,maxlat=18., minlon = -17., maxlon = 60.)
        P_ctl_yrs_awave = area_weighted_avg(P_ctl_yrs[i,:,:], area_array, landmaskxr, 'land', minlat=-18.,maxlat=18., minlon = -17., maxlon = 60.)
        obj_years[i] = (E_ctl_yrs_awave/P_ctl_yrs_awave - (1 + ep_ctl_yrs_awave/P_ctl_yrs_awave - (1 + (ep_ctl_yrs_awave/P_ctl_yrs_awave)**omega)**(1./omega))**2)
    obj[j] = np.sum(obj_years)
    j += 1

omega_opt_AF = rge[np.argmin(np.abs(obj))]


obj = np.empty((len(rge)))

j = 0 


for omega in rge: # parameter sweep
    obj_years = np.empty(((np.shape(ep_ctl_yrs))[0]))
    for i in range((np.shape(ep_ctl_yrs))[0]):
        Ep_estim_ctl_yrs_awave = area_weighted_avg(Ep_estim_ctl_yrs[i,:,:], area_array, landmaskxr, 'land', minlat=-18.,maxlat=18., minlon = -17., maxlon = 60.)
        E_ctl_yrs_awave = area_weighted_avg(E_ctl_yrs[i,:,:], area_array, landmaskxr, 'land',minlat=-18.,maxlat=18., minlon = -17., maxlon = 60.)
        P_ctl_yrs_awave = area_weighted_avg(P_ctl_yrs[i,:,:], area_array, landmaskxr, 'land', minlat=-18.,maxlat=18., minlon = -17., maxlon = 60.)
        obj_years[i] = (E_ctl_yrs_awave/P_ctl_yrs_awave - (1 + Ep_estim_ctl_yrs_awave/P_ctl_yrs_awave - (1 + (Ep_estim_ctl_yrs_awave/P_ctl_yrs_awave)**omega)**(1./omega))**2)
    obj[j] = np.sum(obj_years)
    j += 1

omega_opt_AF_estim = rge[np.argmin(np.abs(obj))]



# print('Amazon $\Delta P$ = '+str(area_weighted_avg(array,area_array,landmaskxr,'land',minlat=-12.,maxlat=8., minlon = -70., maxlon = -45.)))
# print('Africa $\Delta P$ = '+str(area_weighted_avg(array,area_array,landmaskxr,'land',minlat=-18.,maxlat=18., minlon = -17., maxlon = 60.)))


PAM_ctl = precipitation_avg_ctl.where(landmask == 1.).sel(lat = slice(-12., 8.)).sel(lon=slice(-70.,-45.))
EAM_ctl = net_lhe_avg_ctl.where(landmask == 1.).sel(lat = slice(-12., 8.)).sel(lon=slice(-70.,-45.))
epAM_ctl = ep_avg_ctl.where(landmask == 1.).sel(lat = slice(-12., 8.)).sel(lon=slice(-70.,-45.))
epAMveg = epveg_avg.where(landmask == 1.).sel(lat = slice(-12., 8.)).sel(lon=slice(-70.,-45.))
Ep_AM_estim_ctl = Ep_estim_avg_ctl.where(landmask == 1.).sel(lat = slice(-12., 8.)).sel(lon=slice(-70.,-45.))
epAM = ep_avg.where(landmask == 1.).sel(lat = slice(-12., 8.)).sel(lon=slice(-70.,-45.))
PAM = precipitation_avg.where(landmask == 1.).sel(lat = slice(-12., 8.)).sel(lon=slice(-70.,-45.))
EAM = net_lhe_avg.where(landmask == 1.).sel(lat = slice(-12., 8.)).sel(lon=slice(-70.,-45.))
PAMveg = precipitationveg_avg.where(landmask == 1.).sel(lat = slice(-12., 8.)).sel(lon=slice(-70.,-45.))
EAMveg = net_lheveg_avg.where(landmask == 1.).sel(lat = slice(-12., 8.)).sel(lon=slice(-70.,-45.))
fu_ctlAM = (epAM_ctl/PAM_ctl) + 1. - (1. + (epAM_ctl/PAM_ctl)**(omega_opt_AM))**(1./omega_opt_AM)
fu_ctlAM_estim = (Ep_AM_estim_ctl/PAM_ctl) + 1. - (1. + (Ep_AM_estim_ctl/PAM_ctl)**(omega_opt_AM_estim))**(1./omega_opt_AM_estim)

PAF_ctl = precipitation_avg_ctl.where(landmask == 1.).sel(lat = slice(-18., 18.)).sel(lon=slice(-17.,60.))
EAF_ctl = net_lhe_avg_ctl.where(landmask == 1.).sel(lat = slice(-18., 18.)).sel(lon=slice(-17.,60.))
epAF_ctl = ep_avg_ctl.where(landmask == 1.).sel(lat = slice(-18., 18.)).sel(lon=slice(-17.,60.))
Ep_AF_estim_ctl = Ep_estim_avg_ctl.where(landmask == 1.).sel(lat = slice(-18., 18.)).sel(lon=slice(-17.,60.))
epAFveg = epveg_avg.where(landmask == 1.).sel(lat = slice(-18., 18.)).sel(lon=slice(-17.,60.))
epAF = ep_avg.where(landmask == 1.).sel(lat = slice(-18., 18.)).sel(lon=slice(-17.,60.))
PAF = precipitation_avg.where(landmask == 1.).sel(lat = slice(-18., 18.)).sel(lon=slice(-17.,60.))
EAF = net_lhe_avg.where(landmask == 1.).sel(lat = slice(-18., 18.)).sel(lon=slice(-17.,60.))
PAFveg = precipitationveg_avg.where(landmask == 1.).sel(lat = slice(-18., 18.)).sel(lon=slice(-17.,60.))
EAFveg = net_lheveg_avg.where(landmask == 1.).sel(lat = slice(-18., 18.)).sel(lon=slice(-17.,60.))
fu_ctlAF = (epAF_ctl/PAF_ctl) + 1. - (1. + (epAF_ctl/PAF_ctl)**(omega_opt_AF))**(1./omega_opt_AF)
fu_ctlAF_estim = (Ep_AF_estim_ctl/PAF_ctl) + 1. - (1. + (Ep_AF_estim_ctl/PAF_ctl)**(omega_opt_AF_estim))**(1./omega_opt_AF_estim)

med = 18

fig, axes = plt.subplots(1,2, sharex = True, sharey = True, figsize = (25,10))

axes[0].plot((epAM_ctl/PAM_ctl), (EAM_ctl/PAM_ctl), 'b.', label = 'control')
axes[0].plot((Ep_AM_estim_ctl/PAM_ctl), (EAM_ctl/PAM_ctl), 'b*', label = 'control estim')
axes[0].plot((epAM/PAM), (EAM/PAM), 'r.', label = 'perturbed (SB)')
axes[0].plot((epAMveg/PAMveg), (EAMveg/PAMveg), 'g.', label = 'perturbed (VP05)')
axes[0].plot([1.0,1.0],[0.0,1.0],color='dimgray',linewidth=1, linestyle = 'dashed')
axes[0].plot([1.0,10.0],[1.0,1.0],color='dimgray',linewidth=1, linestyle = 'dashed')
axes[0].plot([0.0,1.0],[0.0,1.0],color='r',linewidth=1, linestyle = 'dashed')


axes[0].set_xlim(0.0, 10.0)
axes[0].set_ylim(0.0, 1.5)
axes[0].set_title('Amazon', fontsize = med)
axes[0].set_xlabel('E$_P$/P', fontsize = med)
axes[0].set_ylabel('E/P', fontsize = med)

axes[0].plot((epAM_ctl/PAM_ctl), fu_ctlAM, 'k.', label = 'Fu Eq.')
axes[0].plot((Ep_AM_estim_ctl/PAM_ctl), fu_ctlAM_estim, 'k*', label = 'Fu Eq. estim')

# axes[0].plot((ep/P), fu, 'k*')
# axes[0].plot((epveg/Pveg), fu_veg, 'kd')

axes[0].spines['right'].set_visible(False)
axes[0].spines['top'].set_visible(False)



axes[1].plot((epAF_ctl/PAF_ctl), (EAF_ctl/PAF_ctl), 'b.', label = 'control')
axes[1].plot((Ep_AF_estim_ctl/PAF_ctl), (EAF_ctl/PAF_ctl), 'b*', label = 'control estim')

axes[1].plot((epAF/PAF), (EAF/PAF), 'r.', label = 'perturbed (SB)')
axes[1].plot((epAFveg/PAFveg), (EAFveg/PAFveg), 'g.', label = 'perturbed (VP05)')
axes[1].plot([1.0,1.0],[0.0,1.0],color='dimgray',linewidth=1, linestyle = 'dashed')
axes[1].plot([1.0,10.0],[1.0,1.0],color='dimgray',linewidth=1, linestyle = 'dashed')
axes[1].plot([0.0,1.0],[0.0,1.0],color='r',linewidth=1, linestyle = 'dashed')

axes[1].set_xlim(0.0, 10.0)
axes[1].set_ylim(0.0, 1.5)
axes[1].set_title('Central Africa', fontsize = med)
axes[1].set_xlabel('E$_P$/P', fontsize = med)
axes[1].set_ylabel('E/P', fontsize = med)


axes[1].plot((epAF_ctl/PAF_ctl), fu_ctlAF, 'k.', label = 'Fu Eq.')
axes[1].plot((Ep_AF_estim_ctl/PAF_ctl), fu_ctlAF_estim, 'k*', label = 'Fu Eq. estim')

# axes[1].plot((epAF/PAF), fu, 'k*')
# axes[1].plot((epAFveg/PAFveg), fu_veg, 'kd')

axes[1].spines['right'].set_visible(False)
axes[1].spines['top'].set_visible(False)
axes[1].tick_params(labelsize = med)
axes[0].tick_params(labelsize = med)

fig.legend(fontsize = med, bbox_to_anchor=(0.8,0.4))

fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/'+dire+'/budyko_curve_fueq_ep_vs_ep_estim_'+str(runmin)+'-'+str(runmax)+'.pdf', bbox_inches = 'tight', dpi=400)
fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/'+dire+'/budyko_curve_fueq_ep_vs_ep_estim_'+str(runmin)+'-'+str(runmax)+'.png', bbox_inches = 'tight', dpi=400)



outdir = 'Isca/ISCA_HPC/'+dire

any_configuration_plot_noshift(outdir,runmin,runmax,-90.,90.,(ep_avg_ctl/precipitation_avg_ctl).where((epveg_avg/precipitationveg_avg)<=2.),area_array,'','limits_ctl','energy_limit',landmaskxr, minval = 0., maxval = 2.)
any_configuration_plot_noshift(outdir,runmin,runmax,-90.,90.,(epveg_avg/precipitationveg_avg).where((epveg_avg/precipitationveg_avg)<=2.),area_array,'','limits_veg','energy_limit',landmaskxr, minval = 0., maxval = 2.)
any_configuration_plot_noshift(outdir,runmin,runmax,-90.,90.,(ep_avg/precipitation_avg).where((ep_avg/precipitation_avg)<=2.),area_array,'','limits_pertb','energy_limit',landmaskxr, minval = 0., maxval = 2.)
any_configuration_plot_noshift(outdir,runmin,runmax,-90.,90.,(ep_avg_ctl),area_array,'mm/d','ep_ctl','fromwhite',landmaskxr, minval = 0., maxval = 50.)
any_configuration_plot_noshift(outdir,runmin,runmax,-90.,90., Ep_estim_avg_ctl,area_array,'mm/d','ep_estim_ctl','fromwhite',landmaskxr, minval = 0., maxval = 10.)


fig,axes = plt.subplots(1,2,sharex = True, sharey = True, figsize = (20,10))

axes[0].plot(precipitation_avg_ctl.where(landmask == 1.).sel(lat = slice(-12., 8.)).sel(lon=slice(-70.,-45.)), net_lhe_avg_ctl.where(landmask == 1.).sel(lat = slice(-12., 8.)).sel(lon=slice(-70.,-45.)), 'b.', label = 'SB ctl')
axes[0].plot(precipitation_avg.where(landmask == 1.).sel(lat = slice(-12., 8.)).sel(lon=slice(-70.,-45.)), net_lhe_avg.where(landmask == 1.).sel(lat = slice(-12., 8.)).sel(lon=slice(-70.,-45.)), 'r.', label = 'SB pert')
axes[0].plot(precipitationveg_avg.where(landmask == 1.).sel(lat = slice(-12., 8.)).sel(lon=slice(-70.,-45.)), net_lheveg_avg.where(landmask == 1.).sel(lat = slice(-12., 8.)).sel(lon=slice(-70.,-45.)), 'g.', label = 'VP pert')
# axes[0].plot(precipitation_avg_ctl.where(landmask == 1.).sel(lat = slice(-12., 8.)).sel(lon=slice(-70.,-45.)), ep_avg_ctl.where(landmask == 1.).sel(lat = slice(-12., 8.)).sel(lon=slice(-70.,-45.)), 'k.', label = 'E$_P$ ctl')
# axes[0].plot(precipitation_avg_ctl.where(landmask == 1.).sel(lat = slice(-12., 8.)).sel(lon=slice(-70.,-45.)), epveg_avg.where(landmask == 1.).sel(lat = slice(-12., 8.)).sel(lon=slice(-70.,-45.)), 'k*', label = 'E$_P$ veg')
axes[0].set_xlim(0.,15.)
axes[0].set_ylim(0.,6.)

axes[0].spines['right'].set_visible(False)
axes[0].spines['top'].set_visible(False)

axes[0].set_title('Amazon', fontsize = med)
axes[0].set_ylabel('E (mm/d)', fontsize = med)
axes[0].set_xlabel('P (mm/d)', fontsize = med)


axes[1].plot(precipitation_avg_ctl.where(landmask == 1.).sel(lat = slice(-18., 18.)).sel(lon=slice(-17.,60.)), net_lhe_avg_ctl.where(landmask == 1.).sel(lat = slice(-18., 18.)).sel(lon=slice(-17.,60.)), 'b.', label = 'SB ctl')
axes[1].plot(precipitation_avg.where(landmask == 1.).sel(lat = slice(-18., 18.)).sel(lon=slice(-17.,60.)), net_lhe_avg.where(landmask == 1.).sel(lat = slice(-18., 18.)).sel(lon=slice(-17.,60.)), 'r.', label = 'SB pert')
axes[1].plot(precipitationveg_avg.where(landmask == 1.).sel(lat = slice(-18., 18.)).sel(lon=slice(-17.,60.)), net_lheveg_avg.where(landmask == 1.).sel(lat = slice(-18., 18.)).sel(lon=slice(-17.,60.)), 'g.', label = 'VP pert')
# axes[1].plot(precipitation_avg_ctl.where(landmask == 1.).sel(lat = slice(-18., 18.)).sel(lon=slice(-17.,60.)), ep_avg_ctl.where(landmask == 1.).sel(lat = slice(-18., 18.)).sel(lon=slice(-17.,60.)), 'k.', label = 'E$_P$ ctl')
# axes[1].plot(precipitation_avg_ctl.where(landmask == 1.).sel(lat = slice(-18., 18.)).sel(lon=slice(-17.,60.)), epveg_avg.where(landmask == 1.).sel(lat = slice(-18., 18.)).sel(lon=slice(-17.,60.)), 'k*', label = 'E$_P$ veg')


axes[1].set_xlabel('P (mm/d)', fontsize = med)
axes[1].set_title('Central Africa', fontsize = med)


axes[1].spines['right'].set_visible(False)
axes[1].spines['top'].set_visible(False)
axes[0].tick_params(labelsize = med)
axes[1].tick_params(labelsize = med)

fig.legend(fontsize = med, bbox_to_anchor=(0.8,0.8))

fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/'+dire+'/budyko_curve_E_versus_P_noEP.pdf', bbox_inches = 'tight', dpi=400)
fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/'+dire+'/budyko_curve_E_versus_P_noEP.png', bbox_inches = 'tight', dpi=400)
