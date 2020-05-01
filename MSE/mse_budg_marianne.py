'''19/6/2019 Evaluate vertically integrated MSE budget
'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '/scratch/mp586/Code/MSE')
import gradients as gr, model_constants as mc
from pylab import rcParams
import sh
sys.path.insert(0, '/scratch/mp586/Code/PYCODES')
from plotting_routines_kav7 import *
GFDL_BASE = os.environ['GFDL_BASE']
sys.path.insert(0, os.path.join(GFDL_BASE,'src/extra/python/scripts'))
import cell_area as ca

def mse_budg(run):
 
    # Load data
    data = xr.open_mfdataset('/scratch/mp586/Isca_DATA/ISCA_HPC/withtv/'+run+'/*/atmos_monthly_interp.nc')
    data = data.mean('time')

    data_dp = xr.open_dataset('/scratch/mp586/Isca_DATA/ISCA_HPC/withtv/'+run+'/run0121/atmos_monthly.nc')
    dp = data_dp.phalf.diff('phalf')*100
    dp = xr.DataArray(dp[::-1], coords = [data.pfull.values], dims = ['pfull'])

    # Evaluate upward radiative fluxes at the surface
    data['rflux_surf'] = data.t_surf ** 4. * mc.stefan - data.flux_sw - data.flux_lw
    # Evaluate downward radiative fluxes at the TOA
    data['rflux_toa'] = data.toa_sw - data.olr

    # Evaluate net flux of energy into column
    Fnet = data.flux_lhe + data.flux_t + data.rflux_surf + data.rflux_toa

    uT_eddy = data.ucomp_temp - data.ucomp * data.temp
    vT_eddy = data.vcomp_temp - data.vcomp * data.temp
    wT_eddy = data.omega_temp - data.omega * data.temp
    
    uT_dx_eddy = -1.* mc.cp_air * gr.ddx(uT_eddy)
    vT_dy_eddy = -1. * mc.cp_air * gr.ddy(vT_eddy)
    wT_dp_eddy = -1. * mc.cp_air * gr.ddp(wT_eddy)
    
    uq_eddy = data.sphum_u - data.ucomp * data.sphum
    vq_eddy = data.sphum_v - data.vcomp * data.sphum
    wq_eddy = data.sphum_w - data.omega * data.sphum
    
    uq_dx_eddy = -1. * mc.L * gr.ddx(uq_eddy)
    vq_dy_eddy = -1. * mc.L * gr.ddy(vq_eddy)
    wq_dp_eddy = -1. * mc.L * gr.ddp(wq_eddy)
    
    uz_eddy = data.ucomp_height - data.ucomp * data.height
    vz_eddy = data.vcomp_height - data.vcomp * data.height
    wz_eddy = data.omega_height - data.omega * data.height
    
    uz_dx_eddy = -1. * mc.grav * gr.ddx(uz_eddy)
    vz_dy_eddy = -1. * mc.grav * gr.ddy(vz_eddy)
    wz_dp_eddy = -1. * mc.grav * gr.ddp(wz_eddy)


    div_VT_mean  = data.ucomp * gr.ddx(mc.cp_air * data.temp) + data.vcomp * gr.ddy(mc.cp_air * data.temp, vector = False) + data.omega * gr.ddp(mc.cp_air * data.temp)
    div_Vq_mean  = data.ucomp * gr.ddx(mc.L * data.sphum) + data.vcomp * gr.ddy(mc.L * data.sphum, vector = False) + data.omega * gr.ddp(mc.L * data.sphum)
    div_Vz_mean  = data.ucomp * gr.ddx(mc.grav * data.height) + data.vcomp * gr.ddy(mc.grav * data.height, vector = False) + data.omega * gr.ddp(mc.grav * data.height)

    def column_int(var_in):
        var_int = var_in.sum('pfull')/mc.grav
        return var_int

    # def column_int_bds(var_in,lowbd,upbd):
    #     var_int = var_in.sel(pfull = slice(lowbd,upbd)).sum('pfull')/mc.grav
    #     return var_int
    
    mse_sens_p = (uT_dx_eddy + vT_dy_eddy + wT_dp_eddy - div_VT_mean)*dp
    mse_latent_p = (uq_dx_eddy + vq_dy_eddy + wq_dp_eddy - div_Vq_mean)*dp
    mse_height_p = (uz_dx_eddy + vz_dy_eddy + wz_dp_eddy - div_Vz_mean)*dp 

    # mse_latent_PBL = column_int_bds((uq_dx_eddy + vq_dy_eddy + wq_dp_eddy - div_Vq_mean)*dp, 1000, 750) 

    mse_sens = column_int(mse_sens_p)
    mse_latent = column_int(mse_latent_p)
    mse_height = column_int(mse_height_p)

    mse_tot = mse_sens + mse_latent + mse_height

    # checks 
    # mse = mc.cp_air * data.temp + mc.L * data.sphum + mc.grav * data.height
    # mse_divV = -1. * mse * (gr.ddx(data.ucomp) + gr.ddy(data.vcomp, vector = False) + gr.ddp(data.omega))
    # mse_divV_ci = column_int(mse_divV*dp)
    # divV should be = zero (incompressible fluid!). div_v itself is close to zero, even when column integrated, but because mse is
    # on the order of 10**6 it ends up being a massive term!
    # see discussion with Ruth on slack!

    # mse_all = mse_tot + mse_divV_ci

    # div_hv should in theory be the same as mse_tot, but it's not. adding mse_divV_ci to mse_tot almost equals div_hv.
    # div_hv = mc.cp_air*(gr.ddx(data.ucomp_temp) + gr.ddy(data.vcomp_temp, vector = False) + gr.ddp(data.omega_temp)) + mc.L*(gr.ddx(data.sphum_u) + gr.ddy(data.sphum_v, vector = False) + gr.ddp(data.sphum_w)) + mc.grav*(gr.ddx(data.ucomp_height) + gr.ddy(data.vcomp_height, vector = False) + gr.ddp(data.omega_height))
    # div_hv = - div_hv 
    # div_hv_ci = column_int(div_hv*dp)

    # any_configuration_plot('Isca/ISCA_HPC/withtv/',120,480,-90.,90.,Fnet,area_array,'W/m2','AP_Fnet','PE_scale',landmaskxrAP,minval = -100., maxval = 100.)
    # any_configuration_plot('Isca/ISCA_HPC/withtv/',120,480,-90.,90.,mse_tot,area_array,'W/m2','AP_mse_orig','PE_scale',landmaskxrAP,minval = -100., maxval = 100.)
    # any_configuration_plot('Isca/ISCA_HPC/withtv/',120,480,-90.,90.,Fnet + mse_tot,area_array,'W/m2','AP_Fnet_plus_mse_orig','PE_scale',landmaskxrAP,minval = -50., maxval = 50.)
    # any_configuration_plot('Isca/ISCA_HPC/withtv/',120,480,-90.,90.,mse_all,area_array,'W/m2','AP_mse_allterms','PE_scale',landmaskxrAP,minval = -100., maxval = 100.)
    # any_configuration_plot('Isca/ISCA_HPC/withtv/',120,480,-90.,90.,mse_all - div_hv_ci,area_array,'W/m2','AP_mse_all_minus_div_hv','PE_scale',landmaskxrAP,minval = -100., maxval = 100.)

    # with vs without vertical component: 
    # mse_withdp = mc.cp_air*(gr.ddx(data.ucomp_temp) + gr.ddy(data.vcomp_temp, vector = False) + gr.ddp(data.omega_temp)) + mc.L*(gr.ddx(data.sphum_u) + gr.ddy(data.sphum_v, vector = False) + gr.ddp(data.sphum_w)) + mc.grav*(gr.ddx(data.ucomp_height) + gr.ddy(data.vcomp_height, vector = False) + gr.ddp(data.omega_height))
    # mse_withdp = - mse_withdp

    # mse_nodp = mc.cp_air*(gr.ddx(data.ucomp_temp) + gr.ddy(data.vcomp_temp, vector = False)) + mc.L*(gr.ddx(data.sphum_u) + gr.ddy(data.sphum_v, vector = False)) + mc.grav*(gr.ddx(data.ucomp_height) + gr.ddy(data.vcomp_height, vector = False))
    # mse_nodp = - mse_nodp

    # (xr.DataArray(column_int(mse_withdp*dp)) - xr.DataArray(column_int(mse_nodp*dp))).plot() --> saved a screenshot



    return mse_tot, Fnet, mse_sens, mse_latent, mse_height, mse_sens_p, mse_height_p, mse_latent_p

if __name__ == '__main__':

    run_list = [
    'aquaplanet_frierson_insolation_0qflux_mld20_commitd15c267',
    'aquaplanet_frierson_insolation_0qflux_mld20_plus_2xCO2_spinup_361_commitd15c267',
    'square_South_America_frierson_insolation_lepref1_0qflux_samealbedo_to_01land_samehcp_landocean_commitd15c267',
    'square_South_America_frierson_insolation_newbucket_0qflux_samealbedo_samehcp_landocean_commitd15c267',
    'square_South_America_frierson_insolation_newbucket_0qflux_samealbedo_samehcp_landocean_plus_2xCO2_spinup_361_commitd15c267',
    'square_South_America_frierson_insolation_newbucket_0qflux_samealbedo_to_01land_samehcp_landocean_commitd15c267']

    ds = ['','','','','','']
    ds_precip = ['','','','','','']

    # for i in range(len(ds)):
    #     [mse_tot, Fnet, mse_sens, mse_latent, mse_height, mse_sens_p, mse_height_p, mse_latent_p] = mse_budg(run_list[i])
    #     ds[i] = xr.Dataset({'mse_tot': mse_tot, 'Fnet': Fnet, 'mse_sens': mse_sens, 'mse_latent':mse_latent, 'mse_height':mse_height, 'mse_sens_p': mse_sens_p, 'mse_latent_p':mse_latent_p, 'mse_height_p':mse_height_p},
    #             # coords={'lon': (['lon'], Fnet.lon), 'lat': (['lat'], Fnet.lat)},
    #             # attrs={'units':'W/m2', 'long_name':run_list[i]},
    #             coords={'pfull': (['pfull'], mse_sens_p.pfull), 'lon': (['lon'], Fnet.lon), 'lat': (['lat'], Fnet.lat)},
    #             attrs={'units':'W/m2', 'long_name':run_list[i]}
    #             )
    # for i in range(len(ds)):
    #     ds[i].to_netcdf('/scratch/mp586/Isca_DATA/ISCA_HPC/withtv/MSEresult_files/'+run_list[i]+'_MSE.nc')

    for i in range(len(run_list)):
        ds[i] = xr.open_dataset('/scratch/mp586/Isca_DATA/ISCA_HPC/withtv/MSEresult_files/'+run_list[i]+'_MSE.nc')

    for i in range(len(run_list)):
        ds_precip[i] = xr.open_mfdataset('/scratch/mp586/Isca_DATA/ISCA_HPC/withtv/'+run_list[i]+'/*/atmos_monthly_interp.nc').mean('time')

    APctl = ds[0]
    APco2 = ds[1]
    APpatch = ds[2]
    AMsamectl = ds[3]
    AMsameco2 = ds[4]
    AM01 = ds[5]

    APctl_P = ds_precip[0]
    APco2_P = ds_precip[1]
    APpatch_P = ds_precip[2]
    AMsamectl_P = ds_precip[3]
    AMsameco2_P = ds_precip[4]
    AM01_P = ds_precip[5]

    area_array, dx, dy = ca.cell_area_all(t_res=42,base_dir='/scratch/mp586/Isca/') # added _all because then dx and dy are also returned 
    area_array = xr.DataArray(area_array) # returned in units of m bzw m^2, because radius in cell_area.py is given in metres
    area_array_3D = np.expand_dims(area_array, axis=0)
    area_array_3D = np.repeat(area_array_3D, 40, axis = 0) # to make area_array 3D (pressure, lat, lon)
        
    landfile=Dataset(os.path.join(GFDL_BASE,'input/square_South_America/land.nc'),mode='r')
    landmask=landfile.variables['land_mask'][:]
    landlats=landfile.variables['lat'][:]
    landlons=landfile.variables['lon'][:]
    # for specified lats
    landmaskxrAM=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

    landfile=Dataset(os.path.join(GFDL_BASE,'input/aquaplanet/land.nc'),mode='r')
    landmask=landfile.variables['land_mask'][:]
    landmaskxrAP=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

    landmask = [landmaskxrAP,landmaskxrAP,landmaskxrAP,landmaskxrAM,landmaskxrAM,landmaskxrAM]

    variables = ['mse_tot','mse_sens','mse_latent','mse_height']
    for i in range(len(run_list)):
        # for j in range(len(variables)):
        #     any_configuration_plot('Isca/ISCA_HPC/withtv/'+run_list[i],120,480,-90.,90., -1.*ds[i][variables[j]],area_array,'W/m2',variables[j],'PE_scale',landmask[i],minval = -100., maxval = 100.)
        # any_configuration_plot('Isca/ISCA_HPC/withtv/'+run_list[i],120,480,-90.,90.,ds[i].Fnet + ds[i].mse_tot,area_array,'W/m2','FNet_minus_msetot','PE_scale',landmask[i],minval = -20., maxval = 20.)
        # any_configuration_plot('Isca/ISCA_HPC/withtv/'+run_list[i],120,480,-90.,90.,ds[i].Fnet,area_array,'W/m2','FNet','PE_scale',landmask[i],minval = -100., maxval = 100.)
        any_configuration_plot('Isca/ISCA_HPC/withtv/'+run_list[i],120,480,-90.,90.,-1.*(ds[i].mse_sens + ds[i].mse_height),area_array,'W/m2','mse_sens_plus_height','PE_scale',landmask[i],minval = -100., maxval = 100.)

    # for j in range(len(variables)):
    #     any_configuration_plot('Isca/ISCA_HPC/withtv/',120,480,-90.,90.,-1.*(APco2-APctl)[variables[j]],area_array,'W/m2','APco2_minus_APctl_'+variables[j],'PE_scale',landmask[0],minval = -100., maxval = 100.)
    #     any_configuration_plot('Isca/ISCA_HPC/withtv/',120,480,-90.,90.,-1.*(APpatch-APctl)[variables[j]],area_array,'W/m2','APpatch_minus_APctl_'+variables[j],'PE_scale',landmask[0],minval = -100., maxval = 100.)
    #     any_configuration_plot('Isca/ISCA_HPC/withtv/',120,480,-90.,90.,-1.*(AMsamectl-APctl)[variables[j]],area_array,'W/m2','AMsamectl_minus_APctl_'+variables[j],'PE_scale',landmask[3],minval = -100., maxval = 100.)
    #     any_configuration_plot('Isca/ISCA_HPC/withtv/',120,480,-90.,90.,-1.*(AM01-AMsamectl)[variables[j]],area_array,'W/m2','AM01_minus_AMsamectl_'+variables[j],'PE_scale',landmask[3],minval = -100., maxval = 100.)
    #     any_configuration_plot('Isca/ISCA_HPC/withtv/',120,480,-90.,90.,-1.*(AMsameco2-AMsamectl)[variables[j]],area_array,'W/m2','AMsameco2_minus_AMsamectl_'+variables[j],'PE_scale',landmask[3],minval = -100., maxval = 100.)

    # any_configuration_plot('Isca/ISCA_HPC/withtv/',120,480,-90.,90.,(APco2-APctl).Fnet,area_array,'W/m2','APco2_minus_APctl_Fnet','PE_scale',landmask[0],minval = -100., maxval = 100.)
    # any_configuration_plot('Isca/ISCA_HPC/withtv/',120,480,-90.,90.,(APpatch-APctl).Fnet,area_array,'W/m2','APpatch_minus_APctl_Fnet','PE_scale',landmask[0],minval = -100., maxval = 100.)
    # any_configuration_plot('Isca/ISCA_HPC/withtv/',120,480,-90.,90.,(AMsamectl-APctl).Fnet,area_array,'W/m2','AMsamectl_minus_APctl_Fnet','PE_scale',landmask[3],minval = -100., maxval = 100.)
    # any_configuration_plot('Isca/ISCA_HPC/withtv/',120,480,-90.,90.,(AM01-AMsamectl).Fnet,area_array,'W/m2','AM01_minus_AMsamectl_Fnet','PE_scale',landmask[3],minval = -100., maxval = 100.)
    # any_configuration_plot('Isca/ISCA_HPC/withtv/',120,480,-90.,90.,(AMsameco2-AMsamectl).Fnet,area_array,'W/m2','AMsameco2_minus_AMsamectl_Fnet','PE_scale',landmask[3],minval = -100., maxval = 100.)

    # any_configuration_plot('Isca/ISCA_HPC/withtv/',120,480,-90.,90.,(APco2-APctl).Fnet + (APco2-APctl).mse_tot,area_array,'W/m2','APco2_minus_APctl_Fnet_minus_mse_tot','PE_scale',landmask[0],minval = -20., maxval = 20.)
    # any_configuration_plot('Isca/ISCA_HPC/withtv/',120,480,-90.,90.,(APpatch-APctl).Fnet + (APpatch-APctl).mse_tot,area_array,'W/m2','APpatch_minus_APctl_Fnet_minus_mse_tot','PE_scale',landmask[0],minval = -20., maxval = 20.)
    # any_configuration_plot('Isca/ISCA_HPC/withtv/',120,480,-90.,90.,(AMsamectl-APctl).Fnet + (AMsamectl-APctl).mse_tot,area_array,'W/m2','AMsamectl_minus_APctl_Fnet_minus_mse_tot','PE_scale',landmask[3],minval = -20., maxval = 20.)
    # any_configuration_plot('Isca/ISCA_HPC/withtv/',120,480,-90.,90.,(AM01-AMsamectl).Fnet + (AM01-AMsamectl).mse_tot,area_array,'W/m2','AM01_minus_AMsamectl_Fnet_minus_mse_tot','PE_scale',landmask[3],minval = -20., maxval = 20.)
    # any_configuration_plot('Isca/ISCA_HPC/withtv/',120,480,-90.,90.,(AMsameco2-AMsamectl).Fnet + (AMsameco2-AMsamectl).mse_tot,area_array,'W/m2','AMsameco2_minus_AMsamectl_Fnet_minus_mse_tot','PE_scale',landmask[3],minval = -20., maxval = 20.)

    # any_configuration_plot('Isca/ISCA_HPC/withtv/',120,480,-90.,90.,-1.*((APco2-APctl).mse_sens + (APco2-APctl).mse_height),area_array,'W/m2','APco2_minus_APctl_mse_sens_plus_mse_height','PE_scale',landmask[0],minval = -100., maxval = 100.)
    # any_configuration_plot('Isca/ISCA_HPC/withtv/',120,480,-90.,90.,-1.*((APpatch-APctl).mse_sens + (APpatch-APctl).mse_height),area_array,'W/m2','APpatch_minus_APctl_mse_sens_plus_mse_height','PE_scale',landmask[0],minval = -100., maxval = 100.)
    # any_configuration_plot('Isca/ISCA_HPC/withtv/',120,480,-90.,90.,-1.*((AMsamectl-APctl).mse_sens + (AMsamectl-APctl).mse_height),area_array,'W/m2','AMsamectl_minus_APctl_mse_sens_plus_mse_height','PE_scale',landmask[3],minval = -100., maxval = 100.)
    # any_configuration_plot('Isca/ISCA_HPC/withtv/',120,480,-90.,90.,-1.*((AM01-AMsamectl).mse_sens + (AM01-AMsamectl).mse_height),area_array,'W/m2','AM01_minus_AMsamectl_mse_sens_plus_mse_height','PE_scale',landmask[3],minval = -100., maxval = 100.)
    # any_configuration_plot('Isca/ISCA_HPC/withtv/',120,480,-90.,90.,-1.*((AMsameco2-AMsamectl).mse_sens + (AMsameco2-AMsamectl).mse_height),area_array,'W/m2','AMsameco2_minus_AMsamectl_mse_sens_plus_mse_height','PE_scale',landmask[3],minval = -100., maxval = 100.)



    mses = ['mse_latent_p','mse_sens_p','mse_height_p']

    Y, Z = np.meshgrid(ds[0].lat, ds[0].pfull)
    v = np.linspace(-100.,100.,21)
    fig, axes = plt.subplots(5, 3, sharex = True, sharey = True, figsize = (30,15))

    for i in range(len(mses)):
        cset = axes[0,i].contourf(Y, Z, area_weighted_avg_4D(xr.DataArray(-1.*(APco2-APctl)[mses[i]], coords = [ds[0].pfull,ds[0].lat,ds[0].lon], dims = ['pfull','lat','lon']),area_array_3D,landmaskxrAP,'all_sfcs',minlon = 0., maxlon = 40., axis=2), v, cmap='bwr_r', extend = 'both')
        axes[1,i].contourf(Y, Z, area_weighted_avg_4D(xr.DataArray(-1.*(APpatch-APctl)[mses[i]], coords = [ds[0].pfull,ds[0].lat,ds[0].lon], dims = ['pfull','lat','lon']),area_array_3D,landmaskxrAP,'all_sfcs',minlon = 0., maxlon = 40., axis=2), v, cmap='bwr_r', extend = 'both')
        axes[2,i].contourf(Y, Z, area_weighted_avg_4D(xr.DataArray(-1.*(AMsamectl-APctl)[mses[i]], coords = [ds[0].pfull,ds[0].lat,ds[0].lon], dims = ['pfull','lat','lon']),area_array_3D,landmaskxrAP,'all_sfcs',minlon = 0., maxlon = 40., axis=2), v, cmap='bwr_r', extend = 'both')
        axes[3,i].contourf(Y, Z, area_weighted_avg_4D(xr.DataArray(-1.*(AM01-AMsamectl)[mses[i]], coords = [ds[0].pfull,ds[0].lat,ds[0].lon], dims = ['pfull','lat','lon']),area_array_3D,landmaskxrAP,'all_sfcs',minlon = 0., maxlon = 40., axis=2), v, cmap='bwr_r', extend = 'both')
        axes[4,i].contourf(Y, Z, area_weighted_avg_4D(xr.DataArray(-1.*(AMsameco2-AMsamectl)[mses[i]], coords = [ds[0].pfull,ds[0].lat,ds[0].lon], dims = ['pfull','lat','lon']),area_array_3D,landmaskxrAP,'all_sfcs',minlon = 0., maxlon = 40., axis=2), v, cmap='bwr_r', extend = 'both')

    axes[0,0].set_title('$\Delta$ $div (Lq*v)$', size = 22)
    axes[0,1].set_title('$\Delta$ $div (cT*v)$', size = 22)
    axes[0,2].set_title('$\Delta$ $div (gz*v)$', size = 22)
    axes[0,0].set_ylabel('$AP_{CO2}$', size = 18)
    axes[1,0].set_ylabel('$AP-patch$', size = 18)
    axes[2,0].set_ylabel('$AM-same$', size = 18)
    axes[3,0].set_ylabel('$AM01$', size = 18)
    axes[4,0].set_ylabel('$AM-same_{CO2}$', size = 18)
    axes[0,0].tick_params(labelsize = 18)
    axes[1,0].tick_params(labelsize = 18)
    axes[2,0].tick_params(labelsize = 18)
    axes[3,0].tick_params(labelsize = 18)
    axes[4,0].tick_params(labelsize = 18)
    axes[4,1].tick_params(labelsize = 18)
    axes[4,2].tick_params(labelsize = 18)

    fig.gca().invert_yaxis()
    cbar = fig.colorbar(cset,ax=axes) # ax = axes tells it to take space away from all the subplots. could adjust location by setting ax to axes[0,0] for example. 
    cbar.set_label('$W/m^2$', size = 22)
    cbar.ax.tick_params(labelsize=22)  
    fig.text(0.4, 0.04, 'Latitude ($^{\circ} N)$', ha='center', size = 22)
    fig.text(0.04, 0.5, 'Pressure (hPa)', va='center', rotation='vertical', size = 22)
    fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/mses_vertical_landlons.png', format = 'png', dpi = 400, bbox_inches='tight')
    fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/mses_vertical_landlons.pdf')

    lats = ds[0].lat
    variables = ['mse_tot','mse_latent','Fnet']
    labels = ['$F_{MSE}$', '$F_{Lq}$', '$F_{net}$']
    for i in range(len(run_list)):
        plt.close()
        fig, ax = plt.subplots(1,1,figsize = (25,10))
        for j in range(len(variables)):
            ax.plot(lats, area_weighted_avg(xr.DataArray(-1.*(ds[i][variables[j]]), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[i],'all_sfcs',minlon = 0., maxlon = 40., axis=1), label = labels[j], linewidth = 2)
        ax2 = ax.twinx()
        ax2.plot(lats, area_weighted_avg(xr.DataArray(86400.*ds_precip[i].precipitation, coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[i],'all_sfcs',minlon = 0., maxlon = 40., axis=1), 'b', label = 'P', linewidth = 2)        
        ax.plot(lats, area_weighted_avg(xr.DataArray(-1.*(ds[i].mse_sens + ds[i].mse_height), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[i],'all_sfcs',minlon = 0., maxlon = 40., axis=1), label = '$F_{DSE}$', linewidth = 2)
        ax.legend(fontsize = 20)
        ax.set_ylabel('Divergence ($W/m^2$)', fontsize=18)
        ax2.set_ylabel('Precipitation (mm/d)', fontsize=18, color = 'b')
        ax.set_ylim(-200.,200)
        ax2.set_ylim(0.,8.)
        ax.set_xlabel('Latitude ($^{\circ}$N)', fontsize=16)      
        fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/'+run_list[i]+'/zonavg_energy_terms_land.png')
        fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/'+run_list[i]+'/zonavg_energy_terms_land.pdf')

    dataset_diffs = [APpatch - APctl, AMsamectl - APctl, AMsameco2 - AMsamectl, AM01 - AMsamectl]
    dataset_diffs_names = ['APpatch_minus_APctl', 'AMsamectl_minus_APctl', 'AMsameco2_minus_AMsamectl', 'AM01_minus_AMsamectl']

    dataset_precip_diffs = [APpatch_P - APctl_P, AMsamectl_P - APctl_P, AMsameco2_P - AMsamectl_P, AM01_P - AMsamectl_P]

    for i in range(len(dataset_diffs)):
        plt.close()
        fig, ax = plt.subplots(1,1,figsize = (25,10))
        for j in range(len(variables)):
            ax.plot(lats, area_weighted_avg(xr.DataArray(-1.*(dataset_diffs[i][variables[j]]), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[i],'all_sfcs',minlon = 0., maxlon = 40., axis=1), label = labels[j], linewidth = 2)
        ax2 = ax.twinx()
        ax2.plot(lats, area_weighted_avg(xr.DataArray(86400.*dataset_precip_diffs[i].precipitation, coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[i],'all_sfcs',minlon = 0., maxlon = 40., axis=1), 'b', label = '$P$', linewidth = 2)        
        ax.plot(lats, area_weighted_avg(xr.DataArray(-1.*(dataset_diffs[i].mse_sens + dataset_diffs[i].mse_height), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[i],'all_sfcs',minlon = 0., maxlon = 40., axis=1), label = '$F_{DSE}$', linewidth = 2)
        ax.legend(fontsize = 20)
        ax.set_ylabel('Divergence ($W/m^2$)', fontsize=18)
        ax2.set_ylabel('Precipitation (mm/d)', fontsize=18, color = 'b')
        ax.set_ylim(-50.,50.)
        ax2.set_ylim(-4.,4.)
        ax.set_xlabel('Latitude ($^{\circ}$N)', fontsize=16)
        fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/'+dataset_diffs_names[i]+'_zonavg_energy_terms_land.png')
        fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/'+dataset_diffs_names[i]+'_zonavg_energy_terms_land.png')



    exit()
    
    variables = ['mse_tot','mse_sens','mse_latent','mse_height','Fnet']
    pref = [-1.,-1.,-1.,-1.,1.]
    small = 18 #largefonts 14 # smallfonts 10 # medfonts = 14
    med = 20 #largefonts 18 # smallfonts 14 # medfonts = 16
    lge = 22 #largefonts 22 # smallfonts 18 # medfonts = 20

    v = np.linspace(-100.,100.,21)

    # South America Only


    for j in range(len(variables)):

        plt.close()
        fig, axes = plt.subplots(2,2, figsize = (25,12))
        fig.subplots_adjust(hspace = 0.2, wspace = 0.05)      
        lats = ds[0].lat
        lons = ds[0].lon
        array = pref[j]*(APpatch - APctl)[variables[j]]
        axes[0,0].set_title('(a) AP-patch - AP', size = med)
        m = Basemap(projection='kav7',lon_0=0.,resolution='c', ax = axes[0,0])
        array = xr.DataArray(array,coords=[lats,lons],dims=['lat','lon'])

        array = np.asarray(array)
        array, lons_cyclic = addcyclic(array, lons)
        array,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))
        array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])

        lons = lons_cyclic
        m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0], fontsize=small)
        m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,0], fontsize=small)

        lon, lat = np.meshgrid(lons, lats)
        xi, yi = m(lon, lat)

        cs = m.contourf(xi,yi,array, v, cmap='bwr_r', extend = 'both')


        array = pref[j]*(AMsamectl - APctl)[variables[j]]
        landmask = np.asarray(landmaskxrAM)
        landlons = np.asarray(landmaskxrAM.lon)
        lats = ds[0].lat
        lons = ds[0].lon
        axes[0,1].set_title('(a) AM-same - AP', size = med)
        m = Basemap(projection='kav7',lon_0=0.,resolution='c', ax = axes[0,1])
        array = xr.DataArray(array,coords=[lats,lons],dims=['lat','lon'])

        array = np.asarray(array)
        array, lons_cyclic = addcyclic(array, lons)
        array,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))
        array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])

        lons = lons_cyclic
        m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0], fontsize=small)
        m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,0], fontsize=small)

        lon, lat = np.meshgrid(lons, lats)
        xi, yi = m(lon, lat)

        cs = m.contourf(xi,yi,array, v, cmap='bwr_r', extend = 'both')
        landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))
        landmask, lons_cyclic = addcyclic(landmask, landlons)
        m.contour(xi,yi,landmask, 1, color = 'k', linestyles = 'dashed')

        

        array = pref[j]*(AMsameco2 - AMsamectl)[variables[j]]
        landmask = np.asarray(landmaskxrAM)
        landlons = np.asarray(landmaskxrAM.lon)
        lats = ds[0].lat
        lons = ds[0].lon
        axes[1,0].set_title('(a) AM-same-CO2 - AM-same', size = med)
        m = Basemap(projection='kav7',lon_0=0.,resolution='c', ax = axes[1,0])
        array = xr.DataArray(array,coords=[lats,lons],dims=['lat','lon'])

        array = np.asarray(array)
        array, lons_cyclic = addcyclic(array, lons)
        array,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))
        array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])

        lons = lons_cyclic
        m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0], fontsize=small)
        m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,0], fontsize=small)

        lon, lat = np.meshgrid(lons, lats)
        xi, yi = m(lon, lat)

        cs = m.contourf(xi,yi,array, v, cmap='bwr_r', extend = 'both')
        landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))
        landmask, lons_cyclic = addcyclic(landmask, landlons)
        m.contour(xi,yi,landmask, 1, color = 'k', linestyles = 'solid')

        array = pref[j]*(AM01 - AMsamectl)[variables[j]]
        lats = ds[0].lat
        lons = ds[0].lon
        landmask = np.asarray(landmaskxrAM)
        landlons = np.asarray(landmaskxrAM.lon)
        axes[1,1].set_title('(a) AM-dark - AM-same', size = med)
        m = Basemap(projection='kav7',lon_0=0.,resolution='c', ax = axes[1,1])
        array = xr.DataArray(array,coords=[lats,lons],dims=['lat','lon'])

        array = np.asarray(array)
        array, lons_cyclic = addcyclic(array, lons)
        array,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))
        array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])

        lons = lons_cyclic
        m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0], fontsize=small)
        m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,0], fontsize=small)

        lon, lat = np.meshgrid(lons, lats)
        xi, yi = m(lon, lat)

        cs = m.contourf(xi,yi,array, v, cmap='bwr_r', extend = 'both')
        landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))
        landmask, lons_cyclic = addcyclic(landmask, landlons)
        m.contour(xi,yi,landmask, 1, color = 'k', linestyles = 'solid')

        cbar = fig.colorbar(cs, orientation = 'vertical', ax = axes, shrink = 0.65) 
        cbar.set_label('mm/d', size=med)
        cbar.ax.tick_params(labelsize=med)

        fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/MSE_4cases_'+variables[j]+'.png', bbox_inches = 'tight')
        fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/MSE_4cases_'+variables[j]+'.pdf', bbox_inches = 'tight', dpi = 400)


        