'''19/6/2019 Evaluate vertically integrated MSE budget
'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import sys
sys.path.insert(0, '/scratch/mp586/Code/Graphics/Graphics_for_chapter3/MSE')
import gradients as gr, model_constants as mc
from pylab import rcParams
import sh
from scipy.odr import *
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


    # Gross mosit stability calculation 

    # find level at which stream function has its max. : pm_tropis_lev 

    v_zonavg = data.vcomp.mean('lon')

    # divv_cumsum = np.cumsum(divv*dp/mc.grav, axis = 0)
    # v_cumsum = np.cumsum(v_zonavg*dp/mc.grav, axis = 0)

    c = 2*np.pi*6376.0e3*np.cos(v_zonavg.lat*np.pi/180) / mc.grav
    msf = xr.DataArray((np.cumsum(v_zonavg*dp, axis=v_zonavg.dims.index('pfull')))*c)
    pm_levs = np.argmax(np.abs(msf), axis = 0)
    pm_all_lats = data.pfull[pm_levs]

    # y, z = np.meshgrid(data.lat, data.pfull)
    # y,z
    # plt.close()
    # fig, ax = plt.subplots(1,1)
    # ax.contourf(y, z, msf, cmap='RdBu_r')
    # ax.plot(data.lat, pm_all_lats, 'wx')
    # fig.gca().invert_yaxis()

    area_array, dx, dy = ca.cell_area_all(t_res=42,base_dir='/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/') # added _all because then dx and dy are also returned 
    area_array = xr.DataArray(area_array, coords=[data.lat, data.lon], dims = ['lat','lon'])
    pm_tropics_lev = int(np.round(np.average((pm_levs).sel(lat=slice(-10.,10.)), weights = area_array.sel(lat=slice(-10.,10.))[:,0]))) # has to be a model level, not the exact avg

    mse_v = mc.cp_air * data.vcomp_temp + mc.L * data.sphum_v + mc.grav * data.vcomp_height
    mse = mc.cp_air * data.temp + mc.L * data.sphum + mc.grav * data.height
    divv = (gr.ddx(data.ucomp) + gr.ddy(data.vcomp))
    h_divv  =  mse * divv



    ### actually I think it there shouldnt be minus signs here because of the integration limits TOA to SFC, but I don't really understand why... but there aren't any in the streamfunction calculation which has the same TOA-SFC rather than SFC- TOA integration limits
    ### doesn't matter for calculations below because they are all fractional changes. 
    ll_div = column_int((divv*dp)[:pm_tropics_lev,:,:]) # from pm to surface in Neelin's definition
    hdivv = column_int((h_divv*dp)) # from toa to surface in Neelin's definition

    divv_output = data.div # is this 2D or 3D divergence???
    h_divv_output = mse * divv_output
    ll_div_output = column_int((divv_output*dp)[:pm_tropics_lev,:,:]) # from pm to surface in Neelin's definition
    hdivv_output = column_int((h_divv_output*dp)) # from toa to surface in Neelin's definition

    massflux = column_int((v_zonavg*dp)[pm_tropics_lev:,:]) # from TOA to pm in wei 2018 and Frierson 2007
    mseflux = column_int((mse_v.mean('lon')*dp)) # from toa to surface in wei and frierson 

    pm = data.pfull[pm_tropics_lev]

    return mse_tot, Fnet, mse_sens, mse_latent, mse_height, mse_sens_p, mse_height_p, mse_latent_p, ll_div, massflux, mseflux, hdivv, ll_div_output, hdivv_output, pm

if __name__ == '__main__':

    run_list = [
    'aquaplanet_frierson_insolation_0qflux_mld20_commitd15c267',
    'aquaplanet_frierson_insolation_0qflux_mld20_plus_2xCO2_spinup_361_commitd15c267',
    'square_South_America_frierson_insolation_lepref1_0qflux_samealbedo_to_01land_samehcp_landocean_commitd15c267',
    'square_South_America_frierson_insolation_newbucket_0qflux_samealbedo_samehcp_landocean_commitd15c267',
    'square_South_America_frierson_insolation_newbucket_0qflux_samealbedo_samehcp_landocean_plus_2xCO2_spinup_361_commitd15c267',
    'square_South_America_frierson_insolation_newbucket_0qflux_samealbedo_to_01land_samehcp_landocean_commitd15c267']

    ds = ['','','','','','']

    # for i in range(len(ds)):
    #     [mse_tot, Fnet, mse_sens, mse_latent, mse_height, mse_sens_p, mse_height_p, mse_latent_p, ll_div, massflux, mseflux, hdivv, ll_div_output, hdivv_output, pm] = mse_budg(run_list[i])
    #     ds[i] = xr.Dataset({'mse_tot': mse_tot, 'Fnet': Fnet, 'mse_sens': mse_sens, 'mse_latent':mse_latent, 'mse_height':mse_height, 'll_div':ll_div, 'h_divv':hdivv, 'll_div_op':ll_div_output, 'hdivv_op' : hdivv_output, 'pm': pm, 'massflux':massflux, 'mseflux':mseflux, 'mse_sens_p': mse_sens_p, 'mse_latent_p':mse_latent_p, 'mse_height_p':mse_height_p},
    #             # coords={'lon': (['lon'], Fnet.lon), 'lat': (['lat'], Fnet.lat)},
    #             # attrs={'units':'W/m2', 'long_name':run_list[i]},
    #             coords={'pfull': (['pfull'], mse_sens_p.pfull), 'lon': (['lon'], Fnet.lon), 'lat': (['lat'], Fnet.lat)},
    #             attrs={'units':'W/m2', 'long_name':run_list[i]}
    #             )
    # for i in range(len(ds)):
    #     ds[i].to_netcdf('/scratch/mp586/Isca_DATA/ISCA_HPC/withtv/MSEresult_files/'+run_list[i]+'_MSE.nc')

    for i in range(len(run_list)):
        ds[i] = xr.open_dataset('/scratch/mp586/Isca_DATA/ISCA_HPC/withtv/MSEresult_files/'+run_list[i]+'_MSE.nc')

    ds_precip = ['','','','','','']
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


    dataset_diffs = [APpatch - APctl, AMsamectl - APctl, AMsameco2 - AMsamectl, AM01 - AMsamectl]
    dataset_perts = [APpatch, AMsamectl, AMsameco2, AM01]
    dataset_pert_names = ['APpatch', 'AMsamectl', 'AMsameco2', 'AM01']

    dataset_diffs_names = ['AP-dark-patch_minus_AP-ctl', 'AM-same_minus_AP-ctl', 'AM-same-co2_minus_AM-same', 'AM-dark_minus_AM-same']
    dataset_ctls = [APctl, APctl, AMsamectl, AMsamectl]

    lon = APctl.lon
    lat = APctl.lat.sel(lat=slice(-10.,10.))
    X, Y = np.meshgrid(lon, lat)
    v = np.linspace(-2,2,21)

    frac_llconv_avgs = np.empty((4))
    frac_gms_avgs = np.empty((4))
    frac_fnet_avgs = np.empty((4))
    frac_llconv_sds = np.empty((4))
    frac_gms_sds = np.empty((4))
    frac_fnet_sds = np.empty((4))
    del_llconv_all = np.empty((4, 64, 128))
    del_gms_all  = np.empty((4, 64, 128))
    del_fnet_all  = np.empty((4, 64, 128))

    for i in range(len(dataset_diffs)):
        del_llconv = - (dataset_perts[i].ll_div - dataset_ctls[i].ll_div)
        del_llconv_all[i,:,:] = del_llconv
        # frac_llconv = del_llconv/(-dataset_ctls[i].ll_div)
        frac_llconv = (dataset_perts[i].ll_div/dataset_ctls[i].ll_div) - 1.
        frac_llconv_avg, frac_llconv_sd = area_weighted_avg(xr.DataArray(frac_llconv, coords = [ds[0].lat, ds[0].lon], dims = ['lat','lon']), area_array, landmask[3], 'land', minlat = -5., maxlat = 5., return_sd = True)

        # ############### inoue method
        # gms_ctl = (dataset_ctls[i].mse_tot/(dataset_ctls[i].mse_sens + dataset_ctls[i].mse_height))
        # gms_pert = (dataset_perts[i].mse_tot/(dataset_perts[i].mse_sens + dataset_perts[i].mse_height))
        # del_gms = gms_pert - gms_ctl
        # frac_gms = del_gms/(gms_ctl)
        # frac_gms_avg, frac_gms_sd = area_weighted_avg(xr.DataArray(frac_gms, coords = [ds[0].lat, ds[0].lon], dims = ['lat','lon']), area_array, landmask[3], 'land', minlat = -5., maxlat = 5., return_sd = True)
        # In [52]: frac_gms_sds
        # Out[52]: array([  0.17442484,   1.06834837,   1.2144091 ,  21.84435502])

        ######## neelin method
        gms_ctl = dataset_ctls[i].Fnet/(-dataset_ctls[i].ll_div)
        gms_pert = dataset_perts[i].Fnet/(-dataset_perts[i].ll_div)
        del_gms = gms_pert - gms_ctl
        del_gms_all[i,:,:] = del_gms
        frac_gms = gms_pert/gms_ctl - 1.
        frac_gms_avg, frac_gms_sd = area_weighted_avg(xr.DataArray(frac_gms, coords = [ds[0].lat, ds[0].lon], dims = ['lat','lon']), area_array, landmask[3], 'land', minlat = -5., maxlat = 5., return_sd = True)
        # In [9]: frac_gms_sds
        # Out[9]: array([  0.16292841,  26.16757147,  27.57816263,   0.65813981])
        # In [20]: frac_fnet_sds
        # Out[20]: array([ 0.04274859,  0.09930309,  0.0158936 ,  0.27438932])
        # In [23]: frac_llconv_sds
        # Out[23]: array([  0.45105979,   0.14635751,   3.59865019,  23.56154668])

        # gms_ctl = dataset_ctls[i].Fnet/(-dataset_ctls[i].ll_div_op)
        # gms_pert = dataset_perts[i].Fnet/(-dataset_perts[i].ll_div_op)
        # del_gms = gms_pert - gms_ctl
        # frac_gms = del_gms/(gms_ctl)
        # frac_gms_avg, frac_gms_sd = area_weighted_avg(xr.DataArray(frac_gms, coords = [ds[0].lat, ds[0].lon], dims = ['lat','lon']), area_array, landmask[3], 'land', minlat = -5., maxlat = 5., return_sd = True)


        # # ########## neelin method original
        # gms_ctl = (-dataset_ctls[i].h_divv)/(-dataset_ctls[i].ll_div)
        # gms_pert = (-dataset_perts[i].h_divv)/(-dataset_perts[i].ll_div)
        # del_gms = gms_pert - gms_ctl
        # frac_gms = gms_pert/gms_ctl - 1.
        # frac_gms_avg, frac_gms_sd = area_weighted_avg(xr.DataArray(frac_gms, coords = [ds[0].lat, ds[0].lon], dims = ['lat','lon']), area_array, landmask[3], 'land', minlat = -5., maxlat = 5., return_sd = True)
        # In [3]: frac_gms_sds
        # Out[3]: array([  27.7810138 ,  792.59325928,  130.07566745,    4.93360152])



        del_fnet = (dataset_perts[i].Fnet - dataset_ctls[i].Fnet)
        del_fnet_all[i,:,:] = del_fnet
        frac_fnet = dataset_perts[i].Fnet/(dataset_ctls[i].Fnet) - 1.
        frac_fnet_avg, frac_fnet_sd = area_weighted_avg(xr.DataArray(frac_fnet, coords = [ds[0].lat, ds[0].lon], dims = ['lat','lon']), area_array, landmask[3], 'land', minlat = -5., maxlat = 5., return_sd = True)


        frac_llconv_avgs[i] = frac_llconv_avg
        frac_gms_avgs[i] = frac_gms_avg
        frac_fnet_avgs[i] = frac_fnet_avg
        frac_llconv_sds[i] = frac_llconv_sd
        frac_gms_sds[i] = frac_gms_sd
        frac_fnet_sds[i] = frac_fnet_sd


               
        plt.close()

        fig, axes = plt.subplots(1, 3, figsize = (30,10))
        cs = axes[0].contourf(X, Y, frac_llconv.sel(lat=slice(-10.,10.)), v, cmap = 'BrBG', extend = 'both')
        axes[1].contourf(X, Y, frac_fnet.sel(lat=slice(-10.,10.)), v, cmap = 'BrBG', extend = 'both')
        axes[2].contourf(X, Y, frac_gms.sel(lat=slice(-10.,10.)), v, cmap = 'BrBG', extend = 'both')
        axes[0].set_title('frac_llconv', fontsize = 22)
        axes[1].set_title('frac_Fnet', fontsize = 22)
        axes[2].set_title('frac_gms', fontsize = 22)
        cbar = fig.colorbar(cs, orientation = 'horizontal', ax = axes, shrink = 0.65) 
        fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/frac_gms_neelin_changes_alltropics_'+dataset_diffs_names[i]+'.png', bbox_inches='tight')
        fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/frac_gms_neelin_changes_alltropics_'+dataset_diffs_names[i]+'.pdf', bbox_inches='tight')

        plt.close()

        fig, axes = plt.subplots(2,2, figsize = (20,20))
        cs = axes[0,0].contourf(X, Y, frac_llconv.sel(lat=slice(-10.,10.)), v, cmap = 'BrBG', extend = 'both')
        axes[0,1].contourf(X, Y, frac_fnet.sel(lat=slice(-10.,10.)), v, cmap = 'BrBG', extend = 'both')
        axes[1,0].contourf(X, Y, frac_gms.sel(lat=slice(-10.,10.)), v, cmap = 'BrBG', extend = 'both')
        axes[1,1].contourf(X, Y, (frac_llconv - frac_fnet + frac_gms).sel(lat=slice(-10.,10.)), v, cmap = 'BrBG', extend = 'both')
        axes[0,0].set_title('frac_llconv', fontsize = 22)
        axes[0,1].set_title('frac_Fnet', fontsize = 22)
        axes[1,0].set_title('frac_gms', fontsize = 22)
        axes[1,1].set_title('sum', fontsize = 22)
        cbar = fig.colorbar(cs, orientation = 'horizontal', ax = axes, shrink = 0.65) 
        fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/fracsum_gms_neelin_changes_alltropics_'+dataset_diffs_names[i]+'.png', bbox_inches='tight')
        fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/fracsum_gms_neelin_changes_alltropics_'+dataset_diffs_names[i]+'.pdf', bbox_inches='tight')

        plt.close()
        fig, axes = plt.subplots(1, 3, figsize = (30,10))
        cs = axes[0].contourf(X, Y, (-dataset_perts[i].ll_div).sel(lat=slice(-10.,10.)), cmap = 'BrBG', extend = 'both')
        axes[1].contourf(X, Y, (dataset_perts[i].Fnet).sel(lat=slice(-10.,10.)), cmap = 'BrBG', extend = 'both')
        axes[2].contourf(X, Y, (gms_pert).sel(lat=slice(-10.,10.)), cmap = 'BrBG', extend = 'both')
        axes[0].set_title('llconv', fontsize = 22)
        axes[1].set_title('Fnet', fontsize = 22)
        axes[2].set_title('gms', fontsize = 22)
        fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/gms_terms_neelin_alltropics_'+dataset_pert_names[i]+'.png', bbox_inches='tight')
        fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/gms_terms_neelin_alltropics_'+dataset_pert_names[i]+'.pdf', bbox_inches='tight')
        plt.close()
        fig, axes = plt.subplots(1, 3, figsize = (30,10))
        cs = axes[0].contourf(X, Y, (del_llconv).sel(lat=slice(-10.,10.)), cmap = 'BrBG', extend = 'both')
        axes[1].contourf(X, Y, (del_fnet).sel(lat=slice(-10.,10.)), cmap = 'BrBG', extend = 'both')
        axes[2].contourf(X, Y, (del_gms*10**-9).sel(lat=slice(-10.,10.)), cmap = 'BrBG', extend = 'both')
        axes[0].set_title('llconv', fontsize = 22)
        axes[1].set_title('Fnet', fontsize = 22)
        axes[2].set_title('gms', fontsize = 22)
        fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/delgms_terms_neelin_alltropics_AM01.png', bbox_inches='tight')
        fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/delgms_terms_neelin_alltropics_AM01.pdf', bbox_inches='tight')




    names = ['(a) AP-dark-patch minus AP', '(b) AM-same minus AP', '(c) AM-same-2xCO2 minus AM-same', '(d) AM-dark minus AM-same']
    fig, axes = plt.subplots(2, 2, sharex = True, sharey = True, figsize = (20,10))
    axes[0,0].bar([0.0],frac_llconv_avgs[0], 0.1, align = 'edge', color = ['darkcyan'], label = '- '+r'$\nabla . v_{h,2}$')
    axes[0,0].bar([0.12],frac_gms_avgs[0], 0.1, align = 'edge', color = ['lightcoral'], label = 'GMS')
    axes[0,0].bar([0.24],frac_fnet_avgs[0], 0.1, align = 'edge', color = ['slategray'], label = '$F_{net}$')
    axes[0,1].bar([0.0, 0.12, 0.24],[frac_llconv_avgs[1], frac_gms_avgs[1], frac_fnet_avgs[1]], 0.1, align = 'edge', color = ['darkcyan', 'lightcoral', 'slategray'])
    axes[1,0].bar([0.0, 0.12, 0.24],[frac_llconv_avgs[2], frac_gms_avgs[2], frac_fnet_avgs[2]], 0.1, align = 'edge', color = ['darkcyan', 'lightcoral', 'slategray'])
    axes[1,1].bar([0.0, 0.12, 0.24],[frac_llconv_avgs[3], frac_gms_avgs[3], frac_fnet_avgs[3]], 0.1, align = 'edge', color = ['darkcyan', 'lightcoral', 'slategray'])
    axes[0,0].set_title(names[0], fontsize = 22)
    axes[0,1].set_title(names[1], fontsize = 22)
    axes[1,0].set_title(names[2], fontsize = 22) 
    axes[1,1].set_title(names[3], fontsize = 22) 
    axes[0,0].set_ylim(-5.5,2.)
    plt.rcParams['ytick.labelsize']=22
    axes[0,0].set_xticklabels([])
    axes[0,1].set_xticklabels([])
    axes[1,0].set_xticklabels([])
    axes[1,1].set_xticklabels([])
    axes[0,0].plot(np.linspace(0,0.35,10), np.zeros((10)), 'k')
    axes[0,1].plot(np.linspace(0,0.35,10), np.zeros((10)), 'k')
    axes[1,0].plot(np.linspace(0,0.35,10), np.zeros((10)), 'k')
    axes[1,1].plot(np.linspace(0,0.35,10), np.zeros((10)), 'k')
    axes[0,0].spines['right'].set_visible(False)
    axes[0,0].spines['top'].set_visible(False)
    axes[1,0].spines['right'].set_visible(False)
    axes[1,0].spines['top'].set_visible(False)    
    axes[1,1].spines['right'].set_visible(False)
    axes[1,1].spines['top'].set_visible(False)    
    axes[0,1].spines['right'].set_visible(False)
    axes[0,1].spines['top'].set_visible(False)
    axes[0,0].legend(fontsize = 22, ncol = 3, mode = 'expand', loc = 8)

    fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/frac_gms_neelin_changes_land5s5n_bars_allexps.png', bbox_inches='tight')
    fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/frac_gms_neelin_changes_land5s5n_bars_allexps.pdf', bbox_inches='tight')





    del_llconv_all = xr.DataArray(del_llconv_all, coords = [np.linspace(0,3,4), ds[0].lat, ds[0].lon], dims = ['exp','lat','lon'])
    del_fnet_all = xr.DataArray(del_fnet_all, coords = [np.linspace(0,3,4), ds[0].lat, ds[0].lon], dims = ['exp','lat','lon'])
    del_gms_all = xr.DataArray(del_gms_all*10**-7, coords = [np.linspace(0,3,4), ds[0].lat, ds[0].lon], dims = ['exp','lat','lon'])
    names = ['(a) AP-dark-patch minus AP', '(b) AM-same minus AP', '(c) AM-same-2xCO$_2$ minus AM-same', '(d) AM-dark minus AM-same']
    colors = ['lightseagreen','tan','dodgerblue','cyan']

    small = 22
    med = 28
    lge = 34

    fig, axes = plt.subplots(2, 2, figsize = (30,20))
    plt.rcParams['ytick.labelsize']=med
    plt.rcParams['xtick.labelsize']=med
    # axes[0,0].scatter(del_fnet_all[0,:,:].sel(lat=slice(-10.,10.)).sel(lon=slice(0.,40.)), del_llconv_all[0,:,:].sel(lat=slice(-10.,10.)).sel(lon=slice(0.,40.)), c=colors[0])
    # axes[0,1].scatter(del_fnet_all[1,:,:].sel(lat=slice(-10.,10.)).sel(lon=slice(0.,40.)), del_llconv_all[1,:,:].sel(lat=slice(-10.,10.)).sel(lon=slice(0.,40.)), c=colors[1])
    # axes[1,0].scatter(del_fnet_all[2,:,:].sel(lat=slice(-10.,10.)).sel(lon=slice(0.,40.)), del_llconv_all[2,:,:].sel(lat=slice(-10.,10.)).sel(lon=slice(0.,40.)), c=colors[2])
    # axes[1,1].scatter(del_fnet_all[3,:,:].sel(lat=slice(-10.,10.)).sel(lon=slice(0.,40.)), del_llconv_all[3,:,:].sel(lat=slice(-10.,10.)).sel(lon=slice(0.,40.)), c=colors[3])
    axes[0,0].set_title(names[0], fontsize = med)
    axes[0,1].set_title(names[1], fontsize = med)
    axes[1,0].set_title(names[2], fontsize = med) 
    axes[1,1].set_title(names[3], fontsize = med)
    axes[0,0].spines['right'].set_visible(False)
    axes[0,0].spines['top'].set_visible(False)
    axes[1,0].spines['right'].set_visible(False)
    axes[1,0].spines['top'].set_visible(False)    
    axes[1,1].spines['right'].set_visible(False)
    axes[1,1].spines['top'].set_visible(False)    
    axes[0,1].spines['right'].set_visible(False)
    axes[0,1].spines['top'].set_visible(False)
    fig.text(0.05, 0.5, r'$\Delta (- \nabla . \vec{v}_{h,2}$) (kg/m$^2$s)', va='center', rotation='vertical', fontsize = lge)
    fig.text(0.45, 0.05, r'$\Delta F_{net}$ (W/m$^2$)', va='center', rotation='horizontal', fontsize = lge)
    for i in range(4):
        del_llconv_flat = np.asarray(del_llconv_all[i,:,:].sel(lat=slice(-10.,10.)).sel(lon=slice(0.,40.))).flatten()
        del_fnet_flat = np.asarray(del_fnet_all[i,:,:].sel(lat=slice(-10.,10.)).sel(lon=slice(0.,40.))).flatten()
        del_gms_flat = np.asarray(del_gms_all[i,:,:].sel(lat=slice(-10.,10.)).sel(lon=slice(0.,40.))).flatten()
        mask = ~np.isnan(del_fnet_flat)
        [slope, intercept, r_value, p_value, std_err] = stats.linregress(del_fnet_flat[mask],del_llconv_flat[mask])
        print(p_value)
        if i <= 1:
            axes[0,i].plot(del_fnet_flat, del_fnet_flat*slope + intercept, colors[i], label = 'r = '+"%.2f" % r_value) #+', p = '+"%.2f" % p_value)
            axes[0,i].scatter(del_fnet_flat, del_llconv_flat, c=colors[i])
            axes[0,i].legend(fontsize = med, loc = 4)
        else:
            axes[1,i-2].plot(del_fnet_flat, del_fnet_flat*slope + intercept, colors[i], label = 'r = '+"%.2f" % r_value) #+', p = '+"%.2f" % p_value)
            axes[1,i-2].legend(fontsize = med, loc = 4)
            axes[1,i-2].scatter(del_fnet_flat, del_llconv_flat, c=colors[i])

        # [slope, intercept, nan, nan, nan] = orthoregress(del_fnet_flat[mask],del_llconv_flat[mask])
        # if i <= 1:
        #     axes[0,i].plot(del_fnet_flat, del_fnet_flat*slope + intercept, 'k.')
        # else:
        #     axes[1,i-2].plot(del_fnet_flat, del_fnet_flat*slope + intercept, 'k.')
    fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/frac_gms_neelin_changes_llconv_fnet_correlations_allexps_land2.png', bbox_inches='tight')
    fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/frac_gms_neelin_changes_llconv_fnet_correlations_allexps_land2.pdf', bbox_inches='tight')    

#pvalue land
# 9.32935893403e-06
# 0.0119494214523
# 0.000438748613251
# 9.83208000867e-09

# pvalue all-lons
# 0.0
# 6.67229555769e-81
# 0.00251747840837
# 3.20811083283e-210

# pvalue ocean
# 2.5248880787e-41
# 3.81796291578e-236
# 0.10959520288
# 1.75495593986e-11


    fig, axes = plt.subplots(2, 2, figsize = (30,20))
    plt.rcParams['ytick.labelsize']=med
    plt.rcParams['xtick.labelsize']=med
    axes[0,0].set_title(names[0], fontsize = med)
    axes[0,1].set_title(names[1], fontsize = med)
    axes[1,0].set_title(names[2], fontsize = med) 
    axes[1,1].set_title(names[3], fontsize = med)
    axes[0,0].spines['right'].set_visible(False)
    axes[0,0].spines['top'].set_visible(False)
    axes[1,0].spines['right'].set_visible(False)
    axes[1,0].spines['top'].set_visible(False)    
    axes[1,1].spines['right'].set_visible(False)
    axes[1,1].spines['top'].set_visible(False)    
    axes[0,1].spines['right'].set_visible(False)
    axes[0,1].spines['top'].set_visible(False)
    fig.text(0.05, 0.5, r'$\Delta (- \nabla . \vec{v}_{h,2}$) (kg/m$^2$s)', va='center', rotation='vertical', fontsize = lge)
    fig.text(0.45, 0.05, r'$\Delta P$ (mm/d)', va='center', rotation='horizontal', fontsize = lge)
    for i in range(4):
        del_llconv_flat = np.asarray(del_llconv_all[i,:,:].sel(lat=slice(-10.,10.)).sel(lon=slice(0.,40.))).flatten()
        del_precip_flat = np.asarray((dataset_precip_diffs[i].precipitation*86400).sel(lat=slice(-10.,10.)).sel(lon=slice(0.,40.))).flatten()
        mask = ~np.isnan(del_precip_flat)
        [slope, intercept, r_value, p_value, std_err] = stats.linregress(del_precip_flat[mask],del_llconv_flat[mask])
        print(p_value)
        if i <= 1:
            axes[0,i].plot(del_precip_flat, del_precip_flat*slope + intercept, colors[i], label = 'r = '+"%.2f" % r_value) #+', p = '+"%.2f" % p_value)
            axes[0,i].scatter(del_precip_flat, del_llconv_flat, c=colors[i])
            axes[0,i].legend(fontsize = med, loc = 4)
        else:
            axes[1,i-2].plot(del_precip_flat, del_precip_flat*slope + intercept, colors[i], label = 'r = '+"%.2f" % r_value) #+', p = '+"%.2f" % p_value)
            axes[1,i-2].legend(fontsize = med, loc = 4)
            axes[1,i-2].scatter(del_precip_flat, del_llconv_flat, c=colors[i])

        # [slope, intercept, nan, nan, nan] = orthoregress(del_fnet_flat[mask],del_llconv_flat[mask])
        # if i <= 1:
        #     axes[0,i].plot(del_fnet_flat, del_fnet_flat*slope + intercept, 'k.')
        # else:
        #     axes[1,i-2].plot(del_fnet_flat, del_fnet_flat*slope + intercept, 'k.')
    fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/frac_gms_neelin_changes_llconv_precip_correlations_allexps_land2.png', bbox_inches='tight')
    fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/frac_gms_neelin_changes_llconv_precip_correlations_allexps_land2.pdf', bbox_inches='tight')    


#pval land
# 1.46239366616e-77
# 3.9531026458e-74
# 2.66205228556e-05
# 4.56860136595e-07

# pvals ocean
# 8.41859594802e-35
# 0.0
# 0.141925502703
# 1.13398685668e-05


    fig, axes = plt.subplots(2, 2, figsize = (30,20))
    plt.rcParams['ytick.labelsize']=med
    plt.rcParams['xtick.labelsize']=med
    axes[0,0].set_title(names[0], fontsize = med)
    axes[0,1].set_title(names[1], fontsize = med)
    axes[1,0].set_title(names[2], fontsize = med) 
    axes[1,1].set_title(names[3], fontsize = med)
    axes[0,0].spines['right'].set_visible(False)
    axes[0,0].spines['top'].set_visible(False)
    axes[1,0].spines['right'].set_visible(False)
    axes[1,0].spines['top'].set_visible(False)    
    axes[1,1].spines['right'].set_visible(False)
    axes[1,1].spines['top'].set_visible(False)    
    axes[0,1].spines['right'].set_visible(False)
    axes[0,1].spines['top'].set_visible(False)
    fig.text(0.05, 0.5, r'$\Delta (- \nabla . \vec{v}_{h,2}$) (kg/m$^2$s)', va='center', rotation='vertical', fontsize = lge)
    fig.text(0.45, 0.05, r'$\Delta (\delta h)$ (x10$^7$ J/kg)', va='center', rotation='horizontal', fontsize = lge)
    for i in range(4):
        del_llconv_flat = np.asarray(del_llconv_all[i,:,:].sel(lat=slice(-10.,10.)).sel(lon=slice(0.,40.)).where(del_gms_all[i,:,:].sel(lat=slice(-10.,10.)).sel(lon=slice(0.,40.)) > -0.9)).flatten() # use for no_outliers land .where(del_gms_all[i,:,:].sel(lat=slice(-10.,10.)).sel(lon=slice(0.,40.)) > -0.9), # use for zoom land .where(abs(del_gms_all[i,:,:].sel(lat=slice(-10.,10.)).sel(lon=slice(0.,40.))) < 0.01)
        del_fnet_flat = np.asarray(del_fnet_all[i,:,:].sel(lat=slice(-10.,10.)).sel(lon=slice(0.,40.)).where(del_gms_all[i,:,:].sel(lat=slice(-10.,10.)).sel(lon=slice(0.,40.)) > -0.9)).flatten() # use for zoom ocean .where(abs(del_gms_all[i,:,:].sel(lat=slice(-10.,10.)).sel(lon=slice(41.,359.))) < 10)
        del_gms_flat = np.asarray(del_gms_all[i,:,:].sel(lat=slice(-10.,10.)).sel(lon=slice(0.,40.)).where(del_gms_all[i,:,:].sel(lat=slice(-10.,10.)).sel(lon=slice(0.,40.)) > -0.9)).flatten() 
        mask = ~np.isnan(del_fnet_flat)
        [slope, intercept, r_value, p_value, std_err] = stats.linregress(del_gms_flat[mask],del_llconv_flat[mask])
        print(p_value)
        if i <= 1:
            axes[0,i].plot(del_gms_flat, del_gms_flat*slope + intercept, colors[i], label = 'r = '+"%.2f" % r_value) #+', p = '+"%.2f" % p_value)
            axes[0,i].legend(fontsize = med, loc = 4)
            axes[0,i].scatter(del_gms_flat, del_llconv_flat, c=colors[i])
        else:
            axes[1,i-2].plot(del_gms_flat, del_gms_flat*slope + intercept, colors[i], label = 'r = '+"%.2f" % r_value) #+', p = '+"%.2f" % p_value)
            axes[1,i-2].scatter(del_gms_flat, del_llconv_flat, c=colors[i])
            axes[1,i-2].legend(fontsize = med, loc = 4)


        # [slope, intercept, nan, nan, nan] = orthoregress(del_fnet_flat[mask],del_llconv_flat[mask])
        # if i <= 1:
        #     axes[0,i].plot(del_fnet_flat, del_fnet_flat*slope + intercept, 'k.')
        # else:
        #     axes[1,i-2].plot(del_fnet_flat, del_fnet_flat*slope + intercept, 'k.')
    axes[1,1].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/frac_gms_neelin_changes_llconv_gms_correlations_allexps_no_outlier_land.png', bbox_inches='tight')
    fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/frac_gms_neelin_changes_llconv_gms_correlations_allexps_no_outlier_land.pdf', bbox_inches='tight')    

# pvalue land no outlier
# 0.330285128104
# 0.429013947741
# 0.184366045682
# 0.83921659427

# pvalue land with outlier
# 0.50968282567
# 0.774249451374
# 0.184366045682
# 0.83921659427

# pvalue ocean
# 0.251726831412
# 0.747867981535
# 0.11561272222
# 0.31216053503

#pvalue ocean zoom
# 0.251726831412
# 0.0508639878679
# 0.818952045519
# 0.028853442929



    # frierson/wei method
    for i in range(len(dataset_perts)):
        gms_ctl = dataset_ctls[i].mseflux / dataset_ctls[i].massflux
        gms_pert = dataset_perts[i].mseflux / dataset_perts[i].massflux
        gms_frac = (gms_pert - gms_ctl)/gms_ctl
        v_frac = dataset_diffs[i].massflux / dataset_ctls[i].massflux
        hv_frac = dataset_diffs[i].mseflux / dataset_ctls[i].mseflux
        fig, ax = plt.subplots(1,1, figsize = (10,10))
        ax.plot(lat, gms_frac.sel(lat=slice(-10.,10.)), color = 'lightcoral', label = 'GMS')
        ax.plot(lat, v_frac.sel(lat=slice(-10.,10.)), color = 'darkcyan', label = 'v')
        ax.plot(lat, hv_frac.sel(lat=slice(-10.,10.)), color = 'slategray', label = 'hv')
        ax.legend(fontsize = 22)
        fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/frac_gms_wei_changes_alltropics_'+dataset_diffs_names[i]+'.png', bbox_inches='tight')
        fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/frac_gms_wei_changes_alltropics_'+dataset_diffs_names[i]+'.pdf', bbox_inches='tight')


        # fig, axes = plt.subplots(1, 3, figsize = (30,10))
        # cs = axes[0].contourf(X, Y, ((dataset_ctls[i].Fnet)*del_llconv).sel(lat=slice(-10.,10.)).sel(lon=slice(0.,40.)), v, cmap = 'BrBG', extend = 'both')
        # axes[1].contourf(X, Y, (gms_ctl*frac_fnet).sel(lat=slice(-10.,10.)).sel(lon=slice(0.,40.)), v, cmap = 'BrBG', extend = 'both')
        # axes[2].contourf(X, Y, (del_gms).sel(lat=slice(-10.,10.)).sel(lon=slice(0.,40.)), v, cmap = 'BrBG', extend = 'both')
        # axes[0].set_title('Fnet * del_divv')
        # axes[1].set_title('frac_Fnet * gms')
        # axes[2].set_title('del_gms')
        # cbar = fig.colorbar(cs, orientation = 'horizontal', ax = axes, shrink = 0.65) 
        # fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/gms_changes'+dataset_diffs_names[i]+'.png', bbox_inches='tight')



    # variables = ['ll_div','Fnet']
    # colors = ['k','cyan']
    # labels = ['$v_{h,2}$', '$F_{net}$']
    # names = ['(a) AP-dark-patch minus AP', '(b) AM-same minus AP', '(c) AM-same-2xCO$_2$ minus AM-same', '(d) AM-dark minus AM-same']
    # fig, axes = plt.subplots(2,2, figsize = (20,10))
    # for i in range(len(dataset_diffs)):
    #     for j in range(len(variables)):
    #         if i<=1:
    #             axes[0,i].plot(lats, area_weighted_avg(xr.DataArray((dataset_diffs[i][variables[j]])/APctl[variables[j]], coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[3],'land',minlon = 0., maxlon = 40., axis=1), label = labels[j], color = colors[j], linewidth = 2)
    #             axes[0,i].plot(lats, area_weighted_avg(xr.DataArray((dataset_diffs[i].mse_tot/(dataset_diffs[i].mse_sens + dataset_diffs[i].mse_height))/(APctl.mse_tot/(APctl.mse_sens + APctl.mse_height)), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[3],'land',minlon = 0., maxlon = 40., axis=1), label = '\Delta h', color = 'magenta', linewidth = 2)                
    #             axes[0,i].set_ylim(-100.,100.)
    #             axes[0,i].set_title(names[i], fontsize = 28)                
    #         else:
    #             axes[1,i-2].plot(lats, area_weighted_avg(xr.DataArray((dataset_diffs[i][variables[j]])/AMsamectl[variables[j]], coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[3],'land',minlon = 0., maxlon = 40., axis=1), color = colors[j], linewidth = 2)
    #             axes[1,i-2].plot(lats, area_weighted_avg(xr.DataArray((dataset_diffs[i].Fnet/(-dataset_diffs[i].ll_div))/(AMsamectl.Fnet/(-AMsamectl.ll_div)), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[3],'land',minlon = 0., maxlon = 40., axis=1), color = 'magenta', linewidth = 2)                
    #             axes[1,i-2].set_ylim(-100.,100.)
    #             axes[1,i-2].set_title(names[i], fontsize = 28)              
    #     fig.text(0.06, 0.5, '$\Delta$ Energy Input and Flux Divergence (W/m$^2$)', va='center', rotation='vertical', size = 24)
    #     # fig.text(0.94, 0.5, '$\Delta$ Precipitation (mm/d)', va='center', rotation='vertical', size = 24, color = 'b')
    #     fig.text(0.45, 0.07, 'Latitude ($^{\circ}$N)', va='center', rotation='horizontal', size = 24)
    # axes[0,0].legend(fontsize = 22, ncol = 3, mode = 'expand', loc = 8)
    # # axes[0,1].set_yticklabels([])
    # # axes[1,1].set_yticklabels([])
    # # axes[0,0].set_xticklabels([])
    # # axes[0,1].set_xticklabels([])
    # plt.rcParams['ytick.labelsize']=22
    # plt.rcParams['xtick.labelsize']=22
    # fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/allexps_zonavg_fracgms_terms_land.pdf', bbox_inches = 'tight')
    # fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/allexps_zonavg_fracgms_terms_land.png', bbox_inches = 'tight')



# (dataset_diffs[i].Fnet/(-dataset_diffs[i].ll_div))/(APctl.Fnet/(-APctl.ll_div))




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



    # mses = ['mse_latent_p','mse_sens_p','mse_height_p']

    # Y, Z = np.meshgrid(ds[0].lat, ds[0].pfull)
    # v = np.linspace(-100.,100.,21)
    # fig, axes = plt.subplots(5, 3, sharex = True, sharey = True, figsize = (30,15))

    # for i in range(len(mses)):
    #     cset = axes[0,i].contourf(Y, Z, area_weighted_avg_4D(xr.DataArray(-1.*(APco2-APctl)[mses[i]], coords = [ds[0].pfull,ds[0].lat,ds[0].lon], dims = ['pfull','lat','lon']),area_array_3D,landmaskxrAP,'all_sfcs',minlon = 0., maxlon = 40., axis=2), v, cmap='bwr_r', extend = 'both')
    #     axes[1,i].contourf(Y, Z, area_weighted_avg_4D(xr.DataArray(-1.*(APpatch-APctl)[mses[i]], coords = [ds[0].pfull,ds[0].lat,ds[0].lon], dims = ['pfull','lat','lon']),area_array_3D,landmaskxrAP,'all_sfcs',minlon = 0., maxlon = 40., axis=2), v, cmap='bwr_r', extend = 'both')
    #     axes[2,i].contourf(Y, Z, area_weighted_avg_4D(xr.DataArray(-1.*(AMsamectl-APctl)[mses[i]], coords = [ds[0].pfull,ds[0].lat,ds[0].lon], dims = ['pfull','lat','lon']),area_array_3D,landmaskxrAP,'all_sfcs',minlon = 0., maxlon = 40., axis=2), v, cmap='bwr_r', extend = 'both')
    #     axes[3,i].contourf(Y, Z, area_weighted_avg_4D(xr.DataArray(-1.*(AM01-AMsamectl)[mses[i]], coords = [ds[0].pfull,ds[0].lat,ds[0].lon], dims = ['pfull','lat','lon']),area_array_3D,landmaskxrAP,'all_sfcs',minlon = 0., maxlon = 40., axis=2), v, cmap='bwr_r', extend = 'both')
    #     axes[4,i].contourf(Y, Z, area_weighted_avg_4D(xr.DataArray(-1.*(AMsameco2-AMsamectl)[mses[i]], coords = [ds[0].pfull,ds[0].lat,ds[0].lon], dims = ['pfull','lat','lon']),area_array_3D,landmaskxrAP,'all_sfcs',minlon = 0., maxlon = 40., axis=2), v, cmap='bwr_r', extend = 'both')

    # axes[0,0].set_title('$\Delta$ $div (Lq*v)$', size = 22)
    # axes[0,1].set_title('$\Delta$ $div (cT*v)$', size = 22)
    # axes[0,2].set_title('$\Delta$ $div (gz*v)$', size = 22)
    # axes[0,0].set_ylabel('$AP_{CO2}$', size = 18)
    # axes[1,0].set_ylabel('$AP-patch$', size = 18)
    # axes[2,0].set_ylabel('$AM-same$', size = 18)
    # axes[3,0].set_ylabel('$AM01$', size = 18)
    # axes[4,0].set_ylabel('$AM-same_{CO2}$', size = 18)
    # axes[0,0].tick_params(labelsize = 18)
    # axes[1,0].tick_params(labelsize = 18)
    # axes[2,0].tick_params(labelsize = 18)
    # axes[3,0].tick_params(labelsize = 18)
    # axes[4,0].tick_params(labelsize = 18)
    # axes[4,1].tick_params(labelsize = 18)
    # axes[4,2].tick_params(labelsize = 18)

    # fig.gca().invert_yaxis()
    # cbar = fig.colorbar(cset,ax=axes) # ax = axes tells it to take space away from all the subplots. could adjust location by setting ax to axes[0,0] for example. 
    # cbar.set_label('$W/m^2$', size = 22)
    # cbar.ax.tick_params(labelsize=22)  
    # fig.text(0.4, 0.04, 'Latitude ($^{\circ} N)$', ha='center', size = 22)
    # fig.text(0.04, 0.5, 'Pressure (hPa)', va='center', rotation='vertical', size = 22)
    # fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/mses_vertical_landlons.png', format = 'png', dpi = 400, bbox_inches='tight')
    # fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/mses_vertical_landlons.pdf')


    run_list = ['aquaplanet_frierson_insolation_0qflux_mld20_commitd15c267']
    mses = ['mse_latent_p','mse_sens_p','mse_height_p']
    Y, Z = np.meshgrid(ds[0].lat, ds[0].pfull)
    v = np.linspace(-300.,300.,21)

    for j in range(len(run_list)):
        plt.close()
        fig, axes = plt.subplots(1, 3, sharex = True, sharey = True, figsize = (35,10))
        for i in range(len(mses)):
            cset = axes[i].contourf(Y, Z, area_weighted_avg_4D(xr.DataArray(-1.*(ds[j][mses[i]]), coords = [ds[0].pfull,ds[0].lat,ds[0].lon], dims = ['pfull','lat','lon']),area_array_3D,landmaskxrAP,'all_sfcs',minlon = 0., maxlon = 40., axis=2), v, cmap='bwr', extend = 'both')
        axes[0].set_title('(a) '+r'$\nabla . (Lq \cdot v)$', size = 28)
        axes[1].set_title('(b) '+r'$\nabla . (c_p T \cdot v)$', size = 28)
        axes[2].set_title('(c) '+r'$\nabla . (gz \cdot v)$', size = 28)
        plt.rcParams['ytick.labelsize']=22
        plt.rcParams['xtick.labelsize']=22

        fig.gca().invert_yaxis()
        cbar = fig.colorbar(cset,ax=axes) # ax = axes tells it to take space away from all the subplots. could adjust location by setting ax to axes[0,0] for example. 
        cbar.set_label('Flux Divergence (W/kg)', size = 24)
        cbar.ax.tick_params(labelsize=24)  
        fig.text(0.08, 0.5, 'Pressure (hPa)', va='center', rotation='vertical', fontsize = 28)
        fig.text(0.45, 0.05, 'Latitude ($^{\circ}$N)', va='center', rotation='horizontal', fontsize = 28)

        fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/'+run_list[j]+'/mses_vertical_landlons.png', format = 'png', dpi = 400, bbox_inches='tight')
        fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/'+run_list[j]+'/mses_vertical_landlons.pdf', bbox_inches='tight')


    run_list = ['aquaplanet_frierson_insolation_0qflux_mld20_commitd15c267']
    mses = ['mse_latent_p','mse_sens_p','mse_height_p']
    Y, Z = np.meshgrid(ds[0].lat, ds[0].pfull)
    v = np.linspace(-300.,300.,21)

    for j in range(len(run_list)):
        plt.close()
        fig, axes = plt.subplots(1, 3, sharex = True, sharey = True, figsize = (35,10))
        for i in range(len(mses)):
            cset = axes[i].contourf(Y, Z, area_weighted_avg_4D(xr.DataArray(-1.*(ds[j][mses[i]]), coords = [ds[0].pfull,ds[0].lat,ds[0].lon], dims = ['pfull','lat','lon']),area_array_3D,landmaskxrAP,'all_sfcs', axis=2), v, cmap='bwr', extend = 'both')
        axes[0].set_title('(a) '+r'$\nabla . (Lq \cdot v)$', size = 28)
        axes[1].set_title('(b) '+r'$\nabla . (c_p T \cdot v)$', size = 28)
        axes[2].set_title('(c) '+r'$\nabla . (gz \cdot v)$', size = 28)
        plt.rcParams['ytick.labelsize']=22
        plt.rcParams['xtick.labelsize']=22

        fig.gca().invert_yaxis()
        cbar = fig.colorbar(cset,ax=axes) # ax = axes tells it to take space away from all the subplots. could adjust location by setting ax to axes[0,0] for example. 
        cbar.set_label('Flux Divergence (W/kg)', size = 24)
        cbar.ax.tick_params(labelsize=24)  
        fig.text(0.08, 0.5, 'Pressure (hPa)', va='center', rotation='vertical', fontsize = 28)
        fig.text(0.45, 0.05, 'Latitude ($^{\circ}$N)', va='center', rotation='horizontal', fontsize = 28)

        fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/'+run_list[j]+'/mses_vertical.png', format = 'png', dpi = 400, bbox_inches='tight')
        fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/'+run_list[j]+'/mses_vertical.pdf', bbox_inches='tight')




    run_list = ['aquaplanet_frierson_insolation_0qflux_mld20_commitd15c267']
    mses = ['mse_latent_p','mse_sens_p','mse_height_p']
    Y, Z = np.meshgrid(ds[0].lat.sel(lat=slice(-30.,30.)), ds[0].pfull)
    v = np.linspace(-300.,300.,21)

    for j in range(len(run_list)):
        plt.close()
        fig, axes = plt.subplots(1, 3, sharex = True, sharey = True, figsize = (35,10))
        for i in range(len(mses)):
            zonavg = area_weighted_avg_4D(xr.DataArray(-1.*(ds[j][mses[i]]), coords = [ds[0].pfull,ds[0].lat,ds[0].lon], dims = ['pfull','lat','lon']),area_array_3D,landmaskxrAM,'land',minlon = 0., maxlon = 40., axis=2)
            zonavg = xr.DataArray(zonavg, coords = [ds[0].pfull,ds[0].lat], dims = ['pfull','lat'])
            cset = axes[i].contourf(Y, Z, zonavg.sel(lat=slice(-30.,30.)), v, cmap='bwr', extend = 'both')
        axes[0].set_title('(a) '+r'$\nabla . (Lq \cdot v)$', size = 28)
        axes[1].set_title('(b) '+r'$\nabla . (c_p T \cdot v)$', size = 28)
        axes[2].set_title('(c) '+r'$\nabla . (gz \cdot v)$', size = 28)
        plt.rcParams['ytick.labelsize']=22
        plt.rcParams['xtick.labelsize']=22

        fig.gca().invert_yaxis()
        cbar = fig.colorbar(cset,ax=axes) # ax = axes tells it to take space away from all the subplots. could adjust location by setting ax to axes[0,0] for example. 
        cbar.set_label('Flux Divergence (W/kg)', size = 24)
        cbar.ax.tick_params(labelsize=24)  
        fig.text(0.08, 0.5, 'Pressure (hPa)', va='center', rotation='vertical', fontsize = 28)
        fig.text(0.45, 0.05, 'Latitude ($^{\circ}$N)', va='center', rotation='horizontal', fontsize = 28)

        fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/'+run_list[j]+'/mses_vertical_land_only.png', format = 'png', dpi = 400, bbox_inches='tight')
        fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/'+run_list[j]+'/mses_vertical_land_only.pdf', bbox_inches='tight')


    lats = ds[0].lat
    variables = ['mse_tot','mse_latent','mse_height', 'mse_sens']
    colors = ['k','cyan','grey','grey']
    labels = ['$F_{MSE}$', '$F_{Lq}$', '$F_{gz}$/10', '$F_{c_pT}$/10']
    factor = [-1, -1, -0.1, -0.1]
    style = ['-','-','--',':']
    fig, axes = plt.subplots(1,2,figsize = (25,12))
    plt.rcParams['ytick.labelsize']=24
    plt.rcParams['xtick.labelsize']=24 
    ctl_runs = ['(a) AP', '(b) AM-same']
    for i in range(len(ctl_runs)):
        for j in range(len(variables)):
            axes[i].plot(lats, area_weighted_avg(xr.DataArray(factor[j]*(ds[i][variables[j]]), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[i],'all_sfcs',minlon = 0., maxlon = 360., axis=1), color = colors[j], label = labels[j], linewidth = 2, linestyle = style[j])
        ax2 = axes[i].twinx()
        ax2.plot(lats, area_weighted_avg(xr.DataArray(86400.*ds_precip[i].precipitation, coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[i],'all_sfcs',minlon = 0., maxlon = 360., axis=1), 'blue', label = 'P', linewidth = 2)        
        axes[i].plot(lats, area_weighted_avg(xr.DataArray(-1.*(ds[i].mse_sens + ds[i].mse_height), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[i],'all_sfcs',minlon = 0., maxlon = 360., axis=1), label = '$F_{DSE}$', color = 'grey', linewidth = 2)
        axes[i].plot(lats, area_weighted_avg(xr.DataArray(ds[i].Fnet, coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[i],'all_sfcs',minlon = 0., maxlon = 360., axis=1), label = '$F_{net}$', color = 'magenta', linewidth = 2)
        axes[i].set_ylim(-250.,250)
        ax2.set_ylim(-15.,15.)
        ax2.set_yticklabels([])
        axes[i].set_xlabel('Latitude ($^{\circ}$N)', fontsize=28)      
        axes[i].set_title(ctl_runs[i], fontsize=28)
    axes[0].set_ylabel('Energy Input and Flux Divergence (W/m$^2$)', fontsize=28)
    axes[0].legend(fontsize = 24, loc = 'lower center', ncol = 3)
    axes[1].set_yticklabels([])    
    ax2.set_ylabel('Precipitation (mm/d)', fontsize=24, color = 'b')
    ax2.set_yticklabels(['','','','0','5','10','15'])      
    fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/zonavg_energy_terms_ctls.png', bbox_inches = 'tight')
    fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/zonavg_energy_terms_ctls.pdf', bbox_inches = 'tight')

    fig, axes = plt.subplots(1,2,figsize = (20,10))
    plt.rcParams['ytick.labelsize']=24
    plt.rcParams['xtick.labelsize']=24 
    ctl_runs = ['(a) AP', '(b) AM-same']
    for i in range(len(ctl_runs)):
        for j in range(len(variables)):
            axes[i].plot(lats, area_weighted_avg(xr.DataArray(factor[j]*(ds[i][variables[j]]), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[3],'land',minlon = 0., maxlon = 40., axis=1), color = colors[j], label = labels[j], linewidth = 2, linestyle = style[j])
        ax2 = axes[i].twinx()
        ax2.plot(lats, area_weighted_avg(xr.DataArray(86400.*ds_precip[i].precipitation, coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[3],'land',minlon = 0., maxlon = 40., axis=1), 'blue', label = 'P', linewidth = 2)        
        axes[i].plot(lats, area_weighted_avg(xr.DataArray(-1.*(ds[i].mse_sens + ds[i].mse_height), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[3],'land',minlon = 0., maxlon = 40., axis=1), label = '$F_{DSE}$', color = 'grey', linewidth = 2)
        axes[i].plot(lats, area_weighted_avg(xr.DataArray(ds[i].Fnet, coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[3],'land',minlon = 0., maxlon = 40., axis=1), label = '$F_{net}$', color = 'magenta', linewidth = 2)
        axes[i].set_ylim(-250.,250)
        ax2.set_ylim(-15.,15.)
        ax2.set_yticklabels([])
        axes[i].set_xlabel('Latitude ($^{\circ}$N)', fontsize=28)      
        axes[i].set_title(ctl_runs[i], fontsize=28)
    axes[0].set_ylabel('Energy Input and Flux Divergence (W/m$^2$)', fontsize=28)
    axes[0].legend(fontsize = 24, loc = 'lower center', ncol = 3)
    axes[1].set_yticklabels([])    
    ax2.set_ylabel('Precipitation (mm/d)', fontsize=24, color = 'b')
    ax2.set_yticklabels(['','','','0','5','10','15'])      
    fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/zonavg_energy_terms_land_ctls.png', bbox_inches = 'tight')
    fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/zonavg_energy_terms_land_ctls.pdf', bbox_inches = 'tight')




    dataset_diffs = [APpatch - APctl, AMsamectl - APctl, AMsameco2 - AMsamectl, AM01 - AMsamectl]
    dataset_diffs_names = ['AP-dark-patch_minus_AP-ctl', 'AM-same_minus_AP-ctl', 'AM-same-co2_minus_AM-same', 'AM-dark_minus_AM-same']

    dataset_precip_diffs = [APpatch_P - APctl_P, AMsamectl_P - APctl_P, AMsameco2_P - AMsamectl_P, AM01_P - AMsamectl_P]
    variables = ['mse_tot','mse_latent','Fnet']

    # lats = ds[0].lat
    # variables = ['mse_tot','mse_latent','mse_height', 'mse_sens']
    # colors = ['k','cyan','grey','grey']
    # labels = ['$F_{MSE}$', '$F_{Lq}$', '$F_{gz}$/10', '$F_{c_pT}$/10']
    # factor = [-1, -1, -0.1, -0.1]
    # style = ['-','-','--',':']
    # for i in range(len(run_list)):
    #     plt.close()
    #     fig, ax = plt.subplots(1,1,figsize = (20,10))
    #     for j in range(len(variables)):
    #         ax.plot(lats, area_weighted_avg(xr.DataArray(factor[j]*(ds[i][variables[j]]), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[i],'all_sfcs',minlon = 0., maxlon = 360., axis=1), color = colors[j], label = labels[j], linewidth = 2, linestyle = style[j])
    #     ax2 = ax.twinx()
    #     ax2.plot(lats, area_weighted_avg(xr.DataArray(86400.*ds_precip[i].precipitation, coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[i],'all_sfcs',minlon = 0., maxlon = 360., axis=1), 'b', label = 'P', linewidth = 2)        
    #     ax.plot(lats, area_weighted_avg(xr.DataArray(-1.*(ds[i].mse_sens + ds[i].mse_height), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[i],'all_sfcs',minlon = 0., maxlon = 360., axis=1), label = '$F_{DSE}$', color = 'grey', linewidth = 2)
    #     ax.plot(lats, area_weighted_avg(xr.DataArray(ds[i].Fnet, coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[i],'all_sfcs',minlon = 0., maxlon = 360., axis=1), label = '$F_{net}$', color = 'magenta', linewidth = 2)
    #     ax.legend(fontsize = 24)
    #     ax.set_ylabel('Energy Input and Flux Divergence (W/m$^2$)', fontsize=22)
    #     ax2.set_ylabel('Precipitation (mm/d)', fontsize=22, color = 'b')
    #     ax.set_ylim(-250.,250)
    #     ax2.set_ylim(-15.,15.)
    #     ax2.set_yticklabels(['','','','0','5','10','15'])
    #     ax.set_xlabel('Latitude ($^{\circ}$N)', fontsize=22)      
    #     fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/'+run_list[i]+'/zonavg_energy_terms.png')
    #     fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/'+run_list[i]+'/zonavg_energy_terms.pdf')

    # for i in range(len(dataset_diffs)):
    #     plt.close()
    #     fig, ax = plt.subplots(1,1,figsize = (25,10))
    #     for j in range(len(variables)):
    #         ax.plot(lats, area_weighted_avg(xr.DataArray(-1.*(dataset_diffs[i][variables[j]]), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[i],'all_sfcs',minlon = 0., maxlon = 40., axis=1), label = labels[j], linewidth = 2)
    #     ax2 = ax.twinx()
    #     ax2.plot(lats, area_weighted_avg(xr.DataArray(86400.*dataset_precip_diffs[i].precipitation, coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[i],'all_sfcs',minlon = 0., maxlon = 40., axis=1), 'b', label = '$P$', linewidth = 2)        
    #     ax.plot(lats, area_weighted_avg(xr.DataArray(-1.*(dataset_diffs[i].mse_sens + dataset_diffs[i].mse_height), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[i],'all_sfcs',minlon = 0., maxlon = 40., axis=1), label = '$F_{DSE}$', linewidth = 2)
    #     ax.legend(fontsize = 20)
    #     ax.set_ylabel('Divergence ($W/m^2$)', fontsize=18)
    #     ax2.set_ylabel('Precipitation (mm/d)', fontsize=18, color = 'b')
    #     ax.set_ylim(-50.,50.)
    #     ax2.set_ylim(-4.,4.)
    #     ax.set_xlabel('Latitude ($^{\circ}$N)', fontsize=16)
    #     fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/'+dataset_diffs_names[i]+'_zonavg_energy_terms_land.png')
    #     fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/'+dataset_diffs_names[i]+'_zonavg_energy_terms_land.png')


    variables = ['mse_tot','mse_latent','Fnet']
    colors = ['k','cyan','magenta']
    labels = ['$F_{MSE}$', '$F_{Lq}$', '$F_{net}$']
    factor = [-1.,-1.,1.]
    names = ['(a) AP-dark-patch minus AP', '(b) AM-same minus AP', '(c) AM-same-2xCO$_2$ minus AM-same', '(d) AM-dark minus AM-same']
    fig, axes = plt.subplots(2,2, figsize = (20,10))
    for i in range(len(dataset_diffs)):
        for j in range(len(variables)):
            if i<=1:
                axes[0,i].plot(lats, area_weighted_avg(xr.DataArray(factor[j]*(dataset_diffs[i][variables[j]]), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[i],'all_sfcs',minlon = 0., maxlon = 360., axis=1), label = labels[j], color = colors[j], linewidth = 2)
                ax2 = axes[0,i].twinx()
                axes[0,i].set_ylim(-50.,50.)
                axes[0,i].set_title(names[i], fontsize = 28)                
                ax2.set_ylim(-4.,4.)
                if i == 0:
                    ax2.set_yticklabels([])
            else:
                axes[1,i-2].plot(lats, area_weighted_avg(xr.DataArray(factor[j]*(dataset_diffs[i][variables[j]]), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[i],'all_sfcs',minlon = 0., maxlon = 360., axis=1), color = colors[j], linewidth = 2)
                ax2 = axes[1,i-2].twinx()
                axes[1,i-2].set_ylim(-50.,50.)
                ax2.set_ylim(-4.,4.)
                axes[1,i-2].set_title(names[i], fontsize = 28)
                if i == 2:
                    ax2.set_yticklabels([])                
        ax2.plot(lats, area_weighted_avg(xr.DataArray(86400.*dataset_precip_diffs[i].precipitation, coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[i],'all_sfcs',minlon = 0., maxlon = 360., axis=1), 'b', label = '$P$', linewidth = 2)

        fig.text(0.06, 0.5, '$\Delta$ Energy Input and Flux Divergence (W/m$^2$)', va='center', rotation='vertical', size = 24)
        fig.text(0.94, 0.5, '$\Delta$ Precipitation (mm/d)', va='center', rotation='vertical', size = 24, color = 'b')
        fig.text(0.45, 0.07, 'Latitude ($^{\circ}$N)', va='center', rotation='horizontal', size = 24)
    axes[0,0].plot(lats, area_weighted_avg(xr.DataArray(-1.*(dataset_diffs[0].mse_sens + dataset_diffs[0].mse_height), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[0],'all_sfcs',minlon = 0., maxlon = 360., axis=1), label = '$F_{DSE}$', color='grey',linewidth = 2)
    axes[0,1].plot(lats, area_weighted_avg(xr.DataArray(-1.*(dataset_diffs[1].mse_sens + dataset_diffs[1].mse_height), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[1],'all_sfcs',minlon = 0., maxlon = 360., axis=1), color='grey',linewidth = 2)
    axes[1,0].plot(lats, area_weighted_avg(xr.DataArray(-1.*(dataset_diffs[2].mse_sens + dataset_diffs[2].mse_height), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[2],'all_sfcs',minlon = 0., maxlon = 360., axis=1), color='grey',linewidth = 2)
    axes[1,1].plot(lats, area_weighted_avg(xr.DataArray(-1.*(dataset_diffs[3].mse_sens + dataset_diffs[3].mse_height), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[3],'all_sfcs',minlon = 0., maxlon = 360., axis=1), color='grey',linewidth = 2)
    axes[0,0].legend(fontsize = 22, ncol = 4, mode = 'expand', loc = 8)
    axes[0,1].set_yticklabels([])
    axes[1,1].set_yticklabels([])
    axes[0,0].set_xticklabels([])
    axes[0,1].set_xticklabels([])
    plt.rcParams['ytick.labelsize']=22
    plt.rcParams['xtick.labelsize']=22
    fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/allexps_zonavg_energy_terms.pdf', bbox_inches = 'tight')
    fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/allexps_zonavg_energy_terms.png', bbox_inches = 'tight')


    variables = ['mse_tot','mse_latent','Fnet']
    colors = ['k','cyan','magenta']
    labels = ['$F_{MSE}$', '$F_{Lq}$/10', '$F_{net}$']
    factor = [-1.,-1/10.,1.]
    ylim = [100.,100.,10.,100.]
    ylimp = [8.,8.,3.,8.]
    names = ['(a) AP-dark-patch minus AP', '(b) AM-same minus AP', '(c) AM-same-2xCO$_2$ minus AM-same', '(d) AM-dark minus AM-same']
    fig, axes = plt.subplots(2,2, figsize = (25,12))
    for i in range(len(dataset_diffs)):
        for j in range(len(variables)):
            if i<=1:
                axes[0,i].plot(lats, area_weighted_avg(xr.DataArray(factor[j]*(dataset_diffs[i][variables[j]]), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[3],'land',minlon = 0., maxlon = 40., axis=1), label = labels[j], color = colors[j], linewidth = 2)
                ax2 = axes[0,i].twinx()
                axes[0,i].set_ylim(-ylim[i],ylim[i])
                axes[0,i].set_title(names[i], fontsize = 28)                
                ax2.set_ylim(-ylimp[i],ylimp[i])
                # if i == 0:
                #     ax2.set_yticklabels([])
            else:
                axes[1,i-2].plot(lats, area_weighted_avg(xr.DataArray(factor[j]*(dataset_diffs[i][variables[j]]), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[3],'land',minlon = 0., maxlon = 40., axis=1), color = colors[j], linewidth = 2)
                ax2 = axes[1,i-2].twinx()
                axes[1,i-2].set_ylim(-ylim[i],ylim[i])
                ax2.set_ylim(-ylimp[i],ylimp[i])
                axes[1,i-2].set_title(names[i], fontsize = 28)
                # if i == 2:
                #     ax2.set_yticklabels([])                
        ax2.plot(lats, area_weighted_avg(xr.DataArray(86400.*dataset_precip_diffs[i].precipitation, coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[3],'land',minlon = 0., maxlon = 40., axis=1), 'b', label = '$P$', linewidth = 2)

        fig.text(0.06, 0.5, '$\Delta$ Energy Input and Flux Divergence (W/m$^2$)', va='center', rotation='vertical', size = 24)
        fig.text(0.94, 0.5, '$\Delta$ Precipitation (mm/d)', va='center', rotation='vertical', size = 24, color = 'b')
        fig.text(0.45, 0.05, 'Latitude ($^{\circ}$N)', va='center', rotation='horizontal', size = 24)
    axes[0,0].plot(lats, area_weighted_avg(xr.DataArray(-1/10.*(dataset_diffs[0].mse_sens + dataset_diffs[0].mse_height), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[3],'land',minlon = 0., maxlon = 40., axis=1), label = '$F_{DSE}$/10', color='grey',linewidth = 2)
    axes[0,1].plot(lats, area_weighted_avg(xr.DataArray(-1/10.*(dataset_diffs[1].mse_sens + dataset_diffs[1].mse_height), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[3],'land',minlon = 0., maxlon = 40., axis=1), color='grey',linewidth = 2)
    axes[1,0].plot(lats, area_weighted_avg(xr.DataArray(-1/10.*(dataset_diffs[2].mse_sens + dataset_diffs[2].mse_height), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[3],'land',minlon = 0., maxlon = 40., axis=1), color='grey',linewidth = 2)
    axes[1,1].plot(lats, area_weighted_avg(xr.DataArray(-1/10.*(dataset_diffs[3].mse_sens + dataset_diffs[3].mse_height), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[3],'land',minlon = 0., maxlon = 40., axis=1), color='grey',linewidth = 2)
    axes[0,0].legend(fontsize = 22, ncol = 4, mode = 'expand', loc = 8)
    # axes[0,1].set_yticklabels([])
    # axes[1,1].set_yticklabels([])
    # axes[0,0].set_xticklabels([])
    # axes[0,1].set_xticklabels([])
    plt.rcParams['ytick.labelsize']=22
    plt.rcParams['xtick.labelsize']=22
    fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/allexps_zonavg_energy_terms_land.pdf', bbox_inches = 'tight')
    fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/allexps_zonavg_energy_terms_land.png', bbox_inches = 'tight')









    lats = ds[0].lat    
    variables = ['mse_tot'] #,'mse_sens','mse_latent','mse_height','Fnet']
    pref = [-1.,-1.,-1.,-1.,1.]
    small = 18 #largefonts 14 # smallfonts 10 # medfonts = 14
    med = 24 #largefonts 18 # smallfonts 14 # medfonts = 16
    lge = 28 #largefonts 22 # smallfonts 18 # medfonts = 20

    v = np.linspace(-20.,20.,17)

    # South America Only


    for j in range(len(variables)):

        plt.close()
        fig, axes = plt.subplots(2,2, figsize = (25,12))
        fig.subplots_adjust(hspace = 0.2, wspace = 0.05)      
        lats = ds[0].lat
        lons = ds[0].lon
        array = pref[j]*(APpatch - APctl)[variables[j]]
        axes[0,0].set_title('(a) AP-dark-patch minus AP', size = med)
        m = Basemap(projection='kav7',lon_0=0.,resolution='c', ax = axes[0,0])
        array = xr.DataArray(array,coords=[lats,lons],dims=['lat','lon'])
        landmask = np.asarray(landmaskxrAM)
        landlons = np.asarray(landmaskxrAM.lon)
        array = np.asarray(array)
        array, lons_cyclic = addcyclic(array, lons)
        array,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))
        array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])

        lons = lons_cyclic
        m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0], fontsize=small)
        m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,0], fontsize=small)

        lon, lat = np.meshgrid(lons, lats)
        xi, yi = m(lon, lat)

        cs = m.contourf(xi,yi,array, v, cmap='bwr', extend = 'both')
        landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))
        landmask, lons_cyclic = addcyclic(landmask, landlons)
        m.contour(xi,yi,landmask, 1, color = 'k', linestyles = 'dotted')

        array = pref[j]*(AMsamectl - APctl)[variables[j]]
        landmask = np.asarray(landmaskxrAM)
        landlons = np.asarray(landmaskxrAM.lon)
        lats = ds[0].lat
        lons = ds[0].lon
        axes[0,1].set_title('(b) AM-same minus AP', size = med)
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

        cs = m.contourf(xi,yi,array, v, cmap='bwr', extend = 'both')
        landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))
        landmask, lons_cyclic = addcyclic(landmask, landlons)
        m.contour(xi,yi,landmask, 1, color = 'k', linestyles = 'dashed')

        

        array = pref[j]*(AMsameco2 - AMsamectl)[variables[j]]
        landmask = np.asarray(landmaskxrAM)
        landlons = np.asarray(landmaskxrAM.lon)
        lats = ds[0].lat
        lons = ds[0].lon
        axes[1,0].set_title('(c) AM-same-2xCO2 minus AM-same', size = med)
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

        cs = m.contourf(xi,yi,array, v, cmap='bwr', extend = 'both')
        landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))
        landmask, lons_cyclic = addcyclic(landmask, landlons)
        m.contour(xi,yi,landmask, 1, color = 'k', linestyles = 'solid')

        array = pref[j]*(AM01 - AMsamectl)[variables[j]]
        lats = ds[0].lat
        lons = ds[0].lon
        landmask = np.asarray(landmaskxrAM)
        landlons = np.asarray(landmaskxrAM.lon)
        axes[1,1].set_title('(d) AM-dark minus AM-same', size = med)
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

        cs = m.contourf(xi,yi,array, v, cmap='bwr', extend = 'both')
        landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))
        landmask, lons_cyclic = addcyclic(landmask, landlons)
        m.contour(xi,yi,landmask, 1, color = 'k', linestyles = 'solid')

        cbar = fig.colorbar(cs, orientation = 'vertical', ax = axes, shrink = 0.65) 
        cbar.set_label('Energy Flux Divergence (W/m$^2$)', size=med)
        cbar.ax.tick_params(labelsize=med)

        fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/MSE_4cases_'+variables[j]+'.png', bbox_inches = 'tight')
        fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/MSE_4cases_'+variables[j]+'.pdf', bbox_inches = 'tight', dpi = 400)


        