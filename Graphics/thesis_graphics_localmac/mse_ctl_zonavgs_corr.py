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

    lats = ds[0].lat
    variables = ['mse_tot','mse_latent','mse_height', 'mse_sens']
    colors = ['k','cyan','grey','grey']
    labels = ['$F_{MSE}$', '$F_{Lq}$', '$F_{gz}$/10', '$F_{c_pT}$/10']
    factor = [-1, -1, -0.1, -0.1]
    style = ['-','-','--',':']
    fig, axes = plt.subplots(1,2,figsize = (20,10))
    plt.rcParams['ytick.labelsize']=24
    plt.rcParams['xtick.labelsize']=24 
    ctl_runs = ['(a) AP', '(b) AM-same']
    ds_ctl_here = [APctl, AMsamectl]
    ds_Pctl_here = [APctl_P, AMsamectl_P]
    for i in range(len(ctl_runs)):
        for j in range(len(variables)):
            axes[i].plot(lats, area_weighted_avg(xr.DataArray(factor[j]*(ds_ctl_here[i][variables[j]]), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[i],'all_sfcs',minlon = 0., maxlon = 360., axis=1), color = colors[j], label = labels[j], linewidth = 2, linestyle = style[j])
#        ax2 = axes[i].twinx()
#        ax2.plot(lats, area_weighted_avg(xr.DataArray(86400.*ds_Pctl_here[i].precipitation, coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[i],'all_sfcs',minlon = 0., maxlon = 360., axis=1), 'blue', label = 'P', linewidth = 2)        
        axes[i].plot(lats, area_weighted_avg(xr.DataArray(-1.*(ds_ctl_here[i].mse_sens + ds_ctl_here[i].mse_height), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[i],'all_sfcs',minlon = 0., maxlon = 360., axis=1), label = '$F_{DSE}$', color = 'grey', linewidth = 2)
        axes[i].plot(lats, area_weighted_avg(xr.DataArray(ds_ctl_here[i].Fnet, coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[i],'all_sfcs',minlon = 0., maxlon = 360., axis=1), label = '$F_{net}$', color = 'magenta', linewidth = 2)
        axes[i].set_ylim(-250.,250)
#        ax2.set_ylim(-12.,12.)
#        ax2.set_yticklabels([])
        axes[i].set_xlabel('Latitude ($^{\circ}$N)', fontsize=28)      
        axes[i].set_title(ctl_runs[i], fontsize=28)
    axes[0].set_ylabel('Energy Input and Flux Divergence (W/m$^2$)', fontsize=28)
    axes[0].legend(fontsize = 24, loc = 'lower center', ncol = 3)
    axes[1].set_yticklabels([])    
    axes[0].axhline(y=0., color='dimgray', linewidth = 1)    
    axes[1].axhline(y=0., color='dimgray', linewidth = 1)    
#    ax2.set_ylabel('Precipitation (mm/d)', fontsize=24, color = 'b')
#    ax2.set_yticklabels(['','','','0','4','8','12'])      
    fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/zonavg_energy_terms_ctls_CORRECTED_noP.png', bbox_inches = 'tight')
    fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/zonavg_energy_terms_ctls_CORRECTED_noP.pdf', bbox_inches = 'tight')

    fig, axes = plt.subplots(1,2,figsize = (20,10))
    plt.rcParams['ytick.labelsize']=24
    plt.rcParams['xtick.labelsize']=24 
    ctl_runs = ['(a) AP', '(b) AM-same']
    for i in range(len(ctl_runs)):
        for j in range(len(variables)):
            axes[i].plot(lats, area_weighted_avg(xr.DataArray(factor[j]*(ds_ctl_here[i][variables[j]]), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[3],'land',minlon = 0., maxlon = 40., axis=1), color = colors[j], label = labels[j], linewidth = 2, linestyle = style[j])
#        ax2 = axes[i].twinx()
#        ax2.plot(lats, area_weighted_avg(xr.DataArray(86400.*ds_Pctl_here[i].precipitation, coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[3],'land',minlon = 0., maxlon = 40., axis=1), 'blue', label = 'P', linewidth = 2)        
        axes[i].plot(lats, area_weighted_avg(xr.DataArray(-1.*(ds_ctl_here[i].mse_sens + ds_ctl_here[i].mse_height), coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[3],'land',minlon = 0., maxlon = 40., axis=1), label = '$F_{DSE}$', color = 'grey', linewidth = 2)
        axes[i].plot(lats, area_weighted_avg(xr.DataArray(ds_ctl_here[i].Fnet, coords = [ds[0].lat,ds[0].lon], dims = ['lat','lon']),area_array,landmask[3],'land',minlon = 0., maxlon = 40., axis=1), label = '$F_{net}$', color = 'magenta', linewidth = 2)
        axes[i].set_ylim(-250.,250)
#        ax2.set_ylim(-12.,12.)
#        ax2.set_yticklabels([])
        axes[i].set_xlabel('Latitude ($^{\circ}$N)', fontsize=28)      
        axes[i].set_title(ctl_runs[i], fontsize=28)
    axes[0].set_ylabel('Energy Input and Flux Divergence (W/m$^2$)', fontsize=28)
    axes[0].legend(fontsize = 24, loc = 'lower center', ncol = 3)
    axes[0].axhline(y=0., color='dimgray', linewidth = 1)    
    axes[1].axhline(y=0., color='dimgray', linewidth = 1)    
    axes[1].set_yticklabels([])    
#    ax2.set_ylabel('Precipitation (mm/d)', fontsize=24, color = 'b')
#    ax2.set_yticklabels(['','','','0','4','8','12'])      
    fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/zonavg_energy_terms_land_ctls_CORRECTED_noP.png', bbox_inches = 'tight')
    fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/zonavg_energy_terms_land_ctls_CORRECTED_noP.pdf', bbox_inches = 'tight')



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
    axes[0,0].axhline(y=0., color='dimgray', linewidth = 1)    
    axes[1,0].axhline(y=0., color='dimgray', linewidth = 1) 
    axes[0,1].axhline(y=0., color='dimgray', linewidth = 1)    
    axes[1,1].axhline(y=0., color='dimgray', linewidth = 1)       
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
    axes[0,0].axhline(y=0., color='dimgray', linewidth = 1)    
    axes[1,0].axhline(y=0., color='dimgray', linewidth = 1) 
    axes[0,1].axhline(y=0., color='dimgray', linewidth = 1)    
    axes[1,1].axhline(y=0., color='dimgray', linewidth = 1)    
    plt.rcParams['ytick.labelsize']=22
    plt.rcParams['xtick.labelsize']=22
    fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/allexps_zonavg_energy_terms_land.pdf', bbox_inches = 'tight')
    fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/withtv/allexps_zonavg_energy_terms_land.png', bbox_inches = 'tight')





