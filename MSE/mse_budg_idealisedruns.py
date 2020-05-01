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



def mse_budg(run):
    
    # Load data
    data = xr.open_mfdataset('/scratch/mp586/Isca_DATA/ISCA_HPC/withtv/'+run+'/run01*/atmos_monthly_interp.nc')
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
    
    # Evaluate moist static energy (cpT + Lq + gz) flux into the column by mean flow
    mse = mc.cp_air * data.temp + mc.L * data.sphum + mc.grav * data.height
    
    # Evaluate eddy fluxes
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
    
    eddies = (uT_dx_eddy + vT_dy_eddy + wT_dp_eddy + 
              uq_dx_eddy + vq_dy_eddy + wq_dp_eddy +
              uz_dx_eddy + vz_dy_eddy + wz_dp_eddy)
    
    u_dmsedx_mean = -1. * (data.ucomp * gr.ddx(mse))
    v_dmsedy_mean = -1. * (data.vcomp * gr.ddy(mse, vector=False))
    w_dmsedp_mean = -1. * (data.omega * gr.ddp(mse))
    
    # Evaluate rate of change of internal plus potential energy (cvT + Lq + gz)
    e = ((mc.cp_air - mc.rdgas) * data.temp + mc.L * data.sphum + mc.grav * data.height)
    dedt = gr.ddt(e)
    
    # Take column integrals of dedt, advective fluxes, and mse
    def column_int(var_in):
        var_int = var_in.sum('pfull')/mc.grav
        return var_int
    
    u_dmsedx_mean_ci = column_int(u_dmsedx_mean*dp)
    v_dmsedy_mean_ci = column_int(v_dmsedy_mean*dp)
    w_dmsedp_mean_ci = column_int(w_dmsedp_mean*dp)
    eddies = column_int(eddies*dp)
    dedt = column_int(dedt*dp)
    mse = column_int(mse*dp)
    
    # Approximate residual as eddy flux to close budget
    residual = Fnet - dedt + u_dmsedx_mean + v_dmsedy_mean + w_dmsedp_mean + eddies
    
    return mse, Fnet, u_dmsedx_mean, v_dmsedy_mean, w_dmsedp_mean, eddies, residual



def mse_plots(run, pentads=[32,38,44,50,56], land_mask=None, mse_levels=np.arange(3100.,3221.,20.), mse_cbar = np.arange(3100.,3221.,40.)):
    
    plot_dir = '/scratch/rg419/plots/asymmetry_paper/'
    mkdir = sh.mkdir.bake('-p')
    mkdir(plot_dir)
    
    data = xr.open_dataset('/disco/share/rg419/Data_moist/climatologies/' + run + '.nc')
    
    mse, Fnet, dedt, u_dmsedx, v_dmsedy, w_dmsedp, eddies, residual = mse_budg(run)
    
#    rcParams['figure.figsize'] = 12, 8
    rcParams['figure.figsize'] = 12, 7
    rcParams['font.size'] = 11
    fig, ((ax1, ax2, ax3, ax4, ax5), 
          (ax6, ax7, ax8, ax9, ax10),
          (ax11, ax12, ax13, ax14, ax15),
          (ax16, ax17, ax18, ax19, ax20),
          (ax21, ax22, ax23, ax24, ax25)) = plt.subplots(5, 5, sharex='col', sharey='row')
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12, ax13, ax14, ax15, ax16, ax17, ax18, ax19, ax20, ax21, ax22, ax23, ax24, ax25]
    
    #levels=np.arange(-160.,161.,20.)
    levels=np.arange(-210.,211.,30.)
    
    i=0
    for pentad in pentads:
        f1 = (mse/1000000.).sel(xofyear=pentad).plot.contourf(ax=axes[0+5*i],x='lon',y='lat', add_labels=False, extend='both', levels=mse_levels, add_colorbar=False, cmap='Reds')
        (data.precipitation.sel(xofyear=pentad)*86400.).plot.contour(ax=axes[0+5*i],x='lon',y='lat', add_labels=False, levels=np.arange(0.,16.,5), cmap='Blues')
        f2 = Fnet.sel(xofyear=pentad).plot.contourf(ax=axes[1+5*i],x='lon',y='lat', add_labels=False, extend='both', levels=levels, add_colorbar=False)
        u_dmsedx.sel(xofyear=pentad).plot.contourf(ax=axes[2+5*i],x='lon',y='lat', add_labels=False,  extend='both', levels=levels, add_colorbar=False)
        v_dmsedy.sel(xofyear=pentad).plot.contourf(ax=axes[3+5*i],x='lon',y='lat', add_labels=False,  extend='both', levels=levels, add_colorbar=False)
        w_dmsedp.sel(xofyear=pentad).plot.contourf(ax=axes[4+5*i],x='lon',y='lat', add_labels=False, extend='both', levels=levels, add_colorbar=False)
        #eddies.sel(xofyear=pentad).plot.contourf(ax=axes[4+5*i],x='lon',y='lat', add_labels=False, extend='both', levels=levels, add_colorbar=False)
        #residual.sel(xofyear=pentad).plot.contourf(ax=axes[4+5*i],x='lon',y='lat', add_labels=False, extend='both', levels=levels, add_colorbar=False)
        
        #for ax in axes[0+5*i: 5+5*i]:
        #    (data.precipitation.sel(xofyear=pentad)*86400.).plot.contour(ax=ax,x='lon',y='lat', add_labels=False, levels=[5.,10.], colors='k', alpha=0.4)
        i=i+1
    
    for ax in axes:
        ax.set_ylim(-15.,45.)
        ax.set_yticks(np.arange(-15.,46.,15.))
        ax.set_xlim(90.,245.)
        ax.set_xticks(np.arange(90.,226.,45.))
        ax.grid(True,linestyle=':')
    
    
    if not land_mask==None:
        land = xr.open_dataset(land_mask)
        for ax in axes:
            land.land_mask.plot.contour(x='lon', y='lat', ax=ax, levels=np.arange(-1.,2.,1.), add_labels=False, colors='k')
            land.zsurf.plot.contour(ax=ax, x='lon', y='lat', levels=np.arange(2000.,3001.,1000.), add_labels=False, colors='k', zorder=2)
    
    axes_mse = axes[0:25:5]
    axes_heating = [ax2, ax3, ax4, ax5, ax7, ax8, ax9, ax10, ax12, ax13, ax14, ax15, ax17, ax18, ax19, ax20, ax22, ax23, ax24, ax25]
    
    i=0
    labels = ['a','f','k','p','u']
    timing = ['Pentad 32', 'Pentad 38', 'Pentad 44', 'Pentad 50', 'Pentad 56']
    for ax in axes_mse:
        ax.set_ylabel('Latitude')
        ax.text(52,45, labels[i], weight='bold')
        ax.text(20,-10, timing[i], rotation=90, weight='bold')
        i=i+1
    
    i=0
    labels = ['b','c','d','e','g','h','i','j','l','m','n','o','q','r','s','t','v','w','x','y']
    for ax in axes_heating:
        ax.text(75,45, labels[i], weight='bold')
        i=i+1
    
    for ax in axes[20:25]:
        ax.set_xlabel('Longitude')
        

    ax1.set_title('$\mathbf{\overline{h}}$', fontsize=11, weight='bold')
    ax2.set_title('$\mathbf{F_{net}}$', fontsize=11, weight='bold')
    ax3.set_title('$\mathbf{-\overline{u}\partial \overline{h}/\partial \overline{x}}$', fontsize=11, weight='bold')
    ax4.set_title('$\mathbf{-\overline{v}\partial \overline{h}/\partial \overline{y}}$', fontsize=11, weight='bold')
    #ax5.set_title('Eddies', fontsize=11, weight='bold')
    #ax5.set_title('Residual', fontsize=11, weight='bold')
    ax5.set_title('$\mathbf{-\overline{\omega}\partial \overline{h}/\partial \overline{p}}$', fontsize=11, weight='bold')
    
    plt.subplots_adjust(left=0.085, right=0.985, top=0.96, bottom=0.03, hspace=0.25, wspace=0.2)
    cb1=fig.colorbar(f1, ax=axes_mse, use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.08, aspect=30, shrink=1., ticks=mse_cbar)
    cb1.set_label('MJ/kg')
    cb1=fig.colorbar(f2, ax=axes_heating, use_gridspec=True, orientation = 'horizontal',fraction=0.05, pad=0.08, aspect=60, shrink=0.5)
    cb1.set_label('Wm$^{-2}$')
    
    plt.savefig(plot_dir + 'mse_budg_new_' + run + '.pdf', format='pdf')
    plt.close()


#mse_plots('half_land_nowishe', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc', mse_levels=np.arange(3100.,3301.,20.))
#mse_plots('half_nh_land_tibet_ol8_nowishe', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_nh_land_tibet.nc', mse_levels=np.arange(3100.,3301.,20.), mse_cbar=np.arange(3100.,3301.,50.))

#mse_plots('half_shallow', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc')
#mse_plots('half_10_land', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_10_shallow.nc', mse_levels=np.arange(3100.,3301.,20.), mse_cbar=np.arange(3100.,3301.,50.))

#mse_plots('half_nh_land', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_nh_shallow.nc', mse_levels=np.arange(3100.,3301.,20.), mse_cbar=np.arange(3100.,3301.,50.))
#mse_plots('half_shallow_sn_2.00', pentads=[64,76,88,100,112], land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc')
#mse_plots('half_shallow_sn_2.00', pentads=[64,70,76,82,88], land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc')

mse_plots('half_land', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_shallow.nc')
#mse_plots('half_nh_land_tibet_ol8', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_nh_land_tibet.nc', mse_levels=np.arange(3100.,3301.,20.), mse_cbar=np.arange(3100.,3301.,50.))
mse_plots('half_nh_land_tibet_ol8', land_mask = '/scratch/rg419/Experiments/asym_aquaplanets/input/half_nh_land_tibet.nc', mse_levels=np.arange(3160.,3281.,20.), mse_cbar=np.arange(3160.,3281.,40.))

    