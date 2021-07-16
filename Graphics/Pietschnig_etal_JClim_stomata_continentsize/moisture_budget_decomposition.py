### Parts of this Code have been adapted from Ruth Geen's code ###
### for calculating atmospheric energy and moisture fluxes ###

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

# landfile=Dataset(os.path.join(GFDL_BASE,'input/square_South_America/land.nc'),mode='r')                                                         
# landmask=landfile.variables['land_mask'][:]                             
# landlats=landfile.variables['lat'][:]                                   
# landlons=landfile.variables['lon'][:]                                                                                 
# landmaskxrAM=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it     

# data = xr.open_mfdataset('/scratch/mp586/Isca_DATA/ISCA_HPC/square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commitfe93b9d/*/atmos_monthly_interp.nc')
# data = data.mean('time')

# data_dp = xr.open_dataset('/scratch/mp586/Isca_DATA/ISCA_HPC/square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commitfe93b9d/run0121/atmos_monthly.nc')
# dp = data_dp.phalf.diff('phalf')*100
# dp = xr.DataArray(dp[::-1], coords = [data.pfull.values], dims = ['pfull'])

# PE = data.precipitation - (data.flux_lhe/28.)/86400.

# uq_eddy = data.sphum_u - data.ucomp * data.sphum
# vq_eddy = data.sphum_v - data.vcomp * data.sphum
# wq_eddy = data.sphum_w - data.omega * data.sphum

# uq_dx_eddy = -1. * gr.ddx(uq_eddy)
# vq_dy_eddy = -1. * gr.ddy(vq_eddy)
# wq_dp_eddy = -1. * gr.ddp(wq_eddy)

# div_Vq_mean = data.ucomp * gr.ddx(data.sphum) + data.vcomp * gr.ddy(data.sphum, vector = False) + data.omega * gr.ddp(data.sphum)

# def column_int(var_in):
#     var_int = var_in.sum('pfull')/mc.grav
#     return var_int

# # def column_int_bds(var_in,lowbd,upbd):
# #     var_int = var_in.sel(pfull = slice(lowbd,upbd)).sum('pfull')/mc.grav
# #     return var_int

# uq_p = (uq_dx_eddy + vq_dy_eddy + wq_dp_eddy - div_Vq_mean)*dp

# divqu_ci = column_int(uq_p)

# cons = PE-divqu_ci

# X, Y = np.meshgrid(data.lon, data.lat)
# Yflat = Y.flatten()
# plt.scatter(np.asarray(-divqu_ci.where(landmask == 0.)).flatten(),np.asarray(-PE.where(landmask==0.)).flatten(), c = Yflat, s=20., cmap='Blues_r')
# plt.ylim(-0.0001,0.0001)
# plt.xlim(-0.0001,0.0001)
# plt.scatter(np.asarray(-divqu_ci.where(landmask == 1.)).flatten(),np.asarray(-PE.where(landmask==1.)).flatten(), c = Yflat, s=20., cmap='Greens_r')
# plt.plot([-0.0001,0.0001],[-0.0001,0.0001],'k-')
# plt.xlabel('div.(qu)')
# plt.ylabel('E-P')

# a = plt.contourf(X, Y, cons, np.linspace(-0.00003,0.00003,61), cmap = 'BrBG', extend = 'both')
# plt.colorbar(a, extend = 'both')

# data_ctl = xr.open_mfdataset('/scratch/mp586/Isca_DATA/ISCA_HPC/square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_commitfe93b9d/*/atmos_monthly_interp.nc')
# data_ctl = data_ctl.mean('time')

# data_pert = xr.open_mfdataset('/scratch/mp586/Isca_DATA/ISCA_HPC/square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commitfe93b9d/*/atmos_monthly_interp.nc')
# data_pert = data_pert.mean('time')

# data_diff = data_pert-data_ctl


# dPE = data_diff.precipitation - (data_diff.flux_lhe/28.)/86400.

# dTH = data_ctl.ucomp * gr.ddx(data_diff.sphum) + data_ctl.vcomp * gr.ddy(data_diff.sphum, vector = False) + data_ctl.omega * gr.ddp(data_diff.sphum)
# dTH_ci = column_int(dTH*dp)

# dMCD = data_diff.ucomp * gr.ddx(data_ctl.sphum) + data_diff.vcomp * gr.ddy(data_ctl.sphum, vector = False) + data_diff.omega * gr.ddp(data_ctl.sphum)
# dMCD_ci = column_int(dMCD*dp)

# uq_eddy_ctl = data_ctl.sphum_u - data_ctl.ucomp * data_ctl.sphum
# vq_eddy_ctl = data_ctl.sphum_v - data_ctl.vcomp * data_ctl.sphum
# wq_eddy_ctl = data_ctl.sphum_w - data_ctl.omega * data_ctl.sphum

# uq_dx_eddy_ctl = gr.ddx(uq_eddy_ctl)
# vq_dy_eddy_ctl = gr.ddy(vq_eddy_ctl)
# wq_dp_eddy_ctl = gr.ddp(wq_eddy_ctl)


# uq_eddy_pert = data_pert.sphum_u - data_pert.ucomp * data_pert.sphum
# vq_eddy_pert = data_pert.sphum_v - data_pert.vcomp * data_pert.sphum
# wq_eddy_pert = data_pert.sphum_w - data_pert.omega * data_pert.sphum

# uq_dx_eddy_pert = gr.ddx(uq_eddy_pert)
# vq_dy_eddy_pert = gr.ddy(vq_eddy_pert)
# wq_dp_eddy_pert = gr.ddp(wq_eddy_pert)

# dTE = (uq_dx_eddy_pert - uq_dx_eddy_ctl) + (vq_dy_eddy_pert - vq_dy_eddy_ctl) + (wq_dp_eddy_pert - wq_dp_eddy_ctl)
# dTE_ci = column_int(dTE*dp)

# residual = dPE + dMCD_ci + dTH_ci + dTE_ci
# v = np.linspace(-0.00001,0.00001,41)
# fig, axes = plt.subplots(3,2,figsize=(20,30))
# X, Y = np.meshgrid(data_ctl.lon, data_ctl.lat)
# axes[0,0].contourf(X, Y, - residual, v, cmap = 'BrBG', extend = 'both')
# axes[1,0].contourf(X, Y, dPE, v, cmap = 'RdBu_r', extend = 'both')
# axes[1,1].contourf(X, Y, - dTH_ci, v, cmap = 'RdBu_r', extend = 'both')
# axes[2,0].contourf(X, Y, - dMCD_ci, v, cmap = 'RdBu_r', extend = 'both')
# cb = axes[2,1].contourf(X, Y, - dTE_ci, v, cmap = 'RdBu_r', extend = 'both')

# axes[0,0].set_title('residual')
# axes[1,0].set_title('dPE')
# axes[1,1].set_title('dTH')
# axes[2,1].set_title('dTE')
# axes[2,0].set_title('dMCD')
# fig.colorbar(cb)


### 2D ####

# dTH = data_ctl.ucomp * gr.ddx(data_diff.sphum) + data_ctl.vcomp * gr.ddy(data_diff.sphum, vector = False) + data_diff.sphum * (gr.ddx(data_ctl.ucomp) + gr.ddy(data_ctl.vcomp))
# dTH_ci = column_int(dTH*dp)

# dMCD = data_diff.ucomp * gr.ddx(data_ctl.sphum) + data_diff.vcomp * gr.ddy(data_ctl.sphum, vector = False) + data_ctl.sphum * (gr.ddx(data_diff.ucomp) + gr.ddy(data_diff.vcomp))
# dMCD_ci = column_int(dMCD*dp)

# uq_eddy_ctl = data_ctl.sphum_u - data_ctl.ucomp * data_ctl.sphum
# vq_eddy_ctl = data_ctl.sphum_v - data_ctl.vcomp * data_ctl.sphum

# uq_dx_eddy_ctl = gr.ddx(uq_eddy_ctl)
# vq_dy_eddy_ctl = gr.ddy(vq_eddy_ctl)

# uq_eddy_pert = data_pert.sphum_u - data_pert.ucomp * data_pert.sphum
# vq_eddy_pert = data_pert.sphum_v - data_pert.vcomp * data_pert.sphum

# uq_dx_eddy_pert = gr.ddx(uq_eddy_pert)
# vq_dy_eddy_pert = gr.ddy(vq_eddy_pert)

# dTE = (uq_dx_eddy_pert - uq_dx_eddy_ctl) + (vq_dy_eddy_pert - vq_dy_eddy_ctl)
# dTE_ci = column_int(dTE*dp)

# residual = - (dPE + dMCD_ci + dTH_ci + dTE_ci)
# v = np.linspace(-0.00001,0.00001,41)
# fig, axes = plt.subplots(3,2,figsize=(20,30))
# X, Y = np.meshgrid(data_ctl.lon, data_ctl.lat)
# axes[0,0].contourf(X, Y, -1.*residual, v, cmap = 'RdBu_r', extend = 'both')
# axes[1,0].contourf(X, Y, dPE, v, cmap = 'RdBu_r', extend = 'both')
# axes[1,1].contourf(X, Y, -1.*dTH_ci, v, cmap = 'RdBu_r', extend = 'both')
# axes[2,0].contourf(X, Y, -1.*dMCD_ci, v, cmap = 'RdBu_r', extend = 'both')
# cb = axes[2,1].contourf(X, Y, -1.*dTE_ci, v, cmap = 'RdBu_r', extend = 'both')

# axes[0,0].set_title('residual')
# axes[1,0].set_title('dPE')
# axes[1,1].set_title('dTH')
# axes[2,1].set_title('dTE')
# axes[2,0].set_title('dMCD')
# fig.colorbar(cb)


def moisture_budget(run_ctl, run_pert):

	def column_int(var_in):
		var_int = var_in.sum('pfull')/mc.grav
		return var_int

	data_ctl = xr.open_mfdataset('/scratch/mp586/Isca_DATA/ISCA_HPC/'+run_ctl+'/*/atmos_monthly_interp.nc')
	data_ctl = data_ctl.mean('time')

	data_pert = xr.open_mfdataset('/scratch/mp586/Isca_DATA/ISCA_HPC/'+run_pert+'/*/atmos_monthly_interp.nc')
	data_pert = data_pert.mean('time')

	data_diff = data_pert-data_ctl

	data_dp = xr.open_dataset('/scratch/mp586/Isca_DATA/ISCA_HPC/'+run_ctl+'/run0121/atmos_monthly.nc')
	dp = data_dp.phalf.diff('phalf')*100
	dp = xr.DataArray(dp[::-1], coords = [data_ctl.pfull.values], dims = ['pfull'])

	dPE = data_diff.precipitation - (data_diff.flux_lhe/28.)/86400.

	dTH = data_ctl.ucomp * gr.ddx(data_diff.sphum) + data_ctl.vcomp * gr.ddy(data_diff.sphum, vector = False) + data_ctl.omega * gr.ddp(data_diff.sphum)
	dTH_ci = column_int(dTH*dp)

	dMCD = data_diff.ucomp * gr.ddx(data_ctl.sphum) + data_diff.vcomp * gr.ddy(data_ctl.sphum, vector = False) + data_diff.omega * gr.ddp(data_ctl.sphum)
	dMCD_ci = column_int(dMCD*dp)

	uq_eddy_ctl = data_ctl.sphum_u - data_ctl.ucomp * data_ctl.sphum
	vq_eddy_ctl = data_ctl.sphum_v - data_ctl.vcomp * data_ctl.sphum
	wq_eddy_ctl = data_ctl.sphum_w - data_ctl.omega * data_ctl.sphum

	uq_dx_eddy_ctl = gr.ddx(uq_eddy_ctl)
	vq_dy_eddy_ctl = gr.ddy(vq_eddy_ctl)
	wq_dp_eddy_ctl = gr.ddp(wq_eddy_ctl)


	uq_eddy_pert = data_pert.sphum_u - data_pert.ucomp * data_pert.sphum
	vq_eddy_pert = data_pert.sphum_v - data_pert.vcomp * data_pert.sphum
	wq_eddy_pert = data_pert.sphum_w - data_pert.omega * data_pert.sphum

	uq_dx_eddy_pert = gr.ddx(uq_eddy_pert)
	vq_dy_eddy_pert = gr.ddy(vq_eddy_pert)
	wq_dp_eddy_pert = gr.ddp(wq_eddy_pert)

	dTE = (uq_dx_eddy_pert - uq_dx_eddy_ctl) + (vq_dy_eddy_pert - vq_dy_eddy_ctl) + (wq_dp_eddy_pert - wq_dp_eddy_ctl)
	dTE_ci = column_int(dTE*dp)

	residual = dPE + dMCD_ci + dTH_ci + dTE_ci

	return dPE*86400., dMCD_ci*86400., dTH_ci*86400., dTE_ci*86400., residual*86400. # convert from kg/m2s --> mm/d



if __name__ == '__main__':

	control_sb_dirs = [
	'narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387',
	'narrow_twelve_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387',
	'narrow_twentyfour_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387',
	'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387',
	'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387',
	'squareland_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387'
	]

	vp05_dirs = [
	'narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
	'narrow_twelve_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
	'narrow_twentyfour_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
	'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
	'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
	'squareland_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387'
	]

	simple_bucket_dirs = [
	'narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
	'narrow_twelve_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
	'narrow_twentyfour_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
	'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
	'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
	'squareland_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387'
	]

	ds = ['','','','','','']

	# for i in range(len(ds)):
	# 	[dPE, dMCD_ci, dTH_ci, dTE_ci, residual] = moisture_budget(control_sb_dirs[i], vp05_dirs[i])
	# 	ds[i] = xr.Dataset({'dPE':dPE, 'dMCD_ci':dMCD_ci, 'dTH_ci':dTH_ci,'dTE_ci':dTE_ci,'residual':residual},
	# 			coords={'lon': (['lon'], dPE.lon), 'lat': (['lat'], dPE.lat)},
	# 			attrs={'units':'mm/d', 'long_name':vp05_dirs[i]}
	# 			)

	# for i in range(len(ds)):
	# 	ds[i].to_netcdf('/scratch/mp586/Isca_DATA/ISCA_HPC/moisture_budget_decomp/'+vp05_dirs[i]+'_MBD.nc')

	for i in range(len(simple_bucket_dirs)):
		 ds[i] = xr.open_dataset('/scratch/mp586/Isca_DATA/ISCA_HPC/moisture_budget_decomp/'+simple_bucket_dirs[i]+'_MBD.nc')

	dsvp = ['','','','','','']
	for i in range(len(simple_bucket_dirs)):
		 dsvp[i] = xr.open_dataset('/scratch/mp586/Isca_DATA/ISCA_HPC/moisture_budget_decomp/'+vp05_dirs[i]+'_MBD.nc')	

	landnames = ['narrow_three','narrow_twelve','narrow_twentyfour','square_South_America','square_Africa','squareland']
	avg_dTH = np.empty((len(simple_bucket_dirs),1))
	sd_dTH = np.empty((len(simple_bucket_dirs),1))
	avg_dPE = np.empty((len(simple_bucket_dirs),1))
	sd_dPE = np.empty((len(simple_bucket_dirs),1))
	avg_dTE = np.empty((len(simple_bucket_dirs),1))
	sd_dTE = np.empty((len(simple_bucket_dirs),1))
	avg_dMCD = np.empty((len(simple_bucket_dirs),1))
	sd_dMCD = np.empty((len(simple_bucket_dirs),1))
	avg_res = np.empty((len(simple_bucket_dirs),1))
	sd_res = np.empty((len(simple_bucket_dirs),1))

	avgvp_dTH = np.empty((len(simple_bucket_dirs),1))
	sdvp_dTH = np.empty((len(simple_bucket_dirs),1))
	avgvp_dPE = np.empty((len(simple_bucket_dirs),1))
	sdvp_dPE = np.empty((len(simple_bucket_dirs),1))
	avgvp_dTE = np.empty((len(simple_bucket_dirs),1))
	sdvp_dTE = np.empty((len(simple_bucket_dirs),1))
	avgvp_dMCD = np.empty((len(simple_bucket_dirs),1))
	sdvp_dMCD = np.empty((len(simple_bucket_dirs),1))
	avgvp_res = np.empty((len(simple_bucket_dirs),1))
	sdvp_res = np.empty((len(simple_bucket_dirs),1))



	area_array, dx, dy = ca.cell_area_all(t_res=42,base_dir='/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/') # added _all because then dx and dy are also returned 
	area_array = xr.DataArray(area_array) # returned in units of m bzw m^2, because radius in cell_area.py is given in metres


	minlat = -10.
	maxlat = 10. 

	for i in range(len(simple_bucket_dirs)):
		landfile=Dataset(os.path.join(GFDL_BASE,'input/'+landnames[i]+'/land.nc'),mode='r')                                                         
		landmask=landfile.variables['land_mask'][:]                             
		landlats=landfile.variables['lat'][:]                                   
		landlons=landfile.variables['lon'][:]                                                                                 
		landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it     

		avg_dTH[i], sd_dTH[i] = area_weighted_avg(ds[i].dTH_ci,area_array,landmaskxr,'land', minlat = minlat, maxlat = maxlat,return_sd = True)
		avg_dTE[i], sd_dTE[i] = area_weighted_avg(ds[i].dTE_ci,area_array,landmaskxr,'land', minlat = minlat, maxlat = maxlat,return_sd = True)
		avg_dPE[i], sd_dPE[i] = area_weighted_avg(ds[i].dPE,area_array,landmaskxr,'land', minlat = minlat, maxlat = maxlat,return_sd = True)
		avg_dMCD[i], sd_dMCD[i] = area_weighted_avg(ds[i].dMCD_ci,area_array,landmaskxr,'land', minlat = minlat, maxlat = maxlat,return_sd = True)
		avg_res[i], sd_res[i] = area_weighted_avg(ds[i].residual,area_array,landmaskxr,'land', minlat = minlat, maxlat = maxlat,return_sd = True)

	fig, ax = plt.subplots(1,1, figsize = (10,10))
	ax.plot([3,12,24,40,60,100], avg_dPE, 'v', color = 'salmon', markersize = 10., label = '$\Delta$(P-E)')
	ax.plot([3,12,24,40,60,100], - avg_dTH, 's', color = 'darkred', markersize = 10., label = '$\Delta$TH')
	ax.plot([3,12,24,40,60,100], - avg_dMCD, 'd', color = 'cyan', markersize = 10., label = '$\Delta$MCD')
	ax.plot([3,12,24,40,60,100], - avg_dTE, 'o', color = 'slategrey', markersize = 10., label = '$\Delta$TE')
	ax.plot([3,12,24,40,60,100], - avg_res, 'p', color = 'darkcyan', markersize = 10., label = 'R')
	ax.plot([0,100],[0,0],'k')
	med = 24

	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.tick_params(labelsize = med)
	ax.tick_params(labelsize = med)
	ax.set_xlabel('Continental extent ($^{\circ}$ lon)', fontsize = med)
	ax.set_ylabel('Contribution per Area (mm/d/m$^2$)', fontsize = med)
	# ax.set_ylim([0.,8.])
	# ax.set_xlim(0.,110.)
	ax.set_xticks([3,12,24,40,60,100])
	ax.legend(fontsize = med, loc = 'upper right')
	ax.set_title('Moisture Budget Decomposition', fontsize = med)

	plt.savefig('/scratch/mp586/Code/Graphics/PE_decomp_cont_size_bucket_'+str(minlat)+'-'+str(maxlat)+'N.png', bbox_inches = 'tight', format = 'png', dpi=400)
	plt.savefig('/scratch/mp586/Code/Graphics/PE_decomp_cont_size_bucket_'+str(minlat)+'-'+str(maxlat)+'N.pdf', bbox_inches = 'tight', format = 'pdf', dpi=400)


	for i in range(len(vp05_dirs)):
		landfile=Dataset(os.path.join(GFDL_BASE,'input/'+landnames[i]+'/land.nc'),mode='r')                                                         
		landmask=landfile.variables['land_mask'][:]                             
		landlats=landfile.variables['lat'][:]                                   
		landlons=landfile.variables['lon'][:]                                                                                 
		landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it     

		avgvp_dTH[i], sdvp_dTH[i] = area_weighted_avg(dsvp[i].dTH_ci,area_array,landmaskxr,'land', minlat = minlat, maxlat = maxlat,return_sd = True)
		avgvp_dTE[i], sdvp_dTE[i] = area_weighted_avg(dsvp[i].dTE_ci,area_array,landmaskxr,'land', minlat = minlat, maxlat = maxlat,return_sd = True)
		avgvp_dPE[i], sdvp_dPE[i] = area_weighted_avg(dsvp[i].dPE,area_array,landmaskxr,'land', minlat = minlat, maxlat = maxlat,return_sd = True)
		avgvp_dMCD[i], sdvp_dMCD[i] = area_weighted_avg(dsvp[i].dMCD_ci,area_array,landmaskxr,'land', minlat = minlat, maxlat = maxlat,return_sd = True)
		avgvp_res[i], sdvp_res[i] = area_weighted_avg(dsvp[i].residual,area_array,landmaskxr,'land', minlat = minlat, maxlat = maxlat,return_sd = True)

	fig, ax = plt.subplots(1,1, figsize = (10,10))
	ax.plot([3,12,24,40,60,100], avgvp_dPE, 'v', color = 'salmon', markersize = 10., label = '$\Delta$(P-E)')
	ax.plot([3,12,24,40,60,100], - avgvp_dTH, 's', color = 'darkred', markersize = 10., label = '$\Delta$TH')
	ax.plot([3,12,24,40,60,100], - avgvp_dMCD, 'd', color = 'cyan', markersize = 10., label = '$\Delta$MCD')
	ax.plot([3,12,24,40,60,100], - avgvp_dTE, 'o', color = 'slategrey', markersize = 10., label = '$\Delta$TE')
	ax.plot([3,12,24,40,60,100], - avgvp_res, 'p', color = 'darkcyan', markersize = 10., label = 'R')
	ax.plot([0,100],[0,0],'k')
	med = 24

	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.tick_params(labelsize = med)
	ax.tick_params(labelsize = med)
	ax.set_xlabel('Continental extent ($^{\circ}$ lon)', fontsize = med)
	ax.set_ylabel('Contribution per Area (mm/d/m$^2$)', fontsize = med)
	# ax.set_ylim([0.,8.])
	# ax.set_xlim(0.,110.)
	ax.set_xticks([3,12,24,40,60,100])
	ax.legend(fontsize = med, loc = 'upper right')
	ax.set_title('Moisture Budget Decomposition', fontsize = med)

	plt.savefig('/scratch/mp586/Code/Graphics/PE_decomp_cont_size_vp05_'+str(minlat)+'-'+str(maxlat)+'N.png', bbox_inches = 'tight', format = 'png', dpi=400)
	plt.savefig('/scratch/mp586/Code/Graphics/PE_decomp_cont_size_vp05_'+str(minlat)+'-'+str(maxlat)+'N.pdf', bbox_inches = 'tight', format = 'pdf', dpi=400)

	fig, axes = plt.subplots(1,2, figsize = (20,10))

	axes[0].plot([3,12,24,40,60,100], avg_dPE, 'v', color = 'salmon', markersize = 10., label = '$\Delta$(P-E)')
	axes[0].plot([3,12,24,40,60,100], - avg_dTH, 's', color = 'darkred', markersize = 10., label = '$\Delta$TH')
	axes[0].plot([3,12,24,40,60,100], - avg_dMCD, 'd', color = 'cyan', markersize = 10., label = '$\Delta$MCD')
	axes[0].plot([3,12,24,40,60,100], - avg_dTE, 'o', color = 'slategrey', markersize = 10., label = '$\Delta$TE')
	axes[0].plot([3,12,24,40,60,100], - avg_res, 'p', color = 'darkcyan', markersize = 10., label = 'R')
	axes[0].plot([0,100],[0,0],'k')

	axes[1].plot([3,12,24,40,60,100], avgvp_dPE, 'v', color = 'salmon', markersize = 10., label = '$\Delta$(P-E)')
	axes[1].plot([3,12,24,40,60,100], - avgvp_dTH, 's', color = 'darkred', markersize = 10., label = '$\Delta$TH')
	axes[1].plot([3,12,24,40,60,100], - avgvp_dMCD, 'd', color = 'cyan', markersize = 10., label = '$\Delta$MCD')
	axes[1].plot([3,12,24,40,60,100], - avgvp_dTE, 'o', color = 'slategrey', markersize = 10., label = '$\Delta$TE')
	axes[1].plot([3,12,24,40,60,100], - avgvp_res, 'p', color = 'darkcyan', markersize = 10., label = 'R')
	axes[1].plot([0,100],[0,0],'k')
	med = 24

	axes[0].spines['right'].set_visible(False)
	axes[0].spines['top'].set_visible(False)
	axes[1].spines['right'].set_visible(False)
	axes[1].spines['top'].set_visible(False)
	axes[0].tick_params(labelsize = med)
	axes[0].tick_params(labelsize = med)
	axes[0].set_xlabel('Continental extent ($^{\circ}$ lon)', fontsize = med)
	axes[0].set_ylabel('Contribution per Area (mm/d/m$^2$)', fontsize = med)
	axes[1].tick_params(labelsize = med)
	axes[1].tick_params(labelsize = med)
	axes[0].set_xlabel('Continental extent ($^{\circ}$ lon)', fontsize = med)
	axes[1].set_xlabel('Continental extent ($^{\circ}$ lon)', fontsize = med)

	axes[0].set_ylim([-0.4,1.])
	axes[0].set_xlim(0.,110.)
	axes[1].set_ylim([-.4,1.])
	axes[1].set_xlim(0.,110.)
	axes[0].set_xticks([3,12,24,40,60,100])
	axes[1].legend(fontsize = med, loc = 'upper right')
	axes[0].set_title('(a) bucket', fontsize = med)
	axes[1].set_title('(b) 50%cond', fontsize = med)

	plt.savefig('/scratch/mp586/Code/Graphics/PE_decomp_cont_size_bucket_and_vp05_'+str(minlat)+'-'+str(maxlat)+'N.png', bbox_inches = 'tight', format = 'png', dpi=400)
	plt.savefig('/scratch/mp586/Code/Graphics/PE_decomp_cont_size_bucket_and_vp05_'+str(minlat)+'-'+str(maxlat)+'N.pdf', bbox_inches = 'tight', format = 'pdf', dpi=400)




	field_bucket = [ds[3].dPE, - ds[3].dTH_ci,  - ds[3].dMCD_ci,  - ds[3].dTE_ci, - ds[3].residual]
	field_vp05 = [dsvp[3].dPE, - dsvp[3].dTH_ci,  - dsvp[3].dMCD_ci,  - dsvp[3].dTE_ci, - dsvp[3].residual]

	names = ['$\Delta$(P-E)','$\Delta$TH','$\Delta$MCD','$\Delta$TE','R']

	landfile=Dataset(os.path.join(GFDL_BASE,'input/square_South_America/land.nc'),mode='r')                                                         
	landmask=landfile.variables['land_mask'][:]                             
	landlats=landfile.variables['lat'][:]                                   
	landlons=landfile.variables['lon'][:]                                                                                 
	landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it     
	landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))
	landmask, lons_cyclic = addcyclic(landmask, landlons)

	v = np.linspace(-1.,1.,41)
	fig, axes = plt.subplots(5,2, figsize=(16,10))

	for j in range(2):
		if j == 0:
			field = field_bucket
		else:
			field = field_vp05
		for i in range(len(field)):
			array = field[i]
			axes[i,0].set_ylabel(names[i], size = med)


			#fig = plt.figure()

			lats = ds[3].lat
			lons = ds[3].lon

			m = Basemap(projection='cyl',resolution='c', ax = axes[i,j],llcrnrlat=-40, urcrnrlat=40,llcrnrlon=-180, urcrnrlon=180)
			array = xr.DataArray(array,coords=[lats,lons],dims=['lat','lon'])

			array = np.asarray(array)
			array, lons_cyclic = addcyclic(array, lons)
			array,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))

			array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])

			lons = lons_cyclic
			m.drawparallels(np.arange(-40.,40.,20.))
			m.drawmeridians(np.arange(-180.,180.,60.))


			lon, lat = np.meshgrid(lons, lats)
			xi, yi = m(lon, lat)

			cs = m.contourf(xi,yi,array, v, cmap='bwr_r', extend = 'both')

			if np.any(landmask != 0.):
			    m.contour(xi,yi,landmask, 1)


	axes[0,0].set_title('(a) AM 100%cond', fontsize = med)
	axes[0,1].set_title('(b) AM 50%cond', fontsize = med)

	# plt.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.9, wspace=0.02, hspace=0.02)
	# cb_ax = plt.axes([0.3, 0.05, 0.4, 0.03])
	# cbar = plt.colorbar(cs, cax=cb_ax, orientation = 'horizontal')


	cbar = fig.colorbar(cs, orientation = 'vertical', ax = axes, shrink = 0.5)
	cbar.set_label('mm/d', size = med)
	cbar.ax.tick_params(labelsize=med)

	plt.savefig('/scratch/mp586/Code/Graphics/PE_decomp_AM_maps.png', bbox_inches = 'tight', format = 'png', dpi=400)
	plt.savefig('/scratch/mp586/Code/Graphics/PE_decomp_AM_maps.pdf', bbox_inches = 'tight', format = 'pdf', dpi=400)


