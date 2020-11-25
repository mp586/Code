from netCDF4 import Dataset
import numpy as np
import matplotlib as mpl 
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

area_array, dx, dy = ca.cell_area_all(t_res=42,base_dir='/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/') # added _all because then dx and dy are also returned 
area_array = xr.DataArray(area_array) # returned in units of m bzw m^2, because radius in cell_area.py is given in metres

control_sb_dirs = [
'narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387',
'narrow_six_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387',
'narrow_twelve_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387',
'narrow_twentyfour_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387',
'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387',
'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387',
'squareland_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387'
]

control_vp0_dirs = [
'narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref0_commit7bb4387',
'x',
'x',
'x',
'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref0_commit7bb4387',
'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref0_commit7bb4387',
'x'
]

control_vp02_dirs = [
'x',
'x',
'x',
'x',
'x',
'x',
'x'
]

control_vp05_dirs = [
'narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_commit7bb4387',
'x',
'narrow_twelve_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_commit7bb4387',
'narrow_twentyfour_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_commit7bb4387',
'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_commit7bb4387',
'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_commit7bb4387',
'squareland_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_commit7bb4387'
]

control_vp07_dirs = [
'x',
'x',
'x',
'x',
'x',
'x',
'x'
]



vp0_dirs = [
'narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref0_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'narrow_six_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref0_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'x',
'x',
'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref0_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref0_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'x'
]

vp02_dirs = [
'x',
'x',
'x',
'x',
'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref02_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref02_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'x'
]

vp05_dirs = [
'narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'x',
'narrow_twelve_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'narrow_twentyfour_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'squareland_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387'
]

vp07_dirs = [
'x',
'x',
'x',
'x',
'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref07_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref07_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'x'
]


simple_bucket_dirs = [
'narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'x',
'narrow_twelve_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'narrow_twentyfour_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'squareland_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387'
]




[precipitation_ctl,precipitation_avg_ctl,x,x,x]=seasonal_surface_variable('Isca_DATA/ISCA_HPC/narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387','isca',1,10,'precipitation','mm/d', factor=86400)


lats = precipitation_avg_ctl.lat
lons = precipitation_avg_ctl.lon

landmasks = ['narrow_three','narrow_six','narrow_twelve','narrow_twentyfour','square_South_America','square_Africa','squareland']

landfile=Dataset(os.path.join(GFDL_BASE,'input/'+landmasks[0]+'/land.nc'),mode='r')
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]

landmask_array = np.empty((len(landmasks), len(landlats), len(landlons)))



for i in range(len(landmasks)):
	landfile=Dataset(os.path.join(GFDL_BASE,'input/'+landmasks[i]+'/land.nc'),mode='r')
	landmask_array[i,:,:]=landfile.variables['land_mask'][:]


ctl_dict = {'CV0_ctl': control_vp0_dirs, 
'CV02_ctl': control_vp02_dirs,
'CV05_ctl': control_vp05_dirs, 
'CV07_ctl': control_vp05_dirs, 
'SB_ctl': control_sb_dirs}

ctl_list = ['CV0_ctl','CV02_ctl','CV05_ctl','CV07_ctl','SB_ctl']


pert_dict = {'CV0': vp0_dirs, 
'CV02': vp02_dirs,
'CV05': vp05_dirs, 
'CV07' : vp07_dirs,
'SB': simple_bucket_dirs}

pert_list = ['CV0','CV02','CV05','CV07','SB']
pert_list_names = ['0%cond','20%cond','50%cond','70%cond','bucket']

precip_pert_matrix = np.zeros((len(vp0_dirs),len(pert_list)))
precip_ctl_matrix = np.zeros((len(vp0_dirs),len(pert_list)))

precip_pert_sds = np.zeros((len(vp0_dirs),len(pert_list)))
precip_ctl_sds = np.zeros((len(vp0_dirs),len(pert_list)))

med = 24
# fig = plt.figure(figsize=(20,20))

minlats = [-10.]
maxlats = [10.]
for k in range(len(minlats)):
	minlat = minlats[k]
	maxlat = maxlats[k]
	for i in range(len(vp0_dirs)):
		landmaskxr=xr.DataArray(landmask_array[i,:,:],coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it
		for j in range(len(pert_list)):
			testdir = ctl_dict[ctl_list[j]][i]
			if testdir != 'x':
				testdir = 'Isca_DATA/ISCA_HPC/'+testdir
				[precipitation_ctl,precipitation_avg_ctl,x,x,x]=seasonal_surface_variable(testdir,'isca',121,481,'precipitation','mm/d', factor=86400)
				precip_ctl_matrix[i,j], precip_ctl_sds[i,j] = area_weighted_avg(precipitation_avg_ctl,area_array,landmaskxr,'land', minlat = minlat, maxlat = maxlat,return_sd = True)

			
			testdir = pert_dict[pert_list[j]][i]
			if testdir != 'x':
				testdir = 'Isca_DATA/ISCA_HPC/'+testdir
				[precipitation,precipitation_avg,x,x,x]=seasonal_surface_variable(testdir,'isca',120,480,'precipitation','mm/d', factor=86400)
				precip_pert_matrix[i,j], precip_pert_sds[i,j] = area_weighted_avg(precipitation_avg,area_array,landmaskxr,'land', minlat = minlat, maxlat = maxlat,return_sd = True)


	precip_ctl_matrix[precip_ctl_matrix == 0] = 'nan'
	precip_pert_matrix[precip_pert_matrix == 0] = 'nan'

	# fig, axes = plt.subplots(1,1, figsize = (10,10))

	# for j in range(len(pert_list)):
	# 	axes.plot([3,6,12,24,40,60,100],precip_ctl_matrix[:,j], '.', label = ctl_list[j])
	# 	axes.plot([3,6,12,24,40,60,100],precip_pert_matrix[:,j], '*', label = pert_list[j])
	# axes.spines['right'].set_visible(False)
	# axes.spines['top'].set_visible(False)
	# axes.tick_params(labelsize = med)
	# axes.tick_params(labelsize = med)
	# axes.set_xlabel('Continental extent ($^{\circ}$ lon)', fontsize = med)
	# axes.set_ylabel('$P$ per Area (mm/d/m$^2$)', fontsize = med)
	# fig.legend(fontsize = med, bbox_to_anchor=(0.83,0.83))
	# fig.savefig('/scratch/mp586/Code/Graphics/P_ctl_P_pert_stomata_cont_size_'+str(minlat)+'-'+str(maxlat)+'N.png', bbox_inches = 'tight', format = 'png', dpi=400)
	# plt.close()


	precip_ctl_matrix_del6 = np.delete(precip_ctl_matrix,1,0)
	precip_pert_matrix_del6 = np.delete(precip_pert_matrix,1,0)

	# fig, axes = plt.subplots(1,1, figsize = (10,10))

	# axes.plot([3,12,24,40,60,100],precip_ctl_matrix_del6[:,4], 's', markersize = 10., label = 'SB ctl')
	# axes.plot([3,12,24,40,60,100],precip_pert_matrix_del6[:,4], 'o', markersize = 10., label = 'SB pert')
	# axes.plot([3,12,24,40,60,100],precip_pert_matrix_del6[:,2], 'v', markersize = 10., label = 'CV05 pert')


	# axes.spines['right'].set_visible(False)
	# axes.spines['top'].set_visible(False)
	# axes.tick_params(labelsize = med)
	# axes.tick_params(labelsize = med)
	# axes.set_xlabel('Continental extent ($^{\circ}$ lon)', fontsize = med)
	# axes.set_ylabel('$P$ per Area (mm/d/m$^2$)', fontsize = med)
	# axes.set_ylim(0.,8.)
	# axes.set_xlim(0.,110.)
	# axes.set_xticks([3,12,24,40,60,100])
	# fig.legend(fontsize = med, loc = 'lower left', bbox_to_anchor=(0.1,0.1))
	# fig.savefig('/scratch/mp586/Code/Graphics/P_ctl_P_pert_SBctl_SBpert_VP05_cont_size_'+str(minlat)+'-'+str(maxlat)+'N.png', bbox_inches = 'tight', format = 'png', dpi=400)
	# fig.savefig('/scratch/mp586/Code/Graphics/P_ctl_P_pert_SBctl_SBpert_VP05_cont_size_'+str(minlat)+'-'+str(maxlat)+'N.pdf', bbox_inches = 'tight', format = 'pdf', dpi=400)
	# plt.close()



	warming_only = precip_pert_matrix_del6[:,4] - precip_ctl_matrix_del6[:,4]
	for j in range(len(pert_list) - 1):
		fig, axes = plt.subplots(1,2, sharex = True, figsize = (23,10))
		stomata_only = precip_ctl_matrix_del6[:,j] - precip_ctl_matrix_del6[:,4]
		addition = warming_only + stomata_only
		full = precip_pert_matrix_del6[:,j] - precip_ctl_matrix_del6[:,4]

		axes[0].plot([6,14,25,40,60,100],precip_ctl_matrix_del6[:,4], 's', color = 'slategrey', markersize = 10., label = 'bucket ctl')
		axes[0].plot([6,14,25,40,60,100],precip_pert_matrix_del6[:,4], 'o', color = 'darkred', markersize = 10., label = 'bucket pert')
		axes[0].plot([6,14,25,40,60,100],precip_pert_matrix_del6[:,j], 'v', color = 'salmon', markersize = 10., label = '50%cond pert')


		axes[0].spines['right'].set_visible(False)
		axes[0].spines['top'].set_visible(False)
		axes[0].tick_params(labelsize = med)
		axes[0].tick_params(labelsize = med)
		axes[0].set_xlabel('Continental extent ($^{\circ}$ lon)', fontsize = med)
		axes[0].set_ylabel('P per Area (mm/d/m$^2$)', fontsize = med)
		axes[0].set_ylim([0.,8.])
		axes[0].set_xlim(0.,110.)
		axes[0].set_xticks([6,14,25,40,60,100])
		axes[0].legend(fontsize = med, loc = 'lower left')
		axes[0].set_title('(a) P vs continental extent', fontsize = med)


		axes[1].plot([6,14,25,40,60,100],[0,0,0,0,0,0],'k')
		axes[1].plot([6,14,25,40,60,100],warming_only,'p', color='deepskyblue', markersize = 10.,label = 'bucket (pert - ctl)')
		axes[1].plot([6,14,25,40,60,100],stomata_only,'P',color = 'seagreen',markersize = 10.,label = 'stomata ('+pert_list_names[j]+' ctl - bucket ctl)')
		axes[1].plot([6,14,25,40,60,100],addition,'D', color='lightgreen', markersize = 10.,label = 'bucket + stomata')
		axes[1].plot([6,14,25,40,60,100],full,'^', color='navy', markersize = 10.,label = 'full change ('+pert_list_names[j]+' pert - bucket ctl)')
		axes[1].legend(fontsize = med, loc = 'lower right')

		axes[1].spines['right'].set_visible(False)
		axes[1].spines['top'].set_visible(False)
		axes[1].tick_params(labelsize = med)
		axes[1].tick_params(labelsize = med)
		axes[1].set_ylim([-2., 2.])
		axes[1].set_title('(b) $\Delta$ P decomposition vs continental extent', fontsize = med)
		axes[1].set_xlabel('Continental extent ($^{\circ}$ lon)',fontsize = med)
		axes[1].set_ylabel('$\Delta$ P per Area (mm/d/m$^2$)',fontsize = med)
		fig.savefig('/scratch/mp586/Code/Graphics/P_stomata_'+pert_list[j]+'_v_warming_and_contsize_'+str(minlat)+'-'+str(maxlat)+'N_paper.png', bbox_inches = 'tight', format = 'png', dpi=400)
		fig.savefig('/scratch/mp586/Code/Graphics/P_stomata_'+pert_list[j]+'_v_warming_and_contsize_'+str(minlat)+'-'+str(maxlat)+'N_paper.pdf', bbox_inches = 'tight', format = 'pdf', dpi=400)
		fig.savefig('/scratch/mp586/Code/Graphics/P_stomata_'+pert_list[j]+'_v_warming_and_contsize_'+str(minlat)+'-'+str(maxlat)+'N_paper.eps', bbox_inches = 'tight', format = 'eps', dpi=600)
		plt.close()



		print(pert_list[j])
		[slope, intercept, r_value, p_value, std_err] = stats.linregress([6,14,25,40,60,100],precip_ctl_matrix_del6[:,4])
		print("slope linregress ctl = "+str(slope))
		[slope, intercept, r_value, p_value, std_err] = stats.linregress([6,14,25,40,60,100],precip_pert_matrix_del6[:,4])
		print("slope linregress bucket = "+str(slope))
		[slope, intercept, r_value, p_value, std_err] = stats.linregress([6,14,25,40,60,100],precip_pert_matrix_del6[:,j])
		print("slope linregress veg = "+str(slope))

		[slope_lst, intercept_lst, x, x, x] = orthoregress([6,14,25,40,60,100],precip_ctl_matrix_del6[:,4])
		print("slope totalregress ctl= "+str(slope_lst))
		[slope_lst, intercept_lst, x, x, x] = orthoregress([6,14,25,40,60,100],precip_pert_matrix_del6[:,4])
		print("slope totalregress pert= "+str(slope_lst))
		[slope_lst, intercept_lst, x, x, x] = orthoregress([6,14,25,40,60,100],precip_pert_matrix_del6[:,j])
		print("slope totalregress veg= "+str(slope_lst))

# CV05
# slope linregress ctl = -0.063495873788
# slope linregress bucket = -0.0674657883614
# slope linregress veg = -0.0656376330027
# slope totalregress ctl= -0.0635055054459
# slope totalregress pert= -0.0674703639336
# slope totalregress veg= -0.0656417572586


	for j in range(len(pert_list) - 1):
		fig, axes = plt.subplots(1,1, figsize = (13,10))
		stomata_only = precip_ctl_matrix_del6[:,j] - precip_ctl_matrix_del6[:,4]
		addition = warming_only + stomata_only
		full = precip_pert_matrix_del6[:,j] - precip_ctl_matrix_del6[:,4]

		axes.plot([6,14,25,40,60,100],precip_ctl_matrix_del6[:,4], 's', color = 'slategrey', markersize = 10., label = 'bucket ctl')
		axes.plot([6,14,25,40,60,100],precip_pert_matrix_del6[:,4], 'o', color = 'darkred', markersize = 10., label = 'bucket pert')
		axes.plot([6,14,25,40,60,100],precip_pert_matrix_del6[:,j], 'v', color = 'salmon', markersize = 10., label = '50%cond pert')


		axes.spines['right'].set_visible(False)
		axes.spines['top'].set_visible(False)
		axes.tick_params(labelsize = med)
		axes.tick_params(labelsize = med)
		axes.set_xlabel('Continental extent ($^{\circ}$ lon)', fontsize = med)
		axes.set_ylabel('P per Area (mm/d/m$^2$)', fontsize = med)
		axes.set_ylim([0.,8.])
		axes.set_xlim(0.,110.)
		axes.set_xticks([6,14,25,40,60,100])
		axes.legend(fontsize = med, loc = 'lower left')
#		axes.set_title('P vs continental extent', fontsize = med)


		# axes[1].plot([6,14,25,40,60,100],[0,0,0,0,0,0],'k')
		# axes[1].plot([6,14,25,40,60,100],warming_only,'p', color='deepskyblue', markersize = 10.,label = 'bucket (pert - ctl)')
		# axes[1].plot([6,14,25,40,60,100],stomata_only,'P',color = 'seagreen',markersize = 10.,label = 'stomata ('+pert_list_names[j]+' ctl - bucket ctl)')
		# axes[1].plot([6,14,25,40,60,100],addition,'D', color='lightgreen', markersize = 10.,label = 'bucket + stomata')
		# axes[1].plot([6,14,25,40,60,100],full,'^', color='navy', markersize = 10.,label = 'full change ('+pert_list_names[j]+' pert - bucket ctl)')
		# axes[1].legend(fontsize = med, loc = 'lower right')

		# axes[1].spines['right'].set_visible(False)
		# axes[1].spines['top'].set_visible(False)
		# axes[1].tick_params(labelsize = med)
		# axes[1].tick_params(labelsize = med)
		# axes[1].set_ylim([-2., 2.])
		# axes[1].set_title('(b) $\Delta$ P decomposition vs continental extent', fontsize = med)
		# axes[1].set_xlabel('Continental extent ($^{\circ}$ lon)',fontsize = med)
		# axes[1].set_ylabel('$\Delta$ P per Area (mm/d/m$^2$)',fontsize = med)
		fig.savefig('/scratch/mp586/Code/Graphics/P_stomata_'+pert_list[j]+'_v_warming_and_contsize_'+str(minlat)+'-'+str(maxlat)+'N_paper_onlya.png', bbox_inches = 'tight', format = 'png', dpi=400)
		fig.savefig('/scratch/mp586/Code/Graphics/P_stomata_'+pert_list[j]+'_v_warming_and_contsize_'+str(minlat)+'-'+str(maxlat)+'N_paper_onlya.pdf', bbox_inches = 'tight', format = 'pdf', dpi=400)
		fig.savefig('/scratch/mp586/Code/Graphics/P_stomata_'+pert_list[j]+'_v_warming_and_contsize_'+str(minlat)+'-'+str(maxlat)+'N_paper_onlya.eps', bbox_inches = 'tight', format = 'eps', dpi=600)
		plt.close()


		fig, axes = plt.subplots(1,1, figsize = (13,10))


		axes.plot([6,14,25,40,60,100],[0,0,0,0,0,0],'k')
		axes.plot([6,14,25,40,60,100],warming_only,'p', color='deepskyblue', markersize = 10.,label = 'bucket (pert - ctl)')
		axes.plot([6,14,25,40,60,100],stomata_only,'P',color = 'seagreen',markersize = 10.,label = 'stomata ('+pert_list_names[j]+' ctl - bucket ctl)')
		axes.plot([6,14,25,40,60,100],addition,'D', color='lightgreen', markersize = 10.,label = 'bucket + stomata')
		axes.plot([6,14,25,40,60,100],full,'^', color='navy', markersize = 10.,label = 'full change ('+pert_list_names[j]+' pert - bucket ctl)')
		axes.legend(fontsize = med, loc = 'lower right')

		axes.spines['right'].set_visible(False)
		axes.spines['top'].set_visible(False)
		axes.tick_params(labelsize = med)
		axes.tick_params(labelsize = med)
		axes.set_ylim([-2., 2.])
		axes.set_xlim(0.,110.)
		axes.set_xticks([6,14,25,40,60,100])
#		axes.set_title('$\Delta$ P decomposition vs continental extent', fontsize = med)
		axes.set_xlabel('Continental extent ($^{\circ}$ lon)',fontsize = med)
		axes.set_ylabel('$\Delta$ P per Area (mm/d/m$^2$)',fontsize = med)
		fig.savefig('/scratch/mp586/Code/Graphics/P_stomata_'+pert_list[j]+'_v_warming_and_contsize_'+str(minlat)+'-'+str(maxlat)+'N_paper_onlyb.png', bbox_inches = 'tight', format = 'png', dpi=400)
		fig.savefig('/scratch/mp586/Code/Graphics/P_stomata_'+pert_list[j]+'_v_warming_and_contsize_'+str(minlat)+'-'+str(maxlat)+'N_paper_onlyb.pdf', bbox_inches = 'tight', format = 'pdf', dpi=400)
		fig.savefig('/scratch/mp586/Code/Graphics/P_stomata_'+pert_list[j]+'_v_warming_and_contsize_'+str(minlat)+'-'+str(maxlat)+'N_paper_onlyb.eps', bbox_inches = 'tight', format = 'eps', dpi=600)
		plt.close()



	for j in range(len(pert_list)):
		plt.plot([3,6,12,24,40,60,100],(precip_pert_matrix[:,j] - precip_ctl_matrix[:,4])/precip_ctl_matrix[:,4], '.', label = pert_list[j])
	plt.legend()
	plt.plot([3,6,12,24,40,60,100],[0,0,0,0,0,0,0],'k')
	plt.ylim([-1., 1.])
	plt.xlabel('Continental extent ($^{\circ}$ lon)')
	plt.ylabel('Rel. $\Delta$P per Area (mm/d/m$^2$)')
	plt.savefig('/scratch/mp586/Code/Graphics/P_rel_change_stomata_cont_size_'+str(minlat)+'-'+str(maxlat)+'N.png', bbox_inches = 'tight', format = 'png', dpi=400)
	plt.close()

	for j in range(len(pert_list)):
		plt.plot([3,6,12,24,40,60,100],(precip_pert_matrix[:,j] - precip_ctl_matrix[:,4]), '.', label = pert_list[j])
	plt.legend()
	plt.plot([3,6,12,24,40,60,100],[0,0,0,0,0,0,0],'k')
	plt.ylim([-4., 4.])
	plt.xlabel('Continental extent ($^{\circ}$ lon)')
	plt.ylabel('Abs. $\Delta$P per Area (mm/d/m$^2$)')
	plt.savefig('/scratch/mp586/Code/Graphics/P_change_stomata_cont_size_'+str(minlat)+'-'+str(maxlat)+'N.png', bbox_inches = 'tight', format = 'png', dpi=400)
	plt.close()





