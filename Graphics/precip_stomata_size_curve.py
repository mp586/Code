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


#control dirs 
control_sb_dirs = [
'narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387',
'narrow_six_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387',
'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387',
'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387'
]

control_vp0_dirs = [
'narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref0_commit7bb4387',
'x',
'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref0_commit7bb4387',
'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref0_commit7bb4387'
]

control_vp02_dirs = [
'x',
'x',
'x',
'x'

]
control_vp05_dirs = [
'narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_commit7bb4387',
'x',
'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_commit7bb4387',
'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_commit7bb4387'
]

# perturbed dirs
vp0_dirs = [
'narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref0_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'narrow_six_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref0_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref0_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref0_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387'
]

vp02_dirs = [
'x',
'x',
'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref02_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref02_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387'
]

vp05_dirs = [
'narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'x',
'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref05_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387'
]

simple_bucket_dirs = [
'narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'x',
'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387'
]

[precipitation_ctl,precipitation_avg_ctl,x,x,x]=seasonal_surface_variable('Isca_DATA/ISCA_HPC/narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387','isca',1,10,'precipitation','mm/d', factor=86400)


lats = precipitation_avg_ctl.lat
lons = precipitation_avg_ctl.lon


landmasks = ['narrow_three','narrow_six','square_South_America','square_Africa']

landfile=Dataset(os.path.join(GFDL_BASE,'input/'+landmasks[0]+'/land.nc'),mode='r')
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]

landmask_array = np.empty((4, len(landlats), len(landlons)))


for i in range(len(landmasks)):
	landfile=Dataset(os.path.join(GFDL_BASE,'input/'+landmasks[i]+'/land.nc'),mode='r')
	landmask_array[i,:,:]=landfile.variables['land_mask'][:]


ctl_dict = {'vp0_ctl': control_vp0_dirs, 
'vp02_ctl': control_vp02_dirs,
'vp05_ctl': control_vp05_dirs, 
'sb_ctl': control_sb_dirs}

ctl_list = ['vp0_ctl','vp02_ctl','vp05_ctl','sb_ctl']


pert_dict = {'vp0': vp0_dirs, 
'vp02': vp02_dirs,
'vp05': vp05_dirs, 
'sb': simple_bucket_dirs}

pert_list = ['vp0','vp02','vp05','sb']

precip_pert_matrix = np.empty((len(vp0_dirs),len(pert_list)))
precip_ctl_matrix = np.empty((len(vp0_dirs),len(pert_list)))

# fig = plt.figure(figsize=(20,20))

minlats = [-10.,-10.,0.]
maxlats = [10.,0.,10.]
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
				precip_ctl_matrix[i,j] = area_weighted_avg(precipitation_avg_ctl,area_array,landmaskxr,'land', minlat = minlat, maxlat = maxlat)

			
			testdir = pert_dict[pert_list[j]][i]
			if testdir != 'x':
				testdir = 'Isca_DATA/ISCA_HPC/'+testdir
				[precipitation,precipitation_avg,x,x,x]=seasonal_surface_variable(testdir,'isca',120,480,'precipitation','mm/d', factor=86400)
				precip_pert_matrix[i,j] = area_weighted_avg(precipitation_avg,area_array,landmaskxr,'land', minlat = minlat, maxlat = maxlat)


	precip_ctl_matrix[precip_ctl_matrix == 0] = 'nan'
	precip_pert_matrix[precip_pert_matrix == 0] = 'nan'

	for j in range(len(pert_list)):
		plt.plot([3,6,40,60],precip_pert_matrix[:,j], '*', label = pert_list[j])
		plt.plot([3,6,40,60],precip_ctl_matrix[:,j], '.', label = ctl_list[j])

	plt.legend()
	plt.xlabel('Continental extent ($^{\circ}$ lon)')
	plt.ylabel('Precipitation per Area (mm/d/m$^2$)')
	plt.savefig('/scratch/mp586/Code/Graphics/P_ctl_P_pert_stomata_cont_size_'+str(minlat)+'-'+str(maxlat)+'N.png', bbox_inches = 'tight', format = 'png', dpi=400)
	plt.close()

	for j in range(len(pert_list)):
		plt.plot([3,6,40,60],(precip_pert_matrix[:,j] - precip_ctl_matrix[:,3])/precip_ctl_matrix[:,3], '.', label = pert_list[j])
	plt.legend()
	plt.plot([3,6,40,60],[0,0,0,0],'k')
	plt.ylim([-1., 1.])
	plt.xlabel('Continental extent ($^{\circ}$ lon)')
	plt.ylabel('Rel. $\Delta$P per Area (mm/d/m$^2$)')
	plt.savefig('/scratch/mp586/Code/Graphics/P_rel_change_stomata_cont_size_'+str(minlat)+'-'+str(maxlat)+'N.png', bbox_inches = 'tight', format = 'png', dpi=400)
	plt.close()

	for j in range(len(pert_list)):
		plt.plot([3,6,40,60],(precip_pert_matrix[:,j] - precip_ctl_matrix[:,3]), '.', label = pert_list[j])
	plt.legend()
	plt.plot([3,6,40,60],[0,0,0,0],'k')
	plt.ylim([-4., 4.])
	plt.xlabel('Continental extent ($^{\circ}$ lon)')
	plt.ylabel('Abs. $\Delta$P per Area (mm/d/m$^2$)')
	plt.savefig('/scratch/mp586/Code/Graphics/P_change_stomata_cont_size_'+str(minlat)+'-'+str(maxlat)+'N.png', bbox_inches = 'tight', format = 'png', dpi=400)
	plt.close()



	warming_only = precip_pert_matrix[:,3] - precip_ctl_matrix[:,3]
	for j in range(len(pert_list) - 1):
		stomata_only = precip_ctl_matrix[:,j] - precip_ctl_matrix[:,3]
		addition = warming_only + stomata_only
		full = precip_pert_matrix[:,j] - precip_ctl_matrix[:,3]

		plt.plot([3,6,40,60],[0,0,0,0],'k')
		plt.plot([3,6,40,60],warming_only,'b*', label = 'warming (SB_p - SB_c)')
		plt.plot([3,6,40,60],stomata_only,'g*', label = 'stomata (vp_c - SB_c)')
		plt.plot([3,6,40,60],addition,'r*', label = 'addition')
		plt.plot([3,6,40,60],full,'m*', label = 'full change (vp_p - SB_c)')
		plt.legend()
		plt.ylim([-4., 4.])
		plt.xlabel('Continental extent ($^{\circ}$ lon)')
		plt.ylabel('Abs. $\Delta$P per Area (mm/d/m$^2$)')
		plt.title(''+pert_list[j])
		plt.savefig('/scratch/mp586/Code/Graphics/P_change_linadd_stomata_'+pert_list[j]+'_cont_size_'+str(minlat)+'-'+str(maxlat)+'N.png', bbox_inches = 'tight', format = 'png', dpi=400)
		plt.close()

