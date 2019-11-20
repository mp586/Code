import xarray as xar
import os

base_dir = os.environ['GFDL_DATA']
base_dir = base_dir + '/ISCA_HPC'

# set name of experiment
exp_names_ctl = ['full_continents_newbucket_fullnewbucketqflux_without_pogfix_commit7bb4387',
'full_continents_newbucket_fixedSSTs_zonally_symmetric_commit7bb4387',
'two_continents_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387',
'two_continents_newbucket_fixedSSTs_zonally_symmetric_EqMax_commit7bb4387',
'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387',
'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387',
'full_continents_newbucket_fixedSSTs_zonally_symmetric_EqMax_commit7bb4387',
'two_continents_shortSA_newbucket_fixedSSTs_from_realworld_zonallysymm_commitfa0ee3e']

#month number of first and last months to be included
start_file = 121
end_file = 480

file_name = 'atmos_monthly.nc'


for exp_name in exp_names_ctl:


	#Generate list of all the files we want to include in our dataset
	files = ['/run%04d/' % m for m in range(start_file, end_file+1)]
	final_files = [base_dir+'/'+exp_name+f+file_name for f in files]

	#Open xar dataset object using all of those files
	dataset = xar.open_mfdataset(final_files, decode_times = False)

	# Define a list of variables you wish to DROP
	list_to_drop = ['vor', 'div','sphum_u','sphum_v','sphum_w','bucket_depth_lh','bucket_depth_conv','bucket_depth_cond']

	temp_dataset = dataset

	#Go through this list and create a new dataset called 'temp_dataset' without your dropped variable
	for var in list_to_drop:
	    print(var)
	    temp_dataset = temp_dataset.drop(var)

	#Once we have dropped all the variables you want, then we can output it
	print('dropped')

	#I've chosen to put the single output file into the data directory for that experiment, but it could go anywhere
	temp_dataset.to_netcdf(base_dir+'/one_files/'+exp_name+'.nc')

exp_names_pertb = ['full_continents_newbucket_fullnewbucketqflux_2xCO2_spinup361_without_pogfix_commit7bb4387',
'full_continents_newbucket_fixedSSTs_zonally_symmetric_plus_2pt52K_and_2xCO2_spinup_361_commit7bb4387',
'two_continents_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'two_continents_newbucket_fixedSSTs_zonally_symmetric_EqMax_plus_2pt52K_and_2xCO2_spinup_361_commit7bb4387',
'square_South_America_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'square_Africa_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387',
'full_continents_newbucket_fixedSSTs_zonally_symmetric_EqMax_plus_2pt52K_and_2xCO2_spinup_361_commit7bb4387',
'two_continents_shortSA_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commitfa0ee3e']
#month number of first and last months to be included
start_file = 120
end_file = 479

file_name = 'atmos_monthly.nc'


for exp_name in exp_names_pertb:


	#Generate list of all the files we want to include in our dataset
	files = ['/run%04d/' % m for m in range(start_file, end_file+1)]
	final_files = [base_dir+'/'+exp_name+f+file_name for f in files]

	#Open xar dataset object using all of those files
	dataset = xar.open_mfdataset(final_files, decode_times = False)

	# Define a list of variables you wish to DROP
	list_to_drop = ['vor', 'div','sphum_u','sphum_v','sphum_w','bucket_depth_lh','bucket_depth_conv','bucket_depth_cond']

	temp_dataset = dataset

	#Go through this list and create a new dataset called 'temp_dataset' without your dropped variable
	for var in list_to_drop:
	    print(var)
	    temp_dataset = temp_dataset.drop(var)

	#Once we have dropped all the variables you want, then we can output it
	print('dropped')

	#I've chosen to put the single output file into the data directory for that experiment, but it could go anywhere
	temp_dataset.to_netcdf(base_dir+'/one_files/'+exp_name+'.nc')



exp_name = 'full_continents_lepref07_fixedSSTs_zonally_symmetric_plus_2pt52K_and_2xCO2_spinup_361'

start_file = 120
end_file = 479

file_name = 'atmos_monthly.nc'

#Generate list of all the files we want to include in our dataset
files = ['/run%04d/' % m for m in range(start_file, end_file+1)]
final_files = [base_dir+'/'+exp_name+f+file_name for f in files]

#Open xar dataset object using all of those files
dataset = xar.open_mfdataset(final_files, decode_times = False)

# Define a list of variables you wish to DROP
list_to_drop = ['vor', 'div','sphum_u','sphum_v','sphum_w']

temp_dataset = dataset

#Go through this list and create a new dataset called 'temp_dataset' without your dropped variable
for var in list_to_drop:
    print(var)
    temp_dataset = temp_dataset.drop(var)

#Once we have dropped all the variables you want, then we can output it
print('dropped')

#I've chosen to put the single output file into the data directory for that experiment, but it could go anywhere
temp_dataset.to_netcdf(base_dir+'/one_files/'+exp_name+'.nc')

exp_name = 'full_continents_lepref07_fixedSSTs_zonally_symmetric'

start_file = 121
end_file = 480

file_name = 'atmos_monthly.nc'

# Generate list of all the files we want to include in our dataset
files = ['/run%04d/' % m for m in range(start_file, end_file+1)]
final_files = [base_dir+'/'+exp_name+f+file_name for f in files]

#Open xar dataset object using all of those files
dataset = xar.open_mfdataset(final_files, decode_times = False)

# Define a list of variables you wish to DROP
list_to_drop = ['vor', 'div','sphum_u','sphum_v','sphum_w']

temp_dataset = dataset

#Go through this list and create a new dataset called 'temp_dataset' without your dropped variable
for var in list_to_drop:
    print(var)
    temp_dataset = temp_dataset.drop(var)

#Once we have dropped all the variables you want, then we can output it
print('dropped')

#I've chosen to put the single output file into the data directory for that experiment, but it could go anywhere
temp_dataset.to_netcdf(base_dir+'/one_files/'+exp_name+'.nc')
