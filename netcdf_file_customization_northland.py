import xarray as xar
import os

base_dir = os.environ['GFDL_DATA']
base_dir = base_dir + '/Northland'

# set name of experiment
exp_names_ctl = ['Aqua2_acdc']

#month number of first and last months to be included
start_file = 1
end_file = 600

file_name = 'atmos_monthly_interp.nc'


for exp_name in exp_names_ctl:


	#Generate list of all the files we want to include in our dataset
	files = ['/run%04d/' % m for m in range(start_file, end_file+1)]
	final_files = [base_dir+'/'+exp_name+f+file_name for f in files]

	#Open xar dataset object using all of those files
	dataset = xar.open_mfdataset(final_files, decode_times = False)

	# Define a list of variables you wish to DROP
	list_to_drop = ['vor', 'div','bucket_depth_lh','bucket_depth_conv', 'precipitation', 'bucket_depth', 
	't_surf', 'flux_sw', 'flux_lw', 'toa_sw', 'rrtm_albedo', 'co2', 'ml_heat_cap', 'olr', 'flux_lhe', 'flux_t', 'omega', 'rh', 'zsurf']

	temp_dataset = dataset

	#Go through this list and create a new dataset called 'temp_dataset' without your dropped variable
	for var in list_to_drop:
	    print(var)
	    temp_dataset = temp_dataset.drop(var)

	#Once we have dropped all the variables you want, then we can output it
	print('dropped')

	#I've chosen to put the single output file into the data directory for that experiment, but it could go anywhere
	temp_dataset.to_netcdf(base_dir+'/one_files/'+exp_name+'.nc')

