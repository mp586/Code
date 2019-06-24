import numpy as np
from netCDF4 import Dataset
#### This is really really slow for large files like the actual bucket depth file 'bucket_depth_exp_slightlyhigherevap.txt'

with open('/scratch/mp586/GFDL_DATA/bucket_depth_exp_slightlyhigherevap.txt') as f:
    lines = f.read().splitlines()
    # lines is now a list with each entry representing all bucket depths for the respective timestep
    
    dummy = lines[0].split()
    a = np.empty([np.shape(lines)[0],np.shape(dummy)[0]])
    for i in range (0,np.shape(lines)[0]):
        a[i,:] = lines[i].split()
        i += 1
    dataset = Dataset('/scratch/mp586/GFDL_DATA/bucket_depth_exp_slightlyhigherevap.nc','w',format='NETCDF4_CLASSIC')
    time = dataset.createDimension('time',None)
    latlon = dataset.createDimension('latlon',np.shape(a)[1])
    bucket_depth = dataset.createVariable('bucket_depth',np.dtype(a[0,0]),('time','latlon'))
    bucket_depth = a

    
