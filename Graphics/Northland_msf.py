from netCDF4 import Dataset
import numpy as np
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


# def msf_calc(nc_file,v='vcomp', a=6376.0e3, g=9.8):
#     """Calculate][] the mass streamfunction for the atmosphere.
#     Based on a vertical integral of the meridional wind.
#     Ref: Physics of Climate, Peixoto & Oort, 1992.  p158.
#     Parameters
#     ----------
#     data :  xarray.DataSet
#         GCM output data
#     v : str, optional
#         The name of the meridional flow field in `data`.  Default: 'vcomp'
#     a : float, optional
#         The radius of the planet. Default: Earth 6317km
#     g : float, optional
#         Surface gravity. Default: Earth 9.8m/s^2
#     Returns
#     -------
#     streamfunction : xarray.DataArray
#         The meridional mass streamfunction.

#     Author = James Penn
#     """
#     data = xr.open_dataset(nc_file+'_interp.nc')
#     data_phalf = xr.open_dataset(nc_file+'.nc')
#     vbar = data[v].mean('lon')
#     c = 2*np.pi*a*np.cos(vbar.lat*np.pi/180) / g
# # take a diff of half levels, and assign to pfull coordinates
#     dp = xr.DataArray(data_phalf.phalf[::-1].diff('phalf').values*100, coords=[('pfull', data.pfull)])
#     msf = (np.cumsum(vbar*dp, axis=vbar.dims.index('pfull')))*c
#     # msf[i-runmin] = (np.cumsum(vbar*dp, axis='pfull')*c)
#     # why cumsum and not # (np.sum(vbar*dp, axis = vbar.dims.index('pfull')))*c
#     msf=xr.DataArray(msf[0,:,:]*1e-10,coords=[data.pfull,data.lat],dims=['pfull','lat'])
#     return (msf)


# # In [57]: for i in range(runmin,runmax):
# #     ...:         runnr="{0:04}".format(i)
# #     ...:         data = xr.open_dataset('/scratch/mp586/'+testdir+'/run'+runnr+'/atmos_monthly.nc')
# #     ...:         if i == runmin:
# #     ...:             msf = np.empty_like(data[v].mean('lon'))
# #     ...:         vbar = data[v].mean('lon')
# #     ...:         msf = msf.repeat((runmax-runmin),axis=0)
# #     ...:         c = 2*np.pi*a*np.cos(vbar.lat*np.pi/180) / g
# #     ...: # take a diff of half levels, and assign to pfull coordinates
# #     ...:         dp = xr.DataArray(data.phalf.diff('phalf').values*100, coords=[('pfull', data.pfull)])
# #     ...:         print(i-runmin)
# #     ...:         msf[i-runmin] = (np.cumsum(vbar*dp, axis='pfull')*c)[0]
# #     ...:      

# def plot_msf(msf,units='10^10 kg/s'):

# 	matplotlib.rcParams['contour.negative_linestyle']= 'dashed'

# 	lats = msf.lat
# 	pfull = msf.pfull
	
# 	y, p = np.meshgrid(lats, pfull)

# 	v = np.linspace(-15.,15.,21)
# 	cset1 = plt.contourf(y, p, msf, v, cmap='RdBu_r', extend = 'both')

# 	cont = plt.contour(y,p,msf, 10, colors = 'k', linewidth=1.)
# 	# plt.clabel(cont, inline=2, fmt='%1.1f',fontsize=14)
# 	cbar = plt.colorbar(cset1)
# 	cbar.set_label(units)
# 	plt.xlabel('Latitude N')
# 	plt.ylabel('Pressure (hPa)')
# 	plt.gca().invert_yaxis()

# 	plt.savefig('/scratch/mp586/Code/Graphics/Northland_msf_'+name+'.png')
# 	plt.close()




if __name__ == '__main__':

    # names = ['LandworldLakes_acdc_atmos_mean_5_to_50',    'Northland_dark_acdc_atmos_mean_5_to_50',
    # 'Landworld_acdc_atmos_mean_5_to_50',         'Northland_dry_acdc_atmos_mean_5_to_50',
    # 'Northland_bright_acdc_atmos_mean_5_to_50',  'Northland_empty_acdc_atmos_mean_5_to_50']
    # for name in names:
    # 	msf = msf_calc('/scratch/mp586/Isca_DATA/Northland/'+name)
    # 	plot_msf(msf, name)

    names = ['Northland_bright_acdc', 'Northland_dry_acdc', 'Northland_dark_acdc', 'Aqua2_acdc']

    msf, msf_avg, bright_msf_seasonal, msf_month_avg = mass_streamfunction_interp('Isca_DATA/Northland/Northland_bright_acdc', 'isca', 121, 481)
    msf, msf_avg, dark_msf_seasonal, msf_month_avg = mass_streamfunction_interp('Isca_DATA/Northland/Northland_dark_acdc', 'isca', 121, 481)
    msf, msf_avg, dry_msf_seasonal, msf_month_avg = mass_streamfunction_interp('Isca_DATA/Northland/Northland_dry_acdc', 'isca', 121, 481)
    msf, msf_avg, aqua_msf_seasonal, msf_month_avg = mass_streamfunction_interp('Isca_DATA/Northland/Aqua2_acdc', 'isca', 121, 481)

    plot_streamfunction_seasonal((dark_msf_seasonal - bright_msf_seasonal)*10, 'Northland_dark_minus_bright', 121, 481, units='10^9 kg/s', minval = -420, maxval =420, ctl_pert = 'ctl')
    plot_streamfunction_seasonal((dry_msf_seasonal - bright_msf_seasonal)*10, 'Northland_dry_minus_bright', 121, 481, units='10^9 kg/s', minval = -400, maxval =400, ctl_pert = 'ctl')
    plot_streamfunction_seasonal((dark_msf_seasonal)*10, 'Northland_dark', 121, 481, units='10^9 kg/s', minval = -420, maxval =420, ctl_pert = 'ctl')
    plot_streamfunction_seasonal((aqua_msf_seasonal)*10, 'Aqua2', 121, 481, units='10^9 kg/s', minval = -420, maxval =420, ctl_pert = 'ctl')





