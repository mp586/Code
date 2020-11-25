from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import xarray as xr
import pandas as pd
import os
from matplotlib.patches import Rectangle

import sys
sys.path.insert(0, '/scratch/mp586/Code/PYCODES')
from plotting_routines_kav7 import *
import stats as st

GFDL_BASE = os.environ['GFDL_BASE']
sys.path.insert(0, os.path.join(GFDL_BASE,'src/extra/python/scripts'))
import cell_area as ca



def quivers_2cases(runmin, runmax, plotname, dire, landmaskxr, cbarlabel, ucomp1_avg, ucomp1_avg_ctl, wcomp1_avg, wcomp1_avg_ctl, ucomp2_avg, ucomp2_avg_ctl, wcomp2_avg, wcomp2_avg_ctl, array1, array2, minval, maxval, vertmult, minlat=-10., maxlat=10., cmap = 'BrBG'): 
    small = 24 #largefonts 14 # smallfonts 10 # medfonts = 14
    med = 28 #largefonts 18 # smallfonts 14 # medfonts = 16
    lge = 30 #largefonts 22 # smallfonts 18 # medfonts = 20
    veclen = 10.
    units_numerator = 'm'
    units_denom = 's'
    shiftby = 180. # = 180. --> 0 degrees in the middle, = 105. --> idealized continents overlap realistic ones 

    area_array, dx, dy = ca.cell_area_all(t_res=42,base_dir='/scratch/mp586/Isca/') # added _all because then dx and dy are also returned 
    area_array = xr.DataArray(area_array) # returned in units of m bzw m^2, because radius in cell_area.py is given in metres
    area_array_3D = np.expand_dims(area_array, axis=0)
    area_array_3D = np.repeat(area_array_3D, 40, axis = 0) # to make area_array 3D (pressure, lat, lon)


    v = np.linspace(minval,maxval,21) # , endpoint=True)
    

    fig, axes = plt.subplots(1, 3, sharex = True, sharey = True, figsize = (18,8))


    # panel 1: Only South America 
    uwind = (ucomp1_avg - ucomp1_avg_ctl) # [::-1,:,:]
    wwind = ((wcomp1_avg - wcomp1_avg_ctl)*vertmult) # [::-1,:,:]
    array = array1
    lons = uwind.lon 
    lats = uwind.lat
    pres = wwind.pfull

    wwind = area_weighted_avg_4D(wwind,area_array_3D,landmaskxr,'all_sfcs',minlat = minlat, maxlat = maxlat, axis=1)
    uwind = area_weighted_avg_4D(uwind,area_array_3D,landmaskxr,'all_sfcs',minlat = minlat, maxlat = maxlat, axis=1)
    array = area_weighted_avg_4D(array,area_array_3D,landmaskxr,'all_sfcs',minlat = minlat, maxlat = maxlat, axis=1)

    uwind, lons_cyclic = addcyclic(uwind, lons)
    wwind, lons_cyclic = addcyclic(wwind, lons)

    uwind = np.asarray(uwind)
    wwind = np.asarray(wwind)
    uwind,lons_shift = shiftgrid(np.max(lons_cyclic)-shiftby,uwind,lons_cyclic,start=False,
                   cyclic=np.max(lons_cyclic))
    wwind,lons_shift = shiftgrid(np.max(lons_cyclic)-shiftby,wwind,lons_cyclic,start=False,
                   cyclic=np.max(lons_cyclic))  

    array, lons_cyclic = addcyclic(array, lons)
    array = np.asarray(array)
    array, lons_shift = shiftgrid(np.max(lons_cyclic)-shiftby,array,lons_cyclic,
                     start=False,cyclic=np.max(lons_cyclic))


    # landmask,landlons = shiftgrid(np.max(landlons)-80.,landmask,landlons,start=False,cyclic=np.max(landlons))
    # landmask, landlons = addcyclic(landmask, landlons)

    array = xr.DataArray(array,coords=[pres,lons_shift],dims=['pres','lon'])
    uwind = xr.DataArray(uwind,coords=[pres,lons_shift],dims=['pres','lon'])
    wwind = xr.DataArray(wwind,coords=[pres,lons_shift],dims=['pres','lon'])

    array = array.sel(lon=slice(-50.,50.))
    uwind = uwind.sel(lon=slice(-50.,50.))
    wwind = wwind.sel(lon=slice(-50.,50.))

    lons_shift = array.lon
    X, Z = np.meshgrid(lons_shift, pres)




    cset1 = axes[0].contourf(X, Z, array, v, cmap=cmap, extend = 'both')

    Q = axes[0].quiver(X, Z, uwind, wwind, scale = 25, scale_units = 'inches') # if angle isn't set to 'xy', can't invert yaxis on quivers, but angle 'xy' doesn't plot quivers of (u,u) in 45 degree angle! angle 'uv' which is the default does and that's what I want
    # in order for the vectors to all be the same length on all panels and quiverkey to apply to all of them, set scale and scale_units 
    axes[0].set_title('(a) 6$^{\circ}$ 0%cond', fontsize = med)

    # Africa Only 

    uwind = (ucomp2_avg - ucomp1_avg_ctl) # [::-1,:,:]
    wwind = ((wcomp2_avg - wcomp1_avg_ctl)*vertmult) # [::-1,:,:]
    array = array2
    lons = uwind.lon
    lats = uwind.lat
    pres = wwind.pfull

    wwind = area_weighted_avg_4D(wwind,area_array_3D,landmaskxr,'all_sfcs',minlat = minlat, maxlat = maxlat, axis=1)
    uwind = area_weighted_avg_4D(uwind,area_array_3D,landmaskxr,'all_sfcs',minlat = minlat, maxlat = maxlat, axis=1)
    array = area_weighted_avg_4D(array,area_array_3D,landmaskxr,'all_sfcs',minlat = minlat, maxlat = maxlat, axis=1)

    uwind, lons_cyclic = addcyclic(uwind, lons)
    wwind, lons_cyclic = addcyclic(wwind, lons)

    uwind = np.asarray(uwind)
    wwind = np.asarray(wwind)
    uwind,lons_shift = shiftgrid(np.max(lons_cyclic)-shiftby,uwind,lons_cyclic,start=False,
                   cyclic=np.max(lons_cyclic))
    wwind,lons_shift = shiftgrid(np.max(lons_cyclic)-shiftby,wwind,lons_cyclic,start=False,
                   cyclic=np.max(lons_cyclic))  

    array, lons_cyclic = addcyclic(array, lons)
    array = np.asarray(array)
    array, lons_shift = shiftgrid(np.max(lons_cyclic)-shiftby,array,lons_cyclic,
                     start=False,cyclic=np.max(lons_cyclic))


    # landmask,landlons = shiftgrid(np.max(landlons)-80.,landmask,landlons,start=False,cyclic=np.max(landlons))
    # landmask, landlons = addcyclic(landmask, landlons)

    array = xr.DataArray(array,coords=[pres,lons_shift],dims=['pres','lon'])
    uwind = xr.DataArray(uwind,coords=[pres,lons_shift],dims=['pres','lon'])
    wwind = xr.DataArray(wwind,coords=[pres,lons_shift],dims=['pres','lon'])

    array = array.sel(lon=slice(-50.,50.))
    uwind = uwind.sel(lon=slice(-50.,50.))
    wwind = wwind.sel(lon=slice(-50.,50.))

    lons_shift = array.lon
    X, Z = np.meshgrid(lons_shift, pres)


    cset1 = axes[1].contourf(X, Z, array, v, cmap=cmap, extend = 'both')


    Q = axes[1].quiver(X, Z, uwind, wwind, scale = 25, scale_units = 'inches') # if angle isn't set to 'xy', can't invert yaxis on quivers, but angle 'xy' doesn't plot quivers of (u,u) in 45 degree angle! angle 'uv' which is the default does and that's what I want
    #qk = axes[1].quiverkey(Q, 0.9, 0.9, veclen, str(veclen)+r'$\frac{'+units_numerator+'}{'+units_denom+'}$', labelpos='E', coordinates='figure')
    axes[1].set_title('(b) 6$^{\circ}$ bucket', fontsize = med)

    # Two continents - America 

    uwind = ((ucomp1_avg - ucomp1_avg_ctl) - (ucomp2_avg - ucomp1_avg_ctl)) # [::-1,:,:]
    wwind = (((wcomp1_avg - wcomp1_avg_ctl) - (wcomp2_avg - wcomp1_avg_ctl))*vertmult) # [::-1,:,:]
    array = array1 - array2 
    lons = uwind.lon
    lats = uwind.lat
    pres = wwind.pfull

    wwind = area_weighted_avg_4D(wwind,area_array_3D,landmaskxr,'all_sfcs',minlat = minlat, maxlat = maxlat, axis=1)
    uwind = area_weighted_avg_4D(uwind,area_array_3D,landmaskxr,'all_sfcs',minlat = minlat, maxlat = maxlat, axis=1)
    array = area_weighted_avg_4D(array,area_array_3D,landmaskxr,'all_sfcs',minlat = minlat, maxlat = maxlat, axis=1)

    uwind, lons_cyclic = addcyclic(uwind, lons)
    wwind, lons_cyclic = addcyclic(wwind, lons)

    uwind = np.asarray(uwind)
    wwind = np.asarray(wwind)
    uwind,lons_shift = shiftgrid(np.max(lons_cyclic)-shiftby,uwind,lons_cyclic,start=False,
                   cyclic=np.max(lons_cyclic))
    wwind,lons_shift = shiftgrid(np.max(lons_cyclic)-shiftby,wwind,lons_cyclic,start=False,
                   cyclic=np.max(lons_cyclic))  

    array, lons_cyclic = addcyclic(array, lons)
    array = np.asarray(array)
    array, lons_shift = shiftgrid(np.max(lons_cyclic)-shiftby,array,lons_cyclic,
                     start=False,cyclic=np.max(lons_cyclic))


    # landmask,landlons = shiftgrid(np.max(landlons)-80.,landmask,landlons,start=False,cyclic=np.max(landlons))
    # landmask, landlons = addcyclic(landmask, landlons)

    array = xr.DataArray(array,coords=[pres,lons_shift],dims=['pres','lon'])
    uwind = xr.DataArray(uwind,coords=[pres,lons_shift],dims=['pres','lon'])
    wwind = xr.DataArray(wwind,coords=[pres,lons_shift],dims=['pres','lon'])

    array = array.sel(lon=slice(-50.,50.))
    uwind = uwind.sel(lon=slice(-50.,50.))
    wwind = wwind.sel(lon=slice(-50.,50.))

    lons_shift = array.lon
    X, Z = np.meshgrid(lons_shift, pres)




    cset1 = axes[2].contourf(X, Z, array, v, cmap=cmap, extend = 'both')

    Q = axes[2].quiver(X, Z, uwind, wwind, scale = 25, scale_units = 'inches') # if angle isn't set to 'xy', can't invert yaxis on quivers, but angle 'xy' doesn't plot quivers of (u,u) in 45 degree angle! angle 'uv' which is the default does and that's what I want
#    qk = axes[2].quiverkey(Q, 0.83, 0.87, veclen, str(veclen)+r'$\frac{'+units_numerator+'}{'+units_denom+'}$', labelpos='E', coordinates='figure', fontproperties={'size': med})
    axes[2].set_title('(c) 6$^{\circ}$ 0% - bucket', fontsize = med)

    fig.gca().invert_yaxis()

    cbar = fig.colorbar(cset1,ax=axes)
    cbar.ax.tick_params(labelsize=small)
    cbar.set_label(cbarlabel, size = med)

    axes[0].set_ylabel('Pressure (hPa)', fontsize = med)
    axes[0].tick_params(labelsize = small)
    axes[1].tick_params(labelsize = small)
    axes[2].tick_params(labelsize = small)


    fig.text(0.4, 0.03, 'Longitude ($^{\circ}$E)', va='center', rotation='horizontal', fontsize = lge)


    fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/'+dire+'/w'+str(vertmult)+'_'+plotname+'_'+str(minlat)+'N-'+str(maxlat)+'N_lgefonts_nonshift_paper.eps', dpi=600)
    fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/'+dire+'/w'+str(vertmult)+'_'+plotname+'_'+str(minlat)+'N-'+str(maxlat)+'N_lgefonts_nonshift.svg', bbox_inches='tight')
    fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/'+dire+'/w'+str(vertmult)+'_'+plotname+'_'+str(minlat)+'N-'+str(maxlat)+'N_lgefonts_nonshift.pdf', dpi=400)



ctl_model = 'isca'
if (ctl_model == 'Isca') or (ctl_model == 'isca'): 
    control_model = 'Isca_DATA'
elif (ctl_model == 'gfdl') or (ctl_model == 'GFDL'):
    control_model = 'GFDL_DATA'

HPC = 'yes'
control_dir = 'narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_commit7bb4387'
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    control_dir= control_model + '/ISCA_HPC/' + control_dir
else: 
    control_dir= control_model + '/' + control_dir

#print control_dir
ctl_runmin=121
ctl_runmax=481

model = 'isca'
if (model == 'Isca') or (model == 'isca'): 
    model_data = 'Isca_DATA'
    output_dir1 = 'Isca'
elif (model == 'gfdl') or (model == 'GFDL'):
    model_data = 'GFDL_DATA'
    output_dir1 = ''

HPC = 'yes'
testdir_in1= 'narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_corrected_vegpref0_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387'
dire = testdir_in1
runmin=120
runmax=480
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    exp1_name = 'ISCA_HPC_'+testdir_in1
    testdir = model_data + '/ISCA_HPC/' + testdir_in1
    testdir_in1 = '/ISCA_HPC/' + testdir_in1
else: 
    exp1_name = testdir_in1
    testdir = model_data + '/' + testdir_in1

land = 'narrow_three'
landfile=Dataset(os.path.join(GFDL_BASE,'input/'+land+'/land.nc'),mode='r')

landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]
# for specified lats
landmaskxrSA=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it


area_array, dx, dy = ca.cell_area_all(t_res=42,base_dir='/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/') # added _all because then dx and dy are also returned 
area_array = xr.DataArray(area_array) # returned in units of m bzw m^2, because radius in cell_area.py is given in metres

area_array_3D = np.expand_dims(area_array, axis=0)
area_array_3D = np.repeat(area_array_3D, 40, axis = 0) # to make area_array 3D (pressure, lat, lon)

data = xr.open_mfdataset('/scratch/mp586/'+control_dir+'/*/atmos_monthly_interp.nc')
data = data.mean('time')

omega1_avg_ctl = data.omega
rh1_avg_ctl = data.rh
sphum1_avg_ctl = data.sphum
ucomp1_avg_ctl = data.ucomp
temp1_avg_ctl = data.temp

data = xr.open_mfdataset('/scratch/mp586/'+testdir+'/*/atmos_monthly_interp.nc')
data = data.mean('time')

omega1_avg = data.omega
rh1_avg = data.rh
sphum1_avg = data.sphum
ucomp1_avg = data.ucomp
temp1_avg = data.temp

################ read in data from exp 2 ###############################

model = 'isca'
if (model == 'Isca') or (model == 'isca'): 
    model_data = 'Isca_DATA'
    output_dir1 = 'Isca'
elif (model == 'gfdl') or (model == 'GFDL'):
    model_data = 'GFDL_DATA'
    output_dir1 = ''

HPC = 'yes'
testdir_in1= 'narrow_three_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361_commit7bb4387'
runmin=120
runmax=480
if (HPC == 'Yes') or (HPC == 'yes') or (HPC == 'y'):
    exp1_name = 'ISCA_HPC/'+testdir_in1
    testdir = model_data + '/ISCA_HPC/' + testdir_in1
    testdir_in1 = '/ISCA_HPC/' + testdir_in1
else: 
    exp1_name = testdir_in1
    testdir = model_data + '/' + testdir_in1

land = 'narrow_three'
landfile=Dataset(os.path.join(GFDL_BASE,'input/'+land+'/land.nc'),mode='r')

landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]
# for specified lats
landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

data = xr.open_mfdataset('/scratch/mp586/'+testdir+'/*/atmos_monthly_interp.nc')
data = data.mean('time')

omega2_avg = data.omega
rh2_avg = data.rh
sphum2_avg = data.sphum
ucomp2_avg = data.ucomp
temp2_avg = data.temp


################ read in data from exp 3 ###############################


outdir = 'Isca' + '/' + exp1_name

# vert_horiz_winds(outdir,runmin,runmax,'delta wind',ucomp2_avg_ctl - ucomp2_ctl_zonavg,(omega2_avg_ctl - omega2_ctl_zonavg)*300.,rh2_avg_ctl.sel(lat = slice(-10.,10.)).mean(dim = 'lat'),20.,90.,veclen=5,units_numerator = '', units_denom = '',save = False)


# # don't need to invert vertical axis for RH because I invert yaxis in plot, but that doesnt work for the quivers so need to invert those!
# vert_horiz_winds(outdir,runmin,runmax.'weighted change',(((ucomp2_avg - ucomp2_avg_ctl) - (ucomp1_avg - ucomp1_avg_ctl))/ucomp2_avg.max())[::-1,:,:],(((omega2_avg - omega2_avg_ctl) - (omega1_avg - omega1_avg_ctl))/omega2_avg.max())[::-1,:,:],((rh2_avg - rh2_avg_ctl) - (rh1_avg - rh1_avg_ctl)).sel(lat = slice(-10.,10.)).mean(dim = 'lat'),-10.,10.,veclen=1,units_numerator = 'Pa m', units_denom = 's s',save = False)


# #this needs to produce vectors that are at a 45 degree angle - angle in quiver has to be set to 'uv' (which is the default) or not set at all duh. Wasn't working because I had set angle to 'xy' so that it would invert the yaxis. Instead can just feed the info in reversed omega coords! 
# vert_horiz_winds(outdir,runmin,runmax,'u u',(ucomp2_avg),ucomp2_avg,((rh2_avg).sel(lat = slice(-10.,10.))).mean(dim = 'lat'),0,80.,veclen=10,units_numerator = 'Pa m', units_denom = 's s',save = False)


g = 9.81
Rspec = 287.058
pfull = data.pfull * 100 # convert from hPa to Pa
pfull = np.expand_dims(pfull, axis = 1)
pfull = np.expand_dims(pfull, axis = 2)
pfull = np.repeat(pfull, 64, axis=1)
pfull = np.repeat(pfull, 128, axis = 2)
pres_lev = data.pfull
pfull = xr.DataArray(pfull, coords = [pres_lev, data.lat, data.lon], dims = ['pfull','lat','lon'])
wcomp1_avg_ctl = - (omega1_avg_ctl * temp1_avg_ctl * Rspec)/(pfull * g)
wcomp1_avg = - (omega1_avg * temp1_avg * Rspec)/(pfull * g)
wcomp2_avg = - (omega2_avg * temp2_avg * Rspec)/(pfull * g)

#conversion following https://www.ncl.ucar.edu/Document/Functions/Contributed/omega_to_w.shtml


#vert_horiz_winds(outdir,runmin,runmax,'u w*80',(((ucomp2_avg - ucomp2_avg_ctl) - (ucomp1_avg - ucomp1_avg_ctl)))[::-1,:,:],(((wcomp2_avg - wcomp2_avg_ctl) - (wcomp1_avg - wcomp1_avg_ctl))*80.)[::-1,:,:],((rh2_avg - rh2_avg_ctl) - (rh1_avg - rh1_avg_ctl)).sel(lat = slice(-10.,10.)).mean(dim = 'lat'),-10.,10.,veclen=5,units_numerator = 'm', units_denom = 's',save = False)


# Panel plot with 4 cases: America only, Africa only, Two continents, and Two continents minus America - climate change for all 



quivers_2cases(runmin, runmax, 'ucorr_interp_quivers_2cases_CV0vbucket_deltarh_fct_latweights_', dire, landmaskxr, '$\Delta$ r (%)', ucomp1_avg, ucomp1_avg_ctl, wcomp1_avg, wcomp1_avg_ctl, ucomp2_avg, ucomp1_avg_ctl, wcomp2_avg, wcomp1_avg_ctl, (rh1_avg - rh1_avg_ctl), (rh2_avg - rh1_avg_ctl), minval=-10., maxval=10., vertmult=3000, minlat=-10., maxlat=10.)
quivers_2cases(runmin, runmax, 'ucorr_interp_quivers_2cases_CV0vbucket_deltasphum_fct_latweights_', dire, landmaskxr, '$\Delta$ q (kg/kg)', ucomp1_avg, ucomp1_avg_ctl, wcomp1_avg, wcomp1_avg_ctl, ucomp2_avg, ucomp1_avg_ctl, wcomp2_avg, wcomp1_avg_ctl, (sphum1_avg - sphum1_avg_ctl), (sphum2_avg - sphum1_avg_ctl), minval=-.003, maxval=.003, vertmult=3000, minlat=-10., maxlat=10.)
quivers_2cases(runmin, runmax, 'ucorr_interp_quivers_2cases_CV0vbucket_deltatemp_fct_latweights_', dire, landmaskxr, '$\Delta$ T (K)', ucomp1_avg, ucomp1_avg_ctl, wcomp1_avg, wcomp1_avg_ctl, ucomp2_avg, ucomp1_avg_ctl, wcomp2_avg, wcomp1_avg_ctl, (temp1_avg - temp1_avg_ctl), (temp2_avg - temp1_avg_ctl), minval=-5., maxval=5., vertmult=3000, minlat=-10., maxlat=10., cmap = 'RdBu_r')





