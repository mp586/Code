#NB Don't do stats on an array that contains nans
#Mean over array with nans ignores the nan values so don't need to take them out first.

from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
import xarray as xr
import pandas as pd
import os
from scipy.stats import linregress as linreg
from scipy import odr # regression routines for parameter fits.

import sys
sys.path.insert(0, '/scratch/mp586/Code/PYCODES')
from plotting_routines_kav7 import *
import stats as st

GFDL_BASE = os.environ['GFDL_BASE']
sys.path.insert(0, os.path.join(GFDL_BASE,'src/extra/python/scripts'))
import cell_area as ca


# Get model data 

ctl_model = 'isca'
if (ctl_model == 'Isca') or (ctl_model == 'isca'): 
    control_model = 'Isca_DATA'
elif (ctl_model == 'gfdl') or (ctl_model == 'GFDL'):
    control_model = 'GFDL_DATA'


control_dir= control_model + '/' + 'ISCA_HPC/two_continents_newbucket_fixedSSTs_from_realworld_zonallysymm'
print control_dir
ctl_runmin=121  # Should be a January month for seasonal variables to be correct
ctl_runmax=481

model = 'isca'
if (model == 'Isca') or (model == 'isca'): 
    model_data = 'Isca_DATA'
    output_dir = 'Isca'
elif (model == 'gfdl') or (model == 'GFDL'):
    model_data = 'GFDL_DATA'
    output_dir = ''

testdir= 'ISCA_HPC/two_continents_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361'
runmin=120  # Should be a January month for seasonal variables to be correct
runmax=480


outdir = output_dir + '/' + testdir
testdir = model_data + '/' + testdir

landfile=Dataset(os.path.join(GFDL_BASE,'input/two_continents/land.nc'),mode='r') # landmask both continents

landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]
# for specified lats
landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

landfileWC=Dataset(os.path.join(GFDL_BASE,'input/square_South_America/land.nc'),mode='r')

landmaskWC=landfileWC.variables['land_mask'][:] # land mask west continent
landmaskWCxr=xr.DataArray(landmaskWC,coords=[landlats,landlons],dims=['lat','lon']) 


landfileEC=Dataset(os.path.join(GFDL_BASE,'input/square_Africa/land.nc'),mode='r')

landmaskEC=landfileEC.variables['land_mask'][:] # land mask east continent
landmaskECxr=xr.DataArray(landmaskEC,coords=[landlats,landlons],dims=['lat','lon']) 


area_array = ca.cell_area(t_res=42,base_dir='/scratch/mp586/GFDL_BASE/GFDL_FORK/GFDLmoistModel/')
area_array = xr.DataArray(area_array)


################################################

# Variable definitions 

dthresh = 1. # mm / day
dd = 0.07 # Held and Soden alpha = our d

[tsper,tspermean,tsper_seasonal_avg,tsper_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'t_surf','K')
[evper,evpermean,evper_seasonal_avg,evper_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'flux_lhe','W/m^2',factor = 1./28.) # latent heat flux in mm/d
[prper,prpermean,prper_seasonal_avg,prper_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'precipitation','kg/m2s', factor=86400)
[rhper,rhpermean,rhper_seasonal_avg,rhper_month_avg,time]=seasonal_surface_variable(testdir,model,runmin,runmax,'rh','%',level=39, factor = 1/100.)


[tscon,tsconmean,tscon_seasonal_avg,tscon_month_avg,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'t_surf','K')
[evcon,evconmean,evcon_seasonal_avg,evcon_month_avg,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'flux_lhe','W/m^2',factor = 1./28.) # latent heat flux mm/d
[prcon,prconmean,prcon_seasonal_avg,prcon_month_avg,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'precipitation','kg/m2s', factor=86400)
[rhcon,rhconmean,rhcon_seasonal_avg,rhcon_month_avg,time]=seasonal_surface_variable(control_dir,ctl_model,ctl_runmin,ctl_runmax,'rh','%',level=39, factor = 1/100.)

# annual means 

pmeconm = prconmean - evconmean
dts = tspermean - tsconmean # d Tsurf
dpr = (prpermean - prconmean) # dP
dprw = (dpr.where(np.asarray(prconmean) > dthresh)) # select only those regions where P > threshold in the control
dprtw = dprw.sel(lat=slice(-30.,30.)).mean(dim = ['lat','lon']) # tropical mean dP in regions where Pcontrol > threshold

prpert = prpermean.sel(lat=slice(-30.,30.)).mean(dim = ['lat','lon'])
prcont = prconmean.sel(lat=slice(-30.,30.)).mean(dim = ['lat','lon'])

dpmehs = pmeconm * dd * dts

prrcptw = prpermean.where(np.asarray(prconmean) > dthresh).sel(lat=slice(-30.,30.)).mean(dim = ['lat','lon'])
# tropical mean perturbed P (don't have rcp data!) only where threshold condition is met 
dtstw = dts.where(np.asarray(prconmean) > dthresh).sel(lat=slice(-30.,30.)).mean(dim = ['lat','lon'])

bbw = dprtw/(prrcptw * dtstw)

curly_r = prconmean.where(np.asarray(prconmean) > dthresh).where(landmask == 1.).sel(lat=slice(-30.,30)) / prcont


# Orthogonal distance regression linear models for f and c

def f(B, x): # We will use the ordinary linear model.
    '''Linear function y = m*x + b'''
    return B[0]*x + B[1]


odrfit = odr.Model(f)

nlps = np.sum(landmaskxr.sel(lat=slice(-30.,30.))) # number of landpoints in the tropics. landmask has 1.s for landpoitns and 0. for ocean

rhconl = rhconmean.where(landmask == 1.).where(np.asarray(prconmean) > dthresh).sel(lat=slice(-30.,30.))
rhconl = (np.asarray(rhconl)).flatten()
rhconl = rhconl[~np.isnan(rhconl)]
curly_r = (np.asarray(curly_r)).flatten()
curly_r = curly_r[~np.isnan(curly_r)]

prconl = prconmean.where(np.asarray(prconmean) > dthresh).where(landmask == 1.).sel(lat=slice(-30.,30.))
prconl = np.asarray(prconl).flatten()
prconl = prconl[~np.isnan(prconl)]


evconl = evconmean.where(np.asarray(prconmean) > dthresh).where(landmask == 1.).sel(lat=slice(-30.,30.))
evconl = np.asarray(evconl).flatten()
evconl = evconl[~np.isnan(evconl)]

# First f. (Known here as ff.) f is d curly_r / drhconl. rhconl (control relative humidity) and curly_r have one value per land point above the wetness threshold. Should get something between about 3 and 5, I think...

xstd = np.repeat(rhconl.std(),nlps) # We don't have a lot of information about 
ystd = np.repeat(curly_r.std(),nlps) # errors, so we'll assume that the ratio of standard deviations of the x and y data is the same as the ratio of the errors in x and y. Let me know if that doesn't make any sense...
ffdata = odr.RealData(rhconl,curly_r,sx=xstd,sy=ystd)
fffit = odr.ODR(ffdata, odrfit, beta0=[1., 2.])
ffout = fffit.run()
ff = ffout.beta[0]



# Then c. prconl is control precipitation for each >1 mm/day land point. evconl is evaporation for same.

xstd = np.repeat(prconl.std(),nlps)
ystd = np.repeat(evconl.std(),nlps)
ccdata = odr.RealData(prconl,evconl,sx=xstd,sy=ystd)
ccfit = odr.ODR(ccdata, odrfit, beta0=[1., 2.])
ccout = ccfit.run()
cc = ccout.beta[0]

# Equation 1. prconm is gridded control precipitation. dtstw is tropical mean (i.e it includes ocean!) wet region temperature increase. 
# drhlt is wet region tropical land mean rh change. prcontw is tropical mean (including ocean again) wet region control precipitation

prconm = prconmean
dtstw = (tspermean - tsconmean).where(np.asarray(prconmean) > dthresh).sel(lat=slice(-30.,30.)).mean(dim = ['lat','lon'])
drhlt = (rhpermean - rhconmean).where(np.asarray(prconmean) > dthresh).where(landmask == 1.).sel(lat=slice(-30.,30.)).mean(dim = ['lat','lon'])
prcontw = prconmean.where(np.asarray(prconmean) > dthresh).sel(lat=slice(-30.,30.)).mean(dim = ['lat','lon'])
drh = (rhpermean - rhconmean)

dprlp = ((bbw * prconm * dtstw) + (ff * prcontw * (drh - drhlt))) / (1. - (bbw * dtstw)) 

# Byrne / Chadwick parameterisation of moisture
# ----------------------------------------------

# We will use "dd" as above, but note that both Byrne and Chadwick suggest something more sophisticated.

# Super simple version from bottom of page 2 of iscaprecipshort that does scaling relative to land only:


rhconm = rhconmean
dtslt = (tspermean - tsconmean).where(np.asarray(prconmean) > dthresh).where(landmask == 1.).sel(lat=slice(-30.,30.)).mean(dim = ['lat','lon'])
dts = tspermean - tsconmean

drhbcs = rhconm * dd * (dtslt - dts) # dtslt is tropical wet region land mean temperature change. rhconm is gridded control RH. (Obviously time mean.) 

# Put it in precipitation estimate (equation 3).

dprlpbcm = ((bbw * prconm * dtstw) + (ff * prcontw * drhbcs)) / (1. - (bbw * dtstw))

# Precipitation minus evaporation. Simple:

# dpmelp = (1.-cc)*dprlp # Equation 1 version. (I.e. we know RH.)
# dpmelpbc = (1.-cc) * dprlpbc # with Byrne-Chadwick stuck in.

############ plotting ##############

small = 10 #largefonts 14 # smallfonts 10 # medfonts = 14
med = 18 #largefonts 18 # smallfonts 14 # medfonts = 16
lge = 22 #largefonts 22 # smallfonts 18 # medfonts = 20

v = np.linspace(-2.,2.,21)
nmb_contours = [-2.,1.,2.]


lats=prconmean.lat
lons=prconmean.lon

landlats = np.asarray(landmaskxr.lat)
landlons = np.asarray(landmaskxr.lon)

landmask = np.asarray(landmaskxr)


# Read landmask

# Add rectangles
#    landmask,landlons = shiftgrid(np.max(landlons)-100.,landmask,landlons,start=True,cyclic=np.max(landlons)) # this works when the array shift is commented....
landmask,landlons = shiftgrid(np.max(landlons)-180.,landmask,landlons,start=False,cyclic=np.max(landlons))

landmask, lons_cyclic = addcyclic(landmask, landlons)


fig, axes = plt.subplots(1, 3, figsize = (30,10))

pltnames = ['Pred. 1.','Pred. 2','Actual']
arrays = [dprlpbcm, dprlp, (prpermean - prconmean)]

for i in range(3):
	axes[i].set_title(pltnames[i], size = med)
	#fig = plt.figure()

	m = Basemap(projection='kav7',lon_0=0.,resolution='c', ax = axes[i])

	array = xr.DataArray(arrays[i],coords=[lats,lons],dims=['lat','lon'])

	array = np.asarray(array)
	array, lons_cyclic = addcyclic(array, lons)
	array,lons_cyclic = shiftgrid(np.max(lons_cyclic)-180.,array,lons_cyclic,start=False,cyclic=np.max(lons_cyclic))

	array = xr.DataArray(array,coords=[lats,lons_cyclic],dims=['lat','lon'])

	m.drawparallels(np.arange(-90.,99.,30.),labels=[1,0,0,0], fontsize=small)
	m.drawmeridians(np.arange(-180.,180.,60.),labels=[0,0,0,1], fontsize=small)

	lon, lat = np.meshgrid(lons_cyclic, lats)
	xi, yi = m(lon, lat)

	cs = m.contourf(xi,yi,array, v, cmap='BrBG', extend = 'both')


	if np.any(landmask != 0.):
	    m.contour(xi,yi,landmask, 1)

# Add Colorbar
cbar = fig.colorbar(cs, orientation = 'vertical', ax = axes, shrink = 0.4, fraction = 0.1) 
cbar.set_label('mm/d', size=med)
cbar.ax.tick_params(labelsize=med)

fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/two_continents_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361/bc_hugo_pred_120-480.png', bbox_inches='tight', dpi=100)
fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/two_continents_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361/bc_hugo_pred_120-480.svg', bbox_inches='tight', dpi=100)



#scatter plotts 


landfileWC=Dataset(os.path.join(GFDL_BASE,'input/square_South_America/land.nc'),mode='r')

landmaskWC=landfileWC.variables['land_mask'][:] # land mask west continent
landmaskWCxr=xr.DataArray(landmaskWC,coords=[landlats,landlons],dims=['lat','lon']) 


landfileEC=Dataset(os.path.join(GFDL_BASE,'input/square_Africa/land.nc'),mode='r')

landmaskEC=landfileEC.variables['land_mask'][:] # land mask east continent
landmaskECxr=xr.DataArray(landmaskEC,coords=[landlats,landlons],dims=['lat','lon']) 

dppred1 = dprlp
dppred2 = dprlpbcm

dppred1EC = np.asarray(dppred1.where(landmaskEC == 1.)).flatten()
dppred1EC = dppred1EC[~np.isnan(dppred1EC)]
dppred1WC = np.asarray(dppred1.where(landmaskWC == 1.)).flatten()
dppred1WC = dppred1WC[~np.isnan(dppred1WC)]

dppred2EC = np.asarray(dppred2.where(landmaskEC == 1.)).flatten()
dppred2EC = dppred2EC[~np.isnan(dppred2EC)]
dppred2WC = np.asarray(dppred2.where(landmaskWC == 1.)).flatten()
dppred2WC = dppred2WC[~np.isnan(dppred2WC)]

dpactualEC = np.asarray((prpermean - prconmean).where(landmaskEC == 1.)).flatten()
dpactualEC = dpactualEC[~np.isnan(dpactualEC)]
dpactualWC = np.asarray((prpermean - prconmean).where(landmaskWC == 1.)).flatten()
dpactualWC = dpactualWC[~np.isnan(dpactualWC)]






fig, ax = plt.subplots(1,2,sharex = True, sharey = True, figsize=(25,15))

ax[1].plot(dpactualWC,dppred1WC,'b.',label='West C') # West Island
ax[1].plot(dpactualEC, dppred1EC,'r.',label = 'East C') # East Island
ax[1].set_xlabel('Actual $\Delta P$ (mm/d)')
ax[1].set_ylabel('Predicted $\Delta P$ (mm/day)')
ax[1].legend()
ax[1].set_title('Prediction 2')

ax[0].plot(dpactualWC,dppred2WC,'b.',label='West C') # West Island
ax[0].plot(dpactualEC, dppred2EC,'r.',label = 'East C') # East Island
ax[0].set_xlabel('Actual $\Delta P$ (mm/d)')
ax[0].set_ylabel('Predicted $\Delta P$ (mm/day)')
ax[0].legend()
ax[0].set_title('Prediction 1')


#############################################################################################################################
### LINEAR REGRESIONS ###

[k,dy,r,p,stderr] = linreg(dpactualWC,dppred1WC) # aa = 8.4, dq = -32
x1 = np.linspace(np.min(dpactualWC),np.max(dpactualWC),500)
y = k*x1 + dy
ax[1].plot(x1,y,'b-')
ax[1].annotate('West C: k = '+str("%.2f" % k)+', dy = '+str("%.2f" % dy)+', r_val = '+str("%.2f" % r), xy=(0.05,0.1), xycoords='axes fraction', fontsize = med)

[k,dy,r,p,stderr] = linreg(dpactualEC,dppred1EC) # aa = 8.4, dq = -32
x1 = np.linspace(np.min(dpactualEC),np.max(dpactualEC),500)
y = k*x1 + dy
ax[1].plot(x1,y,'r-')
ax[1].annotate('East C: k = '+str("%.2f" % k)+', dy = '+str("%.2f" % dy)+', r_val = '+str("%.2f" % r), xy=(0.05,0.05), xycoords='axes fraction', fontsize = med)



[k,dy,r,p,stderr] = linreg(dpactualWC,dppred2WC) # aa = 8.4, dq = -32
x1 = np.linspace(np.min(dpactualWC),np.max(dpactualWC),500)
y = k*x1 + dy
ax[0].plot(x1,y,'b-')
ax[0].annotate('West C: k = '+str("%.2f" % k)+', dy = '+str("%.2f" % dy)+', r_val = '+str("%.2f" % r), xy=(0.05,0.1), xycoords='axes fraction', fontsize = med)

[k,dy,r,p,stderr] = linreg(dpactualEC,dppred2EC) # aa = 8.4, dq = -32
x1 = np.linspace(np.min(dpactualEC),np.max(dpactualEC),500)
y = k*x1 + dy
ax[0].plot(x1,y,'r-')
ax[0].annotate('East C: k = '+str("%.2f" % k)+', dy = '+str("%.2f" % dy)+', r_val = '+str("%.2f" % r), xy=(0.05,0.05), xycoords='axes fraction', fontsize = med)

fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/two_continents_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361/bc_hugo_pred_120-480_scatterplot.png', bbox_inches='tight', dpi=100)



# redefine landmask because overwritten above for cyclic -- map plots 
landmask=landfile.variables['land_mask'][:]
landlats=landfile.variables['lat'][:]
landlons=landfile.variables['lon'][:]
# for specified lats
landmaskxr=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it


dpactuall = np.asarray((prpermean - prconmean).where(np.asarray(prconmean) > dthresh).where(landmask == 1.)).flatten()
dtsl = np.asarray((tspermean - tsconmean).where(np.asarray(prconmean) > dthresh).where(landmask == 1.)).flatten()
drhl = np.asarray((rhpermean - rhconmean).where(np.asarray(prconmean) > dthresh).where(landmask == 1.)).flatten()
devl = np.asarray((evpermean - evconmean).where(np.asarray(prconmean) > dthresh).where(landmask == 1.)).flatten()

dpactuall = dpactuall[~np.isnan(dpactuall)]
dtsl = dtsl[~np.isnan(dtsl)]
drhl = drhl[~np.isnan(drhl)]
devl = devl[~np.isnan(devl)]




fig, ax = plt.subplots(1,3, sharex = True, figsize=(30,10))

ax[0].plot(dpactuall,dtsl,'r.')
ax[0].set_xlabel('Actual $\Delta P$ (mm/d)', fontsize = med)
ax[0].set_ylabel('Actual $\Delta T_S$ (K)', fontsize = med)
ax[0].legend()
ax[0].set_title('$\Delta T_S$ vs $\Delta P$', fontsize = med)

ax[1].plot(dpactuall,devl,'g.')
ax[1].set_xlabel('Actual $\Delta P$ (mm/d)', fontsize = med)
ax[1].set_ylabel('Actual $\Delta E$ (mm/day)', fontsize = med)
ax[1].legend()
ax[1].set_title('$\Delta E$ vs $\Delta P$', fontsize = med)

ax[2].plot(dpactuall,drhl,'k.')
ax[2].set_xlabel('Actual $\Delta P$ (mm/d)', fontsize = med)
ax[2].set_ylabel('Actual $\Delta r$ (mm/day)', fontsize = med)
ax[2].legend()
ax[2].set_title('$\Delta r_S$ vs $\Delta P$', fontsize = med)



#############################################################################################################################
### LINEAR REGRESIONS ###

[k,dy,r,p,stderr] = linreg(dpactuall,dtsl) # aa = 8.4, dq = -32
x1 = np.linspace(np.min(dpactuall),np.max(dpactuall),500)
y = k*x1 + dy
ax[0].plot(x1,y,'r-')
ax[0].annotate('r = '+str("%.2f" % r), xy=(0.05,0.03), xycoords='axes fraction', fontsize = med)
# +', k = '+str("%.2f" % k)+', dy = '+str("%.2f" % dy) 

[k,dy,r,p,stderr] = linreg(dpactuall,devl) # aa = 8.4, dq = -32
x1 = np.linspace(np.min(dpactuall),np.max(dpactuall),500)
y = k*x1 + dy
ax[1].plot(x1,y,'g-')
ax[1].annotate('r = '+str("%.2f" % r), xy=(0.65,0.03), xycoords='axes fraction', fontsize = med)


[k,dy,r,p,stderr] = linreg(dpactuall,drhl) # aa = 8.4, dq = -32
x1 = np.linspace(np.min(dpactuall),np.max(dpactuall),500)
y = k*x1 + dy
ax[2].plot(x1,y,'k-')
ax[2].annotate('r = '+str("%.2f" % r), xy=(0.65,0.03), xycoords='axes fraction', fontsize = med)

fig.savefig('/scratch/mp586/Code/Graphics/Isca/ISCA_HPC/two_continents_newbucket_fixedSSTs_from_realworld_zonallysymm_plus_uniform_warming_and_2xCO2_spinup_361/dP_vs_dTs_dE_dRH_120-480_.png', bbox_inches='tight', dpi=100)

