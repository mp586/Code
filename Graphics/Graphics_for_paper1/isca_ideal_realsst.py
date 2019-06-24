#!/usr/bin/env python

# Analysis of Marianne's two continent runs that use "real SSTs" zonally symettrised and that we'd like to use in her paper, using orthogonal distance regression.

# No energy budget stuff. Just the basics of precipitation using the equations that work for CMIP5.


# Must run on baikonur. 

from numpy import *
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import odr
import latweight # My quick latitude weighting routine

# Directories. 

datdir = '/scratch/fhl202/marianne_idealised/'

# Choices. First which land mask?

lmaskdir = 'two_continents_land.nc' # Idealised land mask for two continents.

# Second. Which data?

con_file = Dataset(datdir+'two_continent_bucket_real_SSTs_control.nc') # Control real SST
#con_file = Dataset(datdir+'two_continent_bucket_control.nc') # Control normal

per_file = Dataset(datdir+'two_continent_bucket_real_SSTS_warm+2xco2.nc') # 2xCO2 real SST
#per_file = Dataset(datdir+'two_continent_bucket_2xco2.nc') # 2xCO2 normal


# Fetch netcdf data

lon = con_file.variables['lon']
lat = con_file.variables['lat']
tcon = con_file.variables['t_surf']
prcon = con_file.variables['precipitation']
pcon = con_file.variables['pfull']
speccon = con_file.variables['sphum']
relcon = con_file.variables['rh']
swcon = con_file.variables['flux_sw']
lwdcon = con_file.variables['flux_lw']
lhcon = con_file.variables['flux_lhe']
shcon = con_file.variables['flux_t']
#qfluxcon = con_file.variables['flux_oceanq'] # Have to turn qflux off for "0qflux", but it actually isn't used anywhere in the analysis.
ucon = con_file.variables['ucomp']
vcon = con_file.variables['vcomp']
omegacon = con_file.variables['omega']
lwucon = 5.67*10**(-8) * (array(tcon)**4.) # construct upward longwave flux

tper = per_file.variables['t_surf']
prper = per_file.variables['precipitation']
pper = per_file.variables['pfull']
specper = per_file.variables['sphum']
relper = per_file.variables['rh']
swper = per_file.variables['flux_sw']
lwdper = per_file.variables['flux_lw']
lhper = per_file.variables['flux_lhe']
shper = per_file.variables['flux_t']
#qfluxper = per_file.variables['flux_oceanq']
uper = per_file.variables['ucomp']
vper = per_file.variables['vcomp']
omegaper = per_file.variables['omega']
lwuper = 5.67*10**(-8) * (array(tper)**4.) # construct upward longwave flux

lmask_file = Dataset(lmaskdir) # Land mask
lmask = lmask_file.variables['land_mask']

# Temporal means of the last 20 years for the tropics. (30N - 30S). For 3D variables, we'll take the surface only for now.

tconm = tcon[240:,21:43,:].mean(0)
prconm = prcon[240:,21:43,:].mean(0) * 86400. # Make mm / day.
specconm = speccon[240:,-1,21:43,:].mean(0)
relconm = relcon[240:,-1,21:43,:].mean(0) / 100. # Turn into fraction of 1. Easier.
swconm = swcon[240:,21:43,:].mean(0)
lwdconm = lwdcon[240:,21:43,:].mean(0)
lwuconm = lwucon[240:,21:43,:].mean(0)
lhconm = lhcon[240:,21:43,:].mean(0)
shconm = shcon[240:,21:43,:].mean(0)
uconm = ucon[240:,:,21:43,:].mean(0) # Keep all vertical levels for the winds.
vconm = vcon[240:,:,21:43,:].mean(0)
omegaconm = omegacon[240:,:,21:43,:].mean(0)
#qfluxconm = qfluxcon[240:,21:43,:].mean(0)
satconm = specconm / relconm # bad estimate of saturation specific humidity

tperm = tper[240:,21:43,:].mean(0)
prperm = prper[240:,21:43,:].mean(0) * 86400. # Make mm / day.
specperm = specper[240:,-1,21:43,:].mean(0)
relperm = relper[240:,-1,21:43,:].mean(0) / 100.
swperm = swper[240:,21:43,:].mean(0)
lwdperm = lwdper[240:,21:43,:].mean(0)
lwuperm = lwuper[240:,21:43,:].mean(0)
lhperm = lhper[240:,21:43,:].mean(0)
shperm = shper[240:,21:43,:].mean(0)
uperm = uper[240:,:,21:43,:].mean(0) 
vperm = vper[240:,:,21:43,:].mean(0)
omegaperm = omegaper[240:,:,21:43,:].mean(0)
#qfluxperm = qfluxper[240:,21:43,:].mean(0) # THE SAME AS CONTROL.
satperm = specperm / relperm # bad estimate of saturation specific humidity

lmaskt = lmask[21:43,:] # Tropical land mask


# Land and ocean-only data.

tconl = tconm[where(lmaskt==1.)]
prconl = prconm[where(lmaskt==1.)]
specconl = specconm[where(lmaskt==1.)]
relconl = relconm[where(lmaskt==1.)]
swconl = swconm[where(lmaskt==1.)]
lwdconl = lwdconm[where(lmaskt==1.)]
lwuconl = lwuconm[where(lmaskt==1.)]
lhconl = lhconm[where(lmaskt==1.)]
shconl = shconm[where(lmaskt==1.)]
#qfluxconl = qfluxconm[where(lmaskt==1.)]
satconl = satconm[where(lmaskt==1.)]

tperl = tperm[where(lmaskt==1.)]
prperl = prperm[where(lmaskt==1.)]
specperl = specperm[where(lmaskt==1.)]
relperl = relperm[where(lmaskt==1.)]
swperl = swperm[where(lmaskt==1.)]
lwdperl = lwdperm[where(lmaskt==1.)]
lwuperl = lwuperm[where(lmaskt==1.)]
lhperl = lhperm[where(lmaskt==1.)]
shperl = shperm[where(lmaskt==1.)]
#qfluxperl = qfluxperm[where(lmaskt==1.)]
satperl = satperm[where(lmaskt==1.)]

tcono = tconm[where(lmaskt==0.)]
prcono = prconm[where(lmaskt==0.)]
speccono = specconm[where(lmaskt==0.)]
relcono = relconm[where(lmaskt==0.)]
swcono = swconm[where(lmaskt==0.)]
lwdcono = lwdconm[where(lmaskt==0.)]
lwucono = lwuconm[where(lmaskt==0.)]
lhcono = lhconm[where(lmaskt==0.)]
shcono = shconm[where(lmaskt==0.)]
#qfluxcono = qfluxconm[where(lmaskt==0.)]
satcono = satconm[where(lmaskt==0.)]

tpero = tperm[where(lmaskt==0.)]
prpero = prperm[where(lmaskt==0.)]
specpero = specperm[where(lmaskt==0.)]
relpero = relperm[where(lmaskt==0.)]
swpero = swperm[where(lmaskt==0.)]
lwdpero = lwdperm[where(lmaskt==0.)]
lwupero = lwuperm[where(lmaskt==0.)]
lhpero = lhperm[where(lmaskt==0.)]
shpero = shperm[where(lmaskt==0.)]
#qfluxpero = qfluxperm[where(lmaskt==0.)]
satpero = satperm[where(lmaskt==0.)]

nlps = shape(where(lmaskt==1.))[1] # Number of land points



# Weights for latitude weighting:

wgt = latweight.lw(inlat=lat)


# Useful shorthands.

dttrop = (tperm.mean() - tconm.mean()) # tropical mean temp change
dttropo = (tpero.mean() - tcono.mean()) # tropical ocean mean temp change
dttropl = (tperl.mean() - tconl.mean()) # tropical land mean temp change
dprtrop = (prperm.mean() - prconm.mean()) # tropical mean precipitation change
dspectrop = (specperm.mean() - specconm.mean()) # tropical mean specific humidity change
dspectropl = (specperl.mean() - specconl.mean()) # trop land mean spec change
dt = tperm - tconm # temperature change
dpr = prperm - prconm # precipitation change
dev = (lhperm - lhconm) / 28.94 # evaporation (Think I've got the factor right for Wm^-2 => mm / day.)
drel = relperm - relconm # rh change
dspec = specperm - specconm # specific humidity change
dreltropl = (relperl.mean() - relconl.mean()) # tropical land rel change.


# Linear models for f, c and g
# ----------------------------

def f(B, x): # We will use the ordinary linear model.
    '''Linear function y = m*x + b'''
    return B[0]*x + B[1]

odrfit = odr.Model(f)


# b. (Not a regression estimate.)

tropprecipmean = (prperm.mean() - prconm.mean()) / prconm.mean() # mean percentage precip change in the tropics

bb = tropprecipmean / (tperm.mean() - tconm.mean()) # mm/day per K change in precipitation.


# d

dd = 0.07 # We'll assume a 7%/K increase in saturation specific humidity.

# f

curly_r = prconl / prconm.mean() # Function that links local precipitation to tropical mean precipitation change. (Note we want to consider land data only for the individual points.)

xstd = repeat(relconl.std(),nlps)
ystd = repeat(curly_r.std(),nlps)
ffdata = odr.RealData(relconl,curly_r,sx=xstd,sy=ystd)
fffit = odr.ODR(ffdata, odrfit, beta0=[1., 2.])
ffout = fffit.run()

ff = ffout.beta[0]

# g (the relationship between q and RH). (r = gq).

xstd = repeat(specconl.std(),nlps)
ystd = repeat(relconl.std(),nlps)
ggdata = odr.RealData(specconl,relconl,sx=xstd,sy=ystd)
ggfit = odr.ODR(ggdata, odrfit, beta0=[1., 2.])
ggout = ggfit.run()

gg = ggout.beta[0]



# Predictions

dprlp = ((bb * prconm * dttrop) + (ff * prconm.mean() * (drel - dreltropl))) / (1. - (bb * dttrop)) # LP


# Now BC scaling of moisture. But this is cheating as temperature change tells us most of what we need to know.

drelbc = relconm * dd * (dttropl - dt) 

# BC scaling of SH.

alpha = 0.07 # HS fractional increase in SH per K.

dspecbc = (specconm / specconm.mean()) * dspectrop
#dspecbc = (specconm / specconm.mean()) * alpha * dttrop * specconm # USE H&S STYLE FORMULATION.
dspecbctropl = dspecbc[where(lmaskt==1.)].mean()

# If we assume gg relationship from control, what change in RH do we get with BC?

drelbcgg = (dspecbc * gg) - (dspecbctropl * gg)

drelgg = (dspec * gg) - (dspectropl * gg) # Also using actual humidity. How well does this do?


# Put it in precipitation estimate (equation 3).

dprlpbc = ((bb * prconm * dttrop) + (ff * prconm.mean() * drelbc)) / (1. - (bb * dttrop))

dspecbc = (specconm / specconm.mean()) * dspectrop


# BC scaling of SH.

#dspecbc = (specconm / specconm.mean()) * dspectrop
#dspecbctropl = dspecbc[where(lmaskt==1.)].mean()

# If we assume gg relationship from control, what change in RH do we get with BC?

#drelbcgg = (dspecbc * gg) - (dspectropl * gg)

#drelgg = dspec * gg - (dspectropl * gg) # Also using actual humidity. How well does this do?

# Then make precipitation estimates:

dprlpbcgg = ((bb * prconm * dttrop) + (ff * prconm.mean() * drelbcgg)) / (1. - (bb * dttrop))

dprlpgg = ((bb * prconm * dttrop) + (ff * prconm.mean() * drelgg)) / (1. - (bb * dttrop))


# PLOT
# ----

# # Rainfall. RH known predictions. Quite good.

# plt.clf()
# plt.plot([-1.,2.],[-1.,2.],'k',linewidth=3)
# plt.plot([-1.,2.],[0.,0.],'k--',linewidth=3)
# plt.plot([0.,0.],[-1.,2.],'k--',linewidth=3)
# plt.plot(dprlp[:,0:15],dpr[:,0:15],'b.',alpha=0.75) # West
# plt.plot(dprlp[:,29:50],dpr[:,29:50],'r.',alpha=0.75) # East
# plt.xlabel('$\Delta P$ from equation 1')
# plt.ylabel('$\Delta P$')
# plt.title('$\Delta P$ reproduced from $\Delta r$')
# plt.ylim(-1.,2.)
# plt.xlim(-1.,2.)
# plt.savefig('isca_ideal_realsst_fig1.png')


# plt.clf()
# plt.plot([-1.,2.],[-1.,2.],'k',linewidth=3)
# plt.plot([-1.,2.],[0.,0.],'k--',linewidth=3)
# plt.plot([0.,0.],[-1.,2.],'k--',linewidth=3)
# plt.plot(dprlpbc[:,0:15],dpr[:,0:15],'b.',alpha=0.75) # West
# plt.plot(dprlpbc[:,29:50],dpr[:,29:50],'r.',alpha=0.75) # East
# plt.xlabel('$\Delta P$ from equation 1')
# plt.ylabel('$\Delta P$')
# plt.title('$\Delta P$ reproduced from $\Delta T$')
# plt.ylim(-1.,2.)
# plt.xlim(-1.,2.)
# plt.savefig('isca_ideal_realsst_fig1bc.eps')




# # With BC scaling. Surprisingly good. But because giving it T puts most of the answer in.

# plt.plot(dprlpbc[:,0:15],dpr[:,0:15],'b.',alpha=0.5) # West
# plt.plot(dprlpbc[:,29:50],dpr[:,29:50],'r.',alpha=0.5) # East
# plt.plot([-1.,1.5],[-1.,1.5],'k',linewidth=3)


# # With BC scaling using gg from control.

# plt.plot(dprlpbcgg[:,0:15],dpr[:,0:15],'bx',alpha=0.5) # West
# plt.plot(dprlpbcgg[:,29:50],dpr[:,29:50],'rx',alpha=0.5) # East
# plt.plot([-1.,1.5],[-1.,1.5],'k',linewidth=3)

# # checksum using gg from control but actual Spec changes.

# plt.plot(dprlpgg[:,0:15],dpr[:,0:15],'bx',alpha=0.5) # West
# plt.plot(dprlpgg[:,29:50],dpr[:,29:50],'rx',alpha=0.5) # East
# plt.plot([-1.,1.5],[-1.,1.5],'k',linewidth=3)






# # RH preds analysis.

# plt.plot(drelbc[:,0:15],drel[:,0:15],'b.',alpha=0.5) # West
# plt.plot(drelbc[:,29:50],drel[:,29:50],'r.',alpha=0.5) # East
# plt.plot([-0.1,0.1],[-0.1,0.1],'k',linewidth=3)

# # SH prediction. (Doesn't rely on knowing temperatures...)

# plt.plot(dspecbc[:,0:15],dspec[:,0:15],'b.',alpha=0.5) # West
# plt.plot(dspecbc[:,29:50],dspec[:,29:50],'r.',alpha=0.5) # East
# plt.plot([0.0005,0.003],[0.0005,0.003],'k',linewidth=3)

# # RH pred using gg from the control. First BC.

# plt.plot(drelbcgg[:,0:15],drel[:,0:15],'b.',alpha=0.5) # West
# plt.plot(drelbcgg[:,29:50],drel[:,29:50],'r.',alpha=0.5) # East
# plt.plot([-0.1,0.1],[-0.1,0.1],'k',linewidth=3)

# # Now RH pred using gg with actual spec change as a checksum.

# plt.plot(drelgg[:,0:15],drel[:,0:15],'b.',alpha=0.5) # West
# plt.plot(drelgg[:,29:50],drel[:,29:50],'r.',alpha=0.5) # East
# plt.plot([-0.1,0.1],[-0.1,0.1],'k',linewidth=3)


# # BC RH change predicted if there is no change in temperature., An initial cause if you like...

# plt.plot(dspecbc[:,0:15]/satconm[:,0:15],dspec[:,0:15]/satconm[:,0:15],'b.',alpha=0.5) # West
# plt.plot(dspecbc[:,29:50]/satconm[:,29:50],dspec[:,29:50]/satconm[:,29:50],'r.',alpha=0.5) # East
# plt.plot([0.0,0.15],[0.0,0.15],'k',linewidth=3)

# # RH against no T change RH change

# plt.plot(dspec[:,0:15]/satconm[:,0:15],drel[:,0:15],'b.',alpha=0.5) # West
# plt.plot(dspec[:,29:50]/satconm[:,29:50],drel[:,29:50],'r.',alpha=0.5) # East
# plt.plot([-0.1,0.1],[-0.1,0.1],'k',linewidth=3)

# # RH against BC no T change RH change

# plt.plot(dspecbc[:,0:15]/satconm[:,0:15],drel[:,0:15],'bx',alpha=0.5) # West
# plt.plot(dspecbc[:,29:50]/satconm[:,29:50],drel[:,29:50],'rx',alpha=0.5) # East
# plt.plot([-0.1,0.1],[-0.1,0.1],'k',linewidth=3)


# Figure for Marianne's paper.


plt.clf()
plt.figure(figsize=[6,6])

plt.subplot(2,2,1)
plt.plot([-0.1,0.15],[0.,0.],'k--',linewidth=3)
plt.plot([0.,0.],[-0.1,0.15],'k--',linewidth=3)
plt.plot([-0.1,0.15],[-0.1,0.15],'k',linewidth=3)
plt.plot(drelgg[:,0:15],drel[:,0:15],'b.',alpha=0.75) # West
plt.plot(drelgg[:,29:50],drel[:,29:50],'r.',alpha=0.75) # East
plt.title('(a) $\Delta r$ from $\Delta q$')
plt.xlabel('$\Delta r_{pred}$')
plt.ylabel('$\Delta r$')
plt.xlim(-0.1,0.15)
plt.ylim(-0.1,0.15)

plt.subplot(2,2,2)
plt.plot([0.,0.],[-1.,2.],'k--',linewidth=3)
plt.plot([-1.,2.],[0.,0.],'k--',linewidth=3)
plt.plot([-1.,2.],[-1.,2.],'k',linewidth=3)
plt.plot(dprlpgg[:,0:15],dpr[:,0:15],'b.',alpha=0.75) # West
plt.plot(dprlpgg[:,29:50],dpr[:,29:50],'r.',alpha=0.75) # East
plt.title('(b) $\Delta P$ from $\Delta q$')
plt.xlabel('$\Delta P_{pred}$')
plt.ylabel('$\Delta P$')
plt.xlim(-1.,2.)
plt.ylim(-1.,2.)

plt.subplot(2,2,3)
plt.plot([-0.1,0.15],[0.,0.],'k--',linewidth=3)
plt.plot([0.,0.],[-0.1,0.15],'k--',linewidth=3)
plt.plot([-0.1,0.15],[-0.1,0.15],'k',linewidth=3)
plt.plot(drelbcgg[:,0:15],drel[:,0:15],'b.',alpha=0.75) # West
plt.plot(drelbcgg[:,29:50],drel[:,29:50],'r.',alpha=0.75) # East
plt.title('(c) $\Delta r$ from $\Delta q_M$')
plt.xlabel('$\Delta r_{pred}$')
plt.ylabel('$\Delta r$')
plt.xlim(-0.1,0.15)
plt.ylim(-0.1,0.15)

plt.subplot(2,2,4)
plt.plot([0.,0.],[-1.,2.],'k--',linewidth=3)
plt.plot([-1.,2.],[0.,0.],'k--',linewidth=3)
plt.plot([-1.,2.],[-1.,2.],'k',linewidth=3)
plt.plot(dprlpbcgg[:,0:15],dpr[:,0:15],'b.',alpha=0.75) # West
plt.plot(dprlpbcgg[:,29:50],dpr[:,29:50],'r.',alpha=0.75) # East
plt.title('(d) $\Delta P$ from $\Delta q_M$')
plt.xlabel('$\Delta P_{pred}$')
plt.ylabel('$\Delta P$')
plt.xlim(-1.,2.)
plt.ylim(-1.,2.)

plt.subplots_adjust(top=0.92, bottom=0.1, left=0.15, right=0.95, hspace=0.4,
                    wspace=0.4)

plt.savefig('isca_ideal_realsst_fig2.png')



###### Marianne Plotting version #########
small = 14 #largefonts 14 # smallfonts 10 # medfonts = 14
med = 18 #largefonts 18 # smallfonts 14 # medfonts = 16
lge = 18
markcolor = ['c','m']
markstyle = '.'
marksize = 30.

fig, axes = plt.subplots(3, 2, figsize=(20,25))

axes[0,1].plot([-.8,1.8],[-.8,1.8],color='dimgray',linewidth=1)
axes[0,1].plot([-1.,2.],[0.,0.],color='dimgray',linewidth=1, linestyle = 'dashed')
axes[0,1].plot([0.,0.],[-1.,2.],color='dimgray',linewidth=1, linestyle = 'dashed')
axes[0,1].scatter(dprlp[:,0:15],dpr[:,0:15],c = markcolor[0], marker = markstyle, s = marksize, label = 'America') # West
axes[0,1].scatter(dprlp[:,29:50],dpr[:,29:50],c = markcolor[1], marker = markstyle, s = marksize, label = 'Africa') # East
axes[0,1].set_xlabel('$\Delta P_{pred}$', fontsize = lge)
axes[0,1].set_ylabel('$\Delta P_{actual}$', fontsize = lge)
axes[0,1].set_title('(a) $\Delta P$ from $\Delta r$', fontsize = lge)
axes[0,1].set_ylim(-1.,2.)
axes[0,1].set_xlim(-1.,2.)
axes[0,1].spines['right'].set_visible(False)
axes[0,1].spines['top'].set_visible(False)
axes[0,1].spines['left'].set_visible(False)
axes[0,1].spines['bottom'].set_visible(False)
axes[0,1].legend(fontsize = lge, markerscale = 2)

axes[1,0].plot([-0.1,0.15],[0.,0.],color='dimgray',linewidth=1, linestyle = 'dashed')
axes[1,0].plot([0.,0.],[-0.1,0.15],color='dimgray',linewidth=1, linestyle = 'dashed')
axes[1,0].plot([-0.08,0.13],[-0.08,0.13],color='dimgray',linewidth=1)
axes[1,0].scatter(drelgg[:,0:15],drel[:,0:15],c = markcolor[0], marker = markstyle, s = marksize) # West
axes[1,0].scatter(drelgg[:,29:50],drel[:,29:50],c = markcolor[1], marker = markstyle, s = marksize) # East
axes[1,0].set_title('(b) $\Delta r$ from $\Delta q_S$', fontsize = lge)
axes[1,0].set_xlabel('$\Delta r_{pred}$', fontsize = lge)
axes[1,0].set_ylabel('$\Delta r_{actual}$', fontsize = lge)
axes[1,0].set_xlim(-0.1,0.15)
axes[1,0].set_ylim(-0.1,0.15)
axes[1,0].spines['right'].set_visible(False)
axes[1,0].spines['top'].set_visible(False)
axes[1,0].spines['left'].set_visible(False)
axes[1,0].spines['bottom'].set_visible(False)

axes[1,1].plot([0.,0.],[-1.,2.],color='dimgray',linewidth=1, linestyle = 'dashed')
axes[1,1].plot([-1.,2.],[0.,0.],color='dimgray',linewidth=1, linestyle = 'dashed')
axes[1,1].plot([-.8,1.8],[-.8,1.8],color='dimgray',linewidth=1)
axes[1,1].scatter(dprlpgg[:,0:15],dpr[:,0:15],c = markcolor[0], marker = markstyle, s = marksize) # West
axes[1,1].scatter(dprlpgg[:,29:50],dpr[:,29:50],c = markcolor[1], marker = markstyle, s = marksize) # East
axes[1,1].set_title('(c) $\Delta P$ from $\Delta q_S$', fontsize = lge)
axes[1,1].set_xlabel('$\Delta P_{pred}$', fontsize = lge)
axes[1,1].set_ylabel('$\Delta P_{actual}$', fontsize = lge)
axes[1,1].set_xlim(-1.,2.)
axes[1,1].set_ylim(-1.,2.)
axes[1,1].spines['right'].set_visible(False)
axes[1,1].spines['top'].set_visible(False)
axes[1,1].spines['left'].set_visible(False)
axes[1,1].spines['bottom'].set_visible(False)

axes[2,0].plot([-0.1,0.15],[0.,0.],color='dimgray',linewidth=1, linestyle = 'dashed')
axes[2,0].plot([0.,0.],[-0.1,0.15],color='dimgray',linewidth=1, linestyle = 'dashed')
axes[2,0].plot([-0.08,0.13],[-0.08,0.13],color='dimgray',linewidth=1)
axes[2,0].scatter(drelbcgg[:,0:15],drel[:,0:15],c = markcolor[0], marker = markstyle, s = marksize) # West
axes[2,0].scatter(drelbcgg[:,29:50],drel[:,29:50],c = markcolor[1], marker = markstyle, s = marksize) # East
axes[2,0].set_title('(d) $\Delta r$ from $\Delta q_{TH}$', fontsize = lge)
axes[2,0].set_xlabel('$\Delta r_{pred}$', fontsize = lge)
axes[2,0].set_ylabel('$\Delta r_{actual}$', fontsize = lge)
axes[2,0].set_xlim(-0.1,0.15)
axes[2,0].set_ylim(-0.1,0.15)
axes[2,0].spines['right'].set_visible(False)
axes[2,0].spines['top'].set_visible(False)
axes[2,0].spines['left'].set_visible(False)
axes[2,0].spines['bottom'].set_visible(False)

axes[2,1].plot([0.,0.],[-1.,2.],color='dimgray',linewidth=1, linestyle = 'dashed')
axes[2,1].plot([-1.,2.],[0.,0.],color='dimgray',linewidth=1, linestyle = 'dashed')
axes[2,1].plot([-.8,1.8],[-.8,1.8],color='dimgray',linewidth=1)
axes[2,1].scatter(dprlpbcgg[:,0:15],dpr[:,0:15],c = markcolor[0], marker = markstyle, s = marksize) # West
axes[2,1].scatter(dprlpbcgg[:,29:50],dpr[:,29:50],c = markcolor[1], marker = markstyle, s = marksize) # East
axes[2,1].set_title('(e) $\Delta P$ from $\Delta q_{TH}$', fontsize = lge)
axes[2,1].set_xlabel('$\Delta P_{pred}$', fontsize = lge)
axes[2,1].set_ylabel('$\Delta P_{actual}$', fontsize = lge)
axes[2,1].set_xlim(-1.,2.)
axes[2,1].set_ylim(-1.,2.)
axes[2,1].spines['right'].set_visible(False)
axes[2,1].spines['top'].set_visible(False)
axes[2,1].spines['bottom'].set_visible(False)
axes[2,1].spines['left'].set_visible(False)

fig.delaxes(axes[0,0])

fig.savefig('isca_simple_scaling_plot.png', bbox_inches='tight', dpi=100)
fig.savefig('isca_simple_scaling_plot.pdf', bbox_inches='tight', dpi=400)
fig.savefig('isca_simple_scaling_plot.svg', bbox_inches='tight', dpi=100)