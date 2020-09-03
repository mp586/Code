import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '/scratch/mp586/Code/MSE')
import gradients as gr, model_constants as mc
from pylab import rcParams
import sh
from scipy.odr import *
sys.path.insert(0, '/scratch/mp586/Code/PYCODES')
from plotting_routines_kav7 import *
GFDL_BASE = os.environ['GFDL_BASE']
sys.path.insert(0, os.path.join(GFDL_BASE,'src/extra/python/scripts'))
import cell_area as ca
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
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
    def column_int(var_in):
        var_int = var_in.sum('pfull')/mc.grav
        return var_int
run = '/scratch/mp586/Isca_DATA/ISCA_HPC/withtv/aquaplanet_frierson_insolation_0qflux_mld20_commitd15c267/run0121/atmos_monthly_interp.nc'
    data_dp = xr.open_dataset(run)
    dp = data_dp.phalf.diff('phalf')*100
    dp = xr.DataArray(dp[::-1], coords = [data.pfull.values], dims = ['pfull'])
run = '/scratch/mp586/Isca_DATA/ISCA_HPC/withtv/aquaplanet_frierson_inso
   ...: lation_0qflux_mld20_commitd15c267/run0121/atmos_monthly.nc'
run = '/scratch/mp586/Isca_DATA/ISCA_HPC/withtv/aquaplanet_frierson_insolation_0qflux_mld20_commitd15c267/run0121/atmos_monthly.nc'
    data_dp = xr.open_dataset(run)
    dp = data_dp.phalf.diff('phalf')*100
    dp = xr.DataArray(dp[::-1], coords = [data.pfull.values], dims = ['pfull'])
data = xr.open_mfdataset('/scratch/mp586/Isca_DATA/ISCA_HPC/withtv/aquaplanet_frierson_insolation_0qflux_mld20_commitd15c267/*/atmos_monthly_interp.nc')
data
ciwv = column_int(data.sphum * dp)
ciwv
data.sphum
sphum = data.sphum
civw = column_int(sphum*dp)
civw
dp
    dp = xr.DataArray(dp[::-1], coords = [data.pfull.values], dims = ['pfull'])
civw = column_int(sphum*dp)
civw
civw.values

    area_array, dx, dy = ca.cell_area_all(t_res=42,base_dir='/scratch/mp586/Isca/') # added _all because then dx and dy are also returned 
    area_array = xr.DataArray(area_array) # returned in units of m bzw m^2, because radius in cell_area.py is given in metres
    area_array, dx, dy = ca.cell_area_all(t_res=42,base_dir='/scratch/mp586/Isca/') # added _all because then dx and dy are also returned 
    area_array = xr.DataArray(area_array) # returned in units of m bzw m^2, because radius in cell_area.py is given in metres
    landfile=Dataset(os.path.join(GFDL_BASE,'input/aquaplanet/land.nc'),mode='r')
    landmask=landfile.variables['land_mask'][:]
    landmaskxrAP=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it
    landlats=landfile.variables['lat'][:]
    landlons=landfile.variables['lon'][:]
    landfile=Dataset(os.path.join(GFDL_BASE,'input/aquaplanet/land.nc'),mode='r')
    landmask=landfile.variables['land_mask'][:]
    landmaskxrAP=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it
        frac_llconv_avg, frac_llconv_sd = area_weighted_avg(civw, area_array, landmaskxrAP, 'all_sfcs')
ciwv
civw = column_int(sphum*dp)
civw
area_integral(civw[0,:,:],area_array,landmaskxr,'all_sfcs')
area_integral(civw[0,:,:],area_array,landmaskxrAP,'all_sfcs')
civw[0,:,:].sum()
np.sum(ciwv[0,:,:])
civw
ciwv =xr.DataArray(civw, coords=[civw.lat,ciwv.lon], dims = ['lat','lon'])
ciwv =xr.DataArray(civw, coords=[civw.time,civw.lat,ciwv.lon], dims = ['time','lat','lon'])
civw[0,:,:].sum()
ciwv[0,:,:].sum()
ciwv
civw.values
ciwv
area_array
area_array.dim_0
area_array = xr.DataArray(area_array, coords=[data.lat,data.lon], dims = ['lat','lon'])
area_integral(civw[0,:,:],area_array,landmaskxrAP,'all_sfcs')
(ciwv[0,:,:]*area_array).sum()
area_weighted_avg(ciwv[0,:,:],area_array,landmaskxrAP,'all_sfcs').sum()
area_weighted_avg(ciwv[0,:,:],area_array,landmaskxrAP,'all_sfcs')
area_weighted_avg(ciwv[0,:,:],area_array,landmaskxrAP,'all_sfcs')*area_array.sum()
ciwv_AA = ciwv*area_array
ciwv_AA
ciwv
area_Array
area_array
ciwv_AA
ciwv_AA.sum('lat','lon')
ciwv_AA.sum(['lat','lon'])
timeseries_totatmwv = ciwv_AA.sum(['lat','lon'])
timeseries_totatmwv.plot()
timeseries_totatmwv
plt.plot(timeseries_totatmwv.time,timeseries_totatmwv)
timeseries_totatmwv
time = np.linspace(0,360,1)
time
time = np.linspace(0,360,361)
time
timeseries_totatmwv.time
time
len(time)
time = np.linspace(1,360,361)
time
time = np.linspace(0,360,360)
time
time = np.linspace(1,360,360)
time
plt.plot(time,timeseries_totatmwv)
ciwv
ciwv
civw
ciwv.groupby('time.year').mean('time')
ciwv.groupby('time.season').mean('time')
runmin = 121
runmin = 120
runmin = 121
runmax = 481
time=[np.array(np.linspace(0,(runmax-runmin-1),(runmax-runmin),dtype='datetime64[M]'))]
ciwv
civw
ciwv = xr.DataArray(ciwv,coords = [time, data.lat, data.lon], dims = ['time','lat','lon'])
ciwv
ciwv = xr.DataArray(ciwv.values,coords = [time, data.lat, data.lon], dims = ['time','lat','lon'])
ciwv = xr.DataArray(np.asarray(values)
,coords = [time, data.lat, data.lon], dims = ['time','lat','lon'])
ciwv = xr.DataArray(np.asarray(ciwv),coords = [time, data.lat, data.lon], dims = ['time','lat','lon'])
ciwv
ciwv.coords['time']
ciwv.groupby('time.year')
ciwv.assign_coords('time',time)
ciwv.assign_coords(time=time)
time
ciwv.time
time = xr.DataArray(time, coords = [time], dims = ['time'])
time
time = xr.DataArray(np.asarray(time), coords = [time], dims = ['time'])
time = np.linspace(1,360,360)
plt.plot(time,timeseries_totatmwv)
plt.plot(time,timeseries_totatmwv.rolling(time=12, center=True).mean())
civw
civw.groupby('time.year').mean()
data.sphum.groupby('time.year')
data.time
civw = xr.DataArray(civw.values,coords = [time, data.lat, data.lon], dims = ['time','lat','lon'])
    time=[np.array(np.linspace(0,(runmax-runmin-1),(runmax-runmin),dtype='datetime64[M]'))] # has to start at zero so that it gives jan to dec. otherwise (1,12) goes from feb to jan!
civw = xr.DataArray(civw.values,coords = [time, data.lat, data.lon], dims = ['time','lat','lon'])
data.sphum.groupby('time.year').mean('time')
civw
civw.time
civw.groupby('year')
civw.groupby('time.year')
civw.groupby('year').mean(dim='time')
[sphum_ctl,sphum_avg_ctl,sphum_seasonal_avg_ctl,sphum_month_avg_ctl,time]=seasonal_4D_variable('Isca_DATA/ISCA_HPC/withtv/aquaplanet_frierson_insolation_0qflux_mld20_commitd15c267','isca',1,481,'sphum','kg/kg')
[sphum_ctl,sphum_avg_ctl,sphum_seasonal_avg_ctl,sphum_month_avg_ctl,sphum_annual_avg_ctl,time]=seasonal_4D_variable('Isca_DATA/ISCA_HPC/withtv/aquaplanet_frierson_insolation_0qflux_mld20_commitd15c267','isca',1,481,'sphum','kg/kg')
%hist
sphum_annual_avg_ctl
    dp = xr.DataArray(dp[::-1], coords = [data.pfull.values], dims = ['pres_lev'])
column_int(sphum_annual_avg_ctl*dp)
    def column_int(var_in):
        var_int = var_in.sum('pres_lev')/mc.grav
        return var_int
column_int(sphum_annual_avg_ctl*dp)
sphum_annual_avg_ctl
sphum_annual_avg_ctl*dp
dp
sphum_annual_avg_ctl
sphum_annual_avg_ctl*dp
dp
dp.pres_lev
sphum_annual_avg_ctl.pres_lev
    data_dp = xr.open_dataset('/scratch/mp586/Isca_DATA/ISCA_HPC/withtv/'+run+'/run0121/atmos_monthly.nc')
    dp = data_dp.phalf.diff('phalf')*100
    dp = xr.DataArray(dp, coords = [data.pfull.values], dims = ['pres_lev'])
    data_dp = xr.open_dataset(run)
    dp = data_dp.phalf.diff('phalf')*100
    dp = xr.DataArray(dp, coords = [data.pfull.values], dims = ['pres_lev'])
dp
sphum_ctl.pres_lev
ciwv
data.sphum
data.pfull
dp
sphum_annual_avg_ctl*dp
sphum_annual_avg_ctl
dp
time
ciwv = xr.DataArray(ciwv, coords = [time[0],data.lat,data.lon], dims = ['time','lat','lon'])
    time=[np.array(np.linspace(0,480,361,dtype='datetime64[M]'))]
time
len(time)
len(time[0])
    time=[np.array(np.linspace(1,360,361,dtype='datetime64[M]'))]
len(time[0])
%hist
runmin = 121
runmax = 481
time=[np.array(np.linspace(0,(runmax-runmin-1),(runmax-runmin),dtype='datetime64[M]'))]
ciwv = xr.DataArray(ciwv,coords = [time[0], data.lat, data.lon], dims = ['time','lat','lon'])
ciwv.groupby('time.year')
ciwv.groupby('time.year').mean('time')
time_int = np.linspace(1,360,360)
plt.plot(time_int,timeseries_totatmwv)
plt.plot(time_int,timeseries_totatmwv)
plt.plot(ciwv.groupby('time.year').mean('time'))
ciwv_yrs =  ciwv.groupby('time.year').mean('time')
ciwv_yrs_AA = ciwv_yrs*area_array
timeseries_ciwv_yrs = ciwv_yrs_AA.sum()
plt.plot(np.linspace(1,30,30),timeseries_ciwv_yrs)
timeseries_ciwv_yrs
timeseries_ciwv_yrs = ciwv_yrs_AA.sum(['lat','lon'])
plt.plot(np.linspace(1,30,30),timeseries_ciwv_yrs)
plt.plot(time_int,timeseries_totatmwv)
timeseries_ciwv_yrs = ciwv_yrs_AA.sum(['lat','lon'])
plt.plot(np.linspace(1,30,30),timeseries_ciwv_yrs)
%hist
len(data.time)
runmin-runmax
data
def column_int(var):
    var_int = var_in.sum('pfull')/mc.grav
    return var_int

def check_total_atmospheric_water(data, dp = None, continents = 'aquaplanet'):

    if dp == None:
        data_dp = xr.open_dataset(run)
        dp = data_dp.phalf.diff('phalf')*100
        dp = xr.DataArray(dp[::-1], coords = [data.pfull.values], dims = ['pfull'])

    area_array, dx, dy = ca.cell_area_all(t_res=42,base_dir='/scratch/mp586/Isca/') # added _all because then dx and dy are also returned 
    area_array = xr.DataArray(area_array) # returned in units of m bzw m^2, because radius in cell_area.py is given in metres
    landfile=Dataset(os.path.join(GFDL_BASE,'input/'+continents+'/land.nc'),mode='r')
    landmask=landfile.variables['land_mask'][:]
    landlats=landfile.variables['lat'][:]
    landlons=landfile.variables['lon'][:]
    landmaskxrAP=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

    civw = column_int(sphum*dp)

    time=[np.array(np.linspace(0,len(data.time),len(data.time),dtype='datetime64[M]'))]
    ciwv = xr.DataArray(ciwv,coords = [time[0], data.lat, data.lon], dims = ['time','lat','lon'])
    ciwv_yrs =  ciwv.groupby('time.year').mean('time')
    ciwv_yrs_AA = ciwv_yrs*area_array
    timeseries_ciwv_yrs = ciwv_yrs_AA.sum()
    plt.plot(np.linspace(1,len(ciwv.time),len(ciwv.time)),timeseries_ciwv_yrs)
check_total_atmospheric_water(data)
def column_int(var_in):
    var_int = var_in.sum('pfull')/mc.grav
    return var_int

def check_total_atmospheric_water(data, dp = None, continents = 'aquaplanet'):

    if dp == None:
        data_dp = xr.open_dataset(run)
        dp = data_dp.phalf.diff('phalf')*100
        dp = xr.DataArray(dp[::-1], coords = [data.pfull.values], dims = ['pfull'])

    area_array, dx, dy = ca.cell_area_all(t_res=42,base_dir='/scratch/mp586/Isca/') # added _all because then dx and dy are also returned 
    area_array = xr.DataArray(area_array) # returned in units of m bzw m^2, because radius in cell_area.py is given in metres
    landfile=Dataset(os.path.join(GFDL_BASE,'input/'+continents+'/land.nc'),mode='r')
    landmask=landfile.variables['land_mask'][:]
    landlats=landfile.variables['lat'][:]
    landlons=landfile.variables['lon'][:]
    landmaskxrAP=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

    civw = column_int(sphum*dp)

    time=[np.array(np.linspace(0,len(data.time),len(data.time),dtype='datetime64[M]'))]
    ciwv = xr.DataArray(ciwv,coords = [time[0], data.lat, data.lon], dims = ['time','lat','lon'])
    ciwv_yrs =  ciwv.groupby('time.year').mean('time')
    ciwv_yrs_AA = ciwv_yrs*area_array
    timeseries_ciwv_yrs = ciwv_yrs_AA.sum()
    plt.plot(np.linspace(1,len(ciwv.time),len(ciwv.time)),timeseries_ciwv_yrs)
def column_int(var_in):
    var_int = var_in.sum('pfull')/mc.grav
    return var_int

def check_total_atmospheric_water(data, dp = None, continents = 'aquaplanet'):

    if dp == None:
        data_dp = xr.open_dataset(run)
        dp = data_dp.phalf.diff('phalf')*100
        dp = xr.DataArray(dp[::-1], coords = [data.pfull.values], dims = ['pfull'])

    area_array, dx, dy = ca.cell_area_all(t_res=42,base_dir='/scratch/mp586/Isca/') # added _all because then dx and dy are also returned 
    area_array = xr.DataArray(area_array) # returned in units of m bzw m^2, because radius in cell_area.py is given in metres
    landfile=Dataset(os.path.join(GFDL_BASE,'input/'+continents+'/land.nc'),mode='r')
    landmask=landfile.variables['land_mask'][:]
    landlats=landfile.variables['lat'][:]
    landlons=landfile.variables['lon'][:]
    landmaskxrAP=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

    civw = column_int(sphum*dp)

    time=[np.array(np.linspace(0,len(data.time),len(data.time),dtype='datetime64[M]'))]
    ciwv = xr.DataArray(ciwv,coords = [time[0], data.lat, data.lon], dims = ['time','lat','lon'])
    ciwv_yrs =  ciwv.groupby('time.year').mean('time')
    ciwv_yrs_AA = ciwv_yrs*area_array
    timeseries_ciwv_yrs = ciwv_yrs_AA.sum()
    plt.plot(np.linspace(1,len(ciwv.time),len(ciwv.time)),timeseries_ciwv_yrs)
    plt.show()
check_total_atmospheric_water(data)
def check_total_atmospheric_water(data, dp = None, continents = 'aquaplanet'):

    if dp == None:
        data_dp = xr.open_dataset(run)
        dp = data_dp.phalf.diff('phalf')*100
        dp = xr.DataArray(dp[::-1], coords = [data.pfull.values], dims = ['pfull'])

    area_array, dx, dy = ca.cell_area_all(t_res=42,base_dir='/scratch/mp586/Isca/') # added _all because then dx and dy are also returned 
    area_array = xr.DataArray(area_array) # returned in units of m bzw m^2, because radius in cell_area.py is given in metres
    landfile=Dataset(os.path.join(GFDL_BASE,'input/'+continents+'/land.nc'),mode='r')
    landmask=landfile.variables['land_mask'][:]
    landlats=landfile.variables['lat'][:]
    landlons=landfile.variables['lon'][:]
    landmaskxrAP=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

    ciwv = column_int(sphum*dp)

    time=[np.array(np.linspace(0,len(data.time),len(data.time),dtype='datetime64[M]'))]
    ciwv = xr.DataArray(ciwv,coords = [time[0], data.lat, data.lon], dims = ['time','lat','lon'])
    ciwv_yrs =  ciwv.groupby('time.year').mean('time')
    ciwv_yrs_AA = ciwv_yrs*area_array
    timeseries_ciwv_yrs = ciwv_yrs_AA.sum()
    plt.plot(np.linspace(1,len(ciwv.time),len(ciwv.time)),timeseries_ciwv_yrs)
    plt.show()
check_total_atmospheric_water(data)
def column_int(var_in):
    var_int = var_in.sum('pfull')/mc.grav
    return var_int

def check_total_atmospheric_water(data, dp = None, continents = 'aquaplanet'):

    if dp == None:
        data_dp = xr.open_dataset(run)
        dp = data_dp.phalf.diff('phalf')*100
        dp = xr.DataArray(dp[::-1], coords = [data.pfull.values], dims = ['pfull'])

    area_array, dx, dy = ca.cell_area_all(t_res=42,base_dir='/scratch/mp586/Isca/') # added _all because then dx and dy are also returned 
    area_array = xr.DataArray(area_array) # returned in units of m bzw m^2, because radius in cell_area.py is given in metres
    landfile=Dataset(os.path.join(GFDL_BASE,'input/'+continents+'/land.nc'),mode='r')
    landmask=landfile.variables['land_mask'][:]
    landlats=landfile.variables['lat'][:]
    landlons=landfile.variables['lon'][:]
    landmaskxrAP=xr.DataArray(landmask,coords=[landlats,landlons],dims=['lat','lon']) # need this in order to use .sel(... slice) on it

    ciwv = column_int(data.sphum*dp)

    time=[np.array(np.linspace(0,len(data.time),len(data.time),dtype='datetime64[M]'))]
    ciwv = xr.DataArray(ciwv,coords = [time[0], data.lat, data.lon], dims = ['time','lat','lon'])
    ciwv_yrs =  ciwv.groupby('time.year').mean('time')
    ciwv_yrs_AA = ciwv_yrs*area_array
    timeseries_ciwv_yrs = ciwv_yrs_AA.sum()
    plt.plot(np.linspace(1,len(ciwv.time),len(ciwv.time)),timeseries_ciwv_yrs)
    plt.show()
check_total_atmospheric_water(data)
%hist