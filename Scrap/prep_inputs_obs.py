#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Preparing Inputs for stochastic model for observations, SST

Currently written to run on Astraeus

Created on Mon Mar 17 11:56:39 2025

@author: gliu

"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import time
import sys
import cartopy.crs as ccrs
import glob
import matplotlib as mpl

import pandas as pd

import scipy as sp

#%% Import modules
stormtrack = 0
if stormtrack:
    sys.path.append("/home/glliu/00_Scripts/01_Projects/00_Commons/")
    sys.path.append("/home/glliu/00_Scripts/01_Projects/01_AMV/02_stochmod/stochmod/model/")
    
    # Path to the processed dataset (qnet and ts fields, full, time x lat x lon)
    #datpath =  "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_RCP85/01_PREPROC/"
    datpath =  "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/anom/"
    figpath =  "/home/glliu/02_Figures/01_WeeklyMeetings/20240621/"
    
else:
    sys.path.append("/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/03_Scripts/stochmod/model/")
    sys.path.append("/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/03_Scripts/")

    # Path to the processed dataset (qnet and ts fields, full, time x lat x lon)
    datpath =  "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/"
    figpath =  "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/02_Figures/20220511/"
from amv import proc,viz
import scm


import amv.proc as hf # Update hf with actual hfutils script, most relevant functions
import amv.loaders as dl
#%% User Edits

# Plot Settings
mpl.rcParams['font.family'] = 'Avenir'
proj = ccrs.PlateCarree()
bbplot = [-80, 0, 35, 75]
mons3 = proc.get_monstr()




# =========================
#%% (1) Load ERA5 HFF
# =========================

flxname         = "qnet"
dpath           = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/"
if flxname == "qnet":
    dof   = (2024-1979 + 1 - 2 - 1) * 3
    vname           = "qnet_damping"
    ncname_era5     = "ERA5_qnet_hfdamping_NAtl_1979to2024_ensorem1_detrend1.nc"
else:
    dof   = (2021-1979 + 1 - 2 - 1) * 3
    vname           = "thflx_damping"
    ncname_era5     = "ERA5_thflx_hfdamping_NAtl_1979to2021_ensorem1_detrend1.nc"
ds_era5         = xr.open_dataset(dpath + ncname_era5)

dt              = 3600*24  # *30 #Effective Processing Period of 1 day


lone = ds_era5.lon
late = ds_era5.lat
if flxname == "thflx": # Old version needed to be converted to months
    damping_era5 = ds_era5.thflx_damping / dt * -1
else:
    damping_era5 = ds_era5[vname] * -1

# Copy ERA5
hff_era5 = damping_era5.copy()


#%% Perform significance testing

# Set Significance calculation settings

hff   = damping_era5.copy()
rsst  = ds_era5.sst_autocorr.copy()
rflx  = ds_era5.sst_flx_crosscorr.copy()
setdict = {  # Taken from hfcalc_params
    'ensorem': 1,      # 1=enso removed, 0=not removed
    'ensolag': 1,      # Lag Applied toENSO and Variable before removal
    'monwin': 3,      # Size of month window for HFF calculations
    'detrend': 1,      # Whether or not variable was detrended
    'tails': 2,      # tails for t-test

    'p': 0.05,   # p-value for significance testing
    'sellags': [0,],   # Lags included (indices, so 0=lag1)
    'lagstr': "lag1",  # Name of lag based on sellags
    # Significance test option: 1 (No Mask); 2 (SST autocorr); 3 (SST-FLX crosscorr); 4 (Both), 5 (Replace with SLAB values)
    'method': 4
}

# dof was set above

# Compute and Apply Mask
st = time.time()
dampingmasked, freq_success, sigmask = scm.prep_HF(hff, rsst, rflx,
                                                   setdict['p'], setdict['tails'], dof, setdict['method'],
                                                   returnall=True)  # expects, [month x lag x lat x lon], should generalized with ensemble dimension?
print("Completed significance testing in %.2fs" % (time.time()-st))



dampingout = dampingmasked[:,setdict['sellags'],:,:].squeeze()
dampingout = xr.where(np.isnan(dampingout),0.,dampingout)
#dampingout[np.isnan(dampingout)] = 0


#%% Put into DataArray


# Get Dimensions
mons        = np.arange(1,13,1)
lon         = ds_era5.lon.values
lat         = ds_era5.lat.values

dims        = dict(mon=mons,lat=lat,lon=lon,)
da          = xr.DataArray(dampingout,name='damping',
                    dims=dims,coords=dims)

edict = proc.make_encoding_dict(da)
outpath_damping = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/03_reemergence/01_Data/proc/model_input/damping/"
outname         = outpath_damping + "ERA5_%s_damping_pilot.nc" % flxname
da.to_netcdf(outname,encoding=edict)


#%% Plot the HFF sigtest results

imon        = 9
ilag        = 0
fig, ax, _  = viz.init_orthomap(1, 1, bbplot, figsize=(14, 6))
ax          = viz.add_coast_grid(ax, bbplot, fill_color='lightgray',
                        proj=proj, line_color="dimgray")

pcm = ax.pcolormesh(lon,lat,dampingmasked[imon,ilag,:,:],
                   vmin=-35, vmax= 35,transform=proj,cmap='cmo.balance')

viz.hcbar(pcm)


# ================================
#%% (2) Load the MLD and get kprev
# ================================



dpath_proc = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/NATL_proc_obs/"
ncname     = dpath_proc + "MIMOC_RegridERA5_mld_NAtl_Climatology.nc"
ds_mld     = xr.open_dataset(ncname).load()

ds_mld     = ds_mld.rename(dict(mld='h'))
ds_mld     = ds_mld.rename(dict(month='mon'))

outpath_h  = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/03_reemergence/01_Data/proc/model_input/mld/"
outname    = outpath_h + "MIMOC_regridERA5_h_pilot.nc"
edict      = proc.make_encoding_dict(ds_mld)
ds_mld.to_netcdf(outname,encoding=edict)

#%% Compute Kprev

# Compute kprev for ens-mean mixed layer depth cycle
infunc = lambda x: scm.find_kprev(x,debug=False,returnh=False)
st     = time.time()

kprevall = xr.apply_ufunc(
    infunc, # Pass the function
    ds_mld, # The inputs in order that is expected
    input_core_dims =[['mon'],], # Which dimensions to operate over for each argument... 
    output_core_dims=[['mon'],], # Output Dimension
    vectorize=True, # True to loop over non-core dims
    )
print("Completed kprev calc in %.2fs" % (time.time()-st))


kprevall   = kprevall.transpose('mon','lat','lon').rename(dict(h='kprev'))
edict      = proc.make_encoding_dict(kprevall)
outname    = outpath_h + "MIMOC_regridERA5_kprev_pilot.nc"
kprevall.to_netcdf(outname,encoding=edict)

#%%

nlat = len(lat)
nlon = len(lon)

for a in range(nlat):
    for o in range(nlon):
        dspt = ds_mld.isel(lat=a,lon=o)
        out = scm.find_kprev(dspt.h.data,debug=False,returnh=False)
        
# ================================
#%% (3) Work on Land/Ice mask
# ================================ 

# Load Sea Ice Masks
dpath_ice = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/NATL_proc_obs/"
nc_masks = dpath_ice + "OISST_ice_masks_1981_2020.nc"
ds_masks = xr.open_dataset(nc_masks).load()

# Make Land Mask
land_nc     = dpath_ice + "ERA5_land_mask_1980_and_2020_NATL.nc"
ds_land     = xr.open_dataset(land_nc).load()

ds_land   = xr.where(np.isnan(ds_land.land_mask),1,np.nan,)


#ds_land  = xr.where(np.isnan(ds_mld.h.mean('mon')),np.nan,1)

# Adjust ice mask
ds_ice = ds_masks.mask_mon
ds_ice = xr.where(ds_ice,np.nan,1)

# Make the mask
#mask   = ds_land

#%% ok, I guess we have to regrid the oisst mask
# let's just do this conservatively
import tqdm

#ds_ice = ds_ice.data
ds_ice_arr = ds_ice.data

nlat = len(ds_ice.lat)
nlon = len(ds_ice.lon)

for a in tqdm.tqdm(range(nlat)):
    for o in range(nlon):
        
        dspt = ds_ice.isel(lon=o,lat=a)
        
        lonf = dspt.lon.data
        latf = dspt.lat.data
        
        ds_land_pt = proc.selpt_ds(ds_land,lonf,latf)
        
        if np.isnan(ds_land_pt.data):
            ds_ice[a,o] = 1
            
#%% Save the plot

# Reverse the mask
ds_ice_out = xr.where(ds_ice == 1,np.nan,1)
 # Ok I'm not sure whats going on but I will work twith this later
 
 
 
#%% (4) Do calculations for the forcing


vname  = "QNET"

ncpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/NATL_proc_obs/"
ncname = "ERA5_Fprime_%s_timeseries_%spilotObs_nroll0_NAtl.nc" % (vname,vname)
ds     = xr.open_dataset(ncpath+ncname).load()

# COmpute the monthly standard deviation
dsmon = ds.groupby('time.month').std('time')
dsmon = dsmon.rename(dict(month='mon'))
 
#outname = ""
#ds_ice_arr.plot()
outpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/03_reemergence/01_Data/proc/model_input/forcing/"
edict   = proc.make_encoding_dict(dsmon)
outname = outpath + "ERA5_Fprime_%s_std_pilot.nc" % vname
dsmon.to_netcdf(outname,encoding=edict)

#%% Repeat above but using qnet damping


        
        
        



#%%



#%%

