#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Preprocess ERA-20C Reanalysis


Variables to Investigate


ssr  (Surface net solar radiation)    - surface_net_downward_shortwave_flux ("J m**-2)
str  (Surface net thermal radiation)  - surface_net_upward_longwave_flux (J m**-2)

ssrd (Surface solar radiation downwards)

sshf () - "surface_upward_sensible_heat_flux"
slhf () - "surface_upward_latent_heat_flux"

Created on Fri May  6 16:12:23 2022

@author: gliu
"""

import xarray as xr
import numpy as np
import time
import matplotlib.pyplot as plt

import sys
sys.path.append("/home/glliu/00_Scripts/01_Projects/00_Commons/")
sys.path.append("/home/glliu/00_Scripts/01_Projects/01_AMV/02_stochmod/stochmod/model/")
from amv import proc
import scm


#%% ERA 20C

# Set up the file names and path
#vname_cmip   = ("rsus" ,"rlus" ,"rsds" ,"rlds" ,"hfss" ,"hfls","ts")
vname_era20c = ("ssr"  ,"str"  ,"sshf" ,"slhf" ,"skt")
vname_new    = ("fsns" ,"flns" ,"hfss" ,"hfls" ,"ts")

# Multiplier to make it positive downwards
flipdir      = [1,1,1,1,1]

checkvlims   = ([0,250],[-70,-40],[-70,25],[-200,0],[240,320])

# Get land and sea ice
icename = "ci/moda/ci.mon.mean.nc"

ncnames = ("ssr/moda/ssr.mon.mean.nc",
           "str/moda/str.mon.mean.nc",
           "sshf/moda/sshf.mon.mean.nc",
           "slhf/moda/slhf.mon.mean.nc",
           "skt/moda/skt.mon.mean.nc"
           )

# Mounted to volumes
datpath = "/Volumes/data/reanalysis/era.20c/sfc/"
ncs     = [datpath+nc for nc in ncnames]
outpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/"


lonname = 'longitude'
latname = 'latitude'
tname   = 'time'

# Rename the dims dict
olddims = [latname,lonname]
newdims = ['lat','lon']

dimswap = dict(zip(olddims,newdims))


dataset_name = "era20c"
#%% Set up land/ice mask

icethres  = 0.05
landthres = 0.05

# Make an ice mask (1948-01-01 to 2017-05-01)
ds_ice = xr.open_dataset(datpath+icename) # Load ice variable 
icefrac = ds_ice.ci.values

# Make the ice mask
ntime,nlat,nlon  = icefrac.shape
icemask          = np.zeros((nlat,nlon)) * np.nan # Preallocate NaN Array
icefree = np.array((icefrac >= icethres).sum(0))  # Get Ice Free Points
icemask[np.where(icefree==0)] = 1                 # Ice Free = 1

# Load land mask
landnc   = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/raw/era_landsea_mask.nc"
ds_lnd   = xr.open_dataset(landnc)
landmask = (ds_lnd.lsm.values).squeeze()
landmask[landmask>=landthres]  = np.nan
landmask[landmask<landthres] = 1

# Combine masks
limask   = landmask*icemask
plt.pcolormesh(ds_ice[lonname],ds_ice[latname],limask),plt.colorbar()

# Save landice mask
savename = "%s%s_limask_icefrac%.2f.npy" % (outpath,dataset_name,icethres)
np.save(savename,limask)

#%% Preprocess and check each variable

vid = 0

das = []

for vid in range(len(ncnames)):
    
    # Open Dataset
    ds = xr.open_dataset(ncs[vid])
    
    # Convert Units if needed
    if ds[vname_era20c[vid]].units == 'J m**-2':
        print("Converting from Jm**-2 to Wm**-2")
        ds[vname_era20c[vid]] /= 86400
    
    # Obsolete because of above, but leaving here for reference
    if vid == 0:
        # Compute seconds between each step
        times = ds.time.values
        dt = (times[1:] - times[:-1]).astype('timedelta64[s]').astype(dtype='int')
        dt = np.append(dt,dt[0]) # Append 31 days for last month (dec)
    
    # Apply the land-ice mask
    ds[vname_era20c[vid]] *= limask[None,:,:]
    
    # Rename and save
    ds         = ds.rename_dims(dimswap)
    ds         = ds[vname_era20c[vid]]
    ds         = ds.rename(dimswap)
    savenetcdf = "%s%s_%s.nc" % (outpath,dataset_name,vname_new[vid])
    ds         = ds.rename(vname_new[vid])
    ds.to_netcdf(savenetcdf,
             encoding={vname_new[vid]: {'zlib': True}})
    print("Saving netCDF to %s"% (savenetcdf,))
    
    # Check the plot
    vlims = checkvlims[vid]
    ds.mean('time').plot(cmap='jet',vmin=vlims[0],vmax=vlims[1])
    
    # Append the file
    das.append(ds)
    plt.show()

#%%

# Load data from above, compute qnet
vnames = ["fsns","flns","hfss","hfls"]
ds_qnet = scm.compute_qnet(outpath,dataset_name,vnames=vnames)

# Save output
savenetcdf = "%s%s_%s.nc" % (outpath,dataset_name,"qnet")
ds_qnet.to_netcdf(savenetcdf,
         encoding={"qnet": {'zlib': True}})
print("Saving netCDF to %s"% (savenetcdf,))

# Visualize to make sure
ds_qnet = xr.open_dataset(savenetcdf)
ds_qnet.qnet.mean('time').plot(cmap='jet',vmin=-200,vmax=200)

#%%







# # Set up the file names and path
# #vname_cmip   = ("rsus" ,"rlus" ,"rsds" ,"rlds" ,"hfss" ,"hfls","ts")
# vname_era20c = ("ssr"  ,"str"  ,"sshf" ,"slhf" ,"skt")
# vname_new    = ("fsns" ,"flns" ,"hfss" ,"hfls" ,"ts")

# # Get land and sea ice
# icename = "ci/moda/ci.mon.mean.nc"

# ncnames = ("ssr/moda/ssr.mon.mean.nc",
#            "str/moda/str.mon.mean.nc",
#            "sshf/moda/sshf.mon.mean.nc",
#            "slhf/moda/sshf.mon.mean.nc",
#            "skt/moda/skt.mon.mean.nc"
#            )

# # Mounted to volumes
# datpath = "/Volumes/data/reanalysis/era.20c/sfc/"
# ncs     = [datpath+nc for nc in ncnames]

# outpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/"

#%%


