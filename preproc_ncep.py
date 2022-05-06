#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

-------------------------
Preproc NCEP Reanalysis 1
-------------------------

Preprocess Monthly Means
1948 - 2017.05.01. 2022

rsus (shortwave), upward
rlus (longwave), upward

rsds (shortwave), downward
rlds (longwave), downward


hfss (sensible hflx)
hfls (latent hflx)

To get started, mount the WHOI CMIP5 server by using
"connect to folder" : smb\\cmip5.whoi.edu\data


General steps
- Combine upwelling/downwelling radiation terms
- Calculate Net Heat Flux, downward positive
- Apply Land/Ice Mask to everything
- Save output with variable names

Created on Thu May  5 11:31:53 2022

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
#%%

# Set up the file names and path
vname_cmip = ("rsus" ,"rlus" ,"rsds" ,"rlds" ,"hfss" ,"hfls","ts")
vname_ncep = ("uswrf","ulwrf","dswrf","dlwrf","shtfl","lhtfl","skt")
vname_new  = ("fsns" ,"flns" ,"fsns" ,"flns" ,"hfss","hfls","ts")

ncnames = ('uswrf.sfc/monthly/uswrf.sfc.mon.mean.nc',
           'ulwrf.sfc/monthly/ulwrf.sfc.mon.mean.nc',
           'dswrf.sfc/monthly/dswrf.sfc.mon.mean.nc',
           'dlwrf.sfc/monthly/dlwrf.sfc.mon.mean.nc',
           'shtfl.sfc/monthly/shtfl.sfc.mon.mean.nc',
           'lhtfl.sfc/monthly/lhtfl.sfc.mon.mean.nc',
           'skt.sfc/monthly/skt.sfc.mon.mean.nc'
           )

# Get land and sea ice
icename = "icec.sfc/monthly/icec.sfc.mon.mean.nc"

# Mounted to volumes
datpath = "/Volumes/data/reanalysis/ncep.reanalysis.i/surface_gauss/"
ncs     = [datpath+nc for nc in ncnames]


outpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/"

# Part 1: Applying the Land/Ice Mask ******************************************
#%% Make land-ice mask

icethres = 0.05

# Make an ice mask (1948-01-01 to 2017-05-01)
ds_ice = xr.open_dataset(datpath+icename) # Load ice variable 
icefrac = ds_ice.icec.values              

# Make the ice mask
ntime,nlat,nlon  = icefrac.shape
icemask          = np.zeros((nlat,nlon)) * np.nan # Preallocate NaN Array
icefree = np.array((icefrac >= icethres).sum(0))  # Get Ice Free Points
icemask[np.where(icefree==0)] = 1                 # Ice Free = 1

# Load land mask
landnc   = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/raw/land.sfc.gauss.nc"
ds_lnd   = xr.open_dataset(landnc)
landmask = (ds_lnd.land.values).squeeze()
landmask[landmask==1] = np.nan
landmask[landmask==0] = 1

# Combine masks
limask = landmask*icemask
plt.pcolormesh(ds_ice.lon,ds_ice.lat,limask),plt.colorbar()

# Save landice mask
savename = "%sncep_ncar_limask_icefrac%.2f.npy" % (outpath,icethres)
#savename = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/ncep_ncar_limask360_fract.npy"
np.save(savename,limask)
#%% Lets examine/process surface shortwave...

# Open shortwave datasets
dsu = xr.open_dataset(ncs[0]) # upward
dsd = xr.open_dataset(ncs[2]) # downward

# Calculate NET fsns (downwards positive)
ds_fsns = (dsd.dswrf - dsu.uswrf) * limask[None,:,:]

# Rename and save
savenetcdf = "%sncep_ncar_fsns.nc" % outpath
ds_fsns    = ds_fsns.rename('fsns')
print(ds_fsns)
ds_fsns.to_netcdf(savenetcdf,
         encoding={'fsns': {'zlib': True}})
print("Saving netCDF to %s"% (savenetcdf,))

# Load and check
ds_fsns = xr.open_dataset(savenetcdf)
ds_fsns.fsns.mean('time').plot(cmap='jet',vmin=0,vmax=250)

#%% Do the same for longwave

# Open longwave datasets
dlu = xr.open_dataset(ncs[1]) # upward
dld = xr.open_dataset(ncs[3]) # downward

# Calculate NET fsns (downwards positive)
ds_flns = (dld.dlwrf - dlu.ulwrf) * limask[None,:,:]

# Rename and save
savenetcdf = "%sncep_ncar_flns.nc" % outpath
ds_flns    = ds_flns.rename('flns')
print(ds_flns)
ds_flns.to_netcdf(savenetcdf,
         encoding={'flns': {'zlib': True}})
print("Saving netCDF to %s"% (savenetcdf,))

# Load and check
ds_flns = xr.open_dataset(savenetcdf)
ds_flns.flns.mean('time').plot(cmap='jet',vmin=-70,vmax=-40)

#%% Apply land-ice mask to the rest of the variables

for i in [4,5,6]:
    
    # Apply limask
    ds = xr.open_dataset(ncs[i])
    
    ds *= limask[None,:,:] 
    if i < 6:
        ds *= -1 # flip to downwards positive
    
    # Save netcdf
    savenetcdf = "%sncep_ncar_%s.nc" % (outpath,vname_new[i])
    dsn = ds[vname_ncep[i]].rename(vname_new[i])
    print(ds)
    dsn.to_netcdf(savenetcdf,
             encoding={vname_new[i]: {'zlib': True}})
    print("Saving netCDF to %s"% (savenetcdf,))
    
    ds.close()
    dsn.close()
    
    # Load and check
    ds = xr.open_dataset(savenetcdf)
    if i == 4:
        ds[vname_new[i]].mean('time').plot(cmap='jet',vmin=-70,vmax=25)
    elif i == 5:
        ds[vname_new[i]].mean('time').plot(cmap='jet',vmin=-200,vmax=0)
    elif i == 6:
        ds[vname_new[i]].mean('time').plot(cmap='jet')
    plt.show()
    ds.close()
    

#%% Load and compute Qnet

def compute_qnet(datpath,dataset_name,vnames=None,downwards_positive=True,
                 verbose=True):
    """
    Computes dataarray of Qnet (downwards positive) given shortwave, longwave,
    sensible, and latent heat fluxes. Assumes inputs are positive into the ocean.
    
    Looks for files of the form <datpath><dataset name>_<flux name>.nc

    Parameters
    ----------
    datpath : STR
        Path to data
    dataset_name : STR
        Name of dataset
    vnames : List of STR, optional
        Name of shortwave,longwave,sensible, and latent fluxes. The default is:
        ["fsns","flns","hfss","hfls"]
    downwards_positive : BOOL, optional
        Set to False to convert to upwards positive. The default is True.
    verbose : BOOL, optional
        Print debuggin messages. The default is True.
    
    Returns
    -------
    ds_new : TYPE
        DESCRIPTION.

    """
    
    if vnames is None:
        vnames = ["fsns","flns","hfss","hfls"]
    
    flxes = []
    for v in vnames:
        loadname = "%s%s_%s.nc" % (datpath,dataset_name,v)
        if verbose:
            print("Loading %s!" % loadname)
        ds = xr.open_dataset(loadname)
        flxes.append(ds[v])
        
    # Compute Qnet = FSNS + (FLNS + SHFLX + LHFLX)
    ds_new = flxes[0] + (flxes[1]+flxes[2]+flxes[3])
    ds_new = ds_new.rename("qnet")
    
    if downwards_positive == False:
        ds_new *= -1 # Convert to downwards positive
        
    return ds_new
    

# Load data from above, compute qnet
vnames = ["fsns","flns","hfss","hfls"]
ds_qnet = compute_qnet(outpath,'ncep_ncar',vnames=vnames)

# Save output
savenetcdf = "%sncep_ncar_%s.nc" % (outpath,"qnet")
ds_qnet.to_netcdf(savenetcdf,
         encoding={"qnet": {'zlib': True}})
print("Saving netCDF to %s"% (savenetcdf,))


# Visualize to make sure
ds_qnet = xr.open_dataset(savenetcdf)
ds_qnet.qnet.mean('time').plot(cmap='jet',vmin=-200,vmax=200)

#%%





    



