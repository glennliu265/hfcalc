#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 09:13:18 2022

Preprocesses CESM1-LE for Heat Flux Damping Calculations.
This is RCP85 version (should be written to readjust this)

- Takes TS, LandFrac, IceFrac, and makes SST
- Combines Heat Fluxes
- Prepares Ensemble Average Metric




Based on the following scripts:
    hf1_enavgrm

@author: gliu
"""



import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from tqdm import tqdm

#%% User Edits

datpath = "/vortex/jetstream/climate/data1/yokwon/CESM1_LE/downloaded/atm/proc/tseries/monthly/"
outpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_RCP85/01_PREPROC/"

mnum    = np.concatenate([np.arange(1,36),np.arange(101,106)])

# Part 1 (Land/Ice Mask Creation)
vnames  = ("LANDFRAC","ICEFRAC") # Variables
mthres  = (0.30,0.05) # Mask out if grid ever exceeds this value

# Part 2 ()
maskmode = "enssum"


#%% Functions


# RCP85 Loader
def load_rcp85(vname,N,datpath=None):
    if datpath is None:
        datpath = "/vortex/jetstream/climate/data1/yokwon/CESM1_LE/downloaded/atm/proc/tseries/monthly/"
        
    # Append variable name to path
    vdatpath = "%s%s/" % (datpath,vname)
        
    # Files are split into 2
    if N<34:
        fn1 = "b.e11.BRCP85C5CNBDRD.f09_g16.%03i.cam.h0.%s.200601-208012.nc" % (N,vname)
        fn2 = "b.e11.BRCP85C5CNBDRD.f09_g16.%03i.cam.h0.%s.208101-210012.nc" % (N,vname)
        ds = []
        for fn in [fn1,fn2]:
            dsf = xr.open_dataset(vdatpath + fn)
            ds.append(dsf)
        ds = xr.concat(ds,dim='time')
    else:
        fn1 = "%sb.e11.BRCP85C5CNBDRD.f09_g16.%03i.cam.h0.%s.200601-210012.nc" % (vdatpath,N,vname)
        ds = xr.open_dataset(fn1)
    return ds[vname]

# ----------------------------
#%% Part 1. Make Land/Ice Mask
# ----------------------------
""" Makes land/ice mask. Saves for all ens members separate and all summed."""
nens = len(mnum) 

# Initialize Mask
mask = [] #np.ones((nens,192,288))

# Loop for each ensemble member
for e in tqdm(range(nens)):
    N = mnum[e] # Get ensemble member
    emask = np.ones((192,288)) * np.nan # Preallocate
    for v in range(2):
    
        # Load dataset
        vname = vnames[v]
        ds = load_rcp85(vname,N,datpath=datpath)
        invar = ds.values # [Time x Lat x Lon]
        
        # Mask (set 0 where it ever exceeds, propagate in time)
        inthres       = mthres[v]
        if v == 0: # Landmask
            maskpts       = ((invar <= inthres).prod(0)) # [Lat x Lon]
            emask[maskpts==1] = 1 # 1 means it is ocean point
        elif v == 1:
            maskpts       = ((invar <= inthres).prod(0)) # [Lat x Lon]
            emask[maskpts==0] = np.nan # 0 means it has sea ice
    mask.append(emask.copy())
# Make into array
mask = np.array(mask)  # [ENS x LAT x LON]

# Save all members
savename = "%slandice_mask_rcp85_byens.npy" % (outpath)
np.save(savename,mask)

# Save ensemble sum
mask_enssum = mask.prod(0)
savename = "%slandice_mask_rcp85_ensavg.npy" % (outpath)
np.save(savename,mask_enssum)


# ------------------------------------------------------------
#%% For each variable: Apply LI Mask, Compute Ensemble Average
# ------------------------------------------------------------


usemask = np.load("%slandice_mask_rcp85_ensavg.npy" % (outpath)) # [Lat x Lon]
vnames  = ("TS","FSNS","FLNS","LHFLX","SHFLX")
nvar    = len(vnames)

for e in tqdm(range(nens)):
    
    # ********************************
    N = mnum[e]
    ensavg = np.zeros((2,1140,192,288)) # [Var x Time x Lat x Lon]
    qnet   = np.zeros((1140,192,288)) 
    
    for v,vname in enumerate(vnames):
        
        # Load the data
        ds    = load_rcp85(vname,N,datpath=datpath)
        
        # Apply the mask
        ds_msk = ds * usemask[None,:,:]
        
        
        # Just Save it
        if vname == "TS":
            
            # Add values for ensemble averaging
            ensavg[0,:,:,:] += ds_msk.values
            
            # Save the dataset
            savename = "%sCESM1_rcp85_%s_ens%02i.nc" % (outpath,"ts",e+1)
            ds_msk = ds_msk.rename('ts')
            ds_msk.to_netcdf(savename,encoding={'ts': {'zlib': True}})
            
        else:
            if vname == "FSNS":
                ds_msk *= -1 # Multiple to upwards positive
            qnet += ds_msk.values

    
    # Make/save qnet
    ensavg[1,:,:,:] += qnet.copy()
    coords  = {'time':ds_msk.time,'lat':ds_msk.lat,'lon':ds_msk.lon}
    da = xr.DataArray(qnet,
                dims=coords,
                coords=coords,
                name = 'qnet',
                )
    savename = "%sCESM1_rcp85_%s_ens%02i.nc" % (outpath,"qnet",e+1)
    da.to_netcdf(savename,
             encoding={'qnet': {'zlib': True}})
    
    
# Compute and save ensemble averages
vnames = ['ts','qnet']
    
for v in range(2):
    v_ensavg = ensavg[v,:,:,:]/nens
    
    da = xr.DataArray(v_ensavg,
                dims=coords,
                coords=coords,
                name = vnames[v],
                )
    savename = "%sCESM1_rcp85_%s_ensAVG.nc" % (outpath,vnames[v])
    da.to_netcdf(savename,
             encoding={vnames[v]: {'zlib': True}})

