#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

------------------
Make Land Ice Mask
------------------

Given ICEFRAC and LANDFRAC, make a land-ice mask by excluding all points where 
a certain percentage is exceeded (see [mthres] variable)

Inputs

Monthly LANDFRAC and ICEFRAC (in percent): [<ens> x time x lat x lon360]

Created on Fri Jun 21 13:41:15 2024

@author: gliu
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from tqdm import tqdm
import xarray as xr
import glob

import time
import sys
import scipy as sp

#%% Load Functions for now from amv.proc. but eventually make a helper func [hfutils --> hf]
# Copied stormtrack import from DATA_CODE_CENTRAL.md

# stormtrack
amvpath = "/home/glliu/00_Scripts/01_Projects/00_Commons/" # amv module

sys.path.append(amvpath)

from amv import proc as hf

#%% Functions

#%% User Edits

# Model Output/Dataset Information (searches for files in <datpath>/)
datname         = "cesm1_htr_5degbilinear" #"cesm2_pic" # Dataset
datpath         = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/proc/" # Where to search
vnames          = ("LANDFRAC","ICEFRAC") # Land and Ice fraction variable names
croptime        = True
tstart          = '1920-01-01'
tend            = '2005-12-31'

outpath       = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/masks/"

## CESM2 PiControl
# datname         = "cesm2_pic"
# datpath         = "/stormtrack/data4/glliu/01_Data/CESM2_PiControl/FCM/atm/"
# searchstr       = 
# tstart          = '0200-01-01'
# tend            = '2000-12-31'

# Time Crop
timestr         = 'year%sto%s' % (tstart[:4],tend[:4])

# Land Ice Mask Creation Options
mthres          = (0.30,0.05) # (Land Thres, Ice Thres) Mask out if grid ever exceeds this percentage value 
maskname        = "land%03i_ice%03i" % (mthres[0]*100,mthres[1]*100)


#%% Declare Functions

def make_limask(landfrac,icefrac,landthres,icethres):
    
    # Get Maximum Percentage over time
    maxicethres  = icefrac.ICEFRAC.max('time') # [Lat x Lon]
    maxlandthres = landfrac.LANDFRAC.max('time') # [Lat x Lon]
    
    # Apply Thresholds
    icemask  = xr.where(maxicethres>icethres,np.nan,1)
    landmask = xr.where(maxlandthres>landthres,np.nan,1) 
    
    return landmask.rename('mask'),icemask.rename('mask')

# ----------------------------
#%% Part 1. Make Land/Ice Mask
# ----------------------------
""" Makes land/ice mask. Saves for all ens members separate and all summed."""

# Loop LANDFRAC, ICEFRAC
if datname == "cesm2_pic":
    maskvar = []
    for vv in range(2): 
        
        # Create filename/list and load
        vname     = vnames[vv]
        keepvars  = ["time","lat","lon",vname]
        searchstr = "%s%s/*%s*.nc" % (datpath,vname,vname) # Searches for datpath + *LANDFRAC*.nc
        nclist    = glob.glob(searchstr)
        nclist.sort()
        
        # Drop Unnecessary variables
        ds_all    = xr.open_mfdataset(nclist,concat_dim="time",combine='nested')
        ds_all    = hf.ds_dropvars(ds_all,keepvars)
        ds_all    = hf.fix_febstart(ds_all)
        
        # Load the Data
        st = time.time()
        ds_all = ds_all.load()
        print("Loaded in %.2fs" % (time.time()-st))
        
        # Append to files
        maskvar.append(ds_all)
    
else: # Just do simple search for other variables (currently works with cesm1_htr_5degbilinear)
    maskvar = []
    for vv in range(2):
        
        # Create filename/list and load
        vname     = vnames[vv]
        keepvars  = ["time","lat","lon",vname]
        searchstr = "%s%s*%s*.nc" % (datpath,datname,vname) # Searches for datpath + *LANDFRAC*.nc
        nclist    = glob.glob(searchstr)
        nclist.sort()
        
        # Load and drop unnecessary variables
        ds_all    = xr.open_mfdataset(nclist,concat_dim="ens",combine='nested')
        ds_all    = hf.ds_dropvars(ds_all,keepvars)
        ds_all    = hf.fix_febstart(ds_all)
        
        # Load the Data
        st = time.time()
        ds_all = ds_all.load() # ('ens', 'time', 'lat', 'lon')
        print("Loaded in %.2fs" % (time.time()-st))
        
        # Append to files
        maskvar.append(ds_all)

#%% Make the Mask (Should be universal for all models)



st                  = time.time()
# Crop time if optionis set
if croptime:
    print("Cropping time between %s and %s" % (tstart,tend))
    maskvar = [ds.sel(time=slice(tstart,tend)) for ds in maskvar]
landfrac,icefrac    = maskvar
landthres,icethres  = mthres

# Make the mask
landmask,icemask    = make_limask(landfrac,icefrac,landthres,icethres)
print("Computed mask in %.2fs" % (time.time()-st))


limask              = landmask*icemask
outputs             = [landmask,icemask,limask]
outnames            = ["landmask_%02ip" % (landthres*100),
                       'icemask_%02ip' % (icethres*100),
                       'limask_%sp_%sp' % (landthres,icethres)]
edict               = dict(mask=dict(zlib=True))

for vv in range(3):
    outname = "%s%s_%s.nc" % (outpath,datname,outnames[vv])
    if croptime:
        outname = hf.addstrtoext(outname,"_%s" % timestr,adjust=-1)
    outputs[vv].to_netcdf(outname,encoding=edict)
    
    if 'ens' in list(outputs[vv].dims):
        print("Saving Ensemble Sum for %s" % (outnames[vv]))
        output_enssum = outputs[vv].prod('ens',skipna=False)
        outname_new   = hf.addstrtoext(outname,"_enssum",adjust=-1)
        output_enssum.to_netcdf(outname_new,encoding=edict)
        


