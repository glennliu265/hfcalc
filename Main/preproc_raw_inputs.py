#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Preprocess Raw Inputs For HFF Calculations/Stochastic Model Parameterizations


Given a loader, datasetname, and variables

- Load + concat files
- Standardize the naming + flip longitude (time x lat x lon)
- Fix February Start
- Crop to Time
- Crop to Region ==Make region dict==#ToDo
- Save to `proc`

Input:
    Loader, Dataset Name, Variable Names
Output:
    Variable: [time x lat x lon180]
    Name:
    <datname>_<vname>_<bbox_name>_<timestr>.nc in /proc/ directory
    ex. cesm2_pic_LHFLX_NAtl_0200to2000.nc
    
    
    

Created on Fri Jun 21 13:48:22 2024

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

#%% Import Modules (report with hfutils)

# stormtrack
amvpath = "/home/glliu/00_Scripts/01_Projects/00_Commons/" # amv module
sys.path.append(amvpath)
import amv.proc as hf

#%% Write a loader


def load_cesm2_pic(vname,datpath=None,searchstr=None,debug=True):
    # Function to Load CESM2 PiC Output on stormtrack, Searches for 
    # File is assumed to be at <datpath>/<vname>/*<vname>*.nc
    
    # For Raw Atmospheric Variables
    if vname == "HMXL":
        print("Loading HMXL that has be regridded by [prep_mld_PIC]!")
        # Load variable that has been processed by prep_mld_PIC.py
        if datpath is None:
            datpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/proc/"
        ncname  = "cesm2_pic_HMXL_regrid_bilinear.nc"
        ds_all  = xr.open_dataset(datpath+ncname)
    else:
        if datpath is None: # Indicate data path on stormtrack
            datpath = "/stormtrack/data4/glliu/01_Data/CESM2_PiControl/FCM/atm/"
        keepvars  = ["time","lat","lon",vname]
        if searchstr is None:
            searchstr = "%s%s/*%s*.nc" % (datpath,vname,vname) # Searches for datpath + *LANDFRAC*.nc
        nclist    = glob.glob(searchstr)
        nclist.sort()
        if debug:
            print("Found %i files" % (len(nclist)))
            print("\n%s" % str(nclist))
        
        ds_all    = xr.open_mfdataset(nclist,concat_dim="time",combine='nested')

    return ds_all


#%% Indicate variable information

# Indicate Dataset
datname           = "cesm2_pic"
datpath           = "/stormtrack/data4/glliu/01_Data/CESM2_PiControl/FCM/atm/"
vnames            = ["SHFLX","FSNS","FLNS","TS","LANDFRAC","ICEFRAC"]

# Output Path
outpath          = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/proc/"

# Indicate Time Crop
tstart            =  '0200-01-01' # "2006-01-01" # 
tend              =  '2000-12-31' #"2101-01-01" # 
timestr           =  '%sto%s'  % (tstart[:4],tend[:4]) # ex. 0000to2000

# Indicate Region Crop
bbox_crop         = [-90,20,0,90]
bbox_name         = "NAtl"

# Set Loader
if datname == "cesm2_pic":
    loader = load_cesm2_pic
else:
    print("Loader for %s not written!" % (datname))

#%% Load the output, looping by variable

nvars = len(vnames)

for vv in range(nvars):
    
    vname = vnames[vv]
    
    # Indicate savename <datname>_<vname>_<bbox_name>_<timestr>.nc
    savename       = "%s%s_%s_%s_%s.nc" % (outpath,datname,vname,bbox_name,timestr)
    print("Output will be saved to %s" % savename ) 
    
    # Load the ds
    ds_all   = loader(vname)[vname]
    
    # Format the ds (flip longitude, rename dimensions, time x lat x lon)
    ds_fmt   = hf.format_ds(ds_all)
    
    # Fix february start (mostly for cesm output)
    ds_fmt   = hf.fix_febstart(ds_fmt)
    
    # Crop to time
    ds_tcrop = ds_fmt.sel(time=slice(tstart,tend))
    
    # Crop to region
    ds_reg   = hf.sel_region_xr(ds_tcrop,bbox_crop)
    
    # Load the ds
    st       = time.time()
    ds_reg   = ds_reg.load()
    print("Data loaded in %.2fs" % (time.time()-st))
    
    # Save the ds
    st       = time.time()
    edict   = {vname:dict(zlib=True)}
    ds_reg.to_netcdf(savename,encoding=edict)
    print("Saved %s in %.2fs" % (vname,time.time()-st))
    print("\tOutput saved to %s " % savename)
