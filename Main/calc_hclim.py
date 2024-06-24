#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Compute Seasonal Mixed-Layer Depth Cycle (Global)

- Works with bilinear regridded output from stochmod/prep_MLD_PIC.py
- Currently supports CESM2... might need to rewrite for other cases

Procedure
    - Load the data
    - Fix febstart
    - Restrict to time period (need to test sensitivity)
    - Convert cm --> m
    - Compute mean seasonal cycle
    - Save Global Output


 
 

Created on Fri Jun 21 14:15:58 2024

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

# import amv.hfutils as hf

# -------------
#%% User Edits
# -------------


# Indicate Dataset
datname       = "cesm2_pic"
datpath       = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/proc/"
ncname        = "cesm2_pic_HMXL_regrid_bilinear.nc"
vname         = "HMXL"
                             
# Indicate Time Crop
tstart            =  '0200-01-01' # "2006-01-01" # 
tend              =  '2000-12-31' #"2101-01-01" # 
timestr           =  '%sto%s'  % (tstart[:4],tend[:4]) # ex. 0000to2000


# Indicate savename
savename       = "%s%s_%s_%s_scycle.nc" % (datpath,datname,vname,timestr)
print("Output will be saved to %s" % savename ) 

#%% Load and Preprocess 

# Open ds
ds              = xr.open_dataset(datpath+ncname)[vname]

# Fix February Start
ds              = hf.fix_febstart(ds)

# Crop the Time Period
dscrop          = ds.sel(time=slice(tstart,tend)).load()

# Calculate seasonal mean
ds_scycle       = dscrop.groupby('time.month').mean('time')

# Convert cm --> m/s
ds_scycle_meter = ds_scycle/100

# Save Output
edict = {vname:dict(zlib=True)}
ds_scycle_meter.to_netcdf(savename,encoding=edict)
