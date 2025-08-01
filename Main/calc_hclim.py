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


# Indicate Time Crop
tstart            =  '1920-01-01' #'0200-01-01' # "2006-01-01" # 
tend              =  '2005-12-31' #'2000-12-31' #"2101-01-01" # 
timestr           =  '%sto%s'  % (tstart[:4],tend[:4]) # ex. 0000to2000
preproc_name      = 'cesm1_htr_5degbilinear' #cesm2_regcrop"



# Set to True to save stochastic model output
sm_input          = True # Set to True to prepare for stochastic model
sm_bbox           = [-80,0,0,65]
sm_bbox_name      = "NAtl"

# Location to check if input is in meters
lonf        = -30
latf        = 50
chk_thres   = 1000 # If greater than this, convert to meters.

# Indicate Dataset
if preproc_name == "cesm2_prep_mld_pic":
    # Direct output from prep_MLD_PIC.py for CESM2
    datname       = "cesm2_pic"
    datpath       = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/proc/"
    ncname        = "cesm2_pic_HMXL_regrid_bilinear.nc"
    vname         = "HMXL"
    
    # Indicate savename
    savename       = "%s%s_%s_%s_scycle.nc" % (datpath,datname,vname,timestr)
    
elif preproc_name == "cesm2_regcrop":
    # Output from [preproc_raw_inputs]
    datname       = "cesm2_pic"
    datpath       = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/proc/"
    ncname        = "cesm2_pic_HMXL_NAtl_0200to2000.nc"
    vname         = "HMXL"
    
    # Indicate savename
    savename       = datpath + hf.addstrtoext(ncname,"_scycle",adjust=-1)#"%scesm2_pic_HMXL_NAtl_0200to2000_scycle.nc" % (datpath)
    
elif preproc_name == "cesm1_htr_5degbilinear":
    # Output from [preproc_raw_inputs]
    datname       = preproc_name#"cesm2_pic"
    datpath       = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/proc/"
    ncname        = "cesm1_htr_5degbilinear_HMXL_Global_1920to2005.nc"
    vname         = "HMXL"
    
    # Indicate savename
    savename       = datpath + hf.addstrtoext(ncname,"_scycle",adjust=-1)#"%scesm2_pic_HMXL_NAtl_0200to2000_scycle.nc" % (datpath)

print("Output will be saved to %s" % savename ) 

# Stochastic model input options
vname_out         = "h"
input_path        = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/03_reemergence/proc/model_input/mld/"
#sm_savename       = "%s%s_%s_NAtl_0200to2000.nc" % (input_path,datname,vname)
sm_savename       = "%s%s_%s_%s_1920to2005.nc" % (input_path,datname,vname,sm_bbox_name)
print("Stochastic Model output will be saved to %s" % (sm_savename))

#%% Load and Preprocess 

# Open ds
ds              = xr.open_dataset(datpath+ncname)[vname]

# Preprocess ds
if np.any(ds.lon > 180):
    print("Flipping Lon!")
    #ds              = hf.format_ds(ds) # This doesn't work because of ens dimension
    ds = hf.lon360to180_xr(ds)

# Fix February Start
ds              = hf.fix_febstart(ds)

# Crop the region
ds              = hf.sel_region_xr(ds,sm_bbox)

# Crop the Time Period
dscrop          = ds.sel(time=slice(tstart,tend)).load()

# Calculate seasonal mean
ds_scycle       = dscrop.groupby('time.month').mean('time')

# Convert cm --> m/s
if ds_scycle.mean() > chk_thres:
    print("Input in [Centimeters] since mean is %.2f, converting to meters..." % (ds_scycle.mean().item()))
    ds_scycle_meter = ds_scycle/100
else:
    ds_scycle_meter = ds_scycle

# do a check

# Save Output
edict = {vname:dict(zlib=True)}
ds_scycle_meter.to_netcdf(savename,encoding=edict)

#%% Save Stochastic Model Output

# Reload things
ds_load = xr.open_dataset(savename)[vname].load()

# Crop to Region

# Rename Dimensions
ds_load = ds_load.rename(vname_out)
ds_load = ds_load.rename(dict(month='mon'))
edict   = {'h':dict(zlib=True)}
ds_load.to_netcdf(sm_savename,encoding=edict)




