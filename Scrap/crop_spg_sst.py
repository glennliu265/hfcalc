#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Crop SST Datasets to SPG region

Created on Wed Feb 26 09:03:21 2025

@author: gliu

"""


import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import time
import sys
import cartopy.crs as ccrs
import glob

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

#%% Other Settings

outpath     = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/03_reemergence/01_Data/proc/SPG_Box/"
bbox_spg    = [-80,0,30,70] # SPG Crop for all Datasets

#%%  Preprocessing Section
"""
Lets Load each dataset and crop to this region

Rules

[sst] [(ens) x time x lat x lon180], Units: degree Celsius

"""
#%% Start with OISST

outpath_spg = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/spg_data/"

# OISST
dpath_proc  = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/NATL_proc_obs/proc/"
nc_sst      = dpath_proc + "OISST_sst_NAtl_1982to2020.nc"
ds          = xr.open_dataset(nc_sst).load() # [Time x Lon x Lat]
ds_sst      = ds.sst.transpose('time','lat','lon')
ds_sst      = proc.sel_region_xr(ds_sst,bbox_spg)

outname     = outpath_spg + "OISST_SPG.nc"
edict       = proc.make_encoding_dict(ds_sst)
ds_sst.to_netcdf(outname,encoding=edict)

#%% Next, HadiSST
dpath  = "/Users/gliu/Downloads/02_Research/01_Projects/04_Predict_AMV/03_Scripts/CESM_data/"
ncname = dpath + "hadisst.1870-01-01_2018-12-01.nc"
ds     = xr.open_dataset(ncname)

ds_sst      = ds.sst.transpose('time','lat','lon')
ds_sst      = proc.sel_region_xr(ds_sst,bbox_spg)

outname     = outpath_spg + "HadISST_SPG.nc"
edict       = proc.make_encoding_dict(ds_sst)
ds_sst.to_netcdf(outname,encoding=edict)

#%% Next, ERA20C
dpath  = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/"
ncname = dpath + "era20c_ts.nc"
ds     = xr.open_dataset(ncname).load()

# Flip at long and rename
ds     = proc.format_ds(ds)
ds     = ds.rename({'ts':'sst'})

# Convert to Celsius
if np.any(ds > 200):
    print("Converting to Celsius")
    ds = ds - 273.15

ds_sst      = ds.sst.transpose('time','lat','lon')
ds_sst      = proc.sel_region_xr(ds_sst,bbox_spg)

outname     = outpath_spg + "ERA20C_SPG.nc"
edict       = proc.make_encoding_dict(ds_sst)
ds_sst.to_netcdf(outname,encoding=edict)


#%% Next, CESM2 (works with output from hfcalc/scrap/process_crop_data)

# # Craop this is yearly data :(
# dpath  = "/Users/gliu/MIT Dropbox/Glenn Liu/Glenn Liu’s files/Home/CMIP6_LENS_Datasets/"
# ncname = dpath + "CESM2_sst_NAtl_1850to2014_detrend0_regridNonedeg.nc"

dpath   = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/03_reemergence/01_Data/proc/CESM1/NATL_proc/"
ncname  = dpath + "cesm2_pic_TS_NAtl_0001_2000.nc"
ds      = xr.open_dataset(ncname).load()

# Flip at long and rename
ds     = proc.format_ds(ds)
ds     = ds.rename({'TS':'sst'})


# Convert to Celsius
if np.any(ds > 200):
    print("Converting to Celsius")
    ds = ds - 273.15
    

ds_sst      = ds.sst.transpose('time','lat','lon')
ds_sst      = proc.sel_region_xr(ds_sst,bbox_spg)

outname     = outpath_spg + "cesm2_pic_SPG.nc"
edict       = proc.make_encoding_dict(ds_sst)
ds_sst.to_netcdf(outname,encoding=edict)

# Need to Redo Above...
#%% Do the same for CESM1 PIC

dpath    = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/01_Data/CESM_proc/"
ncname   = dpath + "TS_PIC_FULL.nc"
ds       = xr.open_dataset(ncname)

ds       = proc.lon360to180_xr(ds)
dsreg    = proc.sel_region_xr(ds,bbox_spg)
ts_reg   = dsreg.TS.data # yr x mon x lat x lon

nyr,nmon,nlat,nlon  = ts_reg.shape
ts_reg              = ts_reg.reshape(nyr*nmon,nlat,nlon)

ts_reg = ts_reg - 273.15 # Convert to Celsius

times   = xr.cftime_range(start='0400',periods=1801*12,freq="MS",calendar="noleap")
dims    = {'time':times,'lat':dsreg.lat,'lon':dsreg.lon}
ds_sst  = xr.DataArray(ts_reg,dims=dims,coords=dims,name="sst")

outname     = outpath_spg + "cesm1_pic_SPG.nc"
edict       = proc.make_encoding_dict(ds_sst)
ds_sst.to_netcdf(outname,encoding=edict)




