#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Process AVISO MSLA (downloaded 2025)

Written to crop to the same grid as process_natl_dat

Created on Thu Feb 20 10:40:33 2025

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

#%% Process MSLA Climatology (1993-2022)

dpath  = "/Users/gliu/Downloads/02_Research/01_Projects/03_SeaLevel/01_Data/00_Raw/aviso/ncfiles/msla_h_clim/"
ncfile = dpath + "dt_global_allsat_msla_h_y1993_2022_m%02i.nc"

nclist = [ncfile % (ii+1) for ii in range(12)]

dsall_sla  = xr.open_mfdataset(nclist,combine='nested',concat_dim='time')


bbox_natl = [-100,20,-10,90]

#%% MADT Climatolgoy (1993-2022)
dpath       = "/Users/gliu/Downloads/02_Research/01_Projects/03_SeaLevel/01_Data/00_Raw/aviso/ncfiles/madt_h_clim/"
ncfile      = dpath + "dt_global_allsat_madt_h_y1993_2022_m%02i.nc"
nclist      = [ncfile % (ii+1) for ii in range(12)]

dsall_adt   = xr.open_mfdataset(nclist,combine='nested',concat_dim='time')

def format_aviso(ds,crop=[-100,20,-10,90]):
    
    latname  = "latitude"
    lonname  = "longitude"
    timename = "time"
    dsa      = proc.format_ds(ds,lonname=lonname,latname=latname,timename=timename)
    
    dsareg  = proc.sel_region_xr(dsa,crop)
    return dsareg


dsareg = format_aviso(dsall_adt.adt)

outpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/NATL_proc_obs/proc/"
outname = outpath + "AVISO_adt_NAtl_1993_2022_clim.nc"
edict   = proc.make_encoding_dict(dsareg)

dsareg.to_netcdf(outname,encoding=edict)





