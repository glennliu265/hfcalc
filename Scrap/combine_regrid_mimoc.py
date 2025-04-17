#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Merge MIMOC Climatology and Regrid to XCESMF

Created on Mon Mar  3 13:15:03 2025

@author: gliu
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import time
import sys
import cartopy.crs as ccrs
import glob

#import pandas as pd
#import scipy as sp
import xesmf as xe


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


#%% Load and Merge MIMOC (copied from viz_mldvar.py)

mldpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/01_Data/MIMOC_ML_v2.2_PT_S/"
#testpath = mldpath + "MIMOC_ML_v2.2_PT_S_MLP_month01.nc"
nclist      = glob.glob(mldpath+"*.nc")
nclist.sort()
print(nclist)

# Read in and concatenate by month variable
ds_all = []
for nc in nclist:
    ds = xr.open_dataset(nc)
    print(ds)
    ds_all.append(ds.DEPTH_MIXED_LAYER)

ds_all          = xr.concat(ds_all,dim="month")

lon             = ds.LONGITUDE.data
lat             = ds.LATITUDE.data

ds_all['LAT']   = lat
ds_all['LONG']  = lon

#%% Save Reformatted Mixed-Layer

latname         = "LAT"
lonname         = "LONG"
timename        = 'month'

ds_rfm          = proc.format_ds(ds_all,timename='month',latname=latname,lonname=lonname)

ds_rfm          = ds_rfm.rename({'time':'month'})
ds_rfm['month'] = np.arange(1,13,1) # Change from 0-11 to 1-12

#%% Check Latitude

lons = ds_rfm['lon'].data
print(lons[1:] - lons[:-1])
#ds_rfm['lon'].data[1:] - ds_rfm['lon'].data[:-1]

#%%
outname = ""
outpath = ""

#%% Load ERA5 Grid

figpath    = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/03_reemergence/02_Figures/20250305/"

dpath_proc = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/NATL_proc_obs/"
nc_era5    = "ERA5_sst_NAtl_1979to2021.nc"
ds_era5    = xr.open_dataset(dpath_proc + nc_era5).load()

lat_e = ds_era5.lat
lon_e = ds_era5.lon

# Regrid (Copied form predict_nasst/regrid_ocean_variable_hmxl.py)
method    = 'bilinear'
ds_out    = xr.Dataset({'lat':lat_e,"lon":lon_e})
ds        = ds_rfm

# Initialize Regridder
regridder = xe.Regridder(ds,ds_out,method,periodic=False)

# Regrid
daproc    = regridder(ds) # Need to input dataarray

daproc    = daproc.rename("mld")

edict     = proc.make_encoding_dict(daproc)
outname   = dpath_proc + "MIMOC_RegridERA5_mld_NAtl_Climatology.nc"

daproc.to_netcdf(outname,encoding=edict)


#%% Plot Sample Plot

bbplot  = [-80, 0, 35, 75]
proj    = ccrs.PlateCarree()
imon    = 1
mons3   = proc.get_monstr()
cints   = np.arange(0,550,50)

fig, axs, _ = viz.init_orthomap(1, 2, bbplot, figsize=(14, 6))

for a,ax in enumerate(axs):
    ax = viz.add_coast_grid(ax, bbplot, fill_color='lightgray',
                            proj=proj, line_color="dimgray")
    
    if a == 0:
        title   = "Original (0.5$\degree$)"
        plotvar = ds_rfm.isel(month=imon)
        
    elif a == 1:
        title = "Regrid ERA5 (0.25$\degree$)"
        plotvar = daproc.isel(month=imon)
    
    ax.set_title(title)
    
    pcm = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,transform=proj,vmin=0,vmax=500,cmap='cmo.deep')
    cl  = ax.contour(plotvar.lon,plotvar.lat,plotvar,transform=proj,levels=cints,
                     colors="k",linewidths=0.75)
    ax.clabel(cl)
    
cb = viz.hcbar(pcm,ax=axs.flatten(),pad=0.01)
cb.set_label("%s Mixed Layer Depth [m]" % mons3[imon])

#plt.show()

savename = "%sMIMOC_MLD_ERA5_Regrid_Comparison.png" % figpath
plt.savefig(savename,dpi=150,bbox_inches='tight',transparent=True)

#%% Impact on Yeager Box

bbox_yeager    = [-50,-10,50,60] # Original
dsreg_original = proc.sel_region_xr(ds_rfm,bbox_yeager)
dsreg_regrid   = proc.sel_region_xr(daproc,bbox_yeager)


aavg_original  = proc.area_avg_cosweight(dsreg_original,)
aavg_regrid    = proc.area_avg_cosweight(dsreg_regrid)


fig,ax=viz.init_monplot(1,1,)

ax.plot(mons3,aavg_original,color="k",label="Original")
ax.plot(mons3,aavg_regrid,color="red",ls='dashed',label="Regrid")
ax.legend()
savename = "%sMIMOC_MLD_ERA5_Regrid_Comparison_YeagerBoxAvg.png" % figpath
plt.savefig(savename,dpi=150,bbox_inches='tight',transparent=True)
