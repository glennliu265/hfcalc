#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 11:28:45 2025

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

#%% Functions

def preproc_ds(ds_in,bbin):
    
    # Select Region
    dsreg = proc.sel_region_xr(ds_in,bbin)
    
    # Take Area Weighted Average
    spg_avg = proc.area_avg_cosweight(dsreg,sqrt=False)
    
    # Deseason
    spg_anom = spg_avg.groupby('time.month') - spg_avg.groupby('time.month').mean('time')
    
    # Detrend (simple linear)
    spg_dt   = sp.signal.detrend(spg_anom,0,type='linear')#proc.xrdetrend(spg_anom,timename='time',verbose=False)
    indict   = {'time':ds_in.time}
    spg_dt   = xr.DataArray(spg_dt,dims=indict,coords=indict)
    return spg_dt

def calc_acf_mon(sst,lags):
    
    nmon        = int(len(sst)/12)
    acfs        = np.zeros((12,len(lags))) * np.nan
    
    sst_in_mon  = sst.reshape(nmon,12).T # [mon x year]
    
    for im in range(12):
        acf = proc.calc_lagcovar(sst_in_mon,sst_in_mon,lags,im+1,1,debug=False)
        acfs[im,:] = acf.copy()
        
    return acfs

def get_year_range(ds):
    times  = ds.time
    ystart = str(times[0].data)[:4]
    yend   = str(times[-1].data)[:4]
    ystr   = "%s to %s" % (ystart,yend)
    return ystr

def winter_acf(ts):
    ts_wint = ds.where(ts.time.dt.month.isin([12,1,2,3]),drop=True)
    
    ts_wint_ann = ts_wint.groupby('time.year').mean('time')
    
    
    return ts_wint_ann

   #ts_wint = ts.groupby('time.season')
    
    


#%% Copy the files in

#outpath     = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/03_reemergence/01_Data/proc/SPG_Box/"

dpath       = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/spg_data/"
dset_names  = ("OISST","HadISST","ERA20C","cesm1_pic","cesm2_pic")
cols        = ("cornflowerblue",'salmon','goldenrod','gray','k')

#ystarts = ["1982","1950"]

ndatasets   = len(dset_names)
ncnames     = ["%s%s_SPG.nc" % (dpath,dset_names[nn]) for nn in range(ndatasets)]
ds_all      = [xr.open_dataset(nc).load() for nc in ncnames]
ntimes      = [len(ds.time) for ds in ds_all]


ystrs       = [get_year_range(ds) for ds in ds_all]

lags  =  np.arange(61)
mons3 = proc.get_monstr()

# Paths
figpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/03_reemergence/02_Figures/20250221/"
proc.makedir(figpath)


# Load SSh and Sea Ice fraction for viz
# Load Sea Ice Masks
dpath_ice = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/NATL_proc_obs/"
nc_masks = dpath_ice + "OISST_ice_masks_1981_2020.nc"
ds_masks = xr.open_dataset(nc_masks).load()

# Load AVISO
dpath_aviso = dpath_ice + "proc/"
nc_adt = dpath_aviso + "AVISO_adt_NAtl_1993_2022_clim.nc"
ds_adt = xr.open_dataset(nc_adt).load()
cints_adt = np.arange(-100, 110, 10)


#%% Check OISST Ice Coverage

ds = ds_all[0]

(ds.sst.isel(time=0)==-1.8).plot()




#%% Check Sea Ice Point


itime = 2

ssts_slice = ds.sst.isel(time=itime).data.flatten()

#bins = np.arange(-2.0,26,.5)
bins = np.arange(-2,2.1,0.1)
fig,ax = plt.subplots(1,1,figsize=(8,3))
ax.hist(ssts_slice,bins=bins,edgecolor="w")
ax.grid(True,ls='dotted')

ax.axvline([-1.8],c='k')
title = "SPG SST Distribution (-2 to 2 $\degree$C) @ %s" % (ds.time.isel(time=itime).data)
ax.set_title(title)


#%%

plotvar = ds.sst.isel(time=itime)

proj   = ccrs.PlateCarree()
bbplot = [-80, 0, 35, 75]
bbplot2 = [-70,-10,45,70]
cints  = np.arange(-2,0.2,0.2)
# for imon in range(12):

fig, ax, _ = viz.init_orthomap(1, 1, bbplot2, figsize=(14, 6))
ax = viz.add_coast_grid(ax, bbplot2, fill_color='lightgray',
                        proj=proj, line_color="dimgray")

pcm = plt.pcolormesh(plotvar.lon,plotvar.lat,plotvar < -1.8,transform=proj,cmap='cmo.thermal')
cb  = viz.hcbar(pcm,pad=0.01)

cl = ax.contour(plotvar.lon,plotvar.lat,plotvar,
                levels=cints,transform=proj,colors="w",lw=0.75)
ax.clabel(cl)

title = "SST @ %s" % (ds.time.isel(time=itime).data)
ax.set_title(title)

#%% Plot Western Side

proj    = ccrs.PlateCarree()
bbplot  = [-80, 0, 35, 75]
bbplot2 = [-70,-10,50,70]
cints   = np.arange(-2,0.2,0.2)
# for imon in range(12):

fig, ax, _  = viz.init_orthomap(1, 1, bbplot2, figsize=(14, 6))
ax          = viz.add_coast_grid(ax, bbplot, fill_color='lightgray',
                        proj=proj, line_color="dimgray")

pcm         = plt.pcolormesh(plotvar.lon,plotvar.lat,plotvar,transform=proj,cmap='cmo.thermal')
cb          = viz.hcbar(pcm)

cl          = ax.contour(plotvar.lon,plotvar.lat,plotvar,
                levels=cints,transform=proj,colors="w",lw=0.75)
ax.clabel(cl)


#%% Compute the ACF

