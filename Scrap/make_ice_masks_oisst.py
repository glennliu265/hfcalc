#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Make Sea Ice Masks and visualize Sea Ice Concentration from OISST Data

Uses output from [process_crop_data.py]

Created on Tue Feb 25 18:09:48 2025

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
import matplotlib as mpl
from cmcrameri import cm


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

#%% User Edits

# Plot Settings
mpl.rcParams['font.family']     = 'Avenir'
proj                            = ccrs.PlateCarree()
bbplot                          = [-80,0,35,75]

# Paths
figpath                         = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/03_reemergence/02_Figures/20250221/"
proc.makedir(figpath)

#%% Load and check Sea Ice

# Daily Max Sea Ice
dpath_proc          = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/NATL_proc_obs/proc/"
nc_sicemax          = dpath_proc + "../OISST_icec_max_1981_2020.nc"
ds_sice             = xr.open_dataset(nc_sicemax)
ds_sice_natl         = proc.sel_region_xr(ds_sice,[-100,20,10,90])


# Monthly Mean Sea ice
nc_sice_monmean             = dpath_proc + "../OISST_icec_monmean_1981_2020.nc"
ds_sice_monmean             = xr.open_dataset(nc_sice_monmean).load()
ds_sice_monmean             = proc.lon360to180_xr(ds_sice_monmean)
ds_sice_monmean_natl        = proc.sel_region_xr(ds_sice_monmean,[-100,20,10,90])
ds_sice_monmean_natl_max    = ds_sice_monmean_natl.icec.mean('month').max('year')

# Sea Ice Exceedence Frequency (Daily)
nc_sice_exceed            = dpath_proc + "../OISST_iceexceed005_monmean_1981_2020.nc"
ds_sice_exceed            = xr.open_dataset(nc_sice_exceed).load()
ds_sice_exceed            = proc.lon360to180_xr(ds_sice_exceed)
ds_sice_exceed_natl       = proc.sel_region_xr(ds_sice_exceed,[-100,20,10,90])





#%% Compute and make ice masks

# Sea Ice Exceedence Frequency (Daily)
exceed_daily              = ds_sice_exceed_natl.icec.sum('time')


# Sea Ice Exceedence Frequency (Monthly)
exceed_mon = ds_sice_monmean_natl.icec > 0.05
exceed_mon = exceed_mon.sum('month').sum('year')

# Get Total Counts
max_count_daily           = 365 * len(ds_sice_monmean_natl.year) + len(ds_sice_monmean_natl.year)/4 # Account for leap years
max_count_monthly         = 12  *len (ds_sice_monmean_natl.year)

#
mask_exceed_daily = (exceed_daily/max_count_daily) < 0.05
mask_exceed_mon   = (exceed_mon/max_count_monthly) < 0.05

ds_masks = xr.merge([mask_exceed_daily.rename("mask_daily"),
                     mask_exceed_mon.rename("mask_mon")
                     ])
edict             = proc.make_encoding_dict(ds_masks)
outname           = dpath_proc + "../OISST_ice_masks_1981_2020.nc"
ds_masks.to_netcdf(outname,encoding=edict)

#%% Plot Sea Ice (Daily Max)


# Daily Max
plotvar   = ds_sice_natl.icec.max('time')
# Monthly Mean
#plotvar = ds_sice_monmean_natl.icec.mean('month').max('year')


# Make the PLot
fig,ax,_  = viz.init_orthomap(1,1,bbplot,figsize=(14,6))
ax        = viz.add_coast_grid(ax,bbplot,fill_color='lightgray',proj=proj)


pcm       = ax.pcolormesh(ds_sice_natl.lon,ds_sice_natl.lat,
                          plotvar,transform=proj,cmap='cmo.ice',
                         vmin=0,vmax=1)

# plotvar   = ds_adt.isel(time=imon)
cl        = ax.contour(ds_sice_natl.lon,ds_sice_natl.lat,
                          plotvar,colors="forestgreen",
                      linewidths=1.5,transform=proj,levels=[0.05,])
ax.clabel(cl,fontsize=14)

#ax.set_title("Month %s, Lag %i" % (mons3[imon],ilag+1),fontsize=14)
cb = viz.hcbar(pcm,ax=ax,pad=0.01,fraction=0.045)
cb.set_label("Max Sea Ice Concentration (OISST, 1981 to 2020)",fontsize=14)
cb.ax.tick_params(labelsize=14)


outname= figpath + "OISST_Max_SeaIceConcentration.png"
plt.savefig(outname,dpi=150,bbox_inches='tight')

#%% Compare Daily and Monthly Max


fig,axs,_  = viz.init_orthomap(1,2,bbplot,figsize=(14,6))

for a,ax in enumerate(axs):
    ax        = viz.add_coast_grid(ax,bbplot,fill_color='lightgray',proj=proj)
    
    if a == 0:
        title   = "Daily Resolution"
        plotvar = ds_sice_natl.icec.max('time')
    else:
        title="Monthly Resolution"
        plotvar = ds_sice_monmean_natl_max
        
    
    ax.set_title(title,fontsize=24,y=1.05)
    
    pcm       = ax.pcolormesh(ds_sice_natl.lon,ds_sice_natl.lat,
                              plotvar,transform=proj,cmap='cmo.ice',
                             vmin=0,vmax=1)



    
    
    cl        = ax.contour(ds_sice_natl.lon,ds_sice_natl.lat,
                              plotvar>0.05,colors="forestgreen",
                          linewidths=1.5,transform=proj,levels=[0,1])
    ax.clabel(cl,fontsize=14)

#ax.set_title("Month %s, Lag %i" % (mons3[imon],ilag+1),fontsize=14)
cb = viz.hcbar(pcm,ax=axs.flatten(),pad=0.01,fraction=0.045)
cb.set_label("Max Sea Ice Concentration (OISST, 1981 to 2020)",fontsize=24)
cb.ax.tick_params(labelsize=16)   
    
outname= figpath + "OISST_Max_SeaIceConcentration_DailyvMonthly.png"
plt.savefig(outname,dpi=150,bbox_inches='tight')

#%% Plot Frequency of Exceedence

fig,axs,_  = viz.init_orthomap(1,2,bbplot,figsize=(14,6))

for a,ax in enumerate(axs):
    ax        = viz.add_coast_grid(ax,bbplot,fill_color='lightgray',proj=proj)
    
    if a == 0:
        title     = "Daily Resolution"
        plotvar   = exceed_daily / max_count_daily
        mask_plot = ds_masks.mask_daily
    else:
        title="Monthly Resolution"
        plotvar = exceed_mon / max_count_monthly
        mask_plot = ds_masks.mask_mon
    
    ax.set_title(title,fontsize=24,y=1.05)
    
    pcm       = ax.pcolormesh(ds_sice_natl.lon,ds_sice_natl.lat,
                              plotvar,transform=proj,cmap='cmo.ice',
                             vmin=0,vmax=1,zorder=-1)

    cl        = ax.contour(ds_sice_natl.lon,ds_sice_natl.lat,
                              plotvar>0.05,colors="forestgreen",
                          linewidths=1.5,transform=proj,levels=[0,1])
    ax.clabel(cl,fontsize=14)
    
    
    cl2        = ax.contour(ds_sice_natl.lon,ds_sice_natl.lat,
                              mask_plot,colors="yellow",linestyles='dashed',
                          linewidths=1,transform=proj,levels=[0,1])

#ax.set_title("Month %s, Lag %i" % (mons3[imon],ilag+1),fontsize=14)
cb = viz.hcbar(pcm,ax=axs.flatten(),pad=0.01,fraction=0.045)
cb.set_label("% of time exceeding 5% Ice Concentration",fontsize=24)
cb.ax.tick_params(labelsize=16)   
    
outname= figpath + "OISST_Max_siceExceed.png"
plt.savefig(outname,dpi=150,bbox_inches='tight')

