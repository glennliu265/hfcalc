#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Compare Yeager 2015 Box Across Different Datasets
Works with outptu from crop_spg_sst

(1) Load SST Dataset 
(2) (apply Ice Mask?)
(3) Take Average Over SPG Box
(4) Deseason and Detrend
(5) Compute Monthly ACF

Created on Wed Feb 19 11:53:16 2025

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
    

def crop_dataset(ds,ycrop):
    tstart = "%04i-01-01" % ycrop[0]
    tend   = "%04i-12-31" % ycrop[1]
    dscrop = ds.sel(time=slice(tstart,tend))
    return dscrop

def winter_acf(ts):
    ts_wint = ds.where(ts.time.dt.month.isin([12,1,2,3]),drop=True)
    
    ts_wint_ann = ts_wint.groupby('time.year').mean('time')
    
    
    return ts_wint_ann

   #ts_wint = ts.groupby('time.season')
    
    


#%% Copy the files in


#outpath     = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/03_reemergence/01_Data/proc/SPG_Box/"

dpath       = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/spg_data/"
dset_names  = ("OISST","ERA5","HadISST","ERA20C","cesm1_pic","cesm2_pic")
cols        = ("cornflowerblue","hotpink",'darkmagenta','forestgreen','gray','k')
els         = ('solid','solid','dotted','dotted','dashed','dashed')

#ystarts = ["1982","1950"]

crop_ts = True
ycrops      = ([1982,2020],[1982,2020],
              [1980,2018],[1972,2010],
              [2162,2200],[1962,2000])

ndatasets   = len(dset_names)
ncnames     = ["%s%s_SPG.nc" % (dpath,dset_names[nn]) for nn in range(ndatasets)]
ds_all      = [xr.open_dataset(nc).load() for nc in ncnames]

# Handle some time cropping
if crop_ts:
    ds_all  = [crop_dataset(ds_all[ii],ycrops[ii]) for ii in range(ndatasets)]
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
nc_adt      = dpath_aviso + "AVISO_adt_NAtl_1993_2022_clim.nc"
ds_adt      = xr.open_dataset(nc_adt).load()
cints_adt   = np.arange(-100, 110, 10)


# Load MIMOC
dpath_mimoc = dpath_ice
nc_mld      = dpath_ice + "MIMOC_RegridERA5_mld_NAtl_Climatology.nc"
ds_mld      = xr.open_dataset(nc_mld).load()

# Load Stochastic Model output from `single_point_run`
nc_sm       = "stochastic_model_point.nc"
ds_sm       = xr.open_dataset(dpath + nc_sm).load()



#%% Preprocessing

bbox_yeager = [-50,-10,50,60] # Original
#bbox_yeager = [-30,-10,50,60] # Away from Sea Ice 

bbfn,bbstr  = proc.make_locstring_bbox(bbox_yeager)

spg_anom    = [preproc_ds(ds.sst,bbox_yeager) for ds in ds_all]
spg_acf     = [calc_acf_mon(sst.data,lags) for sst in spg_anom]

#%% Compute Stochastic model acfs

#ssts_sm         = ds_sm.sst
ssts_sm_list    = [ds_sm.sst.isel(ens=ii) for ii in range(10)]

if crop_ts:
    ssts_sm_list = [ds.isel(time=slice(-(38*12)-1,-1)) for ds in ssts_sm_list]

spg_acf_sm = [calc_acf_mon(sst.data,lags) for sst in ssts_sm_list]


#%% Plot The ACF

kmonth  = 1
xtks    = np.arange(0,61,3)
fig,ax  = plt.subplots(1,1,figsize=(12,4.5))

title   = "%s SST ACF \n Box: %s" % (mons3[kmonth],bbstr)
ax,_    = viz.init_acplot(kmonth,xtks,lags,ax=ax,title=None)
ax.set_title(title,fontsize=16)

for nn in range(ndatasets):
    label = "%s (Years %s)" % (dset_names[nn],ystrs[nn])
    ax.plot(lags,spg_acf[nn][kmonth,:],label=label,
            c=cols[nn],linewidth=2.5,ls=els[nn])

# Add Stochastic Model Plots
for e in range(10):
    if e == 0:
        label="Stochastic Model"
    else:
        label=""
    ax.plot(lags,spg_acf_sm[e][kmonth,:],label=label,c='goldenrod',alpha=0.5)

#ax.legend(fontsize=12)

figname = figpath + "ACF_Comparison_Obs_mon%02i.png" % (kmonth+1)

plt.savefig(figname,dpi=150,bbox_inches='tight',transparent=True)

#%% Also Include a Locator

imon   = kmonth
proj   = ccrs.PlateCarree()
bbplot = [-80, 0, 35, 75]

# for imon in range(12):

fig, ax, _ = viz.init_orthomap(1, 1, bbplot, figsize=(14, 6))
ax = viz.add_coast_grid(ax, bbplot, fill_color='lightgray',
                        proj=proj, line_color="dimgray")

# plotvar = thflx_damping.isel(lag=ilag, month=imon)
# pcm = ax.pcolormesh(lon, lat, plotvar, transform=proj, cmap='cmo.balance',
#                     vmin=-55, vmax=55, zorder=-1)


# Plot the Box
viz.plot_box(bbox_yeager,proj=proj,color="purple",linewidth=3,linestyle='dashed')

# # Plot Sea Ice
plotvar = ds_masks.mask_mon
cl = ax.contour(plotvar.lon, plotvar.lat,
                plotvar, colors="cyan",
                linewidths=2, transform=proj, levels=[0, 1], zorder=-1)
ax.clabel(cl, fontsize=12)

plotvar = ds_adt.isel(time=imon)
cl = ax.contour(plotvar.lon, plotvar.lat, plotvar.adt*100, colors="k",
                linewidths=0.75, transform=proj, levels=cints_adt)
ax.clabel(cl)



ax.set_title("Bounding Box: %s" % (bbstr), fontsize=18,y=1.05)
#cb = viz.hcbar(pcm, ax=ax, pad=0.01, fraction=0.045)
# cb.set_label(
#     "Turbulent Heat Flux Feedback (W m$^{-2}$ $\degree$C$^{-1}$)", fontsize=14)
# cb.ax.tick_params(labelsize=14)

outname = figpath + "Bounding_Box_%s.png" % (bbfn)
plt.savefig(outname, dpi=150, bbox_inches='tight', transparent=True)

#%% Plot the actual anomaly timeseries

fig,ax = plt.subplots(1,1,figsize=(8,3.5),
                      constrained_layout=True)


times  = spg_anom[1].time

for nn in range(ndatasets-2):
    ts = spg_anom[nn]
    ax.plot(ts.time,ts,label=dset_names[nn],c=cols[nn],ls=els[nn],lw=2)
ax.legend(ncol=4)

ax.set_xlim([times[0],times[-1]])

outname = figpath + "SPG_Timeseries_%s.png" % (bbfn)
plt.savefig(outname, dpi=150, bbox_inches='tight', transparent=True)

#%% NOTE OLD SCRIPT BELOW =-=========

#%% Do Calculation for OISST --------------------------------------------------

bbox_yeager = [-50,-10,50,60] # Original
bbox_yeager = [-50,-10,50,55] # Shift South
bbox_yeager = [-50,-30,50,60] # Shift West
bbox_yeager = [-30,-10,50,60] # Shift East

bbox_spg    = [-80,0,30,70] # SPG Crop for all Datasets


# Here is the file
dpath_proc  = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/NATL_proc_obs/proc/"
nc_sst      = dpath_proc + "OISST_sst_NAtl_1982to2020.nc"
ds          = xr.open_dataset(nc_sst).load()
spg_ts      = preproc_ds(ds.sst,bbin=bbox_yeager)

spg_ts      = spg_ts.rename("OISST_SST")

outname     = outpath + "OISST_SPG_Box_SST.nc"

spg_ts.to_netcdf(outname)

#%% Compute the autocorrelation function

lags        = np.arange(0,60)
sst_in      = spg_ts.data
acfs        = calc_acf_mon(spg_ts.data,lags)
outname     = outpath + "OISST_SPG_BOX_ACF.npy"

np.save(outname,acfs)

#%% Do the same for CESM1 and CESM2 PIC

# CESM1 PIC
ncname = "TS_PIC_FULL.nc"
ncpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/01_Data/CESM_proc/"
ds     = xr.open_dataset(ncpath + ncname)

#%% Visualize the ACF....
kmonth = 1

xtks   = np.arange(0,61,3)

fig,ax=plt.subplots(1,1,figsize=(12,4.5))
ax,_ = viz.init_acplot(kmonth,xtks,lags,ax=ax)

ax.plot(lags,acfs[kmonth,:],c='navy',label="OISST")
    
#%% Copied from check_hfdamping

dpath         = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/" 
ncname_obs    = "ERA5_RegridOISST_hfdamping_NAtl_1982to2020_ensorem1_detrend1.nc"
ds            = xr.open_dataset(dpath + ncname_obs)

dt            = 3600*24*30


lon           = ds.lon
lat           = ds.lat
thflx_damping = ds.thflx_damping / dt

cov           = ds.cov
autocov       = ds.autocov

cc            = ds.sst_flx_crosscorr
ac            = ds.sst_autocorr



dpath_proc      = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/NATL_proc_obs/proc/"

nc_sst          = dpath_proc + "OISST_sst_NAtl_1982to2020.nc"
nc_flx          = dpath_proc + "ERA5_RegridOISST_thflx_NAtl_1982to2020.nc"
nc_flx_raw      = dpath_proc + "ERA5_thflx_NAtl_1982to2020.nc"

ds_sst          = xr.open_dataset(nc_sst)
ds_flx_regrid   = xr.open_dataset(nc_flx)
ds_flx_raw      = xr.open_dataset(nc_flx_raw)

nc_adt          = dpath_proc + "AVISO_adt_NAtl_1993_2022_clim.nc"
ds_adt          = xr.open_dataset(nc_adt)


#%% Plot Autocorr


bbin = bbox_yeager
ilag = 0
imon = 1
mons3 = proc.get_monstr()

cints_adt = np.arange(-100,110,10)

proj      = ccrs.PlateCarree()
bbplot    = [-80,0,35,75]
fig,ax,_  = viz.init_orthomap(1,1,bbplot,figsize=(14,6))
ax        = viz.add_coast_grid(ax,bbplot,fill_color='k',proj=proj)

plotvar   = ac.isel(lag=ilag,month=imon) # 
pcm       = ax.pcolormesh(lon,lat,plotvar,transform=proj,cmap='cmo.balance',
                         vmin=-1,vmax=1)

plotvar   = ds_adt.isel(time=imon)
cl        = ax.contour(plotvar.lon,plotvar.lat,plotvar.adt*100,colors="k",
                      linewidths=0.75,transform=proj,levels=cints_adt)
ax.clabel(cl)

viz.plot_box(bbin,proj=proj,color="yellow",linewidth=1.5)

ax.set_title("Month %s, Lag %i \n%s" % (mons3[imon],ilag+1,bbin),fontsize=14)
cb = fig.colorbar(pcm,ax=ax,pad=0.01)
cb.set_label("Lag %i SST Autocorrelation" % (ilag+1))


