#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Check ACF

Created on Thu Feb 27 14:23:28 2025

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
dset_names  = ("OISST","HadISST","ERA20C","c-esm1_pic","cesm2_pic")
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


#%% Load the ACF Computed with pointwise_crosscorrelation

dpath2  = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/NATL_proc_obs/"
nc      =  dpath2 + "OISST_1982_2020_ACF_lag00to60_ALL_ensALL.nc"
ds      = xr.open_dataset(nc).load()

#%% Select a point

lonf = -30
latf = 50
imon = 1

acfpt = ds.acf.sel(lon=lonf,lat=latf,method='nearest').isel(ens=0,thres=0)

acfpt.isel(mons=imon).plot()


#%% Get the SST at the same point and plot

lags    = np.arange(61)
sst_pt  = ds_all[0].sst.sel(lon=lonf,lat=latf,method='nearest').load()
ssta    = proc.xrdeseason(sst_pt)
ssta_dt = sp.signal.detrend(ssta)

ntime   = len(ssta_dt)
nyr     = int(ntime/12)


sst_yrmon = ssta_dt.reshape(nyr,12).T # Transpose to month x year

out     = proc.calc_lagcovar(sst_yrmon,sst_yrmon,lags,imon+1,0,debug=False)

#%% compare the 2

fig,ax = plt.subplots(1,1)
ax.plot(lags,acfpt[imon,:],label="Cross Corr Script")
ax.plot(lags,out,label="Manual",ls='dashed')
ax.legend()


#%% Integrate the Lag

acfs = ds.isel(thres=0,ens=0).acf.squeeze()

t2 = proc.calc_T2(acfs,axis=-1,ds=True)


#%% Visualize the T2

imon = 0

bbox_yeager = [-50,-10,50,60] # Original


proj = ccrs.PlateCarree()
bbplot = [-80, 0, 35, 75]

for imon in range(12):
    fig, ax, _ = viz.init_orthomap(1, 1, bbplot, figsize=(14, 6))
    ax = viz.add_coast_grid(ax, bbplot, fill_color='lightgray',
                            proj=proj, line_color="dimgray")
    
    plotvar = t2[:,:,imon]
    lon = ds.lon
    lat = ds.lat
    pcm = ax.pcolormesh(lon, lat, plotvar.T, transform=proj, cmap='cmo.tempo',
                        vmin=0, vmax=24, zorder=-1)
    
    # # Plot Sea Ice
    plotvar = ds_masks.mask_mon
    cl = ax.contour(plotvar.lon, plotvar.lat,
                    plotvar, colors="yellow",
                    linewidths=2, transform=proj, levels=[0, 1], zorder=-1)
    ax.clabel(cl, fontsize=12)
    
    plotvar = ds_adt.isel(time=imon)
    cl = ax.contour(plotvar.lon, plotvar.lat, plotvar.adt*100, colors="k",
                    linewidths=0.75, transform=proj, levels=cints_adt)
    ax.clabel(cl)
    
    ax.set_title("$T^2$ Monthly, %s" % (mons3[imon]), fontsize=16)
    cb = viz.hcbar(pcm, ax=ax, pad=0.01, fraction=0.045)
    cb.set_label(
        "$T^2 (Months)$", fontsize=14)
    cb.ax.tick_params(labelsize=14)
    
    # Plot the Box
    viz.plot_box(bbox_yeager,proj=proj,color="purple",linewidth=3,linestyle='dashed')
    
    outname = figpath + "OISST_T2_lagmax%02i_mon%02i.png" % (acfs.lags.data[-1], imon+1)
    plt.savefig(outname, dpi=150, bbox_inches='tight', transparent=True)


#%% Take the wintertime Mean

mon_sel = [0,1,2]
fig, ax, _ = viz.init_orthomap(1, 1, bbplot, figsize=(14, 6))
ax = viz.add_coast_grid(ax, bbplot, fill_color='lightgray',
                        proj=proj, line_color="dimgray")


plotvar = t2[:,:,mon_sel].mean(-1)
lon = ds.lon
lat = ds.lat
pcm = ax.pcolormesh(lon, lat, plotvar.T, transform=proj, cmap='cmo.tempo',
                    vmin=0, vmax=24, zorder=-1)

# # Plot Sea Ice
plotvar = ds_masks.mask_mon
cl = ax.contour(plotvar.lon, plotvar.lat,
                plotvar, colors="yellow",
                linewidths=2, transform=proj, levels=[0, 1], zorder=-1)
ax.clabel(cl, fontsize=12)

plotvar = ds_adt.isel(time=mon_sel).mean('time')
cl = ax.contour(plotvar.lon, plotvar.lat, plotvar.adt*100, colors="k",
                linewidths=0.75, transform=proj, levels=cints_adt)
ax.clabel(cl)

ax.set_title("$T^2$ Monthly, JFM Mean", fontsize=16)
cb = viz.hcbar(pcm, ax=ax, pad=0.01, fraction=0.045)
cb.set_label(
    "$T^2 (Months)$", fontsize=14)
cb.ax.tick_params(labelsize=14)

# Plot the Box
viz.plot_box(bbox_yeager,proj=proj,color="purple",linewidth=3,linestyle='dashed')

outname = figpath + "OISST_T2_lagmax%02i_wintertime.png" % (acfs.lags.data[-1])
plt.savefig(outname, dpi=150, bbox_inches='tight', transparent=True)



#%% Plot ACFs in the Yeager Box

xtks = np.arange(0,66,6)

lonf = -40
latf = 56
locfn,loctitle = proc.make_locstring(lonf,latf)

imon     = 1
acfs_reg = proc.sel_region_xr(acfs,bbox_yeager)
nlonr,nlatr,_,nlags = acfs_reg.shape

fig,ax = plt.subplots(1,1,figsize=(12.5,4),constrained_layout=True)
ax,_ = viz.init_acplot(imon,xtks,lags,ax=ax,title=None)

for o in range(nlonr):
    for a in range(nlatr):
        
        acf_in = acfs_reg.isel(lon=o,lat=a,mons=imon)
        ax.plot(lags,acf_in,alpha=0.15)


acf_in = acfs_reg.isel(mons=imon).mean(['lat','lon'])
ax.plot(lags,acf_in,color="k",label="Region Average ACF")

acf_in = acfs_reg.isel(mons=imon).sel(lon=lonf,lat=latf,method='nearest')

ax.plot(lags,acf_in,color="yellow",label="ACF (%s)" % loctitle)
        
ax.legend()
    

    
#%% Compute the REI (copied from calc_remidx_general)

def calc_rei(x): return proc.calc_remidx_xr(x, return_rei=True)

#Runs in approx 160 sec

st = time.time()
# Apply looping through basemonth, lon, lat. ('lon', 'lat', 'mon', 'rem_year')
rei_mon = xr.apply_ufunc(
    calc_rei,
    acfs,
    input_core_dims=[['lags']],
    output_core_dims=[['rem_year',]],
    vectorize=True,
)

print("Function applied in in %.2fs" % (time.time()-st))

# Add numbering based on the re-emergence year
rei_mon['rem_year'] = np.arange(1, 1+len(rei_mon.rem_year))


#%% Plot the Re-emergence Index
imon = 0

for imon in range(12):
    fig, ax, _ = viz.init_orthomap(1, 1, bbplot, figsize=(14, 6))
    ax = viz.add_coast_grid(ax, bbplot, fill_color='lightgray',
                            proj=proj, line_color="dimgray")
    
    plotvar = rei_mon.isel(mons=imon,rem_year=0)
    lon     = ds.lon
    lat     = ds.lat
    pcm = ax.pcolormesh(lon, lat, plotvar.T, transform=proj, cmap='cmo.dense',
                        vmin=0, vmax=0.75, zorder=-1)
    
    # # Plot Sea Ice
    plotvar = ds_masks.mask_mon
    cl = ax.contour(plotvar.lon, plotvar.lat,
                    plotvar, colors="yellow",
                    linewidths=2, transform=proj, levels=[0, 1], zorder=-1)
    ax.clabel(cl, fontsize=12)
    
    plotvar = ds_adt.isel(time=imon)
    cl = ax.contour(plotvar.lon, plotvar.lat, plotvar.adt*100, colors="k",
                    linewidths=0.75, transform=proj, levels=cints_adt)
    ax.clabel(cl)
    
    ax.set_title("$Re-emergence Index, %s" % (mons3[imon]), fontsize=16)
    cb = viz.hcbar(pcm, ax=ax, pad=0.01, fraction=0.045)
    cb.set_label(
        "$REI$", fontsize=14)
    cb.ax.tick_params(labelsize=14)
    
    # Plot the Box
    viz.plot_box(bbox_yeager,proj=proj,color="purple",linewidth=3,linestyle='dashed')
    
    outname = figpath + "OISST_REI_lagmax%02i_mon%02i.png" % (acfs.lags.data[-1], imon+1)
    plt.savefig(outname, dpi=150, bbox_inches='tight', transparent=True)

#%% Take Wintertime Avg

mons_sel = [0,1,2]

fig, ax, _ = viz.init_orthomap(1, 1, bbplot, figsize=(14, 6))
ax = viz.add_coast_grid(ax, bbplot, fill_color='lightgray',
                        proj=proj, line_color="dimgray")

plotvar = rei_mon.isel(mons=mons_sel,rem_year=0).mean('mons')
lon     = ds.lon
lat     = ds.lat
pcm = ax.pcolormesh(lon, lat, plotvar.T, transform=proj, cmap='cmo.dense',
                    vmin=0, vmax=0.75, zorder=-1)

# # Plot Sea Ice
plotvar = ds_masks.mask_mon
cl = ax.contour(plotvar.lon, plotvar.lat,
                plotvar, colors="yellow",
                linewidths=2, transform=proj, levels=[0, 1], zorder=-1)
ax.clabel(cl, fontsize=12)

plotvar = ds_adt.isel(time=mons_sel).mean('time')
cl = ax.contour(plotvar.lon, plotvar.lat, plotvar.adt*100, colors="k",
                linewidths=0.75, transform=proj, levels=cints_adt)
ax.clabel(cl)

ax.set_title("$Re-emergence Index, JFM Mean", fontsize=16)
cb = viz.hcbar(pcm, ax=ax, pad=0.01, fraction=0.045)
cb.set_label(
    "$REI$", fontsize=14)
cb.ax.tick_params(labelsize=14)

# Plot the Box
viz.plot_box(bbox_yeager,proj=proj,color="purple",linewidth=3,linestyle='dashed')

outname = figpath + "OISST_REI_lagmax%02i_Wintertime.png" % (acfs.lags.data[-1])
plt.savefig(outname, dpi=150, bbox_inches='tight', transparent=True)

#%% Look at relatioship...

rei_wint = rei_mon.isel(mons=mons_sel,rem_year=0).mean('mons').data.flatten()
t2_wint  = t2[:,:,mon_sel].mean(-1).flatten()


xx,yy = np.meshgrid(rei_mon.lon.data,rei_mon.lat.data)


fig,ax   = plt.subplots(1,1,constrained_layout=True,figsize=(12,12))

sc = ax.scatter(t2_wint,rei_wint,c=yy.flatten(),s=25,alpha=0.15)

sc = ax.scatter(t2_wint,rei_wint,c=yy.flatten(),s=0,alpha=1)

ax.set_xlabel("T2 (Months)",fontsize=16)
ax.set_ylabel("REI Index (Correlation)",fontsize=16)

viz.hcbar(sc,ax=ax)

sc = ax.scatter(t2_wint,rei_wint,c=yy.flatten(),s=25,alpha=0.15)


