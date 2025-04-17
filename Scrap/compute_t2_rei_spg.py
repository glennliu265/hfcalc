#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Copied section from compare_spg_box_sst 

Created on Wed Mar  5 11:12:55 2025

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

import matplotlib as mpl
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

def winter_avg(ds):
    dswint = ds.where(ds.time.dt.month.isin([12,1,2,3]),drop=True)
    dswint = dswint.groupby('time.year').mean('time')
    return dswint


   #ts_wint = ts.groupby('time.season')
    
    


#%% Copy the files in


# Plot Settings
mpl.rcParams['font.family'] = 'Avenir'
proj = ccrs.PlateCarree()
bbplot = [-80, 0, 35, 75]
mons3 = proc.get_monstr()


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


#%% Anomalize Timeseries

# Deseason
ds_anom = [proc.xrdeseason(ds.sst) for ds in ds_all]

# Detrend at each location % (Note, takes just under 1 minute...)
order   = 3
def pointwise_detrend(ds,order):
    
    def point_detrend(ts,order):
        ntime  = len(ts)
        times  = np.arange(ntime)
        if np.any(np.isnan(ts)):
            dt_out = times*np.nan
        else:
            dt_out = proc.polyfit_1d(times,ts,order)
        return dt_out[2]
    detrend_pt = lambda x: point_detrend(x,order)
    
    st = time.time()
    ds_detrended = xr.apply_ufunc(
        detrend_pt,
        ds,
        input_core_dims=[['time']],
        output_core_dims=[['time']],
        vectorize=True,
        )
    print("Detrended in %.2fs" % (time.time()-st))
    return ds_detrended
ds_anom_detrended = [pointwise_detrend(ds,order) for ds in ds_anom]

#%% Next, let's compute the wintertime average, the, Lag 1 Autocorrelation


import tqdm

ds_winter = [winter_avg(ds) for ds in ds_anom_detrended]

def pointwise_lag_manual(ds,lag):
    ds = ds.transpose('year','lat','lon')
    
    def calc_leadlag(ts,lag):
        if np.any(np.isnan(ts)):
            return np.nan
        tslag  = ts[lag:]
        tsbase = ts[:(-lag)]
        try:
            x     = tsbase
            y     = tslag
            xanom = x - x.mean()
            yanom = y - y.mean()
            xstd  = x.std()
            ystd  = y.std()
            
            cov   = ((xanom * yanom).sum())/len(xanom)
            corr  = cov/(xstd*ystd)

            #corr   = np.corrcoef(tsbase,tslag)[0,1]
            
        except:
            return np.nan
        return corr
        
    ntime,nlat,nlon = ds.shape
    corrs           = np.zeros((nlat,nlon))*np.nan
    
    for o in tqdm.tqdm(range(nlon)):
        for a in range(nlat):
            ts = ds.isel(lon=o,lat=a).data
            corrs[a,o] = calc_leadlag(ts,lag)#.copy()
    return corrs
            
# lag = 1
# ds  = ds_winter[0]
# corrs = pointwise_lag_manual(ds,lag)

lag = 1
dswint_r1 = [pointwise_lag_manual(ds,lag) for ds in ds_winter] 

#%%



            
#%% Visualize the Wintertime R1 Patterns
bbox_yeager = [-50,-10,50,60] # Original

nn         = 0
for nn in range(ndatasets):
    plotvar    = dswint_r1[nn]
    plotlon    = ds_winter[nn].lon
    plotlat    = ds_winter[nn].lat
    
    proj       = ccrs.PlateCarree()
    cints      = np.arange(0,1.05,0.05)
    cints_mld  = np.arange(0,1100,100)
    
    bbplot     = [-80, 0, 35, 65]
    fig, ax, _ = viz.init_orthomap(1, 1, bbplot, figsize=(14, 6))
    ax         = viz.add_coast_grid(ax, bbplot, fill_color='lightgray',
                            proj=proj, line_color="dimgray")
    
    pcm        = ax.pcolormesh(plotlon,plotlat,plotvar,
                               vmin=0,vmax=1,transform=proj,cmap='cmo.amp')
    
    cb = viz.hcbar(pcm,)
    cb.set_label("R1 (Wintertime)")
    
    
    # # Plot Sea Ice
    plotvar = ds_masks.mask_mon
    cl = ax.contour(plotvar.lon, plotvar.lat,
                    plotvar, colors="yellow",
                    linewidths=2, transform=proj, levels=[0, 1], zorder=2)
    ax.clabel(cl, fontsize=12)
    
    # Plot the SSH
    plotvar = ds_adt.isel(time=[11,0,1,2]).mean('time')
    cl = ax.contour(plotvar.lon, plotvar.lat, plotvar.adt*100, colors="k",
                    linewidths=0.75, transform=proj, levels=cints_adt)
    ax.clabel(cl)
    
    # Plot MLD
    plotvar= ds_mld.mld.max('month')#isel(month=[11,0,1,2]).mean('month')
    cl = ax.contour(plotvar.lon, plotvar.lat, plotvar, colors="darkblue",
                    linewidths=0.55, transform=proj,levels=cints_mld,linestyles='dotted')
    ax.clabel(cl)
    
    # Plot the Box
    viz.plot_box(bbox_yeager,proj=proj,color="purple",linewidth=3,linestyle='dashed')

    
    ax.set_title("%s (%s)" % (dset_names[nn],ystrs[nn]))
    outname = "%s%s_Wintertime_R1.png" % (figpath,dset_names[nn])
    plt.savefig(outname,dpi=150,bbox_inches='tight',transparent=True)

# =================================
#%% Load + compute OISST REI and compare with MLD (copied from check_acf)


dpath2  = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/NATL_proc_obs/"
nc      =  dpath2 + "OISST_1982_2020_ACF_lag00to60_ALL_ensALL.nc"
ds      = xr.open_dataset(nc).load()
acfs = ds.isel(thres=0,ens=0).acf.squeeze()

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



t2 = proc.calc_T2(acfs,axis=-1,ds=True)


#%% Plot the Re-emergence Index with Month
imon    = 2
plotmax = False

cints_mld  = np.arange(0,1100,50)
plot_range = True

#for imon in range(12):
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

# # plot ADT
# plotvar = ds_adt.isel(time=imon)
# cl = ax.contour(plotvar.lon, plotvar.lat, plotvar.adt*100, colors="k",
#                 linewidths=0.75, transform=proj, levels=cints_adt)
# ax.clabel(cl)


# Plot MLD
if plot_range:
    plotvar = ds_mld.mld.max('month') - ds_mld.mld.min('month')
    mldstr  = "Range"
else:
    plotvar = ds_mld.mld.isel(month=imon)
    mldstr  = ""
    
cl      = ax.contour(plotvar.lon, plotvar.lat, plotvar, colors="k",
                linewidths=0.75, transform=proj,levels=cints_mld,linestyles='solid')
ax.clabel(cl)

ax.set_title("Re-emergence Index, %s" % (mons3[imon]), fontsize=16)
cb = viz.hcbar(pcm, ax=ax, pad=0.01, fraction=0.045)
cb.set_label(
    "$REI$", fontsize=14)
cb.ax.tick_params(labelsize=14)

# Plot the Box
#viz.plot_box(bbox_yeager,proj=proj,color="purple",linewidth=3,linestyle='dashed')

outname = figpath + "OISST_REI_lagmax%02i_MLD%s_mon%02i.png" % (acfs.lags.data[-1], mldstr,imon+1)
plt.savefig(outname, dpi=150, bbox_inches='tight', transparent=True)

#%% Plot Relationships (MLD RANGE vs T2)


mld_offset = ds_mld.mld.data[:,1:,1:].transpose(2,1,0)
mld_range  = mld_offset.max(-1) - mld_offset.min(-1)
mld_max    = mld_offset.max(-1)
t2_wint    = t2[:,:,[0,1,2]].mean(-1)
rei_wint   = rei_mon.isel(mons=[0,1,2],rem_year=1).mean('mons')

xx,yy      = np.meshgrid(rei_mon.lon.data,rei_mon.lat.data)
fig,axs     = plt.subplots(1,2,constrained_layout=True,figsize=(12,6))

#sc         = ax.scatter(t2_wint,mld_range,s=25,alpha=0.15)
ax = axs[0]
sc         = ax.scatter(mld_range,t2_wint,c=xx.flatten(),s=25,alpha=0.15)
ax.set_xlabel("MLD Range (meters)",fontsize=16)
ax.set_ylabel("T2 (Months)",fontsize=16)


ax = axs[1]
sc         = ax.scatter(mld_range,rei_wint,c=xx.flatten(),s=25,alpha=0.15)
ax.set_xlabel("MLD Range (meters)",fontsize=16)
ax.set_ylabel("REI",fontsize=16)

#sc = ax.scatter(t2_wint,rei_wint,c=yy.flatten(),s=25,alpha=0.15)

#sc = ax.scatter(t2_wint,rei_wint,c=yy.flatten(),s=0,alpha=1)

viz.hcbar(sc,ax=axs.flatten())

#sc = ax.scatter(t2_wint,rei_wint,c=yy.flatten(),s=25,alpha=0.15)



#%%

    

# def pointwise_lag(ds,lag):
    
    
    
    
#     calclag = lambda x: calc_leadlag(x,lag)
    
#     st = time.time()
#     ds_lag = xr.apply_ufunc(
#         calclag,
#         ds,
#         input_core_dims=[['year']],
#         output_core_dims=[[],],
#         vectorize=True,
#         )
#     print("Detrended in %.2fs" % (time.time()-st))
#     return ds_lag



# def point_detrend(ts,order):
#     ntime  = len(ts)
#     times  = np.arange(ntime)
#     if np.any(np.isnan(ts)):
#         dt_out = times*np.nan
#     else:
#         dt_out = proc.polyfit_1d(times,ts,order)
#     return dt_out[2]
# detrend_pt = lambda x: point_detrend(x,order)
    
# ds = ds_anom[0]


#%% Preprocessing

bbox_yeager = [-50,-10,50,60] # Original
#bbox_yeager = [-30,-10,50,60] # Away from Sea Ice 

bbfn,bbstr  = proc.make_locstring_bbox(bbox_yeager)

spg_anom    = [preproc_ds(ds.sst,bbox_yeager) for ds in ds_all]
spg_acf     = [calc_acf_mon(sst.data,lags) for sst in spg_anom]
