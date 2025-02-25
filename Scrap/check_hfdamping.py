#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Check Heat Flux Damping from Observations (ERA5, OISST) from calc_hff_general

Created on Wed Feb 19 10:56:59 2025
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

#%%
mpl.rcParams['font.family']     = 'Avenir'
figpath                         = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/03_reemergence/02_Figures/20250221/"
proc.makedir(figpath)
#%% HFF Estimates

dpath         = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/" 
ncname_obs    = "ERA5_RegridOISST_hfdamping_NAtl_1982to2020_ensorem1_detrend1.nc"
ds            = xr.open_dataset(dpath + ncname_obs)

dt            = 3600*24#*30 #Effective Processing Period of 1 day


lon           = ds.lon
lat           = ds.lat
thflx_damping = ds.thflx_damping / dt * -1

cov           = ds.cov
autocov       = ds.autocov

cc            = ds.sst_flx_crosscorr
ac            = ds.sst_autocorr



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
mask_exceed_mon   = (exceed_mon/max_count_daily) < 0.05

ds_masks = xr.merge([mask_exceed_daily.rename("mask_daily"),
                     mask_exceed_mon.rename("mask_mon")
                     ])
edict             = proc.make_encoding_dict(ds_masks)
outname           = dpath_proc + "../OISST_ice_masks_1981_2020.nc"
ds_masks.to_netcdf(outname,encoding=edict)

#%% Plot Sea Ice (Daily Max)
proj      = ccrs.PlateCarree()
bbplot    = [-80,0,35,75]
fig,ax,_  = viz.init_orthomap(1,1,bbplot,figsize=(14,6))
ax        = viz.add_coast_grid(ax,bbplot,fill_color='lightgray',proj=proj)

# Daily Max
plotvar   = ds_sice_natl.icec.max('time')

# Monthly Mean
#plotvar = ds_sice_monmean_natl.icec.mean('month').max('year')

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
        title   = "Daily Resolution"
        plotvar = exceed_daily / max_count_daily
    else:
        title="Monthly Resolution"
        plotvar = exceed_mon / max_count_monthly
        
    
    ax.set_title(title,fontsize=24,y=1.05)
    
    pcm       = ax.pcolormesh(ds_sice_natl.lon,ds_sice_natl.lat,
                              plotvar,transform=proj,cmap='cmo.ice',
                             vmin=0,vmax=1,zorder=-1)



    
    
    cl        = ax.contour(ds_sice_natl.lon,ds_sice_natl.lat,
                              plotvar>0.05,colors="forestgreen",
                          linewidths=1.5,transform=proj,levels=[0,1])
    ax.clabel(cl,fontsize=14)

#ax.set_title("Month %s, Lag %i" % (mons3[imon],ilag+1),fontsize=14)
cb = viz.hcbar(pcm,ax=axs.flatten(),pad=0.01,fraction=0.045)
cb.set_label("% of time exceeding 5% Ice Concentration",fontsize=24)
cb.ax.tick_params(labelsize=16)   
    
outname= figpath + "OISST_Max_siceExceed.png"
plt.savefig(outname,dpi=150,bbox_inches='tight')




#%% Plot HFF Estimates with sea ice

ilag = 0
imon = 1
mons3 = proc.get_monstr()


proj     = ccrs.PlateCarree()
bbplot   = [-80,0,35,75]

#for imon in range(12):
fig,ax,_ = viz.init_orthomap(1,1,bbplot,figsize=(14,6))
ax       = viz.add_coast_grid(ax,bbplot,fill_color='lightgray',proj=proj,line_color="dimgray")

plotvar  = thflx_damping.isel(lag=ilag,month=imon) # 
pcm      = ax.pcolormesh(lon,lat,plotvar,transform=proj,cmap='cmo.balance',
                         vmin=-55,vmax=55,zorder=-1)

# # Plot Sea Ice
# plotvar   = ds_sice_natl.icec.max('time')
# cl        = ax.contour(ds_sice_natl.lon,ds_sice_natl.lat,
#                           plotvar>0.05,colors="cyan",
#                       linewidths=1,transform=proj,levels=[0,],zorder=0)
# ax.clabel(cl,fontsize=12)

plotvar   = ds_adt.isel(time=imon)
cl        = ax.contour(plotvar.lon,plotvar.lat,plotvar.adt*100,colors="k",
                      linewidths=0.75,transform=proj,levels=cints_adt)
ax.clabel(cl)

ax.set_title("Month %s, Lag %i" % (mons3[imon],ilag+1),fontsize=14)
cb = viz.hcbar(pcm,ax=ax,pad=0.01,fraction=0.045)
cb.set_label("Turbulent Heat Flux Feedback (W m$^{-2}$ $\degree$C$^{-1}$)",fontsize=14)
cb.ax.tick_params(labelsize=14)

outname= figpath + "OISST_ERA5_THFF_lag%i_mon%02i.png" % (ilag+1,imon+1)
plt.savefig(outname,dpi=150,bbox_inches='tight')


#%% Check the Data



nc_sst          = dpath_proc + "OISST_sst_NAtl_1982to2020.nc"
nc_flx          = dpath_proc + "ERA5_RegridOISST_thflx_NAtl_1982to2020.nc"
nc_flx_raw      = dpath_proc + "ERA5_thflx_NAtl_1982to2020.nc"

ds_sst          = xr.open_dataset(nc_sst)
ds_flx_regrid   = xr.open_dataset(nc_flx)
ds_flx_raw      = xr.open_dataset(nc_flx_raw)

nc_adt          = dpath_proc + "AVISO_adt_NAtl_1993_2022_clim.nc"
ds_adt          = xr.open_dataset(nc_adt)


#%% Plot Autocorr

ilag = 0
imon = 1
mons3 = proc.get_monstr()

cints_adt = np.arange(-100,110,10)

proj      = ccrs.PlateCarree()
bbplot    = [-80,0,35,75]
fig,ax,_  = viz.init_orthomap(1,1,bbplot,figsize=(14,6))
ax        = viz.add_coast_grid(ax,bbplot,fill_color='k',proj=proj)

plotvar   = ac.isel(lag=ilag,month=imon) # 
pcm       = ax.pcolormesh(lon,lat,plotvar,transform=proj,cmap=cm.acton,
                         vmin=0.30,vmax=1)

plotvar   = ds_adt.isel(time=imon)
cl        = ax.contour(plotvar.lon,plotvar.lat,plotvar.adt*100,colors='maroon',
                      linewidths=0.75,transform=proj,levels=cints_adt)
ax.clabel(cl)

ax.set_title("Month %s, Lag %i" % (mons3[imon],ilag+1),fontsize=14)
cb = fig.colorbar(pcm,ax=ax,pad=0.01)
cb.set_label("Lag %i Correlation" % (ilag+1),fontsize=14)
cb.ax.tick_params(labelsize=14)

outname= figpath + "OISST_ERA5_SST_Autocorr_lag%i_mon%02i.png" % (ilag+1,imon+1)
plt.savefig(outname,dpi=150,bbox_inches='tight',transparent=True)


#%% Examine Numerator and Denominator


imonth = 9
ilag   = 0

fig,axs,_  = viz.init_orthomap(3,1,bbplot,figsize=(14,22))


for ii in range(3):
    ax = axs[ii]
    ax = viz.add_coast_grid(ax,bbplot,fill_color='k',proj=proj)
    
    if ii == 0:
        plotvar = cov.isel(month=imonth,lag=ilag)/dt
        pname   = "Covariance (SST,FLX)"
        
        vlims   = [-45,45]
    elif ii == 1:
        plotvar = autocov.isel(month=imonth,lag=ilag)
        pname   = "Autocovariance (SST)"
        vlims   = [-1,1]
    else:
        plotvar = thflx_damping.isel(month=imonth,lag=ilag)
        pname   = "THFLX Feedback"
        vlims   = [-45,45]
        
        
    pcm      = ax.pcolormesh(lon,lat,plotvar,transform=proj,cmap='cmo.balance',
                             vmin = vlims[0],vmax = vlims[1])
    
    fig.colorbar(pcm,ax=ax,pad=0.01)
    
    ax.set_title(pname,fontsize=16)
plt.suptitle("Month %s, Lag %i" % (mons3[imonth],ilag+1),fontsize=42)




#%% Load CESM1 Estimates

#ds_cesm = 
cpath           = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/CESM-PIC-FULL-Damping/"
nhflx_cesm      = np.load(cpath+ "NHFLX_Damping_monwin3_lags123_ensorem1_lag1_pcs2_monwin3.npz.npy")
cc_cesm         = np.load(cpath+ "NHFLX_Crosscorrelation_monwin3_lags123_ensorem1_lag1_pcs2_monwin3.npz.npy")
ac_cesm         = np.load(cpath+ "SST_Autocorrelation_monwin3_lags123_ensorem1_lag1_pcs2_monwin3.npz.npy")

htrpath         = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/CESM-HTR-RCP85/"
ds_cesm         = xr.open_dataset(htrpath+"CESM1_htr_qnet_damping_ensorem1_detrend1_1920to2005_allens.nc")                      

ds_cesm         = proc.lon360to180_xr(ds_cesm)

qnet_damping    = ds_cesm.qnet_damping
cc_cesm         = ds_cesm.sst_flx_crosscorr
ac_cesm         = ds_cesm.sst_autocorr
cov_cesm        = ds_cesm.cov
autocov_cesm    = ds_cesm.autocov

clon = ds_cesm.lon
clat = ds_cesm.lat

#%% Make similar plot as above, but for cesm

imonth = 1
ilag   = 0
iens   = 0
fig,axs,_  = viz.init_orthomap(3,1,bbplot,figsize=(14,22))

for ii in range(3):
    ax = axs[ii]
    ax = viz.add_coast_grid(ax,bbplot,fill_color='k',proj=proj)
    
    if ii == 0:
        
        plotvar = cov_cesm.isel(month=imonth,lag=ilag,ens=iens)
        pname   = "Covariance (SST,FLX)"
        vlims   = [-10,10]
        
    elif ii == 1:
        
        plotvar = autocov_cesm.isel(month=imonth,lag=ilag,ens=iens)
        pname   = "Autocovariance (SST)"
        vlims   = [-1,1]
        
    else:
        
        plotvar = qnet_damping.isel(month=imonth,lag=ilag,ens=iens)
        pname   = "NHFLX Feedback"
        vlims   = [-45,45]
        
        
    pcm      = ax.pcolormesh(clon,clat,plotvar,transform=proj,cmap='cmo.balance',
                             vmin = vlims[0],vmax = vlims[1])
    
    fig.colorbar(pcm,ax=ax,pad=0.01)
    
    ax.set_title(pname)
plt.suptitle("Month %s, Lag %i, Ens %02i" % (mons3[imonth],ilag+1,iens+1),fontsize=42)





#%%

(ds_flx_regrid.thflx.isel(time=0)/dt).plot()

#%% Take Average Over the Region








