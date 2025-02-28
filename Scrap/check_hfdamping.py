#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Check Heat Flux Damping from Observations (ERA5, OISST) from calc_hff_general

Created on Wed Feb 19 10:56:59 2025
@author: gliu

"""


from amv import proc, viz
import amv.proc as hf  # Update hf with actual hfutils script, most relevant functions
import scm
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


# %% Import modules
stormtrack = 0
if stormtrack:
    sys.path.append("/home/glliu/00_Scripts/01_Projects/00_Commons/")
    sys.path.append(
        "/home/glliu/00_Scripts/01_Projects/01_AMV/02_stochmod/stochmod/model/")

    # Path to the processed dataset (qnet and ts fields, full, time x lat x lon)
    # datpath =  "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_RCP85/01_PREPROC/"
    datpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/anom/"
    figpath = "/home/glliu/02_Figures/01_WeeklyMeetings/20240621/"

    # hfcalc_path = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/03_Scripts/hfcalc/"
else:
    sys.path.append(
        "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/03_Scripts/stochmod/model/")
    sys.path.append(
        "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/03_Scripts/")

    # Path to the processed dataset (qnet and ts fields, full, time x lat x lon)
    datpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/"
    figpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/02_Figures/20220511/"

    # hfcalc_path = "/stormtrack/home/glliu/00_Scripts/01_Projects/01_AMV/01_hfdamping/hfcalc/" # hfcalc module


# Import Hf Calc params
# sys.path.append(hfcalc_path)
# import hfcalc_params as hp


# %% User Edits

# Plot Settings
mpl.rcParams['font.family'] = 'Avenir'
proj = ccrs.PlateCarree()
bbplot = [-80, 0, 35, 75]
mons3 = proc.get_monstr()

# Paths
figpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/03_reemergence/02_Figures/20250221/"
proc.makedir(figpath)

# Load Sea Ice Masks
dpath_ice = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/NATL_proc_obs/"
nc_masks = dpath_ice + "OISST_ice_masks_1981_2020.nc"
ds_masks = xr.open_dataset(nc_masks).load()

# Load AVISO
dpath_aviso = dpath_ice + "proc/"
nc_adt = dpath_aviso + "AVISO_adt_NAtl_1993_2022_clim.nc"
ds_adt = xr.open_dataset(nc_adt).load()


# Other Plotting Selections
cints_adt = np.arange(-100, 110, 10)


# This should be in the hfcalc_params
# +1 due to inclusive, -2 and -1 due to year window, *3 due to month window
dof = (2020-1982 + 1 - 2 - 1) * 3

# %% HFF Estimates

dpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/"
ncname_obs = "ERA5_RegridOISST_hfdamping_NAtl_1982to2020_ensorem1_detrend1.nc"
ds = xr.open_dataset(dpath + ncname_obs)

dt = 3600*24  # *30 #Effective Processing Period of 1 day


lon = ds.lon
lat = ds.lat
thflx_damping = ds.thflx_damping / dt * -1

cov = ds.cov
autocov = ds.autocov

cc = ds.sst_flx_crosscorr
ac = ds.sst_autocorr

# %% Plot HFF Estimates with sea ice

# Select Lag and Month
ilag = 0
imon = 1

proj = ccrs.PlateCarree()
bbplot = [-80, 0, 35, 75]

# for imon in range(12):

fig, ax, _ = viz.init_orthomap(1, 1, bbplot, figsize=(14, 6))
ax = viz.add_coast_grid(ax, bbplot, fill_color='lightgray',
                        proj=proj, line_color="dimgray")

plotvar = thflx_damping.isel(lag=ilag, month=imon)
pcm = ax.pcolormesh(lon, lat, plotvar, transform=proj, cmap='cmo.balance',
                    vmin=-55, vmax=55, zorder=-1)

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

ax.set_title("Month %s, Lag %i" % (mons3[imon], ilag+1), fontsize=16)
cb = viz.hcbar(pcm, ax=ax, pad=0.01, fraction=0.045)
cb.set_label(
    "Turbulent Heat Flux Feedback (W m$^{-2}$ $\degree$C$^{-1}$)", fontsize=14)
cb.ax.tick_params(labelsize=14)

outname = figpath + "OISST_ERA5_THFF_lag%i_mon%02i.png" % (ilag+1, imon+1)
plt.savefig(outname, dpi=150, bbox_inches='tight', transparent=True)


# %% Check the Data


nc_sst = dpath_proc + "OISST_sst_NAtl_1982to2020.nc"
nc_flx = dpath_proc + "ERA5_RegridOISST_thflx_NAtl_1982to2020.nc"
nc_flx_raw = dpath_proc + "ERA5_thflx_NAtl_1982to2020.nc"

ds_sst = xr.open_dataset(nc_sst)
ds_flx_regrid = xr.open_dataset(nc_flx)
ds_flx_raw = xr.open_dataset(nc_flx_raw)


# %% Plot Autocorr

ilag = 0
imon = 1
mons3 = proc.get_monstr()


proj = ccrs.PlateCarree()
bbplot = [-80, 0, 35, 75]
fig, ax, _ = viz.init_orthomap(1, 1, bbplot, figsize=(14, 6))
ax = viz.add_coast_grid(ax, bbplot, fill_color='k', proj=proj)

plotvar = ac.isel(lag=ilag, month=imon)
pcm = ax.pcolormesh(lon, lat, plotvar, transform=proj, cmap=cm.acton,
                    vmin=0.30, vmax=1)

plotvar = ds_adt.isel(time=imon)
cl = ax.contour(plotvar.lon, plotvar.lat, plotvar.adt*100, colors='maroon',
                linewidths=0.75, transform=proj, levels=cints_adt)
ax.clabel(cl)

ax.set_title("Month %s, Lag %i" % (mons3[imon], ilag+1), fontsize=14)
cb = fig.colorbar(pcm, ax=ax, pad=0.01)
cb.set_label("Lag %i Correlation" % (ilag+1), fontsize=14)
cb.ax.tick_params(labelsize=14)

outname = figpath + \
    "OISST_ERA5_SST_Autocorr_lag%i_mon%02i.png" % (ilag+1, imon+1)
plt.savefig(outname, dpi=150, bbox_inches='tight', transparent=True)


# %% Examine Numerator and Denominator


imonth = 9
ilag = 0

fig, axs, _ = viz.init_orthomap(3, 1, bbplot, figsize=(14, 22))


for ii in range(3):
    ax = axs[ii]
    ax = viz.add_coast_grid(ax, bbplot, fill_color='k', proj=proj)

    if ii == 0:
        plotvar = cov.isel(month=imonth, lag=ilag)/dt
        pname = "Covariance (SST,FLX)"

        vlims = [-45, 45]
    elif ii == 1:
        plotvar = autocov.isel(month=imonth, lag=ilag)
        pname = "Autocovariance (SST)"
        vlims = [-1, 1]
    else:
        plotvar = thflx_damping.isel(month=imonth, lag=ilag)
        pname = "THFLX Feedback"
        vlims = [-45, 45]

    pcm = ax.pcolormesh(lon, lat, plotvar, transform=proj, cmap='cmo.balance',
                        vmin=vlims[0], vmax=vlims[1])

    fig.colorbar(pcm, ax=ax, pad=0.01)

    ax.set_title(pname, fontsize=16)
plt.suptitle("Month %s, Lag %i" % (mons3[imonth], ilag+1), fontsize=42)

# %% Load CESM1 Estimates

# ds_cesm =
cpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/CESM-PIC-FULL-Damping/"
nhflx_cesm = np.load(
    cpath + "NHFLX_Damping_monwin3_lags123_ensorem1_lag1_pcs2_monwin3.npz.npy")
cc_cesm = np.load(
    cpath + "NHFLX_Crosscorrelation_monwin3_lags123_ensorem1_lag1_pcs2_monwin3.npz.npy")
ac_cesm = np.load(
    cpath + "SST_Autocorrelation_monwin3_lags123_ensorem1_lag1_pcs2_monwin3.npz.npy")

htrpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/CESM-HTR-RCP85/"
ds_cesm = xr.open_dataset(
    htrpath+"CESM1_htr_qnet_damping_ensorem1_detrend1_1920to2005_allens.nc")

ds_cesm = proc.lon360to180_xr(ds_cesm)

qnet_damping = ds_cesm.qnet_damping
cc_cesm = ds_cesm.sst_flx_crosscorr
ac_cesm = ds_cesm.sst_autocorr
cov_cesm = ds_cesm.cov
autocov_cesm = ds_cesm.autocov

clon = ds_cesm.lon
clat = ds_cesm.lat

# %% Make similar plot as above, but for cesm

imonth = 1
ilag = 0
iens = 0
fig, axs, _ = viz.init_orthomap(3, 1, bbplot, figsize=(14, 22))

for ii in range(3):
    ax = axs[ii]
    ax = viz.add_coast_grid(ax, bbplot, fill_color='k', proj=proj)

    if ii == 0:

        plotvar = cov_cesm.isel(month=imonth, lag=ilag, ens=iens)
        pname = "Covariance (SST,FLX)"
        vlims = [-10, 10]

    elif ii == 1:

        plotvar = autocov_cesm.isel(month=imonth, lag=ilag, ens=iens)
        pname = "Autocovariance (SST)"
        vlims = [-1, 1]

    else:

        plotvar = qnet_damping.isel(month=imonth, lag=ilag, ens=iens)
        pname = "NHFLX Feedback"
        vlims = [-45, 45]

    pcm = ax.pcolormesh(clon, clat, plotvar, transform=proj, cmap='cmo.balance',
                        vmin=vlims[0], vmax=vlims[1])

    fig.colorbar(pcm, ax=ax, pad=0.01)

    ax.set_title(pname)
plt.suptitle("Month %s, Lag %i, Ens %02i" %
             (mons3[imonth], ilag+1, iens+1), fontsize=42)


# %%

(ds_flx_regrid.thflx.isel(time=0)/dt).plot()

# %% Take Average Over the Region

# %% Do significance testing (based on prep_HF_input)

hff = thflx_damping
rsst = ac
rflx = cc
setdict = {  # Taken from hfcalc_params
    'ensorem': 1,      # 1=enso removed, 0=not removed
    'ensolag': 1,      # Lag Applied toENSO and Variable before removal
    'monwin': 3,      # Size of month window for HFF calculations
    'detrend': 1,      # Whether or not variable was detrended
    'tails': 2,      # tails for t-test

    'p': 0.05,   # p-value for significance testing
    'sellags': [0,],   # Lags included (indices, so 0=lag1)
    'lagstr': "lag1",  # Name of lag based on sellags
    # Significance test option: 1 (No Mask); 2 (SST autocorr); 3 (SST-FLX crosscorr); 4 (Both), 5 (Replace with SLAB values)
    'method': 4
}
# dof was set above

# Compute and Apply Mask
st = time.time()
dampingmasked, freq_success, sigmask = scm.prep_HF(hff, rsst, rflx,
                                                   setdict['p'], setdict['tails'], dof, setdict['method'],
                                                   returnall=True)  # expects, [month x lag x lat x lon], should generalized with ensemble dimension?
print("Completed significance testing in %.2fs" % (time.time()-st))

# Remake Mask such that 1 = insignificant
landmask = xr.where(np.isnan(thflx_damping.isel(month=0, lag=0)), np.nan, 1)

sigmask_inv = np.isnan(sigmask.copy())
sigmask_inv = sigmask_inv * landmask.data[None, None, :, :]


# Make Masks for AC and CC Only
_,_,acmask = scm.prep_HF(hff, rsst, rflx,
                        setdict['p'], setdict['tails'], dof, 2,
                        returnall=True)
_,_,ccmask = scm.prep_HF(hff, rsst, rflx,
                        setdict['p'], setdict['tails'], dof, 3,
                        returnall=True)
def flip_mask(mask,landmask):
    mask_inv = np.isnan(mask.copy())
    mask_inv = mask_inv * landmask.data[None, None, :, :]
    return mask_inv

acmask_inv = flip_mask(acmask,landmask)
ccmask_inv = flip_mask(ccmask,landmask)
    

# %% Plot the significance testing results (Copied HFF plot froma bove with sea ice)


# Select Lag and Month
ilag = 2
imon = 7

proj = ccrs.PlateCarree()
bbplot = [-80, 0, 35, 75]

# for imon in range(12):

fig, ax, _ = viz.init_orthomap(1, 1, bbplot, figsize=(14, 6))
ax = viz.add_coast_grid(ax, bbplot, fill_color='lightgray',
                        proj=proj, line_color="dimgray")

plotvar = thflx_damping.isel(lag=ilag, month=imon)
pcm = ax.pcolormesh(lon, lat, plotvar, transform=proj, cmap='cmo.balance',
                    vmin=-55, vmax=55, zorder=-1)

# # Plot Sea Ice
plotvar = ds_masks.mask_mon
cl = ax.contour(plotvar.lon, plotvar.lat,
                plotvar, colors="yellow",
                linewidths=2, transform=proj, levels=[0, 1], zorder=-1)
ax.clabel(cl, fontsize=12)

# Plot the SSH
plotvar = ds_adt.isel(time=imon)
cl = ax.contour(plotvar.lon, plotvar.lat, plotvar.adt*100, colors="k",
                linewidths=0.75, transform=proj, levels=cints_adt)
ax.clabel(cl)

# Plot the dots
# Plot significant points
plotmsk1 = sigmask_inv[imon, ilag, :, :]
# Had to reverse and transpose to lon x lat
viz.plot_mask(plotvar.lon, plotvar.lat, plotmsk1.T, reverse=True, ax=ax, proj=proj, geoaxes=True,
              markersize=.4, color='gray')

ax.set_title("Month %s, Lag %i" % (mons3[imon], ilag+1), fontsize=16)
cb = viz.hcbar(pcm, ax=ax, pad=0.01, fraction=0.045)
cb.set_label(
    "Turbulent Heat Flux Feedback (W m$^{-2}$ $\degree$C$^{-1}$)", fontsize=14)
cb.ax.tick_params(labelsize=14)

outname = figpath + \
    "OISST_ERA5_THFF_lag%i_mon%02i_sigmask.png" % (ilag+1, imon+1)
plt.savefig(outname, dpi=150, bbox_inches='tight', transparent=True)

#%% Examine Breakdown by Denominator and Numerator

imon = 1
ilag = 0

fig, axs, _ = viz.init_orthomap(1, 3, bbplot, figsize=(24, 6.5))

for ii in range(3):
    ax = axs[ii]
    ax = viz.add_coast_grid(ax, bbplot, fill_color='k', proj=proj)
    
    
    if ii == 0:
        plotvar     = cov.isel(month=imon, lag=ilag)/dt
        pname       = "Covariance (SST,FLX)"
        plotmask    = ccmask_inv
        
        vlims = [-45, 45]
        
    elif ii == 1:
        plotvar     = autocov.isel(month=imon, lag=ilag)
        pname       = "Autocovariance (SST)"
        plotmask    = acmask_inv
        
        vlims = [-1, 1]
    else:
        plotvar = thflx_damping.isel(month=imon, lag=ilag)
        pname = "THFLX Feedback"
        plotmask   = sigmask_inv
        
        vlims = [-45, 45]

    pcm = ax.pcolormesh(lon, lat, plotvar, transform=proj, cmap='cmo.balance',
                        vmin=vlims[0], vmax=vlims[1])
    
    # Plot the dots
    # Plot significant points
    plotmsk1 = plotmask[imon, ilag, :, :]
    # Had to reverse and transpose to lon x lat
    viz.plot_mask(plotvar.lon, plotvar.lat, plotmsk1.T, reverse=True, ax=ax, proj=proj, geoaxes=True,
                  markersize=.4, color='gray')

    viz.hcbar(pcm, ax=ax, pad=0.01)

    ax.set_title(pname, fontsize=16)
plt.suptitle("Month %s, Lag %i" % (mons3[imon], ilag+1), fontsize=42)

outname = figpath + \
    "OISST_ERA5_THFF_bycomponent_lag%i_mon%02i_sigmask.png" % (ilag+1, imon+1)
plt.savefig(outname, dpi=150, bbox_inches='tight', transparent=True)