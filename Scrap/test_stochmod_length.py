#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Visualize Differences in ACF based on stochastic model length

Copied upper section of viz_stochmod_output_obs_scrap on 2025.04.08

Created on Tue Apr  8 12:37:47 2025

@author: gliu
"""


import xarray as xr
import numpy as np
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import scipy as sp

import matplotlib.patheffects as PathEffects

import matplotlib.pyplot as plt
import sys
import glob
import os

import tqdm
import time

#from cmcrameri import cm

#%%

# local device

amvpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/03_Scripts/" # amv module
scmpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/03_Scripts/stochmod/model/"

sys.path.append(amvpath)
sys.path.append(scmpath)

from amv import proc,viz
import scm
import amv.loaders as dl
import cvd_utils as cvd

#%% Plotting Inputs

# Set Plotting Options
darkmode = False
if darkmode:
    dfcol = "w"
    bgcol = np.array([15,15,15])/256
    sp_alpha = 0.05
    transparent = True
    plt.style.use('dark_background')
    mpl.rcParams['font.family']     = 'Avenir'
else:
    dfcol = "k"
    bgcol = "w"
    sp_alpha = 0.75
    transparent = False
    plt.style.use('default')

bboxplot    = [-80, 0, 20, 65]
mpl.rcParams['font.family'] = 'Avenir'
mons3       = proc.get_monstr(nletters=3)

fsz_tick    = 18
fsz_axis    = 20
fsz_title   = 16

rhocrit     = proc.ttest_rho(0.05, 2, 86)

proj        = ccrs.PlateCarree()

#%% Load some other things to plot


# Load Sea Ice Masks
dpath_ice = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/NATL_proc_obs/"
nc_masks = dpath_ice + "OISST_ice_masks_1981_2020.nc"
ds_masks = xr.open_dataset(nc_masks).load()

# Load AVISO
dpath_aviso = dpath_ice + "proc/"
nc_adt = dpath_aviso + "AVISO_adt_NAtl_1993_2022_clim.nc"
ds_adt = xr.open_dataset(nc_adt).load()
cints_adt = np.arange(-100, 110, 10)

#%% Set some paths

figpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/03_reemergence/02_Figures/20250409/"



# Analysis Settings
bboxreg        = [-50,-10,50,60]
bbname         = "Yeager 2015"
lags    = np.arange(61)
#%% Load the data


# Load ERA5
dpath_era   = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/NATL_proc_obs/"
ncera       = dpath_era + "ERA5_sst_NAtl_1979to2024.nc" #"ERA5_sst_NAtl_1979to2021.nc"
dsera       = xr.open_dataset(ncera).load().sst


# Load stochastic model data (1kyr run, with entrainment forcing) -------------
sm_output_path = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/03_reemergence/01_Data/sm_experiments/"
expnames       = ["SST_Obs_Pilot_00_Tdcorr0_qnet","SST_Obs_Pilot_00_Tdcorr1_qnet"]
expnames_sm    = ["Stochastic Model (with entrainment forcing)","Stochastic Model (entrainment damping only"]
#ds_all         = [dl.load_smoutput(ex,sm_output_path) for ex in expnames]



# Write a loading function
def load_sm_regavg(expname,sm_output_path,bboxreg):
    
    # Crop to region
    ds    = dl.load_smoutput(expname,sm_output_path,load=False)
    dsreg = proc.sel_region_xr(ds,bboxreg)
    st    = time.time()
    dsreg = dsreg.load()
    print("Loaded region in %.2fs" % (time.time()-st))
    
    # Take the Region Average
    dsreg_avg = proc.area_avg_cosweight(dsreg)
    return dsreg_avg

expname = expnames[0]
dsreg_avg = load_sm_regavg(expname,sm_output_path)

#%% Compute ACF for ERA5

#from 
dsreg_era = proc.sel_region_xr(dsera,bboxreg)
dsreg_era = proc.area_avg_cosweight(dsreg_era)
aavg_anom = proc.xrdeseason(dsreg_era)
aavg_era  = sp.signal.detrend(aavg_anom)

tsm_era   = scm.compute_sm_metrics([aavg_era],lags=lags)
    

#%% Examine sensivity to different amounts of time

# Deseason
dsreg_avg_anom = proc.xrdeseason(dsreg_avg)
aavgs           = dsreg_avg_anom.SST.data # [Run x Year]
aavgs          = sp.signal.detrend(aavgs,axis=0)
aavgs_flat     = aavgs.flatten()


# First, compute the metrics for the whole timeseries...
tsm_stochmod = scm.compute_sm_metrics([aavgs_flat],lags=lags)

#tmaxes_name = ["46 Years (ERA5)","75 Years","100 Years","500 Years","1000 Years",]
#tmaxes      = [552,75*12,100*12,500*12,1000*12]


#%% Sample Years

# Indicate Time Length, Number of Samples
nyrs     = 50
tlen     = nyrs*12

nsamp    = 10000

istarts  = np.arange(0,aavgs_flat.shape[0] - tlen,12) # Note need to sample carefully to start on Jan..
tranges  = []
tseries  = []

for ii in tqdm.tqdm(range(nsamp)):
    #print(istart)
    istart = np.random.choice(istarts,size=1)[0]
    trange = np.arange(istart,istart+tlen)
    
    tranges.append(trange)
    tseries.append(aavgs_flat[trange])
    
st = time.time()
tmetrics = scm.compute_sm_metrics(tseries,lags=lags)
print("Computed Metrics in $.2fs" % (time.time()-st))

#%% Plot ACF for a given Month


istart_cut = 7000*12 # Set a cutoff point to  indicate periods before cutoff with specified color


xtks   = lags[::6]
kmonth = 1
fig,ax = plt.subplots(1,1,figsize=(8,4.5),constrained_layout=True)
ax,_   = viz.init_acplot(kmonth,xtks,lags,ax=ax)

acf_samples = []
for ii in range(nsamp):
    plotvar = tmetrics['acfs'][kmonth][ii]
    label = "t= %i to %i" % (tranges[ii][0],tranges[ii][-1])
    
    if tranges[ii][0] > istart_cut:
        cplot = 'cornflowerblue'
    else:
        cplot = 'firebrick'
    
    ax.plot(lags,plotvar,label="",
            alpha=.5,c=cplot)
    acf_samples.append(plotvar)


plotvar = np.array(acf_samples).mean(0)
ax.plot(lags,plotvar,label="Stochastic Model, %i Sample Average" % nsamp,
        lw=2.5,color='blue')    

plotvar = tsm_era['acfs'][kmonth][0]
ax.plot(lags,plotvar,label="ERA5",lw=2.5,color="k")

plotvar = tsm_stochmod['acfs'][kmonth][0]
ax.plot(lags,plotvar,label="Stochastic Model (10kyr)",lw=2.5,color="orange")

title = "%s ACF, %s\nNumber Years: %i, Number Samples: %i" % (mons3[kmonth],bbname,
                                                    nyrs,nsamp)
ax.set_title(title)
ax.legend()

ax.set_ylim([-0.75,1.25])

figname = '%sSM_Sample_Test_%s_month%0i_nyr%04i_nsamp%i_tr1%i.png' % (figpath,expname,kmonth+1,nyrs,nsamp,tranges[0][0])

plt.savefig(figname,dpi=150)

#%%


# Take T2 of the samples
t2_samples = proc.calc_T2(np.array(acf_samples),axis=1)
plt.hist(t2_samples)

id_min_to_max = np.argsort(t2_samples)
max10 = id_min_to_max[-10:][::-1]
min10 = id_min_to_max[:10]

sstmax10 = []
trmax10 = []
sstmin10 = []
trmin10 = []
for ii in range(10):
    
    trangemax = tranges[max10[ii]]
    sstmax10.append(aavgs_flat[trangemax])
    trmax10.append(trangemax)
    
    trangemin = tranges[min10[ii]]
    sstmin10.append(aavgs_flat[trangemin])
    trmin10.append(trangemin)
    

# Also get variance of the samples
tsarr = np.array(tseries)
tsvar = np.var(tseries,1)
    
#%% PLot the Timepseries

fig,axs= plt.subplots(2,1,constrained_layout=True,figsize=(12.5,4.5),sharex=True)

alph = 1
for ii in range(10):
    
    # Plot max
    ax    = axs[0]
    label = "Rank %i (%i to %i)"  % (ii+1,trmax10[ii][0],trmax10[ii][-1])
    ax.plot(sstmax10[ii],label=label,alpha=alph)
    
    # Plot min
    ax    = axs[1]
    label = "Rank %i (%i to %i)"  % (ii+1,trmin10[ii][0],trmin10[ii][-1])
    ax.plot(sstmin10[ii],label=label,alpha=alph)
for ax in axs:
    #ax.legend()
    ax.set_xlim([0,600])
    
    
#%% Check T2 and Variance Relationship

fig,ax = plt.subplots(1,1,constrained_layout=True)

ax.scatter(tsvar,t2_samples)

ax.set_xlabel("Variance (degC$2$)")
ax.set_ylabel("Feb T2 (Months)")
#ax.set_aspect('equal')
    

#%% Examine visually the characteristics of the most persistent timeseries



# =====================================================================
##%% Test Sampling Above, but split into late and early stochastic model
# =====================================================================
# Just get the first simulation
aavg_flat_subset = aavgs_flat[(1000*12):((1000*12)*2)]

# Split into 100 year chunks
tlen_split = 100*12
chunks     = []
for ii in range(10):
    iirange = [(tlen_split*ii),(tlen_split*(ii+1))]
    chunks.append(aavg_flat_subset[iirange[0]:iirange[1]])

tsm_chunks = scm.compute_sm_metrics(chunks,lags=lags)

#%% Visualize Timeseries Chunks
fig,ax = plt.subplots(1,1,figsize=(12,4.5),constrained_layout=True)
ax,_   = viz.init_acplot(kmonth,xtks,lags,ax=ax)



acf_samples = []
for ii in range(10):
    plotvar = tsm_chunks['acfs'][kmonth][ii]
    iirange = [(tlen_split*ii),(tlen_split*(ii+1))]
    label = "t= %i to %i" % (iirange[0],iirange[1])
    
    alpha = 0.1 + 0.9* ii/10
    cplot = 'cornflowerblue'
    
    
    ax.plot(lags,plotvar,label=label,lw=3,alpha=alpha)
    acf_samples.append(plotvar)


plotvar = np.array(acf_samples).mean(0)
ax.plot(lags,plotvar,label="Stochastic Model, %i Sample Average" % nsamp,
        lw=2.5,color='blue')    

plotvar = tsm_era['acfs'][kmonth][0]
ax.plot(lags,plotvar,label="ERA5",lw=2.5,color="k")

plotvar = tsm_stochmod['acfs'][kmonth][0]
ax.plot(lags,plotvar,label="Stochastic Model (10kyr)",lw=2.5,color="orange")

title = "%s ACF, %s\nNumber Years: %i, Number Samples: %i" % (mons3[kmonth],bbname,
                                                    nyrs,nsamp)
ax.set_title(title)
ax.legend(ncol=5)

ax.set_ylim([-0.75,1.25])

figname = '%sSM_Sample_Test_%s_month%0i_Chunks.png' % (figpath,expname,kmonth+1,)

plt.savefig(figname,dpi=150)


    
    
#%% Increase Sampling Length

# Indicate the number of years
nyrs      = np.array([46,100,150,200,500,1000,5000,10000])
tmaxes    = nyrs*12
aavg_tmax = [aavgs_flat[:tmax] for tmax in tmaxes]

tsm_tmax  = scm.compute_sm_metrics(aavg_tmax,lags=lags)

#%% Make the plots

fig,ax = plt.subplots(1,1,figsize=(12,4.5),constrained_layout=True)
ax,_   = viz.init_acplot(kmonth,xtks,lags,ax=ax)

acf_samples = []
for ii in range(len(nyrs)):
    plotvar = tsm_tmax['acfs'][kmonth][ii]
    
    #iirange = [(tlen_split*ii),(tlen_split*(ii+1))]
    label = "Years 0 to %i" % (nyrs[ii])#"t= %i to %i" % (iirange[0],iirange[1])
    
    alpha = 0.15 + 0.9* ii/len(nyrs)
    cplot = 'cornflowerblue'
    
    ax.plot(lags,plotvar,label=label,lw=3,alpha=alpha)
    acf_samples.append(plotvar)

plotvar = tsm_era['acfs'][kmonth][0]
ax.plot(lags,plotvar,label="ERA5",lw=2.5,color="k")

plotvar = tsm_stochmod['acfs'][kmonth][0]
ax.plot(lags,plotvar,label="Stochastic Model (10kyr)",lw=2.5,color="orange")

title = "%s ACF, %s" % (mons3[kmonth],bbname)
ax.set_title(title)
ax.legend(ncol=5)

ax.set_ylim([-0.75,1.25])

figname = '%sSM_Sample_Test_%s_month%0i_Tmax_Limit.png' % (figpath,expname,kmonth+1,)

plt.savefig(figname,dpi=150)


#%% Make an iterative plot

for ic,ii in tqdm.tqdm((enumerate(range(50,500,10)))):
    tmax    = ii * 12
    aavg_in = aavgs_flat[:tmax]
    tsm_tmax  = scm.compute_sm_metrics([aavg_in],lags=lags)
    
    fig,ax = plt.subplots(1,1,figsize=(12,4.5),constrained_layout=True)
    ax,_   = viz.init_acplot(kmonth,xtks,lags,ax=ax)
    
    # Plot Tmax
    plotvar = tsm_tmax['acfs'][kmonth][0]
    label = "Years 0 to %i" % (ii)#"t= %i to %i" % (iirange[0],iirange[1])  
    ax.plot(lags,plotvar,label=label,lw=3,alpha=alpha)

    plotvar = tsm_era['acfs'][kmonth][0]
    ax.plot(lags,plotvar,label="ERA5",lw=2.5,color="k")

    plotvar = tsm_stochmod['acfs'][kmonth][0]
    ax.plot(lags,plotvar,label="Stochastic Model (10kyr)",lw=2.5,color="orange")

    title = "%s ACF, %s, nyrs = %03i" % (mons3[kmonth],bbname,ii)
    ax.set_title(title)
    ax.legend(ncol=5)

    ax.set_ylim([-0.75,1.25])

    figname = '%sSM_Sample_Test_%s_month%0i_Tmax_Limit_nyr%03i.png' % (figpath,expname,kmonth+1,ii)

    plt.savefig(figname,dpi=150)





    
    
    
