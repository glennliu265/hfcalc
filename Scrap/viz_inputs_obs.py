#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Visualize the stochastiuc model inputs for observations

Created on Wed Mar 26 10:22:26 2025

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

figpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/03_reemergence/02_Figures/20250402/"
proc.makedir(figpath)

#%% Indicate Regions
bbname = ["Yeager 2015",
          "SPG West",
          "SPG East",
          "West REI Maxima",
          "East REI Maxima",
          ]
bbsel = ([-50,-10,50,60],
         [-42,-25,50,60],
         [-35,-10,55,63],
         [-50,-40,50,55],
         [-25,-10,50,55]
         )
bbcol = ('k',
         'orange',
         'violet',
         'navy',
         'firebrick')

bbsty = (
    'solid',
    'dashed',
    'solid',
    'dashed',
    'dotted'
    )
nregs = len(bbname)


#%% Load Inputs

input_path = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/03_reemergence/01_Data/proc/model_input/"

ptype   = [
    "forcing",
    "forcing",
    "damping",
    "damping",
    "mld"
    ]

ncnames = [
    "ERA5_Fprime_THFLX_std_pilot.nc",
    "ERA5_Fprime_QNET_std_pilot.nc",
    "ERA5_thflx_damping_pilot.nc",
    "ERA5_qnet_damping_pilot.nc",
    "MIMOC_regridERA5_h_pilot.nc"
    ]

punits = [
    "W m$^{-2}$",
    "W m$^{-2}$",
    "W m$^{-2} \, \degree$C$^{-1}$",
    "W m$^{-2} \, \degree$C$^{-1}$",
    "m" 
    ]

pnames = [
    "Fprime (THFLX)",
    "Fprime (Qnet)",
    "Damping (THFLX)",
    "Damping (Qnet)",
    "MLD"
    
    ]

pnames_fn = [
    "Fprime_THFLX",
    "Fprime_Qnet",
    "Damping_THFLX",
    "Damping_Qnet",
    "MLD"
    
    ]

vnames = [
    "Fprime",
    "Fprime",
    "damping",
    "damping",
    "h"

    ]

nload     = len(ncnames)

ds_inputs = []

for nn in range(nload):
    
    ncname = "%s%s/%s" % (input_path,ptype[nn],ncnames[nn])
    ds     = xr.open_dataset(ncname).load()[vnames[nn]]
    ds_inputs.append(ds)
    
    
#%%
vv          = 4
imon        = 1
bboxspg    = [-80, 0, 30, 65]

vtype = ptype[vv]
plotvar = ds_inputs[vv].isel(mon=imon).squeeze().T

if vtype == "forcing":
    vlms = [0,150]
    cmap = "inferno"
elif vtype == "damping":
    vlms = [-50,50]
    cmap = "cmo.balance"
elif vtype == "mld":
    vlms = [0,400]
    cmap = 'cmo.deep'


fig,ax  = viz.init_regplot(bboxin=bboxspg)



pcm = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar.T,vmin=vlms[0],vmax=vlms[1],
                    transform=proj,cmap=cmap)

ax.set_title(pnames[vv],fontsize=25)

cb = viz.hcbar(pcm,ax=ax)
cb.ax.tick_params(labelsize=16)


#%% Visualize Regional Averages of the parameters

params_byreg = []
aavg_byreg   = []
for rr in range(nregs):
    #print(rr)
    params_reg = [proc.sel_region_xr(ds,bbsel[rr]) for ds in ds_inputs]
    
    aavg       = [proc.area_avg_cosweight(ds) for ds in params_reg]
    
    params_byreg.append(params_reg)
    aavg_byreg.append(aavg)
    

#%% Plot Area Avg Parameters
vv     = 0


for vv in range(5):
    fig,ax = viz.init_monplot(1,1,figsize=(6,4))
    intype = ptype[vv]
    
    if intype == "forcing":
        ylim = [0,85]
    elif intype == "damping":
        ylim = [-5,25]
    elif intype == "mld":
        ylim = [0,300]
        
    
    for rr in range(nregs):
        plotvar = aavg_byreg[rr][vv].squeeze()
        ax.plot(mons3,plotvar,c=bbcol[rr],label=bbname[rr],lw=2.5)
    ax.legend()
    
    ax.set_title(pnames[vv])
        
    ax.set_ylabel("%s [%s]" % (ptype[vv], punits[vv]))
    ax.set_ylim(ylim)
    
    figname = "%sRegional_Avg_Parameters_%s.png" % (figpath,pnames_fn[vv])
    plt.savefig(figname,dpi=150,bbox_inches='tight')




    
#%% Wintertime values


    




    
#%%  Take seasonal averages


ds_savg = [ds.groupby('time.season').mean('time') for ds in ds_inputs]



