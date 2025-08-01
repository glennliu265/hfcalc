#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Compare Damping Computed from [calc_enso_general] for CESM1 and CESM2 
PiControl Runs

Copied upper section from [viz_hfdamping_CESM1_PIC]




Created on Thu Jun 20 16:28:38 2024

@author: gliu

"""


import sys


from scipy.interpolate import interp1d
from scipy.io import loadmat,savemat
from scipy import signal,stats
from tqdm import tqdm

import cartopy.crs as ccrs
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import calendar as cal

import scm
import time
import cmocean
import copy
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import xarray as xr
import sys
from tqdm import tqdm
import copy
import glob
import time
import cartopy.crs as ccrs
import os


#%% Load Modules

amvpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/03_Scripts/" # amv module
scmpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/03_Scripts/stochmod/model/"

sys.path.append(amvpath)
sys.path.append(scmpath)

from amv import proc,viz
import scm
import amv.loaders as dl
import cvd_utils as cvd


#%% Damping Paths


# CESM2_pic, computed with [calc_enso_general]
dp1          = '/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/hff/'
nc1          = 'cesm2_pic_hfdamping_ensorem1_detrend1_1to2000_crop.nc'
ds_cesm2    = xr.open_dataset(dp1+nc1).load()


# CESM1_pic, computed with ???
dp          = '/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/'
nc          = 'CESM1_PIC-FULL_NHFLX_Damping_lag1.nc'
ds_cesm1    = xr.open_dataset(dp+nc).load()

ds_cesm1 = ds_cesm1.rename(dict(latitude='lat',longitude='lon'))


# Note the values are exactly the same, perhaps something is wrong?

#%% plotting params

figpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/03_reemergence/02_Figures/20240619/"
proc.makedir(figpath)

# Plotting Info
bbplot                      = [-80,0,20,65]
mpl.rcParams['font.family'] = 'Avenir'
proj                        = ccrs.PlateCarree()
mons3                       = proc.get_monstr()

# Font Sizes
fsz_title                   = 20
fsz_ticks                   = 14
fsz_axis                    = 16
fsz_legend                  = 16

#%% Load Sea Ice

#%% Compare it
im = 1

for im in range(12):
    cints   = np.arange(-48,52,4)
    
    inhffs = [ds_cesm2.qnet_damping.isel(lag=0,month=im),ds_cesm1.nhflx_damping.isel(month=im)]
    expnames = ["CESM2","CESM1"]
    
    
    fig,axs,_=viz.init_orthomap(1,2,bbplot,figsize=(12,5.5))
    
    for ii in range(2):
        ax = axs[ii]
        ax = viz.add_coast_grid(ax,bbox=bbplot,fill_color="lightgray",line_color="k")
        
        pv = inhffs[ii]
        #pcm = ax.pcolormesh(pv.lon,pv.lat,pv,transform=proj,cmap='cmo.balance',vmin=-40,vmax=40)
        cf  = ax.contourf(pv.lon,pv.lat,pv,levels=cints,transform=proj,cmap='cmo.balance',)
        cl  = ax.contour(pv.lon,pv.lat,pv,levels=cints,transform=proj,colors="dimgray",linewidths=0.75)
        ax.clabel(cl,cints[::2],fontsize=fsz_ticks)
        #fig.colorbar(pcm,ax=ax)
        
        ax.set_title(expnames[ii],fontsize=fsz_title)
        
    
    cb = viz.hcbar(pcm,ax=axs.flatten(),fraction=0.045)
    savename = "%sHFF_CESM1_CESM2_Lag_Comparison_mon%03i.png" % (figpath,im+1)
    plt.savefig(savename,dpi=150,bbox_inches="tight")
    #hffarr = [hf.values for hf in inhffs]




#%% Plot Compare at a plont

lonf = -30
latf = 50

inhffs = [ds_cesm2.qnet_damping.isel(lag=0),ds_cesm1.nhflx_damping]
hffpt = [hf.sel(lon=lonf,lat=latf,method='nearest') for hf in inhffs]




fig,ax = viz.init_monplot(1,1)

for ii in range(2):
    ax.plot(mons3,hffpt[ii],label=expnames[ii])
ax.legend()

#%%

im = 0

inhffs = [ds_cesm2.qnet_damping.isel(lag=0),ds_cesm1.nhflx_damping]
diffs_pw = inhffs[1] - inhffs[0]


plt.pcolormesh(diffs_pw[im,:,:]),plt.colorbar()

#%% Data Paths

datpath     = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/"
lipath      = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/01_Data/model_input/"
llpath      = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/01_Data/model_input/"
figpath     = '/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/02_Figures/20220824/'
proc.makedir(figpath)

mconfigs    = ["PIC-SLAB","PIC-FULL"]
dofs        = [898 - 1 - 2 - 2, 1898 - 1 - 2 - 2] 

lags        = [1,2,3]
bboxplot    =  [-80,0,5,62]
mons3       = [viz.return_mon_label(m,nletters=3) for m in np.arange(1,13)]

