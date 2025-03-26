#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""'

Debug HF Calc


Created on Mon Feb 12 12:10:33 2024

@author: gliu
"""


import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import sys
from tqdm import tqdm
import copy
import glob
import matplotlib as mpl

#%% Import Custom Modules
amvpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/03_Scripts/" # amv module
scmpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/03_Scripts/stochmod/model/"

sys.path.append(amvpath)
sys.path.append(scmpath)

from amv import proc,viz
import scm
import amv.loaders as dl
import yo_box as ybx


#%%

datpath     = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/03_reemergence/01_Data/proc/CMIP6/proc/"
figpath     = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/02_Figures/20240215/"
proc.makedir(figpath)

#%%


dpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/CESM-HTR-RCP85/SST_Recalc/"
nc    = "htr_hfdamping_ensorem1_detrend1_1920to2005_crop_ens42.nc"
ds    = xr.open_dataset(dpath+nc)


#%%
mons3     = proc.get_monstr()
lonf,latf = 330,50
dspt      = ds.sel(lon=lonf,lat=latf,method='nearest') d
fig,ax = viz.init_monplot(1,1)
for ll in range(3):
    ax.plot(mons3,dspt.qnet_damping.isel(lag=ll),label="Lag %i" % (ll+1))
ax.legend()

#%%

 mons3     = proc.get_monstr()
 lonf,latf = 330,50
 dspt      = ds.sel(lon=lonf,lat=latf,method='nearest')
 fig,ax = viz.init_monplot(1,1)
 for ll in range(3):
     ax.plot(mons3,dspt.autocov.isel(lag=ll),label="Lag %i" % (ll+1))
 ax.legend()
 
 
 
#%% Check on stormtrack

ncnew = "CESM1_htr_ts_ens41.nc"
ncp   = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_HTR/"

ncold = "CESM1_htr_TS_ens41.nc"
ncpold = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_HTR/fluxes_masked/"


dstnew = xr.open_dataset(ncp+ncnew) + 273.15
dstold = xr.open_dataset(ncpold+ncold)

#%%

#dsin = dstnew

dstnew.isel(time=0).ts.plot(),plt.show()

dstold.isel(time=0).TS.plot(),plt.show()

