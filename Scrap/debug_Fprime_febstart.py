#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Debug febstart fix in calc_Fprime

Created on Tue Jul  9 15:53:06 2024

@author: gliu
"""


import numpy as np
import xarray as xr
import sys
import time
import matplotlib.pyplot as plt

#%%
# stormtrack
amvpath = "/home/glliu/00_Scripts/01_Projects/00_Commons/" # amv module
sys.path.append(amvpath)
import amv.proc as hf

#%%

dpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/proc/"
ncnew = "cesm1le_htr_5degbilinear_Fprime_timeseries_cesm1le5degqnet_nroll0_Global.nc"
ncold = "cesm1le_htr_5degbilinear_Fprime_timeseries_cesm1le5degqnet_nroll0_NAtl_old.nc"
ncs   = [ncold,ncnew]


dsall = [xr.open_dataset(dpath + nc).load() for nc in ncs] # (time: 1032, ens: 42, lat: 37, lon: 73)

# Compute seasonal cycle
dsscycle = [ds.groupby('time.month').mean('time') for ds in dsall]

#%% Quickly Check seasonal cycle

names = ["Old","New (Fix Febstart)"]

fig,ax = plt.subplots(1,1,figsize=(12,5))
for ii in range(2):
    #ax.plot(dsscycle[ii].mean('ens').Fprime.data[:,22,22],label=names[ii])
    
    ax.plot(np.nanmean(dsscycle[ii].Fprime.data,1)[:,22,22],label=names[ii])
ax.legend()
plt.show()



#fps  = [ ds]