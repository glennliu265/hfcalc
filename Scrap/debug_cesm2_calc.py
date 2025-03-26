#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The CESM2 HFF is exactly the same as the CESm1 HFF, as calculated through 
[calc_enso_general]. What is going on? Let's try to figure out...

Created on Thu Jun 20 17:20:38 2024

@author: gliu
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import time
import sys
import cartopy.crs as ccrs
import glob

#%%
# stormtrack
amvpath = "/home/glliu/00_Scripts/01_Projects/00_Commons/" # amv module
scmpath = "/home/glliu/00_Scripts/01_Projects/01_AMV/02_stochmod/stochmod/model/" # scm module

sys.path.append(amvpath)
sys.path.append(scmpath)

from amv import proc,viz
import scm
import amv.loaders as dl


#%% Stop forward, maybe from teh raw data

lonf = 330
latf = 50


dp = "/stormtrack/data4/glliu/01_Data/CESM2_PiControl/FCM/atm/TS/"
nc = "b.e21.B1850.f09_g17.CMIP6-piControl.001.cam.h0.TS.040001-049912.nc"

ds_cesm2 = xr.open_dataset(dp+nc).sel(lon=lonf,lat=latf,method='nearest').load()


dp = "/vortex/jetstream/climate/data1/yokwon/CESM1_LE/downloaded/atm/proc/tseries/monthly/TS/"
nc = "b.e11.B1850C5CN.f09_g16.005.cam.h0.TS.040001-049912.nc"
ds_cesm1 = xr.open_dataset(dp+nc).sel(lon=lonf,lat=latf,method='nearest').load()


#%% Get Timeseries


tscesm2 = ds_cesm2.TS

tscesm1 = ds_cesm1.TS


#%%


fig,ax = plt.subplots(1,1)
ax.plot(tscesm2,label='cesm2')
ax.plot(tscesm1,label='cesm1')
ax.legend()
plt.show()

#%% Step backwards, maybe starting from the ENSO removal step

# Load CESM1
dp = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_PIC_SLAB/02_ENSOREM/FULL/"
nc = "ENSOREM_TS_lag1_pcs2_monwin3.npz"
ld = np.load(dp+nc,allow_pickle=True)

tsacesm1 = ld['TS'] # (21576, 192, 288)
lat1     = ld['lat']
lon1     = ld['lon']


klon,klat   = proc.find_latlon(lonf,latf,lon1,lat1)
cesm1pt     = tsacesm1[:,klat,klon]


cesm1pt  = cesm1pt.reshape(1798,12)
 
cesm1pt  = cesm1pt[:-200,:] # [1598,12], drop years 2000-2200
  

# Load CESM2
dp       = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/anom/"
nc       = "cesm2_pic_TS_manom_detrend1_1to2000.nc"
dstsa    = xr.open_dataset(dp+nc).sel(lon=lonf,lat=latf,method='nearest').load()
dstsa    = dstsa.sel(time=slice('0403-01-01','2000-12-31'))
cesm2pt  = dstsa.TS.values.reshape(1598,12)

#%% Plot Again


fig,ax = plt.subplots(1,1)
ax.plot(cesm2pt.flatten(),label='cesm2')
ax.plot(cesm1pt.flatten(),label='cesm1')
ax.legend()
plt.show()








