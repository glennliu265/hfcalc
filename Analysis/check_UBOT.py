#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Examine lower level of U with other variables

Created on Thu Aug 18 11:00:39 2022

@author: gliu
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import glob
import sys
from tqdm import tqdm 

dp1  = "/vortex/jetstream/climate/data1/yokwon/CESM1_LE/downloaded/atm/proc/tseries/monthly/"
dp2  = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/00_Commons/CESM1_LE/"

nc1  = "b.e11.B20TRC5CNBDRD.f09_g16.002.cam.h0.U.192001-200512.nc"
nc2  = "b.e11.B20TRC5CNBDRD.f09_g16.002.cam.h0.V.192001-200512.nc"
nc3  = "b.e11.B20TRC5CNBDRD.f09_g16.002.cam.h0.U10.192001-200512.nc"
nc4  = "b.e11.B20TRC5CNBDRD.f09_g16.002.cam.h0.UBOT.192001-200512.nc"

vns   = ["U","V","U10","UBOT"]
dps   = [dp1,dp1,dp2,dp2]
ncs   = [nc1,nc2,nc3,nc4]


ds_all = []
nvar = len(vns)
for v in range(nvar):
    
    ds = xr.open_dataset("%s%s/%s" % (dps[v],vns[v],ncs[v]))
    
    
    print(ds)
    if vns[v] == "UBOT": # UBOT is called U
        ds_all.append(ds.U)
    elif vns[v] in ["U","V"]:
        ds_all.append(ds[vns[v]].isel(lev=-1)) # Select Bottom Level
        
    else:
        ds_all.append(ds[vns[v]])
    
    
Umod     = np.sqrt(ds_all[0]**2 + ds_all[1]**2)

U_m_UBOT    = ds_all[3]-ds_all[0]
Umod_m_U10  = ds_all[2] - Umod

UUBOT_mean = U_m_UBOT.mean('time')
UmodU10_mean  = Umod_m_U10.mean('time')

