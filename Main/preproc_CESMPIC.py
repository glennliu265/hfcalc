#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Quickly Preprocess CESM1-PIC data to match preproc_ncep_ncar.py
for debugging

Created on Fri May  6 11:48:47 2022

@author: gliu
"""

import numpy as np
import xarray as xr

import sys
stormtrack = 0

if stormtrack:
    sys.path.append("/home/glliu/00_Scripts/01_Projects/00_Commons/")
    sys.path.append("/home/glliu/00_Scripts/01_Projects/01_AMV/02_stochmod/stochmod/model/")
else:
    sys.path.append("/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/03_Scripts/stochmod/model/")
    sys.path.append("/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/03_Scripts/")
from amv import proc,viz
import scm


datpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/01_Data/CESM_proc/"
ncnames  = ["TS_PIC_FULL.nc",'NHFLX_PIC_FULL.nc']
outpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/"

limask = np.load("/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/01_Data/CESM-FULL_landicemask360.npy")




ds_ts = xr.open_dataset(datpath+ncnames[0])



ts                 = ds_ts.TS.values
nyr,nmon,nlat,nlon = ts.shape
ts                 = ts.reshape(nyr*nmon,nlat,nlon)
ts                 *= limask


times = xr.cftime_range(start='0400-01-01',periods=ts.shape[0],freq="MS")

dataset_name = "CESM1_FULL_PIC"
savename     = "%s%s_ts.nc" % (outpath,dataset_name)
proc.numpy_to_da(ts,times,ds_ts.lat.values,ds_ts.lon.values,"ts",savenetcdf=savename)



# Do the same for Qnet
ds_ts = xr.open_dataset(datpath+ncnames[1])

ts                 = ds_ts.NHFLX.values
nyr,nmon,nlat,nlon = ts.shape
ts                 = ts.reshape(nyr*nmon,nlat,nlon)
ts                 *= limask


times = xr.cftime_range(start='0400-01-01',periods=ts.shape[0],freq="MS")

dataset_name = "CESM1_FULL_PIC"
savename     = "%s%s_qnet.nc" % (outpath,dataset_name)
proc.numpy_to_da(ts,times,ds_ts.lat.values,ds_ts.lon.values,"qnet",savenetcdf=savename)






