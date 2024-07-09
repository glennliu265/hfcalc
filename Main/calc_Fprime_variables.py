#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Storing inputs for calc_Fprime here...

Created on Wed Jul  3 11:17:34 2024

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

#%% Set Paths

# Note: For Eprime computation, just switch to the correct flux and damping!
#Eprime     = False # Set to True to Compute E' instead of F'


stormtrack   = 1

#%% CESM2 PiC

# Indicate inputs
datname      = "cesm2_pic"
lensflag     = False # Set to True for lens datasets/to detrend with ensemble average
outvar       = "Fprime"  # "Set to Fprime by default, but LHFLX for Eprime calculations..."
outpath      = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/proc/"
regstr       = "NAtl"

# Mixed Layer Depth --> [h: time x lat x lon180]
mldpath      = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/03_reemergence/proc/model_input/mld/"
mldnc        = "cesm2_pic_HMXL_NAtl_0200to2000.nc"
mldname      = "h"

# Net Heat Flux (Positive Upwards) --> [qnet: time x lat x lon180]
flxpath      = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/proc/"
flxnc        = "cesm2_pic_qnet_NAtl_0200to2000.nc"
flxname      = 'qnet'

# SST
sstpath      = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/proc/"
sstnc        = "cesm2_pic_TS_NAtl_0200to2000.nc"
sstname      = 'TS'

# Damping 
damppath     = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/hff/qnet_damping/"
dampnc       = "cesm2_pic_hfdamping_NAtl_0200to2000_ensorem1_detrend1.nc"
dampname     = 'qnet_damping'
ilag         = 0 # Indicate which lag to select

# Damping Information and roll options
dampstr      = "CESM2PiCqnetDamp"
nroll        = 0 # Amount to roll lbd*T' term
rollstr      = "nroll%0i"  % nroll
convert_wm2  = False # Convert hff to wm2

#%% CESM1 LENS Regridded (5 Deg), Qnet Damping

# Indicate inputs
datname      = "cesm1le_htr_5degbilinear"
lensflag     = True # Set to True for lens datasets/to detrend with ensemble average
outvar       = "Fprime"  # "Set to Fprime by default, but LHFLX for Eprime calculations..."
outpath      = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/proc/"
regstr       = "Global"

# Mixed Layer Depth --> [h: time x lat x lon180]
mldpath      = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/03_reemergence/proc/model_input/mld/"
mldnc        = "cesm1_htr_5degbilinear_HMXL_Global_1920to2005.nc"#"cesm1_htr_5degbilinear_HMXL_Global_1920to2005.nc"
mldname      = "h"

# Net Heat Flux (Positive Upwards) --> [qnet: time x lat x lon180]
flxpath      = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/proc/"
flxnc        = "cesm1_htr_5degbilinear_qnet_Global_1920to2005.nc"
flxname      = 'qnet'

# SST
sstpath      = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/proc/"
sstnc        = "cesm1_htr_5degbilinear_TS_Global_1920to2005.nc"
sstname      = 'TS'

# Damping 
damppath     = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/hff/qnet_damping/"
dampnc       = "cesm1_htr_5degbilinear_hfdamping_Global_1920to2005_ensorem1_detrend1.nc"
dampname     = 'qnet_damping'
ilag         = 0 # Indicate which lag to select

# Damping Information and roll options
dampstr      = "cesm1le5degqnet"
nroll        = 0 # Amount to roll lbd*T' term
rollstr      = "nroll%0i"  % nroll
convert_wm2  = False # Convert hff to wm2

#%% CESM1 LENS Regridded (5 Deg), LHFLX Damping

# Indicate inputs
datname      = "cesm1le_htr_5degbilinear"
lensflag     = True # Set to True for lens datasets/to detrend with ensemble average
outvar       = "Eprime"  # "Set to Fprime by default, but LHFLX for Eprime calculations..."
outpath      = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/proc/"

# Mixed Layer Depth --> [h: time x lat x lon180]
mldpath      = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/03_reemergence/proc/model_input/mld/"
mldnc        = "cesm1_htr_5degbilinear_HMXL_Global_1920to2005.nc"#"cesm1_htr_5degbilinear_HMXL_Global_1920to2005.nc"
mldname      = "h"

# Net Heat Flux (Positive Upwards) --> [qnet: time x lat x lon180]
flxpath      = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/proc/"
flxnc        = "cesm1_htr_5degbilinear_LHFLX_Global_1920to2005.nc"
flxname      = 'LHFLX'

# SST
sstpath      = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/proc/"
sstnc        = "cesm1_htr_5degbilinear_TS_Global_1920to2005.nc"
sstname      = 'TS'

# Damping 
damppath     = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/hff/LHFLX_damping/"
dampnc       = "cesm1_htr_5degbilinear_hfdamping_Global_1920to2005_ensorem1_detrend1.nc"
dampname     = 'LHFLX_damping'
ilag         = 0 # Indicate which lag to select

# Damping Information and roll options
dampstr      = "cesm1le5degLHFLX"
nroll        = 0 # Amount to roll lbd*T' term
rollstr      = "nroll%0i"  % nroll
convert_wm2  = False # Convert hff to wm2

# Conversion Factors
dt          = 3600*24*30
cp0         = 3996
rho         = 1026

