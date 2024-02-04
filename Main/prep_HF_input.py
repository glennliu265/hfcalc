#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

-----------------------------
Make Heat Flux Feedback Input
-----------------------------

Takes a given input (from [combine_damping_CESM1LE.py], or [calc_enso_general], etc)
Flips Longitude and Crops to Region
Applies significance testing, and generates output ready for input into the stochastic model

Relies on parameters indicated in [stochmod_params] (or maybe generate them here and copy em over)
Trying to generalize function in [prep_HF]

Created on Thu Feb  1 13:38:30 2024
@author: gliu

"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob

import tqdm
import time

#%% Import Custom Modules

# Import AMV Calculation
amvpath = "/home/glliu/00_Scripts/01_Projects/00_Commons/" # amv module
sys.path.append(amvpath)
from amv import proc,viz
import amv.loaders as dl

# Import stochastic model scripts
sys.path.append("/home/glliu/00_Scripts/01_Projects/01_AMV/02_stochmod/stochmod/model/")
import scm

# Import Hf Calc params
hfpath  = "/stormtrack/home/glliu/00_Scripts/01_Projects/01_AMV/01_hfdamping/hfcalc/" # hfcalc module 
sys.path.append(hfpath)
import hfcalc_params as hp

#%% Print and indicate parameter set, paths

# Region Crop

bbox_crop = [-90,20,0,90]  # Preprocessing box
dof       = 83#hp.dofs_cesm['HTR_FULL'] #None # DOF to use, can manually set
print("Using %i dofs" % (dof))

# Parameter Set
print("Parameter Sets: " + str(hp.hff_names))
setname = "nomasklag1"
setdict = hp.hff_sets[setname]


# Save a no-test version

# Input (nc file with LHFLX damping, crosscorr, and autocorr with [ens x mon x lag x lat x lon360])
ncname = "CESM1_htr_LHFLX_damping_ensorem1_detrend1_1920to2005_allens.nc"
ncpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_HTR/LHFLX_damping/"
vname  = "LHFLX_damping"

# Output Path
datpath = ""
outpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/03_reemergence/proc/"#"#model_input/damping/" # Note this was giving me permission erorrs


#%% Load HFF

# Load, Flip Lon, Crop to Region
st      = time.time()
ds      = xr.open_dataset(ncpath+ncname)
ds      = proc.lon360to180_xr(ds)
ds_reg  = proc.sel_region_xr(ds,bbox_crop).load()
print("Loaded in %.2fs" % (time.time()-st))

# All of size [ens x mon x lag x lat x lon360]
hff     = ds_reg[vname].values               
rsst    = ds_reg['sst_autocorr'].values      
rflx    = ds_reg['sst_flx_crosscorr'].values

#%% Do a case where no testing is applied



#%% Debugging Plots

# ds.LHFLX_damping.isel(ens=0,lag=0,month=0).plot(),plt.show()
# ds_reg.LHFLX_damping.isel(ens=0,lag=0,month=0).plot(),plt.show()
# ds_reg.sst_autocorr.isel(ens=0,lag=0,month=0).plot(),plt.show()
# ds_reg.sst_flx_crosscorr.isel(ens=0,lag=0,month=0).plot(),plt.show()


# plt.pcolormesh(sigmask[0,0,0,:,:]),plt.show()
# plt.pcolormesh(sigmask[2,0,0,:,:]),plt.show()
# plt.pcolormesh(freq_success[2,0,0,:,:]),plt.show()


# plt.pcolormesh(hff[:,0,0,:,:].mean(0),vmin=-35,vmax=35,cmap="RdBu_r"),plt.colorbar(),plt.show()
# plt.pcolormesh(dampingmasked[:,0,0,:,:].mean(0),vmin=-35,vmax=35,cmap="RdBu_r"),plt.colorbar(),plt.show()
# plt.pcolormesh(dampingmasked[:,0,0,:,:].mean(0)),plt.colorbar(),plt.show()
 
 
# plt.pcolormesh(dampingmasked[2,0,0,:,:]),plt.colorbar(),plt.show()


#plt.pcolormesh(dampingmasked[:,0,0,:,:].mean(0)),plt.show()

#%% Perform Significance Testing and make the mastk

# Compute and Apply Mask
st = time.time()
dampingmasked,freq_success,sigmask=scm.prep_HF(hff,rsst,rflx,
                          setdict['p'],setdict['tails'],dof,setdict['method'],
                          returnall=True) #expects, [month x lag x lat x lon], should generalized with ensemble dimension?
print("Completed significance testing in %.2fs" % (time.time()-st))


# Take Lag
dampinglag = dampingmasked[:,:,setdict['sellags'],:,:].squeeze() # [Ens x Month x Lat x Lon]

# Transpose to [Month x Ens x Lat x Lon]
dampingout = dampinglag.transpose(1,0,2,3)
#%% Transpose to [Month x Ens x Lat x Lon] and prepare to save

# Save output
cdict    = {'mon':np.arange(1,13,1),
              'ens':np.arange(1,43,1),
              'lat':ds_reg.lat.values,
              'lon':ds_reg.lon.values}
daout    = xr.DataArray(dampingout,coords=cdict,dims=cdict,name='damping')
edict    = {"damping":{"zlib":True}}

savename = "%sCESM1_HTR_FULL_%s_%s.nc" % (outpath,vname,setname)
daout.to_netcdf(savename,encoding=edict)

# Save EnsAvg Version
daout_eanvg = daout.mean('ens')
savename = "%sCESM1_HTR_FULL_%s_%s_EnsAvg.nc" % (outpath,vname,setname)
daout_eanvg.to_netcdf(savename,encoding=edict)


#%% Assumes timeseries have no correlation (can eventually implement effective dof)


# Assumes all points have same dof
#corrthres = proc.ttest_rho(setdict['p'],setdict['tails'],dof)





