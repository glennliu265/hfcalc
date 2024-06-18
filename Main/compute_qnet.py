#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Compute Net Heat Flux (Qnet) from the flux components.
Copied form predict_nasst/calc_nhflx_predictor.py
Works with output from:
    - preproc_CESM2_PiControl

Created on Tue Jun  6 23:07:02 2023

@author: gliu

"""

import numpy as np
import xarray as xr
import sys
import time
import matplotlib.pyplot as plt

# ------------
#%% User Edits
# ------------

# Dataset and path information
datname     = "cesm2_pic" # Name of dataset (for input and output)
datpath     = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/anom/"
outpath     = datpath # Output Path
vname_out   = "qnet" # 
outname     = "%s%s_%s_anom.nc" % (datpath,datname,vname_out)
print("File will be saved as %s" % outname)


# NC files (can be anomalized, assume **Positive Upwards**, all units [W/m2])
nc_fsns         = "%s%s_FSNS_anom.nc" % (datpath,datname) # Incoming Shortwave
nc_flns         = "%s%s_FLNS_anom.nc" % (datpath,datname) # Outgoing Longwave
nc_shflx        = "%s%s_SHFLX_anom.nc" % (datpath,datname) # Outgoing Sensible Heat Flux
nc_lhflx        = "%s%s_LHFLX_anom.nc" % (datpath,datname) # Outgoing Latent Heat Flux
nclist          = [nc_fsns,nc_flns,nc_shflx,nc_lhflx]
vnames          = ["FSNS","FLNS","SHFLX","LHFLX"] # Shortwave first! Names of variables in netCDFs


# Glossary of names to check (later, add CESM1/CESM2, CMIP5, and CMIP6 definitions, etc)
vname_cmip = ("rsus" ,"rlus" ,"rsds" ,"rlds" ,"hfss" ,"hfls","ts")
vname_ncep = ("uswrf","ulwrf","dswrf","dlwrf","shtfl","lhtfl","air")
vname_cesm = ("FSNS","FLNS",None,None,"SHFLX","LHFLX","TS")
vname_new  = ("fsns" ,"flns" ,"fsns" ,"flns" ,"hfss","hfls","ts")
    
    

fsns_names      = ["FSNS"]
flns_names      = ["FLNS"]
shflx_names     = ["SHFLX"]
lhflx_names     = ["LHFLX"]
qnet_names      = ["NHFLX","QNET","qnet",'nhflx']

# -------------------------
#%% Add Everything Together
# -------------------------

lonf  = 330
latf  = 50
check = []
for vv in range(4):
    stv = time.time()
    
    # Load the data
    st = time.time()
    vname = vnames[vv]
    ds    = xr.open_dataset(nclist[vv])[vname].load() # Maybe add region crop before load...
    print("Loaded data in %.2fs" % (time.time()-st))
    
    if vname in fsns_names:
        ds      = ds * -1 # Flip to positive upwards
        print("Flipping sign for Incoming Shortwave")
    
    
    if (vv == 0):
        ds_qnet = ds.copy()
    else:
        ds_qnet = ds_qnet + ds
    
    # Check Append
    check.append(ds.sel(lon=lonf,lat=latf,method='nearest').isel(time=1).values.item())
    
    ds.close()
    print("Completed %s in %.2fs" % (vname,time.time()-st))


print("QNET = -FSNS + (FLNS + SHFLX + LHFLX)" )
print("%.2f = %.2f + (%.2f +  %.2f +  %.2f)" % (np.array(check).sum(),check[0],check[1],check[2],check[3]))
print("Calculated value = %.2f" % (ds_qnet.sel(lon=lonf,lat=latf,method='nearest').isel(time=1).values.item()))

#%% Save output

st = time.time()
ds_qnet = ds_qnet.rename("qnet")
edict   = {vname_out:dict(zlib=True)}
ds_qnet.to_netcdf(outname,encoding=edict)
print("Save output in %.2fs" % (time.time()-st))
