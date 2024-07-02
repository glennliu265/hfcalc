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

#procname = "cesm2_pic_anom"#"cesm2_pic_regioncut"
procname = "cesm1_htr_5degbilinear"

#"cesm2_pic_anom"

# Version (1), Anomalized Flux Components from CESM1
if procname == "cesm2_pic_anom":
    # Dataset and path information
    datname     = "cesm2_pic" # Name of dataset (for input and output)
    datpath     = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/proc/"
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

# Version (2), Region crop full flux fields from [preproc_raw_inputs]
elif procname == "cesm2_pic_regioncut":
    # Dataset and path information
    datname     = "cesm2_pic" # Name of dataset (for input and output)
    datpath     = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/proc/"
    outpath     = datpath # Output Path
    vname_out   = "qnet" # 
    
    outname     = "%s%s_%s_NAtl_0200to2000.nc" % (datpath,datname,vname_out)
    print("File will be saved as %s" % outname)
    
    # NC files (can be anomalized, assume **Positive Upwards**, all units [W/m2])
    nc_fsns         = datpath + "cesm2_pic_FSNS_NAtl_0200to2000.nc"  # Incoming Shortwave
    nc_flns         = datpath + "cesm2_pic_FLNS_NAtl_0200to2000.nc"  # Outgoing Longwave
    nc_shflx        = datpath + "cesm2_pic_SHFLX_NAtl_0200to2000.nc" # Outgoing Sensible Heat Flux
    nc_lhflx        = datpath + "cesm2_pic_LHFLX_NAtl_0200to2000.nc" # Outgoing Latent Heat Flux
    nclist          = [nc_fsns,nc_flns,nc_shflx,nc_lhflx]
    vnames          = ["FSNS","FLNS","SHFLX","LHFLX"] # Shortwave first! Names of variables in netCDFs
    
# Version (3), Global data (coarsened) from CESM1
elif procname == "cesm1_htr_5degbilinear":
    
    # Dataset and path information
    datname     = procname # Name of dataset (for input and output)
    datpath     = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/proc/"
    outpath     = datpath # Output Path
    vname_out   = "qnet" # 
    
    outname     = "%s%s_%s_Global_1920to2005.nc" % (datpath,datname,vname_out)
    print("File will be saved as %s" % outname)
    
    # NC files (can be anomalized, assume **Positive Upwards**, all units [W/m2])
    nc_fsns         = datpath + "cesm1_htr_5degbilinear_FSNS_Global_1920to2005.nc"  # Incoming Shortwave
    nc_flns         = datpath + "cesm1_htr_5degbilinear_FLNS_Global_1920to2005.nc"  # Outgoing Longwave
    nc_shflx        = datpath + "cesm1_htr_5degbilinear_SHFLX_Global_1920to2005.nc" # Outgoing Sensible Heat Flux
    nc_lhflx        = datpath + "cesm1_htr_5degbilinear_LHFLX_Global_1920to2005.nc" # Outgoing Latent Heat Flux
    nclist          = [nc_fsns,nc_flns,nc_shflx,nc_lhflx]
    vnames          = ["FSNS","FLNS","SHFLX","LHFLX"] # Shortwave first! Names of variables in netCDFs
 
#_Global_1920to2005
    

#%%
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
    
        
    # Set ensemble flag
    lensflag = False
    if 'ens' in list(ds.dims):
        print("Ensemble dimension detected!")
        lensflag = True
        
    
    if vname in fsns_names:
        ds      = ds * -1 # Flip to positive upwards
        print("Flipping sign for Incoming Shortwave")
    
    
    if (vv == 0):
        ds_qnet = ds.copy()
    else:
        ds_qnet = ds_qnet + ds
    
    # Check Append
    if lensflag:
        check.append(ds.sel(lon=lonf,lat=latf,method='nearest').isel(time=1,ens=0).values.item())
    else:
        check.append(ds.sel(lon=lonf,lat=latf,method='nearest').isel(time=1).values.item())
    
    ds.close()
    print("Completed %s in %.2fs" % (vname,time.time()-st))


print("QNET = -FSNS + (FLNS + SHFLX + LHFLX)" )
print("%.2f = %.2f + (%.2f +  %.2f +  %.2f)" % (np.array(check).sum(),check[0],check[1],check[2],check[3]))
if lensflag:
    print("Calculated value = %.2f" % (ds_qnet.sel(lon=lonf,lat=latf,method='nearest').isel(time=1,ens=0).values.item()))
else:
    print("Calculated value = %.2f" % (ds_qnet.sel(lon=lonf,lat=latf,method='nearest').isel(time=1).values.item()))

#%% Save output

st = time.time()
ds_qnet = ds_qnet.rename("qnet")
edict   = {vname_out:dict(zlib=True)}
ds_qnet.to_netcdf(outname,encoding=edict)
print("Save output in %.2fs" % (time.time()-st))
