#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Combines the output created by calc_enso_general

Created on Mon Jun  6 16:22:19 2022

@author: gliu

"""
import glob
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import sys
import time

#%% Load in the data
st = time.time()

scenario = "htr" # 'rcp85' or 'htr;
debug    = False
vname    = "SHFLX" # [FLNS,FSNS,LHFLX,SHFLX,qnet]

vnames = ("FLNS",)


# Set the datpath depending on the scenario
if scenario == "rcp85":
    timestr = "2006to2100"
    datpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_RCP85/01_PREPROC/"
elif scenario == "htr":
    timestr = "1920to2005"
    datpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_HTR/"


namedict = {
    "LHFLX": "Latent Heat Flux",
    "SHFLX": "Sensible Heat Flux",
    "FSNS" : "Shortwave Flux",
    "FLNS" : "Longwave Flux",
    "qnet" : "Net Heat Flux",
    "NHFLX": "Net Heat Flux",
    "RHFLX": "Radiative Heat Flux",
    "THFLX": "Turbulent Heat Flux"
    }


for vname in vnames:
    
    damping_name = "%s_damping" % vname
    datpath_new = "%s%s/"% (datpath,damping_name) 
    print("Searching for files in %s" % datpath_new)
    
    #%% Load in data (repetitive because last 5 lat values are slightly different)
    
    # Get List of Nc files
    nclist = glob.glob(datpath_new+"*.nc")
    nclist.sort()
    nens   = len(nclist)
    print("Found %i files" % (nens))
    if nens == 0:
        print("Warning! no files found at %s. Exiting." % datpath_new)
        sys.exit()
    if debug:
        print(*nclist,sep="\n")
    
    # Check if combined file already exists
    check_allens = ["allens" in s for s in nclist]
    if np.any(check_allens):
        print('Warning! Already merged files. Check %s. Exiting Script.'%(datpath_new))
        sys.exit()
    
    lbd  = []
    cc   = []
    ac   = []
    cv   = []
    av   = []
    
    for i in range(nens):
        
        ds = xr.open_dataset(nclist[i]) # [Month x Lag x Lat x Lon]
        if debug:
            ds[damping_name].isel(month=0,lag=0).plot(),plt.show()
        
        
        lbd.append(ds[damping_name].values)
        cc.append(ds.sst_flx_crosscorr.values)
        ac.append(ds.sst_autocorr.values)
        
        cv.append(ds.cov.values)
        av.append(ds.autocov.values)
        
        
    
    #%% Remake Dataarray
    
    # Prepare variables for input
    invars = [lbd,cc,ac,cv,av]
    invars = [np.array(v) for v in invars] # [ens x mon x lag x lat x lon]
    
    # Set some attributes (take from calc_enso_general)
    varnames = (damping_name,
                "sst_flx_crosscorr",
                "sst_autocorr",
                "cov",
                "autocov")
    varlnames = (namedict[vname],
                 "SST-Heat Flux Cross Correlation",
                 "SST Autocorrelation",
                 "SST-Heat Flux Covariance",
                 "SST Autocovariance")
    units     = ("W/m2/degC",
                 "Correlation",
                 "Correlation",
                 "W/m2*degC",
                 "degC2")
    dims  = {"ens":np.arange(1,nens+1,1),
             "month":ds.month.values,
             "lag":ds.lag.values,
             "lat":ds.lat.values,
             "lon":ds.lon.values
               }
    
    das = []
    for v,name in enumerate(varnames):
        attr_dict = {'long_name':varlnames[v],
                     'units':units[v]}
        da = xr.DataArray(invars[v],
                    dims=dims,
                    coords=dims,
                    name = name,
                    attrs=attr_dict
                    )
        if v == 0:
            ds = da.to_dataset() # Convert to dataset
        else:
            ds = ds.merge(da) # Merge other datasets
            
        # Append to list if I want to save separate dataarrays
        das.append(ds)
    
    #% Save as netCDF
    # ---------------
    savename = "%sCESM1_%s_%s_ensorem1_detrend1_%s_allens.nc" % (datpath_new,scenario,damping_name,timestr)
    encoding_dict = {name : {'zlib': True} for name in varnames} 
    print("Saving as " + savename)
    ds.to_netcdf(savename,
             encoding=encoding_dict)
    print("Combined %s files in %.2fs"%(vname,time.time()-st))
print("Script ran in %.2fs"%(time.time()-st))



#%%

