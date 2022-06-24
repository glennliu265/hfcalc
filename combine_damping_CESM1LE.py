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


#%% Load in the data

datpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/CESM-HTR-RCP85/"
#%% Load in data (repetitive because last 5 lat values are slightly different)

nclist = glob.glob(datpath+"*.nc")
nclist.sort()
nens = len(nclist)
print(*nclist)

lbd  = []
cc   = []
ac   = []


cv = []
av = []
for i in range(nens):
    ds = xr.open_dataset(nclist[i]) # [Month x Lag x Lat x Lon]
    
    ds.nhflx_damping.isel(month=0,lag=0).plot(),plt.show()
    
    
    lbd.append(ds.nhflx_damping.values)
    cc.append(ds.sst_flx_crosscorr.values)
    ac.append(ds.sst_autocorr.values)
    
    cv.append(ds.cov.values)
    av.append(ds.autocov.values)
    
    

#%% Remake Dataarray

# Prepare variables for input
invars = [lbd,cc,ac,cv,av]
invars = [np.array(v) for v in invars] # [ens x mon x lag x lat x lon]

# Set some attributes (take from calc_enso_general)
varnames = ("nhflx_damping",
            "sst_flx_crosscorr",
            "sst_autocorr",
            "cov",
            "autocov")
varlnames = ("Net Heat Flux Damping",
             "SST-Heat Flux Cross Correlation",
             "SST Autocorrelation",
             "SST-Heat Flux Covariance",
             "SST Autocovariance")
units     = ("W/m2/degC",
             "Correlation",
             "Correlation",
             "W/m2*degC",
             "degC2")
dims  = {"ens":np.arange(1,41,1),
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
savename = "%sCESM1_rcp85_hfdamping_ensorem1_detrend1_2006to2101_allens.nc" % (datpath,)
encoding_dict = {name : {'zlib': True} for name in varnames} 
print("Saving as " + savename)
ds.to_netcdf(savename,
         encoding=encoding_dict)



#%%
