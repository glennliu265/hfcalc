#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 03:59:41 2020

@author: gliu
"""

import xarray as xr
import numpy as np
import glob
import time

import sys
sys.path.append("/home/glliu/00_Scripts/01_Projects/00_Commons/")
sys.path.append("/home/glliu/00_Scripts/01_Projects/01_AMV/02_stochmod/stochmod/model/")
from amv import proc




# ENSO Removal Options
ensorem = 1 # Set to 1 if ENSO was removed
ensolag = 1 # Lag between enso removal and variable
pcrem   = 2 # PCs of enso removed
emonwin = 3 # Month window of ENSO Removal


point  = [-30+360,50] # [lon, lat]
kmonth = 1 
lags = np.arange(0,37,1)

# Paths to use
datpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_PIC_SLAB/02_ENSOREM/" # Output Path
outpath = "/home/glliu/01_Data/00_Scrap/" # Output Path

# Load SST
# Set File Names
ensoexp = "lag%i_pcs%i_monwin%i" % (ensolag,pcrem,emonwin)
filepath = datpath+"ENSOREM_%s_"+ensoexp + ".npz"

# Load the files #[time x lat x lon]
st = time.time()
pz = np.load("%sENSOREM_TS_%s.npz"%(datpath,ensoexp),allow_pickle=True)
TS = pz['TS']
lon = pz['lon']
lat = pz['lat']

# Find values at point
lonf,latf = point
klon,klat = proc.find_latlon(lonf,latf,lon,lat)

# Get value at point
sst = TS[:,klat,klon]

tsmodel = proc.year2mon(sst) # mon x year

# Deseason
tsmodel2 = tsmodel - np.mean(tsmodel,1)[:,None]

# Compute autocorrelation and save data for region
ac = proc.calc_lagcovar(tsmodel2,tsmodel2,lags,kmonth+1,0)

# Save file
locstring = "pointlon%i_lat%i"%(point[0],point[1])
np.save("%sCESM-SLAB_PIC_autocorrelation_%s.npy"%(outpath,locstring),ac)

# # Flux Calculation options
# monwin = 3 # Window of months
# lags   = [1,2,3] # Lags to include
# flux   = 'NHFLX'


#%% Calculate regional averages


# Regional Analysis Settings
bbox_SP = [-60,-15,40,65]
bbox_ST = [-80,-10,20,40]
bbox_TR = [-75,-15,0,20]
bbox_NA = [-80,0 ,0,65]
regions = ("SPG","STG","TRO","NAT")        # Region Names
bboxes = (bbox_SP,bbox_ST,bbox_TR,bbox_NA) # Bounding Boxes


TS = TS.transpose(2,1,0) # Lon Lat Time

# Apply landice mask
limask = np.load(outpath+"landicemask_enssum.npy")
TS *= limask.T[:,:,None]
plt.pcolormesh(lon,lat,TS[:,:,0].T),plt.show()  

# Flip Longitude
lon180,sst = proc.lon360to180(lon,TS)
plt.pcolormesh(lon180,lat,sst[:,:,0].T),plt.show()  


nregions = len(regions)
acs = np.zeros((len(lags),nregions)) # Lag x Region
for r in range(nregions):
    
    sstr,_,_ = proc.sel_region(sst,lon180,lat,bboxes[r])
    tsmodel = np.nanmean(sstr,(0,1))
        
    tsmodel = proc.year2mon(tsmodel) # mon x year
        
    # Deseason
    tsmodel2 = tsmodel - np.mean(tsmodel,1)[:,None]
    
    # Compute autocorrelation and save data for region
    acs[:,r] = proc.calc_lagcovar(tsmodel2,tsmodel2,lags,kmonth+1,0)
np.save("%sCESM-SLAB_PIC_autocorrelation_Regions.npy"%(outpath),acs)


