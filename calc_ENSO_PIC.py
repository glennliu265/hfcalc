#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HF Damping Calculations: Compute ENSO from PIC
Created on Fri Oct 30 12:08:19 2020

b.e11.B1850C5CN.f09_g16.005.cam.h0.TS.040001-049912.nc  b.e11.B1850C5CN.f09_g16.005.cam.h0.TS.100001-109912.nc  b.e11.B1850C5CN.f09_g16.005.cam.h0.TS.160001-169912.nc
b.e11.B1850C5CN.f09_g16.005.cam.h0.TS.050001-059912.nc  b.e11.B1850C5CN.f09_g16.005.cam.h0.TS.110001-119912.nc  b.e11.B1850C5CN.f09_g16.005.cam.h0.TS.170001-179912.nc
b.e11.B1850C5CN.f09_g16.005.cam.h0.TS.060001-069912.nc  b.e11.B1850C5CN.f09_g16.005.cam.h0.TS.120001-129912.nc  b.e11.B1850C5CN.f09_g16.005.cam.h0.TS.180001-189912.nc
b.e11.B1850C5CN.f09_g16.005.cam.h0.TS.070001-079912.nc  b.e11.B1850C5CN.f09_g16.005.cam.h0.TS.130001-139912.nc  b.e11.B1850C5CN.f09_g16.005.cam.h0.TS.190001-199912.nc
b.e11.B1850C5CN.f09_g16.005.cam.h0.TS.080001-089912.nc  b.e11.B1850C5CN.f09_g16.005.cam.h0.TS.140001-149912.nc  b.e11.B1850C5CN.f09_g16.005.cam.h0.TS.200001-209912.nc
b.e11.B1850C5CN.f09_g16.005.cam.h0.TS.090001-099912.nc  b.e11.B1850C5CN.f09_g16.005.cam.h0.TS.150001-159912.nc  b.e11.B1850C5CN.f09_g16.005.cam.h0.TS.210001-220012.nc

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

#%% User Edits

# Mode
mode = 'FULL' # "SLAB or FULL"

# TS [time lat lon]
varkeep = ['TS','time','lat','lon','lev'] 

# PCs to calculate
pcrem = 3# User edited variable!

# Subset data for enso index calculation
bbox = [120, 290, -20, 20]

# Outpath
outpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_PIC_SLAB/ENSO/"
outname = "EOF_ENSO_PIC_%s.npz" % mode

#%% Functions 
def check_ENSO_sign(eofs,pcs,lon,lat,verbose=True):
    """
    checks sign of EOF for ENSO and flips sign by check the sum over
    conventional ENSO boxes (see below for more information)
    
    checks to see if the sum within the region and flips if sum < 0
    
    inputs:
        1) eofs [space (latxlon), PC] EOF spatial pattern from eof_simple
        2) pcs  [time, PC] PC timeseries calculated from eof_simple
        3) lat  [lat] Latitude values
        4) lon  [lon] Longitude values
    
    outputs:
        1) eofs [space, PC]
        2) pcs  [time, PC]
    
    """   
    # Set EOF boxes to check over (west,east,south,north)
    eofbox1 = [190,240,-5,5] #EOF1 sign check Lat/Lon Box (Currently using nino3.4)
    eofbox2 = [190,240,-5,5] #EOF2 (Frankignoul et al. 2011; extend of positive contours offshore S. Am > 0.1)
    eofbox3 = [200,260,5,5] #EOF3 (Frankignoul et al. 2011; narrow band of positive value around equator)
    chkboxes = [eofbox1,eofbox2,eofbox3]
    
    # Find dimensions and separate out space
    nlon = lon.shape[0]
    nlat = lat.shape[0]
    npcs = eofs.shape[1] 
    eofs = eofs.reshape(nlat,nlon,npcs)
    eofs = eofs.transpose(1,0,2) # [lon x lat x npcs]
        
    for n in range(npcs):
        chk = proc.sel_region(eofs,lon,lat,chkboxes[n],reg_sum=1)[n]
        
        if chk < 0:
            if verbose:
                print("Flipping EOF %i because sum is %.2f"%(n,chk))
    
            eofs[:,:,n] *= -1
            pcs[:,n] *= -1
    

    eofs = eofs.transpose(1,0,2).reshape(nlat*nlon,npcs) # Switch back to lat x lon x pcs
    
    return eofs,pcs
    
# Define preprocessing variable
def preprocess(ds,varlist=varkeep):
    """"preprocess dataarray [ds],dropping variables not in [varlist] and 
    selecting surface variables at [lev=-1]"""
    # Drop unwanted dimension
    dsvars = list(ds.variables)
    remvar = [i for i in dsvars if i not in varlist]
    ds = ds.drop(remvar)
    
    # Select the ground level
    ds = ds.isel(lev=-1)
    
    # # Correct first month (Note this isn't working)
    # if ds.time.values[0].month != 1:
    #     startyr = str(ds.time.values[0].year)
    #     endyr = str(ds.time.values[-1].year)
    #     correctedtime = xr.cftime_range(start=startyr,end=endyr,freq="MS",calendar="noleap") 
    #     ds = ds.assign_coords(time=correctedtime) 
    
    return ds

# Get List of nc files for preindustrial control
ncpath = r'/vortex/jetstream/climate/data1/yokwon/CESM1_LE/downloaded/atm/proc/tseries/monthly/TS/'
if mode == 'SLAB':
    ncsearch = 'e.e11.E1850C5CN.f09_g16.001.cam.h0.TS.*.nc'
elif mode == 'FULL':
    ncsearch = 'b.e11.B1850C5CN.f09_g16.005.cam.h0.TS.*.nc'
nclist = glob.glob(ncpath+ncsearch)
nclist.sort()
nclist

# Open dataset
st = time.time()
dsall = xr.open_mfdataset(nclist,concat_dim='time',preprocess=preprocess)
print("Opened in %.2fs"%(time.time()-st))

# # Calculate monthly anomalies (Note, the month files are not grouping properly...)
# st = time.time()
# dsanom = dsall - dsall.groupby('time.month').mean('time')
# print("Monthly Anomalies computed in %.2fs"%(time.time()-st))

# Apply Landice Mask
mask = np.load('/home/glliu/01_Data/00_Scrap/landicemask_enssum.npy')
dsall *= mask[None,:,:]

# Slice to region
dsreg = dsall.sel(lon=slice(bbox[0],bbox[1]),lat=slice(bbox[2],bbox[3]))

# Read out the variables
st = time.time()
ts = dsreg.TS.values
lon = dsreg.lon.values
lat = dsreg.lat.values
times = dsreg.time.values
print("Data loaded in %.2fs"%(time.time()-st))

#%% Calculate ENSO

# Apply Area Weight
_,Y = np.meshgrid(lon,lat)
wgt = np.sqrt(np.cos(np.radians(Y))) # [lat x lon]
ts = ts * wgt[None,:,:]

# Reshape for ENSO calculations
ntime,nlat,nlon = ts.shape 
ts = ts.reshape(ntime,nlat*nlon) # [time x space]
ts = ts.T #[space x time]

# Remove NaN points
okdata,knan,okpts = proc.find_nan(ts,1) # Find Non-Nan Points
oksize = okdata.shape[0]

# Calcuate monthly anomalies
okdata = okdata.reshape(oksize,int(ntime/12),12) # [space x yr x mon]
manom = okdata.mean(1)
tsanom = okdata - manom[:,None,:]
#tsanom = tsanom.reshape(nlat*nlon,ntime)
nyr = tsanom.shape[1]

eofall = np.zeros((nlat*nlon,12,pcrem)) *np.nan# [space x month x pc]
pcall  = np.zeros((nyr,12,pcrem)) *np.nan# [year x month x pc]
varexpall  = np.zeros((12,pcrem)) * np.nan #[month x pc]

# Compute EOF!!
for m in range(12):
    
    # Perform EOF
    st = time.time()
    eofsok,pcs,varexp=proc.eof_simple(tsanom[:,:,m],pcrem,1)
    #print("Performed EOF in %.2fs"%(time.time()-st))

    # Place back into full array
    eofs = np.zeros((nlat*nlon,pcrem)) * np.nan
    eofs[okpts,:] = eofsok   

    # Correct ENSO Signs
    eofs,pcs = check_ENSO_sign(eofs,pcs,lon,lat,verbose=True)
    
    
    # Save variables
    eofall[:,m,:] = eofs.copy()
    pcall[:,m,:] = pcs.copy()
    varexpall[m,:] = varexp.copy()
    
    print("Completed month %i in %.2fs"%(m+1,time.time()-st))

# Replace data back with nan points
# Save Output
st = time.time()
np.savez(outpath+outname,**{
         'eofs': eofall,
         'pcs': pcall,
         'varexp': varexpall,
         'lon': lon,
         'lat':lat,
         'times':times}
        )

print("Data saved in %.2fs"%(time.time()-st))