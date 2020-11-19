#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 23:48:43 2020

Calculate heat flux given

SST and some type of FLX


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

# Paths to use
datpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_PIC_SLAB/02_ENSOREM/" # Output Path
outpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_PIC_SLAB/03_HFCALC/" # Output Path






#%% Functions

def combineflux(fluxname,filepath,downward_positive=True,verbose=False):
    """  
    Combine flux files from CESM to get the following variables...
    (Note expressions below assume downwards_positive, but can set to False to
     flip values). CESM1 Flux values are given as upwards positive.
     1) Net Heat Flux (NHFLX) = FSNS - (LHFLX + SHFLX + FLNS) 
     2) Radiative Heat Flux (RHFLX) = FSNS - FLNS
     3) Turbulent Heat Flux (THFLX) = - (LHFLX + SHFLX)
    
    Parameters
    ----------
    fluxname : STR
        Desired flux. Can be:
            'NHFLX' (Net Heat Flux)
            'RHFLX' (Radiative Heat Flux)
            'THFLX' (Turbulent Heat Flux)
            
    filepaths : STR
        File path and name, with "%s" indicating where
        the flux name is. Note that this is assumed to be
        an npz file.
        
        ex: "/PathtoFile/ENSOREM_%s_HISTORICAL.nc"
    OPTIONAL ---
    downward_positive : BOOL
        Set to TRUE for positive flux into the ocean
        
    verbose : BOOL
        Print added values to check
    
    Returns
    -------
    combineflx : ARRAY
        Numpy Array containing the combined flux

    """

    fluxes = ['FSNS','FLNS','LHFLX','SHFLX']
    
    if fluxname == 'NHFLX':
        inflx = fluxes
    elif fluxname == 'THFLX':
        inflx = fluxes[2:]
    elif fluxname == 'RHFLX':
        inflx = fluxes[:]
    
    i = 0
    for f in inflx:
        
        # Open the file [time x lat x lon]
        flux = np.load(filepath%f,allow_pickle=True)[f]
        
        if f == 'FSNS':
            flux *= -1
        
        # Save value for printing
        savevaladd = flux[0,44,44].copy()
        
        # Append
        if i == 0:
            fluxout = flux.copy()
        else:
            savevalori = fluxout[0,44,44].copy()
            fluxout += flux
            if verbose:
                print("%.2f + %.2f = %.2f"%(savevalori,savevaladd,fluxout[0,44,44].copy()))
        i+= 1
    if downward_positive:
        fluxout *= -1
    return fluxout

def indexwindow(invar,m,monwin,combinetime=False,verbose=False):
    """
    index a specific set of months/years for an odd sliding window
    given the following information (see inputs)
    
    drops the first and last years when a the dec-jan boundary
    is crossed, according to the direction of crossing
    time dimension is thus reduced by 2 overall
    
    inputs:
        1) invar [ARRAY: yr x mon x otherdims] : variable to index
        2) m [int] : index of central month in the window
        3) monwin [int]: total size of moving window of months
        4) combinetime [bool]: set to true to combine mons and years into 1 dimension
    
    output:
        1) varout [ARRAY]
            [yr x mon x otherdims] if combinetime=False
            [time x otherdims] if combinetime=True
    
    """
    
    if monwin > 1:  
        winsize = int(np.floor((monwin-1)/2))
        monid = [m-winsize,m,m+winsize]

    
    varmons = []
    msg = []
    for m in monid:

        if m < 0: # Indexing months from previous year
            
            msg.append("Prev Year")
            varmons.append(invar[:-2,m,:])
            
        elif m > 11: # Indexing months from next year
            msg.append("Next Year")
            varmons.append(invar[2:,m-12,:])
            
        else: # Interior years (drop ends)
            msg.append("Interior Year")
            varmons.append(invar[1:-1,m,:])
    if verbose:
        print("Months are %s with years %s"% (str(monid),str(msg)))       
    # Stack together and combine dims, with time in order
    varout = np.stack(varmons) # [mon x yr x otherdims]
    varout = varout.transpose(1,0,2) # [yr x mon x otherdims]
    if combinetime:
        varout = varout.reshape((varout.shape[0]*varout.shape[1],varout.shape[2])) # combine dims
    return varout
        
    
#%% User inputs

# ENSO Removal Options
ensorem = 1 # Set to 1 if ENSO was removed
ensolag = 1 # Lag between enso removal and variable
pcrem   = 2 # PCs of enso removed
emonwin = 3 # Month window of ENSO Removal

# Flux Calculation options
monwin = 3 # Window of months
lags   = [1,2,3] # Lags to include
flux   = 'NHFLX'


#%%
# Set File Names
ensoexp = "lag%i_pcs%i_monwin%i.npz" % (ensolag,pcrem,emonwin)
filepath = datpath+"ENSOREM_%s_"+ensoexp


# Load the files #[time x lat x lon]
st = time.time()
sst = np.load("%sENSOREM_TS_%s"%(datpath,ensoexp),allow_pickle=True)['TS']
flx = combineflux(flux,filepath,verbose=True)
print("Data loaded in %.2fs"%(time.time()-st))

# Separate out time dimensions, combine space dimensions
ntime,nlat,nlon = sst.shape
sst = sst.reshape(int(ntime/12),12,nlat*nlon)
flx = flx.reshape(sst.shape)

# Calculate heat flux for each lag and each month...
nlag = len(lags)
damping   = np.zeros((12,nlag,nlat*nlon)) # [month, lag, lat, lon]
autocorr  = np.zeros(damping.shape)
crosscorr = np.zeros(damping.shape)

st = time.time()
for l in range(nlag):
    lag = lags[l]
    
    for m in range(12):
        
                
                
        lm = m-lag # Get Lag Month
        
        # Restrict to time ----
        flxmon = indexwindow(flx,m,monwin,combinetime=True,verbose=False)
        sstmon = indexwindow(sst,m,monwin,combinetime=True,verbose=False)
        sstlag = indexwindow(sst,lm,monwin,combinetime=True,verbose=False)
        
        # Compute Correlation Coefficients ----
        crosscorr[m,l,:] = proc.pearsonr_2d(flxmon,sstlag,0) # [space]
        autocorr[m,l,:] = proc.pearsonr_2d(sstmon,sstlag,0) # [space]
        
        # Calculate covariance ----
        cov    = proc.covariance2d(flxmon,sstlag,0)
        autocov = proc.covariance2d(sstmon,sstlag,0)
        
        # Compute damping
        damping[m,l,:] = cov/autocov
        
        print("Completed Month %02i for Lag %s (t = %.2fs)" % (m+1,lag,time.time()-st))
        
        
# Reshape output variables
damping = damping.reshape(12,nlag,nlat,nlon)  
autocorr = autocorr.reshape(damping.shape)
crosscorr = crosscorr.reshape(damping.shape)  
        
# Save output variables
if nlag == 3:
    lagstr = '123'
elif nlag == 2:
    lagstr = "12"
else:
    lagstr = str(lags[0])
    
# Save Damping Parameter
outname1 = "%s%s_Damping_monwin%i_lags%s_ensorem%i_%s.npy" % (outpath,flux,monwin,lagstr,ensorem,ensoexp)
np.save(outname1,damping)

# Save Autocorrelation
outname = "%sSST_Autocorrelation_monwin%i_lags%s_ensorem%i_%s.npy" % (outpath,monwin,lagstr,ensorem,ensoexp)
np.save(outname,autocorr)


# Save Damping Parameter
outname = "%s%s_Crosscorrelation_monwin%i_lags%s_ensorem%i_%s.npy" % (outpath,flux,monwin,lagstr,ensorem,ensoexp)
np.save(outname,crosscorr)

# Sample visualization
import matplotlib.pyplot as plt
plt.pcolormesh(damping[0,0,:,:],vmin=-50,vmax=50,cmap='seismic'),plt.colorbar(),plt.show()  
plt.pcolormesh(autocorr[1,0,:,:],vmin=-1,vmax=1,cmap='seismic'),plt.colorbar(),plt.show()      
    
    

    








