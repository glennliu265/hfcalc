#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 15 15:23:17 2020

Remove the ENSO Component of a timeseries

Currently works with CESM1 PIC (SLAB)

@author: gliu
"""

import xarray as xr
import numpy as np
import glob
import time
import sys
#%% User Edits

# Currently only runs on stormtrack
sys.path.append("/home/glliu/00_Scripts/01_Projects/00_Commons/")
from amv import proc,viz

# Paths
datpath   = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_PIC_SLAB/ENSO/" # Path to ENSO Data
varpath   = "/vortex/jetstream/climate/data1/yokwon/CESM1_LE/downloaded/" # Path to raw variables
outpath   = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_PIC_SLAB/02_ENSOREM/" # Output Path

# Name of the ENSO file
mconfig   = "SLAB"
ensofile  = "EOF_ENSO_PIC_%s.npz" % mconfig # Name of the file containing ENSO indices (.npz)
ensoname  = "pcs" # Name of the ENSO Index (PC) variable

# Variables to process
vnames    = ["SHFLX","LHFLX","FSNS","FLNS"]#["TS"]#,"SHFLX","LHFLX","FSNS","FLNS"]

# Removal choices
ensorem   = False # Set to True to remove enso, False to save with enso present (just consolidating files)
pcrem     = 2 # Number of EOFs to remove
ensolag   = 1 # Lag between Variable and ENSO (Variable lags ENSO by _ months)
monwin    = 3 # Moving window of months to consider (3 months centered around 1)
savepat   = True # Set to true to save ENSO patterns
reduceyr  = True # Reduce time-period to account of lags and year crossings 

# if mconfig =="FULL":
#     outpath = outpath + '/FULL/'


#%% Functions

# Function to index time
def indexwindow(invar,m,monwin,combinetime=False):
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
    for m in monid:

        if m < 0: # Indexing months from previous year
            #print("Prev Year")
            varmons.append(invar[:-2,m,:])
            
        elif m > 11: # Indexing months from next year
            #print("Next Year")
            varmons.append(invar[2:,m-12,:])
            
        else: # Interior years (drop ends)
            #print("Int Year")
            varmons.append(invar[1:-1,m,:])
            
    # Stack together and combine dims, with time in order
    varout = np.stack(varmons) # [mon x yr x otherdims]
    varout = varout.transpose(1,0,2) # [yr x mon x otherdims]
    if combinetime:
        varout = varout.reshape((varout.shape[0]*varout.shape[1],varout.shape[2])) # combine dims
    return varout

#%% Load in the data
allstart = time.time()

outpath += "%s/" % mconfig

# Load ENSO
if ensorem:
    ld     = np.load(datpath+ensofile,allow_pickle=True)
    ensoid = ld[ensoname] # [year x  month x pc]
    
    # Standardize the index (stdev = 1)
    ensoid = ensoid/np.std(ensoid,(0,1))

    # Shift ENSO by the specified enso lag
    ensoid = np.roll(ensoid,ensolag,axis=1) # (ex: Idx 0 (Jan) is now 11)
    
# Reduce years if needed
if reduceyr: # Since enso leads, drop years from end of period
    dropyr = int(np.fix(ensolag/12) + 1)
    
    if ensorem: # Drop years for index if removal is part of the step
        ensoid=ensoid[:-dropyr,:,:] 


for v,vname in enumerate(vnames):
    
    # Set experiment name
    if ensorem:
        outname = "ENSOREM_%s_lag%i_pcs%i_monwin%i" %(vname,ensolag,pcrem,monwin)
    else:
        outname = "ENSOREM0_%s_lag%i" % (vname,ensolag)
        
    # Set path
    ncpath = '%satm/proc/tseries/monthly/%s/' % (varpath,vname)
    if mconfig =="SLAB":
        ncsearch = 'e.e11.E1850C5CN.f09_g16.001.cam.h0.%s.*.nc' % vname
    else:
        ncsearch = 'b.e11.B1850C5CN.f09_g16.005.cam.h0.%s.*.nc' % vname
    
    nclist = glob.glob(ncpath+ncsearch)
    nclist.sort()
    
    if len(nclist)<1:
        print("WARNING NO FILES FOUND!")
        break
    
    # Set variables to keep
    varkeep  = [vname,'time','lat','lon','lev'] 
    
    # Define preprocessing function
    def preprocess(ds,varlist=varkeep):
        """"preprocess dataarray [ds],dropping variables not in [varlist] and 
        selecting surface variables at [lev=-1]"""
        # Drop unwanted dimension
        dsvars = list(ds.variables)
        remvar = [i for i in dsvars if i not in varlist]
        ds = ds.drop(remvar)
        
        # Select the ground level
        ds = ds.isel(lev=-1)      
        return ds
    
    # Open dataset # [time x lat x lon]
    st = time.time()
    dsall = xr.open_mfdataset(nclist,concat_dim='time',preprocess=preprocess,combine='nested')
    #print("Opened in %.2fs"%(time.time()-st))
    
    # Read out the variables (bottleneck step here) [time x lat x lon]
    st = time.time()
    invar = dsall[vname].values
    lon = dsall.lon.values
    lat = dsall.lat.values
    times = dsall.time.values
    print("Data loaded in %.2fs"%(time.time()-st))
    
    # Remove monthly anomalies
    manom,invar = proc.calc_clim(invar,0,returnts=1) # Calculate clim with time in axis 0
    vanom = invar - manom[None,:,:,:]
    
    # Reduce years if option is set
    if reduceyr:
        vanom=vanom[dropyr:,:,:,:]
    nyr,nmon,nlat,nlon = vanom.shape
    
    if ensorem:
        # Reshape to combine spatial dimensions
        vanom  = vanom.reshape(nyr,nmon,nlat*nlon) # [year x mon x space]
        vout   = vanom.copy()
        
        # Variable to save enso pattern
        ensopattern = np.zeros((12,nlat,nlon,pcrem))
        
        # Looping by PC...
        for pc in range(pcrem):
            
            # Looping for each month
            for m in range(12):
                
                # Set up indexing
                if monwin > 1:  
                    winsize = int(np.floor((monwin-1)/2))
                    monid = [m-winsize,m,m+winsize]
                    if monid[2] > 11: # Reduce end month if needed
                        monid[2] -= 12
                else:
                    winsize = 0
                    monid = [m]
                    
                if reduceyr:
                    ensoin = indexwindow(ensoid[:,:,[pc]],m,monwin,combinetime=True).squeeze()
                    varin = indexwindow(vanom,m,monwin,combinetime=True)
                    nyr   = int(ensoin.shape[0]/monwin) # Get new year dimension
                else:
                    # Index corresponding timeseries
                    ensoin = ensoid[:,monid,pc] # [yr * mon]  Lagged ENSO
                    varin = vanom[:,monid,:] # [yr * mon * space] Variable
                    
                    # Re-combine time dimensions
                    ensoin = ensoin.reshape(nyr*monwin)
                    varin = varin.reshape(nyr*monwin,nlat*nlon)
                
                # Regress to obtain coefficients [space]
                varreg,_ = proc.regress_2d(ensoin.squeeze(),varin,nanwarn=1)
                
                # Write to enso pattern
                ensopattern[m,:,:,pc] = varreg.reshape(nlat,nlon).copy()
                
                # Expand and multiply out and take mean for period [space,None]*[None,time] t0o [space x time]
                ensocomp = (varreg[:,None] * ensoin[None,:]).squeeze()
                
                # Separate year and mon and take mean along months
                ensocomp = ensocomp.reshape(nlat*nlon,nyr,monwin).mean(2)
                
                # Remove the enso component for the specified month
                vout[winsize:-winsize,m,:] -= ensocomp.T
                
                print("Removed ENSO Component for %s: PC %02i | Month %02i (t=%.2fs)" % (vname,pc+1,m+1,time.time()-allstart))
                # < End Mon Loop>
            # < End PC Loop>
        
        np.savez("%s%s.npz"%(outpath,outname),**{
            vname:vout[winsize:-winsize,:,:].reshape(nyr*12,nlat,nlon),
            'lon':lon,
            'lat':lat}
                 )
        np.savez("%s%s_ensocomponent.npz"%(outpath,outname),**{
            'ensopattern':ensopattern,
            'lon':lon,
            'lat':lat}
                 )
    else: # Otherwise just consolidate the save the variable
        # if monwin > 1:  
        #     winsize = int(np.floor((monwin-1)/2))
        # else:
        #     winsize=0
        np.savez("%s%s.npz"%(outpath,outname),**{
            vname: vanom[:,:,:].reshape(nyr*12,nlat,nlon),
            #vname:vanom[winsize:-winsize,:,:].reshape(nyr*12,nlat,nlon),
            'lon':lon,
            'lat':lat}
                 )
    # End ENSOREM conditional
        
    print("Completed variable %s (t=%.2fs)" % (vname,time.time()-allstart))
