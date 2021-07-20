#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Attempts to write functions to calculate heat flux
Testing script

Created on Mon Jul 19 16:14:32 2021

@author: gliu
"""

import xarray as xr
import numpy as np
import glob
import time

import sys

import matplotlib.pyplot as plt



#%%
stormtrack = 0
mconfig    = "SLAB_FULL" 



if stormtrack == 1:
    # Module Paths
    sys.path.append("/home/glliu/00_Scripts/01_Projects/00_Commons/")
    sys.path.append("/home/glliu/00_Scripts/01_Projects/01_AMV/02_stochmod/stochmod/model/")
elif stormtrack == 0:
    # Module Paths
    sys.path.append("/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/03_Scripts/stochmod/model/")
    sys.path.append("/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/03_Scripts/")
    
    datpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/"
    outpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/02_Figures/20210722"

    lipath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/01_Data/model_input/"
    llpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/01_Data/model_input/"

from amv import proc
import scm

#mconfig = "SLAB_FULL"
#datpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_PIC_SLAB/02_ENSOREM/%s/" % mconfig
#outpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_PIC_SLAB/03_HFCALC/%s/"  % mconfig


tails = 2
p     = 0.05
mode  = 1

#%% Quick fill-in functions

def load_dampraw(mconfig,datpath):
    inpaths = datpath+"CESM-"+mconfig+"-Damping/"
    damping = np.load(inpaths+"NHFLX_Damping_monwin3_lags123_ensorem1_lag1_pcs2_monwin3.npz.npy")
    rflx = np.load(inpaths+"NHFLX_Crosscorrelation_monwin3_lags123_ensorem1_lag1_pcs2_monwin3.npz.npy")
    rsst = np.load(inpaths+"SST_Autocorrelation_monwin3_lags123_ensorem1_lag1_pcs2_monwin3.npz.npy")
    return damping,rsst,rflx


#%% Test the prep_HF script


# Load some data


# Note this only works locally, need to make equivalent on stormtrack
if stormtrack == 0:
    lon,lat=scm.load_latlon(datpath=llpath,lon360=True)
    lon180,_ = scm.load_latlon(datpath=llpath)
    limask = np.load(lipath+"../landicemask_enssum.npy")



mconfigs =["PIC-SLAB","PIC-FULL"]
dofs     =[898 - 1 - 2 - 2, 1898 - 1 - 2 - 2] 

dampings   = []
rssts      = []
rflxs      = []
dampingfin = []
for i,mcf in enumerate(mconfigs):

    # Load Data
    a,b,c = load_dampraw(mcf,datpath)
    dampings.append(a)
    rssts.append(b)
    rflxs.append(c)
    
    # Apply Masks
    d = scm.prep_HF(a,b,c,p,tails,dofs[i],mode)
    #dampingfin.append(d)
    
    # Postprocess
    dampingw = scm.postprocess_HF(d,limask,[0],lon)
    dampingfin.append(dampingw)

#%% Plot comparison
lonf = -30
latf = 50
klon,klat = proc.find_latlon(lonf,latf,lon180,lat)
klon360,_ = proc.find_latlon(lonf+360,latf,lon,lat)

#%% Load data from CESM-FULL Historical

ds = xr.open_dataset("/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/allens_nhflxdamping_monwin3_sig020_dof082_mode4_lag1.nc")

lbdpt = ds.sel(lon=330,lat=50,method='nearest')
dampfull = lbdpt.NHFLX_Damping.values

#%%



fig,ax = plt.subplots(1,1)
for i in range(2):
    div = 1
    # if i == 1:
    #     div = 2
    ax.plot(dampingfin[i][klon,klat,:]/div,label=mconfigs[i])
    
for i in range(42):
    ax.plot(dampfull[:,i],label="",color='gray',alpha=0.25)
ax.plot(dampfull[:,-1],label="HTR-FULL (member)",color='gray',alpha=0.25)
ax.plot(dampfull.mean(1),label="HTR-FULL (ens-avg)",color='k',alpha=1)

ax.legend()
ax.grid(True,ls='dotted')
ax.set_xlabel("Months")
ax.set_ylabel("$\lambda_a$ (W/m2/degC)")





#Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/Data/CESM-SLAB_FULL-Damping/'
#/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data
#%%









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


def calc_HF(sst,flx,lags,monwin,verbose=True):
    """
    damping,autocorr,crosscorr=calc_HF(sst,flx,lags,monwin,verbose=True)
    Calculates the heat flux damping given SST and FLX anomalies using the
    formula:
        lambda = [SST(t),FLX(t+l)] / [SST(t),SST(t+l)]
    
    
    Inputs
    ------
        1) sst     : ARRAY [year x  x lat x lon] 
            sea surface temperature anomalies
        2) flx     : ARRAY [time x lat x lon]
            heat flux anomalies
        3) lags    : List of INTs
            lags to calculate for (0-N)
        4) monwin  : INT (odd #)
            Moving window of months centered on target month
            (ex. For Jan, monwin=3 is DJF and monwin=1 = J)
        
        --- OPTIONAL ---
        4) verbose : BOOL
            set to true to display print messages
    
    Outputs
    -------     
        1) damping   : ARRAY [month x lag x lat x lon]
            Heat flux damping values
        2) autocorr  : ARRAY [month x lag x lat x lon]
            SST autocorrelation
        3) crosscorr : ARRAY [month x lag x lat x lon]
            SST-FLX cross correlation
    """
    # Reshape variables [time x lat x lon] --> [yr x mon x space]
    ntime,nlat,nlon = sst.shape
    sst = sst.reshape(int(ntime/12),12,nlat*nlon)
    flx = flx.reshape(sst.shape)
    
    # Preallocate
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
            cov     = proc.covariance2d(flxmon,sstlag,0)
            autocov = proc.covariance2d(sstmon,sstlag,0)
            
            # Compute damping
            damping[m,l,:] = cov/autocov
            
            print("Completed Month %02i for Lag %s (t = %.2fs)" % (m+1,lag,time.time()-st))
            
    # Reshape output variables
    damping = damping.reshape(12,nlag,nlat,nlon)  
    autocorr = autocorr.reshape(damping.shape)
    crosscorr = crosscorr.reshape(damping.shape)  
            
    return damping,autocorr,crosscorr

def prep_HF(damping,rsst,rflx,p,tails,dof,mode,returnall=False):
    """
    
    Inputs
    ------
        1) damping   : ARRAY [month x lag x lat x lon]
            Heat flux damping values
        2) autocorr  : ARRAY [month x lag x lat x lon]
            SST autocorrelation
        3) crosscorr : ARRAY [month x lag x lat x lon]
            SST-FLX cross correlation
        4) p : NUMERIC
            p-value
        5) tails : INT
            # of tails for t-test (1 or 2)
        6) dof : INT
            Degrees of freedom
        7) mode: INT
            Apply the following significance testing/masking:
            1 --> No Mask
            2 --> SST autocorrelation based
            3 --> SST-FLX cross correlation based
            4 --> Both 2 and 3
        --- OPTIONAL ---
        8) returnall BOOL
            Set to True to return masks and frequency
    
    Outputs
    -------
        1) dampingmasked [month x lag x lat x lon]
        
    """
    # Determine correlation threshold
    ptilde    = 1-p/tails
    critval   = stats.t.ppf(ptilde,dof)
    corrthres = np.sqrt(1/ ((dof/np.power(critval,2))+1))
    
    # Create Mask
    msst = np.zeros(damping.shape)
    mflx = np.zeros(damping.shape)
    msst[rsst > corrthres] = 1
    mflx[rflx > corrthres] = 1
    if mode == 1:
        mtot = np.ones(damping.shape)     # Total Frequency of successes
        mall = np.copy(mtot)              # Mask that will be applied
        mult = 1
    elif mode == 2:
        mtot = np.copy(msst)
        mall = np.copy(msst)
        mult = 1
    elif mode == 3:
        mtot = np.copy(mflx)
        mall = np.copy(mflx)  
        mult = 1
    elif mode == 4:
        mtot = msst + mflx
        mall = msst * mflx
        mult = 2
    
    # Apply Significance Mask
    dampingmasked = damping * mall
    
    if returnall:
        return dampingmasked,mtot,mult
    return dampingmasked

def postprocess(dampingmasked,limask,sellags,lon):
    
    # Inputs
    ## Dampingmasked [month x lag x lat x lon]
    ## limask [lat x lon]
    
    # Select lags, apply landice mask
    mchoose = dampingmasked[:,sellags,:,:] * limask[None,:,:]
    
    # Flip longiude coordinates ([mon lat lon] --> [lon x lat x mon])
    lon1,dampingw = proc.lon360to180(lon,mchoose.transpose(2,1,0))

    # Multiple by 1 to make positive upwards
    dampingw *= -1
    return dampingw
    
    
    
    
    