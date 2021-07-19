#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Attempts to write functions to calculate heat flux

Created on Mon Jul 19 16:14:32 2021

@author: gliu
"""
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
        1) sst     : ARRAY [time x lat x lon] 
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
    
    # Apply Mask
    dampingmasked = damping * mall
    
    if returnall:
        return dampingmasked,mtot,mult
    return dampingmasked
    
    