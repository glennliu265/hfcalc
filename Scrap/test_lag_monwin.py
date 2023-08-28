#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 01:13:02 2020

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


# Test
mons = np.array(["DJF","JFM","FMA","MAM",
        "AMJ","MJJ","JJA","JAS",
        "ASO","SON","OND","NDJ"])

    
# [year mon otherdims]
arr1 = np.tile(mons[:,None],100)[:,:,None].transpose(1,0,2)
arr2 = arr1.copy()


lag      = 3
lagcross = 0 # Add 1 when lag begins to cross into the new year
monwin = 3
for m in range(12):
    lm = m-3
    basearr = indexwindow(arr1,m,monwin,combinetime=False,verbose=True)
    lagarr  = indexwindow(arr2,lm,monwin,combinetime=False,verbose=True)
    
    print("Base Arr is %s"%str(basearr.shape))
    print("Lag Arr is %s"%(str(lagarr.shape)))
    
