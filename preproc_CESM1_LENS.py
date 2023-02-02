#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 09:13:18 2022

Preprocesses CESM1-LE for Heat Flux Damping Calculations.
This is RCP85 version (should be written to readjust this)

- Takes TS, LandFrac, IceFrac, and makes SST
- Combines Heat Fluxes
- Prepares Ensemble Average Metric
- Added regridding capabilities for predict_amv project (pred_prep toggle)

Based on the following scripts:
    hf1_enavgrm

@author: gliu
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from tqdm import tqdm
import xarray as xr
import xesmf as xe

#%% User Edits

datpath     = "/vortex/jetstream/climate/data1/yokwon/CESM1_LE/downloaded/atm/proc/tseries/monthly/"
regrid      = 3
regrid_step = True # Set to true if regrid indicates the stepsize rather than total dimension size.
pred_prep   = True # Set to true to output to folders for AMV Prediction...
predpath    = "/stormtrack/data3/glliu/01_Data/04_DeepLearning/CESM_data/LENS_other/"
mask_sep    = True

# Part 1 (Land/Ice Mask Creation)
vnames      = ("LANDFRAC","ICEFRAC") # Variables
mthres      = (0.30,0.05) # Mask out if grid ever exceeds this value

# Part 2 ()
maskmode = "enssum"
mconfig  = 'htr' # ['rcp85','htr']

if mconfig == "rcp85":
    outpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_RCP85/01_PREPROC/"
    mnum    = np.concatenate([np.arange(1,36),np.arange(101,106)])
    ntime = 1140
elif mconfig == "htr":
    outpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_HTR/"
    mnum    = np.concatenate([np.arange(1,36),np.arange(101,108)])
    ntime = 1032
#%% Functions

### NOTE: 2023.01.24 >> Moved these functions onto amv.loader, need to rewrite this to load functions from there....ArithmeticError 
# RCP85 Loader
def load_rcp85(vname,N,datpath=None):
    if datpath is None:
        datpath = "/vortex/jetstream/climate/data1/yokwon/CESM1_LE/downloaded/atm/proc/tseries/monthly/"
        
    # Append variable name to path
    vdatpath = "%s%s/" % (datpath,vname)
        
    # Files are split into 2
    if N<34:
        fn1 = "b.e11.BRCP85C5CNBDRD.f09_g16.%03i.cam.h0.%s.200601-208012.nc" % (N,vname)
        fn2 = "b.e11.BRCP85C5CNBDRD.f09_g16.%03i.cam.h0.%s.208101-210012.nc" % (N,vname)
        ds = []
        for fn in [fn1,fn2]:
            dsf = xr.open_dataset(vdatpath + fn)
            ds.append(dsf)
        ds = xr.concat(ds,dim='time')
    else:
        fn1 = "%sb.e11.BRCP85C5CNBDRD.f09_g16.%03i.cam.h0.%s.200601-210012.nc" % (vdatpath,N,vname)
        ds = xr.open_dataset(fn1)
    return ds[vname]

def load_htr(vname,N,datpath=None):
    if datpath is None:
        datpath = "/vortex/jetstream/climate/data1/yokwon/CESM1_LE/downloaded/atm/proc/tseries/monthly/"
    
    # Append variable name to path
    vdatpath = "%s%s/" % (datpath,vname)
    
    # Ensemble 1 has a different time
    if N == 1:
        fn = "%sb.e11.B20TRC5CNBDRD.f09_g16.%03i.cam.h0.%s.185001-200512.nc" % (vdatpath,N,vname)
    else:
        fn = "%sb.e11.B20TRC5CNBDRD.f09_g16.%03i.cam.h0.%s.192001-200512.nc" % (vdatpath,N,vname)
    ds = xr.open_dataset(fn)
    if N == 1:
        ds = ds.sel(time=slice("1920-02-01","2006-01-01"))
    return ds[vname]
    
    
# ----------------------------
#%% Part 1. Make Land/Ice Mask
# ----------------------------
""" Makes land/ice mask. Saves for all ens members separate and all summed."""
nens = len(mnum) 

# Initialize Mask
mask = [] #np.ones((nens,192,288))
if mask_sep:
    lmasks = []
    imasks = []
# Loop for each ensemble member
for e in tqdm(range(nens)):
    N = mnum[e] # Get ensemble member
    emask = np.ones((192,288)) * np.nan # Preallocate
    if mask_sep:
        imask = emask.copy()
        lmask = emask.copy()
    for v in range(2):
    
        # Load dataset
        vname = vnames[v]
        if mconfig =='rcp85':
            ds = load_rcp85(vname,N,datpath=datpath)
        elif mconfig == 'htr':
            ds = load_htr(vname,N,datpath=datpath)
        invar = ds.values # [Time x Lat x Lon]
        
        # Mask (set 0 where it ever exceeds, propagate in time)
        inthres       = mthres[v]
        if v == 0: # Landmask
            maskpts       = ((invar <= inthres).prod(0)) # [Lat x Lon]
            emask[maskpts==1] = 1 # 1 means it is ocean point
            if mask_sep:
                lmask[maskpts==1] = 1
        elif v == 1:
            maskpts       = ((invar <= inthres).prod(0)) # [Lat x Lon]
            emask[maskpts==0] = np.nan # 0 means it has sea ice
            if mask_sep:
                imask[maskpts==0] = np.nan
    # Save masks separately
    if mask_sep:
        imasks.append(imask)
        lmasks.append(lmask)
    mask.append(emask.copy())
# Make into array
mask = np.array(mask)  # [ENS x LAT x LON]

# Save all members
savename = "%slandice_mask_%s_byens.npy" % (outpath,mconfig)
np.save(savename,mask)

# Save ensemble sum
mask_enssum = mask.prod(0)
savename = "%slandice_mask_%s_ensavg.npy" % (outpath,mconfig)
np.save(savename,mask_enssum)


if mask_sep:
    masklists = [lmasks,imasks]
    masknames = ("land","ice")
    for mm in range(2):
        
        # Save all masks
        maskarr  = np.array(masklists[mm])
        savename = "%s%s_mask_%s_byens_regrid%ideg.npy" % (outpath,masknames[mm],"CESM1",regrid)
        np.save(savename,maskarr)
        
        # Save ens sum
        masks_enssum = maskarr.prod(0)
        savename    = "%s%s_mask_%s_ensavg_regrid%ideg.npy" % (outpath,masknames[mm],"CESM1",regrid)
        np.save(savename,masks_enssum)

# ------------------------------------------------------------
#%% For each variable: Apply LI Mask, Compute Ensemble Average
# ------------------------------------------------------------
usemask = np.load("%slandice_mask_%s_ensavg.npy" % (outpath,mconfig)) # [Lat x Lon]

if pred_prep: # Just prep the surface temperature
    vnames    = ("TS",)
    calc_qnet = False
    savepath  = predpath
    apply_limask =False
    print("Saving for predict_amv in %s" % predpath)
else:
    vnames    = ("TS","FSNS","FLNS","LHFLX","SHFLX",)#"FSNS","FLNS","LHFLX","SHFLX")# ("TS","FSNS","FLNS","LHFLX","SHFLX")
    calc_qnet = False # Set to True to compute Qnet
    apply_limask=True
    savepath  = outpath
nvar      = len(vnames)

if calc_qnet:
    nvar    += 1

for e in tqdm(range(nens)):
    
    # ********************************
    N = mnum[e]
    
    for v,vname in enumerate(vnames):
        
        # Load the data
        if mconfig =='rcp85':
            ds = load_rcp85(vname,N,datpath=datpath)
        elif mconfig == 'htr':
            ds = load_htr(vname,N,datpath=datpath)
        
        # Apply the mask
        if apply_limask:
            ds_msk = ds * usemask[None,:,:]
        else:
            ds_msk = ds
        if vname == "FSNS":
            ds_msk *= -1 # Multiple to upwards positive
        
        # Regrid, if option is set (note, need to see if qnet is sensitive to regridding step placement)
        if regrid:
            
            lat = ds.lat
            lon = ds.lon
            if regrid_step:
                lat_out = np.arange(lat[0],lat[-1]+regrid,regrid)
                lon_out = np.arange(-180,180,regrid)
            else:
                lat_out = np.linspace(lat[0],lat[-1],regrid)
                lon_out = np.linspace(-180,180,regrid)
            
            ds_out    = xr.Dataset({'lat': (['lat'], lat_out), 'lon': (['lon'], lon_out) })
            regridder = xe.Regridder(ds_msk, ds_out, 'bilinear')
            ds_msk    = regridder(ds_msk)
            
        
            
        # Preallocate for first ensemble member!!
        if e == 0:
            nlat,nlon = len(ds_msk.lat),len(ds_msk.lon)
            ensavg = np.zeros((nvar,ntime,nlat,nlon)) # [Var x Time x Lat x Lon]
            if calc_qnet:
                qnet   = np.zeros((ntime,nlat,nlon)) 
        
        # Just Save it for surface temperature
        if (vname == "TS") or (calc_qnet==False):
            
            #print("Saving %s variable separately!" % (vname))
            # Add values for ensemble averaging
            ensavg[v,:,:,:] += ds_msk.values
            
            
            # Save the dataset
            if pred_prep:
                if apply_limask:
                    savename = "%sCESM1_%s_%s_regrid%ideg_ens%02i.nc" % (savepath,mconfig,"ts",regrid,e+1)
                else:
                    savename = "%sCESM1_%s_%s_regrid%ideg_ens%02i_nomask.nc" % (savepath,mconfig,"ts",regrid,e+1)
            else:
                savename = "%sCESM1_%s_%s_ens%02i.nc" % (savepath,mconfig,vname,e+1)
            ds_msk = ds_msk.rename("ts") # Rename TS to ts
            ds_msk.to_netcdf(savename,encoding={"ts": {'zlib': True}})
        else:
            qnet += ds_msk.values

    coords  = {'time':ds_msk.time,'lat':ds_msk.lat,'lon':ds_msk.lon}
    
    # Make/save qnet
    if calc_qnet:
        ensavg[-1,:,:,:] += qnet.copy()
        
        da = xr.DataArray(qnet,
                    dims=coords,
                    coords=coords,
                    name = 'qnet',
                    )
        if pred_prep:
            if apply_limask:
                savename = "%sCESM1_%s_%s_regrid%ideg_ens%02i.nc" % (savepath,mconfig,"qnet",regrid,e+1)
            else:
                savename = "%sCESM1_%s_%s_regrid%ideg_ens%02i_nomask.nc" % (savepath,mconfig,"qnet",regrid,e+1)
        else: # Always apply mask for HFF calculations
            savename = "%sCESM1_%s_%s_ens%02i.nc" % (savepath,mconfig,"qnet",e+1)
        da.to_netcdf(savename,
                 encoding={'qnet': {'zlib': True}})
    
# Compute and save ensemble averages
if calc_qnet:
    vnames = ['ts','qnet']
    #vnames.append('qnet')
#vnames = ['ts','qnet']
    
for v in range(len(vnames)):
    if vnames[v] == "TS":
        vnames_save = "ts"
    else:
        vnames_save = vnames[v]
        
    v_ensavg = ensavg[v,:,:,:]/nens
    
    da = xr.DataArray(v_ensavg,
                dims=coords,
                coords=coords,
                name = vnames_save,
                )

    
    if pred_prep:
        if apply_limask:
            savename = "%sCESM1_%s_%s_regrid%ideg_ensAVG.nc" % (savepath,mconfig,vnames_save,regrid)
        else:
            savename = "%sCESM1_%s_%s_regrid%ideg_ensAVG_nomask.nc" % (savepath,mconfig,vnames_save,regrid)
    
    else:
        savename = "%sCESM1_%s_%s_ensAVG.nc" % (savepath,mconfig,vnames_save)
    da.to_netcdf(savename,
             encoding={vnames_save: {'zlib': True}})

