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
regrid      = False
regrid_step = False # Set to true if regrid indicates the stepsize rather than total dimension size.
pred_prep   = False # Set to true to output to folders for AMV Prediction...
predpath    = "/stormtrack/data3/glliu/01_Data/04_DeepLearning/CESM_data/LENS_other/"
mask_sep    = True

# Option to process regridded SST
use_SST     = True # Set to true to use SST regridded by prep_mld script (in stochmod). Otherwise use TS
SSTpath     = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/02_stochmod/SST/"
SSTncs      = ["SST_FULL_HTR_bilinear_num%02i.nc" % i for i in range(42)]

# Part 1 (Land/Ice Mask Creation)
vnames      = ("LANDFRAC","ICEFRAC") # Variables
mthres      = (0.30,1.00) # Mask out if grid ever exceeds this value

maskname  = "land%03i_ice%03i" % (mthres[0]*100,mthres[1]*100)


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


def lon180to360(lon180,var,debug=True):
    """
    Convert Longitude from Degrees West to Degrees East 
    Inputs:
        1. lon180 - array with longitude in degrees west
        2. var    - corresponding variable [lon x lat x time]
        3. autoreshape - BOOL, reshape variable autocmatically if size(var) > 3
    Copied from proc on 2024.02.08
    """
    
    # Reshape to combine dimensions
        
    kw = np.where(lon180 < 0)[0]
    ke = np.where(lon180 >= 0)[0]
    lon360 = np.concatenate((lon180[ke],lon180[kw]+360),0)
    if var is None:
        return lon360
    var = np.concatenate((var[ke,...],var[kw,...]),0)
    

    return lon360,var



def fix_febstart(ds):
    # Copied from preproc_CESM.py on 2022.11.15
    if ds.time.values[0].month != 1:
        print("Warning, first month is %s. Fixing."% ds.time.values[0])
        # Get starting year, must be "YYYY"
        startyr = str(ds.time.values[0].year)
        while len(startyr) < 4:
            startyr = '0' + startyr
        nmon = ds.time.shape[0] # Get number of months
        # Corrected Time
        correctedtime = xr.cftime_range(start=startyr,periods=nmon,freq="MS",calendar="noleap")
        ds = ds.assign_coords(time=correctedtime) 
    return ds

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
        imask = np.ones((192,288))
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
                imask[maskpts==0] = np.nan # All points 1, ice points NaN
    # Save masks separately
    if mask_sep:
        imasks.append(imask)
        lmasks.append(lmask)
    mask.append(emask.copy())
# Make into array
mask = np.array(mask)  # [ENS x LAT x LON]

# Save all members
savename = "%smasks/landice_mask_%s_%s_byens.npy" % (outpath,maskname,mconfig)
np.save(savename,mask)

# Save ensemble sum
mask_enssum = mask.prod(0)
savename = "%smasks/landice_mask_%s_%s_ensavg.npy" % (outpath,maskname,mconfig)
np.save(savename,mask_enssum)

# Save as netCDF
coords  = {'ens':np.arange(1,43,1),'lat':ds.lat.values,'lon':ds.lon.values,}
da_mask = xr.DataArray(mask,coords=coords,dims=coords,name='mask')
savename =  "%smasks/landice_mask_%s_%s_byens.nc" % (outpath,maskname,mconfig)
da_mask.to_netcdf(savename)



if mask_sep:
    
    longlob = ds.lon.values
    latglob = ds.lat.values
    
    masklists = [lmasks,imasks]
    masknames = ("land","ice")
    for mm in range(2):
        
        # Regrid the mask
        if regrid_step:
            
            maskarr = np.array(masklists[mm])
            
            lat_out = np.arange(latglob[0],latglob[-1]+regrid,regrid)
            lon_out = np.arange(-180,180,regrid)
            
            ds_out    = xr.Dataset({'lat': (['lat'], lat_out), 'lon': (['lon'], lon_out) })
            ds_in     = xr.Dataset({'lat': (['lat'], latglob), 'lon': (['lon'], longlob) })
            coordsdict = {"ensemble":np.arange(1,43,1),
                          "lat":latglob,
                          "lon":longlob}
            
            ds_msk    = xr.DataArray(maskarr,
                                     dims=coordsdict,
                                     coords=coordsdict,
                                     name="icemask")# Make into a numpy array
            regridder = xe.Regridder(ds_in, ds_out, 'bilinear')
            ds_msk    = regridder(ds_msk)
        
            masklists[mm] = maskarr
        
        if pred_prep:
            savepath = predpath
        else:
            savepath = outpath
        
        
        # Save all masks
        maskarr  = np.array(masklists[mm])
        
        savename = "%s%s_mask_%s_byens_regrid%ideg.npy" % (savepath,masknames[mm],"CESM1",regrid)
        np.save(savename,maskarr)
        
        # Save ens sum
        masks_enssum = maskarr.prod(0)
        savename    = "%s%s_mask_%s_ensavg_regrid%ideg.npy" % (savepath,masknames[mm],"CESM1",regrid)
        np.save(savename,masks_enssum)

# ------------------------------------------------------------
#%% For each variable: Apply LI Mask, Compute Ensemble Average
# ------------------------------------------------------------

# Load  the mask [Lat x Lon]
savename =  "%smasks/landice_mask_%s_%s_byens.nc" % (outpath,maskname,mconfig)
ds_mask  = xr.open_dataset(savename)

usemask  = ds_mask.mask.values.prod(0) # Lat x Lon

nens      = len(ds_mask.ens)

if pred_prep: # Just prep the surface temperature
    vnames    = ("TS",)
    calc_qnet = False
    savepath  = predpath
    apply_limask =False
    print("Saving for predict_amv in %s" % predpath)
else:
    
    if use_SST:
        vnames    = ("FSNS","FLNS","LHFLX","SHFLX",)#"FSNS","FLNS","LHFLX","SHFLX")# ("TS","FSNS","FLNS","LHFLX","SHFLX")
    else:
        vnames    = ("TS","FSNS","FLNS","LHFLX","SHFLX")
    calc_qnet = True # Set to True to compute Qnet
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
        
        ds = fix_febstart(ds).load()
        
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


#%% If you were stupid and forgot to set calc_qnet to True (dumbass me)


# import glob


# calc_qnet=True # Set it True because you forgot... stupid hack fix to erase later. so tired.

# for e in tqdm(range(nens)):
#     for v,vname in enumerate(vnames):
#         dp = '/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_HTR/'
#         fn = '%sCESM1_htr_%s_ens%02i.nc' % (dp,vname,e+1)
#         ds_msk = xr.open_dataset(fn).load() # Signs have already been flipped...
        
#         # Preallocate for first ensemble member!!
#         if e == 0:
#             nlat,nlon = len(ds_msk.lat),len(ds_msk.lon)
#             ensavgq = np.zeros((1,ntime,nlat,nlon)) # [Var x Time x Lat x Lon]
#             if calc_qnet:
#                 qnet   = np.zeros((ntime,nlat,nlon)) 
        
#         qnet += ds_msk.values
    
#     coords  = {'time':ds_msk.time,'lat':ds_msk.lat,'lon':ds_msk.lon}
    
#     # Make/save qnet
#     ensavgq[0,:,:,:] += qnet.copy()
    
#     da = xr.DataArray(qnet,
#                 dims=coords,
#                 coords=coords,
#                 name = 'qnet',
#                 )
#     if pred_prep:
#         if apply_limask:
#             savename = "%sCESM1_%s_%s_regrid%ideg_ens%02i.nc" % (savepath,mconfig,"qnet",regrid,e+1)
#         else:
#             savename = "%sCESM1_%s_%s_regrid%ideg_ens%02i_nomask.nc" % (savepath,mconfig,"qnet",regrid,e+1)
#     else: # Always apply mask for HFF calculations
#         savename = "%sCESM1_%s_%s_ens%02i.nc" % (savepath,mconfig,"qnet",e+1)
#     da.to_netcdf(savename,
#              encoding={'qnet': {'zlib': True}})



#%% Need to process SST, if it is available

if use_SST:
    nclist = [SSTpath + fn for fn in SSTncs]
    nclist.sort()
    ds_sst_all = []
    
    for e in range(nens):
        
        # Load ds
        ds = xr.open_dataset(nclist[e]) # [Time, z_t, lat, lon]
        ds = fix_febstart(ds) # make sure to fix the first month
        ds = ds.sel(time=slice('1920-01-01','2006-01-01'))
        
        # Flip Longitude (do this the dumb way)
        sst = ds.SST.load().values
        lon180 = ds.lon.values
        lat = ds.lat.values
        sst = sst.squeeze().transpose(2,1,0) # Lon x Lat x Time
        lon360,sst360 = lon180to360(lon180,sst)
        newcoord = dict(time=ds.time.values,lat=lat,lon=lon360)
        da_new = xr.DataArray(sst360.transpose(2,1,0),coords=newcoord,dims=newcoord,name="ts")
        
        edict    = {'ts':{'zlib':True}}
        savepath = outpath
        savename = "%sCESM1_%s_%s_ens%02i.nc" % (savepath,mconfig,'ts',e+1)
        da_new.to_netcdf(savename,encoding=edict)
        
        ds_sst_all.append(da_new)

    # Save the ensemble mean
    ds_concat = xr.concat(ds_sst_all,dim='ens')
    
    # Compute the ensemble mean
    ds_ts_ensavg = ds_concat.mean('ens')
    savename = "%sCESM1_%s_%s_ensAVG.nc" % (savepath,mconfig,'ts')
    ds_ts_ensavg.to_netcdf(savename,encoding={'ts':{'zlib':True}})
    

#%% Save the Ensemble average of each variable

qnet_ensavg = ensavg[-1,:,:,:]/nens
coords      = dict(time=ds_msk.time,lat=ds_msk.lat,lon=ds_msk.lon)
da_ensavg_qnet = xr.DataArray(qnet_ensavg,coords=coords,dims=coords,name="qnet")
savename = "%sCESM1_%s_%s_ensAVG.nc" % (savepath,mconfig,'qnet')
edict    = {'qnet':{'zlib':True}}
da_ensavg_qnet.to_netcdf(savename,encoding=edict)

#%%


vnames_ensavg = ['FSNS', 'FLNS', 'LHFLX', 'SHFLX',"qnet"]

    
for v in range(len(vnames_ensavg)):
    vnames_save = vnames_ensavg[v]
    eavg_out = ensavg[v,...]
    coords=dict(())

#%%

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

#%% Do more repairs...
pth    = '/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_HTR'
nclist = ["CESM1_htr_qnet_ens%02i.nc" % i for i in range(1,43,1)]
ds_all = [xr.open_dataset(nc) for nc in nclist]

dsmerge = xr.concat(ds_all,dim='ens')

CESM1_htr_qnet_ens01.nc