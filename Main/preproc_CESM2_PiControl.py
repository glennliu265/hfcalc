#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Preprocess CESM2 PiControl Data and prepares it for HFF calculations in 
[calc_enso_general]

Inputs:
    datname, datpath, LANDFRAC, ICEFRAC, TS, FLUXES (SHFLX, LHFLX, FSNS, FLNS)
    
Essentially makes:
    ts, qnet (or another flux variable) in [time x lat x lon360]
    land and ice mask (applies to variable)

General steps
- Calculate Net Heat Flux, downward positive
- Apply Land/Ice Mask to everything
- Save output with variable names


# Output (Mask Section)
mask [(ens) x lat x lon360]

based on : [preproc_CESM1_LENS]


Created on Mon Jun 17 09:21:27 2024

@author: gliu
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from tqdm import tqdm
import xarray as xr
import glob

import time
import sys
import scipy as sp

#%% User Edits

datname       = "cesm2_pic"
datpath       = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/proc/" # Path to regridded mixed-layer depth
#bbox_reg      = [-80,0,0,65]
#bbox_reg      = [-100,20,-10,90] # Larger Region 

# Indicate Data PAth
datpath       = "/stormtrack/data4/glliu/01_Data/CESM2_PiControl/FCM/atm/"
mask_sep      = True
outpath       = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/masks/"

# Part 1 (Land/Ice Mask Creation)
vnames        = ("LANDFRAC","ICEFRAC") # Variables
mthres        = (0.30,0.05) # Mask out if grid ever exceeds this value
maskname      = "land%03i_ice%03i" % (mthres[0]*100,mthres[1]*100)

# Make Ice Fraction over specific time
croptime          = True
tstart            =  '0200-01-01' # "2006-01-01" # 
tend              =  '2000-02-01' # "2101-01-01" # 
timestr           =  '%04ito%04i'  % (int(tstart[:4]),int(tend[:4]))# "

mnum          = 1 # Just 1 Ensemble Member for PiControl


# Part 2 (Anomalize + Mask)
anom_varnames = ["TS","SHFLX","LHFLX","FSNS","FLNS"] 
outpath_anom  = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/anom/"



#%% Load Functions for now from amv.proc. but eventually make a helper func [hfutils --> hf]
# Copied stormtrack import from DATA_CODE_CENTRAL.md

# stormtrack
amvpath = "/home/glliu/00_Scripts/01_Projects/00_Commons/" # amv module
scmpath = "/home/glliu/00_Scripts/01_Projects/01_AMV/02_stochmod/stochmod/model/" # scm module

sys.path.append(amvpath)
sys.path.append(scmpath)

#from amv import proc,viz
#import scm
#import amv.loaders as dl

from amv import proc as hf



#%%




def load_cesm2_pic(vname,datpath,searchstr=None):
    # Function to Load CESM2 PiC Output on stormtrack, Searches for 
    # File is assumed to be at <datpath>/<vname>/*<vname>*.nc
    keepvars  = ["time","lat","lon",vname]
    if searchstr is None:
        searchstr = "%s%s/*%s*.nc" % (datpath,vname,vname) # Searches for datpath + *LANDFRAC*.nc
    nclist    = glob.glob(searchstr)
    nclist.sort()
    
    ds_all    = xr.open_mfdataset(nclist,concat_dim="time",combine='nested')
    ds_all    = hf.ds_dropvars(ds_all,keepvars)
    ds_all    = hf.fix_febstart(ds_all)
    
    return ds_all
    


# import hfutils as hf

#%%

# # Part 2 ()
# maskmode = "enssum"
# mconfig  = 'htr' # ['rcp85','htr']

# if mconfig == "rcp85":
#     outpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_RCP85/01_PREPROC/"
#     mnum    = np.concatenate([np.arange(1,36),np.arange(101,106)])
#     ntime = 1140
# elif mconfig == "htr":
#     outpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_HTR/"
#     mnum    = np.concatenate([np.arange(1,36),np.arange(101,108)])
#     ntime = 1032
# #%% Functions

# ### NOTE: 2023.01.24 >> Moved these functions onto amv.loader, need to rewrite this to load functions from there....ArithmeticError 
# # RCP85 Loader
# def load_rcp85(vname,N,datpath=None):
#     if datpath is None:
#         datpath = "/vortex/jetstream/climate/data1/yokwon/CESM1_LE/downloaded/atm/proc/tseries/monthly/"
        
#     # Append variable name to path
#     vdatpath = "%s%s/" % (datpath,vname)
        
#     # Files are split into 2
#     if N<34:
#         fn1 = "b.e11.BRCP85C5CNBDRD.f09_g16.%03i.cam.h0.%s.200601-208012.nc" % (N,vname)
#         fn2 = "b.e11.BRCP85C5CNBDRD.f09_g16.%03i.cam.h0.%s.208101-210012.nc" % (N,vname)
#         ds = []
#         for fn in [fn1,fn2]:
#             dsf = xr.open_dataset(vdatpath + fn)
#             ds.append(dsf)
#         ds = xr.concat(ds,dim='time')
#     else:
#         fn1 = "%sb.e11.BRCP85C5CNBDRD.f09_g16.%03i.cam.h0.%s.200601-210012.nc" % (vdatpath,N,vname)
#         ds = xr.open_dataset(fn1)
#     return ds[vname]

# def load_htr(vname,N,datpath=None):
#     if datpath is None:
#         datpath = "/vortex/jetstream/climate/data1/yokwon/CESM1_LE/downloaded/atm/proc/tseries/monthly/"
    
#     # Append variable name to path
#     vdatpath = "%s%s/" % (datpath,vname)
    
#     # Ensemble 1 has a different time
#     if N == 1:
#         fn = "%sb.e11.B20TRC5CNBDRD.f09_g16.%03i.cam.h0.%s.185001-200512.nc" % (vdatpath,N,vname)
#     else:
#         fn = "%sb.e11.B20TRC5CNBDRD.f09_g16.%03i.cam.h0.%s.192001-200512.nc" % (vdatpath,N,vname)
#     ds = xr.open_dataset(fn)
#     if N == 1:
#         ds = ds.sel(time=slice("1920-02-01","2006-01-01"))
#     return ds[vname]


# def lon180to360(lon180,var,debug=True):
#     """
#     Convert Longitude from Degrees West to Degrees East 
#     Inputs:
#         1. lon180 - array with longitude in degrees west
#         2. var    - corresponding variable [lon x lat x time]
#         3. autoreshape - BOOL, reshape variable autocmatically if size(var) > 3
#     Copied from proc on 2024.02.08
#     """
    
#     # Reshape to combine dimensions
        
#     kw = np.where(lon180 < 0)[0]
#     ke = np.where(lon180 >= 0)[0]
#     lon360 = np.concatenate((lon180[ke],lon180[kw]+360),0)
#     if var is None:
#         return lon360
#     var = np.concatenate((var[ke,...],var[kw,...]),0)
    

#     return lon360,var


def ds_dropvars(ds,keepvars):
    '''Drop variables in ds whose name is not in the list [keepvars]'''
    # Drop unwanted dimension
    dsvars = list(ds.variables)
    remvar = [i for i in dsvars if i not in keepvars]
    ds = ds.drop(remvar)
    return ds

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

def make_limask(landfrac,icefrac,landthres,icethres):
    
    # Get Maximum Percentage over time
    maxicethres  = icefrac.ICEFRAC.max('time') # [Lat x Lon]
    maxlandthres = landfrac.LANDFRAC.max('time') # [Lat x Lon]
    
    # Apply Thresholds
    icemask  = xr.where(maxicethres>icethres,np.nan,1)
    landmask = xr.where(maxlandthres>landthres,np.nan,1) 
    
    return landmask.rename('mask'),icemask.rename('mask')

# ----------------------------
#%% Part 1. Make Land/Ice Mask
# ----------------------------
""" Makes land/ice mask. Saves for all ens members separate and all summed."""

# Loop LANDFRAC, ICEFRAC
maskvar = []
for vv in range(2): 
    
    # Create filename/list and load
    vname     = vnames[vv]
    keepvars  = ["time","lat","lon",vname]
    searchstr = "%s%s/*%s*.nc" % (datpath,vname,vname) # Searches for datpath + *LANDFRAC*.nc
    nclist    = glob.glob(searchstr)
    nclist.sort()
    
    # Drop Unnecessary variables
    ds_all    = xr.open_mfdataset(nclist,concat_dim="time",combine='nested')
    ds_all    = hf.ds_dropvars(ds_all,keepvars)
    ds_all    = hf.fix_febstart(ds_all)
    
    
    # Load the Data
    st = time.time()
    ds_all = ds_all.load()
    print("Loaded in %.2fs" % (time.time()-st))
    
    # Append to files
    maskvar.append(ds_all)


#%% Make the Mask (Should be universal for all models)

st                  = time.time()
# Crop time if optionis set
if croptime:
    print("Cropping time between %s and %s" % (tstart,tend))
    maskvar = [ds.sel(time=slice(tstart,tend)) for ds in maskvar]
landfrac,icefrac    = maskvar
landthres,icethres  = mthres

# Make the mask
landmask,icemask    = make_limask(landfrac,icefrac,landthres,icethres)
print("Computed mask in %.2fs" % (time.time()-st))


limask              = landmask*icemask
outputs             = [landmask,icemask,limask]
outnames            = ["landmask_%02ip" % (landthres*100),
                       'icemask_%02ip' % (icethres*100),
                       'limask_%sp_%sp' % (landthres,icethres)]
edict               = dict(mask=dict(zlib=True))
for vv in range(3):
    outname = "%s%s_%s.nc" % (outpath,datname,outnames[vv])
    if croptime:
        outname = hf.addstrtoext(outname,"_%s" % timestr,adjust=-1)
    outputs[vv].to_netcdf(outname,encoding=edict)

#%% ===========================================================================

# For each variable: Load and apply LI Mask, Anomalize/Detrend, Combine Heat Fluxes
#maskname = '/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/masks/cesm2_pic_limask_0.3p_0.05p_0200to2000.nc'
maskname = None
if maskname is None:
    mask = 1
else:
    mask     = xr.open_dataset(maskname).mask.load()

nvars = len(anom_varnames)
for vv in tqdm(range(nvars)):
    
    # Create filename/list and load
    vname = anom_varnames[vv]
    keepvars  = ["time","lat","lon",vname]
    searchstr = "%s%s/*%s*.nc" % (datpath,vname,vname) # Searches for datpath + *LANDFRAC*.nc
    nclist    = glob.glob(searchstr)
    nclist.sort()
    print("Found %i files for %s" % (len(nclist),vname))
    
    # Drop Unnecessary variables
    ds_all    = xr.open_mfdataset(nclist,concat_dim="time",combine='nested')
    ds_all    = hf.ds_dropvars(ds_all,keepvars)
    ds_all    = hf.fix_febstart(ds_all)
    
    # Load it
    st        = time.time()
    ds_all    = ds_all[vname].load()
    print("Loaded in %.2fs" % (time.time()-st))
    
    # Apply Land Ice Mask
    ds_masked = ds_all * mask
    
    # Remove Seasonal Cycle
    ds_anom   = hf.xrdeseason(ds_masked)
    
    # Remove Trend (Simple Linear)
    
    dt_dict   = hf.detrend_dim(ds_anom.values,0,return_dict=True)# ASSUME TIME in first axis
    
    # Save output
    ds_anom_out = xr.DataArray(dt_dict['detrended_var'],dims=ds_anom.dims,coords=ds_anom.coords,name=vname)
    
    outname   = "%s%s_%s_anom.nc" % (outpath_anom,datname,vname)
    edict     = {vname: dict(zlib=True)}
    ds_anom_out.to_netcdf(outname,encoding=edict)
    
# ------------ 0 ------------ 0 ------------ 0 ------------ 0 ------------ 0 --
#%% Compute Qnet
# ------------ 0 ------------ 0 ------------ 0 ------------ 0 ------------ 0 --

    



# ------------------------------------------------------------------------------
# %% Crop Variables to Selected Region
# ------------------------------------------------------------------------------


bbox_reg = [-100,20,]

# Load in Variables







# #%%


  
    

# # ------------------------------------------------------------
# #%% For each variable: Apply LI Mask, Compute Ensemble Average
# # ------------------------------------------------------------

# # Load  the mask [Lat x Lon]
# savename =  "%smasks/landice_mask_%s_%s_byens.nc" % (outpath,maskname,mconfig)
# ds_mask  = xr.open_dataset(savename)

# usemask  = ds_mask.mask.values.prod(0) # Lat x Lon

# nens      = len(ds_mask.ens)

# if pred_prep: # Just prep the surface temperature
#     vnames    = ("TS",)
#     calc_qnet = False
#     savepath  = predpath
#     apply_limask =False
#     print("Saving for predict_amv in %s" % predpath)
# else:
    
#     if use_SST:
#         vnames    = ("FSNS","FLNS","LHFLX","SHFLX",)#"FSNS","FLNS","LHFLX","SHFLX")# ("TS","FSNS","FLNS","LHFLX","SHFLX")
#     else:
#         vnames    = ("TS","FSNS","FLNS","LHFLX","SHFLX")
#     calc_qnet = True # Set to True to compute Qnet
#     apply_limask=True
#     savepath  = outpath
# nvar      = len(vnames)

# if calc_qnet:
#     nvar    += 1

# for e in tqdm(range(nens)):
    
#     # ********************************
#     N = mnum[e]
    
#     for v,vname in enumerate(vnames):
        
#         # Load the data
#         if mconfig =='rcp85':
#             ds = load_rcp85(vname,N,datpath=datpath)
#         elif mconfig == 'htr':
#             ds = load_htr(vname,N,datpath=datpath)
        
#         ds = fix_febstart(ds).load()
        
#         # Apply the mask
#         if apply_limask:
#             ds_msk = ds * usemask[None,:,:]
#         else:
#             ds_msk = ds
#         if vname == "FSNS":
#             ds_msk *= -1 # Multiple to upwards positive
        
#         # Regrid, if option is set (note, need to see if qnet is sensitive to regridding step placement)
#         if regrid:
            
#             lat = ds.lat
#             lon = ds.lon
#             if regrid_step:
#                 lat_out = np.arange(lat[0],lat[-1]+regrid,regrid)
#                 lon_out = np.arange(-180,180,regrid)
#             else:
#                 lat_out = np.linspace(lat[0],lat[-1],regrid)
#                 lon_out = np.linspace(-180,180,regrid)
            
#             ds_out    = xr.Dataset({'lat': (['lat'], lat_out), 'lon': (['lon'], lon_out) })
#             regridder = xe.Regridder(ds_msk, ds_out, 'bilinear')
#             ds_msk    = regridder(ds_msk)
            
        
            
#         # Preallocate for first ensemble member!!
#         if e == 0:
#             nlat,nlon = len(ds_msk.lat),len(ds_msk.lon)
#             ensavg = np.zeros((nvar,ntime,nlat,nlon)) # [Var x Time x Lat x Lon]
#             if calc_qnet:
#                 qnet   = np.zeros((ntime,nlat,nlon)) 
        
#         # Just Save it for surface temperature
#         if (vname == "TS") or (calc_qnet==False):
            
#             #print("Saving %s variable separately!" % (vname))
#             # Add values for ensemble averaging
#             ensavg[v,:,:,:] += ds_msk.values
            
#             # Save the dataset
#             if pred_prep:
#                 if apply_limask:
#                     savename = "%sCESM1_%s_%s_regrid%ideg_ens%02i.nc" % (savepath,mconfig,"ts",regrid,e+1)
#                 else:
#                     savename = "%sCESM1_%s_%s_regrid%ideg_ens%02i_nomask.nc" % (savepath,mconfig,"ts",regrid,e+1)
#             else:
#                 savename = "%sCESM1_%s_%s_ens%02i.nc" % (savepath,mconfig,vname,e+1)
#             ds_msk = ds_msk.rename("ts") # Rename TS to ts
#             ds_msk.to_netcdf(savename,encoding={"ts": {'zlib': True}})
#         else:
#             qnet += ds_msk.values

#     coords  = {'time':ds_msk.time,'lat':ds_msk.lat,'lon':ds_msk.lon}
    
#     # Make/save qnet
#     if calc_qnet:
#         ensavg[-1,:,:,:] += qnet.copy()
        
#         da = xr.DataArray(qnet,
#                     dims=coords,
#                     coords=coords,
#                     name = 'qnet',
#                     )
#         if pred_prep:
#             if apply_limask:
#                 savename = "%sCESM1_%s_%s_regrid%ideg_ens%02i.nc" % (savepath,mconfig,"qnet",regrid,e+1)
#             else:
#                 savename = "%sCESM1_%s_%s_regrid%ideg_ens%02i_nomask.nc" % (savepath,mconfig,"qnet",regrid,e+1)
#         else: # Always apply mask for HFF calculations
#             savename = "%sCESM1_%s_%s_ens%02i.nc" % (savepath,mconfig,"qnet",e+1)
#         da.to_netcdf(savename,
#                  encoding={'qnet': {'zlib': True}})


# #%% If you were stupid and forgot to set calc_qnet to True (dumbass me)


# # import glob


# # calc_qnet=True # Set it True because you forgot... stupid hack fix to erase later. so tired.

# # for e in tqdm(range(nens)):
# #     for v,vname in enumerate(vnames):
# #         dp = '/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_HTR/'
# #         fn = '%sCESM1_htr_%s_ens%02i.nc' % (dp,vname,e+1)
# #         ds_msk = xr.open_dataset(fn).load() # Signs have already been flipped...
        
# #         # Preallocate for first ensemble member!!
# #         if e == 0:
# #             nlat,nlon = len(ds_msk.lat),len(ds_msk.lon)
# #             ensavgq = np.zeros((1,ntime,nlat,nlon)) # [Var x Time x Lat x Lon]
# #             if calc_qnet:
# #                 qnet   = np.zeros((ntime,nlat,nlon)) 
        
# #         qnet += ds_msk.values
    
# #     coords  = {'time':ds_msk.time,'lat':ds_msk.lat,'lon':ds_msk.lon}
    
# #     # Make/save qnet
# #     ensavgq[0,:,:,:] += qnet.copy()
    
# #     da = xr.DataArray(qnet,
# #                 dims=coords,
# #                 coords=coords,
# #                 name = 'qnet',
# #                 )
# #     if pred_prep:
# #         if apply_limask:
# #             savename = "%sCESM1_%s_%s_regrid%ideg_ens%02i.nc" % (savepath,mconfig,"qnet",regrid,e+1)
# #         else:
# #             savename = "%sCESM1_%s_%s_regrid%ideg_ens%02i_nomask.nc" % (savepath,mconfig,"qnet",regrid,e+1)
# #     else: # Always apply mask for HFF calculations
# #         savename = "%sCESM1_%s_%s_ens%02i.nc" % (savepath,mconfig,"qnet",e+1)
# #     da.to_netcdf(savename,
# #              encoding={'qnet': {'zlib': True}})



# #%% Need to process SST, if it is available

# if use_SST:
#     nclist = [SSTpath + fn for fn in SSTncs]
#     nclist.sort()
#     ds_sst_all = []
    
#     for e in range(nens):
        
#         # Load ds
#         ds = xr.open_dataset(nclist[e]) # [Time, z_t, lat, lon]
#         ds = fix_febstart(ds) # make sure to fix the first month
#         ds = ds.sel(time=slice('1920-01-01','2006-01-01'))
        
#         # Flip Longitude (do this the dumb way)
#         sst = ds.SST.load().values
#         lon180 = ds.lon.values
#         lat = ds.lat.values
#         sst = sst.squeeze().transpose(2,1,0) # Lon x Lat x Time
#         lon360,sst360 = lon180to360(lon180,sst)
#         newcoord = dict(time=ds.time.values,lat=lat,lon=lon360)
#         da_new = xr.DataArray(sst360.transpose(2,1,0),coords=newcoord,dims=newcoord,name="ts")
        
#         edict    = {'ts':{'zlib':True}}
#         savepath = outpath
#         savename = "%sCESM1_%s_%s_ens%02i.nc" % (savepath,mconfig,'ts',e+1)
#         da_new.to_netcdf(savename,encoding=edict)
        
#         ds_sst_all.append(da_new)

#     # Save the ensemble mean
#     ds_concat = xr.concat(ds_sst_all,dim='ens')
    
#     # Compute the ensemble mean
#     ds_ts_ensavg = ds_concat.mean('ens')
#     savename = "%sCESM1_%s_%s_ensAVG.nc" % (savepath,mconfig,'ts')
#     ds_ts_ensavg.to_netcdf(savename,encoding={'ts':{'zlib':True}})
    

# #%% Save the Ensemble average of each variable

# qnet_ensavg = ensavg[-1,:,:,:]/nens
# coords      = dict(time=ds_msk.time,lat=ds_msk.lat,lon=ds_msk.lon)
# da_ensavg_qnet = xr.DataArray(qnet_ensavg,coords=coords,dims=coords,name="qnet")
# savename = "%sCESM1_%s_%s_ensAVG.nc" % (savepath,mconfig,'qnet')
# edict    = {'qnet':{'zlib':True}}
# da_ensavg_qnet.to_netcdf(savename,encoding=edict)

# #%%


# vnames_ensavg = ['FSNS', 'FLNS', 'LHFLX', 'SHFLX',"qnet"]

    
# for v in range(len(vnames_ensavg)):
#     vnames_save = vnames_ensavg[v]
#     eavg_out = ensavg[v,...]
#     coords=dict(())

# #%%

# # Compute and save ensemble averages
# if calc_qnet:
#     vnames = ['ts','qnet']
#     #vnames.append('qnet')
# #vnames = ['ts','qnet']
    
# for v in range(len(vnames)):
#     if vnames[v] == "TS":
#         vnames_save = "ts"
#     else:
#         vnames_save = vnames[v]
    
#     v_ensavg = ensavg[v,:,:,:]/nens
    
#     da = xr.DataArray(v_ensavg,
#                 dims=coords,
#                 coords=coords,
#                 name = vnames_save,
#                 )
    
#     if pred_prep:
#         if apply_limask:
#             savename = "%sCESM1_%s_%s_regrid%ideg_ensAVG.nc" % (savepath,mconfig,vnames_save,regrid)
#         else:
#             savename = "%sCESM1_%s_%s_regrid%ideg_ensAVG_nomask.nc" % (savepath,mconfig,vnames_save,regrid)
    
#     else:
#         savename = "%sCESM1_%s_%s_ensAVG.nc" % (savepath,mconfig,vnames_save)
#     da.to_netcdf(savename,
#              encoding={vnames_save: {'zlib': True}})

# #%% Do more repairs...
# pth    = '/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_HTR'
# nclist = ["CESM1_htr_qnet_ens%02i.nc" % i for i in range(1,43,1)]
# ds_all = [xr.open_dataset(nc) for nc in nclist]

# dsmerge = xr.concat(ds_all,dim='ens')

