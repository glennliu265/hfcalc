#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Preprocess variables for heat flux damping calculations

Adapted from
- preproc_CESM1_LENS.py

Created on Wed Jun 15 15:41:48 2022

@author: gliu
"""

from tqdm import tqdm
import glob

import numpy as np
import xarray as xr
import xesmf as xe
import matplotlib.pyplot as plt

#%% User Edits

modelname = "csiro_mk36_lens"


#"gfdl_esm2m_lens"

datpath   = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/00_Commons/CLIVAR_LE/%s/Amon/" % modelname
outpath   = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_lens/%s/" % modelname

#mnum     = np.concatenate([np.arange(1,36),np.arange(101,106)])

# Part 1 (Land/Ice Mask Creation)
landpath   = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/00_Commons/CLIVAR_LE/%s/fx/sftlf/" % modelname
icepath    = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/00_Commons/CLIVAR_LE/%s/OImon/sic/" % modelname
method     = "bilinear" # regridding method

maskpaths  = (landpath,icepath)
mvnames     = ("sftlf","sic") # Variables
mthres     = (0.30,0.05) # Mask out if grid ever exceeds this value

# Part 2 ()
maskmode   = "enssum"


vname_cmip = ("rsus" ,"rlus" ,"rsds" ,"rlds" ,"hfss" ,"hfls","ts")
vname_new  = ("fsns" ,"flns" ,"fsns" ,"flns" ,"hfss","hfls","ts")


def get_lens_nc(modelname,vname,e,compname="Amon"):
    if modelname == "gfdl_esm2m_lens":
        ncname = "%s_%s_GFDL-ESM2M_historical_rcp85_r%ii1p1_195001-210012.nc" % (vname,compname,e+1)
    elif modelname == "csiro_mk36_lens":
        ncname = "%s_%s_CSIRO-Mk3-6-0_historical_rcp85_r%ii1p1_185001-210012.nc" % (vname,compname,e+1)
    return ncname
    


#%% Get list of files to see total ensemble count

ncsearch = "%s/%s/*.nc" % (datpath,vname_cmip[0])
nclist   = glob.glob(ncsearch)
nclist.sort()
nens = len(nclist)
print(len(nclist))




#%% Functions


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

# ----------------------------
#%% Part 1. Make Land/Ice Mask
# ----------------------------
""" Makes land/ice mask. Saves for all ens members separate and all summed."""
# Get number of ensemble members


# Regrid Ice values




# Initialize Mask
mask = [] #np.ones((nens,192,288))

for e in tqdm(range(nens)):
    
    # Get netCDF names
    icenc  = get_lens_nc(modelname,"sic",e,compname="OImon")
    
    if modelname == "gfdl_esm2m_lens":
        landnc = "sftlf_fx_GFDL-ESM2M_historical_r0i0p0.nc"
    elif modelname == "csiro_mk36_lens":
        
        landpath = outpath 
        landnc   = "barot_mask_csiro_mk36_lens.nc"
        
    else:
        landnc = landnc = get_lens_nc(modelname,"sftlf",e,compname="fx")
    
    # Open dataset
    dsice  = xr.open_dataset(icepath+icenc)
    dsland = xr.open_dataset(landpath+landnc)
    if e == 0: # Get Lat Lon
    
        if modelname == "csiro_mk36_lens": # Lat/Lon in other dataset
            ncpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/00_Commons/CLIVAR_LE/csiro_mk36_lens/Amon/ts/"
            ncname = "ts_Amon_CSIRO-Mk3-6-0_historical_rcp85_r10i1p1_185001-210012.nc"
            dstemp = xr.open_dataset(ncpath+ncname)
            lat = dstemp.lat.values
            lon = dstemp.lon.values
            
        else:
            lat   = dsland.lat.values
            lon   = dsland.lon.values
        
    
    # Regrid ice/land values (based on prep_mld_PIC.py)
    ds_out    = xr.Dataset({'lat':lat,'lon':lon}) # Define new grid (note: seems to support flipping longitude)
    regridder = xe.Regridder(dsice,ds_out,method,periodic=True)
    daproc    = regridder(dsice[mvnames[1]]) # Need to input dataarray
    # Get land and ice fractions
    if modelname == "csiro_mk36_lens": # Regrid landmask as well
        regridder_land = xe.Regridder(dsland,ds_out,method,periodic=True)
        daproc_land    = regridder_land(dsland['landmask'])
        landfrac = daproc_land.values
    else:
        landfrac = dsland.sftlf.values
    icefrac  = daproc.values
    
    if np.any(icefrac>1):
        icefrac /=100 # Convert from % to decimal
    
    # Preallocate
    emask = np.ones((len(lat),len(lon))) * np.nan 
    
    # Make Landmask
    invar   = landfrac
    inthres = mthres[0] 
    if modelname == "csiro_mk36_lens":
        maskpts = ~np.isnan(invar)
    else:
        maskpts = ((invar <= inthres)) # 0 is Land
    emask[maskpts==1] = 1 # 1 means it is ocean point
    
    # Make Ice Mask
    invar   = icefrac
    inthres = mthres[1]
    maskpts       = ((invar <= inthres).prod(0)) # 0 means it has sea ice
    emask[maskpts==0] = np.nan
    
    
    mask.append(emask.copy())

# Make into array
mask = np.array(mask)  # [ENS x LAT x LON]

# Save all members
savename    = "%slandice_mask_%s_byens.npy" % (outpath,modelname)
np.save(savename,mask)

# Save ensemble sum
mask_enssum = mask.prod(0)
savename    = "%slandice_mask_%s_ensavg.npy" % (outpath,modelname)
np.save(savename,mask_enssum)

# Get Time dimension for later
ntime = len(daproc.time)
nlat  = len(lat)
nlon  = len(lon)

# ------------------------------------------------------------
#%% For each variable: Apply LI Mask, Compute Ensemble Average
# ------------------------------------------------------------

usemask = np.load("%slandice_mask_%s_ensavg.npy" % (outpath,modelname)) # [Lat x Lon]
nvar    = len(vname_cmip)

for e in tqdm(range(nens)):
    
    ensavg = np.zeros((2,ntime,nlat,nlon)) # [Var x Time x Lat x Lon]
    qnet   = np.zeros((ntime,nlat,nlon)) 
    
    # Loop for each variable
    for v in range(nvar):
        
        # Read in variable
        vname_in = vname_cmip[v]
        ncname   = get_lens_nc(modelname,vname_in,e,)
        ds       = xr.open_dataset(datpath+vname_in+"/"+ncname)
        invar    = ds[vname_in]
        
        # Apply the mask
        invar    *= usemask[None,:,:]
        
        # Just save it
        if vname_in == "ts":
            
            # Add values for ensemble averaging
            ensavg[0,:,:,:] += invar.values
            
            # Save the dataset
            savename = "%s%s_%s_ens%02i.nc" % (outpath,modelname,"ts",e+1)
            ds_msk = invar.rename('ts')
            ds_msk.to_netcdf(savename,encoding={'ts': {'zlib': True}})
            
        # Compute heat fluxes
        else:
            
            # Flip so that it is upwards positive
            if "d" in vname_in:
                invar *= -1 # Multiply such that downwards positive
            
            # # Both latent and sensible are upwards positive
            # if ("u" in vname_in) or vname in ["hfls","hfss"]:
            #     invar *= -1 # Multiply such that downwards positive
            # Add to qnet
            qnet += invar.values
        # End variable loop
    
    # Make/save qnet
    ensavg[1,:,:,:] += qnet.copy()
    coords  = {'time':invar.time,'lat':invar.lat,'lon':invar.lon}
    da = xr.DataArray(qnet,
                dims=coords,
                coords=coords,
                name = 'qnet',
                )
    savename = "%s%s_%s_ens%02i.nc" % (outpath,modelname,"qnet",e+1)
    da.to_netcdf(savename,
             encoding={'qnet': {'zlib': True}})
    
    
# Compute and save ensemble averages
vnames = ['ts','qnet']

for v in range(2):
    v_ensavg = ensavg[v,:,:,:]/nens
    
    da = xr.DataArray(v_ensavg,
                dims=coords,
                coords=coords,
                name = vnames[v],
                )
    savename = "%s%s_%s_ensAVG.nc" % (outpath,modelname,vnames[v])
    da.to_netcdf(savename,
             encoding={vnames[v]: {'zlib': True}})
