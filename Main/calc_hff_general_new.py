#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate Heat Flux Feedback (General)

Based on calc_enso_general
Need to eventually support lens datasets...
Currently runs on Astraeus. Will modify for more flexibilityin the future

Inputs

    Region/Time-Cropped Flux and TS [time x (ens) x lat x lon180] : from [preproc_raw_inputs], located in datpath + proc/
    ENSO Index (for ENSO removal): from [compute_enso_index], located in datpath + enso/
    
    !! NOTE: Assumes time dimension in ENSO Index and Input Variables are of equal length!!



Update: 2025.07.03: Moved calculation settings to "hfcalc_new_settings_template.py"

Created on Mon Jun 24 14:16:29 2024

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

# =============================================================================
#%% CESM2 PiControl, FOM
# =============================================================================

device             = "stormtrack"
datpath            = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/"

# Dataset Info
calcname            = "CESM2_FOM" 
lonname             = 'lon'
latname             = 'lat'
tname               = 'time' 

# Indicate Input Time Crop (for input)
croptime          =   True
tstart            =  '0200-01-01'
tend              =  '2000-12-31'
timestr           =  '%sto%s'  % (tstart[:4],tend[:4]) # ex. 0000to2000

# Indicate HFF calculation crop
croptime_estimate = False # Cut time right before estimating the heat flux feedback
tcrop_start       =  ''
tcrop_end         =  ''
tcrop_fname       =  ""
if croptime_estimate:
    tcrop_fname   = "_%sto%s" % (tcrop_start[:4].replace('-',''),tcrop_end[:4].replace('-',''))

# Indicate bbox crop information
bbox_name         = "NAtl" # "NAtl

# Variables Information  -----

# SST 
sstname           = "TS"
sstnc             = "CESM2_FOM_TS_NAtl_0200to2000.nc"
sstpath           = "/stormtrack/data4/glliu/01_Data/CESM2_PiControl/proc/NAtl/"

# Heat Flux
flxname           = "SHF"
flxnc             = "CESM2_FOM_SHF_NAtl_0200to2000.nc"
flxpath           = sstpath
varlnames         = "Net Heat Flux Damping"


# Enso Information
ensonc            = "CESM2_FOM_ENSO_detrend1_pcs3_0200to2000.nc"

# Additional Information
lensflag          = False#True
nens              = 1 #42

#%% Import Packages
if device == "Astraeus":

    # local device (currently set to run on Astraeus, customize later)
    amvpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/03_Scripts/" # amv module
    scmpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/03_Scripts/stochmod/model/"

elif device == "stormtrack":
    amvpath = "/home/glliu/00_Scripts/01_Projects/00_Commons/" # amv module
    scmpath = "/home/glliu/00_Scripts/01_Projects/01_AMV/02_stochmod/stochmod/model/" # scm module

sys.path.append(amvpath)
sys.path.append(scmpath)

from amv import proc,viz
import amv.proc as hf
import scm
import amv.loaders as dl


#%% Options for each calculation step

# Step 1 (Preprocessing)
detrend           = 1  # Detrend Method, 0 = remove ens avg, 1 = linear
    
# Step 2 (Remove ENSO)
pcrem             = 3    # PCs calculated
ensolag           = 1    # Lag between ENSO month and response month in NATL
reduceyr          = True # Drop years due to ENSO lag
monwin            = 3    # Window of months to consider

# Step 3 (HFF Calculation)
ensorem           = True # Set to False to skip ENSO removal step

#%% Additional Toggles

overwrite         = True # Set to True to ovewrite
debug             = True # Set to true for debugging plots


#%% Set up paths

# Set Paths
anompath = datpath + "anom/"
ensopath = datpath + "enso/"
hffpath  = datpath + "hff/"
maskpath = datpath + "masks/"
procpath = datpath + "proc/"

hf.makedir(anompath)
hf.makedir(ensopath)
hf.makedir(hffpath)
hf.makedir(maskpath)
hf.makedir(procpath)



#%% Main Script Start

st_script = time.time()

# Set up Calculations
vnames_in = [sstname,flxname]
ncs_in    = [sstpath+sstnc,flxpath+flxnc]

# # Load the data
# ds_all    = [xr.open_dataset(ncs_in[vv])[vnames_in[vv]].load() for vv in range(2)]

# # Format the data array
# ds_all    = [hf.format_ds(ds,latname=latname,lonname=lonname,timename=tname) for ds in ds_all]

# -------------------------------------------------------------------------
#%% Part 1: Preprocess Variables (Anomalize, Detrend, Flip latitude if needed)
# -------------------------------------------------------------------------

# in: sst, flux with [lonname,latname,timename]
# out" sstanome: [time x lat x lon]
print("(Step 1.) Calculating Anomalies...")
def preprocess_ds(ds,tstart,tend,detrend):
    # Fix February Start (for CESM1 Output)
    try:
        ds = hf.fix_febstart(ds)
    except:
        print("\tWarning: Skipping Febstart check due to issue with fix_febstart")
    
    # Slice to time period of interest
    if croptime:
        ds = ds.sel(time=slice(tstart,tend),drop=True)
        
    # Deseason ----------------
    dsa   = hf.xrdeseason(ds)
    
    # Detrend ----------------
    if detrend == 0:
        dsadt = dsa - dsa.mean('ens')
    else:
        dsadt = hf.xrdetrend(dsa,verbose=False)
    dsadt = dsadt.transpose('time','lat','lon')
    return dsadt 


def check_exist(fn,overwrite):
    query         = glob.glob(savename_anom)
    if (len(query) > 0) and (overwrite is False):
        return True
    return False

# Preprocess SST and Flux
ds_anoms = []
for vv in range(2):
    
    vname         = vnames_in[vv]
    savename_anom = "%s%s_%s_manom_%s_%s_detrend%0i.nc" % (anompath,calcname,vname,bbox_name,timestr,detrend)
    # Check to see if it exists
    check = check_exist(savename_anom,overwrite)
    if check:
        print("\tSkipping anomalize step... File found: %s"  % savename_anom)
        dsa = xr.open_dataset(savename_anom)[vnames_in[vv]].load()
        
    else:
        
        # Load Dataset and format
        ds_in = xr.open_dataset(ncs_in[vv])[vnames_in[vv]].load()
        ds_in = hf.format_ds(ds_in,latname=latname,lonname=lonname,timename=tname)
        
        
        
        dsa   = preprocess_ds(ds_in,tstart,tend,detrend)
        edict = {vname:{'zlib':True}}
        dsa.to_netcdf(savename_anom,encoding=edict)
        print("\t\tSaved File to %s..." % savename_anom)
    
    ds_anoms.append(dsa)
print("\n\n")

#%% Remove ENSO
print("----------------------------------------------------------------------")
print("(Step 2.) Removing ENSO...")
# Load ENSO Index
savename_ensoid = ensopath + ensonc
query   = glob.glob(savename_ensoid)
if len(query) < 1:
    print("%s not found.\nPlease run [compute_enso_index.py]" % savename_ensoid)
    exit
ds_enso = xr.open_dataset(savename_ensoid).load()
ensoid  = ds_enso.pcs

# Loop by each variable and remove ENSO
ds_ensorems = []
for vv in range(2):
    vname            = vnames_in[vv]
    
    # Check if (ENSO index file) already exists, and skip if so.
    # ex. cesm2_pic_qnet_NAtl_0200to2000_detrend1_ENSOrem_lag1_pcs3_monwin3.nc
    savename_ensorem = "%s%s_%s_%s_%s_detrend%i_ENSOrem_lag%i_pcs%i_monwin%i.nc" % (anompath,calcname,vname,
                                                                                bbox_name,timestr,
                                                                                detrend,ensolag,pcrem,monwin)
    
    check = check_exist(savename_ensorem,overwrite)
    if check:
        print("\tSkipping ENSO removal step... File found: %s"  % savename_ensorem)
        da_ensorem = xr.open_dataset(savename_ensorem)[vnames_in[vv]].load()
    else:
        da = ds_anoms[vv]
        
        # Calculate ENSO
        st        = time.time()
        invar     = da.data
        lon       = da.lon.data
        lat       = da.lat.data
        times     = da.time.data
                
        # Check if ENSO timeseries size matches the variable
        ntime_da   = len(times)
        ntime_enso = len(ensoid.year.data) * 12
        if ntime_da != ntime_enso:
            print("ERROR: time length of variable (%i) != length of ENSO Index (%i)" % (ntime_da,ntime_enso))
            print("\tExiting calculation...")
            exit
        
                
        # Remove ENSO
        vout,ensopattern,times = hf.remove_enso(invar,ensoid,ensolag,monwin,reduceyr=reduceyr,times=times)
        

        # Save output variables
        da_ensorem = hf.numpy_to_da(vout,times,lat,lon,vname,savenetcdf=savename_ensorem)
        
        # Save ENSO component
        savename_ensocomp = "%s%s_%s_detrend%i_ENSOcmp_lag%i_pcs%i_monwin%i_%s.nc" % (ensopath,calcname,vname,detrend,ensolag,pcrem,monwin,timestr)
        coords_ensocomp   = dict(mon=np.arange(1,13,1),lat=lat,lon=lon,pc=np.arange(1,pcrem+1))
        da_ensocomp = xr.DataArray(ensopattern,coords=coords_ensocomp,dims=coords_ensocomp,name=vname)
        edict = {vname:{'zlib':True}}
        da_ensocomp.to_netcdf(savename_ensocomp,encoding=edict)
    
    ds_ensorems.append(da_ensorem)
    print("Completed variable %s (t=%.2fs)" % (vname,time.time()-st))

print("\n\n")

#%% Calculate Net Heat Flux Feedback
print("----------------------------------------------------------------------")

if ensorem:
    dsanoms_in = ds_ensorems
else:
    dsanoms_in = ds_anoms

invars    = []
invars_da = []
for vv in range(2):
    
    vname = vnames_in[vv]
    ds    = dsanoms_in[vv]
    
    
    # Crop time if option is set
    if croptime_estimate:
        ds = ds.sel(time=slice(tcrop_start,tcrop_end),drop=True)
        
        
    # Get Dimensions
    loadvar         = ds.data
    lat             = ds.lat.data
    lon             = ds.lon.data
    times            = ds.time.data
    ntime,nlat,nlon = loadvar.shape
    loadvar         = loadvar.reshape(int(ntime/12),12,nlat,nlon)
    
    invars.append(loadvar)
    invars_da.append(ds)
    
#% Calculate heat flux
sst,flx = invars
damping,autocorr,crosscorr,autocov,cov = hf.calc_HF(sst,flx,[1,2,3],3,verbose=True,posatm=True,return_cov=True)
 
 
# Save heat flux (from hfdamping_mat2nc.py)
# ----------------------------------------
outvars  = [damping,crosscorr,autocorr,cov,autocov]
savename = "%s%s_%s_damping_%s_%s_ensorem%i_detrend%i.nc" % (hffpath,calcname,flxname,
                                                                       bbox_name,timestr,
                                                                       ensorem,detrend)
if croptime_estimate:
    savename = hf.addstrtoext(savename,"_%scrop" % tcrop_fname,adjust=-1)
                 

# Save Output
dims     = {'month'  :np.arange(1,13,1),
              "lag"  :np.arange(1,4,1),
              "lat"  :lat,
              "lon"  :lon}

# Set some attributes
varnames = ("%s_damping" % flxname,
            "sst_flx_crosscorr",
            "sst_autocorr",
            "cov",
            "autocov")
varlnames = ("Net Heat Flux Damping",
              "SST-Heat Flux Cross Correlation",
              "SST Autocorrelation",
              "SST-Heat Flux Covariance",
              "SST Autocovariance")
units     = ("W/m2/degC",
              "Correlation",
              "Correlation",
              "W/m2*degC",
              "degC2")

das = []
for v,name in enumerate(varnames):

    attr_dict = {'long_name':varlnames[v],
                  'units':units[v]}
    da = xr.DataArray(outvars[v],
                dims=dims,
                coords=dims,
                name = name,
                attrs=attr_dict
                )
    if v == 0:
        ds = da.to_dataset() # Convert to dataset
    else:
        ds = ds.merge(da) # Merge other datasets
        
    # Append to list if I want to save separate dataarrays
    das.append(ds)

#% Save as netCDF
# ---------------
st = time.time()
encoding_dict = {name : {'zlib': True} for name in varnames} 
print("Saving as " + savename)
ds.to_netcdf(savename,
          encoding=encoding_dict)








