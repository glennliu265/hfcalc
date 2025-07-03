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
#%%

device             = "stormtrack"
datpath            = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/"

# CESM2 PiControl, SOM
calcname            = "CESM2_SOM" 
lonname             = 'lon'
latname             = 'lat'
tname               = 'time' 

# Indicate Input Time Crop (for input)
croptime          =   True
tstart            =  '0060-01-01'
tend              =  '0360-12-31'
timestr           =  '%sto%s'  % (tstart[:4],tend[:4]) # ex. 0000to2000

# Indicate HFF calculation crop
croptime_estimate = False # Cut time right before estimating the heat flux feedback
tstart            =  ''
tend              =  ''
tcrop_fname       =  ""
if croptime_estimate:
    tcrop_fname   = "_%sto%s" % (tstart[:4].replace('-',''),tend[:4].replace('-',''))

# Indicate bbox crop information
bbox_name         = "NAtl" # "NAtl

# Variables Information  -----

# SST 
sstname           = "TS"
sstnc             = "OAFLUX_sst_1984to2007.nc"
sstpath           = "/Users/gliu/Globus_File_Transfer/Reanalysis/OAFLUX/"

# Heat Flux
flxname           = "qnet"
flxnc             = "OAFLUX_qnet_1984to2007.nc"
flxpath           = sstpath
varlnames         = "Net Heat Flux Damping"


# Enso Information
ensonc            = "OAFLUX_ENSO_detrend1_pcs3_1984to2007.nc"

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

#%% Main Script END

#++++ Draft Below

 # --------------------------------
 #%% Compute the heat flux feedback
 # --------------------------------
 
 # Load inputs with variables removed
 invars = []
 invars_da = []
 for vv,v in enumerate(vnames_in):
     dataset_name = dataset_names[vv]
     if ensorem:
         savename_anom_ld = "%s%s_%s_%s_%s_detrend%i_ENSOrem_lag%i_pcs%i_monwin%i.nc" % (anompath,dataset_name,v,
                                                                                     bbox_name,timestr,
                                                                                     detrend,ensolag,pcrem,monwin)
     
     else:
         savename_anom_ld = "%s%s_%s_manom_%s_%s_detrend%0i.nc" % (anompath,dataset_name,v,bbox_name,timestr,detrend)
     if lensflag:
         savename_anom_ld = hf.addstrtoext(savename_anom_ld,"_ens%02i"%(ensnum),adjust=-1)
     ds       = xr.open_dataset(savename_anom_ld)
     
     lat = ds.lat.values
     lon = ds.lon.values
     
     # Crop time if option is set
     if croptime_estimate:
         ds = ds.sel(time=slice(tcrop_start,tcrop_end),drop=True)
     
     loadvar         = ds[v].values
     ntime,nlat,nlon = loadvar.shape
     loadvar         = loadvar.reshape(int(ntime/12),12,nlat,nlon)
     
     invars.append(loadvar)
     invars_da.append(ds[v])
 
 #% Calculate heat flux
 sst,flx = invars
 damping,autocorr,crosscorr,autocov,cov = hf.calc_HF(sst,flx,[1,2,3],3,verbose=True,posatm=True,return_cov=True)
 
 # Save heat flux (from hfdamping_mat2nc.py)
 # ----------------------------------------
 outvars  = [damping,crosscorr,autocorr,cov,autocov]
 datpath_out = "%s/%s_damping/" % (hffpath,v)
 hf.makedir(datpath_out)
 savename = "%s%s_hfdamping_%s_%s_ensorem%i_detrend%i.nc" % (datpath_out,dataset_name,
                                                                       bbox_name,timestr,
                                                                       ensorem,detrend)
 if croptime_estimate:
     savename = hf.addstrtoext(savename,"_%scrop" % tcrop_fname,adjust=-1)
     
 if lensflag:
     savename = hf.addstrtoext(savename,"_ens%02i"%(ensnum),adjust=-1)
 dims     = {'month'  :np.arange(1,13,1),
               "lag"  :np.arange(1,4,1),
               "lat"  :lat,
               "lon"  :lon}
 
 # Set some attributes
 varnames = ("%s_damping" % v,
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
#%%                                                                          
    
                                                                                    
#     dsa = ds_anoms[vv]
    
    
# #"%s%s_ENSO_detrend%i_pcs%i_%s.nc" % (ensopath,calcname,detrend,pcrem,timestr)

        
# #%%

        

# # -------------------------------------------------------------------------
# #%% Part 3: Calculate the ENSO Component, and remove it.
# # (based on remove_ENSO_PIC.py)

# # ------------------- -------- General Portion --------------------------------
# """
# Loads enso file containing ensoid
# Loads variable file containing anomalized, masked variable

# Regresses and removes enso component at specified lag, etc
# """
# allstart = time.time()

# # Load ENSO Index
# #dataset_name    = dataset_names[0] # Note, this is the SST
# dataset_name = dataset_names[0]
# savename_ensoid = "%s%s_ENSO_detrend%i_pcs%i_%s.nc" % (ensopath,dataset_name,detrend,pcrem,timestr)
# # if lensflag:
# #     savename_ensoid = hf.addstrtoext(savename_ensoid,"_ens%02i"%(ensnum),adjust=0)
# query   = glob.glob(savename_ensoid)
# if len(query) < 1:
#     print("%s not found.\nPlease run [compute_enso_index.py]" % savename_ensoid)
# ld      = xr.open_dataset(savename_ensoid)#np.load(savename_ensoid,allow_pickle=True)
# ensoid  = ld.pcs#['pcs'] # [year x  month x pc]
# if lensflag:
#     ensoid = ensoid.sel(ens=ensnum)

# # Loop by each variable and remove ENSO
# for vv,v in enumerate(vnames_in):
    
#     # Load Target variable
#     dataset_name = dataset_names[vv]
#     savename_anom = "%s%s_%s_manom_%s_%s_detrend%0i.nc" % (anompath,dataset_name,v,bbox_name,timestr,detrend)
#     if lensflag:
#         savename_anom = hf.addstrtoext(savename_anom,"_ens%02i"%(ensnum),adjust=-1)
#     da = xr.open_dataset(savename_anom)
    
    # Check if (ENSO index file) already exists, and skip if so.
    # ex. cesm2_pic_qnet_NAtl_0200to2000_detrend1_ENSOrem_lag1_pcs3_monwin3.nc
    savename_ensorem = "%s%s_%s_%s_%s_detrend%i_ENSOrem_lag%i_pcs%i_monwin%i.nc" % (anompath,dataset_name,v,
                                                                                bbox_name,timestr,
                                                                                detrend,ensolag,pcrem,monwin)
#     if lensflag:
#         savename_ensorem = hf.addstrtoext(savename_ensorem,"_ens%02i"%(ensnum),adjust=-1)
#     query    = glob.glob(savename_ensorem)
#     if (len(query) < 1) or (overwrite == True):
#         # Read out the variables # [time x lat x lon]
#         st        = time.time()
#         invar     = da[v].values
#         lon       = da[lonname].values
#         lat       = da[latname].values
#         times     = da[tname].values
        
#         # Check if ENSO timeseries size matches the variable
#         ntime_da   = len(times)
#         ntime_enso = ensoid.shape[0] * ensoid.shape[1]
#         if ntime_da != ntime_enso:
#             print("Warning time length of variable (%i) != length of ENSO Index (%i)" % (ntime_da,ntime_enso))
#             print("\tCalculation will fail...")
#             # print("\tCrop will be attempted...")
#             # tstart_da = da.time[0].data
#             # if type(tstart_da) == np.ndarray:
#             #     ystart_da = tstart_da.astype('datetime64[Y]').astype(int) + 1970
                
            
        
#         # Remove ENSO
#         vout,ensopattern,times = hf.remove_enso(invar,ensoid,ensolag,monwin,reduceyr=reduceyr,times=times)
        
#         # Save output variables
#         #savename = "%senso/%s_%s_detrend%i_ENSOrem_lag%i_pcs%i_monwin%i_%s.nc" % (datpath,dataset_name,v,detrend,ensolag,pcrem,monwin,timestr)
#         # if lensflag:
#         #     savename = proc.addstrtoext(savename,"_ens%02i"%(ensnum),adjust=-1)
#         da = hf.numpy_to_da(vout,times,lat,lon,v,savenetcdf=savename_ensorem)
        
#         # Save ENSO component
#         savename_ensocomp = "%s%s_%s_detrend%i_ENSOcmp_lag%i_pcs%i_monwin%i_%s.npz" % (ensopath,dataset_name,v,detrend,ensolag,pcrem,monwin,timestr)
#         if lensflag:
#             savename_ensocomp = hf.addstrtoext(savename_ensocomp,"_ens%02i"%(ensnum),adjust=0)
#         np.savez(savename_ensocomp,**{
#             'ensopattern':ensopattern,
#             'lon':lon,
#             'lat':lat}
#                   )
        
#         print("Completed variable %s (t=%.2fs)" % (v,time.time()-allstart))
#         # End Skip
#     else:
#         print("Skipping. Found existing file: %s" % (str(query)))


    






#%%

st_script = time.time()

for ensnum in np.arange(1,nens+1):
#for ensnum in np.arange(1,nens+1):
    
    # -------------------------------------------------------------------------
    # Part 1: Preprocess Variables (Anomalize, Detrend, Flip latitude if needed)
    # -------------------------------------------------------------------------
    
    """
    
    IN  : ncfile, <dataset_name>_<vname>.nc 
        Contains variable [time x lat x lon360] where li-mask has been applied.
        Also contains lat,lon,time.
    
    OUT : ncfile, <dataset_name>_<vname>_manom_detrend#.nc
        Save as above, but anomalized, detrended and latitude corrected.
        
    ex: ncep_ncar_ts.nc  -->  ncep_ncar_ts_manom_detrend1.nc
    
    """
    
    for vv,v in enumerate(vnames_in):
        
        dataset_name = dataset_names[vv]

            
        # Load the variable processed by [preproc_raw_inputs], [time x lat x lon180]
        ncname = "%s%s_%s_%s_%s.nc" % (procpath,dataset_name,v,bbox_name,timestr)
        da     = xr.open_dataset(ncname) # [Time x Lat x Lon]
        
        # Below section should already be done
        
        
        
        

        if lensflag:
            # Fix ens numbering
            if da.ens.data[0] == 0:
                print("Adjusting ensemble numbering")
                da['ens'] = da.ens.data+1
            
            # Compute Ensemble average for calculations later
            eavg_fname  = hf.addstrtoext(ncname,"_ensavg",adjust=-1)
            query = glob.glob(eavg_fname)
            if (len(query) < 1) or (overwrite == True):
                print("Computing ensemble average for %s..." % v)
                ensavg      = da.mean('ens')
                edict       = hf.make_encoding_dict(ensavg)
                ensavg.to_netcdf(eavg_fname,encoding=edict)
            
            
            
            # Select just 1 ensemble member
            da     = da.sel(ens=ensnum)
        
        
        
        query = glob.glob(savename_anom)
        
        if (len(query) < 1) or (overwrite == True):
            
            # Read out the other variables # [time x lat x lon]
            # -------------------------------------------------
            st    = time.time()
            invar = da[v].values
            lon   = da[lonname].values
            lat   = da[latname].values
            times = da[tname].values
            
            print("Data loaded in %.2fs"%(time.time()-st))
            
            # For LENs case, remove ensavg
            # ----------------------------
            if lensflag and detrend:
                
                if dataset_name in ["rcp85", "htr"]:
                    eavg_fname  = "%sCESM1_%s_%s_ensAVG.nc" % (datpath,dataset_name,v)
                else:
                    eavg_fname  = hf.addstrtoext(ncname,"_ensavg",adjust=-1)
                print("Removing ensemble average, loading from %s" % eavg_fname)
                ensavg      = xr.open_dataset(eavg_fname).load()
                ensavg      = ensavg.sel(time=slice(tstart,tend),drop=True)
                ensavg      = ensavg[v].values
                invar       = invar - ensavg
            else:
                # Make sure it is in time x lat x lon (move this earlier)
                invar = da[v].transpose('time','lat','lon').values
            # elif detrend: # Simple Linear Detrend
            
            #     print("Detrending by removing linear fit")
            #     ds_anom   = da[v].transpose('time','lat','lon')
                
            #     # Simple Linear Detrend
            #     dt_dict   = hf.detrend_dim(ds_anom.values,0,return_dict=True)# ASSUME TIME in first axis
                
                
            #     invar = dt_dict['detrended_var']
            #     # Put back into DataArray
            #     #da = xr.DataArray(dt_dict['detrended_var'],dims=ds_anom.dims,coords=ds_anom.coords,name=vname)

                
                
                
            # Remove monthly anomalies
            # ------------------------
            nmon,nlat,nlon = invar.shape
            manom,invar = hf.calc_clim(invar,0,returnts=1) # Calculate clim with time in axis 0
            vanom = invar - manom[None,:,:,:]
            vanom = vanom.reshape(nmon,nlat,nlon) # Reshape back to [time x lat x lon]
            
            # Flip latitude
            if lat[0] > lat[-1]: # If latitude is decreasing...
                lat   = np.flip(lat)
                vanom = np.flip(vanom,axis=1)
            
            # Detrend the variable (taken from calc_amv_hadisst.py)
            # ----------------------------------------------------
            if lensflag is False:
                start= time.time()
                indata = vanom.reshape(nmon,nlat*nlon).T # Transpose to [Space x Time]
                okdata,knan,okpts = hf.find_nan(indata,1)
                x = np.arange(0,nmon,1)
                if detrend == 0:
                    # Compute global weighted average
                    glomean = hf.area_avg(vanom.transpose(2,1,0),[0,360,-90,90],lon,lat,1)
                    
                    # Regress back to the original data to get the global component
                    beta,b=hf.regress_2d(glomean,okdata,nanwarn=0)
                    
                    # Subtract this from the original data
                    okdt = okdata - beta[:,None]
                else:
                    # Polynomial Detrend
                    okdt,model = hf.detrend_poly(x,okdata,detrend)
                    if debug:
                        fig,ax=plt.subplots(1,1)
                        ax.scatter(x,okdata[44,:],label='raw')
                        ax.plot(x,model[44,:],label='fit')
                        ax.scatter(x,okdt[:,44],label='dt')
                        ax.set_title("Visualize Detrending Method %i"% detrend)
                        #okdt = okdt.T
                        
                data_dt = np.zeros((nmon,nlat*nlon)) * np.nan
                data_dt[:,okpts] = okdt
                data_dt = data_dt.reshape(nmon,nlat,nlon) # Back to [time x lat x lon]
                
                # Save the output to same path...
                da      = hf.numpy_to_da(data_dt,times,lat,lon,v,savenetcdf=savename_anom)
            else:
                # Ensemble average was removed earlier...
                data_dt = vanom
                
                # Save detrended option, if set
                # -----------------------------
                da = hf.numpy_to_da(data_dt,times,lat,lon,v,savenetcdf=savename_anom)
            # End Skip
        else:
            print("Skipping. Found existing file: %s" % (str(query)))
    
    # -------------------------------------------------------------------------
    #%% Part 3: Calculate the ENSO Component, and remove it.
    # (based on remove_ENSO_PIC.py)
    
    # ------------------- -------- General Portion --------------------------------
    """
    Loads enso file containing ensoid
    Loads variable file containing anomalized, masked variable
    
    Regresses and removes enso component at specified lag, etc
    """
    allstart = time.time()
    
    # Load ENSO Index
    #dataset_name    = dataset_names[0] # Note, this is the SST
    dataset_name = dataset_names[0]
    savename_ensoid = "%s%s_ENSO_detrend%i_pcs%i_%s.nc" % (ensopath,dataset_name,detrend,pcrem,timestr)
    # if lensflag:
    #     savename_ensoid = hf.addstrtoext(savename_ensoid,"_ens%02i"%(ensnum),adjust=0)
    query   = glob.glob(savename_ensoid)
    if len(query) < 1:
        print("%s not found.\nPlease run [compute_enso_index.py]" % savename_ensoid)
    ld      = xr.open_dataset(savename_ensoid)#np.load(savename_ensoid,allow_pickle=True)
    ensoid  = ld.pcs#['pcs'] # [year x  month x pc]
    if lensflag:
        ensoid = ensoid.sel(ens=ensnum)
    
    # Loop by each variable and remove ENSO
    for vv,v in enumerate(vnames_in):
        
        # Load Target variable
        dataset_name = dataset_names[vv]
        savename_anom = "%s%s_%s_manom_%s_%s_detrend%0i.nc" % (anompath,dataset_name,v,bbox_name,timestr,detrend)
        if lensflag:
            savename_anom = hf.addstrtoext(savename_anom,"_ens%02i"%(ensnum),adjust=-1)
        da = xr.open_dataset(savename_anom)
        
        # Check if (ENSO index file) already exists, and skip if so.
        # ex. cesm2_pic_qnet_NAtl_0200to2000_detrend1_ENSOrem_lag1_pcs3_monwin3.nc
        savename_ensorem = "%s%s_%s_%s_%s_detrend%i_ENSOrem_lag%i_pcs%i_monwin%i.nc" % (anompath,dataset_name,v,
                                                                                    bbox_name,timestr,
                                                                                    detrend,ensolag,pcrem,monwin)
        if lensflag:
            savename_ensorem = hf.addstrtoext(savename_ensorem,"_ens%02i"%(ensnum),adjust=-1)
        query    = glob.glob(savename_ensorem)
        if (len(query) < 1) or (overwrite == True):
            # Read out the variables # [time x lat x lon]
            st        = time.time()
            invar     = da[v].values
            lon       = da[lonname].values
            lat       = da[latname].values
            times     = da[tname].values
            
            # Check if ENSO timeseries size matches the variable
            ntime_da   = len(times)
            ntime_enso = ensoid.shape[0] * ensoid.shape[1]
            if ntime_da != ntime_enso:
                print("Warning time length of variable (%i) != length of ENSO Index (%i)" % (ntime_da,ntime_enso))
                print("\tCalculation will fail...")
                # print("\tCrop will be attempted...")
                # tstart_da = da.time[0].data
                # if type(tstart_da) == np.ndarray:
                #     ystart_da = tstart_da.astype('datetime64[Y]').astype(int) + 1970
                    
                
            
            # Remove ENSO
            vout,ensopattern,times = hf.remove_enso(invar,ensoid,ensolag,monwin,reduceyr=reduceyr,times=times)
            
            # Save output variables
            #savename = "%senso/%s_%s_detrend%i_ENSOrem_lag%i_pcs%i_monwin%i_%s.nc" % (datpath,dataset_name,v,detrend,ensolag,pcrem,monwin,timestr)
            # if lensflag:
            #     savename = proc.addstrtoext(savename,"_ens%02i"%(ensnum),adjust=-1)
            da = hf.numpy_to_da(vout,times,lat,lon,v,savenetcdf=savename_ensorem)
            
            # Save ENSO component
            savename_ensocomp = "%s%s_%s_detrend%i_ENSOcmp_lag%i_pcs%i_monwin%i_%s.npz" % (ensopath,dataset_name,v,detrend,ensolag,pcrem,monwin,timestr)
            if lensflag:
                savename_ensocomp = hf.addstrtoext(savename_ensocomp,"_ens%02i"%(ensnum),adjust=0)
            np.savez(savename_ensocomp,**{
                'ensopattern':ensopattern,
                'lon':lon,
                'lat':lat}
                      )
            
            print("Completed variable %s (t=%.2fs)" % (v,time.time()-allstart))
            # End Skip
        else:
            print("Skipping. Found existing file: %s" % (str(query)))
    
    # --------------------------------
    #%% Compute the heat flux feedback
    # --------------------------------
    
    # Load inputs with variables removed
    invars = []
    invars_da = []
    for vv,v in enumerate(vnames_in):
        dataset_name = dataset_names[vv]
        if ensorem:
            savename_anom_ld = "%s%s_%s_%s_%s_detrend%i_ENSOrem_lag%i_pcs%i_monwin%i.nc" % (anompath,dataset_name,v,
                                                                                        bbox_name,timestr,
                                                                                        detrend,ensolag,pcrem,monwin)
        
        else:
            savename_anom_ld = "%s%s_%s_manom_%s_%s_detrend%0i.nc" % (anompath,dataset_name,v,bbox_name,timestr,detrend)
        if lensflag:
            savename_anom_ld = hf.addstrtoext(savename_anom_ld,"_ens%02i"%(ensnum),adjust=-1)
        ds       = xr.open_dataset(savename_anom_ld)
        
        lat = ds.lat.values
        lon = ds.lon.values
        
        # Crop time if option is set
        if croptime_estimate:
            ds = ds.sel(time=slice(tcrop_start,tcrop_end),drop=True)
        
        loadvar         = ds[v].values
        ntime,nlat,nlon = loadvar.shape
        loadvar         = loadvar.reshape(int(ntime/12),12,nlat,nlon)
        
        invars.append(loadvar)
        invars_da.append(ds[v])
    
    #% Calculate heat flux
    sst,flx = invars
    damping,autocorr,crosscorr,autocov,cov = hf.calc_HF(sst,flx,[1,2,3],3,verbose=True,posatm=True,return_cov=True)
    
    # Save heat flux (from hfdamping_mat2nc.py)
    # ----------------------------------------
    outvars  = [damping,crosscorr,autocorr,cov,autocov]
    datpath_out = "%s/%s_damping/" % (hffpath,v)
    hf.makedir(datpath_out)
    savename = "%s%s_hfdamping_%s_%s_ensorem%i_detrend%i.nc" % (datpath_out,dataset_name,
                                                                          bbox_name,timestr,
                                                                          ensorem,detrend)
    if croptime_estimate:
        savename = hf.addstrtoext(savename,"_%scrop" % tcrop_fname,adjust=-1)
        
    if lensflag:
        savename = hf.addstrtoext(savename,"_ens%02i"%(ensnum),adjust=-1)
    dims     = {'month'  :np.arange(1,13,1),
                  "lag"  :np.arange(1,4,1),
                  "lat"  :lat,
                  "lon"  :lon}
    
    # Set some attributes
    varnames = ("%s_damping" % v,
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
    print("Saved in %.2fs" % (time.time()-st))
print("Script Ran to Completion in %.2fs"%(time.time()-st_script))
#%%

if lensflag:
    ds_all = []
    for e in range(nens):
        savename = "%s%s_hfdamping_%s_%s_ensorem%i_detrend%i_ens%02i.nc" % (datpath_out,dataset_name,
                                                                              bbox_name,timestr,
                                                                              ensorem,detrend,e+1)
        
        ds = xr.open_dataset(savename).load()
        ds_all.append(ds.copy())
    
    ds_all = xr.concat(ds_all,dim='ens')
    edict  = hf.make_encoding_dict(ds_all)
    
    savename = "%s%s_hfdamping_%s_%s_ensorem%i_detrend%i.nc" % (datpath_out,dataset_name,
                                                                          bbox_name,timestr,
                                                                          ensorem,detrend)
    ds_all.to_netcdf(savename,encoding=edict)
    print("Saved combined output to %s" % savename)
    
    
        



# #%% Delete the Following Probably (from old script)

# #vnames_in            = [sstname,flxname]

# # vnames_in           = ['sst','qnet']#['sst','thflx'] #['TS','LHFLX'] # ["qnet","fsns","flns","lhflx","shflx"] #"TS" for historical data
# # ncs_in              = ['','']
# #dataset_names       = ["OAFLUX","OAFLUX"]#['OISST',"ERA5_RegridOISST"] #'cesm1_htr_5degbilinear'#'cesm2_pic'#'rcp85'

# # Note: Currently hard coded (unfortunately) to read the first dataset name as sst, second as flx
# #hflx_name           = "ERA" # set to None if dataset_name is same as hflx_name

# # # For these datasets, loop for each ensemble member...
# # lens_datasets     = ['htr','rcp85','gfdl_esm2m_lens','csiro_mk36_lens','canesm2_lens']


# # # Determine number of ensemble members
# # dataset_name = dataset_names[0]
# # if dataset_name == 'rcp85':
# #     nens = 40
# # elif dataset_name in ('gfdl_esm2m_lens', "csiro_mk36_lens"):
# #     nens = 30
# # elif dataset_name == 'canesm2_lens':
# #     nens = 50
# # elif dataset_name == 'htr': # CESM1 Historical
# #     nens = 42
# # elif "cesm1_htr" in dataset_name:
# #     nens = 42
# # else:
# #     nens = 1








