#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Compute the ENSO index for ENSO removal...
Decided to separate this from calc_enso_general
to reduce the size of datasets, particularly for hi-res output

Main objective is to just get the INDEX...
Regression and removal will occur in the actual preprocessing script...

Script Sections
(1) Load data, apply land ice mask
(2) Detrend and Deseason
(3) Compute the ENSO indices


Output (netcdf form):
    
    eofs    : [lat x lon x month x pc]  - EOF Patterns
    pcs     : [year x mon x pc]         - Principle Component Timeseries
    varexp  : [mon x pc]                - Variance Explained


Copied upper section of [calc_enso_general]

Created on Tue Jun 18 08:39:34 2024

@author: gliu
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import time
import sys
import cartopy.crs as ccrs
import glob

import pandas as pd


#%% Import modules
stormtrack = 1
if stormtrack:
    sys.path.append("/home/glliu/00_Scripts/01_Projects/00_Commons/")
    sys.path.append("/home/glliu/00_Scripts/01_Projects/01_AMV/02_stochmod/stochmod/model/")
    
    # Path to the processed dataset (qnet and ts fields, full, time x lat x lon)
    #datpath =  "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_RCP85/01_PREPROC/"
    datpath =  "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/anom/"
    figpath =  "/home/glliu/02_Figures/01_WeeklyMeetings/20240621/"
    
else:
    sys.path.append("/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/03_Scripts/stochmod/model/")
    sys.path.append("/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/03_Scripts/")

    # Path to the processed dataset (qnet and ts fields, full, time x lat x lon)
    datpath =  "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/"
    figpath =  "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/02_Figures/20220511/"
from amv import proc,viz
import scm


import amv.proc as hf # Update hf with actual hfutils script, most relevant functions

#%% ENSO Calculation and Cropping Options

# Select time crop (prior to preprocessing)
croptime          = True # Cut the time prior to detrending, EOF, etc
tstart            =  '1979-01-01' #'1920-01-01'#'0001-01-01' # "2006-01-01" # 
tend              =  '2021-12-31' #'2005-12-31'#'2000-02-01' # "2101-01-01" # 
timestr           = "%sto%s" % (tstart[:4],tend[:4])

# ENSO Parameters
pcrem             = 3                   # PCs to calculate
bbox              = [120, 290, -20, 20] # ENSO Bounding Box

# Toggles and Options
overwrite        = True # Set to True to overwrite existing output...
save_netcdf      = True # Set true to save netcdf version, false to save npz
debug            = True # Debug toggle

#%% Dataset option (load full TS variable in [time x lat x lon360])
# Example provided below here is for CESM1

# Data Information
# dataset_name        = 'cesm1_htr_5degbilinear'#"cesm2_pic"
# datpath             = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/proc/"#"/stormtrack/data4/glliu/01_Data/CESM2_PiControl/FCM/atm/"
# vname               = "TS"
# lonname             = "lon"
# latname             = "lat"
# timename            = "time"
# concat_dim          = "ens" #"time"
# keepvars            = [timename,latname,lonname,vname]
# ensnum              = 1 # Irrelevant for now, need to add ensemble support...
# detrend             = 1 # 1 to remove linear trend 
maskpath            = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/masks/"
maskname            = 'cesm1_htr_5degbilinear_limask_0.3p_0.05p_year1920to2005_enssum.nc'#"cesm2_pic_limask_0.3p_0.05p.nc"



# NOAA OISST
dataset_name        = 'OISST'#"cesm2_pic"
datpath             = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/NATL_proc_obs/"
vname               = "sst"
lonname             = "lon"
latname             = "lat"
timename            = "time"
concat_dim          = None#"time"
keepvars            = [timename,latname,lonname,vname]
ensnum              = 1 # Irrelevant for now, need to add ensemble support...
detrend             = 1 # 1 to remove linear trend 
outpath             = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/"


# ERA5 SST
dataset_name        = 'ERA5'#"cesm2_pic"
datpath             = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/03_reemergence/data/NATL_proc_obs/"
vname               = "sst"
lonname             = "lon"
latname             = "lat"
timename            = "time"
concat_dim          = None#"time"
keepvars            = [timename,latname,lonname,vname]
ensnum              = 1 # Irrelevant for now, need to add ensemble support...
detrend             = 1 # 1 to remove linear trend 
outpath             = ""

# Mask Information (first run a maskmaker script/section such as that in preproc_CESM2_PiControl.py)
maskpath            = None
maskname            = None

# Output Path (Checks for an "enso" folder)
outpath             = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/"


#%%

# 1A. Load Variable ----------

# Create filename/list and load 
if dataset_name == "cesm2_pic":
    searchstr = "%s%s/*%s*.nc" % (datpath,vname,vname) # Searches for datpath + *LANDFRAC*.nc
elif dataset_name == "OISST": # Just grab tropical pacific
    searchstr = "%s%s*%s*TropicalPacific.nc" % (datpath,dataset_name,vname) # Searches for datpath + *LANDFRAC*.nc
elif dataset_name == "ERA5": # Tropical Pacific Box
    searchstr = "%s%s*%s*TropicalPacific*.nc" % (datpath,dataset_name,vname) # Searches for datpath + *LANDFRAC*.nc
else:
    searchstr = "%s%s*%s*.nc" % (datpath,dataset_name,vname) # Searches for datpath + dataset_name*LANDFRAC*.nc"
nclist    = glob.glob(searchstr)
nclist.sort()
print("Found %i files for %s" % (len(nclist),vname))

# Drop Unnecessary variables
if concat_dim is None or len(nclist) == 1: # Assume no concantenation is needed
    ds_all    = xr.open_dataset(nclist[0])
else:
    print("Concatenating files by dim <%s>" % (concat_dim))
    ds_all    = xr.open_mfdataset(nclist,concat_dim=concat_dim,combine='nested')
ds_all    = hf.ds_dropvars(ds_all,keepvars)
try:
    ds_all    = hf.fix_febstart(ds_all)
except:
    print("Warning, Time is not in datetime... converting")
    print("First timestep (pre conversion) is %s" % (ds_all[timename][0]))
    timeconv = pd.to_datetime(ds_all.time.data)
    ds_all[timename] = timeconv
    print("First timestep (post conversion) is %s" % (ds_all[timename][0]))
    
    #ds_all['time2'] = pd.to_datetime(ds_all.time.data)
    print("Warning: February Start Fix was not implemented")
    #print("First timestep is %s" % (ds_all[timename][0]))

# Load it
st        = time.time()
ds_all    = ds_all[vname]#.load()
print("Loaded in %.2fs" % (time.time()-st))

# 1B. Load and Apply Mask
if (maskpath is None) or (maskname is None):
    print("No mask will be applied")
else:
    st = time.time()
    mask = xr.open_dataset(maskpath+maskname).mask.load()
    ds_all = ds_all * mask
    print("Mask applied in %.2fs" % (time.time()-st))


# Set ensemble flag
lensflag = False
if 'ens' in list(ds_all.dims):
    print("Ensemble dimension detected!")
    lensflag = True

#%% 2. Preprocessing

"""

From this point on you should have:
    
    TS [time x lat x lon360]... Land-Ice Mask Applied

"""

# Select the Box and load
st    = time.time()
dsreg = ds_all.sel(lon=slice(bbox[0],bbox[1]),lat=slice(bbox[2],bbox[3]))
dsreg = dsreg.load()
print("Output Loaded in %.2fs" % (time.time()-st))

# Remove seasonal cycle
#ds_anom  = hf.xrdeseason(dsreg) # hf.xrdeseason is not working?
ds_anom  = proc.xrdeseason(dsreg,check_mon=False)

# Remove Trend if option is set
if detrend:
    
    if 'ens' in list(ds_anom.dims):
        print("Detrending by removing ensemble mean")
        # Detrend by removing ensemble average
        da = ds_anom - ds_anom.mean('ens')
        
        da = da.transpose('ens','time','lat','lon')
        
    else:
        print("Detrending by removing linear fit")
        ds_anom   = ds_anom.transpose('time','lat','lon')
        
        # Simple Linear Detrend
        dt_dict   = hf.detrend_dim(ds_anom.values,0,return_dict=True)# ASSUME TIME in first axis
        
        # Put back into DataArray
        da = xr.DataArray(dt_dict['detrended_var'],dims=ds_anom.dims,coords=ds_anom.coords,name=vname)

else:
    da = ds_anom.copy()

print("Data preprocessed in %.2fs" % (time.time()-st))



#%% Part 3, Compute ENSO Indices (copied from calc_ENSO_general)

# ------------------- -------- General Portion --------------------------------

"""

IN : ncfile, <dataset_name>_<vname>_manom_detrend#.nc
    Anomalized, detrended ts with landice masked applied

OUT : npz file <dataset_name>_ENSO_detrend#_pcs#.npz
    PC File containing:
        eofall (ENSO EOF Patterns)          [lon x lat x month x pc]
        pcall  (ENSO principle components)  [(ens) x time x month x pc]
        varexpall (ENSO variance explained) [month x pc]]
        lon,lat,time,ensobbox variables

ex: ncep_ncar_ts_manom_detrend1.nc --> ncep_ncar_ENSO_detrend1_pcs3.npz

"""

st = time.time()

# # Open the dataset
# savename = "%s%s_%s_manom_detrend%i_%s.nc" % (datpath,dataset_name,vnames_in[0],detrend,timestr)
# if lensflag:
#     savename = proc.addstrtoext(savename,"_ens%02i"%(ensnum),adjust=-1)
# da = xr.open_dataset(savename)

# Slice to region
#da = da.sel(lon=slice(bbox[0],bbox[1]),lat=slice(bbox[2],bbox[3]))

# Check if ENSO has already been calculated and skip if so
proc.makedir("%senso/"% datpath) 
savename = "%senso/%s_ENSO_detrend%i_pcs%i_%s.nc" % (outpath,dataset_name,detrend,pcrem,timestr)
# if lensflag:
#     savename = proc.addstrtoext(savename,"_ens%02i"%(ensnum),adjust=-1)
query = glob.glob(savename)
if (len(query) < 1) or (overwrite == True):
    
    # Read out the variables # [(ens) x time x lat x lon]
    st        = time.time()
    invar     = da.values
    lon       = da[lonname].values
    lat       = da[latname].values
    times     = da[timename].values
    print("Data loaded in %.2fs"%(time.time()-st))
    
    if lensflag:
        nens = len(da.ens)
        
        eofall_ens      = []
        pcall_ens       = []
        varexpall_ens   = []
        for e in range(nens):
            invar_ens              = invar[e,...]
            # Portion Below is taken from calc_ENSO_PIC.py VV ***********
            eofall,pcall,varexpall = scm.calc_enso(invar_ens,lon,lat,pcrem,bbox=bbox)
            
            eofall_ens.append(eofall.copy())
            pcall_ens.append(pcall.copy())
            varexpall_ens.append(varexpall.copy())
        
        eofall    = np.array(eofall_ens)
        pcall     = np.array(pcall_ens)
        varexpall = np.array(varexpall_ens)
        
    else:
        # Portion Below is taken from calc_ENSO_PIC.py VV ***********
        eofall,pcall,varexpall = scm.calc_enso(invar,lon,lat,pcrem,bbox=bbox)
    
    # Sanity Check
    if debug:
        im = 0
        ip = 0
        proj = ccrs.PlateCarree(central_longitude=180)
        fig,ax = plt.subplots(1,1,subplot_kw={'projection':proj})
        ax = viz.add_coast_grid(ax,bbox=bbox)
        if lensflag:
            plotvar = eofall[0,:,:,im,ip]
        else:
            plotvar = eofall[:,:,im,ip]
            
        pcm = ax.pcolormesh(lon,lat,plotvar,vmin=-1,vmax=1,
                            cmap='cmo.balance',transform=ccrs.PlateCarree())
        cb = fig.colorbar(pcm,ax=ax,orientation='horizontal',fraction=0.055,pad=0.1)
        cb.set_label("SST Anomaly ($\degree C \sigma_{ENSO}^{-1}$)")
        ax.set_title("EOF %i, Month %i\n Variance Explained: %.2f" % (ip+1,im+1,varexpall[im,ip]*100)+"%")
    
    if save_netcdf:
        
        mons    = np.arange(1,13,1)
        
        years   = np.arange(int(len(times)/12))
        pcnums  = np.arange(1,pcrem+1)
        
        # Make Dictionary
        coords_eofs   = dict(lat=lat,lon=lon,month=mons,pc=pcnums) # 
        coords_pcs    = dict(year=years,month=mons,pc=pcnums)
        coords_varexp = dict(month=mons,pc=pcnums)
        if lensflag:
            ens     = np.arange(1,nens+1,1)
            # Unpack and repack dict to append item to start # https://www.geeksforgeeks.org/python-append-items-at-beginning-of-dictionary/
            coords_eofs,coords_pcs,coords_varexp = [{**{'ens':ens},**dd} for dd in [coords_eofs,coords_pcs,coords_varexp]]
        
        
        da_eofs       = xr.DataArray(eofall,coords=coords_eofs,dims=coords_eofs,name='eofs')
        da_pcs        = xr.DataArray(pcall,coords=coords_pcs,dims=coords_pcs,name='pcs')
        da_varexp     = xr.DataArray(varexpall,coords=coords_varexp,dims=coords_varexp,name='varexp')
        
        # Merge everything
        da_out        = xr.merge([da_eofs,da_pcs,da_varexp])
        
        # Add Additional Variables
        da_out['time']      = times
        da_out['enso_bbox'] = bbox
        
        edict = proc.make_encoding_dict(da_out)
        da_out.to_netcdf(savename,encoding=edict)
        
    else:
        
        # Save Output
        np.savez(savename,**{
                 'eofs': eofall, # [(ens) x lon x lat x month x pc]
                 'pcs': pcall,   # [Year, Month, PC]
                 'varexp': varexpall,
                 'lon': lon,
                 'lat':lat,
                 'times':times,
                 'enso_bbox':bbox}
                )
    print("Data saved in %.2fs"%(time.time()-st))
else:
    print("Skipping. Found existing file: %s" % (str(query)))
# End Skip

#%% Plot the ENSO Index






