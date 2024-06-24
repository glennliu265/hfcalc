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
croptime          = False # Cut the time prior to detrending, EOF, etc
tstart            =  '0001-01-01' # "2006-01-01" # 
tend              =  '2000-02-01' #"2101-01-01" # 
timestr           = "%sto%s" % (tstart[:4],tend[:4])

# ENSO Parameters
pcrem             = 3                   # PCs to calculate
bbox              = [120, 290, -20, 20] # ENSO Bounding Box

# Toggles and Options
overwrite        = False # Set to True to overwrite existing output...
save_netcdf      = True # Set true to save netcdf version, false to save npz
debug            = True # Debug toggle

#%% Dataset option (load full TS variable in [time x lat x lon360])
# Example provided below here is for CESM1

# Data Information
dataset_name        = "cesm2_pic"
datpath             = "/stormtrack/data4/glliu/01_Data/CESM2_PiControl/FCM/atm/"
vname               = "TS"
lonname             = "lon"
latname             = "lat"
timename            = "time"
keepvars            = [timename,latname,lonname,vname]
lensflag            = False
ensnum              = 1 # Irrelevant for now, need to add ensemble support...
detrend             = 1 # 1 to remove linear trend 

# Mask Information (first run a maskmaker script/section such as that in preproc_CESM2_PiControl.py)
maskpath            = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/masks/"
maskname            = "cesm2_pic_limask_0.3p_0.05p.nc"

# Output Path (Checks for an "enso" folder)
outpath             = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/"

#%%

# 1A. Load Variable ----------

# Create filename/list and load 
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
ds_all    = ds_all[vname]#.load()
print("Loaded in %.2fs" % (time.time()-st))


# 1B. Load and Apply Mask
st = time.time()
mask = xr.open_dataset(maskpath+maskname).mask.load()
ds_all = ds_all * mask
print("Mask applied in %.2fs" % (time.time()-st))


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
ds_anom  = hf.xrdeseason(dsreg)
ds_anom   = ds_anom.transpose('time','lat','lon')

# Remove Trend if option is set
if detrend:
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
        pcall  (ENSO principle components)  [time x month x pc]
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
if lensflag:
    savename = proc.addstrtoext(savename,"_ens%02i"%(ensnum),adjust=0)
query = glob.glob(savename)
if (len(query) < 1) or (overwrite == True):
    
    # Read out the variables # [time x lat x lon]
    st        = time.time()
    invar     = da.values
    lon       = da[lonname].values
    lat       = da[latname].values
    times     = da[timename].values
    print("Data loaded in %.2fs"%(time.time()-st))
    
    # Portion Below is taken from calc_ENSO_PIC.py VV ***********
    eofall,pcall,varexpall = scm.calc_enso(invar,lon,lat,pcrem,bbox=bbox)
    
    # Sanity Check
    if debug:
        im = 0
        ip = 0
        proj = ccrs.PlateCarree(central_longitude=180)
        fig,ax = plt.subplots(1,1,subplot_kw={'projection':proj})
        ax = viz.add_coast_grid(ax,bbox=bbox)
        pcm = ax.pcolormesh(lon,lat,eofall[:,:,im,ip],vmin=-1,vmax=1,
                            cmap='cmo.balance',transform=ccrs.PlateCarree())
        cb = fig.colorbar(pcm,ax=ax,orientation='horizontal',fraction=0.055,pad=0.1)
        cb.set_label("SST Anomaly ($\degree C \sigma_{ENSO}^{-1}$)")
        ax.set_title("EOF %i, Month %i\n Variance Explained: %.2f" % (ip+1,im+1,varexpall[im,ip]*100)+"%")
    
    if save_netcdf:
        
        mons    = np.arange(1,13,1)
        years   = np.arange(pcall.shape[0])
        pcnums  = np.arange(1,pcrem+1)
        
        coords_eofs   = dict(lat=lat,lon=lon,month=mons,pc=pcnums) # 
        da_eofs       = xr.DataArray(eofall,coords=coords_eofs,dims=coords_eofs,name='eofs')
        
        coords_pcs    = dict(year=years,month=mons,pc=pcnums)
        da_pcs        = xr.DataArray(pcall,coords=coords_pcs,dims=coords_pcs,name='pcs')
        
        coords_varexp = dict(month=mons,pc=pcnums)
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
                 'eofs': eofall, # [lon x lat x month x pc]
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







