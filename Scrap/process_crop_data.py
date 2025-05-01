#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Process and Crop the data for heat flux feedback calculations for the #SMIO project

Created on Tue Feb 18 13:19:28 2025

@author: gliu
"""

import os
import sys
import glob

import xarray as xr
import numpy as np
import cartopy.crs as ccrs
from tqdm import tqdm


#%% Module Loader

# stormtrack
amvpath = "/home/glliu/00_Scripts/01_Projects/00_Commons/" # amv module
scmpath = "/home/glliu/00_Scripts/01_Projects/01_AMV/02_stochmod/stochmod/model/" # scm module

sys.path.append(amvpath)
sys.path.append(scmpath)

from amv import proc,viz
import scm
import amv.loaders as dl

#%%


natl_box = [-100,20,-10,90]
enso_box = [120, 290, -20, 20] # Get ENSO Box (SST only)
outpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/03_reemergence/data/NATL_proc_obs/"




# =============================================
#%% Get ERA5, Monthly Surface Latent Heat Flux
# =============================================

# Get List of NetCDF files
vname    = 'slhf'
latname  = 'latitude'
lonname  = 'longitude'
timename = 'time'
dpath    = "/mnt/CMIP6/data/era5/reanalysis/single-levels/monthly-means/surface_latent_heat_flux/"
ncsearch = "%s*.nc" % dpath

vname_out = "lhflx"
outname   = outpath + "ERA5_lhflx_1979_2022_NATL.nc"

# Get List of Files, open dataset
nclist      = glob.glob(ncsearch)
nclist.sort()
dsall       = xr.open_mfdataset(nclist)

# Flip longitude
dsall180    = proc.format_ds(dsall,lonname=lonname,latname=latname)
dsall_natl  = proc.sel_region_xr(dsall180,natl_box)

# Rename and save
dsall_natl  = dsall_natl.rename({vname:vname_out})
edict       = proc.make_encoding_dict(dsall_natl)
dsall_natl.to_netcdf(outname,encoding=edict)

# # Correct silly issue where I named it wrong (shflx --> lhflx)
# oldname = outpath + "ERA5_shflx_1979_2022_NATL.nc"
# newname = oldname.replace('shflx','lhflx')
# dsload  = xr.open_dataset(oldname).load()
# dsload  = dsload.rename({'shflx':'lhflx'})
# edict       = proc.make_encoding_dict(dsload)
# dsload.to_netcdf(newname,encoding=edict)

# =============================================
#%% Do same thing for Monthly Sensible latent heat flux
# =============================================

vname       = "sshf"
vname_out   = "shflx"
dpath       = "/mnt/CMIP6/data/era5/reanalysis/single-levels/monthly-means/surface_sensible_heat_flux/"
ncsearch    = "%s*.nc" % dpath

outname   = outpath + "ERA5_%s_1979_2022_NATL.nc" % vname_out

# Get List of Files, open dataset
nclist      = glob.glob(ncsearch)
nclist.sort()
dsall       = xr.open_mfdataset(nclist)

# Flip longitude
dsall180    = proc.format_ds(dsall,lonname=lonname,latname=latname)
dsall_natl  = proc.sel_region_xr(dsall180,natl_box)

# Rename and save
dsall_natl  = dsall_natl.rename({vname:vname_out})
edict       = proc.make_encoding_dict(dsall_natl)
dsall_natl.to_netcdf(outname,encoding=edict)

# =============================================
#%% now repeat for ERA5 SST
# =============================================

vname           = "sst"
vname_out       = "sst"
latname         = 'latitude'
lonname         = 'longitude'
timename        = 'time'

dpath           = "/mnt/CMIP6/data/era5/reanalysis/single-levels/monthly-means/sea_surface_temperature/"
ncsearch        = "%s*.nc" % dpath
#outname         = outpath + "ERA5_%s_1979_2022_NATL.nc" % vname_out #Note this is the old naming convention, modified below to match HFF preprocessing
#outname         = outpath + #"ERA5_sst_NAtl_1979to2021.nc"
outpath2    = "/stormtrack/data4/glliu/01_Data/Reanalysis/ERA5/"
outname     = "sst_1979_2024.nc"

# Get List of Files, open dataset
nclist          = glob.glob(ncsearch)
nclist.sort()
nclist          = nclist[:-1] # Remove 2022 since it only goes up until Sept
dsall           = xr.open_mfdataset(nclist)

# Flip longitude
dsall180    = proc.format_ds(dsall,lonname=lonname,latname=latname)
dsall_natl  = proc.sel_region_xr(dsall180,natl_box)

# Rename and save
dsall_natl  = dsall_natl.rename({vname:vname_out})
edict       = proc.make_encoding_dict(dsall_natl)
dsall_natl.to_netcdf(outname,encoding=edict)

# Also Crop the tropical Pacific
dsall360  = proc.format_ds(dsall,lonname=lonname,latname=latname,lon180=False)
dsall_trop = proc.sel_region_xr(dsall360,enso_box)
edict      = proc.make_encoding_dict(dsall_trop)
outname    = outpath + "ERA5_sst_TropicalPacific_1979to2021.nc"
dsall_trop.to_netcdf(outname,encoding=edict)


# Also save global version
outname_global  = outpath + "ERA5_sst_Global_1979to2024.nc"
edict           = proc.make_encoding_dict(dsall180)
dsall180.to_netcdf(outname_global,encoding=edict)
# =============================================
#%% Now do for OISST
# =============================================

dpath      = "/vortex/jetstream/climate/data/yokwon/NOAA_OI_0.25deg_v2.1/processed/monthly/"
ncname     = "NOAA_OI_025deg_v2_1_AVHRR_monthly_1982_2020.nc"

lonname    = 'lon'
latname    = 'lat'



dsall      = xr.open_dataset(dpath+ ncname)

# Take ENSO
dsall_trop = proc.sel_region_xr(dsall,enso_box)
edict      = proc.make_encoding_dict(dsall_trop)
outname    = outpath + "OISST_SST_1982_2020_TropicalPacific.nc"
dsall_trop.to_netcdf(outname,encoding=edict)


# Take the North Atlantic Box
dsall180    = proc.format_ds(dsall.sst,lonname=lonname,latname=latname)
dsall_natl  = proc.sel_region_xr(dsall180,natl_box)
edict       = proc.make_encoding_dict(dsall_natl)
outname    = outpath + "OISST_SST_1982_2020_NATL.nc"
dsall_natl.to_netcdf(outname,encoding=edict)

# Need to do 2 crops (tropical pacific)
# North Atlantic
# =============================================
#%% Load the data above and reprocess further (based on specs in calc_hff_general)
# (Compute THFLX, Rename Files)
# =============================================

# I.E. Compute THFLX ----------------------------------------------------------
dpath       = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/03_reemergence/data/NATL_proc_obs/"
ds_lhflx    = xr.open_dataset(dpath + "ERA5_lhflx_1979_2022_NATL.nc").load()
ds_shflx    = xr.open_dataset(dpath + "ERA5_shflx_1979_2022_NATL.nc").load()

# Add LHFLX and SHFLX for THFLX, then resave
ds_thflx            = ds_lhflx.lhflx + ds_shflx.shflx
ds_thflx            = ds_thflx.rename('thflx')
ds_thflx            = ds_thflx.sel(time=slice('1982-01-01','2020-12-31'))
edict               = proc.make_encoding_dict(ds_thflx)
outname_thflx_new   = dpath + "ERA5_thflx_NAtl_1982to2020.nc"
ds_thflx.to_netcdf(outname_thflx_new,encoding=edict)


# Save version leading up to 2021
ds_thflx            = ds_lhflx.lhflx + ds_shflx.shflx
ds_thflx            = ds_thflx.rename('thflx')
ds_thflx            = ds_thflx.sel(time=slice('1979-01-01','2021-12-31'))
edict               = proc.make_encoding_dict(ds_thflx)
outname_thflx_new   = dpath + "ERA5_thflx_NAtl_1979to2021.nc"
ds_thflx.to_netcdf(outname_thflx_new,encoding=edict)

# Resave SST (OISST) with new naming conventions ------------------------------
ds_sst      = xr.open_dataset(dpath + "OISST_SST_1982_2020_NATL.nc").load()
outname_sst_new = dpath + "OISST_sst_NAtl_1982to2020.nc"
ds_sst = ds_sst.sel(time=slice('1982-01-01','2020-12-31'))
edict               = proc.make_encoding_dict(ds_sst)
ds_sst.to_netcdf(outname_sst_new,encoding=edict)

# Check lat lon
dx_oisst = ds_sst.lon.data[1:] - ds_sst.lon.data[:-1]
print(dx_oisst)
dx_era5  = ds_lhflx.lon.data[1:] - ds_lhflx.lon.data[:-1]
print(dx_era5)

# =============================================
#%% Regrid using ERA5 using xesmf
# =============================================
import xesmf as xe
import matplotlib.pyplot as plt

# Reload
dpath       = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/03_reemergence/data/NATL_proc_obs/proc/"
nc_sst      = dpath + "OISST_sst_NAtl_1982to2020.nc"
nc_thflx    = dpath + "ERA5_thflx_NAtl_1982to2020.nc"

ds_sst      = xr.open_dataset(nc_sst).load()
ds_flx      = xr.open_dataset(nc_thflx).load()


# Get Lat and Lon
lon_e = ds_flx.lon
lat_e = ds_flx.lat
lon_o = ds_sst.lon
lat_o = ds_sst.lat

# Visualize Grid ------- <0> -------

xxe,yye = np.meshgrid(lon_e,lat_e) # ERA5

xxo,yyo = np.meshgrid(lon_o,lat_o) # OISST

fig,ax = plt.subplots(1,1,figsize=(16,12),subplot_kw={'projection':ccrs.PlateCarree()})
ax.coastlines()
ax.set_extent([-60,-20,45,65])
ax.scatter(xxe,yye,c='r',s=0.2,marker="x")
ax.scatter(xxo,yyo,c='b',s=0.2,marker="d")
plt.show()

# ------- <0> ------- 

# Regrid (Copied form predict_nasst/regrid_ocean_variable_hmxl.py)
method    = 'bilinear'
ds_out    = xr.Dataset({'lat':lat_o,"lon":lon_o})
ds        = ds_flx#.to_dataset()

# Initialize Regridder
regridder = xe.Regridder(ds,ds_out,method,periodic=False)

# Regrid
daproc    = regridder(ds) # Need to input dataarray

edict     = proc.make_encoding_dict(daproc)
outname   = dpath + "ERA5_RegridOISST_thflx_NAtl_1982to2020.nc"
daproc.to_netcdf(outname,encoding=edict)

# =============================================
#%% Get Sea Ice Concentration (daily?), OISST
# =============================================
dpath   = "/vortex/jetstream/climate/data/yokwon/NOAA_OI_0.25deg_v2.1/downloaded/"
ncstr   = dpath + "icec.day.mean.*.nc"

nclist  = glob.glob(ncstr)

years   = np.arange(1981,2021)
nyrs    = len(nclist)

maxice     = []
meanice    = []
icefreq    = []

for yy in tqdm(range(nyrs)):
    
    year   = years[yy]
    ncname = dpath + "icec.day.mean.%s.nc" %  year
    ds     = xr.open_dataset(ncname).load()
    
    
    # dsmax  = ds.icec.max('time')
    # dsmax['year'] = year
    # maxice.append(dsmax.copy())  
    
    # Take Monthly Means
    dsmean         = ds.icec.groupby('time.month').mean('time')
    dsmean['year'] = year
    meanice.append(dsmean.copy())  
    
    # Get Frequency of Exceedence
    ds_exceed     = ds > 0.05
    ds_exceed     = ds_exceed.icec.sum('time')
    ds_exceed['year'] = year
    icefreq.append(ds_exceed.copy())
    
# # Save Max
# ds_cat      = xr.concat(maxice,dim='year')
# ds_cat      = proc.format_ds(ds_cat,timename='year')
# dpath       = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/03_reemergence/data/NATL_proc_obs/proc/"
# outname     = outpath + "OISST_icec_max_1981_2020.nc"
# edict       = proc.make_encoding_dict(ds_cat)
# ds_cat.to_netcdf(outname,encoding=edict)

# Save Mon Mean
ds_cat      = xr.concat(meanice,dim='year')
#ds_cat      = proc.format_ds(ds_cat,timename='year')
dpath       = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/03_reemergence/data/NATL_proc_obs/proc/"
outname     = outpath + "OISST_icec_monmean_1981_2020.nc"
edict       = proc.make_encoding_dict(ds_cat)
ds_cat.to_netcdf(outname,encoding=edict)

# Save IceFreq
ds_cat      = xr.concat(icefreq,dim='year')
ds_cat      = proc.format_ds(ds_cat,timename='year')
dpath       = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/03_reemergence/data/NATL_proc_obs/proc/"
outname     = outpath + "OISST_iceexceed005_monmean_1981_2020.nc"
edict       = proc.make_encoding_dict(ds_cat)
ds_cat.to_netcdf(outname,encoding=edict)

dsmax = ds.icec.max('time')

# =============================================
#%% Get CESM2 (copied from preproc_CESM2_PiControl)
# =============================================

datname       = "cesm2_pic"

datpath       = "/stormtrack/data4/glliu/01_Data/CESM2_PiControl/FCM/atm/" # Path to regridded mixed-layer depth


def load_cesm2_pic(vname,datpath,searchstr=None):
    # Function to Load CESM2 PiC Output on stormtrack, Searches for 
    # File is assumed to be at <datpath>/<vname>/*<vname>*.nc
    keepvars  = ["time","lat","lon",vname]
    if searchstr is None:
        searchstr = "%s%s/*%s*.nc" % (datpath,vname,vname) # Searches for datpath + *LANDFRAC*.nc
    nclist    = glob.glob(searchstr)
    nclist.sort()
    
    ds_all    = xr.open_mfdataset(nclist,concat_dim="time",combine='nested')
    ds_all    = proc.ds_dropvars(ds_all,keepvars)
    ds_all    = proc.fix_febstart(ds_all)
    
    return ds_all


ds_all = load_cesm2_pic("TS",datpath)

vname     = "TS"
keepvars  = ["time","lat","lon",vname]
searchstr = "%s%s/*%s*.nc" % (datpath,vname,vname) # Searches for datpath + *LANDFRAC*.nc
nclist    = glob.glob(searchstr)
nclist.sort()
print("Found %i files for %s" % (len(nclist),vname))

# Drop Unnecessary variables
ds_all    = xr.open_mfdataset(nclist,concat_dim="time",combine='nested')
ds_all    = proc.ds_dropvars(ds_all,keepvars)
ds_all    = proc.fix_febstart(ds_all)

# Load it
import time
st        = time.time()
ds_all    = ds_all[vname].load()
print("Loaded in %.2fs" % (time.time()-st))

# Reformat DS
ds_all = proc.format_ds(ds_all)

# Crop to NATL Region
ds_reg = proc.sel_region_xr(ds_all,natl_box)

outpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/03_reemergence/proc/"
outname = outpath+"cesm2_pic_TS_NAtl_0001_2000.nc"
edict   = proc.make_encoding_dict(ds_reg)
ds_reg.to_netcdf(outname,encoding=edict)

print(outname)


#%% Load and process land sea mask 
# **Note something seemed to be up with stormtrack so I moved this locally to work with (see section below)


# Get List of NetCDF files
vname    = 'lsm'
latname  = 'latitude'
lonname  = 'longitude'
timename = 'time'
dpath    = "/mnt/CMIP6/data/era5/reanalysis/single-levels/monthly-means/land_sea_mask/"
ncsearch = "%s*.nc" % dpath

vname_out = "mask"
outname   = outpath + "ERA5_land_sea_mask_1979_2022_NATL.nc"

# Get List of Files, open dataset
nclist      = glob.glob(ncsearch)
nclist.sort()
dsall       = xr.open_mfdataset(nclist)

# Flip longitude
dsall180    = proc.format_ds(dsall,lonname=lonname,latname=latname)
dsall_natl  = proc.sel_region_xr(dsall180,natl_box)

# Rename and save
dsall_natl  = dsall_natl.rename({vname:vname_out})
edict       = proc.make_encoding_dict(dsall_natl)
dsall_natl.to_netcdf(outname,encoding=edict)


#%% Local Processing

dpath2 = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/"
dslist = ["era5_1980_land_sea_mask.nc","era5_2020_land_sea_mask.nc"]

ds1 = xr.open_dataset(dpath2 + dslist[0]).load()
ds2 = xr.open_dataset(dpath2 + dslist[1]).load()
#dsall = xr.open_dataset(dpath2 + dslist[ii] for ii in range(2)).load()


ds_in       = [ds1,ds2]
dsall180    = [proc.format_ds(ds,lonname=lonname,latname=latname) for ds in ds_in]
dsall_natl  = [proc.sel_region_xr(ds,natl_box) for ds in dsall180]



mask_80 = xr.where(dsall_natl[0].lsm.max('time') > 0,np.nan,1)
mask_20 = xr.where(dsall_natl[1].lsm.max('time') > 0,np.nan,1)

mask_all = mask_80.data * mask_20.data
lat      = dsall_natl[0].lat.data
lon      = dsall_natl[0].lon.data


dims  = dict(lat=lat,lon=lon)
daout = xr.DataArray(mask_all,dims=dims,coords=dims,
                     name='land_mask')

edict       = proc.make_encoding_dict(daout)
outpath     = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/NATL_proc_obs/"
outname     = outpath + "ERA5_land_mask_1980_and_2020_NATL.nc"
daout.to_netcdf(outname,encoding=edict)