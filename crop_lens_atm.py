#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Crop CESM-LE Atmospheric

Crop atmospheric variables for vertical gradient analysis in CESM1-LE

Created on Thu Aug 11 09:42:54 2022

@author: gliu
"""

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import glob
import sys
from tqdm import tqdm 
#%%

vname                 = "T"
z_ref                 = 850 # in mb
exp                   = "rcp85" # [rcp85 or htr]
loopens               = np.arange(34,40) ##np.range(len(mnum))

# Time period to cut to
croptime              =  True # Time period slice
if exp == 'htr':
    tstart            =  '1920-01-01'  # "2006-01-01"
    tend              =  '2006-01-01'  #"2100-12-01"
elif exp == 'rcp85':
    tstart            =  "2006-01-01"
    tend              =  "2100-12-01"

# Import visualization modules
sys.path.append("/home/glliu/00_Scripts/01_Projects/00_Commons/")
sys.path.append("/home/glliu/00_Scripts/01_Projects/01_AMV/02_stochmod/stochmod/model/")
from amv import proc,viz

# Plotting Parameters
bboxplot_glob = [-180,180,-65,75]
proj          = ccrs.PlateCarree(central_longitude=0)

# Outpath
outpath     = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/atmvar/"
outpath_var = outpath + "%s/" % vname
proc.makedir(outpath_var)

#%% Functions (taken from preproc_CESM1_LENS.py)

# RCP85 Loader
def load_rcp85(vname,N,datpath=None,append=True):
    if datpath is None:
        datpath = "/vortex/jetstream/climate/data1/yokwon/CESM1_LE/downloaded/atm/proc/tseries/monthly/"
    
    # Append variable name to path
    if append: # Add vname to datpath
        vdatpath = "%s%s/" % (datpath,vname)
    else: # Otherwise, datpath is just the datpath
        vdatpath=datpath
        
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

def load_cesmlens(vname,N,exp,datpath=None,append=True,):
    if datpath is None:
        datpath = "/vortex/jetstream/climate/data1/yokwon/CESM1_LE/downloaded/atm/proc/tseries/monthly/"
    
    # Append variable name to path
    if append: # Add vname to datpath
        vdatpath = "%s%s/" % (datpath,vname)
    else: # Otherwise, datpath is just the datpath
        vdatpath=datpath
    
    if exp == 'rcp85':
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
    else:
        # Ensemble 1 has a different time
        if N == 1:
            fn = "%sb.e11.B20TRC5CNBDRD.f09_g16.%03i.cam.h0.%s.185001-200512.nc" % (vdatpath,N,vname)
        else:
            fn = "%sb.e11.B20TRC5CNBDRD.f09_g16.%03i.cam.h0.%s.192001-200512.nc" % (vdatpath,N,vname)
        ds = xr.open_dataset(fn)
        if N == 1:
            ds = ds.sel(time=slice("1920-02-01","2006-01-01"))
    return ds

#%% Set File Paths
if "Q" in vname:
    datpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/00_Commons/CESM1_LE/%s/" % vname
if "T" in vname:
    datpath_all = "/vortex/jetstream/climate/data1/yokwon/CESM1_LE/downloaded/atm/proc/tseries/monthly/%s/" % vname
    datpath_34  = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/00_Commons/CESM1_LE/%s/" % vname
else:
    datpath = "/vortex/jetstream/climate/data1/yokwon/CESM1_LE/downloaded/atm/proc/tseries/monthly/%s/" % vname

#%% Get Ensemble Membrs and time length 

if exp == "rcp85":
    mnum    = np.concatenate([np.arange(1,36),np.arange(101,106)])
    ntime = 1140
elif exp == "htr":
    mnum    = np.concatenate([np.arange(1,36),np.arange(101,108)])
    ntime = 1032

# # Make netCDF string search for each scenario
# if exp == "htr":
#     ncsearch = "b.e11.B20TRC5CNBDRD.f09_g16.*.cam.h0.%s.*.nc" % vname
# elif exp == 'rcp85':
#     ncsearch = "b.e11.BRCP85C5CNBDRD.f09_g16.*.cam.h0.%s.*.nc" % vname

# # Get file names
# # Load in the data
# if vname == "T": # Ens 34 was missing for T, redirect path.
#     nclist34  = glob.glob(datpath_34 + ncsearch)
#     nclistall = glob.glob(datpath_all + ncsearch)
#     nclist    = nclist34+nclistall
# else:å
#     globby = datpath+ ncsearch
#     nclist =glob.glob(globby)
# nclist = [nc for nc in nclist if "OIC" not in nc]
# nclist.sort()
# print("Found %i items" % (len(nclist)))
# # Note this part isn't even really needed...
#%% For each item...


# Test out with e = 0
#e = 0

for e in tqdm(loopens):
    N = mnum[e]
    
    
    
    # Load in the data
    if "T" in vname:
        if N in [34,35]:
            datpath = datpath_34
        else:
            datpath = datpath_all
    ds = load_cesmlens(vname,N,exp,datpath=datpath,append=False)
    
    # Crop to Time Period and Atmospheric Levels (grab surface and reference)
    z_sfc     = ds.lev.isel(lev=-1).values
    ds_slice  = ds[vname].sel(lev=slice(z_ref,z_sfc),time=slice(tstart,tend))
    
    # Slice to the desired time period
    ds_sfc    = ds[vname].sel(lev=z_sfc,method='nearest')
    ds_ref    = ds[vname].sel(lev=z_ref,method='nearest')
    ds_diff   = ds_sfc - ds_ref
    
    # Recreate Attribute Dictionary
    attr_dict= ds[vname].attrs
    attr_dict['z_ref'] = ds_ref.lev.values
    attr_dict['z_sfc'] = ds_sfc.lev.values
    ds_diff.attrs = attr_dict
    
    outnc = "%s%sdiff_%s_%sto%s_ens%02i.nc" % (outpath_var,vname,exp,tstart[:4],tend[:4],e+1)
    ds_diff.to_netcdf(outnc,
             encoding={vname: {'zlib': True}})
    



# Plot the mean pattern for the selected period
# title = "$U_{%.2f}$ - $U_{%.2f}$ \n %s to %s Mean" % (ds_sfc.lev.values,
#                                            ds_ref.lev.values,
#                                            tstart,
#                                            tend)

# fig,ax = plt.subplots(1,1,figsize=(10,4.5),facecolor='white',constrained_layout=True,
#                        subplot_kw={'projection':proj})

# #pcm = ds_diff.mean('time').plot(ax=ax)
# pcm = ax.pcolormesh(ds.lon,ds.lat,ds_diff.mean('time'),cmap='cmo.balance',transform=proj)
# ax = viz.add_coast_grid(ax,proj=proj,bbox=bboxplot_glob,
#                         fill_color='k',ignore_error=True)
# ax.set_title(title)
# cb = fig.colorbar(pcm,ax=ax)
# cb.set_label("%s Differences (%s) " % (ds.U.attrs['long_name'],ds.U.attrs['units']))
# plt.show()


#ds_ref = ds.U.sel(lev=z_ref,method='nearest')





