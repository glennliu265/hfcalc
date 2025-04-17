#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Compute Qnet in ERA5 reanalysis

Works with output downloaded using Copernicus CDI, and placed
into stormtrack in data4/ (see dpath)


Created on Tue Mar 25 14:36:37 2025

@author: gliu
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import tqdm as tqdm
import sys
import cartopy.crs as ccrs

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

# Indicate the variable names
vnames = ["slhf","sshf","ssr","str"]
dpath  = "/stormtrack/data4/glliu/01_Data/Reanalysis/ERA5/"
ystr   = "1979_2024"

# Load netcdf
ds_all = []
for vv in tqdm.tqdm(range(4)):
    ncname = "%s%s_%s.nc" % (dpath,vnames[vv],ystr)
    ds = xr.open_dataset(ncname).load()
    ds_all.append(ds)

#%% Format the DataFrame

latname         = 'latitude'
lonname         = 'longitude'
timename        = 'valid_time'

# Flip the Longitude
ds_all180       = [proc.format_ds(ds,lonname=lonname,latname=latname,timename=timename) for ds in ds_all]

# Select the North Atlantic Bounding Box
natl_box        = [-100,20,-10,90]
ds_all_natl     = [proc.sel_region_xr(ds,natl_box) for ds in ds_all180]

#%% Quick Sanity Check


t       = 36
dtdaily = 60*60*24 

fig,axs = plt.subplots(1,4,constrained_layout=True,figsize=(16,4.5),
                       subplot_kw={'projection':ccrs.PlateCarree()})

for vv in range(4):
    
    ax      = axs[vv]
    vname   = vnames[vv]
    
    plotvar = ds_all_natl[vv][vname].isel(time=t)/dtdaily
    
    pcm     = ax.pcolormesh(plotvar.lon,plotvar.lat,plotvar,
                            vmin=-250,vmax=250,cmap='cmo.balance')
    viz.hcbar(pcm,ax=ax)
    ax.set_title(vname)
plt.suptitle(plotvar.time.data)
plt.show()


"""
Note. Checking the winter values, it seems that all values are positive into the ocean 

(i.e. SLHF, SSHF, and STR are all negative in winter, while
 SSR is positive)
"""

#%% Compute a few values
# Start with qnet

da_all  = xr.merge([ds_all_natl[vv][vnames[vv]]/dtdaily for vv in range(4)])

flxes   = [da_all[vnames[vv]].data for vv in range(4)]
qnet    = np.stack(flxes).sum(0)
ds      = da_all[vnames[0]]
dims    = dict(time=ds.time,lat=ds.lat,lon=ds.lon) #da_all[vnames[0]].dims
da_qnet = xr.DataArray(qnet,dims=dims,coords=dims,name='qnet')

#%% Output files

outpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/03_reemergence/data/NATL_proc_obs/proc/"

# Save Net Fluxes
edict   = proc.make_encoding_dict(da_qnet)
ncname  = outpath + "ERA5_qnet_NAtl_1979to2024.nc"
da_qnet.to_netcdf(ncname,encoding=edict)

# Save all fluxes
edict   = proc.make_encoding_dict(da_all)
ncname  = outpath + "ERA5_All_Fluxes_NAtl_1979to2024.nc"
da_all.to_netcdf(ncname,encoding=edict)

# ============================================================================
#%% Reload and Save SST
# ============================================================================

# Load the data
dpath           = "/stormtrack/data4/glliu/01_Data/Reanalysis/ERA5/"
ystr            = "1979_2024"
vname           = "sst"
ncname          = "%s%s_%s.nc" % (dpath,vname,ystr)
ds              = xr.open_dataset(ncname).load()

latname         = 'latitude'
lonname         = 'longitude'
timename        = 'valid_time'

# Grab ENSO region
enso_box        = [120, 290, -20, 20] # Get ENSO Box (SST only)
ds360           = proc.format_ds(ds.sst,lonname=lonname,latname=latname,timename=timename,lon180=False)
dsall_trop      = proc.sel_region_xr(ds360,enso_box)
edict           = proc.make_encoding_dict(dsall_trop)
outname         = outpath + "ERA5_sst_TropicalPacific_1979to2024.nc"
dsall_trop.to_netcdf(outname,encoding=edict)

# Grab NAtl region
ds180           = proc.format_ds(ds.sst,lonname=lonname,latname=latname,timename=timename,lon180=True)
ds_natl         = proc.sel_region_xr(ds180,natl_box)
ncname          = outpath + "ERA5_sst_NAtl_1979to2024.nc"
ds_natl.to_netcdf(ncname,encoding=edict)


#%% Resave pre-satelite era chunk for analysis

outpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/03_reemergence/data/NATL_proc_obs/proc/"
dpath           = "/stormtrack/data4/glliu/01_Data/Reanalysis/ERA5/"

ncname = dpath + "sst_1940_1978.nc"


ds              = xr.open_dataset(ncname).load()

latname         = 'latitude'
lonname         = 'longitude'
timename        = 'valid_time'

# Grab NAtl region
ds180           = proc.format_ds(ds.sst,lonname=lonname,latname=latname,timename=timename,lon180=True)
ds_natl         = proc.sel_region_xr(ds180,natl_box)
edict           = proc.make_encoding_dict(ds_natl)
ncname          = outpath + "ERA5_sst_NAtl_1940to1978.nc"
ds_natl.to_netcdf(ncname,encoding=edict)



#%%



#ds_flip      = [proc.check_flx(ds_all_natl[vv],flxname=vnames[vv]) for vv in range(4)]
# Check the sign