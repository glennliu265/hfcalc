#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Regrid CESM1 LENs to a specified resolution

Copy the following:
    - regrid_reanalysis_cesm1.py
    - regrid_ocean_variable.py
    - preproc_CESM1_LENS
    
Output
    variable [ens x time x lat x lon360]

Created on Tue Jul  2 09:16:01 2024

@author: gliu

"""

import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from tqdm import tqdm
import xarray as xr
import xesmf as xe

#%% Load custom modules

amvpath = "/home/glliu/00_Scripts/01_Projects/00_Commons/" # amv module
#scmpath = "/home/glliu/00_Scripts/01_Projects/01_AMV/02_stochmod/stochmod/model/" # scm module

sys.path.append(amvpath)
#sys.path.append(scmpath)

from amv import proc,viz
#import scm
import amv.loaders as dl

#%% User Edits

# Indicate variables to process
vnames  = ['HMXL',]#["LANDFRAC","ICEFRAC","SHFLX","LHFLX","FSNS","FLNS"]
new_res = 5 # New resolution (in degrees)
method  = 'bilinear'

# Indicate paths and names
datname = "cesm1_htr_%ideg%s" % (new_res,method)
outpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/proc/"

#%% Get other info and set up grid

# Other Information
ocn_vars    = ["HMXL","XMXL","SSS","SALT","TEMP","BSF","SSH"] # Common Ocean Variables
mnum        = np.concatenate([np.arange(1,36),np.arange(101,108)]) # Ensemble Numbers
nens        = len(mnum)
nvars       = len(vnames)

# Set up new grid
# NTS: maybe use xe.util.grid_global(new_res, new_res) <-- this gives cell centers, need to check effect
lon_new     = np.arange(0,360+new_res,new_res)
lat_new     = np.arange(-90,90+new_res,new_res)
ds_newgrid  = xr.Dataset({'lat': (['lat'], lat_new), 'lon': (['lon'], lon_new) })

#%%


# Looping for each variable
for vv in range(nvars):
    vname   = vnames[vv]
    
    if vname in ocn_vars:
        atm = False
    else:
        atm = True
    
    # Loopin for ensemble members (atmospheric)
    dsvar_all = []
    for e in tqdm(range(nens)):
        N       = mnum[e]
        
        # Load DS
        dsvar   = dl.load_htr(vname,N,atm=atm).load()
        
        # Additional step for ocean variables
        if atm is False:
            dsvar = dsvar.rename({"TLONG": "lon", "TLAT": "lat"})
            dsvar = dsvar.to_dataset()
            dsvar.transpose('time','nlat','nlon')
        else:
            dsvar.transpose('time','lat','lon')
        
        # Set up regridder and perform regridding
        regridder       = xe.Regridder(dsvar, ds_newgrid, method, periodic=True)
        dsvar_regrid    = regridder( dsvar )
        if atm is False:
            dsvar_regrid = dsvar_regrid[vname]
        
        # Append and continue
        dsvar_all.append(dsvar_regrid)
    
    # Concatenate along the ensemble dimension # (ens: 42, time: 1032, lat: 37, lon: 73)
    dsvar_all = xr.concat(dsvar_all,dim='ens')
    
    # Save output
    edict     = proc.make_encoding_dict(dsvar_all)
    savename  = "%s%s_%s.nc" % (outpath,datname,vname)
    dsvar_all.to_netcdf(savename,encoding=edict)
    print("Saved regridded %s to:\n %s" % (vname,savename))



