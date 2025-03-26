#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Double Check/Compare ENSO indices computed with [calc_enso_general]
and [compute_enso_index]

Copied upper section of the latter

Created on Tue Jun 18 11:58:11 2024

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

# Output Path (Checks for an "enso" folder)
outpath             = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/"

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



#%% Load ENSO Indices

ncnew    = outpath + "enso/cesm2_pic_ENSO_detrend1_pcs3_0001to2000.nc" 
ds_new = xr.open_dataset(ncnew).load()

nc_old = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/anom/enso/" + "cesm2_pic_ENSO_detrend1_pcs3_1to2000.npz"
ld_old = np.load(nc_old,allow_pickle=True)
print(ld_old.files)


#%% Check ENSO indices

im          = 1
nn          = 0

pcs         = [ds_new['pcs'].isel(month=im,pc=nn),ld_old['pcs'][:,im,nn]]
plotnames   = ["new Script","old Script"]

fig,ax      = plt.subplots(1,1)
for ii in range(2):
    ax.plot(pcs[ii],label=plotnames[ii])
ax.legend()
plt.show()

print(np.nanmax(np.abs(pcs[0].values - pcs[1])))

#%% Check ENSO Patterns

im = 0
ip = 0

eofsall     = [ds_new['eofs'].values,ld_old['eofs']]
varexpalls  = [ds_new['varexp'].values,ld_old['varexp']]
proj        = ccrs.PlateCarree(central_longitude=180)

lon         = ld_old['lon']
lat         = ld_old['lat']

fig,axs     = plt.subplots(2,1,subplot_kw={'projection':proj})

for ii in range(2):
    ax          = axs[ii]
    eofall      = eofsall[ii]
    varexpall   =  varexpalls[ii]
    
    ax          = viz.add_coast_grid(ax,bbox=bbox)
    pcm = ax.pcolormesh(lon,lat,eofall[:,:,im,ip],vmin=-1,vmax=1,
                        cmap='cmo.balance',transform=ccrs.PlateCarree())
    cb          = fig.colorbar(pcm,ax=ax,orientation='horizontal',fraction=0.055,pad=0.1)
    cb.set_label("SST Anomaly ($\degree C \sigma_{ENSO}^{-1}$)")
    ax.set_title("%s EOF %i, Month %i\n Variance Explained: %.2f" % (plotnames[ii],ip+1,im+1,varexpall[im,ip]*100)+"%")

plt.show()

#%% 


