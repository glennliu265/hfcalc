#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Visualize Heat Flux Feedback
----------------------------
Calculated from different reanalysis datasets
in calc_enso_general.py (need to really rename this...)

Created on Wed May 11 10:04:55 2022

@author: gliu
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import time
import sys
import cartopy.crs as ccrs

import cmocean as cmo
from tqdm import tqdm

#%% Import modules
stormtrack = 0

if stormtrack:
    sys.path.append("/home/glliu/00_Scripts/01_Projects/00_Commons/")
    sys.path.append("/home/glliu/00_Scripts/01_Projects/01_AMV/02_stochmod/stochmod/model/")
else:
    sys.path.append("/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/03_Scripts/stochmod/model/")
    sys.path.append("/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/03_Scripts/")
from amv import proc,viz
import scm


# List of NCs
datpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/"
nclist  = ("CESM1_FULL_PIC_hfdamping_ensorem1_detrend1.nc",
           "era20c_hfdamping_ensorem1_detrend3.nc",
           "noaa_20cr_v2_hfdamping_ensorem1_detrend3.nc",
           "ncep_ncar_hfdamping_ensorem1_detrend3.nc"
          )


nclist  = ("CESM1_FULL_PIC_hfdamping_ensorem1_detrend1.nc",
           'era20c_hfdamping_ensorem1_detrend1_1948to2010.nc',
           'noaa_20cr_v2_hfdamping_ensorem1_detrend1_1948to2014.nc',
           )

ncnames = ("CESM1-PiC",
           "ERA-20C Reanalysis",
           "NOAA-CIRES 20CR v2c",
           "NCEP-NCAR Reanalysis I"
           )

figpath  = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/02_Figures/20220602/"
proc.makedir(figpath)

# Plotting params
proj     = ccrs.PlateCarree(central_longitude=0)
bboxplot = [-80,0,10,62]
cints    = np.arange(-50,55,5)

vname    = 'nhflx_damping'
mons3    = proc.get_monstr(3)



#%% Load in the files

das      = []
for nc in nclist:
    da = xr.open_dataset(datpath+nc)
    das.append(da[vname]*-1)
    print(da)
    

#%% Make the plot for a specific month,lag

ilag = 0
imon = 0

for imon in tqdm(range(12)):
    fig,axs = plt.subplots(1,len(das),constrained_layout=True,figsize=(12,4),
                           subplot_kw={'projection':proj})
    
    for d in range(len(das)):
        
        ax =axs[d]
        
        blabel = [0,0,0,1]
        if d == 0:
            blabel[0] = 1
        
        ax = viz.add_coast_grid(ax,bbox=bboxplot,fill_color='gray',blabels=blabel)
        pcm = das[d].isel(month=imon,lag=ilag).plot.contourf(levels=cints,
                                                       ax=ax,
                                                       cmap=cmo.cm.balance,
                                                       add_colorbar=False,
                                                       title=False,extend='both')
        
        ax.set_title(ncnames[d])
        if d == 0:
            ax.text(-0.35, 0.55, "%s" % (mons3[imon]), va='bottom', ha='center',rotation='horizontal',
                    rotation_mode='anchor',transform=ax.transAxes)
    fig.colorbar(pcm,ax=axs.flatten(),fraction=0.009,pad=0.01)
    plt.savefig("%sHFF_Comparison_lag%i_mon%02i" % (figpath,ilag,imon),bbox_inches='tight',dpi=150)
        
        