#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare Ice Fraction in CESM1 ans CESM2

Copied section from reemergence/viz_icefrac

Created on Mon Jun 24 13:23:24 2024

@author: gliu

"""


import xarray as xr
import numpy as np
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs

import matplotlib.pyplot as plt
import sys
import glob
import os

import tqdm
import time

# ----------------------------------
# %% Import custom modules and paths
# ----------------------------------

# Import re-eergemce parameters

# Indicate the Machine!
machine = "Astraeus"

# First Load the Parameter File
cwd = os.getcwd()
sys.path.append(cwd+ "/..")
import reemergence_params as rparams

# Paths and Load Modules
pathdict = rparams.machine_paths[machine]

sys.path.append(pathdict['amvpath'])
sys.path.append(pathdict['scmpath'])

# Set needed paths
figpath     = pathdict['figpath']
input_path  = pathdict['input_path']
output_path = pathdict['output_path']
procpath    = pathdict['procpath']
rawpath     = pathdict['raw_path']

#%% Import Custom Modules

# Import AMV Calculation
from amv import proc,viz
import amv.loaders as dl

# Import stochastic model scripts
import scm

#%% PLotting Parameters

bboxplot                    = [-80,0,10,65]
mpl.rcParams['font.family'] = 'Avenir'
mons3                       = proc.get_monstr(nletters=3)
fsz_tick                    = 18
fsz_axis                    = 14
fsz_title                   = 16
rhocrit                     = proc.ttest_rho(0.05,2,86)

proj                        = ccrs.PlateCarree()
bboxice = [-70,-10,55,70]

#%% Load the Ice Fractions



ncnames  = ["cesm2_pic_ICEFRAC_NAtl_0001to2000.nc","CESM1LE_ICEFRAC_NAtl_19200101_20050101_bilinear.nc"]
datnames = ["CESM2 PiControl", "CESM1 LENs"]


ds_all   = [xr.open_dataset(rawpath + nc).load() for nc in ncnames]



#%% Part (1), see if there is a difference in ice mask if different periods are
# used in CESM2

ystart              = 300

for ystart in [80,90]:
        
    ds_in               = ds_all[0].ICEFRAC
    
    
    icefrac_max_full    = ds_in.max('time')
    
    icefrac_tcrop       = ds_in.sel(time=slice('%04i-01-01' % ystart,'2000-12-31'))
    icefrac_max_crop    = icefrac_tcrop.max('time')
    
    
    
    plot_frac           = [icefrac_max_full,icefrac_max_crop]
    icemasks            = [xr.where(ds>0.05,0,1) for ds in plot_frac]
    plotnames           = ["Full Time Period (year 0-2000)","Year %04i Onwards" % ystart]
    
    #% Plot the Difference -------
    
    
    
    fig,axs,_ = viz.init_orthomap(1,2,bboxice,figsize=(14,4.5),constrained_layout=True)
    
    for ii in range(2):
        ax       = axs[ii]
        #im       = sel_mons[ii]
        #rei_sss  = rei_byvar[1].isel(mon=im).mean('yr').mean('ens')
        
        ax       = viz.add_coast_grid(ax,bbox=bboxice,fill_color='lightgray')
        
        # Plot Ice Concentration
        pv       = plot_frac[ii]#icecycle.isel(month=im).mean('ensemble')
        pcm      = ax.pcolormesh(pv.lon,pv.lat,pv,cmap='cmo.ice',transform=proj,zorder=-1,
                                 vmin=0,vmax=1)
        
        # Pot ice Edge
        icemask = icemasks[ii]
        ax.contour(icemask.lon,icemask.lat,icemask,colors="red",linewidths=1.5,linestyles='dotted',
                   levels=[0,1],zorder=1,transform=proj,label="Ice Mask Edge")
        
        
        ax.set_title(plotnames[ii],fontsize=fsz_title)
        
    
    
    cb = viz.hcbar(pcm,ax=axs.flatten(),fraction=0.045)
    cb.set_label("Ice Fraction (%)",fontsize=fsz_axis)
    plt.suptitle("Max Ice Fraction and Ice Mask in CESM2",fontsize=fsz_title)
    
    
    savename = "%sIceMaskComparison_CESM2_crop_y%04i.png" % (figpath,ystart)
    plt.savefig(savename,dpi=150,bbox_inches='tight')



