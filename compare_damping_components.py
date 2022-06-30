#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare Components of Damping Across CESM1, Historical and RCP85

Copied script from 

Created on Tue Jun 28 16:46:10 2022

@author: gliu
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
import cmocean
import sys
from tqdm.notebook import trange, tqdm
import glob

# Interactive Plotting
import hvplot.xarray
import panel as pn
import panel.widgets as pnw

#%% Import my modules
sys.path.append("/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/03_Scripts/stochmod/model/")
sys.path.append("/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/03_Scripts/")

from amv import proc,viz
import scm

import importlib
importlib.reload(viz)
#%%

# Variables 
vnames    = ("FLNS","FSNS","LHFLX","SHFLX")
scenarios = ("htr","rcp85")
dofs      = (2005-1920+1,2100-2006+1)

# Significance Test Parameters
p        = 0.20     # p-value
tails    = 2        # two-tailed or one-tailed test...

# Add paths to the files
datpath   = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/CESM-HTR-RCP85/"
figpath   = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/02_Figures/20220629/"
proc.makedir(figpath)
#%% Merge and Begin

def load_damping(vname,scenario,datpath=None):
    
    """Loads data which has been processed by calc_enso_general.py
    and combine_damping_CESM1LE.py
    Returns each variable as  -- [ens x lag x mon x lat x lon]
    """""
    
    # Set Datapath and NetCDF Name
    if datpath is None:
        datpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/CESM-HTR-RCP85/"
    if scenario == "rcp85":
        ncname = "%sCESM1_rcp85_%s_damping_ensorem1_detrend1_2006to2100_allens.nc" % (datpath,vname)
    elif scenario == "htr":
        ncname = "%sCESM1_htr_%s_damping_ensorem1_detrend1_1920to2005_allens.nc" % (datpath,vname)
        
    # Open dataset
    ds           = xr.open_dataset(ncname)
    damping_name = "%s_damping" % vname
    
    # Load Numpy Arrays
    damping = ds[damping_name].values     # Damping [ens x lag x mon x lat x lon]
    rflx    = ds.sst_flx_crosscorr.values
    rsst    = ds.sst_autocorr.values
    
    outvars = [damping,rflx,rsst]
    
    # Adjust dimensions
    # Tranpose [ens x (mon) x (lag) x lat x lon] to [ens x lag x mon x lat x lon]
    if damping.shape[1] == 12:
        print("Transposing Array of shape %s" % str(damping.shape))
        outvars = [outvar.transpose(0,2,1,3,4) for outvar in outvars]
    return damping,rflx,rsst


def make_damping_mask(rflx,rsst,dof,p=0.05,tails=1):
    """
    Takes input of the form: [ens x lag x mon x lat x lon]
    """
    
    # Calculate critical threshold
    rhocrit = proc.ttest_rho(p,tails,dof)
    print("The critical rho value is %f" % rhocrit)
    
    # Make the masks
    msst = rsst > rhocrit
    mflx = rflx > rhocrit
    mall = msst * mflx 
    outvars = [rhocrit,mall,mflx,msst]
    
    return outvars
    
    
#%% 

vdampings  = {}
ccs        = {}
acs        = {}

mflxs      = {}
malls      = {}
rhocrits   = {}
mssts      = {}

for vname in tqdm(vnames):
    
    dampings = []
    rflxs    = []
    rssts    = []
    
    ccmask  = []
    allmask = []
    acmask  = []
    rhoall  = []
    
    for mc in range(2):
        scenario = scenarios[mc]
        outvars = load_damping(vname,scenario,datpath=datpath)
        
        dampings.append(outvars[0])
        rflxs.append(outvars[1])
        rssts.append(outvars[2])
        
        outvars2 = make_damping_mask(outvars[1],outvars[2],dofs[mc])
        rhoall.append(outvars2[0])
        allmask.append(outvars2[1])
        ccmask.append(outvars2[2])
        acmask.append(outvars2[3])
        
        
    
    vdampings[vname] = dampings
    ccs[vname]       = rflxs
    acs[vname]       = rssts
    
    malls[vname]    = allmask
    mssts[vname]    = acmask
    mflxs[vname]    = ccmask
    rhocrits[vname] = rhoall

# Load some additional values
testnc  = glob.glob(datpath+"*allens*.nc")[0]
ds      = xr.open_dataset(testnc)
lon     = ds.lon.values
lat     = ds.lat.values
mon     = ds.month.values
ens     = ds.ens.values
lag     = ds.lag.values
#%% Compute Ensemble Average Differences, Annual Mean

imon = np.arange(0,12,1)
ilag = 0

nens,nmon,nlag,nlat,nlon = vdampings[vname][0].shape
vensavg = np.zeros((4,2,nlat,nlon)) * np.nan
vensstd = np.zeros((4,2,nlat,nlon)) * np.nan
for v in range(4):
    vname = vnames[v]
    for mc in range(2):
        vensavg[v,mc,:,:] = np.mean(vdampings[vname][mc][:,imon,ilag,:,:],(0,1))
        vensstd[v,mc,:,:] = vdampings[vname][mc][:,imon,ilag,:,:].mean(1).std(0)
    print(vname)
    
#%% Visualize the total difference

proj          = ccrs.PlateCarree(central_longitude=0)
bboxplot_glob = [-180,180,-65,75]
vmaxmin       = 30
use_contour   = True
add_clab    = False

cints = np.arange(-20,21,1)

fig,ax = plt.subplots(1,1,figsize=(12,4),facecolor='white',constrained_layout=True,
                       subplot_kw={'projection':proj})
ax = viz.add_coast_grid(ax,proj=proj,bbox=bboxplot_glob,
                        fill_color='k',ignore_error=True)
plotvar = (vensavg[:,1,:,:] - vensavg[:,0,:,:]).sum(0)

if use_contour:
    pcm = ax.contourf(lon,lat,plotvar,levels=cints,cmap='cmo.balance')
else:
    if vmaxmin is None:
        vmaxmin = 2*np.nanstd(plotvar.flatten())
    pcm = ax.pcolormesh(lon,lat,plotvar,vmin=-vmaxmin,vmax=vmaxmin,cmap='cmo.balance')
cl = ax.contour(lon,lat,plotvar,levels=[0,],linecolors="k",linewidths=0.5)
if add_clab:
    ax.clabel(cl)
fig.colorbar(pcm,ax=ax,fraction=0.04,orientation='horizontal')
plt.savefig("%sCESM1LE_HFFDifferences_Summed_AnnAvg.png"% (figpath),dpi=150)
#%% Visualize some annual differences, but for each component flux

fig,axs = plt.subplots(2,2,figsize=(16,6.5),facecolor='white',constrained_layout=True,
                       subplot_kw={'projection':proj})

use_contour    = True

# 1 - Differences
# 2 - Percentage Diff RCP85/HTR
# 3 - Log Ratio Log(RCP85/HTR)
# 4 - Spatially Scaled Difference ((RCP85-HTR)/HTR)
viz_mode = 3

for v in range(4):
    ax = axs.flatten()[v]
    
    blabel=[0,0,0,0]
    if v%2 == 0:
        blabel[0] = 1
    if v>1:
        blabel[-1] = 1
    
    ax = viz.add_coast_grid(ax,proj=proj,bbox=bboxplot_glob,blabels=blabel,
                            fill_color='k',ignore_error=True)
    # Add Subplot Labels
    ax      = viz.label_sp("$\lambda_{%s}$" % vnames[v],ax=ax,
                      labelstyle="%s",alpha=0.75,usenumber=True,fontsize=20)
    if viz_mode == 1:
        plotvar = vensavg[v,1,:,:] - vensavg[v,0,:,:]
        cints   = np.arange(-4.5,4.75,0.25)
        title = "Annual Average Heat Flux Differences by Component (RCP85 - HTR)"
        plotval = 0
    elif viz_mode == 2:
        plotvar = np.abs(vensavg[v,1,:,:] / vensavg[v,0,:,:])
        cints   = np.arange(0,2.1,0.1)
        title = "Annual Average Heat Flux Ratio by Component [ RCP85/HTR ]"
        plotval = 1
    elif viz_mode == 3:
        plotvar = np.log(np.abs(vensavg[v,1,:,:] / vensavg[v,0,:,:]))
        title = "Annual Average Heat Flux Log Ratio by Component [ log(RCP85/HTR) ]"
        cints   = np.arange(-3,3.1,0.1)
        plotval = 0
    elif viz_mode == 4:
        plotvar = (vensavg[v,1,:,:] - vensavg[v,0,:,:]) / vensavg[v,0,:,:]
        cints   = np.arange(0,2.1,0.1)
        title = "Annual Average Heat Flux Scaled Differences by Component [(RCP85 - HTR)/HTR]"
        plotval = 1

    if use_contour:
        pcm = ax.contourf(lon,lat,plotvar,levels=cints,cmap='cmo.balance',extend='both')
        cl = ax.contour(lon,lat,plotvar,levels=[plotval,],linecolors="k",linewidths=0.5)
    else:
        vmaxmin = 2*np.nanstd(plotvar.flatten())
        pcm = ax.pcolormesh(lon,lat,plotvar,vmin=-vmaxmin,vmax=vmaxmin,cmap='cmo.balance')
        fig.colorbar(pcm,ax=ax,fraction=0.04,orientation='horizontal')

if use_contour:
    fig.colorbar(pcm,ax=axs.flatten(),fraction=0.015,pad=0.01)
    
plt.suptitle(title,fontsize=20)
plt.savefig("%sCESM1LE_HFFDifferences_ByComponent_AnnAvg_vizmode%i.png"% (figpath,viz_mode),dpi=150)
#%% Visualize the Difference for a component

il = 0
v  = 3 # Chose the variable
im = np.arange(0,12,1)

mode = "AVG" # {"AVG" or "1 Standard Deviation"}
add_clab = False

vname = vnames[v]

if mode == "AVG":
    titlestr = "Average"
    cmap = 'cmo.balance'

elif mode == "STD":
    titlestr = "1$\sigma$"
    cmap = 'inferno'
    
    
cintdict = {}
cintdict_diff = {}
if vname == "FLNS":
    
    cintdict['STD']      = np.arange(-10,10.5,.5)#**Needs setting
    cintdict_diff['STD'] = np.arange(-6,6.5,.5)  #**
    
    cintdict['AVG']      = np.arange(-10,10.5,0.5)
    cintdict_diff['AVG'] = np.arange(-5.,5.25,0.25)

elif vname == "FSNS":
    cintdict['STD']      = np.arange(-10,10.5,.5)
    cintdict_diff['STD'] = np.arange(-6,6.5,.5)
    
    cintdict['AVG']      = np.arange(-20,21,1)
    cintdict_diff['AVG'] = np.arange(-5.,5.25,0.25)
elif vname == "LHFLX":
    
    cintdict['STD']      = np.arange(-10,10.5,.5)
    cintdict_diff['STD'] = np.arange(-6,6.5,.5)
    
    cintdict['AVG']      = np.arange(-25,26,1)
    cintdict_diff['AVG'] = np.arange(-5.,5.25,0.25)
elif vname == "SHFLX":
    
    cintdict['STD']      = np.arange(-10,10.5,.5)
    cintdict_diff['STD'] = np.arange(-6,6.5,.5)
    
    cintdict['AVG']      = np.arange(-10,10.5,0.5)
    cintdict_diff['AVG'] = np.arange(-5.,5.25,0.25)

bboxplot_glob = [-180,180,-65,75]
proj          = ccrs.PlateCarree(central_longitude=0)
mcstrs        = ["HTR","RCP85"]

fig,axs = plt.subplots(3,1,figsize=(16,10),facecolor='white',constrained_layout=True,
                       subplot_kw={'projection':proj})

for mc in range(3):
    
    ax = axs.flatten()[mc]
    
    blabel=[1,0,0,0]
    if mc == 2:
        blabel[-1] = 1
    ax = viz.add_coast_grid(ax,proj=proj,bbox=bboxplot_glob,blabels=blabel,
                            fill_color='k',ignore_error=True)
    
    if mc == 2:
        mcstr = "RCP85 - HTR"
        cints = cintdict_diff[mode]
        if mode == "AVG":
            
            plotvar = vensavg[v,1,:,:] - vensavg[v,0,:,:]
        elif mode == "STD":
            cmap = 'cmo.balance'
            plotvar = vensstd[v,1,:,:] - vensstd[v,0,:,:]
        plotmsk = np.prod((malls[vname][1][:,il,:,:,:] * malls[vname][0][:-2,il,:,:,:]),(0,1))
    else:
        cints = cintdict[mode]
        mcstr = mcstrs[mc]
        if mode == "AVG":
            plotvar = vensavg[v,mc,:,:]
        elif mode == "STD":
            plotvar = vensstd[v,mc,:,:]
        plotmsk = np.prod(malls[vname][mc][:,il,:,:,:],(0,1))
    
    
    plotvar,lon1 = add_cyclic_point(plotvar,coord=lon)
    pcm = ax.contourf(lon1,lat,plotvar,levels=cints,cmap=cmap,extend='both')
    cl = ax.contour(lon1,lat,plotvar,levels=[0,],linecolors="k",linewidths=0.5)
    if add_clab:
        ax.clabel(cl)
        
    # Plot significant points
    viz.plot_mask(lon,lat,plotmsk.T,reverse=True,ax=ax,markersize=.2,
                  color='gray',proj=proj,geoaxes=True)
    
    
    # Add Subplot Labels
    ax = viz.label_sp(mcstr,ax=ax,labelstyle="%s",alpha=0.75,usenumber=True,fontsize=16)
    
    if mc == 1:
        fig.colorbar(pcm,ax=axs[:2],orientation='vertical',fraction=0.015,pad=0.01)
    elif mc == 2:
        fig.colorbar(pcm,ax=axs[2],orientation='vertical',fraction=0.015,pad=0.01)

plt.suptitle("%s %s Feedback (Lag %i, $Wm^{-2}K^{-1}$)" % ("Annual Average",vname,il+1),y=1.03,x=.67,fontsize=16)
plt.savefig("%s%s_Damping_HTRvRCP85_%s_%s_lag%i.png" % (figpath,vname,mode,"AnnAvg",il+1),
            dpi=200,bbox_inches='tight',transparent=False)

#%% Plot BSF



# BSF
datpath_bsf = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/01_Data/CESM_proc/"
bsf_ds = xr.open_dataset(datpath_bsf+"BSF_FULL_PIC_bilinear.nc")
#ds_reg = bsf_ds.BSF.sel(lon=slice(bboxplot[0],bboxplot[1]),lat=slice(bboxplot[2],bboxplot[-1]))
ds_reg   = bsf_ds.mean('time')
bsf      = ds_reg.BSF.values
bsf_mean = bsf.T
lon180 = ds_reg.lon.values





# Flip back to 360
_,bsf_mean_360 = proc.lon180to360(lon180,bsf_mean)


#%% BSF + LHFLX 

v = 2


cints = np.arange(-5,5.25,.25)

fig,ax  = plt.subplots(1,1,figsize=(12,4),facecolor='white',constrained_layout=True,
                       subplot_kw={'projection':proj})
ax      = viz.add_coast_grid(ax,proj=proj,bbox=bboxplot_glob,
                        fill_color='k',ignore_error=True)
plotvar = (vensavg[v,1,:,:] - vensavg[v,0,:,:])

if use_contour:
    pcm = ax.contourf(lon,lat,plotvar,levels=cints,cmap='cmo.balance',extend='both')
else:
    if vmaxmin is None:
        vmaxmin = 2*np.nanstd(plotvar.flatten())
    pcm = ax.pcolormesh(lon,lat,plotvar,vmin=-vmaxmin,vmax=vmaxmin,cmap='cmo.balance')
    
# Add BSF contours
cl = ax.contour(lon,lat,bsf_mean_360.T,
                levels=np.arange(-30,35,5),colors="k",linewidths=0.5)
if add_clab:
    ax.clabel(cl)
fig.colorbar(pcm,ax=ax,fraction=0.04,orientation='horizontal')
plt.savefig("%sCESM1LE_%sHFFDifferences_Summed_AnnAvg_BSF.png"% (figpath,vnames[v]),dpi=150)

# %% Plot all ensemble members for a given variable
v  = 2
im = 0 
il = 0
mc = "diff" # 0, 1, or "diff"
vname = vnames[v]
if mc == "diff":
    mcstr = "HTRvRCP85"
    cints = np.arange(-5,5.25,.25)
elif mc == 0:
    mcstr = "HTR"
    cints = np.arange(-45,47.5,2.5)
elif mc == 1:
    mcstr = "RCP85"
    cints = np.arange(-45,47.5,2.5)

fig,axs = plt.subplots(5,8,figsize=(24,7),facecolor='white',constrained_layout=True,
                       subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)})

for e in tqdm(range(40)):

    ax = axs.flatten()[e]
    
    #ax.coastlines()
    #ax.set_extent(bboxplot)
    blabel=[0,0,0,0]
    if e > 31:
        blabel[-1] = 1
    if e%8 == 0:
        blabel[0] = 1
    ax = viz.add_coast_grid(ax,proj=proj,blabels=blabel,
                            fill_color='k',ignore_error=True)
    
    #vdampings[vname][mc][e,i,il,:,:]
    if mc == "diff":
        plotvar = vdampings[vname][1][e,im,il,:,:] - vdampings[vname][0][e,im,il,:,:]
        plotmsk = malls[vname][1][e,im,il,:,:] * malls[vname][0][e,im,il,:,:]
    else:
         plotvar = vdampings[vname][mc][e,im,il,:,:]
         plotmsk = malls[vname][mc][e,im,il,:,:]
    #    plotmsk = malls[mc][e,il,im,:,:] 

    # # Flip Longitude
    # lon1,plotvar1 = proc.lon360to180(lon,plotvar.T)
    # _,plotmsk1    = proc.lon360to180(lon,(plotmsk*limask).T)

    pcm = ax.contourf(lon,lat,plotvar,levels=cints,cmap='cmo.balance',extend='both')


    # # Plot significant points
    #viz.plot_mask(lon,lat,plotmsk1,reverse=False,ax=ax,markersize=.25,color='gray')
    
    # Add Subplot Labels
    ax = viz.label_sp(e+1,ax=ax,labelstyle="Ens%s",alpha=0.75,usenumber=True)
fig.colorbar(pcm,ax=axs.flatten(),orientation='horizontal',fraction=0.045,pad=0.01)
#plt.suptitle("%s Lag %i Net Heat Flux Feedback ($Wm^{-2}K^{-1}$)" % (months[im],il+1),y=1.01,fontsize=16)
plt.savefig("%s%s_Damping_%s_month%02i_lag%i.png" % (figpath,vname,mcstr,im+1,il+1),
            dpi=200,bbox_inches='tight',transparent=False)


#%% Loop for each scenario






    
    


#%%

def load_damping(vname,scenario,datpath=None):
    

    
    dampings = []
    rflxs    = []
    rssts    = []
    rhocrits = []
    mssts    = []
    mflxs    = []
    malls    = []

    for mc in range(2):
        
        # Load in the files as a dataset
        ds = xr.open_dataset(datpaths[mc]+ncnames[mc])
    
    
        # Optionally load in the numpy arrays
        damping = ds.nhflx_damping.values     # Damping [ens x lag x mon x lat x lon]
        rflx    = ds.sst_flx_crosscorr.values
        rsst    = ds.sst_autocorr.values
        
        if mc == 1: # Tranpose [ens x (mon) x (lag) x lat x lon] to [ens x lag x mon x lat x lon]
            damping = damping.transpose(0,2,1,3,4)
            rflx    = rflx.transpose(0,2,1,3,4)
            rsst    = rsst.transpose(0,2,1,3,4)
            
        # Append
        dampings.append(damping)
        rflxs.append(rflx)
        rssts.append(rsst)
        
        if mc == 0:
            lon     = ds.longitude.values
            lat     = ds.latitude.values
            mon     = ds.month.values
            ens     = ds.ensemble.values
            lag     = ds.lag_month.values
            
            #% Make land mask
            limask = np.ones((192,288))
            limask[np.isnan((np.sum(damping,(0,1,2))))] = np.nan
            plt.pcolormesh(limask),plt.colorbar()
        
        #%Significance Test
        
        # Calculate critical threshold
        rhocrit = proc.ttest_rho(p,tails,dofs[mc])
        print("The critical rho value is %f" % rhocrit)
        rhocrits.append(rhocrit)
        
        # Make the masks
        msst = rsst > rhocrit
        mflx = rflx > rhocrit
        mall = msst * mflx 
        
        mssts.append(msst)
        mflxs.append(mflx)
        malls.append(mall)
