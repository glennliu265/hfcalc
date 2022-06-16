#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Compare damping LENS

Compare calculated heat flux feedback between different large ensembles
Copied from compare_damping_CESM1LE.py


Created on Thu Jun 16 15:48:54 2022

@author: gliu
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cmocean
import sys
from tqdm.notebook import trange, tqdm

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

def plot_hfdamping(e,plotlag,plotmon,lon,lat,damping,msk,ax=None,cints=None):
    
    if cints is None:
        cints = np.arange(-50,55,5)
    if ax is None:
        ax = plt.gca()
    
    # Select what to plot
    plotvar = damping[e,plotlag,plotmon,:,:]
    plotmsk = msk[e,plotlag,plotmon,:,:]
    
    # Flip Longitude
    lon1,plotvar1 = proc.lon360to180(lon,plotvar.T)
    _,plotmsk1    = proc.lon360to180(lon,plotmsk.T)
    
    # Plot the contours
    pcm = ax.contourf(lon1,lat,plotvar1.T,levels=cints,cmap='cmo.balance',extend='both')
    
    # Plot significant points
    viz.plot_mask(lon1,lat,plotmsk1,reverse=False,ax=ax,markersize=.75,color='gray')
    
    return ax,pcm

#%% User Edits

# Set the File Names


# Name of dataset
dataset_names = (
    "CESM1-LE_HTR",
    "CESM1-LE_RCP85",
    "GFDL_ESM2M_LENS",
    "CSIRO_MK36_LENS",
    )

# Name of netCDF
ncnames = (
    "CESM1-LE_NHFLX_Damping_raw.nc",
    "CESM1_rcp85_hfdamping_ensorem1_detrend1_2006to2101_allens.nc",
    "gfdl_esm2m_lens_hfdamping_ensorem1_detrend1_2006to2101_allens.nc",
    "csiro_mk36_lens_hfdamping_ensorem1_detrend1_2006to2101_allens.nc"
    )

# Path to ata
datpaths = (
    "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/",
    "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/CESM-HTR-RCP85/",
    "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/lens/",
    "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/lens/"
    ) 

# Degrees of freedom
dofs = (
        82,
        (2101-2006) -4 +1,
        (2101-1950) -4 +1,
        (2101-1850) -4 +1,
        )

yr_rngs = ("1920-2005",
           "2006-2100",
           "1950-2100",
           "1850-2100")



# Paths, Names
figpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/02_Figures/20220622/"
proc.makedir(figpath)

# Significance Test Parameters
p       = 0.20     # p-value
tails   = 2        # two-tailed or one-tailed test...

# Plotting Parameters
bboxplot =  [-80,0,5,65]
months   = [viz.return_mon_label(m+1,nletters='all') for m in range(12)]
print(months)

#%% Merge and Begin
mcname    = dataset_names
ndatasets = len(mcname)

dampings = []
rflxs    = []
rssts    = []
rhocrits = []
mssts    = []
mflxs    = []
malls    = []

lons     = []
lats     = []

enss     = []

for mc in tqdm(range(ndatasets)):
    
    # Load in the files as a dataset
    ds = xr.open_dataset(datpaths[mc]+ncnames[mc])


    # Optionally load in the numpy arrays
    damping = ds.nhflx_damping.values     # Damping [ens x lag x mon x lat x lon]
    rflx    = ds.sst_flx_crosscorr.values
    rsst    = ds.sst_autocorr.values
    
    if mc > 0: # Tranpose [ens x (mon) x (lag) x lat x lon] to [ens x lag x mon x lat x lon]
        damping = damping.transpose(0,2,1,3,4)
        rflx    = rflx.transpose(0,2,1,3,4)
        rsst    = rsst.transpose(0,2,1,3,4)
    
    # Append
    dampings.append(damping)
    rflxs.append(rflx)
    rssts.append(rsst)
    if mc < 1:
        lons.append(ds.longitude.values)
        lats.append(ds.latitude.values)
        enss.append(ds.ensemble.values)
    else:
        lons.append(ds.lon.values)
        lats.append(ds.lat.values)
        enss.append(ds.ens.values)
    
    
    if mc == 0:
        #lon     = ds.longitude.values
        #lat     = ds.latitude.values
        mon     = ds.month.values
        #ens     = ds.ensemble.values
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
    
#%% Select Ensemble and Month and Plot

ie = 0 # Member
im = 0 # Month
il = 0 # Lag

fsz = 16

cints = np.arange(-45,47.5,2.5)
cmap  = 'cmo.balance'

proj    = ccrs.PlateCarree(central_longitude=0)

for im in tqdm(range(12)):
    fig,axs = plt.subplots(2,3,subplot_kw={'projection':proj},figsize=(16,8),
                           constrained_layout=True)
    
    for mc in range(2):
        
        for il in range(3):
            
            ax = axs[mc,il]
            
            # Get Mask and Variable
            plotvar = dampings[mc][ie,il,im,:,:]
            plotmsk = malls[mc][ie,il,im,:,:]
            
            #Labeling, Stuffs
            blabel = [0,0,0,0]
            if mc == 1:
                blabel[-1]=1
            else:
                ax.set_title("Lag %i"%(il+1),fontsize=fsz)
            if il == 0:
                blabel[0] =1
                
                ax.text(-0.18, 0.45, '%s'% (mcname[mc]), va='bottom', ha='center',
                    rotation='horizontal', rotation_mode='anchor',
                    transform=ax.transAxes,fontsize=14)
            ax = viz.add_coast_grid(ax,bbox=bboxplot,proj=proj,fill_color='gray',
                                    blabels=blabel)
            
            # Plot it
            lon1,plotvar1 = proc.lon360to180(lon,plotvar.T)
            _,plotmsk1    = proc.lon360to180(lon,plotmsk.T)
            
            pcm = ax.contourf(lon1,lat,plotvar1.T,levels=cints,
                              extend='both',cmap=cmap)
            
            viz.plot_mask(lon1,lat,plotmsk1,reverse=False,ax=ax,markersize=.75,color='gray')
            
    cb = fig.colorbar(pcm,ax=axs.flatten(),fraction=0.025,pad=0.01)
    cb.set_label("Heat Flux Feedback ($Wm^{-2} \degree C^{-1}$)",fontsize=14)
    plt.suptitle("Heat Flux Feedback for ENS %i, Month: %s"%(ie+1,months[im]),fontsize=24)
    
    savename = "%sHFF_Comparison_HTRvRCP85_ens%02i_mon%02i.png" % (figpath,ie+1,im+1)
    plt.savefig(savename,dpi=150,bbox_inches='tight')
        
#%% Examine Intermember Variability in Differences for a given month...

im = 0 
il = 0

mc = 1 # 0, 1, or "diff"

if mc == "diff":
    mcstr = "HTRvRCP85"
    cints = np.arange(-30,37.5,2.5)
elif mc == 0:
    mcstr = "HTR"
    cints = np.arange(-45,47.5,2.5)
elif mc == 1:
    mcstr = "RCP85"
    cints = np.arange(-45,47.5,2.5)


fig,axs = plt.subplots(5,8,figsize=(16,10),facecolor='white',constrained_layout=True,
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
    ax = viz.add_coast_grid(ax,bbox=bboxplot,proj=proj,blabels=blabel,
                            fill_color='gray',ignore_error=True)
    
    
    if mc == "diff":
        plotvar = dampings[1][e,il,im,:,:] - dampings[0][e,il,im,:,:]
        plotmsk = malls[1][e,il,im,:,:] * malls[0][e,il,im,:,:]
    else:
        plotvar = dampings[mc][e,il,im,:,:] 
        plotmsk = malls[mc][e,il,im,:,:] 

    # Flip Longitude
    lon1,plotvar1 = proc.lon360to180(lon,plotvar.T)
    _,plotmsk1    = proc.lon360to180(lon,(plotmsk*limask).T)

    pcm = ax.contourf(lon1,lat,plotvar1.T,levels=cints,cmap='cmo.balance',extend='both')


    # Plot significant points
    viz.plot_mask(lon1,lat,plotmsk1,reverse=False,ax=ax,markersize=.25,color='gray')
    
    # Add Subplot Labels
    ax = viz.label_sp(e+1,ax=ax,labelstyle="Ens%s",alpha=0.75,usenumber=True)
fig.colorbar(pcm,ax=axs.flatten(),orientation='horizontal',fraction=0.025,pad=0.01)
plt.suptitle("%s Lag %i Net Heat Flux Feedback ($Wm^{-2}K^{-1}$)" % (months[im],il+1),y=1.01,fontsize=16)
plt.savefig("%sNHFLX_Damping_%s_month%02i_lag%i.png" % (figpath,mcstr,im+1,il+1),
            dpi=200,bbox_inches='tight',transparent=False)


#%% For each of (4) models, plot the ensemble and annual average

[print(d.shape) for d in dampings]

imon = 0
ilag = 0
mode = "STD" # "STD" or "AVG"

proj = ccrs.PlateCarree(central_longitude=0)


if mode == "AVG":
    cints = np.arange(-45,47.5,2.5)
    cmap  = 'cmo.balance'
    cblbl = "Heat Flux Feedback ($Wm^{-2}$)" 
elif mode == "STD":
    cints = np.arange(0,26,1)
    cmap  = "inferno"
    cblbl = "1$\sigma$ Heat Flux Feedback ($Wm^{-2}$)" 
    
    
#cints = np.arange(-45,47.5,2.5)


fig,axs = plt.subplots(1,4,figsize=(14,4),
                       constrained_layout=True,
                       subplot_kw={'projection':proj})
for mc in range(4):
    
    ax = axs.flatten()[mc]
    blabel=[0,0,0,1]
    if mc == 0:
        blabel[0] = 1
    ax = viz.add_coast_grid(ax,bbox=bboxplot,proj=proj,blabels=blabel,
                            fill_color='gray',ignore_error=True)
    
    if mode == "AVG":
        plotvar = dampings[mc][:,ilag,imon,:,:].mean(0) 
        plotmsk = malls[mc][:,ilag,imon,:,:].mean(0) 
    elif mode == "STD":
        plotvar = dampings[mc][:,ilag,imon,:,:].std(0) 
        plotmsk = malls[mc][:,ilag,imon,:,:].std(0) 
        

    # Flip Longitude
    lon1,plotvar1 = proc.lon360to180(lons[mc],plotvar.T)
    _,plotmsk1    = proc.lon360to180(lons[mc],(plotmsk).T)
    
    pcm = ax.contourf(lon1,lats[mc],plotvar1.T,levels=cints,cmap=cmap,extend='both')
    
    ax.set_title("%s (%s)" % (dataset_names[mc],yr_rngs[mc]))
cb = fig.colorbar(pcm,ax=axs.flatten(),fraction=0.015,pad=0.01)
cb.set_label(cblbl)
plt.suptitle("Ensemble Average Heat Flux Feedback for %s, Lag %i" % (months[imon],ilag+1),y=0.90)
plt.savefig("%sEnsAvg_HFF_mon%02i_lag%i_ens%s.png" % (figpath,imon+1,ilag+1,mode),dpi=150,bbox_inches='tight')

#%% Plot for Clivar LENS, 30-member cases


mc   = 3
imon = 0
ilag = 0

fig,axs = plt.subplots(5,6,figsize=(14,8),
                       constrained_layout=True,
                       subplot_kw={'projection':proj})

for e in range(enss[mc][-1]):
    
    ax = axs.flatten()[e]
    
    
    blabel=[0,0,0,0]
    if e%5 == 0:
        blabel[0] = 1
    if e>24:
        blabel[-1] = 1
    
    ax = viz.add_coast_grid(ax,bbox=bboxplot,proj=proj,blabels=blabel,
                            fill_color='gray',ignore_error=True)
    
#     if mode == "AVG":
#         plotvar = dampings[mc][:,ilag,imon,:,:].mean(0) 
#         plotmsk = malls[mc][:,ilag,imon,:,:].mean(0) 
#     elif mode == "STD":
#         plotvar = dampings[mc][:,ilag,imon,:,:].std(0) 
#         plotmsk = malls[mc][:,ilag,imon,:,:].std(0) 
        

#     # Flip Longitude
#     lon1,plotvar1 = proc.lon360to180(lons[mc],plotvar.T)
#     _,plotmsk1    = proc.lon360to180(lons[mc],(plotmsk).T)
    
#     pcm = ax.contourf(lon1,lats[mc],plotvar1.T,levels=cints,cmap=cmap,extend='both')
    
#     ax.set_title("%s (%s)" % (dataset_names[mc],yr_rngs[mc]))
# cb = fig.colorbar(pcm,ax=axs.flatten(),fraction=0.015,pad=0.01)
# cb.set_label(cblbl)
plt.suptitle("%s Heat Flux Feedback (%s [%s], Lag %i)" % (months[imon],dataset_names[mc],yr_rngs[mc], ilag+1),y=0.90)
# plt.savefig("%s%s_HFF_mon%02i_lag%i_ens%s.png" % (figpath,dataset_names[mc],imon+1,ilag+1,mode),dpi=150,bbox_inches='tight')