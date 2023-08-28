#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
- - - - - - - - - - - - - - - - - - - - - - |
||||||||||||||||||||||| ||||||||| ========= |
-----------------------
Compare Damping CESM1LE ||||||||| --------- |
-----------------------
||||||||||||||||||||||| ||||||||| ========= |
- - - - - - - - - - - - - - - - - - - - - - |

Compare Damping Between the CESM1 LENS Historical and RCP85 Periods

Scavenged code from
    - viz_hfdamping_CESM1LE.py
    
Created on Mon Jun  6 18:05:35 2022

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
use_htr_old = False # Use old HTR (not the one calculated by calc_enso_general)
expname     =  None # [None,30y_70to00,50y_firstlast]
vname       = "LHFLX"

if expname is None:
    
    # Historical
    if use_htr_old:
        #vname = 'nhflx'
        datpath   = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/"
        ncname    = "CESM1-LE_NHFLX_Damping_raw.nc"
        dof       = 82
    else:
        #vname = "qnet"
        datpath   = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/CESM-HTR-RCP85/"
        ncname    = "CESM1_htr_%s_damping_ensorem1_detrend1_1920to2005_allens.nc" % vname
        dof       = 82
    
    # RCP85
    datpath85 = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/CESM-HTR-RCP85/"
    if vname == "qnet":
        ncname85  = "CESM1_rcp85_hfdamping_ensorem1_detrend1_2006to2101_allens.nc"
    else:
        ncname85  = "CESM1_rcp85_%s_damping_ensorem1_detrend1_2006to2101_allens.nc" % vname
    dof85       = (2101-2006) - 4 + 1
    
elif expname is '30y_70to00':
    
    datpath   = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/CESM-HTR-RCP85/htr/"
    datpath85 = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/CESM-HTR-RCP85/rcp85/"
    
    dof       = 30
    dof85     = 30
    
    datpath   = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/CESM-HTR-RCP85/30y_70to00/%s_damping/htr/" % vname
    datpath85 = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/CESM-HTR-RCP85/30y_70to00/%s_damping/rcp85/" % vname
    
    ncname    = "CESM1_htr_%s_damping_ensorem1_detrend1_1970to2000_allens.nc" % vname
    ncname85  = "CESM1_rcp85_%s_damping_ensorem1_detrend1_2070to2100_allens.nc" % vname
    
# Paths, Names
figpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/02_Figures/20220815/"
proc.makedir(figpath)

# Significance Test Parameters
p        = 0.20     # p-value
tails    = 2        # two-tailed or one-tailed test...

# Plotting Parameters
bboxplot =  [-80,0,5,65]
months   = [viz.return_mon_label(m+1,nletters='all') for m in range(12)]
print(months)

#%% Merge and Begin
mcname   = ["HTR","RCP85"]
datpaths = [datpath,datpath85]
ncnames  = [ncname,ncname85]
dofs     = [dof,dof85]



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
    damping_name = "%s_damping" % vname
    if vname == "qnet":
        if use_htr_old or mc==1 and expname is None:
            damping_name = "nhflx_damping"
    
    # Load the values
    if (use_htr_old) or ((mc == 1) and (expname is None)):
        damping = ds[damping_name].values     # Damping [ens x lag x mon x lat x lon]
    else:
        damping = ds[damping_name].values
    rflx    = ds.sst_flx_crosscorr.values
    rsst    = ds.sst_autocorr.values
    
    if damping.shape[1] == 12: # Tranpose [ens x (mon) x (lag) x lat x lon] to [ens x lag x mon x lat x lon]
        damping = damping.transpose(0,2,1,3,4)
        rflx    = rflx.transpose(0,2,1,3,4)
        rsst    = rsst.transpose(0,2,1,3,4)
        
    # Append
    dampings.append(damping)
    rflxs.append(rflx)
    rssts.append(rsst)
    
    if mc == 1:
        lon     = ds.lon.values
        lat     = ds.lat.values
        mon     = ds.month.values
        ens     = ds.ens.values
        lag     = ds.lag.values
        
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
    
    savename = "%s%sHFF_Comparison_HTRvRCP85_ens%02i_mon%02i.png" % (figpath,vname,ie+1,im+1)
    plt.savefig(savename,dpi=150,bbox_inches='tight')
        
#%% Examine Intermember Variability in Differences for a given month...

im = 0 
il = 0
mc = "diff" # 0, 1, or "diff"

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

# -------------------------
# %% Ensemble Average Plots
# -------------------------
#%% Examine Ensemble Average (Annual Average), Global

il = 0
im = np.arange(0,12,1)

mode = "AVG" # {"AVG" or "1 Standard Deviation"}
add_clab = False

if mode == "AVG":
    titlestr = "Average"
    cmap = 'cmo.balance'

elif mode == "STD":
    titlestr = "1$\sigma$"
    cmap = 'inferno'


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
        if mode == "AVG":
            cints = np.arange(-20,21,1)
            plotvar = np.mean(dampings[1][:,il,im,:,:] - dampings[0][:-2,il,im,:,:],(0,1))
        elif mode == "STD":
            cints = np.arange(-6,6.5,.5)
            cmap = 'cmo.balance'
            plotvar = dampings[1][:,il,im,:,:].mean(1).std(0) - dampings[0][:-2,il,im,:,:].mean(1).std(0)
        plotmsk = np.prod((malls[1][:,il,im,:,:] * malls[0][:-2,il,im,:,:]),(0,1))

    else:
        mcstr = mcstrs[mc]
        if mode == "AVG":
            cints = np.arange(-45,47.5,2.5)
            plotvar = np.mean(dampings[mc][:,il,im,:,:],(0,1))
        elif mode == "STD":
            cints = np.arange(0,12.5,0.5)
            plotvar = dampings[mc][:,il,im,:,:].mean(1).std(0)
        plotmsk = np.prod(malls[mc][:,il,im,:,:],(0,1))
    
    
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
plt.savefig("%s%s_Damping_HTRvRCP85_%s_%s_lag%i_%s.png" % (figpath,vname,mode,"AnnAvg",il+1,expname),
            dpi=200,bbox_inches='tight',transparent=False)


#%% Visualize just the Qnet damping differences, but contour style
# matching the differences in vertical atmospheric variables

fig,ax  = plt.subplots(1,1,figsize=(12,6.5),facecolor='white',
                       subplot_kw={'projection':proj})

cints   = np.arange(-20,22,2)

plotvar = np.mean(dampings[1][:,il,im,:,:] - dampings[0][:-2,il,im,:,:],(0,1))
plotmsk = np.prod((malls[1][:,il,im,:,:] * malls[0][:-2,il,im,:,:]),(0,1))
plotvar,lon1 = add_cyclic_point(plotvar,coord=lon)
pcm     = ax.contourf(lon1,lat,plotvar,levels=cints,cmap=cmap,extend='both')
cl      = ax.contour(lon1,lat,plotvar,levels=cints,colors='k',linewidths=0.5)
ax.clabel(cl,levels=cints[::2])

ax      = viz.add_coast_grid(ax,proj=proj,bbox=bboxplot_glob,fill_color='k',ignore_error=True)

ax.set_title("Net Heat Flux Feedback Differences (RCP8.5 - Historical; Lag 1, Ens and Ann. Avg)")
cb      = fig.colorbar(pcm,ax=ax,fraction=0.018)
cb.set_label("%s Difference (%s) " % ("Net Heat Flux Feedback Difference","$W m^{-2}$"))
plt.savefig("%sNHFLX_Damping_HTRvRCP85_AVG_%s_lag%i_%s_contour.png" % (figpath,"AnnAvg",il+1,expname),
            dpi=200,bbox_inches='tight',transparent=False)

#%% Look at seasonal mean differences


hfdiffs      = dampings[1][:,il,:,:,:] - dampings[0][:-2,il,:,:,:] # [Ens x Mon x Lat x Lon]
hfdiffs      = hfdiffs.mean(0)
savgs,seastr = proc.calc_savg(hfdiffs,return_str=True,axis=0)

#%% Make the plot
cints = np.arange(-20,21,1)
fig,axs = plt.subplots(2,2,figsize=(16,6.5),facecolor='white',constrained_layout=True,
                       subplot_kw={'projection':proj})

for s in range(4):
    
    ax = axs.flatten()[s]
    
    blabel=[0,0,0,0]
    if s%2 == 0:
        blabel[0] = 1
    if s>1:
        blabel[-1] = 1
    ax = viz.add_coast_grid(ax,proj=proj,bbox=bboxplot_glob,blabels=blabel,
                            fill_color='k',ignore_error=True)
    
    plotvar      = savgs[s]
    plotvar,lon1 = add_cyclic_point(plotvar,coord=lon)
    pcm          = ax.contourf(lon1,lat,plotvar,levels=cints,cmap='cmo.balance',extend='both')
    cl           = ax.contour(lon1,lat,plotvar,levels=[0,],linecolors="k",linewidths=0.5)
    #ax.clabel(cl,inline_spacing=5)
    
    # Add Subplot Labels
    ax           = viz.label_sp(seastr[s],ax=ax,labelstyle="%s",alpha=0.75,usenumber=True,fontsize=25)

fig.colorbar(pcm,ax=axs.flatten(),orientation='horizontal',fraction=0.045,pad=0.01)
plt.suptitle("Seasonally-Averaged Net Heat Flux Feedback Differences(Lag %i, $Wm^{-2}K^{-1}$) \n RCP85 - HTR" % (il+1),y=1.10,fontsize=16)
plt.savefig("%sNHFLX_Damping_DIFFERENCES_SAVG_lag%i_%s.png" % (figpath,il+1,expname),
            dpi=200,bbox_inches='tight',transparent=False)
    