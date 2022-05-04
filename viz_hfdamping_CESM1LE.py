#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Visualize Heat Flux Damping CESM1-LE

Python Script Version of 
- viz_hfdamping_CESM1LE.ipynb

Uses netCDF output processed from
- hfdamping_mat2nc.py

Includes the following plots
 - Pattern Correlation by Ensemble Member (HFF vs. Correlation)
 - HFF/Correlation Composites of top N events by Pattern Correlation

Created on Tue May  3 12:40:05 2022

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

# Paths, Names
datpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/"
figpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/02_Figures/20220502/"
ncname  = "CESM1-LE_NHFLX_Damping_raw.nc"
proc.makedir(figpath)

# Significance Test Parameters
dof     = 82       # Degrees of Freedom for Significaxnce Testing
p       = 0.20     # p-value
tails   = 2        # two-tailed or one-tailed test...

# Plotting Parameters
bboxplot =  [-80,0,5,60]
months   = [viz.return_mon_label(m+1,nletters='all') for m in range(12)]
print(months)

#%% Load in files

# Load in the files as a dataset
ds = xr.open_dataset(datpath+ncname)

# Optionally load in the numpy arrays
damping = ds.nhflx_damping.values     # Damping [ens x lag x mon x lat x lon]
rflx    = ds.sst_flx_crosscorr.values
rsst    = ds.sst_autocorr.values
lon     = ds.longitude.values
lat     = ds.latitude.values
mon     = ds.month.values
ens     = ds.ensemble.values
lag     = ds.lag_month.values


#%%Significance Test

# Calculate critical threshold

rhocrit = proc.ttest_rho(p,tails,dof)
print("The critical rho value is %f" % rhocrit)

#%% Make land mask
limask = np.ones((192,288))
limask[np.isnan((np.sum(damping,(0,1,2))))] = np.nan
plt.pcolormesh(limask),plt.colorbar()

# Make the masks
msst = rsst > rhocrit
mflx = rflx > rhocrit
mall = msst * mflx 

msst.shape,mflx.shape

#%%
# Make 42-ensemble plots
plotlag = 0
plotmon = 0

for m in range(12):
    plotmon=m
    fig,axs = plt.subplots(7,6,figsize=(16,16),facecolor='white',
                           subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)})
    cints = np.arange(-50,55,5)
    for e in tqdm(range(42)):

        ax = axs.flatten()[e]
        ax.coastlines()
        ax.set_extent(bboxplot)

        plotvar = damping[e,plotlag,plotmon,:,:]
        plotmsk = mall[e,plotlag,plotmon,:,:]

        # Flip Longitude
        lon1,plotvar1 = proc.lon360to180(lon,plotvar.T)
        _,plotmsk1    = proc.lon360to180(lon,(plotmsk*limask).T)
        #_,plotmsk1    = proc.lon360to180(lon,limask.T)

        pcm = ax.contourf(lon1,lat,plotvar1.T,levels=cints,cmap='cmo.balance',extend='both')


        # Plot significant points
        viz.plot_mask(lon1,lat,plotmsk1,reverse=False,ax=ax,markersize=.55,color='gray')

        # Add Subplot Labels
        ax = viz.label_sp(e+1,ax=ax,labelstyle="Ens%s",alpha=0.75,usenumber=True)
    fig.colorbar(pcm,ax=axs.flatten(),orientation='horizontal',fraction=0.025,pad=0.01)
    plt.suptitle("%s Lag %i Net Heat Flux Feedback ($Wm^{-2}K^{-1}$)" % (months[plotmon],plotlag+1),y=0.90,fontsize=16)
    plt.savefig("%sNHFLX_Damping_CESM1LE_mall_month%02i_lag%i.png" % (figpath,plotmon+1,plotlag+1),
                dpi=200,bbox_inches='tight',transparent=False)

#%% Do the Same thing, but for the cross/autocorrelations

vlims  = [-.5,1]




cints_all = [np.arange(-.5,.55,0.05),np.arange(-1,1.1,.1)]


invars = [rflx,rsst]

rnames = ["crosscorr","autocorr"]
rnames_fancy = ("NHFLX-SST Cross-correlation","SST Autocorrelation")
mskin  = [mflx,msst]

for v in range(2):
    for m in range(12):
        
        plotmon=m
        fig,axs = plt.subplots(7,6,figsize=(16,16),facecolor='white',
                               subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)})
        cints = np.arange(-50,55,5)
        
        for e in tqdm(range(42)):

            ax = axs.flatten()[e]
            ax.coastlines()
            ax.set_extent(bboxplot)
            
            plotvar = invars[v][e,plotlag,plotmon,:,:]
            plotmsk = mskin[v][e,plotlag,plotmon,:,:]
            cints   = cints_all[v]

            # Flip Longitude
            lon1,plotvar1 = proc.lon360to180(lon,plotvar.T)
            _,plotmsk1    = proc.lon360to180(lon,(plotmsk*limask).T)
            #_,plotmsk1    = proc.lon360to180(lon,limask.T)

            pcm = ax.contourf(lon1,lat,plotvar1.T,levels=cints,cmap='cmo.balance',extend='both')


            # Plot significant points
            viz.plot_mask(lon1,lat,plotmsk1,reverse=False,ax=ax,markersize=.55,color='gray')

            # Add Subplot Labels
            ax = viz.label_sp(e+1,ax=ax,labelstyle="Ens%s",alpha=0.75,usenumber=True)
        fig.colorbar(pcm,ax=axs.flatten(),orientation='horizontal',fraction=0.025,pad=0.01)
        plt.suptitle("%s Lag %i %s)" % (months[plotmon],plotlag+1,rnames_fancy[v]),y=0.90,fontsize=16)
        plt.savefig("%sNHFLX_SST_%s_CESM1LE_mall_month%02i_lag%i.png" % (figpath,rnames[v],plotmon+1,plotlag+1),
                    dpi=200,bbox_inches='tight',transparent=False)
    
#%% Look at the stdev in correlation for each one


plotnames_fancy = ("$\sigma_{Heat \, Flux \, Feedback}$",
                   "$\sigma_{Cross-correlation}$",
                   "$\sigma_{Autocorrelation}$")
invar_raw = [damping,rflx,rsst]
cints_all = (np.arange(0,27.5,2.5),np.arange(0,0.11,0.01),np.arange(0,0.055,0.005))
cmaps     = ('inferno','cmo.thermal','copper')


ilag      = 0 # Select the Lag

# Loop for each month
for imon in tqdm(range(12)):
    invars = [v[:,ilag,imon,:,:].std(0) for v in invar_raw]
    
    fig,axs = plt.subplots(1,3,figsize=(16,4),facecolor='white',
                           subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)})
    
    for i in range(3):
        
        ax      = axs[i]
        plotvar = invars[i]
        
        blabel = [0,0,0,1]
        if i == 0: # Set labels for left 
            blabel[0]  = 1
            ax.text(-0.28, 0.35, '42-member \n Stdev. \n (%s, \n Lag %i)'% (months[imon],ilag+1), va='bottom', ha='center',
                rotation='horizontal', rotation_mode='anchor',
                transform=ax.transAxes,fontsize=14)
        
        ax = viz.add_coast_grid(ax,bbox=bboxplot,blabels=blabel,fill_color='gray',ignore_error=True)
        
        lon1,plotvar1 = proc.lon360to180(lon,plotvar.T)
        
        pcm = ax.contourf(lon1,lat,plotvar1.T,levels=cints_all[i],
                          extend='both',cmap=cmaps[i])
        ax.set_title(plotnames_fancy[i])
        fig.colorbar(pcm,ax=ax,fraction=0.045,pad=0.1,orientation='horizontal')
        
        plt.savefig("%sNHFLX_SST_Damping_Stdev_CESM1LE_mall_month%02i_lag%i.png" % (figpath,imon+1,ilag+1),
                    dpi=200,bbox_inches='tight',transparent=False)
        
#%% Calculate the pattern correlation

ilag = 0
imon = 0

invar_raw = [damping,rflx,rsst]
invars = [v[:,ilag,imon,:,:].std(0) for v in invar_raw]

# Calculate Pattern correlation
patcorr = np.zeros((2,12,42)) * np.nan # [variable,month,ens]
for imon in tqdm(range(12)):
    for e in range(42):
        invars = [v[e,ilag,imon,:,:] for v in invar_raw]
        # Calculate for cross correlation
        patcorr[0,imon,e] = proc.patterncorr(invars[0],invars[1])
        patcorr[1,imon,e] = proc.patterncorr(invars[0],1/invars[2])
#%% Plot pattern correlations

imon = 0
usebar = False

for imon in tqdm(range(12)):
    fig,ax = plt.subplots(1,1,figsize=(12,4))
    if usebar:
        ax.bar(ens,patcorr[0,imon,:],label="$\lambda_a$ to $CrossCorr$",color="b")
        ax.bar(ens,patcorr[1,imon,:],label="$\lambda_a$ to $AutoCorr^{-1}$",color="r",alpha=0.5)
        ax.set_xlim([0,43])
    
    else:
        
        ax.plot(ens,patcorr[0,imon,:],label="$\lambda_a$:$CrossCorr$",color="b",marker="o")
        ax.plot(ens,patcorr[1,imon,:],label="$\lambda_a$:$AutoCorr^{-1}$",color="r",marker="+")
        ax.set_xlim([1,42])
    
    
    ax.set_ylim([0,1])
    ax.legend(ncol=2)
    ax.set_xticks(ens)
    ax.set_title("%s Pattern Correlation by Ensemble Member" % (months[imon]))
    ax.set_xlabel("Ensemble Member")
    ax.set_ylabel("Pattern Correlation")
    ax.grid(True,ls='dotted')
    
    
    savename = "%sPatternCorr_ByEns_month%i_lag%i.png" % (figpath,imon+1,ilag+1) 
    plt.savefig(savename,dpi=150,bbox_inches='tight')

#%% Plot 3 Panel for select ensemble member and month, with pattern correlation


imon = 2
ilag = 0
e    = 17

invar_raw = [damping,rflx,1/rsst]
msk_raw   = [mall,mflx,msst]

plotnames_fancy = ("Heat Flux Feedback ($\lambda_a$ : $Wm^{-2} \degree C^{-1}$)",
                   "$Q_{net}$-$SST$ Cross-correlation" + "\n Pattern Correlation = %.2f" % (patcorr[0,imon,e]),
                   "($SST$ Autocorrelation)$^{-1}$"+ "\n Pattern Correlation = %.2f" % (patcorr[1,imon,e]))

cints_all = (np.arange(-50,55,5),np.arange(-1,1.1,0.1),np.arange(-2,2.1,0.1))


cmaps     = ('cmo.balance','cmo.balance','cmo.balance')

invars   = [v[e,ilag,imon] for v in invar_raw]
msks_all = [v[e,ilag,imon] for v in msk_raw] 

fig,axs = plt.subplots(1,3,figsize=(16,4),facecolor='white',
                       subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)})


for i in range(3):
    
    ax      = axs[i]
    plotvar = invars[i]
    plotmsk = msks_all[i]
    
    blabel = [0,0,0,1]
    if i == 0: # Set labels for left 
        blabel[0]  = 1
        ax.text(-0.25, 0.35, '%s\n Lag %i\n Ens %i'% (months[imon],ilag+1,e+1), va='bottom', ha='center',
            rotation='horizontal', rotation_mode='anchor',
            transform=ax.transAxes,fontsize=14)
    
    ax = viz.add_coast_grid(ax,bbox=bboxplot,blabels=blabel,fill_color='gray',ignore_error=True)
    
    
    lon1,plotvar1 = proc.lon360to180(lon,plotvar.T)
    _,plotmsk1    = proc.lon360to180(lon,plotmsk.T)
    pcm = ax.contourf(lon1,lat,plotvar1.T,levels=cints_all[i],
                      extend='both',cmap=cmaps[i])
    
    
    
    viz.plot_mask(lon1,lat,plotmsk1,reverse=False,ax=ax,markersize=.35,color='gray')
    
    ax.set_title(plotnames_fancy[i])
    fig.colorbar(pcm,ax=ax,fraction=0.045,pad=0.1,orientation='horizontal')
    #plt.suptitle("J")
    
    plt.savefig("%sNHFLX_SST_Damping_CESM1LE_mall_month%02i_lag%i_ens%02i.png" % (figpath,imon+1,ilag+1,e+1),
                dpi=200,bbox_inches='tight',transparent=False)
    
    
#%% Make Composites of the extreme cases


patcorr_diff = patcorr[0,...] - patcorr[1,...] # [flx - sst]


#%% Make the composites

# Composite top N amount
topN =100

invar_raw = [damping,rflx,1/rsst]
nens,nlag,nmon,nlat,nlon = invar_raw[0].shape
invars = [v[:,ilag,:,:,:].reshape(nens*nmon,nlat,nlon) for v in invar_raw] # Flatten indices

# Transpose to ens x mon to match above, then flatten
patcorr_diff_flat = (patcorr_diff.T).reshape(nens*nmon) 

# Get sorting indices (smallest to largest)
idmin2max = np.argsort(patcorr_diff_flat,axis=0)

# Make composites
invars_max = [v[idmin2max[-topN:],:,:].mean(0) for v in invars]
invars_min = [v[idmin2max[:topN],:,:].mean(0) for v in invars]

# get makers for topN
masktop = np.zeros((nens*nmon))
maskbot = np.zeros((nens*nmon))
masktop[idmin2max[-topN:]] = 1
maskbot[idmin2max[:topN]] = 1
masktop = (masktop.reshape(nens,nmon)).T
maskbot = (maskbot.reshape(nens,nmon)).T

#%% Plot differences, with identified events

fig,ax = plt.subplots(1,1,figsize=(12,4))

pcm = ax.pcolormesh(ens,np.arange(1,13,1),patcorr_diff,shading='nearest',
                    vmin=-.4,vmax=.4,cmap='cmo.balance')
ax.set_aspect('equal')
ax.grid(True,ls='dotted')
cb = fig.colorbar(pcm,ax=ax,fraction=0.025,pad=0.01,orientation='vertical')



viz.plot_mask(ens,np.arange(1,13,1),masktop.T,reverse=True,ax=ax,markersize=10,color='yellow',marker="+")
viz.plot_mask(ens,np.arange(1,13,1),maskbot.T,reverse=True,ax=ax,markersize=10,color='yellow',marker="x")



ax.set_xticks(ens)

ax.set_yticks(np.arange(1,13,1))
ax.set_yticklabels(months)

ax.set_xlabel("Ensemble")
ax.set_title("Pattern Correlation Difference (Cross-correlation - Autocorrelation)",y=1.01)
plt.savefig("%sPattern_Corr_Differences_lag%i_top%i.png" % (figpath,ilag+1,topN),dpi=200,bbox_inches='tight',transparent=False)


#%% Plot composites
cints_all       = (np.arange(-50,55,5),np.arange(-1,1.1,0.1),np.arange(-2,2.1,0.1))
cmaps           = ('cmo.balance','cmo.balance','cmo.balance')
vcat            = ["Top","Bottom"]
plotnames_fancy = ("Heat Flux Feedback ($\lambda_a$ : $Wm^{-2} \degree C^{-1}$)",
                   "$Q_{net}$-$SST$ Cross-correlation",
                   "($SST$ Autocorrelation)$^{-1}$")


fig,axs = plt.subplots(2,3,figsize=(16,8),facecolor='white',constrained_layout=True,
                       subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)})

for i,e in enumerate([invars_max,invars_min]):
    
    for l in range(3):  
        
        ax     = axs[i,l]
        blabel = [0,0,0,0]
        
        if l == 0: # Set labels for left 
            blabel[0]  = 1
            
            ax.text(-0.30, 0.35, '%s %i \n Composite' % (vcat[i],topN), va='bottom', ha='center',
                rotation='horizontal', rotation_mode='anchor',
                transform=ax.transAxes,fontsize=14)
        
        if i == 0:
            ax.set_title(plotnames_fancy[l])
        
        if i == 1: # Set labels for bottom
            blabel[-1] = 1
            
            
        ax  = viz.add_coast_grid(ax,bbox=bboxplot,blabels=blabel,fill_color='gray')
        
        # Plotting Section
        plotvar = e[l]
        lon1,plotvar1 = proc.lon360to180(lon,plotvar.T)
        pcm = ax.contourf(lon1,lat,plotvar1.T,levels=cints_all[l],
                          extend='both',cmap=cmaps[l])
        
        if i == 1:
            fig.colorbar(pcm,ax=ax,fraction=0.055,pad=0.1,orientation='horizontal')
            

plt.savefig("%sNHFLX_Damping_Correlation_CESM1LE_Composites_top%i.png" % (figpath,topN)
            ,dpi=200,bbox_inches='tight',transparent=False)
    

