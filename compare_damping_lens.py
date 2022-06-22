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
    "CANESM2_LENS",
    "CSIRO_MK36_LENS_HTR"
    )

# Name of netCDF
ncnames = (
    "CESM1-LE_NHFLX_Damping_raw.nc",
    "CESM1_rcp85_hfdamping_ensorem1_detrend1_2006to2101_allens.nc",
    "gfdl_esm2m_lens_hfdamping_ensorem1_detrend1_1950to2101_allens.nc",
    "csiro_mk36_lens_hfdamping_ensorem1_detrend1_1850to2101_allens.nc",
    "canesm2_lens_hfdamping_ensorem1_detrend1_1950to2101_allens.nc",
    "csiro_mk36_lens_hfdamping_ensorem1_detrend1_1920to2005_allens.nc"
    )

# Path to ata
datpaths = (
    "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/",
    "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/CESM-HTR-RCP85/",
    "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/lens/",
    "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/lens/",
    "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/lens/",
    "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/lens/"
    ) 

# Degrees of freedom
dofs = (
        82,
        (2101-2006) -4 +1,
        (2101-1950) -4 +1,
        (2101-1850) -4 +1,
        (2101-1950) -4 +1,
        82
        )

yr_rngs = ("1920-2005",
           "2006-2100",
           "1950-2100",
           "1850-2100",
           "1950-2100",
           "1920-2005")

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

covs     = []
autocovs = []

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
        
        # Load covariances
        autocov = ds.autocov.values.transpose(0,2,1,3,4)
        cov     = ds.cov.values.transpose(0,2,1,3,4)
        covs.append(cov)
        autocovs.append(autocov)
    else: # Autocov and Cov not saved for CESM1-LENS
        covs.append(np.zeros(rflx.shape))
        autocovs.append(np.zeros(rflx.shape))
        
    
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
            lon1,plotvar1 = proc.lon360to180(lons[mc],plotvar.T)
            _,plotmsk1    = proc.lon360to180(lons[mc],plotmsk.T)
            
            pcm = ax.contourf(lon1,lats[mc],plotvar1.T,levels=cints,
                              extend='both',cmap=cmap)
            
            viz.plot_mask(lon1,lats[mc],plotmsk1,reverse=False,ax=ax,markersize=.75,color='gray')
            
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
    lon1,plotvar1 = proc.lon360to180(lons[mc],plotvar.T)
    _,plotmsk1    = proc.lon360to180(lons[mc],(plotmsk*limask).T)

    pcm = ax.contourf(lon1,lats[mc],plotvar1.T,levels=cints,cmap='cmo.balance',extend='both')


    # Plot significant points
    viz.plot_mask(lon1,lats[mc],plotmsk1,reverse=False,ax=ax,markersize=.25,color='gray')
    
    # Add Subplot Labels
    ax = viz.label_sp(e+1,ax=ax,labelstyle="Ens%s",alpha=0.75,usenumber=True)
fig.colorbar(pcm,ax=axs.flatten(),orientation='horizontal',fraction=0.025,pad=0.01)
plt.suptitle("%s Lag %i Net Heat Flux Feedback ($Wm^{-2}K^{-1}$)" % (months[im],il+1),y=1.01,fontsize=16)
plt.savefig("%sNHFLX_Damping_%s_month%02i_lag%i.png" % (figpath,mcstr,im+1,il+1),
            dpi=200,bbox_inches='tight',transparent=False)

#%% Make function to determine plotvar_

    



#%% For each of (4) models, plot the ensemble and annual average

[print(d.shape) for d in dampings]

imon = 0
ilag = 0
mode = "STD" # "STD" or "AVG"

plotvar_name = "AUTOCOV" # [HFF, COV, AUTOCOV]

proj = ccrs.PlateCarree(central_longitude=0)


# ======================
if plotvar_name == "HFF":
    invars = dampings
    if mode == "AVG":
        cints = np.arange(-45,47.5,2.5)
        cblbl = "Heat Flux Feedback ($Wm^{-2}$)"
    elif mode == "STD":
        cints = np.arange(0,26,1)
        cblbl = "1$\sigma$ Heat Flux Feedback ($Wm^{-2}$)"
    
elif plotvar_name == "COV":
    invars = covs
    if mode == "AVG":
        cints = np.arange(-30,32,2)
        cblbl = "Cov($Q_{net}'$,$T'$) ($Wm^{-2} \degree C$)"
    elif mode == "STD":
        cints = np.arange(0,10.5,0.5)
        cblbl = "1$\sigma$ Cov($Q_{net}'$,$T'$) ($Wm^{-2} \degree C$)"

elif plotvar_name == "AUTOCOV":
    invars = autocovs
    if mode == "AVG":
        cints = np.arange(-3,3.1,.1)
        cblbl = "Autocov($T'$) ($\degree C^{2}$)"
    elif mode == "STD":
        cints = np.arange(0,1.05,.05)
        cblbl = "1$\sigma$ Autocov($T'$) ($\degree C^{2}$)"
# ======================

#cints = np.arange(-45,47.5,2.5)

fig,axs = plt.subplots(1,5,figsize=(18,4),
                       constrained_layout=True,
                       subplot_kw={'projection':proj})
for mc in range(5):
    
    ax = axs.flatten()[mc]
    blabel=[0,0,0,1]
    if mc == 0:
        blabel[0] = 1
    ax = viz.add_coast_grid(ax,bbox=bboxplot,proj=proj,blabels=blabel,
                            fill_color='gray',ignore_error=True)
    
    
    
    if mode == "AVG":
        cmap  = 'cmo.balance'
        if mc == 2:
            plotvar = invars[mc][1:,ilag,imon,:,:].mean(0) 
            plotmsk = malls[mc][1:,ilag,imon,:,:].mean(0)
        else:
            plotvar = invars[mc][:,ilag,imon,:,:].mean(0) 
            plotmsk = malls[mc][:,ilag,imon,:,:].mean(0)
        
    elif mode == "STD":
        cmap  = "inferno"
        if mc == 2:
            plotvar = invars[mc][1:,ilag,imon,:,:].std(0) 
            plotmsk = malls[mc][1:,ilag,imon,:,:].std(0) 
        else:
            plotvar = invars[mc][:,ilag,imon,:,:].std(0) 
            plotmsk = malls[mc][:,ilag,imon,:,:].std(0) 
    
    # Flip Longitude
    lon1,plotvar1 = proc.lon360to180(lons[mc],plotvar.T)
    _,plotmsk1    = proc.lon360to180(lons[mc],(plotmsk).T)
    
    pcm = ax.contourf(lon1,lats[mc],plotvar1.T,levels=cints,cmap=cmap,extend='both')
    
    ax.set_title("%s (%s)" % (dataset_names[mc],yr_rngs[mc]))
cb = fig.colorbar(pcm,ax=axs.flatten(),fraction=0.015,pad=0.01)
cb.set_label(cblbl)
plt.suptitle("Ensemble %s %s for %s, Lag %i" % (mode,plotvar_name,months[imon],ilag+1),y=0.90)
plt.savefig("%sEnsAvg_%s_mon%02i_lag%i_ens%s.png" % (figpath,plotvar_name,imon+1,ilag+1,mode),dpi=150,bbox_inches='tight')

#%% Plot for Clivar LENS, 30-member cases

mc    = 5
imon  = 0
ilag  = 0

#cints = np.arange(-45,47.5,2.5)

#cblbl = "Heat Flux Feedback ($Wm^{-2}$)" 
    
plotvar_name = "AUTOCOV"
mode         = "AVG"
cmap         = 'cmo.balance'

# ======================
if plotvar_name == "HFF":
    invars = dampings
    if mode == "AVG":
        cints = np.arange(-45,47.5,2.5)
        cblbl = "Heat Flux Feedback ($Wm^{-2}$)"
    elif mode == "STD":
        cints = np.arange(0,26,1)
        cblbl = "1$\sigma$ Heat Flux Feedback ($Wm^{-2}$)"
    
elif plotvar_name == "COV":
    invars = covs
    if mode == "AVG":
        cints = np.arange(-30,32,2)
        cblbl = "Cov($Q_{net}'$,$T'$) ($Wm^{-2} \degree C$)"
    elif mode == "STD":
        cints = np.arange(0,10.5,0.5)
        cblbl = "1$\sigma$ Cov($Q_{net}'$,$T'$) ($Wm^{-2} \degree C$)"

elif plotvar_name == "AUTOCOV":
    invars = autocovs
    if mode == "AVG":
        cints = np.arange(-3,3.1,.1)
        cblbl = "Autocov($T'$) ($\degree C^{2}$)"
    elif mode == "STD":
        cints = np.arange(0,1.05,.05)
        cblbl = "1$\sigma$ Autocov($T'$) ($\degree C^{2}$)"
# ======================

fig,axs = plt.subplots(5,6,figsize=(14,8),
                       constrained_layout=True,
                       subplot_kw={'projection':proj})

for e in tqdm(range(enss[mc][-1])):
    
    ax = axs.flatten()[e]
    
    
    blabel=[0,0,0,0]
    if e%6 == 0:
        blabel[0] = 1
    if e>23:
        blabel[-1] = 1
    
    ax = viz.add_coast_grid(ax,bbox=bboxplot,proj=proj,blabels=blabel,
                            fill_color='gray',ignore_error=True)
    
    
    # Add Subplot Labels
    ax = viz.label_sp(e+1,ax=ax,labelstyle="Ens%s",alpha=0.75,usenumber=True)
    

    plotvar = invars[mc][e,ilag,imon,:,:]
    plotmsk = malls[mc][e,ilag,imon,:,:]

        

    # Flip Longitude
    lon1,plotvar1 = proc.lon360to180(lons[mc],plotvar.T)
    _,plotmsk1    = proc.lon360to180(lons[mc],(plotmsk).T)
    
    pcm = ax.contourf(lon1,lats[mc],plotvar1.T,levels=cints,cmap=cmap,extend='both')
    
    #ax.set_title("%s (%s)" % (dataset_names[mc],yr_rngs[mc]))
cb = fig.colorbar(pcm,ax=axs.flatten(),fraction=0.015,pad=0.01)
cb.set_label(cblbl)
plt.suptitle("%s %s (%s [%s], Lag %i)" % (months[imon],plotvar_name,dataset_names[mc],yr_rngs[mc], ilag+1),y=1.05)
plt.savefig("%s%s_%s_mon%02i_lag%i_ens%s.png" % (figpath,plotvar_name,dataset_names[mc],imon+1,ilag+1,mode),dpi=150,bbox_inches='tight')


#%% Look at CANESM (50 Members)


mc   = 4
imon = 0
ilag = 0
cints = np.arange(-45,47.5,2.5)
cmap  = 'cmo.balance'
cblbl = "Heat Flux Feedback ($Wm^{-2}$)" 


fig,axs = plt.subplots(5,10,figsize=(20,10),
                       constrained_layout=True,
                       subplot_kw={'projection':proj})

for e in tqdm(range(enss[mc][-1])):
    
    ax = axs.flatten()[e]
    
    
    blabel=[0,0,0,0]
    if e%10 == 0:
        blabel[0] = 1
    if e>39:
        blabel[-1] = 1
    
    ax = viz.add_coast_grid(ax,bbox=bboxplot,proj=proj,blabels=blabel,
                            fill_color='gray',ignore_error=True)
    
    
    # Add Subplot Labels
    ax = viz.label_sp(e+1,ax=ax,labelstyle="Ens%s",alpha=0.75,usenumber=True)
    

    plotvar = dampings[mc][e,ilag,imon,:,:]
    plotmsk = malls[mc][e,ilag,imon,:,:]

        

    # Flip Longitude
    lon1,plotvar1 = proc.lon360to180(lons[mc],plotvar.T)
    _,plotmsk1    = proc.lon360to180(lons[mc],(plotmsk).T)
    
    pcm = ax.contourf(lon1,lats[mc],plotvar1.T,levels=cints,cmap=cmap,extend='both')
    
    #ax.set_title("%s (%s)" % (dataset_names[mc],yr_rngs[mc]))
cb = fig.colorbar(pcm,ax=axs.flatten(),fraction=0.015,pad=0.01)
cb.set_label(cblbl)
plt.suptitle("%s Heat Flux Feedback (%s [%s], Lag %i)" % (months[imon],dataset_names[mc],yr_rngs[mc], ilag+1),y=1.01)
plt.savefig("%s%s_HFF_mon%02i_lag%i_ens%s.png" % (figpath,dataset_names[mc],imon+1,ilag+1,mode),dpi=150,bbox_inches='tight')


#%% Explicitly Evaluate Effect of Limiting the Data Period
# On intermember variability

[print(d.shape) for d in dampings]

imon = 0
ilag = 0
mode = "STD" # "STD" or "AVG"

plotids = [0,3,5]

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


fig,axs = plt.subplots(1,3,figsize=(14,4),
                       constrained_layout=True,
                       subplot_kw={'projection':proj})
for a,mc in enumerate(plotids):
    
    ax = axs.flatten()[a]
    blabel=[0,0,0,1]
    if mc == a:
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
plt.suptitle("Ensemble %s Heat Flux Feedback for %s, Lag %i" % (mode,months[imon],ilag+1),y=0.99)
plt.savefig("%sEnsAvg_HFF_mon%02i_lag%i_ens%s_select3panel.png" % (figpath,imon+1,ilag+1,mode),dpi=150,bbox_inches='tight')


#%% Make Scatterplot of Cov vs. Autocov for point with highest intermem. var for each

imon = 0
ilag = 0

lonf = -72+360
latf = 38

data_markers = ("x","+","d","o","v","1")


locfn,loctitle = proc.make_locstring(lonf,latf)

fig,ax = plt.subplots(1,1,constrained_layout=True)

for mc in range(5):
    print(mc)
    #kmax = np.nanargmax(dampings[mc][:,ilag,imon,:,:].std(0).flatten())
    klon,klat = proc.find_latlon(lonf,latf,lons[mc],lats[mc])
    
    plot_acv = autocovs[mc][:,ilag,imon,klat,klon] # [ens]
    plot_cov = covs[mc][:,ilag,imon,klat,klon] # [ens]
    
    
    label= "%s ($\sigma^2_{\lambda}$=%.3f)" % (dataset_names[mc],dampings[mc][:,ilag,imon,klat,klon].std(0))
    
    
    ax.scatter(plot_acv,plot_cov,label=label,alpha=0.6,marker=data_markers[mc])
    
ax.legend(fontsize=8)
ax.grid(True,ls="dotted")
ax.set_xlabel("Autocovariance")
ax.set_ylabel("Covariance")
ax.set_title("%s Covariance vs. Autocovariance @ %s, %s" % (months[imon],loctitle,ilag+1))
savename = "%sScatter_Cov_v_AutoCov_%s_mon%02d_lag%i.png" % (figpath,locfn,imon+1,ilag+1)
plt.savefig(savename,dpi=150,bbox_inches='tight')
    
    #klon,klat = np.unravel_index(kmax,shape=(lats[mc].shape[0],lons[mc].shape[0],))