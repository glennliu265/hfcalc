#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Visualize Heat Flux Feedback for CESM1 Preindustrial Control,
Slab and Full Simulations

Adapted/Copied from 
- calc_HF_func.py
- viz_hfdamping_CESM1LE.ipynb
    

Created on Tue Feb 22 16:11:42 2022
@author: gliu
"""

import sys
sys.path.append("/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/03_Scripts/")
sys.path.append("/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/03_Scripts/stochmod/model/")
from amv import proc,viz
import yo_box as ybx

from scipy.interpolate import interp1d
from scipy.io import loadmat,savemat
from scipy import signal,stats
from tqdm import tqdm

import cartopy.crs as ccrs
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import calendar as cal

import scm
import time
import cmocean
import copy

#%% Data Paths

datpath     = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/"
lipath      = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/01_Data/model_input/"
llpath      = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/01_Data/model_input/"
figpath     = '/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/02_Figures/20220824/'
proc.makedir(figpath)

mconfigs    = ["PIC-SLAB","PIC-FULL"]
dofs        = [898 - 1 - 2 - 2, 1898 - 1 - 2 - 2] 

lags        = [1,2,3]
bboxplot    =  [-80,0,5,62]
mons3       = [viz.return_mon_label(m,nletters=3) for m in np.arange(1,13)]


pubready    = True
#%% Functions

def load_dampraw(mconfig,datpath):
    inpaths = datpath+"CESM-"+mconfig+"-Damping/"
    damping = np.load(inpaths+"NHFLX_Damping_monwin3_lags123_ensorem1_lag1_pcs2_monwin3.npz.npy")
    rflx    = np.load(inpaths+"NHFLX_Crosscorrelation_monwin3_lags123_ensorem1_lag1_pcs2_monwin3.npz.npy")
    rsst    = np.load(inpaths+"SST_Autocorrelation_monwin3_lags123_ensorem1_lag1_pcs2_monwin3.npz.npy")
    return damping,rsst,rflx

#%% # Load the data (taken from calc_HF_func)

# Load some data
#limaskname = "limask180_FULL-HTR.npy" 
# Note this only works locally, need to make equivalent on stormtrack
lon,lat=scm.load_latlon(datpath=llpath,lon360=True)
lon180,_ = scm.load_latlon(datpath=llpath)
limask   = np.load(lipath+"../landicemask_enssum.npy")
limask180    = np.load(lipath+"limask180_FULL-HTR.npy")

dampings   = [] # [mon x lag x lat x lon360]
rssts      = []
rflxs      = []
for i,mcf in enumerate(mconfigs):

    # Load Data
    a,b,c = load_dampraw(mcf,datpath)
    
    # Check sign
    if np.nansum(np.sign(c)) < 0:
        print("WARNING! sst-flx correlation is mostly negative, sign will be flipped")
        c*=-1
    
    dampings.append(a)#*limask[None,None,:,:]) # Apply Landice Mask
    rssts.append(b)
    rflxs.append(c) #Flux must be positive into atm

#%% Calculate critical threshold
p        = 0.05 #0.05
tails    = 2
rhocrits = []

nmon,nlag,nlat,nlon = rssts[0].shape
rmasks    = np.zeros(rssts[0].shape + (2,3,)) * np.nan # [mon x lag x lat x lon x mconfig x masktype]
masktype = ("SST","SST-FLX","BOTH")

masktype_fancy = ("$SST$ Autocorrelation",
                  "$SST-Q_{net}$ Cross-correlation",
                  "Both"
                  )
# mssts = []
# mflxs = []
# malls = []
for i,dof in enumerate(dofs):
    
    # Get critical rho value
    rhocrit = proc.ttest_rho(p,tails,dof)
    rhocrits.append(rhocrit)
    print("The critical rho value for %s is %f" %(mconfigs[i],rhocrit))
    
    # Make the masks
    rmasks[...,i,0] = rssts[i] > rhocrit
    rmasks[...,i,1] = rflxs[i] > rhocrit
    rmasks[...,i,2] = rmasks[...,i,0] * rmasks[...,i,1]
    
    # rsst = rssts[i]
    # rflx = rflxs[i]
    # msst = rsst > rhocrit
    # mflx = rflx > rhocrit
    # mall = msst * mflx 
    # mssts.append(msst)
    # mflxs.append(mflx)
    # malls.append(malls)
#%% Function (modified from .ipynb)

def plot_hfdamping(plotlag,plotmon,lon,lat,damping,msk,ax=None,cints=None,select_id=True):
    
    # damping : [mon x lag x lat x lon]
    # msk     : [mon x lag x lat x lon]
    
    
    if cints is None:
        cints = np.arange(-50,55,5)
    if ax is None:
        ax = plt.gca()
    
    # Select what to plot
    if select_id:
        plotvar = damping[plotmon,plotlag,:,:]
        plotmsk = msk[plotmon,plotlag,:,:]
    else:
        plotvar = damping
        plotmsk = msk
    
    # Flip Longitude
    lon1,plotvar1 = proc.lon360to180(lon,plotvar.T)
    _,plotmsk1    = proc.lon360to180(lon,plotmsk.T)
    
    # Plot the contours
    pcm = ax.contourf(lon1,lat,plotvar1.T,levels=cints,cmap='cmo.balance',extend='both')
    
    # Plot significant points
    viz.plot_mask(lon1,lat,plotmsk1,reverse=False,ax=ax,markersize=.75,color='gray')
    
    return ax,pcm

#%% Visualize the mask

# Visualiz the mask (cross correlation)
rhovalues = [rssts,rflxs]
fig,ax = plt.subplots(1,1,figsize=(10,8),
                       subplot_kw={'projection':ccrs.PlateCarree()})
cints = np.arange(-1,1.1,.1)

plotlag = 1
plotmon = 0
plotmcf = 1
imsk    = 1 # 0 = SST, 1 = Crosscorr

ax.coastlines()
#ax.set_extent(bboxplot)

plotvar  =  rhovalues[imsk][plotmcf][plotlag,plotmon,:,:]
plotmsk  =  rmasks[plotmon,plotlag,:,:,plotmcf,imsk]


# Flip Longitude
lon1,plotvar1 = proc.lon360to180(lon,plotvar.T)
_,plotmsk1    = proc.lon360to180(lon,plotmsk.T)
#_,plotmsk1    = proc.lon360to180(lon,limask.T)

plotvar1 *=limask180
plotmsk1 *= limask180

pcm = ax.contourf(lon1,lat,plotvar1.T,levels=cints,cmap='cmo.balance',extend='both')

# Plot significant points
viz.plot_mask(lon1,lat,plotmsk1,reverse=False,ax=ax,markersize=.75)

# Add Subplot Labels
#ax = viz.label_sp(e+1,ax=ax,labelstyle="Ens%s",alpha=0.75,usenumber=True)
fig.colorbar(pcm,ax=ax,orientation='horizontal',fraction=0.036,pad=.05)
ax.set_title("%s (%s) \n %s; Lag %i " % (masktype_fancy[imsk],mconfigs[plotmcf],mons3[plotmon],lags[plotlag]))

#%% Seasonal Plot of Correlation
plotlag  = 0
plotmcf  = 0
imsk     = 0

outstr     = "lag%s_%s" % (lags[plotlag],mconfigs[plotmcf])
titlestr   = "%s (%s ,Lag: %i)" % (masktype_fancy[imsk],mconfigs[plotmcf],lags[plotlag])


if imsk == 0:
    cints_corr = np.arange(0,1.05,.05)
else:
    cints_corr = np.arange(-1,1.1,.1)

loopmon = np.concatenate([[11,],np.arange(0,11,1)])
fig,axs = plt.subplots(4,3,figsize=(16,16),facecolor='white',constrained_layout=True,
                    subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)})


for i in tqdm(range(12)):
    
    ax     = axs.flatten()[i]
    im     = loopmon[i]
    
    blabel = [0,0,0,0]
    if i%3 == 0:
        blabel[0] = 1
    if i>8:
        blabel[-1] = 1
    
    # Add coastline, Month Title
    ax     = viz.add_coast_grid(ax,bboxplot,blabels=blabel)
    ax.set_title(mons3[im])
    
    # Get plotting variable
    plotvar  =  rhovalues[imsk][plotmcf][im,plotlag,:,:]
    plotmsk  =  rmasks[im,plotlag,:,:,plotmcf,imsk]
    
    
    # Flip Longitude
    lon1,plotvar1 = proc.lon360to180(lon,plotvar.T)
    _,plotmsk1    = proc.lon360to180(lon,plotmsk.T)
    #_,plotmsk1    = proc.lon360to180(lon,limask.T)

    plotvar1 *=limask180
    plotmsk1 *= limask180

    pcm = ax.contourf(lon1,lat,plotvar1.T,levels=cints_corr,cmap='cmo.balance',extend='both')


    # Plot significant points
    viz.plot_mask(lon1,lat,plotmsk1,reverse=False,ax=ax,markersize=.75)
    #ax,pcm = plot_hfdamping(plotlag,im,lon,lat,dampings[plotmcf],msk,ax=ax,cints=cints)
    
cb = fig.colorbar(pcm,ax=axs.flatten(),orientation='horizontal',fraction=0.025)
cb.set_label("Correlation")
plt.suptitle(titlestr,fontsize=16)

plt.savefig("%s%s_corr_%s.png"%(figpath,masktype[imsk],outstr),dpi=150)
#%% || . || . . || . || Try some plots (Damping with selected Mask)


plotlag  = 0
plotmcf  = 0
plotmsk  = 2

outstr   = "lag%s_%s_MASK%s" % (lags[plotlag],mconfigs[plotmcf],masktype[plotmsk])
titlestr = "Heat Flux Damping (%s ,Lag: %i, Mask=%s)" % (mconfigs[plotmcf],lags[plotlag],masktype[plotmsk])


cints = np.arange(-50,55,5)

loopmon = np.concatenate([[11,],np.arange(0,11,1)])
fig,axs = plt.subplots(4,3,figsize=(16,16),facecolor='white',constrained_layout=True,
                    subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)})


for i in tqdm(range(12)):
    
    ax     = axs.flatten()[i]
    im     = loopmon[i]
    
    blabel = [0,0,0,0]
    if i%3 == 0:
        blabel[0] = 1
    if i>8:
        blabel[-1] = 1
    
    ax     = viz.add_coast_grid(ax,bboxplot,blabels=blabel)
    ax.set_title(mons3[im])
    
    
    msk    = rmasks[:,:,:,:,plotmcf,plotmsk]
    ax,pcm = plot_hfdamping(plotlag,im,lon,lat,dampings[plotmcf]*-1,msk,ax=ax,cints=cints)
    

cb = fig.colorbar(pcm,ax=axs.flatten(),orientation='horizontal',fraction=0.025)
cb.set_label("Heat Flux Damping ($W m^{-2} K^{-1}$)")
plt.suptitle(titlestr,fontsize=16)
plt.savefig("%sDamping_%s.png"%(figpath,outstr),dpi=150)
#ax = viz.label_sp(0,ax=ax,labelstyle="%s"%mons3[im],alpha=0.75,usenumber=True)


#%% Plot comparison of HFF with each lag (Ann Avg)


plotlag  = 1
plotmcf  = 1
plotmsk  = 2
im       = 0
avgflag  = 'avg' # Can also use 'avg','max','min'
keep_sign = False

outstr   = "lagdiff_%s_MASK%s" % (mconfigs[plotmcf],masktype[plotmsk])
titlestr = "Heat Flux Damping (%s ,Lag: %i, Mask=%s)" % (mconfigs[plotmcf],lags[plotlag],masktype[plotmsk])


tautitles = (r"$\tau_{1}$",r"$\tau_{2}$ - $\tau_{1}$",r"$\tau_{3}$ - $\tau_{1}$")

cints    = np.arange(-50,55,5)
cintdiff = np.arange(-10,11,1)

loopmon = np.concatenate([[11,],np.arange(0,11,1)])

fig,axs = plt.subplots(1,3,figsize=(16,6),facecolor='white',constrained_layout=True,
                    subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)})


for i in range(3):
    
    ax     = axs.flatten()[i]
    blabel = [0,0,0,1]
    if i == 0:
        blabel[0] = 1    
    ax     = viz.add_coast_grid(ax,bboxplot,blabels=blabel,fill_color="gray")
    ax.set_title(tautitles[i])
    
    if avgflag:
        im = np.arange(0,12,1)
    # if not isinstance(im,int):
    #     avgflag = copy.deepcopy(im) # Copy Operation
    #     im      = np.arange(0,12,1)
    # else:
    #     avgflag = False
        

    if i == 0:
        plothf = dampings[plotmcf][im,0,:,:] # [mon x lag x lat x lon]
        msk    = rmasks[im,0,:,:,plotmcf,plotmsk]
        cint_in = cints
    elif i == 1:
        plothf = dampings[plotmcf][im,1,:,:] - dampings[plotmcf][im,0,:,:]
        msk    = rmasks[im,0,:,:,plotmcf,plotmsk] * rmasks[im,1,:,:,plotmcf,plotmsk]
        cint_in = cintdiff
    elif i == 2:
        plothf = dampings[plotmcf][im,2,:,:] - dampings[plotmcf][im,0,:,:]
        msk    = rmasks[im,0,:,:,plotmcf,plotmsk] * rmasks[im,2,:,:,plotmcf,plotmsk]
        cint_in = cintdiff
        
    if avgflag:
        if avgflag == "avg":
            plothf = plothf.mean(0)
        elif avgflag == "max":
            plothf  = proc.maxabs(plothf,axis=0,keep_sign=keep_sign)
        elif avgflag == "min":
            plothf  = proc.minabs(plothf,axis=0,keep_sign=keep_sign)
        msk = msk.prod(0)
            
    ax,pcm = plot_hfdamping(plotlag,im,lon,lat,plothf*-1*limask,msk*limask,ax=ax,cints=cint_in,select_id=False)
    
    if i == 0:
        cb = fig.colorbar(pcm,ax=axs[0],orientation='horizontal',fraction=0.045)
    elif i == 2:
        cb = fig.colorbar(pcm,ax=axs[1:].flatten(),orientation='horizontal',fraction=0.045)

if avgflag:
    savename= r"%s Heat Flux Differences ($Wm^{-2}\degree C^{-1}$)" % (avgflag.upper())
else:    
    savename= r"%s Heat Flux Differences ($Wm^{-2}\degree C^{-1}$)" % (mons3[im])
plt.suptitle(savename,y=.85,fontsize=14)
plt.savefig("%sDamping_Differences_%s.png"%(figpath,outstr),dpi=150,bbox_inches='tight')

# for i in tqdm(range(12)):
    
#     ax     = axs.flatten()[i]
#     im     = loopmon[i]
    
#     blabel = [0,0,0,0]
#     if i%3 == 0:
#         blabel[0] = 1
#     if i>8:
#         blabel[-1] = 1
    
#     ax     = viz.add_coast_grid(ax,bboxplot,blabels=blabel)
#     ax.set_title(mons3[im])
    
    
#     msk    = rmasks[:,:,:,:,plotmcf,plotmsk]
#     ax,pcm = plot_hfdamping(plotlag,im,lon,lat,dampings[plotmcf],msk,ax=ax,cints=cints)
    

# cb = fig.colorbar(pcm,ax=axs.flatten(),orientation='horizontal',fraction=0.025)
# cb.set_label("Heat Flux Damping ($W m^{-2} K^{-1}$)")
# plt.suptitle(titlestr,fontsize=16)
#plt.savefig("%sDamping_%s.png"%(figpath,outstr),dpi=150)
# #ax = viz.label_sp(0,ax=ax,labelstyle="%s"%mons3[im],alpha=0.75,usenumber=True)

#%% Calculate seasonal averages



sids = ([11,0,1],[2,3,4],[5,6,7],[8,9,10])


ilag  = 0
savgs = []
smasks_mcf = []
for imcf in range(2):
    
    # Compute seasonal averages
    savg,monstrs=proc.calc_savg(dampings[imcf][:,ilag,:,:],return_str=True,axis=0)
    savgs.append(savg)
    
    smasks_m = []
    # Compute seasonally average masks by model
    for s in range(4):
        sid = sids[s]
        smask = np.prod(rmasks[sid,ilag,...,imcf,2],(0)) #[mon x lat x lon x mconfig]
        smasks_m.append(smask)
    smasks_mcf.append(smasks_m)
    
    
    


smasks = []
for s in range(4):
    sid = sids[s]
    
    # Make seasonal masks, where rmasks = [mon x lag x lat x lon x mconfig x masktype]
    # Multiple along month and mconfig option
    smask = np.prod(rmasks[sid,ilag,...,2],(0,-1)) #[mon x lat x lon x mconfig]
    smasks.append(smask)

# ------------------------
#%% Plot Seasonal Averages
# ------------------------
"""
# Stochastic model Draft 5.0 Figure (Copied from viz_inputs_point.py)
# Updated For Revision 1 08242022, New Bounding Box
"""

notitle   = True
plot_mask = True


proj = ccrs.PlateCarree(central_longitude=0)
fig,axs = plt.subplots(2,2,figsize=(9,7.5),constrained_layout=True,
                      subplot_kw={'projection':proj})

sp_id       = 0
bboxplot    = [-80,0,9,65]
cints       = np.arange(-24,26,2)
snamesl     = ('Winter (DJF)','Spring (MAM)','Summer (JJA)','Fall (SON)')

for s,ax in enumerate(axs.flatten()):
    
    
    blabel = [0,0,0,0]
    if s%2 == 0:
        blabel[0] = 1
    if s>1:
        blabel[-1]=1    
    
    plotvar = (savgs[1][s]-savgs[0][s]) * -1 * limask
    
    print(s)
    pcm = ax.contourf(lon,lat,plotvar,levels=cints,cmap='cmo.balance',extend='both')
    #cl  = ax.contour(lon,lat,plotvar,levels=cints,colors='k',linewidths=0.5)
    #ax.clabel(cl,cints[::2],fontsize=8)
    
    if plot_mask:
        viz.plot_mask(lon,lat,(smasks[s]*limask).T,reverse=False,
                      ax=ax,markersize=1,color='gray',proj=proj,geoaxes=True)
        
    ax = viz.add_coast_grid(ax=ax,bbox=bboxplot,fill_color='gray',blabels=blabel,proj=proj)
    ax.set_title(snamesl[s],fontweight='bold')
    #fig.colorbar(pcm,ax=ax)
    ax = viz.label_sp(sp_id,ax=ax,fontsize=14,fig=fig,labelstyle="(%s)",case='lower',alpha=.75)
    sp_id += 1
    
cb = fig.colorbar(pcm,ax=axs.flatten(),orientation='horizontal',fraction=0.030,pad=0.02)
cb.set_label("Atmospheric Heat Flux Feedback ($Wm^{-2}K^{-1}$)")
if notitle is False:
    plt.suptitle("Seasonal Mean Mixed Layer Depth Differences in meters \n (CESM1 - WOA 1994)",fontsize=14,y=.94)

if pubready:
    plt.savefig("%sFig12_CESM_HFF_Differences.png" %(figpath),dpi=900,bbox_inches='tight')
else:
    plt.savefig("%sHFLX_Differences_FULL-SLAB_Savg.png" %(figpath),dpi=150,bbox_inches='tight')

# -----------------------------------------------------------
#%% Plot HFF (SLAB vs FULL), seasonal averages, over a REGION
# -----------------------------------------------------------

# Copied from viz_inputs_point.py

viz_bbox  = [-80,-10,0,25] # Tropics
plot_mask = True

snamesl = ('Winter (DJF)','Spring (MAM)','Summer (JJA)','Fall (SON)')
fig,axs = plt.subplots(2,4,figsize=(18,6),constrained_layout=True,
                      subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)})

sp_id = 0
for imodel in range(2):
    for s in range(4):
        
        ax = axs[imodel,s]
        blabel = [0,0,0,0]
        if s == 0:
            blabel[0] = 1
            
        if imodel==1:
            blabel[-1]=1    
        
        
        plotvar = savgs[imodel][s] * -1 * limask
        print(s)
        pcm = ax.contourf(lon,lat,plotvar,levels=cints,cmap='cmo.balance',extend='both')
        cl  = ax.contour(lon,lat,plotvar,levels=cints,colors='k',linewidths=0.5)
        ax.clabel(cl,cints[::2],fontsize=8)
        
        
        if plot_mask:
            viz.plot_mask(lon,lat,(smasks_mcf[imodel][s]*limask).T,reverse=False,ax=ax,markersize=1,color='gray')
        ax = viz.add_coast_grid(ax=ax,bbox=viz_bbox,fill_color='gray',blabels=blabel)
        if imodel == 0:
            ax.set_title(snamesl[s],fontweight='bold')
        
        # Label Subplot
        ax = viz.label_sp(sp_id,ax=ax,fontsize=14,fig=fig,labelstyle="(%s)",case='lower',alpha=.75)
        sp_id += 1
        
        # Add Extra Label
        if s ==0:
            txt = viz.add_ylabel(mconfigs[imodel],ax=ax,x=-0.15,y=0.5)
        

cb = fig.colorbar(pcm,ax=axs.flatten(),orientation='horizontal',fraction=0.030,pad=0.02)
cb.set_label("Atmospheric Heat Flux Feedback ($Wm^{-2}K^{-1}$)")

# if notitle is False:
#     plt.suptitle("Seasonal Mean Mixed Layer Depth Differences in meters \n (CESM1 - WOA 1994)",fontsize=14,y=.94)
plt.savefig("%sHFLX_values_FULL-SLAB_Savg_regional.png" %(figpath),dpi=150,bbox_inches='tight')
# ------------------------------
#%% Check Conditions at 1 point
# ------------------------------
lonf      = -30+360
latf      = 50
klon,klat = proc.find_latlon(lonf,latf,lon,lat)
flip      = -1
locstring  = "Lon %.2f ; Lat %.2f" % (lon[klon]-360,lat[klat])
locfstring = "Lon_%i_Lat%i" % (lon[klon]-360,lat[klat])
plotmean    = True

# Plot heat flux feedback
fig,axs = plt.subplots(2,1,figsize=(8,8),sharey=True)
dampingpts = []
for imcf in range(2):
    ax = axs.flatten()[imcf]
    dampingpts.append(dampings[imcf][:,:,klat,klon])
    for ilag in range(3):
        ax.plot(mons3,dampingpts[imcf][:,ilag]*-1,label="Lag %i"% (ilag+1),lw=2.5)
    if plotmean:
        ax.plot(mons3,dampingpts[imcf].mean(-1)*-1,label="Lag Avg.",color='k',lw=1,ls='dashed')
    #ax.axhline(rhocrits[imcf],ls='dotted',color='gray')
    
    ax.set_title(mconfigs[imcf])
    ax.grid(True,ls='dotted')
    ax.set_xlim([0,11])
    ax.legend()
    ax.set_ylabel("Heat Flux Damping ($Wm^{-2}K^{-1}$)")
    
plt.suptitle("Heat Flux Feedback Estimates @ %s" % (locstring),y=0.95,fontsize=16)
plt.savefig("%sDamping_LagAvg_%s_plotmean%i.png"%(figpath,locfstring,plotmean),dpi=100)

#%% Plot the correlation
fig,axs = plt.subplots(2,2,figsize=(16,8),sharex=True)
rvars = [rssts,rflxs]
corrnames = [""]
for v in range(2):
    for imcf in range(2):

        ax = axs[v,imcf]
        for ilag in range(3):
            ax.plot(mons3,rvars[v][imcf][:,ilag,klat,klon],label="Lag %i"% (ilag+1),lw=2.5)
        ax.axhline(rhocrits[imcf],color='k',ls='dashed')
        ax.set_title(mconfigs[imcf])
        ax.grid(True,ls='dotted')
        ax.set_xlim([0,11])
        ax.legend()
        if imcf==0:
            ax.set_ylabel(masktype_fancy[v])
        
    plt.suptitle("Correlation for HFF Estimation @ %s" % (locstring),y=0.95,fontsize=16)
    plt.savefig("%sCorrelationg_LagAvg_%s.png"%(figpath,locfstring),dpi=100)

#%% Compute old damping with new damping


dampingold = dampings[1] # [mon lag lat lon]

newpath    = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/CESM1_FULL_PIC_hfdamping_ensorem1.nc"
dampingnew = xr.open_dataset(newpath).nhflx_damping.values 


# Plot the differences
cints   = np.arange(-5,5.025,0.025)
loopmon = np.concatenate([[11,],np.arange(0,11,1)])
fig,axs = plt.subplots(4,3,figsize=(16,16),facecolor='white',constrained_layout=True,
                    subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)})

for i in tqdm(range(12)):
    
    ax     = axs.flatten()[i]
    im     = loopmon[i]
    
    blabel = [0,0,0,0]
    if i%3 == 0:
        blabel[0] = 1
    if i>8:
        blabel[-1] = 1
    
    ax     = viz.add_coast_grid(ax,bboxplot,blabels=blabel)
    ax.set_title(mons3[im])
    
    plotvar = dampingnew[im,0,:,:] - dampingold[im,0,:,:]
    lon1,plotvar1 = proc.lon360to180(lon,plotvar.T[:,:,None])
    pcm = ax.contourf(lon1,lat,plotvar1.squeeze().T,levels=cints,
                      cmap='cmo.balance',extend='both')
    
    cl  = ax.contour(lon1,lat,plotvar1.squeeze().T,levels=cints,
                      colors="k",linewidths=0.75)
    ax.clabel(cl,cints[::2])
    
    #msk    = rmasks[:,:,:,:,plotmcf,plotmsk]
    #ax,pcm = plot_hfdamping(plotlag,im,lon,lat,dampings[plotmcf]*-1,msk,ax=ax,cints=cints)
    

cb = fig.colorbar(pcm,ax=axs.flatten(),orientation='horizontal',fraction=0.025)
cb.set_label("Heat Flux Damping Differences ($W m^{-2} K^{-1}$)")
plt.suptitle("New Damping - Old Damping",fontsize=16)
plt.savefig("%sDamping_Differences_%s.png"%(figpath,outstr),dpi=150)


hfdiff = dampingnew-dampingold



#%%
