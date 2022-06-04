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
figpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/02_Figures/20220609/"
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

vlims     = [-.5,1]
cints_all = [np.arange(-.5,.55,0.05),np.arange(-1,1.1,.1)]


invars    = [rflx,rsst]

rnames       = ["crosscorr","autocorr"]
rnames_fancy = ("NHFLX-SST Cross-correlation","SST Autocorrelation")
mskin        = [mflx,msst]

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
    

#%% Load and make comparisons with calculated timescale


# Load Data
npz = True

if npz:
    # Open and Load
    ncname ="/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/03_reemergence/01_Data/proc/HTR-FULL_SST_autocorrelation_thres0_lag00to60.npz"
    ld = np.load(ncname,allow_pickle=True)
    counts     = ld['class_count']
    acs        = ld['acs']
    cfs        = ld['cfs']
    thresholds = ld['lon']
    lons       = ld['lon']
    lats       = ld['lat']
    lags       = ld['lags']
    threshlabs = ld['threslabs']
    
    # Calculate T2
    t2         = 1 + 2 * np.trapz(acs**2,x=lags,axis=-1) # [lon x lat x ens x mon x thres]
    nlon,nlat,nens,nmon,nthres = t2.shape
    
    # transpose to order # [ens x thres x mon x lat x lon]
    t2 = t2.transpose(2,4,3,1,0)
    
else:
    ncname ="/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/03_reemergence/01_Data/proc/HTR-FULL_SST_autocorrelation_thres0_lag00to60.nc"
    ds     = xr.open_dataset(ncname)
    dssum  = ds.isel(thres=2).sum('lag')
    lons   = dssum.lon.values
    lats   = dssum.lat.values
    sstac  = dssum.SST.values # [ens x mon x lats x lons]
    nens,nmon,nlatr,nlonr=sstac.shape






#%% Prepare to Cut damping variables to the region
bboxreg = [lons[0],lons[-1],lats[0],lats[-1]]

# Reshape damping (also flip longitude, and cut region)
nens,nlag,nmon,nlat,nlon = damping.shape
damping_rs      = damping.reshape(np.array(damping.shape[:-2]).prod(),nlat,nlon) # Combine dims
damping_rs      = damping_rs.transpose(2,1,0) # Flip to [lon x lat x otherdim]

lon1,damping_rs180 = proc.lon360to180(lon,damping_rs) # Flip longitude
dampingr,lonr,latr, = proc.sel_region(damping_rs180,lon1,lat,bboxreg) # Cut Region

dampingr        = dampingr.transpose(2,1,0) # Flip to [otherdim x lat x lon]
dampingr        = dampingr.reshape(damping.shape[:-2]+(len(latr),len(lonr))) # Uncombine dims


# Now compute the pattern correlation
dampingr_in = dampingr[:,0,:,:,:] # Select first lag # [ens x mon x lat x lon]
patcorr_t2 = np.zeros((nens,nmon)) * np.nan
for e in tqdm(range(nens)):
    for m in range(12):
        patcorr_t2[e,m] = proc.patterncorr(dampingr_in[e,m,:,:],sstac[e,m,:,:])

# -------------------------------------
#%% Plot Intermember Variability (Stdev)

ithres = -1
ilag   = 0
imon   = 11
annmean = True

plotpoint = [-72,39] # Set to None to not plot a point, otherwise lat,lon


if annmean is False:
    monstr =  str(imon+1)
else:
    monstr = "Ann. Mean"

bboxplot  = [-80,0,5,65]
if annmean:
    
    t2_std  = np.nanstd(t2[:,ithres,:,:,:],0).mean(0) # [Lat x Lon]
    lbd_std = np.nanstd(dampingr[:,ilag,:,:,:],0).mean(0) 
else:
    t2_std  = np.nanstd(t2[:,ithres,imon,:,:],0) # [Lat x Lon]
    lbd_std = np.nanstd(dampingr[:,ilag,imon,:,:],0) 
    
levels  = np.arange(0,7.5,.5)


lbd_levels = np.arange(0,33,3)

fig,axs = plt.subplots(1,2,subplot_kw={'projection':ccrs.PlateCarree()},
                      constrained_layout=True,figsize=(12,4))

# -------------------- (T2 Stedev)
ax     = axs[0]
ax     = viz.add_coast_grid(ax,bbox=bboxplot,fill_color='gray')
pcm    = ax.contourf(lons,lats,t2_std,levels=levels,cmap='cmo.thermal',extend='both')
cl     = ax.contour(lons,lats,t2_std,levels=levels,colors="w",linewidths=0.35)
ax.clabel(cl,levels=levels[::4])
cb = fig.colorbar(pcm,ax=ax)
cb.set_label("$\sigma_{T_2}$ (Months)")
ax.set_title("$T_2$ for SST")

# -------------------- (Lbd Stdv)
ax     = axs[1]
ax     = viz.add_coast_grid(ax,bbox=bboxplot,fill_color='gray')
pcm1    = ax.contourf(lons,lats,lbd_std,levels=lbd_levels,cmap='cmo.thermal',extend='both')
cl     = ax.contour(lons,lats,lbd_std,levels=lbd_levels,colors="w",linewidths=0.35)
ax.clabel(cl,levels=lbd_levels[::4])
cb1 = fig.colorbar(pcm1,ax=ax)
cb1.set_label("$\sigma_{\lambda_a}$ ($Wm^{-2}\degree C^{-1}$)")
ax.set_title("$\lambda_a$")


if plotpoint is not None:
    for ax in axs:
        ax.plot(plotpoint[0],plotpoint[1],marker="x",color='k',markersize=30)

plt.suptitle("Intermember Variability in CESM1 Historical (42-members), Lag %i, Month: %s"% (ilag+1,monstr))
savename = "%sIntermemVar_T2_SST_CESM1LENS_HTR_lag%i_imon%s_ithres%i.png" % (figpath,ilag+1,monstr,ithres)
if annmean:
    savename = proc.addstrtoext(savename,"annmean")
plt.savefig(savename,dpi=150,bbox_inches='tight')

# ------------------------------------------------
#%% Do the same, but for cross and autocorrelation
# ------------------------------------------------
ilag   = 0
imon   = 11
annmean = True

plotpoint = [] # Set to None to not plot a point


if annmean is False:
    monstr =  str(imon+1)
else:
    monstr = "Ann. Mean"

corrlevels = np.arange(0,0.105,0.005)

fig,axs = plt.subplots(1,2,subplot_kw={'projection':ccrs.PlateCarree()},
                      constrained_layout=True,figsize=(12,4))

# -------------------- (T2 Stedev)
if annmean:
    plotvar = np.nanstd(rflx[:,ilag,:,:,:],0).mean(0)
else:
    plotvar = np.nanstd(rflx[:,ilag,imon,:,:],0)

ax     = axs[0]
ax     = viz.add_coast_grid(ax,bbox=bboxplot,fill_color='gray')
pcm    = ax.contourf(lon,lat,plotvar,levels=corrlevels,cmap='cmo.thermal',extend='both')
cl     = ax.contour(lon,lat,plotvar,levels=corrlevels,colors="w",linewidths=0.35)
ax.clabel(cl,levels=corrlevels[::4])
cb = fig.colorbar(pcm,ax=ax)
#cb.set_label("$\sigma_{T_2}$ (Months)")
ax.set_title("SST-FLX Cross Correlation")

# -------------------- (Lbd Stdv)
if annmean:
    plotvar = np.nanstd(rsst[:,ilag,:,:,:],0).mean(0)
else:
    plotvar = np.nanstd(rsst[:,ilag,imon,:,:],0)
ax     = axs[1]
ax     = viz.add_coast_grid(ax,bbox=bboxplot,fill_color='gray')
pcm1    = ax.contourf(lon,lat,plotvar,levels=corrlevels,cmap='cmo.thermal',extend='both')
cl     = ax.contour(lon,lat,plotvar,levels=corrlevels,colors="w",linewidths=0.35)
ax.clabel(cl,levels=corrlevels[::4])
cb1 = fig.colorbar(pcm1,ax=ax)
# cb1.set_label("$\sigma_{\lambda_a}$ ($Wm^{-2}\degree C^{-1}$)")
ax.set_title("SST Autocorrelation")



plt.suptitle("Intermember Variability in CESM1 Historical (42-members), Lag %i, Month: %s"% (ilag+1,monstr))
savename = "%sIntermemVar_correlations_CESM1LENS_HTR_lag%i_mon%s_thres%i.png" % (figpath,ilag+1,monstr,ithres)
if annmean:
    savename = proc.addstrtoext(savename,"annmean")
plt.savefig(savename,dpi=150,bbox_inches='tight')


#%% Examamine actual values at a point

klon,klat = proc.find_latlon(plotpoint[0],plotpoint[1],lons,lats)

locfn,locstr = proc.make_locstring(plotpoint[0],plotpoint[1])

fig,axs    = plt.subplots(3,1,figsize=(12,8))


# Get variables to plot
selvars  = [dampingr,rflx,rsst]
plotvars = [v[:,ilag,:,klat,klon].mean(-1) for v in selvars]

# Set names, colors
bcols = ["cornflowerblue","goldenrod","orchid"]
bnams = ["Heat FLux Feedback","Cross-Correlation","Autocorrelation"]
ylims = ([0,110],[.15,.30],[.75,.88])



for v in range(3):
    ax = axs[v]
    plotvar = plotvars[v]
    
    ax.bar(np.arange(1,43),plotvar,color=bcols[v],alpha=0.75)
    ax.axhline(np.mean(plotvar),ls='solid',color="k",lw=0.55)
    
    
    ax.set_xticks(np.arange(1,43))
    ax.set_ylabel(bnams[v])
    ax.set_ylim(ylims[v])
    ax.grid(True,ls='dotted')

plt.suptitle("CESM1-LENS, Values at %s"%(locstr),y=.93,fontsize=14)


# # Plot the Damping value
# ax = axs[0]
# lbd_bp    = dampingr[:,ilag,:,klat,klon].mean(-1)
# ax.bar(np.arange(1,43),lbd_bp)
# ax.set_xticks(np.arange(1,43))
# ax.set_ylabel("Heat Flux Feedback")

# # Plot Cross Correlation
# ax = axs[1]
# plotvar    = rflx[:,ilag,:,klat,klon].mean(-1)
# ax.bar(np.arange(1,43),plotvar,color='r')
# ax.set_xticks(np.arange(1,43))
# ax.set_ylabel("Cross-Correlation")
# ax.set_ylim(.15,.30)

# # Plot Cross Correlation
# ax = axs[2]
# plotvar    = rsst[:,ilag,:,klat,klon].mean(-1)
# ax.bar(np.arange(1,43),plotvar,color='violet')
# ax.set_xticks(np.arange(1,43))
# ax.set_ylabel("Autocorrelation")
# ax.set_ylim(.75,.88)

plt.savefig("%sHFF_values_at_%s_CESM1LENS.png"%(figpath,locfn),dpi=150,bbox_inches='tight')
#%% Visualize pcolor of pattern correlation for each month

fig,ax = plt.subplots(1,1,figsize=(12,4))

pcm = ax.pcolormesh(ens,np.arange(1,13,1),patcorr_t2.T,shading='nearest',
                    vmin=-.8,vmax=.8,cmap='cmo.balance')

ax.set_aspect('equal')
ax.grid(True,ls='dotted')
cb = fig.colorbar(pcm,ax=ax,fraction=0.025,pad=0.01,orientation='vertical')

ax.set_xticks(ens)

ax.set_yticks(np.arange(1,13,1))
ax.set_yticklabels(months)

ax.set_xlabel("Ensemble")
ax.set_ylabel("Month (HFF) or Reference Month for Lag ($T_2$)")
ax.set_title("Pattern Correlation (Heat Flux Feedback and SST $T_2$ Timescale)",y=1.01)
plt.savefig("%sPattern_Corr_T2_HFF_lag%i_mon%i_thres%i.png" % (figpath,ilag,imon,ithres),dpi=200,bbox_inches='tight',transparent=False)


#%% Make some scatterplots

# Overall scatterplot (Ok, it's not making much sense....)
fig,ax = plt.subplots(1,1)
sc = ax.scatter(dampingr_in.flatten(),sstac.flatten(),
                alpha=0.3)
ax.set_ylim([0,25])
ax.set_xlim([-50,50])

# -----------------------------------------------------------------------------
#%% Observe behavior over specific regions

# From SM Stylesheet
bbox_SP     = [-60,-15,40,65]
bbox_ST     = [-80,-10,20,40]
bbox_TR     = [-75,-15,10,20]
bbox_NA     = [-80,0 ,0,65]
bbox_NA_new = [-80,0,10,65]
bbox_ST_w   = [-80,-40,20,40]
bbox_ST_e   = [-40,-10,20,40]
regions     = ("SPG","STG","TRO","NAT","NNAT","STGe","STGw")        # Region Names
bboxes      = (bbox_SP,bbox_ST,bbox_TR,bbox_NA,bbox_NA_new,bbox_ST_e,bbox_ST_w) # Bounding Boxes

bbcol       = ["Blue","purple","Yellow","Black","Black","magenta","red"] # >> Need to Update

# Preallocate
nreg    = len(bboxes)
t2_reg  = np.zeros((nens,nmon,nreg))
hff_reg = t2_reg.copy()

# Reshape to lon x lat x otherdims
dampin = dampingr_in.reshape(nens*12,nlatr,nlonr).transpose(2,1,0)
t2in   = sstac.reshape(nens*12,nlatr,nlonr).transpose(2,1,0)


dampreg = []
t2reg   = []

for b,bbox in tqdm(enumerate(bboxes)):
    
    # Save regional average
    dampavg = proc.sel_region(dampin,lonr,latr,bbox,reg_avg=1,awgt=1)
    hff_reg[:,:,b] = dampavg.reshape(nens,12)
    
    # Save the regional values
    rdamp,_,_ = proc.sel_region(dampin,lonr,latr,bbox)
    dampreg.append(rdamp.reshape(rdamp.shape[:2]+(nens,12)))
    
    # Save avg (T2)
    t2avg = proc.sel_region(t2in,lonr,latr,bbox,reg_avg=1,awgt=1)
    t2_reg[:,:,b] = t2avg.reshape(nens,12)
    
    # Save values (T2)
    rt2,_,_ = proc.sel_region(t2in,lonr,latr,bbox)
    t2reg.append(rt2.reshape(rt2.shape[:2]+(nens,12)))
    
    
#%% Repeat pattern correlation for each region

patcorr_t2_reg = np.zeros((nens,nmon,nreg)) * np.nan # [ens x month x region]
for r in tqdm(range(nreg)):
    for e in range(nens):
        for m in range(12):
            dmp = dampreg[r][:,:,e,m]
            t2  = t2reg[r][:,:,e,m]

            patcorr_t2_reg[e,m,r] = proc.patterncorr(dmp,t2)
        
#%% Pcolor plot for each region (Ensemble vs. Reference Month)

for r in range(nreg):
    
    fig,ax = plt.subplots(1,1,figsize=(12,4))
    
    pcm = ax.pcolormesh(ens,np.arange(1,13,1),patcorr_t2_reg[:,:,r].T,shading='nearest',
                        vmin=-.8,vmax=.8,cmap='cmo.balance')
    
    ax.set_aspect('equal')
    ax.grid(True,ls='dotted')
    cb = fig.colorbar(pcm,ax=ax,fraction=0.025,pad=0.01,orientation='vertical')

    ax.set_xticks(ens)
    ax.set_yticks(np.arange(1,13,1))
    ax.set_yticklabels(months)

    ax.set_xlabel("Ensemble")
    ax.set_ylabel("Month (HFF) or Reference Month for Lag ($T_2$)")
    ax.set_title("%s Pattern Correlation (Heat Flux Feedback and SST $T_2$ Timescale)" % (regions[r]),y=1.01)
    plt.savefig("%sPattern_Corr_T2_HFF_lag%i_region%s.png" % (figpath,1,regions[r]),dpi=200,bbox_inches='tight',transparent=False)


#%% Calculate topN composites for each region

topN      = 100
rid_sel   = [0,2,4,5,6]

invar_raw = [dampingr_in,sstac] # (42, 12, 69, 65)
invars    = [v.reshape(nens*nmon,nlatr,nlonr) for v in invar_raw] # Flatten indices : (nens*nmon x lat x lon)
invar_names = ["Damping","T2"]

t2_composites = np.zeros((nlatr,nlonr,nreg,2)) *np.nan # [lat x lon x region x type (0=min,1=max)]
hf_composites = t2_composites.copy()
pcolmarkers  = np.zeros((nmon,nens,nreg,2)) * np.nan


kminmax = np.zeros((nreg,topN,2)) * np.nan #[region x rank x type (0=min,1=max)] 

for r in tqdm(range(nreg)):
    
    # Get patcorr for region and flatten
    patcorr_in        = patcorr_t2_reg[...,r]
    patcorr_diff_flat = (patcorr_in).reshape(nens*nmon) 
    
    # Get sorting indices (smallest to largest)
    idmin2max = np.argsort(patcorr_diff_flat,axis=0)
    
    # Make composites
    invars_max = [v[idmin2max[-topN:],:,:].mean(0) for v in invars]
    invars_min = [v[idmin2max[:topN],:,:].mean(0) for v in invars]
    
    # Save indices of the ens/mon
    kminmax[r,:,1] = idmin2max[-topN:]
    kminmax[r,:,0] = idmin2max[:topN]
    
    # Read out to array
    t2_composites[:,:,r,0] = invars_min[1].copy()
    t2_composites[:,:,r,1] = invars_max[1].copy()
    hf_composites[:,:,r,0] = invars_min[0].copy()
    hf_composites[:,:,r,1] = invars_max[0].copy()
    
    # get makers for topN
    masktop = np.zeros((nens*nmon))
    maskbot = np.zeros((nens*nmon))
    masktop[idmin2max[-topN:]] = 1
    maskbot[idmin2max[:topN]] = 1
    masktop = (masktop.reshape(nens,nmon)).T
    maskbot = (maskbot.reshape(nens,nmon)).T
    
    pcolmarkers[:,:,r,0] = maskbot.copy()
    pcolmarkers[:,:,r,1] = masktop.copy()



#%% Test with index unraveling (so flattening is not needed anytime)

# Try unravelling the index
kunravel = np.unravel_index(kminmax[r,:,0].astype('int'),shape=(nens,nmon))
v1 = invars[0].copy()
v2 = invars[0].reshape(nens,nmon,nlatr,nlonr)
test1 = v1[kminmax[r,:,0].astype('int'),:,:] # [event (ens, mon) x lat x lon]
test2 = v2[kunravel[0],kunravel[1],:,:]
print(np.nanmax(np.abs(test1-test2).flatten()))

#%% Plot Composites for each region

vnames_in = ("Heat Flux Feedback ($\lambda_a$: $Wm^{-2}\degree C^{-1}$)","$T_2$ (Months)")
cints_all = (np.arange(-45,47.5,2.5),np.arange(0,16.5,0.5))
cmaps_all = ('cmo.balance','inferno')

for r in range(nreg):
    fig = plt.figure(constrained_layout=False, facecolor='w',figsize=(12,14))
    
    gs = fig.add_gridspec(nrows=3, ncols=6, left=.02, right=1,
                          hspace=.2, wspace=0.15)
    
    # Plot Pattern Correlation Grid
    ax0 = fig.add_subplot(gs[0, :])
    ax  = ax0
    ax.set_aspect('equal')
    ax.set_xticks(ens)
    ax.set_yticks(np.arange(1,13,1))
    ax.set_yticklabels(months)
    ax.set_xlabel("Ensemble")
    ax.set_ylabel("Reference Month")
    ax.grid(True,ls='dotted')
    pcm = ax.pcolormesh(ens,np.arange(1,13,1),patcorr_t2_reg[...,r].T,shading='nearest',
                        vmin=-.8,vmax=.8,cmap='cmo.balance')
    cb = fig.colorbar(pcm,ax=ax,fraction=0.025,pad=0.01,orientation='vertical')
    ax.set_title("%s Pattern Correlation (Heat Flux Feedback and $T_2$)"%(regions[r]),
                 fontweight='bold',fontsize=16,)
    
    viz.plot_mask(ens,np.arange(1,13,1),pcolmarkers[:,:,r,1].T,reverse=True,ax=ax,markersize=10,color='yellow',marker="+")
    viz.plot_mask(ens,np.arange(1,13,1),pcolmarkers[:,:,r,0].T,reverse=True,ax=ax,markersize=10,color='yellow',marker="x")

    
    # Plot Each Pattern
    ax1 = fig.add_subplot(gs[1, :3],projection=ccrs.PlateCarree())
    ax2 = fig.add_subplot(gs[1, 3:],projection=ccrs.PlateCarree())
    
    ax3 = fig.add_subplot(gs[2, :3],projection=ccrs.PlateCarree())
    ax4 = fig.add_subplot(gs[2, 3:],projection=ccrs.PlateCarree())
    
    plotvars = [hf_composites,t2_composites]
    
    for j in range(2):
        
        if j == 0: # Start with 0 = Bottom, 1 = Top
            axin  = [ax3,ax4]
            axislab = "Bottom %i \n Composite" % topN
            
        else:
            axin = [ax1,ax2]
            axislab = "Top %i \n Composite" % topN
            
        for a,ax in enumerate(axin):
            
            blabel = [0,0,0,1]
            if a ==0:
                blabel[0] = 1
                ax.text(-0.26, 0.5, axislab, va='bottom', ha='center',
                    rotation='horizontal', rotation_mode='anchor',
                    transform=ax.transAxes,fontsize=14)
            ax  = viz.add_coast_grid(ax,bbox=bboxplot,blabels=blabel,fill_color='gray')
            
            plotvar = plotvars[a][:,:,r,j] # [variable][region and max/min]
            cint_in = cints_all[a]
            
            cf = ax.contourf(lonr,latr,plotvar,
                             cmap=cmaps_all[a],levels=cint_in,extend='both',zorder=-1)
            
            if j == 0:
                ax.set_title(vnames_in[a])
            if j == 1:
                cb = fig.colorbar(cf,ax=ax,orientation='vertical',fraction=0.045)
                
    plt.savefig("%sPattern_Corr_T2_HFF_Composites_lag%i_region%s_topN%i.png" % (figpath,1,regions[r],topN),dpi=200,bbox_inches='tight',transparent=False)


#%% Just plot the top 5 composites


mmname = ["Bot. %i"%topN,"Top %i"%topN]
vnames = []
for a in range(2):
    
    cint_in = cints_all[a]
    invar   = invars[a]
    fig,axs = plt.subplots(2,topN,subplot_kw={'projection':ccrs.PlateCarree()},figsize=(16,4),constrained_layout=True)
    
    for r in tqdm(range(nreg)):
        
        for mm in range(2):
            
            for ik in range(topN):
                
                ax        = axs[mm-1,ik] # Plot top (1) first
                
                # Get index, Name
                k         = int(kminmax[r,ik,mm])
                kens,kmon = np.unravel_index(k,shape=(nens,nmon))
                
                
                blabel = [0,0,0,0]
                if mm == 1:
                    blabel[-1] = 1
                if ik == 0:
                    blabel[0] = 1
                    
                    ax.text(-0.28, 0.35, '%s'% (mmname[mm]), va='bottom', ha='center',
                        rotation='horizontal', rotation_mode='anchor',
                        transform=ax.transAxes,fontsize=14)
                
                plotvar = invar[k,:,:]
                
                ax  = viz.add_coast_grid(ax,bbox=bboxplot,blabels=blabel,fill_color='gray')
                ax.set_title("Ens %i, Mon %i" % (kens+1,kmon+1))
                
                cf = ax.contourf(lonr,latr,plotvar,
                                 cmap=cmaps_all[a],levels=cint_in,extend='both',zorder=-1)
        
        
        plt.savefig("%sPattern_Corr_T2_HFF_Indiv_lag%i_region%s_topN%i_%s.png" % (figpath,1,regions[r],topN,invar_names[a]),dpi=200,bbox_inches='tight',transparent=False)
            
            
            
            
            
        


#%% Plot Relationship in each region for all months and fluxes

fig,ax = plt.subplots(1,1)

scs = []
for b in [0,1,2,4,5,6]:
    
    sc = ax.scatter(hff_reg[:,:,b].flatten(),t2_reg[:,:,b].flatten(),c=bbcol[b],alpha=0.4,label=regions[b])
    scs.append(sc)

ax.legend()
ax.set_ylim([0,25])
ax.set_xlim([-10,50])

ax.set_xlabel("Heat Flux Feedback ($Wm^{-1}\degree C^{-1}$)")
ax.set_ylabel("$T_2$ (Months)")
    
#%% On separate subplots


#fig,axs = viz.init_2rowodd(3,proj=None,figsize=(12,8))

fig,axs = plt.subplots(2,3,figsize=(10,6),constrained_layout=True)

for i,b in enumerate([0,1,2,4,5,6]):
    
    #ax = axs[i]
    ax = axs.flatten()[i]
    
    ax.scatter(hff_reg[:,:,b].flatten(),t2_reg[:,:,b].flatten(),c=bbcol[b],alpha=0.6,label=regions[b],marker="x")
    
    ax.set_ylim([0,12])
    ax.set_xlim([-10,50])
    
    if i > 2:
        ax.set_xlabel("Heat Flux Feedback ($Wm^{-1}\degree C^{-1}$)")
    if i in [0,3]:
        ax.set_ylabel("$T_2$ (Months)")
    
    ax.set_title(regions[b])
    ax.grid(True,ls="dotted")
plt.savefig("%sRegionally_avged_hff_T2_relationship.png" % (figpath),dpi=150,bbox_inches='tight')

    

#%%






