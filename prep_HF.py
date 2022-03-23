#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Prepare Heat Flux Damping for input into the stochastic model
Takes the output from the calc_HF.py script.

Created on Fri Nov 20 12:09:01 2020

@author: gliu
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import cmocean
import cartopy.crs as ccrs
import cartopy

import sys
sys.path.append("/home/glliu/00_Scripts/01_Projects/00_Commons/")
#sys.path.append("/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/03_Scripts/")
from amv import proc,viz
from matplotlib import gridspec

from scipy import stats
from scipy.io import loadmat,savemat
import cartopy.feature as cfeature

#%%


datpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_PIC_SLAB/03_HFCALC/FULL/"
outpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_PIC_SLAB/03_HFCALC/proc/"
#datpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/SLAB_PIC/"
#outpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/02_Figures/Weekly_Meetings/"
lipath = "/home/glliu/01_Data/00_Scrap/landicemask_enssum.npy"
llpath = "/home/glliu/01_Data/CESM1_LATLON.mat"
#lipath = "%slandicemask_enssum.npy" % datpath # Land ice mask (.npy file)
#llpath = "%s../CESM1_LATLON.mat" % datpath # Lat/Lon File

# Experiment Parameters ----
vname = 'NHFLX'
lags = [1,2,3] # Lags to include
lagstr = "1-2" # Labeling for plots (DONT FORGET TO CHANGE!)

# Significance Testing Results
mode  = 4 # (1) No mask (2) SST only (3) Flx only (4) Both
p     = 0.05
tails = 2
dof   = 1898 - 1 - 2 - 2


# ENSO Removal Options
ensorem = 1 # Set to 1 if ENSO was removed
ensolag = 1 # Lag between enso removal and variable
pcrem   = 2 # PCs of enso removed
emonwin = 3 # Month window of ENSO Removal

# Flux Calculation options
monwin  = 3 # Window of months
lagname = '123' # Lags to include
flux    = 'NHFLX'
mconfig = "FULL"

# Save options
savevar = 1 # option to save variables

# Plotting Options
plotfigs = False
bbox = [280-360, 360-360, 0, 70]

ensoexp = "lag%i_pcs%i_monwin%i" % (ensolag,pcrem,emonwin)
expin   = "monwin%i_lags%s_ensorem%i_%s" % (monwin,lagname,ensorem,ensoexp)
testres = ""
#NHFLX_Damping_monwin3_lags123_ensorem1_lag1_pcs2_monwin3.npy

#%% Functions

def calc_corrthres(p,tails,dof):

    ptilde    = 1-p/tails
    critval   = stats.t.ppf(ptilde,dof)
    corrthres = np.sqrt(1/ ((dof/np.power(critval,2))+1))
    return corrthres


#%%
# Load in the damping and correlation files
damping = np.load("%s%s_Damping_%s.npz.npy"% (datpath,vname,expin))  # [mon lag lat lon] 
rsst    = np.load("%sSST_Autocorrelation_%s.npz.npy"% (datpath,expin))  # [mon lag lat lon] 
rflx    = np.load("%s%s_Crosscorrelation_%s.npz.npy"% (datpath,vname,expin))  # [mon lag lat lon] 

# Load land ice mask and lat/lon
limask = np.load(lipath) # [lat x lon]

# Load lat/lon
lm = loadmat(llpath)
lon = lm['LON'].squeeze()
lat = lm['LAT'].squeeze()

#%% Significance Testing

#Calcuate the correlation threshold
corrthres = calc_corrthres(p,tails,dof)

# Create and mask
msst = np.zeros(damping.shape)
mflx = np.zeros(damping.shape)
msst[rsst > corrthres] = 1
mflx[rflx > corrthres] = 1

if mode == 1:
    mtot = np.ones(damping.shape)     # Total Frequency of successes
    mall = np.copy(mtot)              # Mask that will be applied
    mult = 1
elif mode == 2:
    mtot = np.copy(msst)
    mall = np.copy(msst)
    mult = 1
elif mode == 3:
    mtot = np.copy(mflx)
    mall = np.copy(mflx)  
    mult = 1
elif mode == 4:
    mtot = msst + mflx
    mall = msst * mflx
    mult = 2

# Apply mask
dampingt = damping * mall

dampingmasked = dampingt.copy()
#%% Further Processing

# Average over selected lags
ilags = [l-1 for l in lags] # get indices
dampingt = dampingt[:,ilags,:,:].mean(1).squeeze()



lagstr = ""
for l in lags:
    lagstr += str(l)
    


# Apply Land/Ice Mask
dampingm = dampingt * limask[None,:,:]

# Flip longiude coordinates ([mon lat lon] --> [lon x lat x mon])
lon1,dampingw = proc.lon360to180(lon,dampingm.transpose(2,1,0))

# Multiple by 1 to make positive upwards
dampingw *= -1

# Save Result
if savevar == 1:
    # Save variables, matching matlab script output
    savedict = {'damping':dampingw, 'LON1': lon1, 'LAT':lat }
    outname = "%s%s_PIC_%sdamping_monwin%i_sig%03d_dof%03d_mode%i.mat" % (datpath,mconfig,flux.lower(),monwin,p*100,dof,mode)
    savemat(outname,savedict)
    print("Saved data to %s" % outname)


if plotfigs:
    #%% Vizualize the results ----
    
    #% Just plot the damping values
    fig,ax = plt.subplots(1,1,subplot_kw={'projection':ccrs.PlateCarree()},figsize=(5,4))
    cint = np.arange(-50,55,5)
    ax = viz.init_map(bbox,ax=ax)
    pcm = ax.contourf(lon1,lat,np.mean(dampingw,2).T,cint,cmap=cmocean.cm.balance)
    cl = ax.contour(lon1,lat,np.mean(dampingw,2).T,cint,colors="k",linewidths = 0.5)
    ax.clabel(cl,fmt="%i",fontsize=8)
    ax.add_feature(cfeature.LAND,color='k')
    ax.set_title(r"Ann. Mean $\lambda_{a,%s}$ (Lags %s)" % (vname,lagstr)+ "\n"+r"p = %.2f | $\rho$ > %.2f " % (p,corrthres),fontsize=12)
    plt.colorbar(pcm,ax=ax,orientation='horizontal',fraction=0.040, pad=0.05)
    plt.savefig(outpath+"SLAB_PIC_%s_Damping__mode%i_monwin%i_lags%s_sig%03d.png"%(vname,mode,monwin,lagstr,p*100),dpi=200)
    
    #%% Just plot specific month .....
    m = 1
    fig,ax = plt.subplots(1,1,subplot_kw={'projection':ccrs.PlateCarree()},figsize=(5,5))
    cint = np.arange(-50,55,5)
    ax = viz.init_map(bbox,ax=ax)
    pcm = ax.contourf(lon1,lat,dampingw[:,:,m-1].T,cint,cmap=cmocean.cm.balance)
    cl = ax.contour(lon1,lat,dampingw[:,:,m-1].T,cint,colors='k',linewidths = 0.5)
    ax.clabel(cl,fmt="%i",fontsize=12)
    ax.add_feature(cfeature.LAND,color='gray')
    ax.set_title(r"Month: %i $\lambda_{a,%s}$ (Lags %s)" % (m,vname,lagstr)+ "\n"+r"p = %.2f | $\rho$ > %.2f " % (p,corrthres),fontsize=12)
    plt.colorbar(pcm,ax=ax,orientation='horizontal',fraction=0.040, pad=0.05)
    plt.savefig(outpath+"SLAB_PIC_%s_Damping_month%i_mode%i_monwin%i_lags%s_sig%03d.png"%(vname,m,mode,monwin,lagstr,p*100),dpi=200)
    #%% Make Plot of success and damping values
    
    fig,axs = plt.subplots(1,2,subplot_kw={'projection':ccrs.PlateCarree()},figsize=(6,4))
    
    bbox = [-60,10,50,75]
    
    ax = axs[0]
    cint = np.arange(0,1.1,.1)
    ax = viz.init_map(bbox,ax=ax)
    pcm = ax.contourf(lon,lat,mfreq.T/maxscore,cint,cmap=cmap)
    cl = ax.contour(lon,lat,mfreq.T/maxscore,cint,colors="k",linewidths = 0.5)
    ax.clabel(cl,np.arange(0,1.2,0.2),fmt="%.1f")
    ax.set_title("% of Sig. Values\n"+ r"Mode = %i, Max = %i " % (mode,maxscore))
    plt.colorbar(pcm,ax=ax,orientation="horizontal")
    #plt.savefig(outpath+"%s_SigPts_monwin%i_lags12_sig%03d.png"%(flux,monwin,p*100),dpi=200)
    
    
    
    ax = axs[1]
    cint = np.arange(-50,55,5)
    ax = viz.init_map(bbox,ax=ax)
    pcm = ax.contourf(lon,lat,np.nanmean(dampseason,2).T,cint,cmap=cmocean.cm.balance)
    cl = ax.contour(lon,lat,np.nanmean(dampseason,2).T,cint,colors="k",linewidths = 0.5)
    ax.clabel(cl,fmt="%i")
    ax.set_title("%s Damping (Ann, Lag, Ens Avg)\n" % flux+ r"| p = %.2f | $\rho$ > %.2f " % (p,corrthres))
    plt.colorbar(pcm,ax=ax,orientation="horizontal")
    plt.savefig(outpath+"%s_Damping_and_SigPts_mode%i_monwin%i_lags12_sig%03d.png"%(flux,mode,monwin,p*100),dpi=200)

#% Just plot the damping values
fig,ax = plt.subplots(1,1,subplot_kw={'projection':ccrs.PlateCarree()},figsize=(5,4))
cint = np.arange(-50,55,5)
ax = viz.init_map(bbox,ax=ax)
pcm = ax.contourf(lon1,lat,np.mean(dampingw,2).T,cint,cmap=cmocean.cm.balance)
cl = ax.contour(lon1,lat,np.mean(dampingw,2).T,cint,colors="k",linewidths = 0.5)
ax.clabel(cl,fmt="%i",fontsize=10)
ax.add_feature(cfeature.LAND,color='gray')
ax.set_title(r"CESM-SLAB Annual Mean $\lambda_{a,%s}$ (Avg. Lags %s)" % (vname,lagstr)+ "\n"+r"DOF= %i | p = %.2f | R > %.2f " % (dof,p,corrthres),fontsize=12)
plt.colorbar(pcm,ax=ax,orientation='horizontal',fraction=0.05, pad=0.05)
plt.savefig(outpath+"SLAB_PIC_%s_Damping_mode%i_monwin%i_lags%s_sig%03d.png"%(vname,mode,monwin,lagstr,p*100),dpi=200)

#%% Just plot specific month .....
m = 1
fig,ax = plt.subplots(1,1,subplot_kw={'projection':ccrs.PlateCarree()},figsize=(5,4))
cint = np.arange(-70,75,5)
ax = viz.init_map(bbox,ax=ax)
pcm = ax.contourf(lon1,lat,dampingw[:,:,m-1].T,cint,cmap=cmocean.cm.balance)
cl = ax.contour(lon1,lat,dampingw[:,:,m-1].T,cint,colors='k',linewidths = 0.5)
ax.clabel(cl,fmt="%i",fontsize=10)
ax.add_feature(cfeature.LAND,color='gray')
ax.set_title(r"CESM-SLAB %s $\lambda_{a,%s}$ (Lags %s)" % (viz.return_mon_label(m),vname,lagstr)+ "\n"+r"DOF= %i | p = %.2f | R > %.2f " % (dof,p,corrthres),fontsize=12)
plt.colorbar(pcm,ax=ax,orientation='horizontal',fraction=0.040, pad=0.05)
plt.savefig(outpath+"SLAB_PIC_%s_Damping_month%i_mode%i_monwin%i_lags%s_sig%03d.png"%(vname,m,mode,monwin,lagstr,p*100),dpi=200)
#%% Make Plot of success and damping values

fig,axs = plt.subplots(1,2,subplot_kw={'projection':ccrs.PlateCarree()},figsize=(6,4))

bbox = [-60,10,50,75]

ax = axs[0]
cint = np.arange(0,1.1,.1)
ax = viz.init_map(bbox,ax=ax)
pcm = ax.contourf(lon,lat,mfreq.T/maxscore,cint,cmap=cmap)
cl = ax.contour(lon,lat,mfreq.T/maxscore,cint,colors="k",linewidths = 0.5)
ax.clabel(cl,np.arange(0,1.2,0.2),fmt="%.1f")
ax.set_title("% of Sig. Values\n"+ r"Mode = %i, Max = %i " % (mode,maxscore))
plt.colorbar(pcm,ax=ax,orientation="horizontal")
#plt.savefig(outpath+"%s_SigPts_monwin%i_lags12_sig%03d.png"%(flux,monwin,p*100),dpi=200)



ax = axs[1]
cint = np.arange(-50,55,5)
ax = viz.init_map(bbox,ax=ax)
pcm = ax.contourf(lon,lat,np.nanmean(dampseason,2).T,cint,cmap=cmocean.cm.balance)
cl = ax.contour(lon,lat,np.nanmean(dampseason,2).T,cint,colors="k",linewidths = 0.5)
ax.clabel(cl,fmt="%i")
ax.set_title("%s Damping (Ann, Lag, Ens Avg)\n" % flux+ r"| p = %.2f | $\rho$ > %.2f " % (p,corrthres))
plt.colorbar(pcm,ax=ax,orientation="horizontal")
plt.savefig(outpath+"%s_Damping_and_SigPts_mode%i_monwin%i_lags12_sig%03d.png"%(flux,mode,monwin,p*100),dpi=200)

#%% Plot differences in damping pattern

lonf = 330
latf = 50

klon,klat = proc.find_latlon(lonf,latf,lon,lat)

damppt = dampingmasked[:,:,klat,klon]


fig,ax = plt.subplots(1,1)
for i in range(3):
    ax.plot(damppt[:,i],label="Lag %i"% (i+1))
ax.legend()




#dampingw


