#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Visualize ENSO calculated from calc_ENSO_PIC

@author: gliu
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import cmocean
import cartopy.crs as ccrs
import cartopy

import sys
sys.path.append("/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/03_Scripts/")
from amv import proc,viz
#from matplotlib import gridspec
import cartopy.feature as cfeature


#%% User Edits

datpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/01_Data/ENSO_PIC/"

# "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/SLAB_PIC/"
npznames = ("EOF_ENSO_PIC_FULL.npz","EOF_ENSO_PIC_SLAB.npz")
outpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/02_Figures/Weekly_Meetings/"

# Subset data for enso index calculation
bbox = [120, 290, -25, 25]



for i in range(2):
    npzname = npznames[i]
    
    # Load data
    ld = np.load(datpath+npzname,allow_pickle=True)
    
    eofs   = ld['eofs']
    pcs    = ld['pcs']
    varexp = ld['varexp']
    lon    = ld['lon']
    lat    = ld['lat']
    times  = ld['times']
#%%

# Find dimensions and separate out space
nlon = lon.shape[0]
nlat = lat.shape[0]
npcs = eofs.shape[2] 

# Reshape EOFs variable
eofs = eofs.reshape(nlat,nlon,12,npcs)

# # Conver to 180
# eofs = eofs.transpose(1,0,2)
# lon,eofs=proc.lon360to180(lon,eofs)
# eofs = eofs.transpose(1,0,2)


# Plot!
cmap = cmocean.cm.balance
cint = np.arange(-1,1.1,.1)
n = 2
m = 0





#fig,axs = plt.subplots(1,1,gridspec_kw={'height_ratios':[3,1]},subplot_kw={'projection': ccrs.PlateCarree(central_longitude=180)})
fig,axs = plt.subplots(1,1,subplot_kw={'projection': ccrs.PlateCarree(central_longitude=180)})

ax = axs
ax = viz.init_map(bbox,ax=ax)
pcm = ax.contourf(lon,lat,eofs[:,:,m,n],levels=cint,cmap=cmap,transform = ccrs.PlateCarree())
fig.colorbar(pcm,ax=ax)
ax.set_title("EOF %i, Month %i, Variance Explained %.2f"%(n+1,m+1,varexp[m,n]*100) + "%")



# ENSO Plot
fig,axs = plt.subplots(1,1,figsize=(5,4),subplot_kw={'projection': ccrs.PlateCarree(central_longitude=180)})
ax = axs
ax = viz.init_map(bbox,ax=ax)
pcm = ax.contourf(lon,lat,eofs[:,:,m,n],levels=cint,cmap=cmap,transform = ccrs.PlateCarree())
fig.colorbar(pcm,ax=ax,fraction=0.05,orientation='horizontal')
ax.set_title("EOF %i, Month %i, Variance Explained %.2f"%(n+1,m+1,varexp[m,n]*100) + "%")
#plt.tight_layout()
plt.savefig("%sENSO_Index_PIC_SLAB_EOF%i_Month%i.png"%(outpath,n+1,m+1),dpi=200)






fig,axs = plt.subplots(1,1,figsize=(5,2))
ax = axs
ax.plot(pcs[:,m,n])
ax.set_title("PC %i, Month %i, Variance Explained %.2f"%(n+1,m+1,varexp[m,n]*100) + "%")
plt.savefig("%sPC_Index_PIC_SLAB_EOF%i_Month%i.png"%(outpath,n+1,m+1),dpi=200)











fig,ax = plt.subplots(1,1,subplot_kw={'projection': ccrs.PlateCarree(central_longitude=180)})







crs0 = ccrs.PlateCarree(central_longitude=180)
ax = plt.axes(projection=crs0)
ax.set_extent(bbox,ccrs.PlateCarree())
ax.add_feature(cartopy.feature.LAND)
ax.add_feature(cartopy.feature.OCEAN)
ax.coastlines('50m')
#ax.set_extent(bbox)
pcm = ax.pcolormesh(lon,lat,eofs[:,:,0],transform=ccrs.PlateCarree())



#%% Load the ENSO Component plots for specific variables
vname = "FSNS"
npzname = "ENSOREM_%s_lag1_pcs2_monwin3_ensocomponent.npz"%vname
ld2 = np.load(datpath+npzname,allow_pickle=True)
pat = ld2['ensopattern']
latf = ld2['lat']
lonf = ld2['lon']

imon = 0
ipc = 0

plt.pcolormesh(lonf,latf,pat[imon,:,:,ipc],cmap=cmocean.cm.balance),plt.colorbar()



fig,axs = plt.subplots(1,1,gridspec_kw={'height_ratios':[3,1]},subplot_kw={'projection': ccrs.PlateCarree(central_longitude=180)})
ax = axs
ax = viz.init_map(bbox,ax=ax)
pcm = ax.contourf(lon,lat,eofs[:,:,m,n],levels=cint,cmap=cmap,transform = ccrs.PlateCarree())
fig.colorbar(pcm,ax=ax)
ax.set_title("EOF %i, Month %i, Variance Explained %.2f"%(n+1,m+1,varexp[m,n]*100) + "%")