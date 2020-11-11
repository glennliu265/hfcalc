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
from matplotlib import gridspec



#%% User Edits

datpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/SLAB_PIC/"
npzname = "EOF_ENSO_PIC_SLAB.npz"

# Subset data for enso index calculation
bbox = [120, 290, -40, 40]




# Load data
ld = np.load(datpath+npzname,allow_pickle=True)
eofs = ld['eofs']
pcs = ld['pcs']
varexp=ld['varexp']
lon = ld['lon']
lat = ld['lat']
times = ld['times']


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
n = 0
m = 0

fig,axs = plt.subplots(2,1,gridspec_kw={'height_ratios':[3,1]},subplot_kw={'projection': ccrs.PlateCarree(central_longitude=180)})
ax = axs[0]
ax = viz.init_map(bbox,ax=ax)
pcm = ax.contourf(lon,lat,eofs[:,:,m,n],levels=cint,cmap=cmap,transform = ccrs.PlateCarree())
fig.colorbar(pcm,ax=ax)
ax.set_title("EOF %i, Month %i, Variance Explained %.2f"%(n+1,m+1,varexp[m,n]*100) + "%")

ax = axs[1]
ax.plot(pcs[:,n])
ax.set_title("PC %i"%(n+1))




fig,ax = plt.subplots(1,1,subplot_kw={'projection': ccrs.PlateCarree(central_longitude=180)})







crs0 = ccrs.PlateCarree(central_longitude=180)
ax = plt.axes(projection=crs0)
ax.set_extent(bbox,ccrs.PlateCarree())
ax.add_feature(cartopy.feature.LAND)
ax.add_feature(cartopy.feature.OCEAN)
ax.coastlines('50m')
#ax.set_extent(bbox)
pcm = ax.pcolormesh(lon,lat,eofs[:,:,0],transform=ccrs.PlateCarree())

