#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Check Land Ice Masking

dummy script. see stochmod/analysis/check_landicemask.ipynb!

Created on Thu May 18 21:01:38 2023

@author: gliu
"""

import numpy as np
import cartopy as crs
import glob
import xarray as xr
import time

machine   = "stormtrack"
atmpath   = "/vortex/jetstream/climate/data1/yokwon/CESM1_LE/downloaded/atm/proc/tseries/monthly/"


#b.e11.B1850C5CN.f09_g16.005.cam.h0.ICEFRAC.040001-049912.nc 
#e.e11.E1850C5CN.f09_g16.001.cam.h0.ICEFRAC.010101-019912.nc

nclist_full = glob.glob(atmpath+"ICEFRAC/b.e11.B1850*.nc")
nclist_slab = glob.glob(atmpath+"ICEFRAC/e.e11.E1850*.nc")
nclist_slab.sort()
nclist_full.sort()

nclist_slab = [nc for nc in nclist_slab if "RAMP" not in nc]
nclist_full = [nc for nc in nclist_full if "RAMP" not in nc]

st = time.time()
dsfull = xr.open_mfdataset(nclist_slab,combine="nested",concat_dim="time").ICEFRAC.load()
dsslab = xr.open_mfdataset(nclist_slab,combine="nested",concat_dim="time").ICEFRAC.load()
print("Loaded data in %.2fs "% (time.time()-st))

icethres = 0.05
iceslab = dsslab.values
icefull = dsfull.values

# MAKE THE MASK
fullmask = (icefull > 0.05).prod(0)
slabmask  = (iceslab > 0.05).prod(0)
#%% Save the ice fraction
savepath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/02_stochmod/ICEFRAC/"
encoding = {"ICEFRAC": {'zlib': True}}
names = ["slab","full"]
for m,msk in enumerate([dsslab,dsfull]):
    savename = "%sICEFRAC_PIC_%s.nc" % (savepath,names[m].upper())
    msk.to_netcdf(savename,encoding=encoding)
#%%
outpath = "/stormtrack/home/glliu/full_mask_test.npy"
np.save(outpath,fullmask)
outpath = "/stormtrack/home/glliu/slab_mask_test.npy"
np.save(outpath,slabmask)


np.save("")

lon = dsslab.lon
lat = dsslab.lat
bbox = [-80,0,20,65]

def plotmask(mask,lon,lat,bbox):

    fig,ax = plt.subplots(1,1,constrained_layout=True,figsize=(8,3),
                           subplot_kw={'projection':ccrs.PlateCarree()})
    
    pcm = ax.pcolormesh(lon,lat,mask,vmin=0,vmax=1,cmap="inferno")
    ax.set_extent(bbox)
    ax.coastlines()
    #viz.add_coast_grid(ax,bbox=bbox_spg,fill_color="gray",blabels=blabel)
    fig.colorbar(pcm,ax=ax,orientation='horizontal',fraction=0.025,pad=0.01)
    plt.show()
    return None

plotmask(fullmask,lon,lat,bbox)#,plt.show()

savename = "%sLand_Ice_Mask_Comparison.png" % figpath
plt.savefig(savename,dpi=150,bbox_inches="tight")


#%% on local device




# Select a region

# if machine == "stormtrack":
#     datpath   = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/02_stochmod/Model_Data/model_input/"
#     slabname  = "SLAB_landicemask360.npy"
#     fullname  = "FULL_landicemask360.npy"
#     htrname   = "limask180_FULL-HTR.npy"
#     htrname_1 = "pacific_limask_180global.npy" 
# else:
#     datpath = ""
#     slabname = "CESM-SLAB_landicemask360.npy"
#     fullname = "CESM-FULL_landicemask360.npy"
# mconfigs = ("slab","full","htr","htr-pac")
# npnames  = [slabname,fullname,htrname,htrname_1]


# masks = []
# for ii in range(4):
#     masks.append(np.load(datpath+npnames[ii],allow_pickle=True))
    
    
