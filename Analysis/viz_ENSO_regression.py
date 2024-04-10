#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Visualize ENSO regression from output of
calc_enso_general

Copied upper section of compare_CESM1LE_damping on 2024.04.10


Created on Wed Apr 10 13:16:08 2024

@author: gliu
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
import cmocean
import sys
import tqdm

import matplotlib as mpl

#%% Import my modules
sys.path.append("/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/03_Scripts/stochmod/model/")
sys.path.append("/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/03_Scripts/")

from amv import proc,viz
import scm

import importlib
importlib.reload(viz)

#%% User Edits

datpath     = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/CESM-HTR-RCP85/enso/"
searchstr   = "htr_%s_detrend1_ENSOcmp_lag1_pcs3_monwin3_1920to2005_ens%02i.npz" 
vnames      = ["ts","qnet"]
vnames_long = ["SST","$Q_{net}$"]
vunits      = ["$\degree$C","$W m^{-2}$"]

save_da     = True
figpath     = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/03_reemergence/02_Figures/20240411/"

#%% Indicate variables and ensembles
nens  = 42
vv    = 0
ds_all = []
for vv in range(2):
    vname = vnames[vv]
    
    
    # Load everything
    ld_all = []
    for e in range(nens):
        fn = datpath + searchstr % (vname,e+1)   
        ld = np.load(fn,allow_pickle=True)       #  (12, 192, 288, 3)
        ld_all.append(ld['ensopattern'])
    
    lon = ld['lon']
    lat = ld['lat']
    pcs = np.arange(1,4,1)
    mons = np.arange(1,13,1)
    ens  = np.arange(1,43,1)
    
    # Make in DataArray
    pats_all = np.array(ld_all) # (42, 12, 192, 288, 3)
    coords   = dict(ens=ens,mon=mons,lat=lat,lon=lon,pc=pcs)
    da_enso  = xr.DataArray(pats_all,coords=coords,dims=coords,name=vname)
    da_enso = proc.lon360to180_xr(da_enso)
    if save_da:
        savename = fn[:-6] + "ALL.nc"
        print("Saving merged patterns in %s" % savename)
        edict    = {vname:{'zlib':True}}
        da_enso.to_netcdf(savename,encoding=edict)
    
    ds_all.append(da_enso)
    
    
#%% Plotting Params

mons3 = proc.get_monstr()
proj  = ccrs.PlateCarree()
mpl.rcParams['font.family'] = 'JetBrains Mono'
bboxplot                    = [-80,0,20,65]


fsz_title=20
fsz_axis = 14

        
#%% Subset to region for plotting


bbsel = [-80,0,0,65]
ds_reg = [proc.sel_region_xr(ds,bbsel) for ds in ds_all]

lonr,latr = ds_reg[0].lon.values,ds_reg[0].lat.values


#%% Visualize the patterns (Ensemble Average)
im=0
nn=0


for vv in range(2):
    if vv == 0:
        vmax = .25
    elif vv == 1:
        vmax = 30
    pmesh    = True
    
    
    for nn in range(3):
        
        for im in tqdm.tqdm(range(12)):
            seltitle = "%s ENSO Pattern \n%s, PC %02i, Ens Avg." % (vnames_long[vv],mons3[im],nn+1)
            
            
            vname  = vnames[vv]
            dsplot = ds_reg[vv].isel(mon=im,pc=nn).mean('ens') # [Lat x Lon]
            
            fig,ax,_ = viz.init_orthomap(1,1,bboxplot,figsize=(10,6.5))
            
            ax = viz.add_coast_grid(ax,bboxplot,fill_color="k")
            if pmesh:
                pcm = ax.pcolormesh(lonr,latr,dsplot,transform=proj,cmap='cmo.balance',vmin=-vmax,vmax=vmax)
            else:
                pcm = ax.contourf(lonr,latr,dsplot,transform=proj,cmap='cmo.balance')
            ax.set_title(seltitle,fontsize=fsz_title)
            cb = fig.colorbar(pcm,ax=ax,orientation='horizontal',fraction=0.05,pad=0.01)
            cb.set_label("%s per 1 Std Dev. PC %02i" % (vunits[vv],nn+1),fontsize=fsz_axis)
            
            savename = "%sENSORegr_CESM1htr_%s_EnsMean_mode%i_mon%02i.png" % (figpath,vnames[vv],nn+1,im+1)
            
            plt.savefig(savename,dpi=150,bbox_inches='tight')
            plt.close()

#%% Visualize Annual mean as well for a single mode

nn = 0
lonsel=-58
latsel=45

vmaxes = [.25,15]
fig,axs,_ = viz.init_orthomap(1,2,bboxplot,figsize=(10,6.5))


vmeans = []
for vv in range(2):
    
    ax = axs[vv]
    ax = viz.add_coast_grid(ax,bboxplot,fill_color="k")
    ax.set_title(vnames_long[vv],fontsize=fsz_axis)
    vmax = vmaxes[vv]
    
    dsplot = ds_reg[vv].isel(pc=nn).mean('ens').mean('mon') # [Lat x Lon]
    if pmesh:
        pcm = ax.pcolormesh(lonr,latr,dsplot,transform=proj,cmap='cmo.balance',vmin=-vmax,vmax=vmax)
    else:
        pcm = ax.contourf(lonr,latr,dsplot,transform=proj,cmap='cmo.balance')
    
    ax.plot(lonsel,latsel,marker="+",markersize=25,transform=proj,c='yellow')
    
    cb = fig.colorbar(pcm,ax=ax,orientation='horizontal',fraction=0.05,pad=0.01)
    cb.set_label("%s per 1 Std Dev. PC %02i" % (vunits[vv],nn+1),fontsize=fsz_axis)
    savename = "%sENSORegr_CESM1htr_EnsMean_mode%i_AnnMean.png" % (figpath,nn+1)
    
    plt.savefig(savename,dpi=150,bbox_inches='tight')
    vmeans.append(dsplot)

#%% Check how the select point compares
fig,axs = plt.subplots(1,2,figsize=(10,4.5))

xls = [.2,15]

for vv in range(2):
    ax = axs[vv]
    ax.set_title(vnames_long[vv])
    plotvar = vmeans[vv].values.flatten()
    bw      = proc.calc_binwidth(plotvar)
    bins    = np.arange(np.nanmin(plotvar),np.nanmax(plotvar)+bw,bw)
    ax.hist(plotvar,bins=bins,edgecolor='w',lw=0.75,color='dimgray',alpha=0.85,label="Histogram")
    
    ax = viz.add_ticks(ax)
    
    ax.set_ylabel("Frequency")
    
    ax.set_xlim([-xls[vv],xls[vv]])
    ax.set_xlabel("%s" % (vunits[vv]))
    
    ax.axvline([vmeans[vv].sel(lon=lonsel,lat=latsel,method='nearest').values],color="goldenrod",
               lw=3,label="Lon: %.2f, Lat: %.2f" % (lonsel,latsel))
    ax.axvline([vmeans[vv].sel(lon=-30,lat=50,method='nearest').values],color="magenta",lw=3
               ,label="Lon %.2f, Lat%.2f" % (-30,50))
    ax.legend(ncol=1,fontsize=8)

plt.suptitle("Regression Coefficients for ENSO PC %i \nEns Mean, Ann. Avg." % (nn+1))
savename = "%sENSORegr_CESM1htr_EnsMean_mode%i_AnnMean_Histogram.png" % (figpath,nn+1)
plt.savefig(savename,dpi=150,bbox_inches='tight')