#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Compare Heat Flux Feedbacks between different estimates

Created on Sun Jun 22 11:11:08 2025

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

#%% Set up project directory paths

projpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/"
figpath  = "/Users/gliu/Downloads/02_Research/01_Projects/05_SMIO/02_Figures/20250625/"

proc.makedir(figpath)

# Set Paths
anompath = projpath + "anom/"
ensopath = projpath + "enso/"
hffpath  = projpath + "hff/"
maskpath = projpath + "masks/"
procpath = projpath + "proc/"

# Check for directories
dirs     = [anompath,ensopath,hffpath,maskpath,procpath]
[proc.makedir(path) for path in dirs]

#%% Indicate file names, load

ncs    = ['OAFLUX_qnet_damping_Global_1984to2007_ensorem1_detrend1.nc',
       'ERA5_qnet_hfdamping_NAtl_1979to2024_ensorem1_detrend1.nc',]
ds_all = [xr.open_dataset(hffpath + nc).load() for nc in ncs]

#%% Check values over bounxing box

bbox_spgne = [-40,-15,52,62]
dsreg      = [proc.sel_region_xr(ds,bbox_spgne) for ds in ds_all]

qnet_dampings = [ds.qnet_damping for ds in dsreg]

dof_era5   = (2021-1979 + 1 - 2 - 1) * 3
dof_oaflux = (2007-1984 + 1 - 2 - 1) * 3

dofs       = [dof_oaflux,dof_era5]

expnames_long = ["OAFLUX (1984-2007)","ERA5 (1979-2024)"]
expcols = ['purple','gray',]

#%% Calculate the Significance of each case (copied from smio/analysis/check_sig_heatflux.py)
# Set Significance calculation settings
"""
Some Options for Signifiance Testing

pilot      : Same as was used for the SSS Paper
noPositive : Just set positive HFF to zero
p10        : Use p = 0.10
p20        : Use p = 0.20
AConly     : Test Autocorrelation Only

"""
signame ="p10" #"p20"# "noPositive"#"noPositive"#"pilot" #"noPositive" # 
print("Significance Testing Option is: %s" % (signame))

setdict = {  # Taken from hfcalc_params
    'ensorem': 1,      # 1=enso removed, 0=not removed
    'ensolag': 1,      # Lag Applied toENSO and Variable before removal
    'monwin': 3,      # Size of month window for HFF calculations
    'detrend': 1,      # Whether or not variable was detrended
    'tails': 2,      # tails for t-test

    'p': 0.05,   # p-value for significance testing
    'sellags': [0,],   # Lags included (indices, so 0=lag1)
    'lagstr': "lag1",  # Name of lag based on sellags
    # Significance test option: 1 (No Mask); 2 (SST autocorr); 3 (SST-FLX crosscorr); 4 (Both), 5 (Replace with SLAB values),6, zero out negative values
    'method': 4
}

if signame == "pilot":
    setdict['method'] = 4 # Apply Significance Testing to Both
elif signame == "noPositive":
    setdict['method'] = 6 # Apply Significance Testing to Both
elif signame == "p10":
    setdict['p'] = 0.10
elif signame == "p20":
    setdict['p'] = 0.20
elif signame == "AConly":
    setdict['method'] = 2


sigmasks = []
for ii in range(2):
    hff  = dsreg[ii].qnet_damping
    rsst = dsreg[ii].sst_autocorr
    rflx = dsreg[ii].sst_flx_crosscorr
    
    dof     = dofs[ii]
    
    # Compute and Apply Mask
    #st = time.time()
    dampingmasked, freq_success, sigmask = scm.prep_HF(hff, rsst, rflx,
                                                       setdict['p'], setdict['tails'], dof, setdict['method'],
                                                       returnall=True)  # expects, [month x lag x lat x lon], should generalized with ensemble dimension?
    #print("Completed significance testing in %.2fs" % (time.time()-st))
    
    dampingout = dampingmasked[:,setdict['sellags'],:,:].squeeze()
    dampingout = xr.where(np.isnan(dampingout),0.,dampingout)
    
    
    sigmasks.append(sigmask)
    
    

#%% Quick Plot of the patterns

fsz_title   = 22
fsz_axis    = 18
fsz_tick    = 14
proj        = ccrs.PlateCarree()

mons3       = proc.get_monstr()
imon        = 0
ilag        = 0
cints       = np.arange(-80,85,5)

for imon in range(12):
    fig,axs = plt.subplots(1,2,figsize=(12,4.5),constrained_layout=True)
    

    for ii in range(2):
        ax      = axs[ii]
        plotvar = qnet_dampings[ii].isel(month=imon,lag=ilag) * -1
        pcm     = ax.contourf(plotvar.lon,plotvar.lat,plotvar,levels=cints,cmap='cmo.balance',extend='both')
        ax.set_title(expnames_long[ii],fontsize=fsz_axis)
        
        cl = ax.contour(plotvar.lon,plotvar.lat,plotvar,levels=cints,colors="w",linewidths=0.55)
        clbl = ax.clabel(cl,fontsize=fsz_tick)
        viz.add_fontborder(clbl,c="dimgray",w=2)
        
        sigmask = sigmasks[ii]
        is_insig = np.where(np.isnan(sigmask[imon,ilag,:,:]),1,0)#xr.where()
        viz.plot_mask(plotvar.lon,plotvar.lat,is_insig.T,reverse=True,marker="x",
                      geoaxes=False,ax=ax,color='k',markersize=0.2)
        
        
    cb= viz.hcbar(pcm,ax=axs.flatten(),fontsize=fsz_tick,fraction=0.065)
    cb.set_label(r"Heat Flux Damping [$W m^{-2} degree C^{-1} $]",fontsize=fsz_tick)
    plt.suptitle("Month %s, Lag %i" % (mons3[imon],ilag+1),fontsize=fsz_title)
    
    
    figname = "%sHFF_Comparison_OAFLUX_ERA5_lag%02i_mon%02i.png" % (figpath,ilag,imon,)
    plt.savefig(figname,dpi=150,bbox_inches='tight')

#%% Examine the monthly progression averaged over the months

def pointwise_stdev(ds):
    var_arr = ds.transpose('lat','lon','month').data
    nlat,nlon,ntime = var_arr.shape
    var_arr = var_arr.reshape(nlat*nlon,ntime)
    return np.nanstd(var_arr,0)

for ilag in range(3):
    fig,ax = viz.init_monplot(1,1,figsize=(8,4.5))
    
    ax.set_title("Region Mean (Line) and 1$\sigma$ (Shading) Damping Estimate\nLag %02i" % (ilag+1))
    for ii in range(2):
        plotvar = dsreg[ii].isel(lag=ilag) * -1 
        
        mu      = proc.area_avg_cosweight(plotvar)
        sigma   = pointwise_stdev(plotvar)
        
        
        
        ax.plot(mons3,mu,label=expnames_long[ii],c=expcols[ii],lw=2.5,marker="o")
        ax.fill_between(mons3,mu-sigma,mu+sigma,alpha=0.1,color=expcols[ii],zorder=-1)
    ax.legend()
    ax.set_ylim([-60,60])
    
    figname = "%sHFF_Comparison_OAFLUX_ERA5_Ravg_lag%02i.png" % (figpath,ilag)
    plt.savefig(figname,dpi=150,bbox_inches='tight')

#%% Determine the Significance?
    
