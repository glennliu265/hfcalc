#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Attempts to write functions to calculate heat flux
Testing script

- Includes seasonal damping plots


Created on Mon Jul 19 16:14:32 2021

@author: gliu
"""

import xarray as xr
import numpy as np
import glob
import time
import cmocean
from tqdm import tqdm

import sys

import cartopy.crs as ccrs
import matplotlib.pyplot as plt

#%% Set working location
stormtrack = 0

if stormtrack == 1:
    # Module Paths
    sys.path.append("/home/glliu/00_Scripts/01_Projects/00_Commons/")
    sys.path.append("/home/glliu/00_Scripts/01_Projects/01_AMV/02_stochmod/stochmod/model/")
elif stormtrack == 0:
    # Module Paths
    sys.path.append("/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/03_Scripts/stochmod/model/")
    sys.path.append("/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/03_Scripts/")
    
    datpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/"
    outpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/02_Figures/20211021/"

    lipath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/01_Data/model_input/"
    llpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/01_Data/model_input/"

from amv import proc,viz
import scm

#mconfig = "SLAB_FULL"
#datpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_PIC_SLAB/02_ENSOREM/%s/" % mconfig
#outpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_PIC_SLAB/03_HFCALC/%s/"  % mconfig


tails      = 2
p          = 0.05
mode       = 4 # 1 = No Mask; 2 = SST autocorr; 3 --> SST-FLX cross corr; 4 = Both
sellags    = [0]
maskval    = 0 # Set the masking value
saveoutput = True

lagstr = "lag"
for i in sellags:
    lagstr += str(sellags[i]+1)
mons3       = ('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
bboxplot = [-100,20,0,80]
#%% Quick fill-in functions

def load_dampraw(mconfig,datpath):
    inpaths = datpath+"CESM-"+mconfig+"-Damping/"
    damping = np.load(inpaths+"NHFLX_Damping_monwin3_lags123_ensorem1_lag1_pcs2_monwin3.npz.npy")
    rflx    = np.load(inpaths+"NHFLX_Crosscorrelation_monwin3_lags123_ensorem1_lag1_pcs2_monwin3.npz.npy")
    rsst    = np.load(inpaths+"SST_Autocorrelation_monwin3_lags123_ensorem1_lag1_pcs2_monwin3.npz.npy")
    return damping,rsst,rflx


#%% Test the prep_HF script
from scipy import stats

# Load some data
# Note this only works locally, need to make equivalent on stormtrack
if stormtrack == 0:
    lon,lat=scm.load_latlon(datpath=llpath,lon360=True)
    lon180,_ = scm.load_latlon(datpath=llpath)
    limask = np.load(lipath+"../landicemask_enssum.npy")

mconfigs =["PIC-SLAB","PIC-FULL"]
dofs     =[898 - 1 - 2 - 2, 1898 - 1 - 2 - 2] 

dampings   = []
rssts      = []
rflxs      = []
dampingfin = []
dunmask    = [] # Reshaped, unmasked damping values
rthres     = [] # Correlation thresholds
mtots      = []
malls      = []
for i,mcf in enumerate(mconfigs):

    # Load Data
    a,b,c = load_dampraw(mcf,datpath)
    dampings.append(a)
    rssts.append(b)
    rflxs.append(c) #Flux must be positive into atm
    
    # Apply Masks
    d,mtot,mall  = scm.prep_HF(a,b,c,p,tails,dofs[i],mode,maskval=0,returnall=True)
    
    
    d2 = scm.prep_HF(a,b,c,p,tails,dofs[i],1,maskval=0)
    
    # Calculate and save correlation threshold
    dof = dofs[i]
    ptilde    = 1-p/tails
    critval   = stats.t.ppf(ptilde,dof)
    corrthres = np.sqrt(1/ ((dof/np.power(critval,2))+1))
    rthres.append(corrthres)
    
    #dampingfin.append(d)
    
    # Postprocess
    dampingw = scm.postprocess_HF(d,limask,sellags,lon)
    dampingfin.append(dampingw)
    
    # Append unmasked data
    d2w = scm.postprocess_HF(d2,limask,sellags,lon)
    dunmask.append(d2w)
    
    # Append mask info
    mall = scm.postprocess_HF(mall,limask,sellags,lon,pos_upward=False)
    mtots.append(mtot)
    malls.append(mall)
    

#%% For Mode 5, replace insignificant values in FULL with those from SLAB
debug = False

if mode == 5:
    idfull  = malls[1] == 0 # Insignificant points in FULL
    print("Replacing %i points"%(np.sum(idfull.flatten())))
    dampingfin[1][idfull] = dampingfin[0][idfull] # Substitute Corresponding points
    malls[1][idfull] = malls[0][idfull] # Substitute for mask as well
    
    
    if debug:
        plt.pcolormesh(idfull[...,0].T),plt.colorbar()
        test = dampingfin[1].copy()
        test[idfull] = dampingfin[0][idfull] # Replace with values
        diff = dampingfin[1] - test # Take Difference
        plt.pcolormesh(diff[:,:,0].T),plt.colorbar()
    
    
    
    



#%% Save the Output
# -----------------
if saveoutput:
    # Save Full PIC
    picout = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/01_Data/model_input/"
    outname = "%sFULL_PIC_NHFLX_Damping_monwin3_sig005_dof%i_mode%i_%s.npy" % (picout,dofs[1],mode,lagstr)
    np.save(outname,dampingfin[1])
    print("Saved to "+outname)
    
    # Save SLAB PIC
    outname = "%sSLAB_PIC_NHFLX_Damping_monwin3_sig005_dof%i_mode%i_%s.npy" % (picout,dofs[0],mode,lagstr)
    np.save(outname,dampingfin[0])
    print("Saved to "+outname)
    
    # Save the masks as well
    outname1 = "%sSLAB_PIC_NHFLX_Damping_monwin3_sig005_dof%i_mode%i_mask_%s.npy" % (picout,dofs[0],mode,lagstr)
    outname0 = "%sFULL_PIC_NHFLX_Damping_monwin3_sig005_dof%i_mode%i_mask_%s.npy" % (picout,dofs[1],mode,lagstr)
    outnames = [outname1,outname0]
    for i in range(2):
        np.save(outnames[i],malls[i])
        print("Saved to "+outnames[i])
        
    
    
    
#%% Select point for comparison
# -----------------------------
lonf = -30
latf = 50
klon,klat = proc.find_latlon(lonf,latf,lon180,lat)
klon360,_ = proc.find_latlon(lonf+360,latf,lon,lat)


#%% Make some plots (Annual Heat Flux without the mask)
# -----------------------------------------------------
bboxplot = [-85,0,0,65]
proj = ccrs.PlateCarree()
fig,axs = plt.subplots(1,2,subplot_kw={'projection':proj})

clvls = np.arange(-40,41,1)

for i in range(2):
    ax = axs[i]
    ax = viz.add_coast_grid(ax,bbox=bboxplot)
    
    plotvar = dunmask[i].mean(-1)
    
    cf = ax.contourf(lon180,lat,plotvar.T,
                     levels=clvls,cmap=cmocean.cm.balance,transform=proj,
                     extend='both')
    ax.contour(lon180,lat,plotvar.T,levels=[0,1],colors='k')
    ax.set_title(mconfigs[i])
fig.colorbar(cf,ax=axs.flatten(),orientation='horizontal')

#%% Try to plot the cross and autocorrelation, monthly
# ----------------------------------------------------
model = 1
lagid = 1
proj = ccrs.PlateCarree()
fig,axs = plt.subplots(4,3,subplot_kw={'projection':proj},figsize=(12,12))
clvl = np.arange(-1,1.05,.05)

for i in tqdm(range(12)):
    ax = axs.flatten()[i]
    ax = viz.add_coast_grid(ax,bbox=bboxplot)
    
    plotvar = rssts[model][i,lagid,:,:]
    
    cf = ax.contourf(lon,lat,plotvar,cmap=cmocean.cm.balance,transform=proj,levels=clvl)
    ax.contour(lon,lat,plotvar,levels=[-1,rthres[model]],colors='k')
    ax.set_title(mons3[i])
plt.suptitle("SST Autocorrelation (Lag %i, %s)" % (lagid+1,mconfigs[model]) ,y=0.90,fontsize=14)
fig.colorbar(cf,ax=axs.flatten(),orientation='vertical',fraction=0.046)

plt.savefig("%sSST_Autocorrelation_%s_Lag%i"%(outpath,mconfigs[model],lagid+1),bbox_inches='tight')

#%

proj = ccrs.PlateCarree()
fig,axs = plt.subplots(4,3,subplot_kw={'projection':proj},figsize=(12,12))
clvl = np.arange(-1,1.05,.05)

for i in tqdm(range(12)):
    ax = axs.flatten()[i]
    ax = viz.add_coast_grid(ax,bbox=bboxplot)
    
    plotvar = rflxs[model][i,lagid,:,:]
    
    cf = ax.contourf(lon,lat,plotvar,cmap=cmocean.cm.balance,transform=proj,levels=clvl)
    ax.contour(lon,lat,plotvar,levels=[-1,rthres[model]],colors='k')
    ax.set_title(mons3[i])

fig.colorbar(cf,ax=axs.flatten(),orientation='vertical',fraction=0.046)

plt.suptitle("SST-FLX Cross-correlation (Lag %i, %s)" % (lagid+1,mconfigs[model]) ,y=0.90,fontsize=14)

plt.savefig("%sSST-FLX_Crosscorrelation_%s_Lag%i"%(outpath,mconfigs[model],lagid+1),bbox_inches='tight')


#%% Check where damping is actually positive (with no mask)
# ---------------------------------------------------------

imod = 1

# Plot Un-masked damping locations
# ----------------------------------------------------
dmp_in = dunmask[imod].copy()
viz.qv_seasonal(lon180,lat,dmp_in,bbox=[-80,0,0,60],anom=True,vmax=50)

# Plot final damping
# ----------------------------------------------------
dmp_in = dampingfin[imod].copy()
viz.qv_seasonal(lon180,lat,dmp_in,bbox=[-80,0,0,60],anom=True,vmax=50)

# Plot locations where the damping values are negative
# ----------------------------------------------------
dmask = (dmp_in<0)
viz.qv_seasonal(lon180,lat,dmask,bbox=[-80,0,0,60],cmap="jet")


# Check which regions fail the test 
# ---------------------------------
viz.qv_seasonal(lon180,lat,malls[imod],bbox=[-80,0,0,60],cmap="jet")

#%% Plot Output Damping for each month with stipling
# --------------------------------------------------

# Select Model
imod     = 0
dmp_in   = dunmask[imod].copy()
sigmask  = malls[imod]

oshift = (lon180[1]-lon180[0])/2
ashift = (lat[1]-lat[0])/2

# Bounidng Boxes, Colorbars,labels
bbox=[-80,0,0,60]
clvl,clb = viz.return_clevels(60, 5,lstep=10)
mons3 = proc.get_monstr(3)

# Plot
fig,axs = plt.subplots(4,3,subplot_kw={'projection':ccrs.PlateCarree()},figsize=(12,14))
imloop = np.roll(np.arange(0,12,1),1)
for i,im in tqdm(enumerate(imloop)):
    
    ax = axs.flatten()[i]
    
    blabel = [0,0,0,0]
    if i%3 == 0:
        blabel[0] = 1
    if i>8:
        blabel[3] = 1
    
    ax = viz.add_coast_grid(ax,bbox=bbox,fill_color="black",blabels=blabel)
    cf = ax.contourf(lon180,lat,dmp_in[:,:,im].T,levels=clvl,cmap=cmocean.cm.balance,extend='both')
    cl = ax.contour(lon180,lat,dmp_in[:,:,im].T,levels=clvl,colors="k",linewidths=0.25)
    ax.clabel(cl,levels=clb,fontsize=10)
    ax.set_title(mons3[im])
    
    # Add Stippling
    smap = viz.plot_mask(lon180,lat,sigmask[:,:,im],ax=ax
                          ,reverse=False,color="k",marker=".",markersize=2)

plt.suptitle("Heat Flux Damping (%s, Lag 1)"% (mconfigs[imod]),y=.92,fontsize=14)
cb = fig.colorbar(cf,ax=axs.flatten(),orientation='horizontal',fraction=0.025,pad=0.05)
cb.set_label("$\lambda _a \, (W m^{-2} \degree C^{-1})$")

plt.savefig("%sHeatFluxDamping_withSig_%s.png" % (outpath,mconfigs[imod]),bbox_inches='tight',dpi=150)



#%% Load data from CESM-FULL Historical

ds = xr.open_dataset("/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/allens_nhflxdamping_monwin3_sig020_dof082_mode4_lag1.nc")

lbdpt = ds.sel(lon=330,lat=50,method='nearest')
dampfull = lbdpt.NHFLX_Damping.values

# Load all data in for basinwide comparison
dampfullall = ds.NHFLX_Damping.values # Lat x Lon360 x Mon x Ens

# Flip the longitude around
nlat,nlon,nmon,nens = dampfullall.shape
dampfullall = dampfullall.reshape(nlat,nlon,nmon*nens).transpose(1,0,2)
_,dampfullall = proc.lon360to180(lon,dampfullall)
dampfullall = dampfullall.reshape(nlon,nlat,nmon,nens)

#%% Plot the comparison at the selected point

fig,ax = plt.subplots(1,1)
for i in range(42):
    ax.plot(dampfull[:,i],label="",color='gray',alpha=0.25)
ax.plot(dampfull[:,-1],label="HTR-FULL (member)",color='gray',alpha=0.25)
ax.plot(dampfull.mean(1),label="HTR-FULL (ens-avg)",color='k',alpha=1)

for i in range(2):
    div = 1
    # if i == 1:
    #     div = 2
    ax.plot(dampingfin[i][klon,klat,:]/div,label=mconfigs[i])

ax.legend()
ax.grid(True,ls='dotted')
ax.set_xlabel("Months")
ax.set_ylabel("$\lambda_a$ (W/m2/degC)")
plt.savefig("%sLambda_Comparison.png"%outpath,dpi=200,bbox_inches='tight')

#%% Plot basinwide comparisons

im = 0

# Set up plotting variable
plotlbd = dampingfin.copy()
dampfull_ensavg = dampfullall.mean(3)
plotlbd.append(dampfull_ensavg)

bboxplot        = [-90,0,0,75]
clevels = np.arange(-30,32,2)
clims = [clevels[0],clevels[-1]]

#
for im in tqdm(range(12)):

    # Plot the comparison
    fig,axs = plt.subplots(1,2,figsize=(10,4),subplot_kw={'projection': ccrs.PlateCarree()})
    for i in range(2):
        ax = axs[i]
        ax.set_title("%s minus HTR-FULL (Ens-Avg)"% (mconfigs[i]))
        ax = viz.add_coast_grid(ax,bbox=bboxplot)
        plotdiff = dampingfin[i][...,im]-dampfull_ensavg[...,im]
        pcm = ax.pcolormesh(lon180,lat,plotdiff.T,cmap='RdBu_r',vmin=clims[0],vmax=clims[-1])
    plt.suptitle("Differences in $\lambda_a$ (W/m2/degC/sec) for %s"% (mons3[im]))
    fig.colorbar(pcm,ax=axs.ravel().tolist())
    plt.savefig("%sLambda_Comparison_%iLags_Mon%02i.png"% (outpath,len(sellags),im+1),dpi=200,bbox_inches='tight')

    
    
    # Plot the actual values
    fig,axs = plt.subplots(1,3,figsize=(15,4),subplot_kw={'projection': ccrs.PlateCarree()})
    for i in range(3):
        ax = axs[i]
        if i == 2:
            ax.set_title("HTR-FULL")
        else:
            ax.set_title(mconfigs[i])
        ax = viz.add_coast_grid(ax,bbox=bboxplot)
        pcm = ax.pcolormesh(lon180,lat,plotlbd[i][...,im].T,cmap='RdBu_r',vmin=clims[0],vmax=clims[-1])
    fig.colorbar(pcm,ax=axs.ravel().tolist())
    plt.suptitle("$\lambda_a$ (W/m2/degC/sec) for %s"% (mons3[im]))
    plt.savefig("%sLambda_Values_%iLags_Mon%02i.png"% (outpath,len(sellags),im+1),dpi=200,bbox_inches='tight')

#%% Plot the annual average

clevels = np.arange(-15,17,2)
clims = [clevels[0],clevels[-1]]

# Plot the comparison
fig,axs = plt.subplots(1,3,figsize=(14,4),subplot_kw={'projection': ccrs.PlateCarree()})

for i in range(2):
    ax = axs[i]
    ax.set_title("%s minus HTR-FULL (Ens-Avg)"% (mconfigs[i]))
    ax = viz.add_coast_grid(ax,bbox=bboxplot)
    plotdiff = dampingfin[i].mean(-1)-dampfull_ensavg.mean(-1)
    
    #pcm = ax.pcolormesh(lon180,lat,plotdiff.T,cmap='RdBu_r',vmin=clims[0],vmax=clims[-1])
    pcm = ax.contourf(lon180,lat,plotdiff.T,cmap='RdBu_r',levels=clevels)
    ax.contour(lon180,lat,plotdiff.T,levels=[0],colors="k",linewidths=0.75)
    
ax = axs[2]
ax.set_title("PIC-SLAB minus PIC-FULL")
ax = viz.add_coast_grid(ax,bbox=bboxplot)
plotdiff = dampingfin[0].mean(-1)-dampingfin[1].mean(-1)
#pcm = ax.pcolormesh(lon180,lat,plotdiff.T,cmap='RdBu_r',vmin=clims[0],vmax=clims[-1])
pcm = ax.contourf(lon180,lat,plotdiff.T,cmap='RdBu_r',levels=clevels)
ax.contour(lon180,lat,plotdiff.T,levels=[0],colors="k",linewidths=0.75)

plt.suptitle("Differences in $\lambda_a$ (W/m2/degC/sec)",fontsize=14,y=0.98)
fig.colorbar(pcm,ax=axs.ravel().tolist())
plt.savefig("%sLambda_Comparison_%iLags_AnnAvg.png"% (outpath,len(sellags)),dpi=200,bbox_inches='tight')

#%% Quickly load MLDs

datpath2 = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/01_Data/model_input/"
h        = np.load(datpath2+"FULL_PIC_HMXL_hclim.npy") # Climatological MLD

#%% Calculate the monthly averages for CESM-SLAB

m = 1
mconfig = mconfigs[m]
print("Plot Mavg for %s"%mconfig)

sids   = [[11,0,1],[2,3,4],[5,6,7],[8,9,10]]
snames = ["Winter (DJF)","Spring (MAM)","Summer (JJA)","Fall (SON)"] 
# Calculate seasonal averages
dampingseas = np.zeros([nlon,nlat,4])*np.nan
hseas = dampingseas.copy()
for i in range (4):
    dampingseas[:,:,i] = dampingfin[m][...,sids[i]].mean(-1)
    hseas[:,:,i] = h[...,sids[i]].mean(-1)



#%% Plot the monthly averages for CESM-SLAB with MLD

clvls = np.arange(-60,65,5)
mldctr = np.arange(0,1200,200)

# Plot damping
fig,axs = plt.subplots(2,2,figsize=(8,6),subplot_kw={'projection': ccrs.PlateCarree()})
for i in range(4):
    
    if i == 0:
        blabel=[1,0,0,0]
    elif i == 1:
        blabel=[0,0,0,0]
    elif i == 2:
        blabel=[1,0,0,1]
    elif i == 3:
        blabel=[0,0,0,1]
    
    ax = axs.flatten()[i]
    ax = viz.add_coast_grid(ax,bbox=bboxplot,blabels=blabel)
    ax.set_title(snames[i])
    
    cf = ax.contourf(lon180,lat,dampingseas[...,i].T,levels=clvls,cmap=cmocean.cm.balance)
    cl = ax.contour(lon180,lat,hseas[...,i].T,levels=mldctr,colors="k",linewidths=.75)
    ax.clabel(cl,levels=[200,1000])
cb = fig.colorbar(cf,ax=axs.flatten(),fraction=0.045)
cb.set_label("$\lambda_a$ ($Wm^{-2}C^{-1}$)")
plt.savefig("%sSeasonal_Damping_CESM-SLAB.png"%outpath,dpi=150,bbox_inches='tight')

#%% 

clvls = np.arange(-60,65,5)
mldctr = np.arange(0,1200,200)

# Plot damping
fig,axs = plt.subplots(2,2,figsize=(8,6),subplot_kw={'projection': ccrs.PlateCarree()})
for i in range(4):
    
    if i == 0:
        blabel=[1,0,0,0]
    elif i == 1:
        blabel=[0,0,0,0]
    elif i == 2:
        blabel=[1,0,0,1]
    elif i == 3:
        blabel=[0,0,0,1]
    
    ax = axs.flatten()[i]
    ax = viz.add_coast_grid(ax,bbox=bboxplot,blabels=blabel)
    ax.set_title(snames[i])
    
    cf = ax.contourf(lon180,lat,dampingseas[...,i].T,levels=clvls,cmap=cmocean.cm.balance)
    cl = ax.contour(lon180,lat,hseas[...,i].T,levels=mldctr,colors="k",linewidths=.75)
    ax.clabel(cl,levels=[200,1000])
cb = fig.colorbar(cf,ax=axs.flatten(),fraction=0.045)
cb.set_label("$\lambda_a$ ($Wm^{-2}C^{-1}$)")
plt.savefig("%sSeasonal_Damping_CESM-FULL.png"%outpath,dpi=150,bbox_inches='tight')


#%% Postprocess the Correlaations

# Determine correlation threshold
corrthress = []
from scipy import stats
ptilde    = 1-p/tails
for dof in dofs:
    critval   = stats.t.ppf(ptilde,dof)
    corrthres = np.sqrt(1/ ((dof/np.power(critval,2))+1))
    corrthress.append(corrthres)


# First preprocess each of the rflx (takem from inside the prep_HF function)
rflxf = []
rsstf = []
for i in range(2):
    rflx = rflxs[i]
    rsst = rssts[i]
    if np.nansum(np.sign(rflx)) < 0:
        print("WARNING! sst-flx correlation is mostly negative, sign will be flipped")
        rflx*=-1
    
    rflxw = scm.postprocess_HF(rflx,limask,sellags,lon)
    rsstw = scm.postprocess_HF(rsst,limask,sellags,lon)
    rflxw *= -1 # Flip back
    rsstw *= -1
    
    rflxf.append(rflxw)
    rsstf.append(rsstw)
#%% Make correlation mask


YY,XX= np.meshgrid(lat,lon180)

rflxmasks = []
rsstmasks = []
for i in range(2):
    
    id_pass = rflxf[i].mean(-1) > corrthress[i]
    rflxmasks.append([YY[id_pass],[XX[id_pass]]])
    
    id_pass = rsstf[i].mean(-1) > corrthress[i]
    rsstmasks.append([YY[id_pass],[XX[id_pass]]])
#%% Plot the correlations

#clims = [0,1]
clvls = np.arange(0,.55,0.05)
# Plot the comparison
fig,axs = plt.subplots(1,2,figsize=(10,4),subplot_kw={'projection': ccrs.PlateCarree()})

for i in range(2):
    ax = axs[i]
    ax.set_title("%s"% (mconfigs[i]))
    ax = viz.add_coast_grid(ax,bbox=bboxplot)
    
    cf = ax.contourf(lon180,lat,rflxf[i].mean(-1).T,cmap='bone_r',levels=clvls)
    
    #ax.scatter(rflxmasks[i][1],rflxmasks[i][0],marker=".",color="k")
    
    #cl = ax.contour(lon180,lat,rflxf[i].mean(-1).T,colors="k",levels=clvls[::2])
    #ax.clabel(cl)


    #pcm =  ax.pcolormesh(lon180,lat,rflxf[i].mean(-1).T,cmap='RdBu_r',vmin=clims[0],vmax=clims[-1])
    
    #plotdiff = dampingfin[i][...,im]-dampfull_ensavg[...,im]
    
    #pcm = ax.pcolormesh(lon180,lat,plotdiff.T,cmap='RdBu_r',vmin=clims[0],vmax=clims[-1])

plt.suptitle("SST-NHFLX Cross Correlation (Lag 1)",fontsize=14,y=0.98)
fig.colorbar(cf,ax=axs.ravel().tolist())
plt.savefig("%sCrossCorr_%iLags.png"% (outpath,len(sellags)),dpi=200,bbox_inches='tight')

#%% PlotAutocorrelations

clevels = np.arange(0.5,1.05,0.05)

fig,axs = plt.subplots(1,2,figsize=(10,4),subplot_kw={'projection': ccrs.PlateCarree()})

for i in range(2):
    ax = axs[i]
    ax.set_title("%s"% (mconfigs[i]))
    ax = viz.add_coast_grid(ax,bbox=bboxplot)
    
    cf = ax.contourf(lon180,lat,rsstf[i].mean(-1).T,cmap='Purples',levels=clevels)
    
    #ax.scatter(rflxmasks[i][1],rflxmasks[i][0],marker=".",color="k")
    
    # for o in range(len(lon)):
    #     for a in range(len(lat)):
    #         if rflxf[i][o,a,:].mean() > corrthress[i]:
    #             ax.scatter(lon180[o],lat[a],marker="x",color="k")
    
    #for o in range(len())
    #for pt in range(len(rflxf[i].flatten())):
        #if rflxf[i].flatten()[pt] > corrthress[i]:
        
        
    
    #cl = ax.contour(lon180,lat,rflxf[i].mean(-1).T,colors="k",levels=clvls[::2])
    #ax.clabel(cl)


    #pcm =  ax.pcolormesh(lon180,lat,rflxf[i].mean(-1).T,cmap='RdBu_r',vmin=clims[0],vmax=clims[-1])
    
    #plotdiff = dampingfin[i][...,im]-dampfull_ensavg[...,im]
    
    #pcm = ax.pcolormesh(lon180,lat,plotdiff.T,cmap='RdBu_r',vmin=clims[0],vmax=clims[-1])

plt.suptitle("SST Autocorrelation (Lag 1)",fontsize=14,y=0.98)
fig.colorbar(cf,ax=axs.ravel().tolist())
plt.savefig("%sAutocorr_%iLags.png"% (outpath,len(sellags)),dpi=200,bbox_inches='tight')


#Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/Data/CESM-SLAB_FULL-Damping/'
#/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data
#%%

def indexwindow(invar,m,monwin,combinetime=False,verbose=False):
    """
    index a specific set of months/years for an odd sliding window
    given the following information (see inputs)
    
    drops the first and last years when a the dec-jan boundary
    is crossed, according to the direction of crossing
    time dimension is thus reduced by 2 overall
    
    inputs:
        1) invar [ARRAY: yr x mon x otherdims] : variable to index
        2) m [int] : index of central month in the window
        3) monwin [int]: total size of moving window of months
        4) combinetime [bool]: set to true to combine mons and years into 1 dimension
    
    output:
        1) varout [ARRAY]
            [yr x mon x otherdims] if combinetime=False
            [time x otherdims] if combinetime=True
    
    """
    
    if monwin > 1:  
        winsize = int(np.floor((monwin-1)/2))
        monid = [m-winsize,m,m+winsize]

    
    varmons = []
    msg = []
    for m in monid:

        if m < 0: # Indexing months from previous year
            
            msg.append("Prev Year")
            varmons.append(invar[:-2,m,:])
            
        elif m > 11: # Indexing months from next year
            msg.append("Next Year")
            varmons.append(invar[2:,m-12,:])
            
        else: # Interior years (drop ends)
            msg.append("Interior Year")
            varmons.append(invar[1:-1,m,:])
    if verbose:
        print("Months are %s with years %s"% (str(monid),str(msg)))       
    # Stack together and combine dims, with time in order
    varout = np.stack(varmons) # [mon x yr x otherdims]
    varout = varout.transpose(1,0,2) # [yr x mon x otherdims]
    if combinetime:
        varout = varout.reshape((varout.shape[0]*varout.shape[1],varout.shape[2])) # combine dims
    return varout


def calc_HF(sst,flx,lags,monwin,verbose=True):
    """
    damping,autocorr,crosscorr=calc_HF(sst,flx,lags,monwin,verbose=True)
    Calculates the heat flux damping given SST and FLX anomalies using the
    formula:
        lambda = [SST(t),FLX(t+l)] / [SST(t),SST(t+l)]
    
    
    Inputs
    ------
        1) sst     : ARRAY [year x  x lat x lon] 
            sea surface temperature anomalies
        2) flx     : ARRAY [time x lat x lon]
            heat flux anomalies
        3) lags    : List of INTs
            lags to calculate for (0-N)
        4) monwin  : INT (odd #)
            Moving window of months centered on target month
            (ex. For Jan, monwin=3 is DJF and monwin=1 = J)
        
        --- OPTIONAL ---
        4) verbose : BOOL
            set to true to display print messages
    
    Outputs
    -------     
        1) damping   : ARRAY [month x lag x lat x lon]
            Heat flux damping values
        2) autocorr  : ARRAY [month x lag x lat x lon]
            SST autocorrelation
        3) crosscorr : ARRAY [month x lag x lat x lon]
            SST-FLX cross correlation
    """
    # Reshape variables [time x lat x lon] --> [yr x mon x space]
    ntime,nlat,nlon = sst.shape
    sst = sst.reshape(int(ntime/12),12,nlat*nlon)
    flx = flx.reshape(sst.shape)
    
    # Preallocate
    nlag = len(lags)
    damping   = np.zeros((12,nlag,nlat*nlon)) # [month, lag, lat, lon]
    autocorr  = np.zeros(damping.shape)
    crosscorr = np.zeros(damping.shape)
    
    st = time.time()
    for l in range(nlag):
        lag = lags[l]
        for m in range(12):
            lm = m-lag # Get Lag Month
            
            # Restrict to time ----
            flxmon = indexwindow(flx,m,monwin,combinetime=True,verbose=False)
            sstmon = indexwindow(sst,m,monwin,combinetime=True,verbose=False)
            sstlag = indexwindow(sst,lm,monwin,combinetime=True,verbose=False)
            
            # Compute Correlation Coefficients ----
            crosscorr[m,l,:] = proc.pearsonr_2d(flxmon,sstlag,0) # [space]
            autocorr[m,l,:] = proc.pearsonr_2d(sstmon,sstlag,0) # [space]
            
            # Calculate covariance ----
            cov     = proc.covariance2d(flxmon,sstlag,0)
            autocov = proc.covariance2d(sstmon,sstlag,0)
            
            # Compute damping
            damping[m,l,:] = cov/autocov
            
            print("Completed Month %02i for Lag %s (t = %.2fs)" % (m+1,lag,time.time()-st))
            
    # Reshape output variables
    damping = damping.reshape(12,nlag,nlat,nlon)  
    autocorr = autocorr.reshape(damping.shape)
    crosscorr = crosscorr.reshape(damping.shape)  
            
    return damping,autocorr,crosscorr

def prep_HF(damping,rsst,rflx,p,tails,dof,mode,returnall=False):
    """
    
    Inputs
    ------
        1) damping   : ARRAY [month x lag x lat x lon]
            Heat flux damping values
        2) autocorr  : ARRAY [month x lag x lat x lon]
            SST autocorrelation
        3) crosscorr : ARRAY [month x lag x lat x lon]
            SST-FLX cross correlation
        4) p : NUMERIC
            p-value
        5) tails : INT
            # of tails for t-test (1 or 2)
        6) dof : INT
            Degrees of freedom
        7) mode: INT
            Apply the following significance testing/masking:
            1 --> No Mask
            2 --> SST autocorrelation based
            3 --> SST-FLX cross correlation based
            4 --> Both 2 and 3
        --- OPTIONAL ---
        8) returnall BOOL
            Set to True to return masks and frequency
    
    Outputs
    -------
        1) dampingmasked [month x lag x lat x lon]
        
    """
    # Determine correlation threshold
    ptilde    = 1-p/tails
    critval   = stats.t.ppf(ptilde,dof)
    corrthres = np.sqrt(1/ ((dof/np.power(critval,2))+1))
    
    # Create Mask
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
    
    # Apply Significance Mask
    dampingmasked = damping * mall
    
    if returnall:
        return dampingmasked,mtot,mult
    return dampingmasked

def postprocess(dampingmasked,limask,sellags,lon):
    
    # Inputs
    ## Dampingmasked [month x lag x lat x lon]
    ## limask [lat x lon]
    
    # Select lags, apply landice mask
    mchoose = dampingmasked[:,sellags,:,:] * limask[None,:,:]
    
    # Flip longiude coordinates ([mon lat lon] --> [lon x lat x mon])
    lon1,dampingw = proc.lon360to180(lon,mchoose.transpose(2,1,0))

    # Multiple by 1 to make positive upwards
    dampingw *= -1
    return dampingw
    
    
    
    
    