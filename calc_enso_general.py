#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Compute ENSO component, and remove it.

Works with output from preproc_ncep.py, but will work to generalize it
This includes the flux data (non-anomalized )

Created on Thu May  5 17:38:51 2022

@author: gliu
"""

#%%

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import time
import sys

#%% Import modules
stormtrack = 0

if stormtrack:
    sys.path.append("/home/glliu/00_Scripts/01_Projects/00_Commons/")
    sys.path.append("/home/glliu/00_Scripts/01_Projects/01_AMV/02_stochmod/stochmod/model/")
else:
    sys.path.append("/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/03_Scripts/stochmod/model/")
    sys.path.append("/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/03_Scripts/")
from amv import proc
import scm
#%% Other edits

datpath =  "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/"
#%%

# Part 2: Preprocess of HFF Estimation ****************************************
# Lets start by doing this just for ts and qnet


# Other Options
debug=False

# Set names
lonname = 'lon'
latname = 'lat'
tname   = 'time'

# Select time crop
tstart  = '1948-01-01'
tend    = '2016-12-31'

# Detrend Method
detrend = 1

# Open dataarray
vnames_in    = ['ts','qnet']
dataset_name ='ncep_ncar'

das = []
for v in vnames_in:
    
    # Open the dataset, slice to time period of interest
    da = xr.open_dataset(datpath+"%s_%s.nc" % (dataset_name,v))
    da = da.sel(time=slice(tstart,tend))
    das.append(da)
    
    # Read out the variables # [time x lat x lon]
    st    = time.time()
    invar = da[v].values
    lon   = da[lonname].values
    lat   = da[latname].values
    times = da[tname].values
    print("Data loaded in %.2fs"%(time.time()-st))
    
    # Remove monthly anomalies
    # ------------------------
    nmon,nlat,nlon = invar.shape
    manom,invar = proc.calc_clim(invar,0,returnts=1) # Calculate clim with time in axis 0
    vanom = invar - manom[None,:,:,:]
    vanom = vanom.reshape(nmon,nlat,nlon) # Reshape back to [time x lat x lon]
    
    # Flip latitude
    if lat[0] > lat[-1]: # If latitude is decreasing...
        lat   = np.flip(lat)
        vanom = np.flip(vanom,axis=1)
    
    # Detrend the variable (take from calc_amv_hadisst.py)
    # ----------------------------------------------------
    start= time.time()
    indata = vanom.reshape(nmon,nlat*nlon).T # Transpose to [Space x Time]
    okdata,knan,okpts = proc.find_nan(indata,1)
    x = np.arange(0,nmon,1)
    if detrend == 0:
        # Compute global weighted average
        glomean = proc.area_avg(vanom.transpose(2,1,0),[0,360,-90,90],lon,lat,1)
        
        # Regress back to the original data to get the global component
        beta,b=proc.regress_2d(glomean,okdata)
        
        # Subtract this from the original data
        okdt = okdata - beta[:,None]
    else:
        # Polynomial Detrend
        okdt,model = proc.detrend_poly(x,okdata,detrend)
        if debug:
            fig,ax=plt.subplots(1,1)
            ax.scatter(x,okdata[44,:],label='raw')
            ax.plot(x,model[44,:],label='fit')
            ax.scatter(x,okdt[:,44],label='dt')
            ax.set_title("Visualize Detrending Method %i"% detrend)
            okdt = okdt.T
            
    data_dt = np.zeros((nmon,nlat*nlon)) * np.nan
    data_dt[:,okpts] = okdt
    data_dt = data_dt.reshape(nmon,nlat,nlon) # Back to [time x lat x lon]
    
    # Save detrended option, if set
    # -----------------------------
    savename = "%s%s_%s_manom_detrend%i.nc" % (datpath,dataset_name,v,detrend)
    da = proc.numpy_to_da(data_dt,times,lat,lon,v,savenetcdf=savename)

    
#%% Part 3, Compute ENSO Indices

# ENSO Parameters
pcrem   = 3 # PCs to calculate
bbox    = [120, 290, -20, 20]
detrend = 1 



# ------------------- -------- General Portion --------------------------------
"""
Looks for <dataset name>_ts_manom_detrend<N>.nc 
containing ts anomalies [time lat lon] with land/ice masked out
"""
st = time.time()

# Open the dataset
savename = "%s%s_ts_manom_detrend%i.nc" % (datpath,dataset_name,detrend)
da = xr.open_dataset(savename)

# Slice to region
da = da.sel(lon=slice(bbox[0],bbox[1]),lat=slice(bbox[2],bbox[3]))

# Read out the variables # [time x lat x lon]
st        = time.time()
invar     = da[v].values
lon       = da[lonname].values
lat       = da[latname].values
times     = da[tname].values
print("Data loaded in %.2fs"%(time.time()-st))

# Portion Below is taken from calc_ENSO_PIC.py VV ***********
eofall,pcall,varexpall = scm.calc_enso(invar,lon,lat,pcrem,bbox=bbox)

# Sanity Check
if debug:
    im = 0
    ip = 0
    proj = ccrs.PlateCarree(central_longitude=180)
    fig,ax = plt.subplots(1,1,subplot_kw={'projection':proj})
    ax = viz.add_coast_grid(ax,bbox=bbox)
    pcm = ax.pcolormesh(lon,lat,eofall[:,:,im,ip],vmin=-1,vmax=1,
                        cmap='cmo.balance',transform=ccrs.PlateCarree())
    cb = fig.colorbar(pcm,ax=ax,orientation='horizontal',fraction=0.055,pad=0.1)
    cb.set_label("SST Anomaly ($\degree C \sigma_{ENSO}^{-1}$)")
    ax.set_title("EOF %i, Month %i\n Variance Explained: %.2f" % (ip+1,im+1,varexpall[im,ip]*100)+"%")

# Save Output
savename = "%senso/%s_ENSO_detrend%i_pcs%i.npz" % (datpath,dataset_name,detrend,pcrem)
np.savez(savename,**{
         'eofs': eofall, # [lon x lat x month x pc]
         'pcs': pcall,   # [Year, Month, PC]
         'varexp': varexpall,
         'lon': lon,
         'lat':lat,
         'times':times,
         'enso_bbox':bbox}
        )
print("Data saved in %.2fs"%(time.time()-st))

#%% Part 4: Calculate the ENSO Component, and remove it.
# (based on remove_ENSO_PIC.py)

ensolag  = 1 # Lag between ENSO month and response month in NATL
reduceyr = True # Drop years due to ENSO lag
monwin   = 3

v        = "ts"

# ------------------- -------- General Portion --------------------------------
"""
Loads enso file containing ensoid
Loads variable file containing anomalized, masked variable

Regresses and removes enso component at specified lag, etc
"""
allstart = time.time()

# Load ENSO
savename = "%senso/%s_ENSO_detrend%i_pcs%i.npz" % (datpath,dataset_name,detrend,pcrem)
ld = np.load(savename,allow_pickle=True)
ensoid = ld['pcs'] # [year x  month x pc]


for v in vnames_in:
    # Load Target variable
    savename = "%s%s_%s_manom_detrend%i.nc" % (datpath,dataset_name,v,detrend)
    da = xr.open_dataset(savename)
    
    # Read out the variables # [time x lat x lon]
    st        = time.time()
    invar     = da[v].values
    lon       = da[lonname].values
    lat       = da[latname].values
    times     = da[tname].values
    
    # Remove ENSO
    vout,ensopattern,times = scm.remove_enso(invar,ensoid,ensolag,monwin,reduceyr=reduceyr,times=times)
    
    # Save output variables
    savename = "%senso/%s_%s_detrend%i_ENSOrem_lag%i_pcs%i_monwin%i.nc" % (datpath,dataset_name,v,detrend,ensolag,pcrem,monwin)
    da = proc.numpy_to_da(vout,times,lat,lon,v,savenetcdf=savename)
    
    # Save ENSO component
    savename = "%senso/%s_%s_detrend%i_ENSOcmp_lag%i_pcs%i_monwin%i.npz" % (datpath,dataset_name,v,detrend,ensolag,pcrem,monwin)
    np.savez(savename,**{
        'ensopattern':ensopattern,
        'lon':lon,
        'lat':lat}
             )
    
    print("Completed variable %s (t=%.2fs)" % (v,time.time()-allstart))


#%% Compute the heat flux feedback

ensorem = False

# Load inputs with variables removed
invars = []
for v in vnames_in:
    
    if ensorem:
        savename = "%senso/%s_%s_detrend%i_ENSOrem_lag%i_pcs%i_monwin%i.nc" % (datpath,dataset_name,v,detrend,ensolag,pcrem,monwin)
    else:
        savename = "%s%s_%s_manom_detrend%i.nc" % (datpath,dataset_name,v,detrend)
    ds       = xr.open_dataset(savename)
    
    loadvar = ds[v].values
    ntime,nlat,nlon = loadvar.shape
    loadvar = loadvar.reshape(int(ntime/12),12,nlat,nlon)
    
    invars.append(loadvar)

#% Calculate heat flux
sst,flx = invars
damping,autocorr,crosscorr = scm.calc_HF(sst,flx,[1,2,3],3,verbose=True,posatm=True)




# Save heat flux
#savename = "%s%s_hfdamping_ensorem%i.npz" % (datpath,dataset_name,ensorem)
#coordsdict = 


# Save 
if debug:
    
    im = 3
    il = 0
    
    proj = ccrs.PlateCarree()
    fig,axs = plt.subplots(4,3,subplot_kw={'projection':proj},
                           figsize=(12,12),constrained_layout=True)
    
    for im in range(12):
        
        ax = axs.flatten()[im]
        plotvar = damping[im,il,:,:]
        lon1,plotvar1 = proc.lon360to180(lon,(plotvar.T)[...,None])
        
        blabel=[0,0,0,0]
        if im%3 == 0:
            blabel[0] = 1
        if im>8:
            blabel[-1] = 1
        
        ax = viz.add_coast_grid(ax,bbox=[-100,0,-10,65],fill_color='gray',
                                blabels=blabel,ignore_error=True)
        pcm = ax.contourf(lon1,lat,plotvar1.squeeze().T*-1,levels = np.arange(-45,50,5),extend='both',
                            cmap='cmo.balance')
    
    cb = fig.colorbar(pcm,ax=axs.flatten(),orientation='horizontal',fraction=0.035,pad=0.05)
    cb.set_label("$\lambda_a$ : $W m^{2} \lambda_a$ ($\degree C ^{-1}$)")
    plt.suptitle("Heat Flux Damping For %s \n Enso Removed: %s | Lag: %i" % (dataset_name,ensorem,il+1))

    
    

# Load sst
sst = np.load()


#%%





# Flip to desired dimensions

# Compute monthly anomalies

# Remove a trend from each point



# Compute ENSO

# Calculate Covariances, Autocovariances

# Compute heat flux feedback