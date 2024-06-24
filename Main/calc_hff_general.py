#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate Heat Flux Feedback (General)

Based on calc_enso_general
Need to eventually support lens datasets...
Currently runs on Astraeus. Will modify for more flexibilityin the future

Inputs


    Region/Time-Cropped Flux and TS [time x (ens) x lat x lon180] : from [preproc_raw_inputs], located in datpath + proc/
    ENSO Index (for ENSO removal): from [compute_enso_index], located in datpath + enso/


Created on Mon Jun 24 14:16:29 2024

@author: gliu
"""


import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from tqdm import tqdm
import xarray as xr
import glob

import time
import sys
import scipy as sp

#%% Import Modules (report with hfutils)

# stormtrack
amvpath = "/home/glliu/00_Scripts/01_Projects/00_Commons/" # amv module
sys.path.append(amvpath)
import amv.proc as hf

# -------------
#%% User Edits
# -------------

# Set Path and options
datpath           = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/"
overwrite         = False # Set to True to ovewrite
debug             = True # Set to true for debugging plots

# Indicate Time Crop (for input)
croptime          = True
tstart            =  '0200-01-01' # "2006-01-01" # 
tend              =  '2000-12-31' #"2101-01-01" # 
timestr           =  '%sto%s'  % (tstart[:4],tend[:4]) # ex. 0000to2000

# Select time crop (for the estimate)
croptime_estimate = False # Cut time right before estimating the heat flux feedback
tcrop_start       = "1970-01-01"#'1920-01-01' '2070-01-01'#
tcrop_end         = "1999-12-31"#'1970-01-01' '2099-12-31'#
tcrop_fname       = ""
if croptime_estimate:
    tcrop_fname      = "_%sto%s" % (tcrop_start.replace('-',''),tcrop_end.replace('-',''))
    
# Indicate bbox crop information
bbox_name         = "NAtl"

# Variables and Dataset Name
vnames_in         = ['TS','qnet'] # ["qnet","fsns","flns","lhflx","shflx"] #"TS" for historical data
dataset_name      = 'cesm2_pic'#'rcp85'
ensnum            = 1

# These should be unnecessary after preproc_raw_inputs
lonname = 'lon'
latname = 'lat'
tname   = 'time' 

# For these datasets, loop for each ensemble member...
lens_datasets     = ['htr','rcp85','gfdl_esm2m_lens','csiro_mk36_lens','canesm2_lens']

# Determine number of ensemble members
if dataset_name == 'rcp85':
    nens = 40
elif dataset_name in ('gfdl_esm2m_lens', "csiro_mk36_lens"):
    nens = 30
elif dataset_name == 'canesm2_lens':
    nens = 50
elif dataset_name == 'htr': # CESM1 Historical
    nens = 42
else:
    nens = 1
    
    
#%% Options for each step

# Step 1 (Preprocessing)
detrend           = 1  # Detrend Method
    
# Step 2 (Remove ENSO)
pcrem             = 3    # PCs calculated
ensolag           = 1    # Lag between ENSO month and response month in NATL
reduceyr          = True # Drop years due to ENSO lag
monwin            = 3    # Window of months to consider

# Step 3 (HFF Calculation)
ensorem           = True # Set to False to skip ENSO removal step

#%% Set up paths

# Set Paths
anompath = datpath + "anom/"
ensopath = datpath + "enso/"
hffpath  = datpath + "hff/"
maskpath = datpath + "masks/"
procpath = datpath + "proc/"

#%%

st_script = time.time()

for ensnum in np.arange(1,nens+1):
#for ensnum in np.arange(1,nens+1):
    
    # -------------------------------------------------------------------------
    # Part 1: Preprocess Variables (Anomalize, Detrend, Flip latitude if needed)
    # -------------------------------------------------------------------------
    
    """
    
    IN  : ncfile, <dataset_name>_<vname>.nc 
        Contains variable [time x lat x lon360] where li-mask has been applied.
        Also contains lat,lon,time.
    
    OUT : ncfile, <dataset_name>_<vname>_manom_detrend#.nc
        Save as above, but anomalized, detrended and latitude corrected.
        
    ex: ncep_ncar_ts.nc  -->  ncep_ncar_ts_manom_detrend1.nc
    
    """
    
    # Set lensflag
    lensflag = False
    if dataset_name in lens_datasets:
        lensflag = True

    for v in vnames_in:
        
        # Load the variable processed by [preproc_raw_inputs], [time x lat x lon180]
        ncname = "%s%s_%s_%s_%s.nc" % (procpath,dataset_name,v,bbox_name,timestr)
        da     = xr.open_dataset(ncname) # [Time x Lat x Lon]
        
        if lensflag:
            da = da.sel(ens=ensnum)
            
        # Below section should already be done
        # # Fix February Start
        # da = hf.fix_febstart(da)
        
        # # Slice to time period of interest
        # if croptime:
        #     da = da.sel(time=slice(tstart,tend),drop=True)
        
        
        # Set Save Name (in /anom/ folder)
        savename_anom = "%s%s_%s_manom_%s_%s_detrend%0i.nc" % (anompath,dataset_name,v,bbox_name,timestr,detrend)
        if lensflag:
            savename = hf.addstrtoext(savename_anom,"_ens%02i"%(ensnum),adjust=-1)
        query = glob.glob(savename_anom)
        
        if (len(query) < 1) or (overwrite == True):
            
            # Read out the other variables # [time x lat x lon]
            # -------------------------------------------------
            st    = time.time()
            invar = da[v].values
            lon   = da[lonname].values
            lat   = da[latname].values
            times = da[tname].values
            
            print("Data loaded in %.2fs"%(time.time()-st))
            
            # For LENs case, remove ensavg
            # ----------------------------
            if lensflag:
                if dataset_name in ["rcp85", "htr"]:
                    eavg_fname  = "%sCESM1_%s_%s_ensAVG.nc" % (datpath,dataset_name,v)
                else:
                    eavg_fname  = "%s%s_%s_ensAVG.nc" % (datpath,dataset_name,v)
                ensavg      = xr.open_dataset(eavg_fname)
                ensavg      = ensavg.sel(time=slice(tstart,tend),drop=True)
                ensavg      = ensavg[v].values
                
                invar = invar - ensavg
                
            # Remove monthly anomalies
            # ------------------------
            nmon,nlat,nlon = invar.shape
            manom,invar = hf.calc_clim(invar,0,returnts=1) # Calculate clim with time in axis 0
            vanom = invar - manom[None,:,:,:]
            vanom = vanom.reshape(nmon,nlat,nlon) # Reshape back to [time x lat x lon]
        
            # Flip latitude
            if lat[0] > lat[-1]: # If latitude is decreasing...
                lat   = np.flip(lat)
                vanom = np.flip(vanom,axis=1)
            
            # Detrend the variable (taken from calc_amv_hadisst.py)
            # ----------------------------------------------------
            if lensflag is False:
                start= time.time()
                indata = vanom.reshape(nmon,nlat*nlon).T # Transpose to [Space x Time]
                okdata,knan,okpts = hf.find_nan(indata,1)
                x = np.arange(0,nmon,1)
                if detrend == 0:
                    # Compute global weighted average
                    glomean = hf.area_avg(vanom.transpose(2,1,0),[0,360,-90,90],lon,lat,1)
                    
                    # Regress back to the original data to get the global component
                    beta,b=hf.regress_2d(glomean,okdata,nanwarn=0)
                    
                    # Subtract this from the original data
                    okdt = okdata - beta[:,None]
                else:
                    # Polynomial Detrend
                    okdt,model = hf.detrend_poly(x,okdata,detrend)
                    if debug:
                        fig,ax=plt.subplots(1,1)
                        ax.scatter(x,okdata[44,:],label='raw')
                        ax.plot(x,model[44,:],label='fit')
                        ax.scatter(x,okdt[:,44],label='dt')
                        ax.set_title("Visualize Detrending Method %i"% detrend)
                        #okdt = okdt.T
                        
                data_dt = np.zeros((nmon,nlat*nlon)) * np.nan
                data_dt[:,okpts] = okdt
                data_dt = data_dt.reshape(nmon,nlat,nlon) # Back to [time x lat x lon]
                
                # Save the output to same path...
                da      = hf.numpy_to_da(data_dt,times,lat,lon,v,savenetcdf=savename_anom)
            else:
                data_dt = vanom
                
                # Save detrended option, if set
                # -----------------------------
                
                da = hf.numpy_to_da(data_dt,times,lat,lon,v,savenetcdf=savename_anom)
            # End Skip
        else:
            print("Skipping. Found existing file: %s" % (str(query)))
    
    # -------------------------------------------------------------------------
    #%% Part 3: Calculate the ENSO Component, and remove it.
    # (based on remove_ENSO_PIC.py)
    
    # ------------------- -------- General Portion --------------------------------
    """
    Loads enso file containing ensoid
    Loads variable file containing anomalized, masked variable
    
    Regresses and removes enso component at specified lag, etc
    """
    allstart = time.time()
    
    # Load ENSO Index
    savename_ensoid = "%s%s_ENSO_detrend%i_pcs%i_%s.npz" % (ensopath,dataset_name,detrend,pcrem,timestr)
    if lensflag:
        savename_ensoid = hf.addstrtoext(savename_ensoid,"_ens%02i"%(ensnum),adjust=0)
    query   = glob.glob(savename_ensoid)
    if len(query) < 1:
        print("%s not found.\nPlease run [compute_enso_index.py]" % savename_ensoid)
    ld      = np.load(savename_ensoid,allow_pickle=True)
    ensoid  = ld['pcs'] # [year x  month x pc]
    
    # Loop by each variable and remove ENSO
    for v in vnames_in:
        
        # Load Target variable
        savename = "%s%s_%s_manom_detrend%i_%s.nc" % (anompath,dataset_name,v,detrend,timestr)
        if lensflag:
            savename = hf.addstrtoext(savename,"_ens%02i"%(ensnum),adjust=-1)
        da = xr.open_dataset(savename)
        
        # Check if (ENSO index file) already exists, and skip if so.
        # ex. cesm2_pic_qnet_NAtl_0200to2000_detrend1_ENSOrem_lag1_pcs3_monwin3.nc
        savename_ensorem = "%s%s_%s_%s_%s_detrend%i_ENSOrem_lag%i_pcs%i_monwin%i.nc" % (anompath,dataset_name,v,
                                                                                   bbox_name,timestr,
                                                                                   detrend,ensolag,pcrem,monwin)
        if lensflag:
            savename_ensorem = proc.addstrtoext(savename_ensorem,"_ens%02i"%(ensnum),adjust=-1)
        query    = glob.glob(savename_ensorem)
        if (len(query) < 1) or (overwrite == True):
            # Read out the variables # [time x lat x lon]
            st        = time.time()
            invar     = da[v].values
            lon       = da[lonname].values
            lat       = da[latname].values
            times     = da[tname].values
            
            # Remove ENSO
            vout,ensopattern,times = hf.remove_enso(invar,ensoid,ensolag,monwin,reduceyr=reduceyr,times=times)
            
            # Save output variables
            #savename = "%senso/%s_%s_detrend%i_ENSOrem_lag%i_pcs%i_monwin%i_%s.nc" % (datpath,dataset_name,v,detrend,ensolag,pcrem,monwin,timestr)
            # if lensflag:
            #     savename = proc.addstrtoext(savename,"_ens%02i"%(ensnum),adjust=-1)
            da = proc.numpy_to_da(vout,times,lat,lon,v,savenetcdf=savename_ensorem)
            
            # Save ENSO component
            savename_ensocomp = "%s%s_%s_detrend%i_ENSOcmp_lag%i_pcs%i_monwin%i_%s.npz" % (ensopath,dataset_name,v,detrend,ensolag,pcrem,monwin,timestr)
            if lensflag:
                savename = proc.addstrtoext(savename_ensocomp,"_ens%02i"%(ensnum),adjust=0)
            np.savez(savename_ensocomp,**{
                'ensopattern':ensopattern,
                'lon':lon,
                'lat':lat}
                     )
            
            print("Completed variable %s (t=%.2fs)" % (v,time.time()-allstart))
            # End Skip
        else:
            print("Skipping. Found existing file: %s" % (str(query)))
    
    # --------------------------------
    #%% Compute the heat flux feedback
    # --------------------------------
    
    # Load inputs with variables removed
    invars = []
    for v in vnames_in:
        
        if ensorem:
            savename_anom_ld = savename_ensorem
        else:
            savename_anom_ld = savename_anom
        if lensflag:
            savename = proc.addstrtoext(savename_anom_ld,"_ens%02i"%(ensnum),adjust=-1)
        ds       = xr.open_dataset(savename_anom_ld)
        
        lat = ds.lat.values
        lon = ds.lon.values
        
        # Crop time if option is set
        if croptime_estimate:
            ds = ds.sel(time=slice(tcrop_start,tcrop_end),drop=True)
        
        loadvar = ds[v].values
        ntime,nlat,nlon = loadvar.shape
        loadvar = loadvar.reshape(int(ntime/12),12,nlat,nlon)
        
        invars.append(loadvar)
    
    #% Calculate heat flux
    sst,flx = invars
    damping,autocorr,crosscorr,autocov,cov = scm.calc_HF(sst,flx,[1,2,3],3,verbose=True,posatm=True,return_cov=True)
    
    # Save heat flux (from hfdamping_mat2nc.py)
    # ----------------------------------------
    outvars  = [damping,crosscorr,autocorr,cov,autocov]
    datpath_out = "%s/useSST/%s_damping/" % (hffpath,v)
    proc.makedir(datpath_out)
    savename = "%s%s_hfdamping_ensorem%i_detrend%i_%s_%scrop.nc" % (datpath_out,dataset_name,ensorem,detrend,timestr,tcrop_fname)
    if lensflag:
        savename = proc.addstrtoext(savename,"_ens%02i"%(ensnum),adjust=-1)
    dims     = {'month'  :np.arange(1,13,1),
                  "lag"  :np.arange(1,4,1),
                  "lat"  :lat,
                  "lon"  :lon}
    
    # Set some attributes
    varnames = ("%s_damping" % v,
                "sst_flx_crosscorr",
                "sst_autocorr",
                "cov",
                "autocov")
    varlnames = ("Net Heat Flux Damping",
                 "SST-Heat Flux Cross Correlation",
                 "SST Autocorrelation",
                 "SST-Heat Flux Covariance",
                 "SST Autocovariance")
    units     = ("W/m2/degC",
                 "Correlation",
                 "Correlation",
                 "W/m2*degC",
                 "degC2")
    
    das = []
    for v,name in enumerate(varnames):

        attr_dict = {'long_name':varlnames[v],
                     'units':units[v]}
        da = xr.DataArray(outvars[v],
                    dims=dims,
                    coords=dims,
                    name = name,
                    attrs=attr_dict
                    )
        if v == 0:
            ds = da.to_dataset() # Convert to dataset
        else:
            ds = ds.merge(da) # Merge other datasets
            
        # Append to list if I want to save separate dataarrays
        das.append(ds)
    
    #% Save as netCDF
    # ---------------
    st = time.time()
    encoding_dict = {name : {'zlib': True} for name in varnames} 
    print("Saving as " + savename)
    ds.to_netcdf(savename,
             encoding=encoding_dict)
    print("Saved in %.2fs" % (time.time()-st))
    
    #%% Save 
    if debug: # Plot seasonal cycle
        il = 0
        proj = ccrs.PlateCarree()
        fig,axs = plt.subplots(4,3,subplot_kw={'projection':proj},
                               figsize=(12,12),constrained_layout=True)
        
        
        plotmon = np.roll(np.arange(0,12),1)
        
        for im in range(12):
            
            monid = plotmon[im]
            
            ax = axs.flatten()[im]
            plotvar = damping[monid,il,:,:]
            lon1,plotvar1 = proc.lon360to180(lon,(plotvar.T)[...,None])
            
            blabel=[0,0,0,0]
            if im%3 == 0:
                blabel[0] = 1
            if im>8:
                blabel[-1] = 1
            
            ax = viz.add_coast_grid(ax,bbox=[-80,0,-10,62],fill_color='gray',
                                    blabels=blabel,ignore_error=True)
            pcm = ax.contourf(lon1,lat,plotvar1.squeeze().T*-1,levels = np.arange(-50,55,5),extend='both',
                                cmap='cmo.balance')
            
            viz.label_sp(monid+1,usenumber=True,alpha=0.7,ax=ax,labelstyle="mon%s")
        
        cb = fig.colorbar(pcm,ax=axs.flatten(),orientation='horizontal',fraction=0.035,pad=0.05)
        cb.set_label("$\lambda_a$ : $W m^{2} \lambda_a$ ($\degree C ^{-1}$)")
        plt.suptitle("Heat Flux Damping For %s \n Enso Removed: %s | Lag: %i" % (dataset_name,ensorem,il+1))
        plt.savefig("%sNHFLX_damping_lag%i_%s_detrend%i_%s.png" % (figpath,il+1,dataset_name,detrend,timestr),dpi=150)
print("Script Ran to Completion in %.2fs"%(time.time()-st_script))
#%%


