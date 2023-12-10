#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute ENSO component, and remove it.
Also compute the heat flux feedback.

Works with output from preproc_ncep.py, but will work to generalize it
This includes the flux data (non-anomalized )

Plots:
    - Plots for each month for a given simulation

Created on Thu May  5 17:38:51 2022

@author: gliu
"""

#%%

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import time
import sys
import cartopy.crs as ccrs
import glob

#%% Import modules
stormtrack = 1
if stormtrack:
    sys.path.append("/home/glliu/00_Scripts/01_Projects/00_Commons/")
    sys.path.append("/home/glliu/00_Scripts/01_Projects/01_AMV/02_stochmod/stochmod/model/")
    
    # Path to the processed dataset (qnet and ts fields, full, time x lat x lon)
    datpath =  "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_RCP85/01_PREPROC/"
    figpath =  "/home/glliu/02_Figures/01_WeeklyMeetings/20220622"
else:
    sys.path.append("/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/03_Scripts/stochmod/model/")
    sys.path.append("/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/03_Scripts/")

    # Path to the processed dataset (qnet and ts fields, full, time x lat x lon)
    datpath =  "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/"
    figpath =  "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/02_Figures/20220511/"
from amv import proc,viz
import scm


#%% User Settings

proc.makedir(figpath)

# Part 1 (Preprocessing) ------------------------------------------
overwrite         = False # Skip the file if it already exists

# Select time crop (prior to preprocessing)
croptime          = True # Cut the time prior to detrending, EOF, etc
tstart            =  '1920-01-01' # "2006-01-01" # 
tend              =  '2006-01-01' #"2101-01-01" # 

# Select time crop (for the estimate)
croptime_estimate = True # Cut time right before estimating the heat flux feedback
tcrop_start       = "1970-01-01"#'1920-01-01' '2070-01-01'#
tcrop_end         = "1999-12-31"#'1970-01-01' '2099-12-31'#
tcrop_fname       = ""
if croptime_estimate:
    tcrop_fname      = "_%sto%s" % (tcrop_start.replace('-',''),tcrop_end.replace('-',''))

# Detrend Method
detrend           = 1 

# Variables and Dataset Name
vnames_in         = ['TS','FLNS'] # ["qnet","fsns","flns","lhflx","shflx"] #"TS" for historical data
dataset_name      = 'htr'#'rcp85'
ensnum            = 1

lens_datasets     = ['htr','rcp85','gfdl_esm2m_lens','csiro_mk36_lens','canesm2_lens']
#"csiro_mk36_lens"
#'CESM1_FULL_PIC'
#'ncep_ncar'
#'era20c'

# Determine number of ensemble members
if dataset_name == 'rcp85':
    nens = 40
elif dataset_name in ('gfdl_esm2m_lens', "csiro_mk36_lens"):
    nens = 30
elif dataset_name == 'canesm2_lens':
    nens = 50
elif dataset_name == 'htr':
    nens = 42
else:
    nens = 1

# Set coordinate names for dataset
lonname  = 'lon'
latname  = 'lat'
tname    = 'time'

# Part 2 (ENSO Index Calculation) ---------------------------------

# ENSO Parameters
pcrem    = 3                   # PCs to calculate
bbox     = [120, 290, -20, 20] # ENSO Bounding Box

# Part 3 (ENSO Removal) -------------------------------------------
ensolag  = 1    # Lag between ENSO month and response month in NATL
reduceyr = True # Drop years due to ENSO lag
monwin   = 3    # Window of months to consider


# Part 4 (HFF Calculations) -------------------------------------------
ensorem  = True

# Toggles
debug    = False # Print Figures, statements for debugging

#%% Main Body

st_script = time.time()

for ensnum in np.arange(1,nens+1):
#for ensnum in np.arange(1,nens+1):
    
    # Part 1: Preprocess Variables (Anomalize, Detrend, Flip latitude if needed)
    
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
            
        # Open the dataset, slice to time period of interest
        if lensflag:
            if dataset_name in ["rcp85", "htr"]:
                if dataset_name == "rcp85":
                    datpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_RCP85/01_PREPROC/"
                elif dataset_name == "htr":
                    datpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_HTR/"
                da = xr.open_dataset("%sCESM1_%s_%s_ens%02i.nc" % (datpath,dataset_name,v,ensnum))
            else:
                datpath = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_lens/%s/" % dataset_name
                da = xr.open_dataset("%s%s_%s_ens%02i.nc" % (datpath,dataset_name,v,ensnum))
        else:
            da = xr.open_dataset(datpath+"%s_%s.nc" % (dataset_name,v))
        
        if croptime:
            da = da.sel(time=slice(tstart,tend),drop=True)
        
        # Check time, and skip file if it already exists
        # ----------------------------------------------
        times   = da[tname].values
        timesyr = times.astype('datetime64[Y]').astype(int) +1970
        timestr = "%ito%i" % (timesyr[0],timesyr[-1])
        # Set Save Name
        savename = "%s%s_%s_manom_detrend%i_%s.nc" % (datpath,dataset_name,v,detrend,timestr)
        if lensflag:
            savename = proc.addstrtoext(savename,"_ens%02i"%(ensnum),adjust=-1)
        query = glob.glob(savename)
        if (len(query) < 1) or (overwrite == True):
            
            # Read out the other variables # [time x lat x lon]
            # -------------------------------------------------
            st    = time.time()
            invar = da[v].values
            lon   = da[lonname].values
            lat   = da[latname].values
            
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
            manom,invar = proc.calc_clim(invar,0,returnts=1) # Calculate clim with time in axis 0
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
                okdata,knan,okpts = proc.find_nan(indata,1)
                x = np.arange(0,nmon,1)
                if detrend == 0:
                    # Compute global weighted average
                    glomean = proc.area_avg(vanom.transpose(2,1,0),[0,360,-90,90],lon,lat,1)
                    
                    # Regress back to the original data to get the global component
                    beta,b=proc.regress_2d(glomean,okdata,nanwarn=0)
                    
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
                        #okdt = okdt.T
                        
                data_dt = np.zeros((nmon,nlat*nlon)) * np.nan
                data_dt[:,okpts] = okdt
                data_dt = data_dt.reshape(nmon,nlat,nlon) # Back to [time x lat x lon]
            else:
                data_dt = vanom
                
                # Save detrended option, if set
                # -----------------------------
                
                da = proc.numpy_to_da(data_dt,times,lat,lon,v,savenetcdf=savename)
            # End Skip
        else:
            print("Skipping. Found existing file: %s" % (str(query)))
        
    #%% Part 2, Compute ENSO Indices
    
    # ------------------- -------- General Portion --------------------------------
    
    """
    
    IN : ncfile, <dataset_name>_<vname>_manom_detrend#.nc
        Anomalized, detrended ts with landice masked applied
    
    OUT : npz file <dataset_name>_ENSO_detrend#_pcs#.npz
        PC File containing:
            eofall (ENSO EOF Patterns)          [lon x lat x month x pc]
            pcall  (ENSO principle components)  [time x month x pc]
            varexpall (ENSO variance explained) [month x pc]]
            lon,lat,time,ensobbox variables
    
    ex: ncep_ncar_ts_manom_detrend1.nc --> ncep_ncar_ENSO_detrend1_pcs3.npz
    
    """
    
    st = time.time()
    
    # Open the dataset
    savename = "%s%s_%s_manom_detrend%i_%s.nc" % (datpath,dataset_name,vnames_in[0],detrend,timestr)
    if lensflag:
        savename = proc.addstrtoext(savename,"_ens%02i"%(ensnum),adjust=-1)
    da = xr.open_dataset(savename)
    
    # Slice to region
    da = da.sel(lon=slice(bbox[0],bbox[1]),lat=slice(bbox[2],bbox[3]))
    
    
    # Check if ENSO has already been calculated and skip if so
    proc.makedir("%senso/"% datpath) 
    savename = "%senso/%s_ENSO_detrend%i_pcs%i_%s.npz" % (datpath,dataset_name,detrend,pcrem,timestr)
    if lensflag:
        savename = proc.addstrtoext(savename,"_ens%02i"%(ensnum),adjust=0)
    query = glob.glob(savename)
    if (len(query) < 1) or (overwrite == True):
        # Read out the variables # [time x lat x lon]
        st        = time.time()
        invar     = da[vnames_in[0]].values
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
    else:
        print("Skipping. Found existing file: %s" % (str(query)))
    # End Skip
        
    #%% Part 3: Calculate the ENSO Component, and remove it.
    # (based on remove_ENSO_PIC.py)
    
    # ------------------- -------- General Portion --------------------------------
    """
    Loads enso file containing ensoid
    Loads variable file containing anomalized, masked variable
    
    Regresses and removes enso component at specified lag, etc
    """
    allstart = time.time()
    
    # Load ENSO
    savename = "%senso/%s_ENSO_detrend%i_pcs%i_%s.npz" % (datpath,dataset_name,detrend,pcrem,timestr)
    if lensflag:
        savename = proc.addstrtoext(savename,"_ens%02i"%(ensnum),adjust=0)
    ld = np.load(savename,allow_pickle=True)
    ensoid = ld['pcs'] # [year x  month x pc]
    
    for v in vnames_in:
        
        # Load Target variable
        savename = "%s%s_%s_manom_detrend%i_%s.nc" % (datpath,dataset_name,v,detrend,timestr)
        if lensflag:
            savename = proc.addstrtoext(savename,"_ens%02i"%(ensnum),adjust=-1)
        da = xr.open_dataset(savename)
        
        # Check if (ENSO index file) already exists, and skip if so.
        savename = "%senso/%s_%s_detrend%i_ENSOrem_lag%i_pcs%i_monwin%i_%s.nc" % (datpath,dataset_name,v,detrend,ensolag,pcrem,monwin,timestr)
        if lensflag:
            savename = proc.addstrtoext(savename,"_ens%02i"%(ensnum),adjust=-1)
        query    = glob.glob(savename)
        if (len(query) < 1) or (overwrite == True):
            # Read out the variables # [time x lat x lon]
            st        = time.time()
            invar     = da[v].values
            lon       = da[lonname].values
            lat       = da[latname].values
            times     = da[tname].values
            
            # Remove ENSO
            vout,ensopattern,times = scm.remove_enso(invar,ensoid,ensolag,monwin,reduceyr=reduceyr,times=times)
            
            # Save output variables
            #savename = "%senso/%s_%s_detrend%i_ENSOrem_lag%i_pcs%i_monwin%i_%s.nc" % (datpath,dataset_name,v,detrend,ensolag,pcrem,monwin,timestr)
            # if lensflag:
            #     savename = proc.addstrtoext(savename,"_ens%02i"%(ensnum),adjust=-1)
            da = proc.numpy_to_da(vout,times,lat,lon,v,savenetcdf=savename)
            
            # Save ENSO component
            savename = "%senso/%s_%s_detrend%i_ENSOcmp_lag%i_pcs%i_monwin%i_%s.npz" % (datpath,dataset_name,v,detrend,ensolag,pcrem,monwin,timestr)
            if lensflag:
                savename = proc.addstrtoext(savename,"_ens%02i"%(ensnum),adjust=0)
            np.savez(savename,**{
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
            savename = "%senso/%s_%s_detrend%i_ENSOrem_lag%i_pcs%i_monwin%i_%s.nc" % (datpath,dataset_name,v,detrend,ensolag,pcrem,monwin,timestr)
        else:
            savename = "%s%s_%s_manom_detrend%i_%s.nc" % (datpath,dataset_name,v,detrend,timestr)
        if lensflag:
            savename = proc.addstrtoext(savename,"_ens%02i"%(ensnum),adjust=-1)
        ds       = xr.open_dataset(savename)
        
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
    datpath_out = "%s%s_damping/" % (datpath,v)
    proc.makedir(datpath_out)
    savename = "%s%s_hfdamping_ensorem%i_detrend%i_%s_%scrop.nc" % (datpath_out,dataset_name,ensorem,detrend,timestr,tcrop_fname)
    if lensflag:
        savename = proc.addstrtoext(savename,"_ens%02i"%(ensnum),adjust=-1)
    dims     = {'month':np.arange(1,13,1),
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

# Flip to desired dimensions

# Compute monthly anomalies

# Remove a trend from each point

# Compute ENSO

# Calculate Covariances, Autocovariances

# Compute heat flux feedback
