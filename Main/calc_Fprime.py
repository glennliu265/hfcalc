#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 10:19:50 2024

calculate_Fprime.py
========================

Copied from calc_Fprime_lens.py
Works with output processed by [preproc_raw_inputs]


Given Damping, MLD, SST , and Qnet, compute Fprime where:
    Qnet = F' + lbd*T
So:
    F' = Qnet - lbd*T 
    
Does the same for E', where Qnet is replaced by QLHFLX
    E' = Qlhflx - lbd*T'

Where T and Qnet are not anomalized. Written to run on stormtracl currently


This script will be used by NHFLX_EOF_monthly.

Inputs:
------------------------

    varname : dims                              - units                 - processing script
    SST     : (ensemble, time, lat, lon)        [degC]                  ????
    qnet    : (ensemble, time, lat, lon)        [W/m2]                  
    h       : (month, ens, lat, lon)            [meters]                calc_hclim.py, prep_mld_PIC.py
    damping : (mon, ens, lat, lon)              [degC/W/m2] OR [1/mon]  calc_enso_general.py
    
Outputs: 
------------------------

    varname : dims                              - units 
    Fprime  : (time, ens, lat, lon)             [W/m2]

Output File Name: "%sCESM1_HTR_FULL_Fprime_timeseries_%s_%s_NAtl.nc" % (rawpath1,dampstr,rollstr)

What does this script do?
------------------------
    1) Load, deseasonalize, and detrend qnet/TS
    2) Load (and optionally convert) HFF
    3) Tile HFF and compute Fprime
    4) Save output

Script History
------------------------
 - Moved from NHFLX_EOF_monthly_lens.py on 2024.02.12
 - Copied Fprime calculation step from preproc_sm_inputs_SSS
 - Created on Mon Feb 12 15:40:17 2024

@author: gliu
"""


import numpy as np
import xarray as xr
import sys
import time
import matplotlib.pyplot as plt

#%%
# stormtrack
amvpath = "/home/glliu/00_Scripts/01_Projects/00_Commons/" # amv module
sys.path.append(amvpath)
import amv.proc as hf



#%% Set Paths

# Note: For Eprime computation, just switch to the correct flux and damping!
#Eprime     = False # Set to True to Compute E' instead of F'

stormtrack   = 1

# # Path to variables processed by prep_data_byvariable_monthly, Output will be saved to rawpath1
# if stormtrack:
#     rawpath1 = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/03_reemergence/proc/CESM1/NATL_proc/"
#     dpath    = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/03_reemergence/proc/model_input/damping/"
#     mldpath  = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/03_reemergence/proc/model_input/mld/"
# else:
#     rawpath1 = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/03_reemergence/01_Data/proc/CESM1/NATL_proc/"
#     mldpath  = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/03_reemergence/01_Data/proc/model_input/mld/"
#     dpath    = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/03_reemergence/01_Data/proc/model_input/damping/"

# Indicate inputs
datname      = "cesm1le_htr_5degbilinear"
lensflag     = True # Set to True for lens datasets/to detrend with ensemble average
outvar       = "Fprime"  # "Set to Fprime by default, but LHFLX for Eprime calculations..."
outpath      = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/proc/"

# Mixed Layer Depth --> [h: time x lat x lon180]
mldpath      = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/03_reemergence/proc/model_input/mld/"
mldnc        = "cesm1_htr_5degbilinear_HMXL_Global_1920to2005.nc"#"cesm1_htr_5degbilinear_HMXL_Global_1920to2005.nc"
mldname      = "h"

# Net Heat Flux (Positive Upwards) --> [qnet: time x lat x lon180]
flxpath      = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/proc/"
flxnc        = "cesm1_htr_5degbilinear_qnet_Global_1920to2005.nc"
flxname      = 'qnet'

# SST
sstpath      = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/proc/"
sstnc        = "cesm1_htr_5degbilinear_TS_Global_1920to2005.nc"
sstname      = 'TS'

# Damping 
damppath     = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/hff/qnet_damping/"
dampnc       = "cesm1_htr_5degbilinear_hfdamping_Global_1920to2005_ensorem1_detrend1.nc"
dampname     = 'qnet_damping'
ilag         = 0 # Indicate which lag to select

# Damping Information and roll options
dampstr      = "cesm1le5degqnet"
nroll        = 0 # Amount to roll lbd*T' term
rollstr      = "nroll%0i"  % nroll
convert_wm2  = False # Convert hff to wm2

# Conversion Factors
dt          = 3600*24*30
cp0         = 3996
rho         = 1026


#%% Functions
def format_ds_mon(ds):
    if "ensemble" in list(ds.dims):
        print("Renaming 'ensemble' --> ens")
        ds = ds.rename({'ensemble':'ens'})
    if "month" in list(ds.dims):
        print("Renaming 'month' --> mon")
        ds = ds.rename({'month':'mon'})
    return ds

# -----------------------------------------------------------------------------
#%% Part 1: Load, Deseasonalize, Detrend qnet and SST
# -----------------------------------------------------------------------------
# Note this was copied from preproc_sm_inputs_SSS.py
st       = time.time()

# Load TS, flux and preprocess -------------------------
# if Eprime:
#     print("Loading LHFLX to compute E'")
#     flxname = "LHFLX"
# else:
#     print("Loading Q_net to compute F'")
#     flxname = "qnet"

# Load SST and Flux
st      = time.time()
ds_sst  = xr.open_dataset(sstpath+sstnc)[sstname].load() # Load SST
ds_flx  = xr.open_dataset(flxpath+flxnc)[flxname].load() # Load qnet
print("Loaded Flux and SST in %.2fs" % (time.time()-st))

# Preprocess 
ds_load = [ds_sst,ds_flx]
# if ds_sst.shape != ds_flx.shape:
#     print("Resizing variables")
#     ds_load = hf.resize_ds(ds_load) # (make sure they are the same size)

# Anomalize
ds_anom  = [hf.xrdeseason(ds) for ds in ds_load]

# Detrend
if lensflag:
    print("Detrending by removing ensemble mean")
    ds_dt    = [ds-ds.mean('ens') for ds in ds_anom] # [ens x time x lat x lon]
    ds_dt    = [ds.transpose('time','ens','lat','lon') for ds in ds_dt]
else:
    print("Applying Simple Linear Detrend")
    ds_dt    = [hf.xrdetrend(ds) for ds in ds_anom]
    ds_dt    = [ds.transpose('time','lat','lon') for ds in ds_dt]
    ds_dt    = [format_ds_mon(ds) for ds in ds_dt]

# -----------------------------------------------------------------------------
#%% Part 2: Load Damping/MLD and Convert HFF
# -----------------------------------------------------------------------------

# Load HFF from [calc_hff_general.py]
dshff    = xr.open_dataset(damppath + dampnc)[dampname]      # [mon x (lag) x (ens) x lat x lon]
if "lag" in list(dshff.dims):
    dshff = dshff.isel(lag=ilag)
    print("Selecting lag %i for heat flux feedback" % dshff.lag)


if lensflag:
    dshff = dshff.transpose('month','ens','lat','lon')

# Load mixed layer depth for conversion from [calc_hclim.py] # [mon x lat x lon]
ds_mld                  = xr.open_dataset(mldpath + mldnc)[mldname]


dshff,ds_mld            = [format_ds_mon(ds) for ds in [dshff,ds_mld]]

# Double Check the Sizes
ds_in                   = ds_dt + [dshff,ds_mld]
ds_in                   = hf.resize_ds(ds_in)
sst,qnet,dshff,ds_mld   = ds_in

# Convert HFF (1/mon to W/m2 per degC) if needed
if convert_wm2:
    print("Converting to Wm2")
    dshff = dshff * (rho*cp0*ds_mld) / dt  * -1 #need to do a check for - value!!
else:
    dshff = dshff

# Load output to numpy
hff     = dshff.values
sst     = sst.values
qnet    = qnet.values

if lensflag is False: # Add singleton dimension for ensemble
    qnet    = qnet[:,None,:,:]
    hff     = hff[:,None,:,:]
    sst     = sst[:,None,:,:]

# -----------------------------------------------------------------------------
#%% Part 3: Tile heat flux feedback and make Fprime 
# -----------------------------------------------------------------------------
# Get Dimensions
ntime,nens,nlat,nlon        = qnet.shape # Check sizes and get dimensions for tiling
ntimeh,nensh,nlath,nlonh    = hff.shape

nyrs                        = int(ntime/12)

hfftile                     = np.tile(hff.transpose(1,2,3,0),nyrs) # [nens,nlat,nlon,ntime]
hfftile                     = hfftile.transpose(3,0,1,2) #   # [nens,nlat,nlon,ntime]

# Check plt.pcolormesh(hfftile[0,0,:,:]-hfftile[12,0,:,:]),plt.colorbar(),plt.show()

#% Calculate F'
Fprime       = qnet - hfftile*np.roll(sst,nroll,axis=0) # Minus is the correct way to go
#Fprime_minus = qnet - hfftile*np.roll(sst,nroll)

# -----------------------------------------------------------------------------
#%% Part 4: Save Fprime output (full timeseries) (Optional)
# -----------------------------------------------------------------------------
if lensflag:
    coords   = dict(time=ds_dt[0].time.values,ens=dshff.ens.values,lat=dshff.lat.values,lon=dshff.lon.values)
else:
    coords   = dict(time=ds_dt[0].time.values,lat=dshff.lat.values,lon=dshff.lon.values)
daf      = xr.DataArray(Fprime.squeeze(),coords=coords,dims=coords,name=outvar)

savename = "%s%s_%s_timeseries_%s_%s_NAtl.nc" % (outpath,datname,outvar,dampstr,rollstr)

edict    = {outvar:{'zlib':True}}
daf.to_netcdf(savename,encoding=edict)
print("Script ran to completion in %.2fs" % (time.time()-st))



#%% Debugging stuff
debug = True

if debug:
    import yo_box as ybx
    from amv import proc
    scmpath = "/home/glliu/00_Scripts/01_Projects/01_AMV/02_stochmod/stochmod/model/" # scm module
    sys.path.append(scmpath)
    import scm
    
    lonf = -30
    latf = 50
    
    # Get SST and Qnet
    dspt        = [proc.selpt_ds(ds,lonf,latf,) for ds in ds_dt]
    sstpt,qnetpt= dspt[0].values,dspt[1].values
    hffpt       = proc.selpt_ds(dshff,lonf,latf) 
    
    # Tile and Compute
    nyrs          = int(len(ds_dt[0].time)/12)
    if lensflag:
        hffpttile     = np.array([np.tile(hffpt.values[:,e][:],nyrs).flatten() for e in range(nens)]).T
    else:
        hffpttile     = np.tile(hffpt,nyrs)
        
    hffpttile = hffpttile * -1
    
    
    Fp        = qnetpt + hffpttile*np.roll(sstpt,nroll)
    Fpminus   = qnetpt - hffpttile*np.roll(sstpt,nroll)
    
    Fpminusr1 = qnetpt - hffpttile*np.roll(sstpt,-2,axis=0)
    
    
    #%% Look at spectra
    
    # original test
    #ints    = [sstpt,qnetpt,Fp,Fpminus,Fpminusr1]
    #labs    = ["SST","Qnet","Fprime Plus","Fprime Minus","Fprime Minus Roll 1"]
    # calculated test
    ints     = [daf.sel(lon=lonf,lat=latf,method='nearest').data,
                ds_dt[1].sel(lon=lonf,lat=latf,method='nearest').data]
    
    labs     = ['Fprime',"qnet"]
    
    
    
    #output  = scm.quick_spectra(ints,)
    
    nsmooth  = 150
    pct      = 0.1
    dtin     = 3600*24*30
    
    specvars = []
    nvars    = len(ints)
    for vv in range(nvars):
        tsens    = ints[vv]
        if lensflag:
            tsens    = [tsens[:,e] for e in range(nens)]
        else:
            tsens    = [tsens]
        specout  = scm.quick_spectrum(tsens, nsmooth, pct, dt=dtin,return_dict=True,
                                     make_arr=True)
        specvars.append(specout)
    
    #%% PLot the ensemble mean
    dtplot  = dtin
    fig,axs = plt.subplots(nvars,1,constrained_layout=True,figsize=(12,10))
    
    for vv in range(nvars):
        ax = axs[0]
        if lensflag:
            plotspec = specvars[vv]['specs'].mean(0)/dtplot
            plotfreq = specvars[vv]['freqs'].mean(0) * dtplot
        else:
            plotspec = specvars[vv]['specs'].squeeze()/dtplot
            plotfreq = specvars[vv]['freqs'].squeeze()* dtplot
        
        ax.plot(plotfreq,plotspec,label=labs[vv])
        
        ax.axhline(np.var(ints[vv])/(plotfreq[-1] - plotfreq[0]))
        ax.legend()
        
    plt.show()
    
    
    
    
    


# #%%

# #%%
# # Indicate Search String for qnet/SST files ------d
# ncstr1   = "CESM1LE_%s_NAtl_19200101_20050101_bilinear.nc"

# # Indicate Mixed-Layer Depth File
# mldnc    = "%sCESM1_HTR_FULL_HMXL_NAtl.nc" % mldpath

# # Fprime Calculation Options
# nroll    = 0


# # Damping Options ----------
# dampstr = "LHFLXnomasklag1" # Damping String  (see below, "load damping of choice")
# """
# Current List of Damping Strings
# -- Name             -- ncfile                                               -- Description
# "nomasklag1"        "CESM1_HTR_FULL_qnet_damping_nomasklag1.nc"             Default Qnet Damping as calculated from covariance-based method.
# "Expfitlbda123"     "CESM1_HTR_FULL_Expfit_lbda_damping_lagsfit123.nc"      Exp Fit to SST - Expfit to SSS; Mean of Lags 1,2,3
# "ExpfitSST123"      "CESM1_HTR_FULL_Expfit_SST_damping_lagsfit123.nc"       Exp Fit to SST (total); Mean of Lags 1,2,3
# "LHFLXnomasklag1"   "CESM1_HTR_FULL_LHFLX_damping_nomasklag1_EnsAvg.nc"     Default LHFLX Damping as calculated from covariance-based method
# """
# if dampstr == "Expfitlbda123":
#     convert_wm2=True
#     hff_nc   = "CESM1_HTR_FULL_Expfit_lbda_damping_lagsfit123.nc"
# elif dampstr == "nomasklag1":
#     convert_wm2=False
#     hff_nc = "CESM1_HTR_FULL_qnet_damping_nomasklag1.nc"
# elif dampstr == "ExpfitSST123":
#     convert_wm2=True
#     hff_nc   = "CESM1_HTR_FULL_Expfit_SST_damping_lagsfit123.nc"#"CESM1_HTR_FULL_qnet_damping_nomasklag1.nc"
# elif dampstr == "LHFLXnomasklag1":
#     convert_wm2= False
#     hff_nc   = "CESM1_HTR_FULL_LHFLX_damping_nomasklag1.nc"

# else:
#     print("Invalid dampstr, currently not supported...")
    
# # Conversion Factors
# dt  = 3600*24*30
# cp0 = 3996
# rho = 1026




# -----------------------------------------------------------------------------
#%% Part 5. Check Whitening at a point
# -----------------------------------------------------------------------------

if debug:
    
    
    
    #% Here's a debug section to check the power spectra and whitening using Fprime
    lonf  = -30
    latf  = 50
    nroll = 0
    
    # Get SST and Qnet
    dspt        = [proc.selpt_ds(ds,lonf,latf,) for ds in ds_dt]
    sstpt,qnetpt= dspt[0].SST.values,dspt[1][flxname].values
    hffpt       = proc.selpt_ds(dshff,lonf,latf)
    
    # Tile and Compute
    nyrs        = int(len(ds_dt[0].time)/12)
    hffpttile     = np.array([np.tile(hffpt.values[:,e][:],nyrs).flatten() for e in range(nens)]).T
    
    
    Fp        = qnetpt + hffpttile*np.roll(sstpt,nroll)
    Fpminus   = qnetpt - hffpttile*np.roll(sstpt,nroll)
    
    Fpminusr1 = qnetpt - hffpttile*np.roll(sstpt,-2,axis=0)


        
    
#%% Look at ACF

tsmetrics = []
nvars    = len(ints)
for vv in range(nvars):
    tsens    = ints[vv]
    tsens    = [tsens[:,e] for e in range(nens)]
    tsm = scm.compute_sm_metrics(tsens)
    tsmetrics.append(tsm)
    
    
#%% Look at ACF

fig,axs = plt.subplots(nvars,1,constrained_layout=True,figsize=(12,10))
lags = np.arange(37)
for vv in range(nvars):
    
    ax,_ = viz.init_acplot(7,lags,lags,ax=ax)
    ax = axs[vv]
    plotspec = np.array(tsmetrics[vv]['acfs'][7]).mean(0)
    plotfreq = np.arange(37)
    ax.plot(plotfreq,plotspec,label=labs[vv])
    
    #ax.axhline(np.var(ints[vv])/(plotfreq[-1] - plotfreq[0]))
    ax.legend()
    