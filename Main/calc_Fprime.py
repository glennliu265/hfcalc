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

Eprime     = False # Set to True to Compute E' instead of F'

stormtrack = 0

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
datname = "cesm2_pic"

# Mixed Layer Depth --> [h: time x lat x lon180]
mldpath      = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/03_reemergence/proc/model_input/mld/"
mldnc        = "cesm2_pic_HMXL_NAtl_0200to2000.nc"
mldname      = "h"

# Net Heat Flux (Positive Upwards) --> [qnet: time x lat x lon180]
flxpath      = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/proc/"
flxnc        = "cesm2_pic_qnet_NAtl_0200to2000.nc"
flxname      = 'qnet'

# SST
sstpath      = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/proc/"
sstnc        = "cesm2_pic_TS_NAtl_0200to2000.nc"
sstname      = 'TS'

# Damping 
damppath     = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/hff/qnet_damping/"
dampnc       = "cesm2_pic_hfdamping_NAtl_0200to2000_ensorem1_detrend1.nc"
dampname     = 'qnet_damping'
ilag         = 0 # Indicate which lag to select

# Damping Information and roll options
dampstr     = "CESM2PiCqnetDamp"
nroll       = 0 # Amount to roll lbd*T' term
rollstr     = "nroll%0i"  % nroll
convert_wm2 = True # Convert hff to wm2

# Conversion Factors
dt  = 3600*24*30
cp0 = 3996
rho = 1026

lensflag    = True # Set to True to detrend with ensemble average

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
    ds_dt    = [ds-ds.mean('ensemble') for ds in ds_anom] # [ens x time x lat x lon]
    ds_dt    = [ds.transpose('time','ensemble','lat','lon') for ds in ds_dt]
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
    
# Load mixed layer depth for conversion from [calc_hclim.py] # [mon x lat x lon]
ds_mld   = xr.open_dataset(mldpath + mldnc)[mldname]


# Double Check the Sizes
ds_in   = ds_dt + [dshff,ds_mld]
ds_in   = hf.resize_ds(ds_in)
sst,qnet,dshff,ds_mld = ds_in

# Convert HFF (1/mon to W/m2 per degC) if needed
if convert_wm2:
    print("Converting to Wm2")
    dshff = dshff * (rho*cp0*ds_mld) / dt  *-1 #need to do a check for - value!!
else:
    dshff= dshff

# Load output to numpy
hff     = dshff.values
sst     = sst.values
qnet    = qnet.values

#%%

#%%
# Indicate Search String for qnet/SST files ------d
ncstr1   = "CESM1LE_%s_NAtl_19200101_20050101_bilinear.nc"

# Indicate Mixed-Layer Depth File
mldnc    = "%sCESM1_HTR_FULL_HMXL_NAtl.nc" % mldpath

# Fprime Calculation Options
nroll    = 0


# Damping Options ----------
dampstr = "LHFLXnomasklag1" # Damping String  (see below, "load damping of choice")
"""
Current List of Damping Strings
-- Name             -- ncfile                                               -- Description
"nomasklag1"        "CESM1_HTR_FULL_qnet_damping_nomasklag1.nc"             Default Qnet Damping as calculated from covariance-based method.
"Expfitlbda123"     "CESM1_HTR_FULL_Expfit_lbda_damping_lagsfit123.nc"      Exp Fit to SST - Expfit to SSS; Mean of Lags 1,2,3
"ExpfitSST123"      "CESM1_HTR_FULL_Expfit_SST_damping_lagsfit123.nc"       Exp Fit to SST (total); Mean of Lags 1,2,3
"LHFLXnomasklag1"   "CESM1_HTR_FULL_LHFLX_damping_nomasklag1_EnsAvg.nc"     Default LHFLX Damping as calculated from covariance-based method
"""
if dampstr == "Expfitlbda123":
    convert_wm2=True
    hff_nc   = "CESM1_HTR_FULL_Expfit_lbda_damping_lagsfit123.nc"
elif dampstr == "nomasklag1":
    convert_wm2=False
    hff_nc = "CESM1_HTR_FULL_qnet_damping_nomasklag1.nc"
elif dampstr == "ExpfitSST123":
    convert_wm2=True
    hff_nc   = "CESM1_HTR_FULL_Expfit_SST_damping_lagsfit123.nc"#"CESM1_HTR_FULL_qnet_damping_nomasklag1.nc"
elif dampstr == "LHFLXnomasklag1":
    convert_wm2= False
    hff_nc   = "CESM1_HTR_FULL_LHFLX_damping_nomasklag1.nc"

else:
    print("Invalid dampstr, currently not supported...")
    
# Conversion Factors
dt  = 3600*24*30
cp0 = 3996
rho = 1026



# -----------------------------------------------------------------------------
#%% Part 2: Load and Convert Damping
# -----------------------------------------------------------------------------


# Load HFF
dshff    = xr.open_dataset(damppath+dampnc)[dampname] # [mon x ens x lat x lon]

# Load mixed layer depth for conversion
ds_mld   = xr.open_dataset(mldnc)[mldname]

# Check sizes, make sure they are all the same...
# if dampstr is not None: # Not sure why, but it seems that the hff default is wrongly cropped
#     ds_list = ds_dt + [dshff,ds_mld]
#     ds_rsz  = proc.resize_ds(ds_list)
#     ds_dt = ds_rsz[:2]
#     dshff = ds_rsz[2]
#     ds_mld = ds_rsz[3]

# Convert HFF (1/mon to W/m2 per degC) if needed
if convert_wm2:

    dshff = dshff * (rho*cp0*ds_mld.h) / dt  * -1 #need to do a check for - value!!
else:
    dshff= dshff.damping

# Load output to numpy
hff     = dshff.values
sst     = ds_dt[0].SST.values
qnet    = ds_dt[1][flxname].values

# -----------------------------------------------------------------------------
#%% Part 3: Tile heat flux feedback and make Fprime 
# -----------------------------------------------------------------------------
ntime,nens,nlat,nlon        = qnet.shape # Check sizes and get dimensions for tiling
ntimeh,nensh,nlath,nlonh    = hff.shape
nyrs                        = int(ntime/12)
hfftile                     = np.tile(hff.transpose(1,2,3,0),nyrs)
hfftile                     = hfftile.transpose(3,0,1,2)
# Check plt.pcolormesh(hfftile[0,0,:,:]-hfftile[12,0,:,:]),plt.colorbar(),plt.show()

#% Calculate F'
Fprime       = qnet - hfftile*np.roll(sst,nroll) # Minus is the correct way to go
#Fprime_minus = qnet - hfftile*np.roll(sst,nroll)

# -----------------------------------------------------------------------------
#%% Part 4: Save Fprime output (full timeseries) (Optional)
# -----------------------------------------------------------------------------
if Eprime:
    outvar = "LHFLX"
else:
    outvar = "Fprime"

coords   = dict(time=ds_dt[0].time.values,ens=dshff.ens.values,lat=dshff.lat.values,lon=dshff.lon.values)
daf      = xr.DataArray(Fprime,coords=coords,dims=coords,name=outvar)
savename = "%sCESM1_HTR_FULL_%s_timeseries_%s_%s_NAtl.nc" % (rawpath1,"Eprime",dampstr,rollstr)
edict    = {outvar:{'zlib':True}}
daf.to_netcdf(savename,encoding=edict)
print("Script ran to completion in %.2fs" % (time.time()-st))

# -----------------------------------------------------------------------------
#%% Part 5. Check Whitening at a point
# -----------------------------------------------------------------------------

#% Here's a debug section to check the power spectra and whitening using Fprime
lonf  = -30
latf  = 50
nroll = 0

# Get SST and Qnet
dspt = [proc.selpt_ds(ds,lonf,latf,) for ds in ds_dt]
sstpt,qnetpt= dspt[0].SST.values,dspt[1][flxname].values
hffpt       = proc.selpt_ds(dshff,lonf,latf)

# Tile and Compute
nyrs        = int(len(ds_dt[0].time)/12)
hffpttile     = np.array([np.tile(hffpt.values[:,e][:],nyrs).flatten() for e in range(nens)]).T


Fp        = qnetpt + hffpttile*np.roll(sstpt,nroll)
Fpminus   = qnetpt - hffpttile*np.roll(sstpt,nroll)

Fpminusr1 = qnetpt - hffpttile*np.roll(sstpt,-2,axis=0)

#%% Look at spectra

ints    = [sstpt,qnetpt,Fp,Fpminus,Fpminusr1]
labs    = ["SST","Qnet","Fprime Plus","Fprime Minus","Fprime Minus Roll 1"]
#output  = scm.quick_spectra(ints,)

nsmooth  = 5
pct      = 0.1
dtin     = 3600*24*30

specvars = []
nvars    = len(ints)
for vv in range(nvars):
    tsens    = ints[vv]
    tsens    = [tsens[:,e] for e in range(nens)]
    specout  = scm.quick_spectrum(tsens, nsmooth, pct, dt=dtin,return_dict=True,
                                 make_arr=True)
    specvars.append(specout)
    
#%% PLot the ensemble mean
dtplot  = dtin
fig,axs = plt.subplots(nvars,1,constrained_layout=True,figsize=(12,10))

for vv in range(nvars):
    ax = axs[vv]
    plotspec = specvars[vv]['specs'].mean(0)/dtplot
    plotfreq = specvars[vv]['freqs'].mean(0) * dtplot
    ax.plot(plotfreq,plotspec,label=labs[vv])
    
    ax.axhline(np.var(ints[vv])/(plotfreq[-1] - plotfreq[0]))
    ax.legend()
    
    
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
    