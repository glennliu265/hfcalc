#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Estimate Stochastic Model Parameters at a single point (or area average)

Inputs:
    SST
    FLX
    ENSO Index
    MLD



Created on Mon Mar  3 15:49:05 2025

@author: gliu

"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import time
import sys
import cartopy.crs as ccrs
import glob

import pandas as pd

import scipy as sp

#%% Import modules
stormtrack = 0
if stormtrack:
    sys.path.append("/home/glliu/00_Scripts/01_Projects/00_Commons/")
    sys.path.append("/home/glliu/00_Scripts/01_Projects/01_AMV/02_stochmod/stochmod/model/")
    
    # Path to the processed dataset (qnet and ts fields, full, time x lat x lon)
    #datpath =  "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_RCP85/01_PREPROC/"
    datpath =  "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/anom/"
    figpath =  "/home/glliu/02_Figures/01_WeeklyMeetings/20240621/"
    
else:
    sys.path.append("/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/03_Scripts/stochmod/model/")
    sys.path.append("/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/03_Scripts/")

    # Path to the processed dataset (qnet and ts fields, full, time x lat x lon)
    datpath =  "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/"
    figpath =  "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/02_Figures/20220511/"
from amv import proc,viz
import scm


import amv.proc as hf # Update hf with actual hfutils script, most relevant functions
import amv.loaders as dl

#%% Load everything



dpath       = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/reanalysis/proc/NATL_proc_obs/"
nc_sst      = dpath + "ERA5_sst_NAtl_1979to2021.nc" #proc/"
nc_thflx    = dpath + "ERA5_thflx_NAtl_1979to2021.nc"
nc_enso     = dpath + "enso/ERA5_ENSO_detrend1_pcs3_1979to2021.nc"
nc_mld      = dpath + "MIMOC_RegridERA5_mld_NAtl_Climatology.nc"


ds_sst      = xr.open_dataset(nc_sst).load()
ds_thflx    = xr.open_dataset(nc_thflx).load()
ds_enso     = xr.open_dataset(nc_enso).load()
ds_mld      = xr.open_dataset(nc_mld).load()

ds_all      = xr.merge([ds_sst,ds_thflx,ds_enso,ds_mld])

#%%

# Limit to region and take area average
bbox_yeager = [-50,-10,50,60] # Original
ds_all_reg  = proc.sel_region_xr(ds_all,bbox_yeager)
ds_regavg   = proc.area_avg_cosweight(ds_all_reg)
dtmon       = 60*60*24

#%% Step (1): Compute HFF


#%% 1.A Preprocess

# Deseason
ds_vars = xr.merge([ds_regavg.sst,ds_regavg.thflx])
dsa     = proc.xrdeseason(ds_vars)

# Detrend
order  = 3
sst    = dsa.sst.data
flx    = dsa.thflx.data
times  = np.arange(len(sst))
invars = [sst,flx]
outdt = [proc.polyfit_1d(times,iv,order) for iv in invars]

invars_detrended = [ii[2] for ii in outdt ] # 3rd argument

# Visualize Detrend
fig,axs = plt.subplots(2,1,constrained_layout=True,figsize=(12,4.5))
ax      = axs[0]
ax.plot(times,sst,color="k",lw=.75,label="Raw")
ax.plot(times,outdt[0][1],color="orange",lw=.75,label="%i-order fit"%order)
ax.plot(times,outdt[0][2],color="blue",lw=.75,label="Detrended")

ax      = axs[1]
ax.plot(times,flx,color="k",lw=.75,label="Raw")
ax.plot(times,outdt[1][1],color="orange",lw=.75,label="%i-order fit"%order)
ax.plot(times,outdt[1][2],color="blue",lw=.75,label="Detrended",ls='dashed')

#%% 1.B Regress out ENSO (skip this step)

ensoid = ds_all.pcs.data # [Yr x Mon x PC]
ensolag = 1
monwin  = 3

# vout,ensopattern
enso_removed = [scm.remove_enso(vv[:,None,None],ensoid,ensolag,monwin) for vv in invars_detrended]
invar_noenso = [enso_removed[ii][0] for ii in range(2)]

#%% 1.C Compute the HFF

#damping,autocorr,crosscorr,autocov,cov = scm.calc_HF(sst,flx,[1,2,3],3,verbose=True,posatm=True,return_cov=True)

nyrs                  = int(invar_noenso[0].shape[0]/12)
invar_noenso_rs       = [vv.reshape(nyrs,12,1,1) for vv in invar_noenso]
sst_noenso,flx_noenso = invar_noenso_rs
# damping,autocorr,crosscorr,autocov,cov
hff_out               = scm.calc_HF(sst_noenso,flx_noenso,[1,2,3],3,verbose=True,posatm=True,
                                    return_cov=True,return_dict=True)

dt                    = 60*60*24
damping               = hff_out['damping'] / dt * -1


lbd_a = damping[:,0,:,:].squeeze() # Just Take Lag 1


# ------------ <0> ------------
#%% Check the Significance
# ------------ <0> ------------

hff     = damping
rsst    = hff_out['autocorr']
rflx    = hff_out['crosscorr']
p       = 0.05
tails   = 2
dof     = (2021-1979 + 1 - 2 - 1) * 3
method  = 4

# Compute and Apply Mask
st = time.time()
dampingmasked, freq_success, sigmask = scm.prep_HF(hff, rsst, rflx,
                                                   p, tails, dof, method,
                                                   returnall=True)  # expects, [month x lag x lat x lon], should generalized with ensemble dimension?

_,_,acmask = scm.prep_HF(hff, rsst, rflx,
                        p, tails, dof, 2,
                        returnall=True)  # expects, [month x lag x lat x lon], should generalized with ensemble dimension?

_,_,ccmask = scm.prep_HF(hff, rsst, rflx,
                        p, tails, dof, 3,
                        returnall=True)  # expects, [month x lag x lat x lon], should generalized with ensemble dimension?


print("Completed significance testing in %.2fs" % (time.time()-st))
#%%

mons3    = proc.get_monstr()

fsz_axis = 14
fig,axs  = viz.init_monplot(3,3,constrained_layout=True,figsize=(14,6.5))

for ii in range(3):
    
    if ii == 0:
        
        vname = "Autocorrelation"
        invar = rsst.squeeze()
        inmask = acmask.squeeze() #  Mon x Lag
        ylim  = [0.2,1]
        
    elif ii == 1:
        
        vname = "Crosscorrelation"
        invar = rflx.squeeze()
        inmask = ccmask.squeeze()
        ylim  = [-0.2,0.5]
        
    else:
        
        vname = "Heat Flux Feedback"
        invar  = hff.squeeze()
        inmask = sigmask.squeeze()
        ylim  = [-35,35]
        
    for ll in range(3):
        
        ax = axs[ii,ll]
        
        
        if ii == 0:
            
            ax.set_title("Lag %i" % (ll + 1),fontsize=fsz_axis)
            
        if ll == 0:
            
            viz.add_ylabel(vname,ax=ax,fontsize=fsz_axis)
            
        plotvar = invar[:,ll]
        plotmsk = inmask[:,ll]
        #plotmsk = np.where(np.isnan(plotmsk),0,1)
        #plotmsk[np.isnan(plotmsk)] = 0
        ax.set_ylim(ylim)
        #ax.axhline([0],ls='dashed',color='k',lw=0.75)
        
        
        ax.plot(mons3,plotvar)
        ax.plot(mons3,plotvar*plotmsk,marker="o",ls='none',markersize=10)
        

# =================================
#%% 2. Compute Fprime
# =================================
# Get Flx and SST
flx = invars_detrended[1]/dt * -1 # Time
sst = invars_detrended[0] # 
mld = ds_regavg.mld

# Tile Damping
#lbd_a_conv = lbd_a * mld 
lbd_a_tile = np.tile(lbd_a,int(len(flx)/12)) #* -1

# Compute Fprime

"""

Qnet = F' - lbd*T' 

Rationale: 
    - Qnet and F' are both downwards positive (as in ERA5)
    - A given +T' should lead to flux upwards (or negative, reducing F')

F' = Qnet + lbd*T'

"""

nroll  = 0 
fprime = flx + lbd_a_tile * np.roll(sst,nroll)


#% Get amplitude of Fprime
fprime_amp = fprime.reshape(int(len(flx)/12),12).std(0)

# ------------ <0> ------------
#%% Debug by checking spectra
# ------------ <0> ------------
"""
Notes: 
    
    - it seems that the spectra is a bit whiter when using nroll = -1


"""
input_spec = [sst,flx,fprime]

nsmooth = [50,]*3
dtmon   = 60*60*24*30
outspec = scm.quick_spectrum(input_spec,nsmooth,0.10,dt=dtmon,
                             return_dict=True)

vnames = ["SST","FLX","Fprime"]
fig,ax = plt.subplots(1,1,figsize=(12,3.5))

for ii in range(3):
    
    
    plotfreq = outspec['freqs'][ii] * dtmon
    plotspec = outspec['specs'][ii] / dtmon
    
    ax.plot(plotfreq,plotspec,label=vnames[ii])
    
ax.legend()


# =================================
#%% Store stochastic model parameters
# =================================

coordsmon = {'mon':np.arange(1,13)}

da_fprime = xr.DataArray(fprime_amp,name="Fprime",
                         dims=coordsmon,coords=coordsmon)
da_lbda   = xr.DataArray(lbd_a,name="lbd_a",
                         dims=coordsmon,coords=coordsmon)
da_mld    = xr.DataArray(mld,name='h',
                         dims=coordsmon,coords=coordsmon)

ds_inputs = xr.merge([da_fprime,da_lbda,da_mld])

def convert_inputs_simple(ds_inputs,rho=1026,cp=3850):
    # For monthly 1D inputs (lbd_a,Fprime,h)...
    
    # Compute Conversion Factor
    ndays = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
    dtmon = ndays * 60*60*24
    conversion_factor = 1 / (rho * cp * ds_inputs.h) * dtmon
    
    # Convert [W/m2/degC] to [1/mon], or [W/m2] to [degC/mon]
    ds_inputs['lbd_a_conv']     = ds_inputs['lbd_a'] * conversion_factor
    ds_inputs['Fprime_conv']    = ds_inputs['Fprime'] * conversion_factor
    
    return ds_inputs


ds_inputs_conv = convert_inputs_simple(ds_inputs)

# =================================
#%% Run the entraining stochastic model
# =================================

nyr_sim  = 1e5
nsim     = 10

def run_entrain_1d_xr(ds_inputs_conv,noise_ts):
    
    h_in        = ds_inputs_conv.h.data[None,None,:]
    lbd_a       = ds_inputs_conv.lbd_a_conv.data[None,None,:]
    F           = (np.tile(ds_inputs_conv.Fprime_conv.data,int(len(noise_ts)/12)) * noise_ts)[None,None,:]
    kprev       = scm.find_kprev(ds_inputs_conv.h.data,returnh=False)[None,None,:]
    
    entrain_out = scm.integrate_entrain(h_in,kprev,lbd_a,F,T0=None,return_dict=True)
    
    return entrain_out

out_ssts = np.zeros((nsim,int(nyr_sim)*12)) * np.nan
for nn in range(nsim):

    noise_ts = np.random.normal(0,1,int(nyr_sim*12))
    entrain_out = run_entrain_1d_xr(ds_inputs_conv,noise_ts)
    sst_sim     = entrain_out['T'].squeeze()
    
    out_ssts[nn,:] = sst_sim.copy()
    

times_sim   = xr.cftime_range(start='0000',periods=int(nyr_sim*12),freq="MS",calendar="noleap")
dims        = {'ens':np.arange(1,nsim+1,1) ,'time':times_sim}
da_sim_out  = xr.DataArray(out_ssts,dims=dims,coords=dims,name='sst')
edict       = proc.make_encoding_dict(da_sim_out)
outpath_spg = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/spg_data/"
outname     = outpath_spg + "stochastic_model_point.nc"
da_sim_out.to_netcdf(outname,encoding=edict)


zd
#%% Lets compare the SSTs by looking at a few Metrics
# ------------ <0> ------------
#proc.

sst = out_ssts[0,:]
in_ssts = [sst,proc.deseason(sst_sim)]

expmarkers      = ['d',"o"]
expcolors       = ['k','goldenrod']
#expcolors_sig   = ['gray','']
expnames    = ['ERA5 SST','Stochastic Model']

# Compute the ACF
nsmooth     = [10,500]
metrics_out = scm.compute_sm_metrics(in_ssts,nsmooth=nsmooth)

# ------------ <0> ------------    
#%% Plot the ACF
# ------------ <0> ------------

lags   = np.arange(37)
xtks   = np.arange(0,39,3)
kmonth = 1

fig,ax = plt.subplots(1,1,figsize=(8,4.5),constrained_layout=True)
ax,_   = viz.init_acplot(kmonth,xtks,lags,ax=ax,title="")

ax.plot(lags,metrics_out['acfs'][kmonth][0],label="ERA5 SST",lw=2.5,marker=expmarkers[0],c=expcolors[0])
ax.plot(lags,metrics_out['acfs'][kmonth][1],label="Stochastic Model",lw=2.5,marker=expmarkers[1],c=expcolors[1])
ax.legend()

ax.set_ylabel("Correlation with %s SST Anomalies" % (mons3[kmonth]))
ax.set_xlabel("Lag from %s [Months]" % (mons3[kmonth]))
#ax.set_title("SST Autocorrelation (%s)" % (mons3[kmonth]),fontsize=20)




#%% Plot the Monthly Variance

mons3  = proc.get_monstr()
fig,ax = viz.init_monplot(1,1,figsize=(6,3.5))

ax.plot(mons3,metrics_out['monvars'][0],label="ERA5 SST",lw=2.5,marker=expmarkers[0],c=expcolors[0])
ax.plot(mons3,metrics_out['monvars'][1],label="Stochastic Model",lw=2.5,marker=expmarkers[1],c=expcolors[1])
ax.legend()

ax.set_ylim([0.1,0.4])
ax.set_ylabel("Monthly Variance [$\degree$C$^2$]")

#%% Plot the Spectra (need to recheck how to do this)

xper        = np.array([40,10,5,1,0.5])
xper_ticks  = 1 / (xper*12)

dtmon_fix   = 60*60*24*30
fig,ax      = plt.subplots(1,1,figsize=(8,4.5),constrained_layout=True)

for ii in range(2):
    plotspec = metrics_out['specs'][ii] / dtmon_fix
    plotfreq = metrics_out['freqs'][ii] * dtmon_fix
    
    CCs = metrics_out['CCs'][ii] / dtmon_fix
    
    
    ax.loglog(plotfreq,plotspec,lw=.75,
            c=expcolors[ii],label=expnames[ii])
    
    
    ax.loglog(plotfreq,CCs[:,1],c=expcolors[ii],ls='dashed',lw=0.55)

ax.set_xlim([1/1000,0.5])
ax.axvline([1/(6)],label="",ls='dotted',c='gray')
ax.axvline([1/(12)],label="",ls='dotted',c='gray')
ax.axvline([1/(5*12)],label="",ls='dotted',c='gray')
ax.axvline([1/(10*12)],label="",ls='dotted',c='gray')
ax.axvline([1/(40*12)],label="",ls='dotted',c='gray')
ax.legend()
ax.set_xlabel("Frequency (1/Month)",fontsize=14)
ax.set_ylabel("Power [$\degree C ^2 cycle \, per \, mon$]")

ax2 = ax.twiny()
ax2.set_xlim([1/1000,0.5])
ax2.set_xscale('log')
ax2.set_xticks(xper_ticks,labels=xper)
ax2.set_xlabel("Period (Years)",fontsize=14)

