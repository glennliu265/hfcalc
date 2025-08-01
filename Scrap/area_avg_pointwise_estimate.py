#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

For a given experiment, take the area-average parameters.
Copied upper section of run_SSS_basinwide

Created on Tue Apr 22 15:15:42 2025

@author: gliu

"""


import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import os

import tqdm
import time

# ----------------------------------
#%% Import custom modules and paths
# ----------------------------------

# Indicate the Machine!
machine = "Astraeus"

# First Load the Parameter File
sys.path.append("../")
import reemergence_params as rparams

# Paths and Load Modules
pathdict   = rparams.machine_paths[machine]

sys.path.append(pathdict['amvpath'])
sys.path.append(pathdict['scmpath'])
from amv import proc,viz
import amv.loaders as dl
import scm
import amv.loaders as dl
import yo_box as ybx

# Set needed paths
input_path  = pathdict['input_path']
output_path = pathdict['output_path']
# procpath    = pathdict['procpath']
# figpath     = pathdict['figpath']
# proc.makedir(figpath)

#%% 


expname     = "SST_Obs_Pilot_00_Tdcorr0_qnet"

expparams   = {
    'varname'           : "SST",
    'bbox_sim'          : [-80,0,20,65],
    'nyrs'              : 1000,
    'runids'            : ["run%02i" % i for i in np.arange(0,10,1)],
    'runid_path'        : None, # If not None, load a runid from another directory
    'Fprime'            : "ERA5_Fprime_QNET_std_pilot.nc",
    'PRECTOT'           : None,
    'LHFLX'             : None,
    'h'                 : "MIMOC_regridERA5_h_pilot.nc",
    'lbd_d'             : None,
    'Sbar'              : None,
    'beta'              : None, # If None, just compute entrainment damping
    'kprev'             : "MIMOC_regridERA5_kprev_pilot.nc",
    'lbd_a'             : "ERA5_qnet_damping_pilot.nc", # NEEDS TO BE CONVERTED TO 1/Mon !!!
    'Qek'               : None, # Now in degC/sec
    'convert_Fprime'    : True,
    'convert_lbd_a'     : True, 
    'convert_PRECTOT'   : False,
    'convert_LHFLX'     : False,
    'froll'             : 0,
    'mroll'             : 0,
    'droll'             : 0,
    'halfmode'          : False,
    "entrain"           : True,
    "eof_forcing"       : False, # CHECK THIS
    "Td_corr"           : False, # Set to True if lbd_d is provided as a correlation, rather than 1/months
    "lbd_e"             : None, # Relevant for SSS
    "Tforce"            : None, # Relevant for SSS
    "correct_Qek"       : False, # Set to True if correction factor to Qek was calculated
    "convert_Qek"       : False, # Set to True if Qek is in W/m2 (True for old SST forcing...) False if in psu/sec or degC/sec (for new scripts)
    }


# Constants
dt    = 3600*24*30 # Timestep [s]
cp    = 3850       # 
rho   = 1026       # Density [kg/m3]
B     = 0.2        # Bowen Ratio, from Frankignoul et al 1998
L     = 2.5e6      # Specific Heat of Evaporation [J/kg], from SSS model document

debug = False

print("==========================")
print("Now Running Experiment: %s" % expname)
print("==========================")

#%% Check and Load Params

print("\tLoading inputs for %s" % expname)

# Apply patch to expdict
expparams = scm.patch_expparams(expparams)

# First, Check if there is EOF-based forcing (remove this if I eventually redo it)
if expparams['eof_forcing']:
    print("\t\tEOF Forcing Detected.")
    eof_flag = True
else:
    print("\t\tEOF Forcing will not be used.")
    eof_flag = False

inputs,inputs_ds,inputs_type,params_vv = scm.load_params(expparams,input_path)



#%% Detect and Process Missing Inputs

_,nlat,nlon=inputs['h'].shape

# Get number of modes
if eof_flag:
    if expparams['varname'] == "SST":
        nmode = inputs['Fprime'].shape[0]
    elif expparams['varname'] == "SSS":
        nmode = inputs['LHFLX'].shape[0]


nr = 0
# Apply roll/shift to seasonal cycle
ninputs = len(inputs)
froll = expparams['froll']
droll = expparams['droll']
mroll = expparams['mroll']
for ni in range(ninputs):
    
    pname = list(inputs.keys())[ni]
    ptype = inputs_type[pname]
    
    if ptype == "mld":
        rollback = mroll
    elif ptype == "forcing":
        rollback = froll
    elif ptype == "damping":
        rollback = droll
    else:
        print("Warning, Parameter Type not Identified. No roll performed.")
        rollback = 0
    
    if rollback != 0:
        print("Rolling %s back by %i" % (pname,rollback))
        
        if eof_flag and len(inputs[pname].shape) > 3:
            rollaxis=1 # Switch to 1st dim to avoid mode dimension
        else:
            rollaxis=0
        inputs[pname] = roll_input(inputs[pname],rollback,axis=rollaxis,halfmode=expparams['halfmode'])

# # Do Unit Conversions ---
inputs_convert  = scm.convert_inputs(expparams,inputs,dt=dt,rho=rho,L=L,cp=cp,return_sep=True)
alpha    = inputs_convert['alpha'] # Amplitude of the forcing
Dconvert = inputs_convert['lbd_a'] # Converted Damping

#%% Save the output in the same style

bbox_spgne  = [-40,-15,52,62] # SPGNE
bbsel       = bbox_spgne
bbname      = "SPGNE"

# Get Damping
mon          = np.arange(1,13,1)
lat          = inputs_ds['lbd_a'].lat
lon          = inputs_ds['lbd_a'].lon
coords_lbda  = dict(mon=mon,lat=lat,lon=lon)
ds_lbda_conv = xr.DataArray(Dconvert,coords=coords_lbda,dims=coords_lbda,name="lbd_a_conv")
ds_lbda      = inputs_ds['lbd_a'].rename('lbd_a')

ds_h  = inputs_ds['h']

# Get Forcing
ds_F     = inputs_ds['Fprime'].squeeze()
ds_Fconv = xr.DataArray(alpha,coords=coords_lbda,dims=coords_lbda,name="Fprime_conv")


inputs_all  = xr.merge([ds_lbda,ds_lbda_conv,ds_h,ds_F,ds_Fconv])
inputs_reg  = proc.sel_region_xr(inputs_all,bbsel)
inputs_aavg = proc.area_avg_cosweight(inputs_reg)

# Save Output (Area Average)
outpath      = "/Users/gliu/Downloads/02_Research/01_Projects/05_SMIO/01_Data/sm_point_run/"
outname_aavg = "%s%s_Method2_inputs.nc" % (outpath,bbname)
edict        = proc.make_encoding_dict(inputs_aavg)
inputs_aavg.to_netcdf(outname_aavg,encoding=edict)

# Save Output (Pointwise)
outname_reg = "%s%s_Pointwise_Parameters.nc" % (outpath,bbname)
edict        = proc.make_encoding_dict(inputs_reg)
inputs_reg.to_netcdf(outname_reg,encoding=edict)


#%% Run Simulation (copied from single_point_estimate)

# =================================
#%% Run the entraining stochastic model
# =================================

nyr_sim  = 1e5
nsim     = 10

simname = "%s_Method2" % bbname
ds_inputs_conv = inputs_aavg

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
outname     = "%s%s_output.nc" % (outpath,simname)
da_sim_out.to_netcdf(outname,encoding=edict)
