#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculation Inputs for calc_hff_general_new.py

Created on Thu Jul  3 17:06:23 2025

@author: gliu

"""

# =============================================================================
# OAFLUX (Qnet, 1948 to 2007 Global)
# =============================================================================
device             = "Astraeus"
datpath            = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/"

# OAFLUX Test Calculation (Qnet, 1948 to 2007, Global) ========================
calcname            = "OAFLUX" 
lonname             = 'lon'
latname             = 'lat'
tname               = 'time' 

# Indicate Input Time Crop (for input)
croptime          =   True
tstart            =  '1984-01-01'
tend              =  '2007-12-31'
timestr           =  '%sto%s'  % (tstart[:4],tend[:4]) # ex. 0000to2000

# Indicate HFF calculation crop
croptime_estimate = False # Cut time right before estimating the heat flux feedback
tstart            =  '1984-01-01'
tend              =  '2007-12-31'
tcrop_fname       = ""
if croptime_estimate:
    tcrop_fname   = "_%sto%s" % (tstart[:4].replace('-',''),tend[:4].replace('-',''))

# Indicate bbox crop information
bbox_name         = "Global" # "NAtl

# Variables Information  -----

# SST 
sstname   = "sst"
sstnc     = "OAFLUX_sst_1984to2007.nc"
sstpath   = "/Users/gliu/Globus_File_Transfer/Reanalysis/OAFLUX/"

# Heat Flux
flxname           = "qnet"
flxnc             = "OAFLUX_qnet_1984to2007.nc"
flxpath           = sstpath
varlnames         = "Net Heat Flux Damping"


# Enso Information
ensonc            = "OAFLUX_ENSO_detrend1_pcs3_1984to2007.nc"



# =============================================================================
#%% CESM2 PiControl, SOM
# =============================================================================

device             = "stormtrack"
datpath            = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/"

# Dataset Info
calcname            = "CESM2_SOM" 
lonname             = 'lon'
latname             = 'lat'
tname               = 'time' 

# Indicate Input Time Crop (for input)
croptime          =   True
tstart            =  '0060-01-01'
tend              =  '0360-12-31'
timestr           =  '%sto%s'  % (tstart[:4],tend[:4]) # ex. 0000to2000

# Indicate HFF calculation crop
croptime_estimate = False # Cut time right before estimating the heat flux feedback
tcrop_start       =  ''
tcrop_end         =  ''
tcrop_fname       =  ""
if croptime_estimate:
    tcrop_fname   = "_%sto%s" % (tcrop_start[:4].replace('-',''),tcrop_end[:4].replace('-',''))

# Indicate bbox crop information
bbox_name         = "NAtl" # "NAtl

# Variables Information  -----

# SST 
sstname           = "TS"
sstnc             = "CESM2_SOM_TS_NAtl_0060to0360.nc"
sstpath           = "/stormtrack/data4/glliu/01_Data/CESM2_PiControl/proc/NAtl/"

# Heat Flux
flxname           = "SHF"
flxnc             = "CESM2_SOM_SHF_NAtl_0060to0360.nc"
flxpath           = sstpath
varlnames         = "Net Heat Flux Damping"


# Enso Information
ensonc            = "CESM2_SOM_ENSO_detrend1_pcs3_0060to0360.nc"

# Additional Information
lensflag          = False#True
nens              = 1 #42


# =============================================================================
#%% CESM2 PiControl, MCOM
# =============================================================================

device             = "stormtrack"
datpath            = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/"

# Dataset Info
calcname            = "CESM2_MCOM" 
lonname             = 'lon'
latname             = 'lat'
tname               = 'time' 

# Indicate Input Time Crop (for input)
croptime          =   True
tstart            =  '0100-01-01'
tend              =  '0500-12-31'
timestr           =  '%sto%s'  % (tstart[:4],tend[:4]) # ex. 0000to2000

# Indicate HFF calculation crop
croptime_estimate = False # Cut time right before estimating the heat flux feedback
tcrop_start       =  ''
tcrop_end         =  ''
tcrop_fname       =  ""
if croptime_estimate:
    tcrop_fname   = "_%sto%s" % (tcrop_start[:4].replace('-',''),tcrop_end[:4].replace('-',''))

# Indicate bbox crop information
bbox_name         = "NAtl" # "NAtl

# Variables Information  -----

# SST 
sstname           = "TS"
sstnc             = "CESM2_MCOM_TS_NAtl_0100to0500.nc"
sstpath           = "/stormtrack/data4/glliu/01_Data/CESM2_PiControl/proc/NAtl/"

# Heat Flux
flxname           = "SHF"
flxnc             = "CESM2_MCOM_SHF_NAtl_0100to0500.nc"
flxpath           = sstpath
varlnames         = "Net Heat Flux Damping"


# Enso Information
ensonc            = "CESM2_MCOM_ENSO_detrend1_pcs3_0100to0500.nc"

# Additional Information
lensflag          = False#True
nens              = 1 #42


# =============================================================================
#%% CESM2 PiControl, FOM
# =============================================================================

device             = "stormtrack"
datpath            = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/output/"

# Dataset Info
calcname            = "CESM2_FOM" 
lonname             = 'lon'
latname             = 'lat'
tname               = 'time' 

# Indicate Input Time Crop (for input)
croptime          =   True
tstart            =  '0200-01-01'
tend              =  '2000-12-31'
timestr           =  '%sto%s'  % (tstart[:4],tend[:4]) # ex. 0000to2000

# Indicate HFF calculation crop
croptime_estimate = False # Cut time right before estimating the heat flux feedback
tcrop_start       =  ''
tcrop_end         =  ''
tcrop_fname       =  ""
if croptime_estimate:
    tcrop_fname   = "_%sto%s" % (tcrop_start[:4].replace('-',''),tcrop_end[:4].replace('-',''))

# Indicate bbox crop information
bbox_name         = "NAtl" # "NAtl

# Variables Information  -----

# SST 
sstname           = "TS"
sstnc             = "CESM2_FOM_TS_NAtl_0200to2000.nc"
sstpath           = "/stormtrack/data4/glliu/01_Data/CESM2_PiControl/proc/NAtl/"

# Heat Flux
flxname           = "SHF"
flxnc             = "CESM2_FOM_SHF_NAtl_0200to2000.nc"
flxpath           = sstpath
varlnames         = "Net Heat Flux Damping"


# Enso Information
ensonc            = "CESM2_FOM_ENSO_detrend1_pcs3_0200to2000.nc"

# Additional Information
lensflag          = False#True
nens              = 1 #42