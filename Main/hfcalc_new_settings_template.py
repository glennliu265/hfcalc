#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculation Inputs for calc_hff_general_new.py

Created on Thu Jul  3 17:06:23 2025

@author: gliu

"""

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
