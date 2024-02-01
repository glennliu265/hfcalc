#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Heat Flux Feedback Calculation Parameters

Common parameters used in heat flux feedback calculations, as well 
as shortform names for preprocessing choices.

Created on Thu Feb  1 13:44:49 2024

@author: gliu
"""

import numpy as np

#%% HFF Significance Testing / Preprocessing Names

# Default run for stochastic model
hname1  = "smpaper" # Run Used for Stochastic Model Paper (PiControl, Slab Replacement)
hparam1 = {
    'ensorem' : 1,      # 1=enso removed, 0=not removed
    'ensolag' : 1,      # Lag Applied toENSO and Variable before removal
    'monwin'  : 3,      # Size of month window for HFF calculations
    'detrend' : 1,      # Whether or not variable was detrended
    'tails'   : 2,      # tails for t-test
    'p'       : 0.05,   # p-value for significance testing
    'sellags' : [0,],   # Lags included (indices, so 0=lag1)
    'lagstr'  : "lag1", # Name of lag based on sellags
    'method'  : 5       # Significance test option: 1 (No Mask); 2 (SST autocorr); 3 (SST-FLX crosscorr); 4 (Both), 5 (Replace with SLAB values)
    }

hname2  = "default" # Default run, same as smpaper but with no SLAB replacement
hparam2 = {
    'ensorem' : 1,      # 1=enso removed, 0=not removed
    'ensolag' : 1,      # Lag Applied toENSO and Variable before removal
    'monwin'  : 3,      # Size of month window for HFF calculations
    'detrend' : 1,      # Whether or not variable was detrended
    'tails'   : 2,      # tails for t-test
    'p'       : 0.05,   # p-value for significance testing
    'sellags' : [0,],   # Lags included (indices, so 0=lag1)
    'lagstr'  : "lag1", # Name of lag based on sellags
    'method'  : 4       # Significance test option: 1 (No Mask); 2 (SST autocorr); 3 (SST-FLX crosscorr); 4 (Both), 5 (Replace with SLAB values)
    }

# Combine and make the dictionaries
hff_names = [hname1,hname2,]
hff_dicts = [hparam1,hparam2,]
hff_sets  = dict(zip(hff_names,hff_dicts))

#%%

