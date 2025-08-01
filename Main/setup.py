#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Set Up Script for Heat Flux Feedback Calculator
(Tried to write this, but stipp WIP)

Created on Fri Jun 20 14:59:50 2025

@author: gliu

"""

import os


def makedir(expdir):
    """
    Check if "expdir" exists, and creates a directory if it doesn't

    Parameters
    ----------
    expdir : TYPE
        DESCRIPTION.

    """
    checkdir = os.path.isdir(expdir)
    if not checkdir:
        print(expdir + " Not Found! \n\tCreating Directory...")
        os.makedirs(expdir)
    else:
        print(expdir+" was found!")
    return None



#%% Still trying to decide how to organize it
# Currently, user provides the dataset paths...
#     Calculation will run and output... 



# # Project Path
projpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/"



# rawpath  = ""
# Calculation name
#calcname = "OAFLUX_1984_2007_qnet"
#projpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/"




#%%




#%% Make the Project Directories (takem from calc_hff_general_new)



# Set Paths
anompath = projpath + "anom/"
ensopath = projpath + "enso/"
hffpath  = projpath + "hff/"
maskpath = projpath + "masks/"
procpath = projpath + "proc/"


# Check for directories
dirs     = [anompath,ensopath,hffpath,maskpath,procpath]
[proc.makedir(path) for path in dirs]


# hf.makedir(anompath)
# hf.makedir(ensopath)
# hf.makedir(hffpath)
# hf.makedir(maskpath)
# hf.makedir(procpath)



# # Main Directories
# rawpath  = projpath + "raw/"                    # Raw Files (SST, Fluxes)
# ensopath = projpath + "enso/"                   # ENSO Indices
# procpath = projpath + "enso_removed/"           # Variables, with ENSO removed
# hffpath  = projpath + "heat_flux_feedback/"     # Heat Flux Feedback Estimates

# dirs = [rawpath,ensopath,procpath,hffpath]
# [proc.makedir(path) for path in dirs]


