#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 23:56:07 2020

@author: gliu
"""
from scipy.io import loadmat
import numpy as np
import time
import matplotlib.pyplot as plt
from scipy import stats


fluxes    = ("nhflx",)

monwin    = 3 # Smoothing window of months
lonremap  = 1 # Remap Longitude 360 -> -180
ensorem   = 1 # 1 if enso was removed
mode      = 3 # Sig Testing (0=notest,1=sst,2=flux,3=both)

# --------------------
# Significance Testing
# --------------------
p     = 0.05
tails = 2
dof   = 82


# -------
#  Paths
# -------
# Path to hfdamping matfiles
datpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/"
llpath  = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/01_Data/"

# ----------------
# Plotting Options
# ----------------

# Geographic Range
bbox = [-100,20,-25,75]


# Colorbar
caxischoose = (-50,50)

# Plot Point Toggle
plotpt =1
if plotpt == 1:
    lonf = -72
    latf = 35

# Preset Labels
if monwin == 3:
    monlab = ('DJF','JFM','FMA',
        'MAM','AMJ','MJJ',
        'JJA','JAS','ASO',
        'SON','OND','NDJ')
elif monwin == 1:
    monlab = ('Jan','Feb','Mar',
    'Apr','May','Jun',
    'Jul','Aug','Sep',
    'Oct','Nov','Dec')



# ****************************************************************************
# Script Start

# -------------------------------------------------
# Perform advance calculations based on user inputs
# -------------------------------------------------
    
# Calculate Correlation Threshold for significance testing
critval   = stats.t.ppf(1-p/tails,dof)
corrthres = np.sqrt(1/((dof/critval**2)+1))



# -------------------
# Load in data
# -------------------

# Lat/Lon 
llmat = loadmat(llpath+'CESM1_LATLON.mat')
lat = llmat['LAT']
lon = llmat['LON']

# SST Autocorrelation Coefficients
smat = loadmat("%s%s_rho_ensorem%i_monwin%i.mat" %(datpath,'SST',ensorem,monwin))
rsst = smat['rsst'] # [ LON(288) x LAT(192) x ENSN(42) x MON(12) x LAG(3)]


# Loop by Damping Variable
for flux in fluxes:
    
    
    
    # Load in damping variable
    dmat = loadmat("%s%s_damping_ensorem%i_monwin%i.mat" %(datpath,flux,ensorem,monwin))
    
    # [ LON(288) x LAT(192) x ENSN(42) x MON(12) x LAG(3)]
    damping = dmat['damping']
    
    # Load in damping-sst cross correlation coefficients
    cmat = loadmat("%s%s_rho_ensorem%i_monwin%i.mat" %(datpath,flux,ensorem,monwin))
    rflx = cmat['rho']
    
    # -----------------------------------
    # Make Masks for Significance Testing
    # -----------------------------------
    
    # Use Threshold to create masks
    if mode != 0:
        cmask = rflx >= corrthres
        smask = rsst >= corrthres
    
    # Mode 1: Apply SST autocorrelation test only
    if mode == 1:
        
        vmask = np.multiply(damping,smask)
    
    # Mode 2: Apply SST-FLX cross correlation testing
    elif mode == 2:
        
        vmask = np.multiply(damping,cmask)
    
    elif mode == 3:
        
        vmask = np.multiply(damping,np.multiply(smask,cmask))
    
    # ----------------
    # Remap Variables
    # ----------------
    
    
    
