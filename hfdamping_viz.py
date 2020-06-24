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

from cartopy import config
import cartopy.feature as cfeature
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cmocean
import matplotlib.ticker as mticker
from cartopy.util import add_cyclic_point



def combine_dims(invar,targdim,combine):
    """ Function to bring target dimension to front and combine all other 
        dimensions. Can also do the reverse, if combine = 0.
        
        Inputs:
            1) invar:   [n-dim array] variable to reshape
            2) targdim: [scalar or nd array]
                            combine = 1, target dim to pull to front
                            combine = 0, dimensions to reshape to
            3) combine: [bool] 1 to combine, 0 to uncombine
        
        Dependencies:
            numpy as np
            
            
        Outputs:
            
            1) outvar: combined variable
            2) rearr:  array before reshaping
        
    """
    
    # Get tuple of variable shape
    varshape = invar.shape
    
    if combine == 1:
        
        # Get location of target dimension
        tdimloc = varshape.index(targdim)
        otherdims = [n for n in varshape if n != targdim]
        
        # Move target dim to front if it is not there already
        if tdimloc != 0:
            
            # total length (needed?)
            ndims   = len(varshape)
            
            # Get rearranged dimension order
            dimsbefore = np.arange(0,tdimloc,1)
            dimsafter  = np.arange(tdimloc+1,ndims,1)
            rearr      = np.hstack([tdimloc,dimsbefore,dimsafter])
            
            # Transpose arrays to bring dim to front
            invar = np.transpose(invar,rearr)
            
            # Save dimension sizes
            varshapenew = invar.shape
            
        else:
            varshapenew = varshape
                
        # Combine all other dimensions
        outvar = np.reshape(invar,(targdim,np.prod(otherdims)))
        
    # To recombine the variable...
    elif combine == 0:
        
        # Reshape the variable according to specifications,
        # Assuming output has been preserved
        outvar = np.reshape(invar,targdim)
        varshapenew = invar.shape
        
        
    return outvar,varshapenew
            
            


fluxes    = ("nhflx",)
monlist   = np.arange(1,13,1)
mnum      = np.arange(1,43,1)
lags      = np.arange(1,3,1)
monwin    = 3 # Smoothing window of months
lonremap  = 1 # Remap Longitude 360 -> -180
ensorem   = 1 # 1 if enso was removed
mode      = 3 # Sig Testing (0=notest,1=sst,2=flux,3=both)

plot_ensavg = 1
plotens  = 0



# --------------------
# Significance Testing
# --------------------
p     = 0.20
tails = 2
dof   = 82


# -------
#  Paths
# -------
# Path to hfdamping matfiles
datpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/01_Data/"
llpath  = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/01_Data/"
outpath = "/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/02_stochmod/02_Figures/20200629/"

# ----------------
# Plotting Options
# ----------------

# Geographic Range
bbox = [-100,20,-25,75]
cmap = cmocean.cm.balance
cint = np.arange(-50,60,10)

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
lat = np.squeeze(llmat['LAT'])
lon = np.squeeze(llmat['LON'])

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
    else:
        vmask = np.copy(damping)
    
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
    # Rearrange longitude and variables from degE to -180:180
    lonremap = 0
    if lonremap == 1:
        
        # Find western/eastern longitudes
        kxw = np.where(lon>180)
        kxe = np.where(lon<=180)
          
        # Concatenate eastern and western longitudes
        lon = np.hstack((lon[kxw]-360,lon[kxe]))
        
        # Permute variable to find the longitude dimension
        vshape = vmask.shape
        
        # Combine non longituded dimensions
        outvar,beforeshape=combine_dims(vmask,len(lon),1)
        
        # Make corresponding change to variable
        varnew = np.squeeze(np.concatenate((outvar[kxw,:],outvar[kxe,:]),1))
        
        # Reshape variable back to original form
        vmask,_ = combine_dims(varnew,beforeshape,0)
    
    # ---------------------------------
    # Compute Ensemble and lag averages
    # ---------------------------------
    
    # Permute month to 3rd axis [lonxlatxmonxensxlag]
    v_ensavg = np.transpose(vmask,(0,1,3,2,4))
    while len(v_ensavg.shape) > 3:
        v_ensavg = np.nanmean(v_ensavg,axis=3)

    
    
    # ------------------
    # Begin plotting....
    # ------------------
    plt.close('all')
    
    # Loop by month
    for m in range(1,len(monlist)+1):
        
        
        
        # -------------------
        # Ensemble Average Plot
        # -------------------
        if plot_ensavg == 1:
            
            # Select data for that month
            vplot = np.transpose(v_ensavg[:,:,m-1],(1,0))
            
            fig = plt.figure()
            
            cmap.set_bad(color='g')
            
            def makemapcontour(lon,lat,var,bbox,cmap,cint):
                
                # Add cyclic point
                var,lon = add_cyclic_point(var,coord=lon)
                
                ax = plt.axes(projection=ccrs.PlateCarree())
                ax.set_extent(bbox)
                
                ax.add_feature(cfeature.COASTLINE,facecolor='k')
                #ax.background_patch.set_fill(False)
                
                # Contour fill
                cs = ax.contourf(lon,lat,var,
                           cint,
                           cmap=cmap,
                           transform=ccrs.PlateCarree())
                
                
                # Negative contours
                cln = ax.contour(lon,lat,var,
                           cint[cint<0],
                           linestyles='dashed',
                           colors='k',
                           linewidths=0.5,
                           transform=ccrs.PlateCarree())
                
                # Positive Contours
                clp = ax.contour(lon,lat,var,
                           cint[cint>=0],
                           colors='k',
                           linewidths=0.5,
                           transform=ccrs.PlateCarree())    
                  
                
                                
                
                
                gl = ax.gridlines(draw_labels=True,linewidth=0.75,color='gray',linestyle='--')

                gl.xlabels_top = gl.ylabels_right = False
                gl.xformatter = LONGITUDE_FORMATTER
                gl.yformatter = LATITUDE_FORMATTER
                  
                return ax,cs,gl,cln,clp
            
            
            
            ax,cs,gl,cln,clp=makemapcontour(lon,lat,vplot,bbox,cmap,cint)
            gl.xlocator = mticker.FixedLocator(np.arange(-120,45,30))
            #gl.ylocator = mticker.FixedLocator(np.linspace(bbox[2],bbox[3],8))
            ax.background_patch.set_facecolor('yellow')
            fig.colorbar(cs)
            
            
            #gl.xlocator = mticker.FixedLocator([np.arange(-90,30,30)])
            
            titlestr = "Ensemble Average %s Feedback for %s " % (str.upper(flux),monlab[m-1])
            ax.set_title(titlestr,fontsize=14)
            
            
            # Save figure
            figname = "%s_Ensavg_m%02d_monwin%i.png" %(str.upper(flux),m,monwin)
            plt.savefig(outpath+figname, bbox_inches="tight",dpi=200)
            
            
            
            
            
        
        
        
        
        
        # -------------------
        # All ensemble plots
        # -------------------
        if plotens == 1:
            
            f1,axs = plt.subplots(6,7,figsize = (14,14))
            
            for e in range(len(mnum)):
                
                # Restrict variable to ensemble member and lag
                vplot = vmask[:,:,e,m-1,lags-1]
                
                # Take average until left with latxlon
                while len(vplot.shape) > 2:
                    vplot = np.nanmean(vplot,axis=2)
                
                vplot= np.transpose(vplot,(1,0))
                
                # Add plot to figure
                # http://unidata.github.io/awips2/python/awips-grids-and-cartopy/
                plt.subplot(6,7,e+1)
                
              
                
                
                
                axs[e+1],cs= makemap(lon,lat,vplot,bbox,cmap)
                
                #cb = plt.colorbar(cs,shrink=1,orientation='vertical')
                #cs.set_clim(caxischoose)
               # cbar = ax.colorbar()

                
                
                
                
                
            
            

            
            
    
        
        
    
