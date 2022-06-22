#!/usr/bin/env python3# -*- coding: utf-8 -*-"""Combines the output created by calc_enso_general for general large ensemblesCopied combine_damping_CESM1LE.pyCreated on Thu Jun 16 16:22:19 2022@author: gliu"""import globimport numpy as npimport xarray as xrimport matplotlib.pyplot as pltimport sys#%% Load in the datadataset_name = "csiro_mk36_lens" #"csiro_mk36_lens" #"gfdl_esm2m_lens"datpath      = "/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_lens/%s/hfdamping_qnet_1920to2005/" % dataset_namedebug        = False#%% Load in data (repetitive because last 5 lat values are slightly different)nclist       = glob.glob(datpath+"*.nc")nclist.sort()nens = len(nclist)print(*nclist)timestr      = nclist[0][-19:-9] #*Manually determined--need to make more flexicheck_allens = ["allens" in s for s in nclist]if np.any(check_allens):    print('Warning! Already merged files. Check %s. Exiting Script.'%(datpath))    sys.exit()lbd  = []cc   = []ac   = []cv   = []av   = []for i in range(nens):    ds = xr.open_dataset(nclist[i]) # [Month x Lag x Lat x Lon]        if debug:        ds.nhflx_damping.isel(month=0,lag=0).plot(),plt.show()        lbd.append(ds.nhflx_damping.values)    cc.append(ds.sst_flx_crosscorr.values)    ac.append(ds.sst_autocorr.values)    cv.append(ds.cov.values)    av.append(ds.autocov.values)    #%% Remake Dataarray# Prepare variables for inputinvars = [lbd,cc,ac,cv,av]invars = [np.array(v) for v in invars] # [ens x mon x lag x lat x lon]# Set some attributes (take from calc_enso_general)varnames = ("nhflx_damping",            "sst_flx_crosscorr",            "sst_autocorr",            "cov",            "autocov")varlnames = ("Net Heat Flux Damping",             "SST-Heat Flux Cross Correlation",             "SST Autocorrelation",             "SST-Heat Flux Covariance",             "SST Autocovariance")units     = ("W/m2/degC",             "Correlation",             "Correlation",             "W/m2*degC",             "degC2")dims  = {"ens":np.arange(1,nens+1,1),         "month":ds.month.values,         "lag":ds.lag.values,         "lat":ds.lat.values,         "lon":ds.lon.values           }das = []for v,name in enumerate(varnames):    attr_dict = {'long_name':varlnames[v],                 'units':units[v]}    da = xr.DataArray(invars[v],                dims=dims,                coords=dims,                name = name,                attrs=attr_dict                )    if v == 0:        ds = da.to_dataset() # Convert to dataset    else:        ds = ds.merge(da) # Merge other datasets            # Append to list if I want to save separate dataarrays    das.append(ds)#% Save as netCDF# ---------------savename = "%s%s_hfdamping_ensorem1_detrend1_%s_allens.nc" % (datpath,dataset_name,timestr)encoding_dict = {name : {'zlib': True} for name in varnames} print("Saving as " + savename)ds.to_netcdf(savename,         encoding=encoding_dict)