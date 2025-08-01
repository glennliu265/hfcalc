#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 20 15:15:53 2025

@author: gliu
"""




calcname = "OAFLUX_1984_2007"
rawpath  = "/Users/gliu/Globus_File_Transfer/Reanalysis/OAFLUX/"
flxname  = "qnet"
sstname  = "sst"
trop_pac = False # Set to True if Tropical Pacific crop is already provided
lensflag = False # True to look for ensemble members

tcrop    = True
ystart   = 1984
yend     = 2007

print("Initiating calculation for: %s" % (calcname))
print("\tRaw NetCDFs: %s" % rawpath)


