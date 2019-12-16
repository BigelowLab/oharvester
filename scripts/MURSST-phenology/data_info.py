# Copyright 2018, by the California Institute of Technology. ALL RIGHTS RESERVED.
# United States Government Sponsorship acknowledged. Any commercial use must be negotiated 
# with the Office of Technology Transfer at the California Institute of Technology.
#
# This software may be subject to U.S. export control laws. By accepting this software, 
# the user agrees to comply with all applicable U.S. export laws and regulations. 
# User has the responsibility to obtain export licenses, or other export authority 
# as may be required before exporting such information to foreign countries or providing access to foreign persons.

#!  /usr/bin/env python3
#

"""
Compute and plot the phenology index 
"""

from netCDF4 import Dataset
import numpy as np
import sys, os
import glob
import calendar
from datetime import datetime, date, timedelta
from time import gmtime, strftime
import time
import pickle
import math
import numpy.polynomial.polynomial as poly
from optparse import OptionParser
import datetime
import calendar
from pathlib import Path

################################################################################################   
def createdir(dirname):
  if not os.path.exists(dirname):
    os.makedirs(dirname)

def podaac_dataset_shorname(dataset_name):
    if dataset_name.find('MUR') != -1:
       return 'MUR-JPL-L4-GLOB-v4.1'
    elif dataset_name.find('CMC') != -1:
       return 'CMC0.2deg-CMC-L4-GLOB-v2.0' 
    elif dataset_name.find('NCEI') != -1:
       return 'AVHRR_OI-NCEI-L4-GLOB-v2.0'
    elif dataset_name.find('MODIS') != -1:
       return 'MODIS_Aqua_L3_CHLA_Daily_4km_V2014.0_R'
    else:
       return null

def latloninfo(rootdir, shortname):
    ncin = Dataset(rootdir+"/"+shortname+"_info.nc", 'r')
    lons = ncin.variables['lon'][:]
    lats = ncin.variables['lat'][:]
    ncin.close()
    return lats, lons

def latloninfo_seawifs(rootdir, startY):
    # read sample file to get dataset infro
    bfound = False
    for j in range(1, 367):
       for filename in glob.glob(rootdir+"/"+"{0:0>4}".format(startY)+"/"+"{0:0>3}".format(j)+"/S"+"{0:0>4}".format(startY)+"{0:0>3}".format(j)+'*.nc'):
         ncin = Dataset(filename, 'r')
         lons = ncin.variables['lon'][:]
         lats = ncin.variables['lat'][:]
         ncin.close()
         bfound = True
       if bfound == True:
         break
    return lats, lons

def getboxinfo(lat_boxsize, lon_boxsize, lat_min, lat_max, lon_min, lon_max, lats, lons):

  nlat = lats.size
  nlon = lons.size

  ### calculate how many boxes
  nx = int(math.ceil((lat_max-lat_min)/lat_boxsize))
  ny = int(math.ceil((lon_max-lon_min)/lon_boxsize))

  ### created box index
  box_lat_start = np.empty(nx)
  box_lat_end = np.empty(nx)
  box_lon_start = np.empty(ny)
  box_lon_end = np.empty(ny)

  for i in range(0, nx):
    lat_bottom = lat_min + i*lat_boxsize
    lat_top    = lat_min + (i+1)*lat_boxsize
    aindex = np.where( (lats >= lat_bottom) & (lats <= lat_top) )
    box_lat_start[i] = aindex[0][0]
    box_lat_end[i] = aindex[0][-1]
  for i in range(0, ny):
    lon_left = lon_min + i*lon_boxsize
    lon_right = lon_min + (i+1)*lon_boxsize
    aindex = np.where( (lons >= lon_left) & (lons <= lon_right) )
    box_lon_start[i] = aindex[0][0]
    box_lon_end[i] = aindex[0][-1]

  return nx, ny, box_lat_start, box_lat_end, box_lon_start, box_lon_end

