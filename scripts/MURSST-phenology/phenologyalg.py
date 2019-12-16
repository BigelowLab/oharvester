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
import configparser
import datetime
import calendar
from sklearn import linear_model, datasets
from pathlib import Path

import data_info

################################################################################################   
def moving_average(a, n=15):
    """Moving Average"""
    ret = np.cumsum(a.filled(0))
    ret[n:] = ret[n:] - ret[:-n]
    counts = np.cumsum(~a.mask)
    counts[n:] = counts[n:] - counts[:-n]
    ret[~a.mask] /= counts[~a.mask]
    ret[a.mask] = np.nan
    return ret

def runningMeanFast(x, N):
    """Running mean using numpy internal convolve routine"""
    #return np.convolve(x, np.ones((N,))/N, mode='full')[(N-1):]
    s = np.r_[x[N-1:0:-1],x,x[-1:-N:-1]]
    y=np.convolve(s, np.ones((N,))/N, mode='valid')
    return y[N/2:len(y)-N/2]

def parseoptions():
  usage = "Usage: %prog [options]"
  parser = OptionParser(usage)
  parser.add_option("-c", "--config", help="config filename", dest="config")

  # Parse command line arguments
  (options, args) = parser.parse_args()

  # print help if no arguments are given
  if(len(sys.argv) == 1):
    parser.print_help()
    exit(-1)

  if options.config == None:
    print('\nConfig filename is required !\nProgram will exit now !\n')
    parser.print_help()
    exit(-1)

  return( options )

def climatology(lats, lons, startY, endY, spring_thresh1, spring_thresh2, summer_thresh_offset, variablename, outdir, outfilename, dirname):

  ny = endY - startY + 1

  # define variables
  sst_all        = np.empty((lats.size, lons.size, 366))
  sst_all_number = np.empty((lats.size, lons.size, 366))
  data_year_all  = np.empty((ny, 12, 31, lats.size, lons.size))
  diff_sst_all        = np.empty((lats.size, lons.size, 366))

  data_year_all[:] = np.nan

  nlat = lats.size
  nlon = lons.size

  for i in range(startY, endY+1):
    print("Processing Year = " + str(i))
    nd = 0
    for j in range(1, 367):
       for filename in glob.glob(dirname+'/'+str(i)+'/'+"{0:0>3}".format(j)+'/*.nc'):
         ncfile = filename.rsplit( "/")[ -1 ]
         ncin = Dataset(filename, 'r')
         asst = ncin.variables[variablename][:]
         sst = asst.squeeze()
         #sst = asst[0,:,:]
         lons = ncin.variables['lon'][:]
         lats = ncin.variables['lat'][:]
         amask = ncin.variables['mask'][:]
         mask = amask[0,:,:]
         #get global attributes
         atts = ncin.__dict__
         for att,val in atts.items():
           if att.find('start_time') != -1: 
             start_time = val
             continue
         ncin.close()
         #aindex = ( (sst > -250.0) & (sst < 350.0) & (mask == 1) )
         #sst_all[aindex,j-1] = sst_all[aindex,j-1] + sst[aindex]
         #sst_all_number[aindex,j-1] = sst_all_number[aindex,j-1] + 1
         aindex = np.where( (sst > -250.0) & (sst < 350.0) & (mask == 1) )
         sst_all[aindex[0][:],aindex[1][:],j-1] = sst_all[aindex[0][:],aindex[1][:],j-1] + sst[aindex[0][:],aindex[1][:]]
         sst_all_number[aindex[0][:],aindex[1][:],j-1] = sst_all_number[aindex[0][:],aindex[1][:],j-1] + 1
         data_year_all[i-startY,int(start_time[4:6])-1,int(start_time[6:8])-1, :,:] = sst[0:nlat,0:nlon]

  # calculate climatology
  sst_all = sst_all/sst_all_number
  diff_sst_all[:,:,1:] = np.diff(sst_all)

  # trend calculation
  data_rate=np.empty((12, nlat,nlon))
  atemp=np.empty((ny, 28, 12))
  data=np.empty((12,ny))

  for k in range(12):
    for j in range(nlat):
      for i in range(nlon):
         atemp[:,:,k] = data_year_all[:,k,0:28,j,i]
         b = atemp < 0.0
         if np.count_nonzero(b) > 0:
           atemp[b]=np.nan
         if np.count_nonzero(b) == atemp.size:
           data_rate[k,j,i] = np.nan 
         else:
           data[k,:] = np.nanmean(atemp[:,:,k], axis=1)
           xx = np.arange(ny)
           yy = np.squeeze(data[k,:])
           #z = np.polyfit(range(ny), data[k,:], 1)
           z = np.polyfit(xx[np.invert(np.isnan(yy))], yy[np.invert(np.isnan(yy))], 1)
           data_rate[k,j,i] = z[0]

  ### create output directory
  data_info.createdir(outdir)

  #output into netCDF file
  fid = Dataset(outdir+'/'+outfilename,'w')

  # Define the dimensions
  nlat = fid.createDimension('lat', nlat) # Unlimited
  nlon = fid.createDimension('lon', nlon) # Unlimited
  time = fid.createDimension('time', 366) # Unlimited
  month = fid.createDimension('month', 12) # Unlimited
  nyear = fid.createDimension('year', ny) # Unlimited

  nc_var = fid.createVariable('lat', 'f8',('lat'),zlib=True)
  fid.variables['lat'][:] = lats
  fid.variables['lat'].standard_name='latitude'
  fid.variables['lat'].long_name='latitude'
  fid.variables['lat'].axis='Y'
  fid.variables['lat'].units='degrees_north'

  nc_var = fid.createVariable('lon', 'f8',('lon'),zlib=True)
  fid.variables['lon'][:] = lons
  fid.variables['lon'].standard_name='longitude'
  fid.variables['lon'].long_name='longitude'
  fid.variables['lon'].units='degrees_east'
  fid.variables['lon'].axis='X'

  nc_var = fid.createVariable('time', 'i4',('time'),zlib=True)
  fid.variables['time'][:] = range(1,367)
  fid.variables['time'].standard_name='time'
  fid.variables['time'].long_name='day number'
  fid.variables['time'].units='day number'
  fid.variables['time'].axis='T'

  nc_var = fid.createVariable('year', 'i4',('year'),zlib=True)
  fid.variables['year'][:] = range(startY,endY+1)
  fid.variables['year'].standard_name='Year'

  nc_var = fid.createVariable('sst_climatology', 'f8',('lat','lon','time'),zlib=True)
  fid.variables['sst_climatology'][:] = sst_all
  fid.variables['sst_climatology'].standard_name='SST Average'
  fid.variables['sst_climatology'].units='kelvin'

  nc_var = fid.createVariable('diff_sst_climatology', 'f8',('lat','lon','time'),zlib=True)
  fid.variables['diff_sst_climatology'][:] = diff_sst_all
  fid.variables['diff_sst_climatology'].standard_name='SST Different with previous day'
  fid.variables['diff_sst_climatology'].units='kelvin'

  nc_var = fid.createVariable('months', 'f8',('month'))

  nc_var = fid.createVariable('data_rate', 'f8',('month', 'lat', 'lon'))
  fid.variables['data_rate'][:] = data_rate
  fid.variables['data_rate'].units='Kelvin/Yr'
  fid.variables['data_rate'].comment='analysed_sst'

  fid.cdm_data_type         = "Grid"
  fid.close()
  
  return

def climatology_chlor(lats, lons, startY, endY, latmin, latmax, lonmin, lonmax, variablename, outdir, outfilename, dirname):
  """Climatology Calculation for Chlorophyll concentration

  The code calculates the average over all the years of Chlorophyll 
  concentration for each day, and the output file is saved in the netCDF
  format.
  """

  ny = endY - startY + 1

  nlat = lats.size
  nlon = lons.size

  # define variables
  sst_all        = np.empty((nlat, nlon, 366))
  sst_all_number = np.empty((nlat, nlon, 366))
  data_year_all  = np.empty((ny, 12, 31, nlat, nlon))

  data_year_all[:] = np.nan

  for i in range(startY, endY+1):
    print("Processing Year = " + str(i))
    nd = 0
    range_end = 367
    if(calendar.isleap(i)):
      range_end = 367
    else:
      range_end = 366
    for j in range(1, range_end):
       for filename in glob.glob(dirname+"/"+"{0:0>4}".format(i)+"/"+"{0:0>3}".format(j)+"/S"+"{0:0>4}".format(i)+"{0:0>3}".format(j)+'*.nc'):
         ncfile = filename.rsplit( "/")[ -1 ]
         ncin = Dataset(filename, 'r')
         asst = ncin.variables[variablename][:]
         sst = asst[::-1,:]
         #get global attributes
         atts = ncin.__dict__
         for att,val in atts.items():
           if att.find('start_time') != -1: 
             start_time = val
             continue
         ncin.close()
         aindex = np.where( (sst > 0.0) )
         sst_all[aindex[0][:],aindex[1][:],j-1] = sst_all[aindex[0][:],aindex[1][:],j-1] + sst[aindex[0][:],aindex[1][:]]
         sst_all_number[aindex[0][:],aindex[1][:],j-1] = sst_all_number[aindex[0][:],aindex[1][:],j-1] + 1
         #data_year_all[i-startY,int(start_time[4:6])-1,int(start_time[6:8])-1, :,:] = sst[0:nlat,0:nlon]

  # calculate climatology
  sst_all = sst_all/sst_all_number

  # trend calculation
  data_rate=np.empty((12, nlat,nlon))
  atemp=np.empty((ny, 28, 12))
  data=np.empty((12,ny))

  """
  for k in range(12):
    for j in range(nlat):
      for i in range(nlon):
         atemp[:,:,k] = data_year_all[:,k,0:28,j,i]
         b = atemp < 0.0
         if np.count_nonzero(b) > 0:
           atemp[b]=np.nan
         data[k,:] = np.nanmean(atemp[:,:,k], axis=1)
         z = np.polyfit(range(ny), data[k,:], 1)
         data_rate[k,j,i] = z[0]
  """

  ### create output directory
  data_info.createdir(outdir)

  #output into netCDF file
  fid = Dataset(outdir+'/'+outfilename,'w')

  # Define the dimensions
  nlat = fid.createDimension('lat', nlat) # Unlimited
  nlon = fid.createDimension('lon', nlon) # Unlimited
  time = fid.createDimension('time', 366) # Unlimited
  month = fid.createDimension('month', 12) # Unlimited
  nyear = fid.createDimension('year', ny) # Unlimited

  nc_var = fid.createVariable('lat', 'f8',('lat'),zlib=True)
  fid.variables['lat'][:] = lats
  fid.variables['lat'].standard_name='latitude'
  fid.variables['lat'].long_name='latitude'
  fid.variables['lat'].axis='Y'
  fid.variables['lat'].units='degrees_north'

  nc_var = fid.createVariable('lon', 'f8',('lon'),zlib=True)
  fid.variables['lon'][:] = lons
  fid.variables['lon'].standard_name='longitude'
  fid.variables['lon'].long_name='longitude'
  fid.variables['lon'].units='degrees_east'
  fid.variables['lon'].axis='X'

  nc_var = fid.createVariable('time', 'i4',('time'),zlib=True)
  fid.variables['time'][:] = range(1,367)
  fid.variables['time'].standard_name='time'
  fid.variables['time'].long_name='day number'
  fid.variables['time'].units='day number'
  fid.variables['time'].axis='T'

  nc_var = fid.createVariable('year', 'i4',('year'),zlib=True)
  fid.variables['year'][:] = range(startY,endY+1)
  fid.variables['year'].standard_name='Year'

  nc_var = fid.createVariable('chlor_climatology', 'f8',('lat','lon','time'),zlib=True)
  fid.variables['chlor_climatology'][:] = sst_all
  fid.variables['chlor_climatology'].standard_name='Chlor Average'
  fid.variables['chlor_climatology'].units='mg m^-3'

  nc_var = fid.createVariable('chlor_num_days', 'i4',('lat','lon','time'),zlib=True)
  fid.variables['chlor_num_days'][:] = sst_all_number

  nc_var = fid.createVariable('months', 'f8',('month'))

  fid.cdm_data_type         = "Grid"
  fid.close()
  
  return

def metric_spring_start(lats, lons, startY, endY, spring_thresh1, spring_thresh2, variablename, outdir, outfilename, clim_filename, dirname):

  ny = endY - startY + 1

  # define max and min matrics
  day_min = np.empty((lats.size, lons.size, ny))
  day_max = np.empty((lats.size, lons.size, ny))
  data_min = np.empty((lats.size, lons.size, ny))
  data_max = np.empty((lats.size, lons.size, ny))

  # define spring start matrics
  day_spring_start1 = np.empty((lats.size, lons.size, ny))
  day_spring_start2 = np.empty((lats.size, lons.size, ny))
  data_year      = np.empty((lats.size, lons.size, 366, ny))

  # define spring start matrics trend
  day_spring_start_trend1 = np.empty((lats.size, lons.size))
  day_spring_start_trend2 = np.empty((lats.size, lons.size))

  # initilization with NAN
  data_year[:] = np.nan
  day_spring_start1[:] = np.nan
  day_spring_start2[:] = np.nan
  day_spring_start_trend1[:] = np.nan
  day_spring_start_trend2[:] = np.nan

  nlat = lats.size
  nlon = lons.size

  # Read in climatology
  ncin = Dataset(clim_filename, 'r')
  diff_clim = ncin.variables['diff_sst_climatology'][:]
  ncin.close()

  day_turn = np.empty([lats.size, lons.size])
  day_turn[:,:] = np.nan

  for i in range(0, lats.size):
    for j in range(0, lons.size):
      xrange = np.arange(366)
      yrange = diff_clim[i,j,:]
      xtemp = xrange[np.invert(np.isnan(yrange))]
      if xtemp.size > 2:
        z = np.polyfit(xrange[np.invert(np.isnan(yrange))], yrange[np.invert(np.isnan(yrange))], 6)
        polynomial = np.poly1d(z)
        ys = polynomial(np.arange(366))

        for k in range(0, 365):
          if(ys[k] > 0):
            day_turn[i, j] = k + 1
            break
      else:
        day_turn[i, j] = np.nan

  """
  fig= plt.figure(figsize=(6,6), dpi=300)
  clevs = np.linspace(1, 100, 101)
  cmap=plt.get_cmap("jet")
  plt.contourf(day_turn, clevs, cmap=cmap)
  plt.show()
  print(lats[22])
  print(lons[40])
  """

  for i in range(startY, endY+1):
    print("Processing Year = " + str(i))
    nd = 0
    for j in range(1, 367):
       for filename in glob.glob(dirname+'/'+str(i)+'/'+"{0:0>3}".format(j)+'/*.nc'):
         if os.path.isfile(filename):
           ncfile = filename.rsplit( "/")[ -1 ]
           ncin = Dataset(filename, 'r')
           asst = ncin.variables[variablename][:]
           sst = asst.squeeze()
           #sst = asst[0,:,:]
           lons = ncin.variables['lon'][:]
           lats = ncin.variables['lat'][:]
           amask = ncin.variables['mask'][:]
           mask = amask[0,:,:]
           #get global attributes
           atts = ncin.__dict__
           for att,val in atts.items():
             if att.find('start_time') != -1: 
               start_time = val
               continue
           ncin.close()
           aindex = np.where( (sst > -250.0) & (sst < 350.0) & (mask == 1) )
           data_year[aindex[0][:],aindex[1][:],j-1,i-startY] = sst[aindex[0][:],aindex[1][:]]

    # data max/date and min/date calculation
    for i1 in range(0, lats.size): 
      for j1 in range(0, lons.size):
           atemp = data_year[i1,j1,:,i-startY]
           atemp1 = atemp.squeeze()
           mx = np.ma.masked_array(atemp1,np.isnan(atemp1))
           atemp = moving_average(mx)
           if np.count_nonzero(np.isnan(atemp)) == atemp.size:
             day_min[i1,j1,i-startY]  = np.nan 
             data_min[i1,j1,i-startY] = np.nan
             day_max[i1,j1,i-startY]  = np.nan
             data_max[i1,j1,i-startY] = np.nan
           else: 
             day_min[i1,j1,i-startY]  = np.nanargmin(atemp)
             if day_min[i1,j1,i-startY] > 300:
               day_min[i1,j1,i-startY] = np.nan
             data_min[i1,j1,i-startY] = np.nanmin(atemp)
             day_max[i1,j1,i-startY]  = np.nanargmax(atemp)
             data_max[i1,j1,i-startY] = np.nanmax(atemp)

    # spring start matrics calculation
    for i1 in range(0, lats.size): 
      for j1 in range(0, lons.size):
           atemp = data_year[i1,j1,:,i-startY]
           atemp = atemp.squeeze()
           mx = np.ma.masked_array(atemp,np.isnan(atemp))
           atemp = moving_average(mx)
           if np.count_nonzero(np.isnan(atemp)) <= atemp.size and ~np.isnan(day_min[i1,j1,i-startY]):
             dmin = int(day_min[i1,j1,i-startY])
             datamin = data_min[i1,j1,i-startY]

             bfind1 = False
             for k in range(int(day_turn[i1, j1]), 300):
               if atemp[k] > spring_thresh1:
                 bfind1 = True 
                 """
                 for i2 in range(1,11):
                   if atemp[k+i2] < spring_thresh1:
                     bfind1 = False
                 """
                 #if bfind1 == True:
                 day_spring_start1[i1,j1,i-startY] = k
                 break
             if bfind1 == False:
               day_spring_start1[i1,j1,i-startY] = np.nan 

             bfind2 = False
             for k in range(int(day_turn[i1, j1]), 300):
               if atemp[k] > spring_thresh2 and k > 0:
                 bfind2 = True 
                 """
                 for i2 in range(1,11):
                   if atemp[k+i2] < spring_thresh2:
                     bfind2 = False
                 """
                 #if bfind2 == True:
                 day_spring_start2[i1,j1,i-startY] = k
                 break
             if bfind2 == False:
               day_spring_start2[i1,j1,i-startY] = np.nan 

  # Calculate the spring start day trend 
  x = np.arange(ny)
  for i1 in range(0, lats.size): 
    for j1 in range(0, lons.size):
      y = day_spring_start1[i1,j1,:].squeeze()
      idx = np.isfinite(y)
      x1 = x[idx]
      if x1.size > ny/2:
        z = np.polyfit(x[idx], y[idx], 1)
        day_spring_start_trend1[i1,j1] = z[0]
        """
        if np.abs(z[0]) < 1.0:
          day_spring_start_trend1[i1,j1] = z[0]
        else:
          print(str(i1) + ", "+ str(j1))
          y1 = y[idx]
          X = np.empty(shape=(x1.size,1))
          X[:,0] = x1
          if (i1 == 22) and (j1 == 40):
           for kk in range(0,36):
             plt.plot(data_year[i1,j1,:,kk])
             plt.plot([0, 365], [281, 281], color='r', linewidth=3)
             plt.xlim([0, 365])
             plt.xlabel("Day #")
             plt.ylabel("SST (K)")
           plt.show()
          ransac = linear_model.RANSACRegressor()
          ransac.fit(X, y1)
          day_spring_start_trend1[i1,j1] = ransac.estimator_.coef_[0]
        """
      else:
        day_spring_start_trend1[i1,j1] = np.nan

      y = day_spring_start2[i1,j1,:].squeeze()
      idx = np.isfinite(y)
      x1 = x[idx]
      if x1.size > ny/2:
        z = np.polyfit(x[idx], y[idx], 1)
        day_spring_start_trend2[i1,j1] = z[0]
        """  
        if np.abs(z[0]) < 1.0:
          day_spring_start_trend2[i1,j1] = z[0]
        else:
          y1 = y[idx]
          X = np.empty(shape=(x1.size,1))
          X[:,0] = x1
          print(x1)
          print(y1)
          ransac = linear_model.RANSACRegressor()
          ransac.fit(X, y1)
          day_spring_start_trend2[i1,j1] = ransac.estimator_.coef_[0]
        """  
      else:
        day_spring_start_trend2[i1,j1] = np.nan


  # Update the variables with _Fill_Value (-9999)
  day_min[np.isnan(day_min)] = -9999
  day_max[np.isnan(day_max)] = -9999
  data_min[np.isnan(data_min)] = -9999
  data_max[np.isnan(data_max)] = -9999
  day_spring_start1[np.isnan(day_spring_start1)] = -9999
  day_spring_start_trend1[np.isnan(day_spring_start_trend1)] = -9999
  day_spring_start2[np.isnan(day_spring_start2)] = -9999
  day_spring_start_trend2[np.isnan(day_spring_start_trend2)] = -9999

  ### create output directory
  data_info.createdir(outdir)

  #output into netCDF file
  fid = Dataset(outdir+'/'+outfilename,'w')

  # Define the dimensions
  nlat = fid.createDimension('lat', nlat) # Unlimited
  nlon = fid.createDimension('lon', nlon) # Unlimited
  nyear = fid.createDimension('year', ny) # Unlimited

  nc_var = fid.createVariable('lat', 'f8',('lat'),zlib=True)
  fid.variables['lat'][:] = lats
  fid.variables['lat'].standard_name='latitude'
  fid.variables['lat'].long_name='latitude'
  fid.variables['lat'].axis='Y'
  fid.variables['lat'].units='degrees_north'

  nc_var = fid.createVariable('lon', 'f8',('lon'),zlib=True)
  fid.variables['lon'][:] = lons
  fid.variables['lon'].standard_name='longitude'
  fid.variables['lon'].long_name='longitude'
  fid.variables['lon'].units='degrees_east'
  fid.variables['lon'].axis='X'

  nc_var = fid.createVariable('year', 'i4',('year'),zlib=True)
  fid.variables['year'][:] = range(startY,endY+1)
  fid.variables['year'].standard_name='Year'

  nc_var = fid.createVariable('day_min', 'i4',('lat','lon','year'),fill_value=-9999,zlib=True)
  fid.variables['day_min'][:] = day_min
  fid.variables['day_min'].standard_name='coldest day of the year'
  fid.variables['day_min'].units='Day number in a year'

  nc_var = fid.createVariable('day_max', 'i4',('lat','lon','year'),fill_value=-9999,zlib=True)
  fid.variables['day_max'][:] = day_max
  fid.variables['day_max'].standard_name='warmest day of the year'
  fid.variables['day_max'].units='Day number in a year'

  nc_var = fid.createVariable('data_min', 'f8',('lat','lon','year'),fill_value=-9999,zlib=True)
  fid.variables['data_min'][:] = data_min
  fid.variables['data_min'].standard_name='coldest data of the year'
  fid.variables['data_min'].units=''

  nc_var = fid.createVariable('data_max', 'f8',('lat','lon','year'),fill_value=-9999,zlib=True)
  fid.variables['data_max'][:] = data_max
  fid.variables['data_max'].standard_name='warmest data of the year'
  fid.variables['data_max'].units=''

  nc_var = fid.createVariable('day_spring1', 'i4',('lat','lon','year'),fill_value=-9999,zlib=True)
  fid.variables['day_spring1'][:] = day_spring_start1
  fid.variables['day_spring1'].standard_name='spring start day of the year'
  fid.variables['day_spring1'].units='Day number in a year'

  nc_var = fid.createVariable('day_spring2', 'i4',('lat','lon','year'),fill_value=-9999,zlib=True)
  fid.variables['day_spring2'][:] = day_spring_start2
  fid.variables['day_spring2'].standard_name='spring start day of the year'
  fid.variables['day_spring2'].units='Day number in a year'

  nc_var = fid.createVariable('day_spring_trend1', 'f8',('lat','lon'),fill_value=-9999,zlib=True)
  fid.variables['day_spring_trend1'][:] = day_spring_start_trend1
  fid.variables['day_spring_trend1'].standard_name='spring start day trend per year'
  fid.variables['day_spring_trend1'].units='Day number per year'

  nc_var = fid.createVariable('day_spring_trend2', 'f8',('lat','lon'),fill_value=-9999,zlib=True)
  fid.variables['day_spring_trend2'][:] = day_spring_start_trend2
  fid.variables['day_spring_trend2'].standard_name='spring start day trend per year'
  fid.variables['day_spring_trend2'].units='Day number per year'

  fid.cdm_data_type         = "Grid"
  fid.close()
  
  return

def metric_summer_start_end(lats, lons, startY, endY, summer_thresh_offset, variablename, outdir, outfilename, clim_filename, dirname):

  ny = endY - startY + 1

  # define max and min matrics
  day_min = np.empty((lats.size, lons.size, ny))
  day_max = np.empty((lats.size, lons.size, ny))
  data_min = np.empty((lats.size, lons.size, ny))
  data_max = np.empty((lats.size, lons.size, ny))

  # define summer start/end matrics
  day_summer_start = np.empty((lats.size, lons.size, ny))
  day_summer_end = np.empty((lats.size, lons.size, ny))
  data_summer_max_all = np.empty((lats.size, lons.size, ny))
  data_summer_max = np.empty((lats.size, lons.size))
  day_summer_max = np.empty((lats.size, lons.size))

  sst_all        = np.empty((lats.size, lons.size, 366))
  sst_all_number = np.empty((lats.size, lons.size, 366))
  data_year_all  = np.empty((ny, 12, 31, lats.size, lons.size))
  data_year      = np.empty((lats.size, lons.size, 366, ny))

  # define summer start/end matrics trend
  day_summer_start_trend = np.empty((lats.size, lons.size))
  day_summer_end_trend = np.empty((lats.size, lons.size))

  # initilization with NAN
  data_year_all[:] = np.nan
  data_year[:] = np.nan
  day_summer_start[:] = np.nan
  day_summer_end[:] = np.nan
  data_summer_max[:] = np.nan
  data_summer_max_all[:] = np.nan
  day_summer_max[:] = np.nan
  day_summer_start_trend[:] = np.nan
  day_summer_end_trend[:] = np.nan

  nlat = lats.size
  nlon = lons.size

  # Read in climatology
  ncin = Dataset(clim_filename, 'r')
  diff_clim = ncin.variables['diff_sst_climatology'][:]
  ncin.close()

  day_turn = np.empty([lats.size, lons.size])
  day_turn[:,:] = np.nan

  for i in range(0, lats.size):
    for j in range(0, lons.size):
      xrange = np.arange(366)
      yrange = diff_clim[i,j,:].squeeze()

      idx = np.argwhere(np.isfinite(yrange))
      if idx.size > len(xrange)/2:
        z = np.polyfit(xrange[idx[:,0].squeeze()], yrange[idx[:,0].squeeze()], 6)
        polynomial = np.poly1d(z)
        ys = polynomial(range(366))

        for k in range(0, 365):
          if(ys[k] > 0):
            day_turn[i, j] = k + 1
            break

  for i in range(startY, endY+1):
    print("Processing Year = " + str(i))
    nd = 0
    for j in range(1, 367):
       for filename in glob.glob(dirname+'/'+str(i)+'/'+"{0:0>3}".format(j)+'/*.nc'):
         ncfile = filename.rsplit( "/")[ -1 ]
         ncin = Dataset(filename, 'r')
         asst = ncin.variables[variablename][:]
         sst = asst.squeeze()
         #sst = asst[0,:,:]
         lons = ncin.variables['lon'][:]
         lats = ncin.variables['lat'][:]
         amask = ncin.variables['mask'][:]
         mask = amask[0,:,:]
         #get global attributes
         atts = ncin.__dict__
         for att,val in atts.items():
           if att.find('start_time') != -1: 
             start_time = val
             continue
         ncin.close()
         aindex = np.where( (sst > -250.0) & (sst < 350.0) & (mask == 1) )
         data_year[aindex[0][:],aindex[1][:],j-1,i-startY] = sst[aindex[0][:],aindex[1][:]]

    # data max/date and min/date calculation
    for i1 in range(0, lats.size): 
      for j1 in range(0, lons.size):
           atemp = data_year[i1,j1,:,i-startY]
           atemp1 = atemp.squeeze()
           mx = np.ma.masked_array(atemp1,np.isnan(atemp1))
           atemp = moving_average(mx)
           if np.count_nonzero(np.isnan(atemp)) == atemp.size:
             day_min[i1,j1,i-startY]  = np.nan 
             data_min[i1,j1,i-startY] = np.nan
             day_max[i1,j1,i-startY]  = np.nan
             data_max[i1,j1,i-startY] = np.nan
           else: 
             day_min[i1,j1,i-startY]  = np.nanargmin(atemp)
             if day_min[i1,j1,i-startY] > 300:
               day_min[i1,j1,i-startY] = np.nan
             data_min[i1,j1,i-startY] = np.nanmin(atemp)
             day_max[i1,j1,i-startY]  = np.nanargmax(atemp)
             data_max[i1,j1,i-startY] = np.nanmax(atemp)

    # summer start/end matrics calculation
    # Step 1: find summer maximum for each year 
    for i1 in range(0, lats.size): 
      for j1 in range(0, lons.size):
           atemp = data_year[i1,j1,150:150+90,i-startY]
           atemp = atemp.squeeze()
           mx = np.ma.masked_array(atemp,np.isnan(atemp))
           atemp = moving_average(mx)
           if np.count_nonzero(np.isnan(atemp)) == atemp.size:
             data_summer_max_all[i1,j1,i-startY]  = np.nan 
           else: 
             data_summer_max_all[i1,j1,i-startY] = np.nanmax(atemp)

  # summer start/end matrics calculation
  # Step 2: summer lowest maximum calculation
  for i1 in range(0, lats.size): 
    for j1 in range(0, lons.size):
       atemp = data_summer_max_all[i1,j1,:]
       atemp = atemp.squeeze()
       if np.count_nonzero(np.isnan(atemp)) == atemp.size:
         data_summer_max[i1,j1]  = np.nan 
         day_summer_max[i1,j1]  = np.nan 
       else: 
         data_summer_max[i1,j1] = np.nanmin(atemp)
         day_summer_max[i1,j1] = np.nanargmin(atemp)+startY

  # summer start/end matrics calculation
  # Step 3: summer start/end calculation
  for i in range(startY, endY+1):
    for i1 in range(0, lats.size): 
      for j1 in range(0, lons.size):
           atemp = data_year[i1,j1,:,i-startY]
           atemp = atemp.squeeze()
           mx = np.ma.masked_array(atemp,np.isnan(atemp))
           atemp = moving_average(mx)
           """
           if np.count_nonzero(np.isnan(atemp)) == atemp.size:
             day_summer_start[i1,j1,i-startY]  = np.nan 
           elif np.isnan(day_turn[i1, j1]):
             day_summer_start[i1,j1,i-startY]  = np.nan 
           else: 
           """
           if np.count_nonzero(np.isnan(atemp)) <= atemp.size and ~np.isnan(day_min[i1,j1,i-startY]):
             for k in range(int(day_turn[i1, j1]), atemp.size):
               if atemp[k] >= data_summer_max[i1,j1]+summer_thresh_offset:
                 day_summer_start[i1,j1,i-startY] = k
                 break
             for k in range(np.nanargmax(atemp), atemp.size):
               if atemp[k] <= data_summer_max[i1,j1]+summer_thresh_offset:
                 day_summer_end[i1,j1,i-startY] = k
                 break

  # Calculate the summer start/end day trend 
  x = np.arange(ny)
  for i1 in range(0, lats.size): 
    for j1 in range(0, lons.size):
      """
      z = np.polyfit(range(ny), day_summer_start[i1,j1,:], 1)
      day_summer_start_trend[i1,j1] = z[0]
      z = np.polyfit(range(ny), day_summer_end[i1,j1,:], 1)
      day_summer_end_trend[i1,j1] = z[0]
      """
      y = day_summer_start[i1,j1,:].squeeze()
      idx = np.isfinite(y)
      x1 = x[idx]
      if x1.size > ny/2:
        z = np.polyfit(x[idx], y[idx], 1)
        day_summer_start_trend[i1,j1] = z[0]
      else:
        day_summer_start_trend[i1,j1] = np.nan

      y = day_summer_end[i1,j1,:].squeeze()
      idx = np.isfinite(y)
      x1 = x[idx]
      if x1.size > ny/2:
        z = np.polyfit(x[idx], y[idx], 1)
        day_summer_end_trend[i1,j1] = z[0]
      else:
        day_summer_end_trend[i1,j1] = np.nan

  # Update the variables with _Fill_Value (-9999)
  day_min[np.isnan(day_min)] = -9999
  day_max[np.isnan(day_max)] = -9999
  data_min[np.isnan(data_min)] = -9999
  data_max[np.isnan(data_max)] = -9999
  day_summer_start[np.isnan(day_summer_start)] = -9999
  day_summer_end[np.isnan(day_summer_end)] = -9999
  day_summer_start_trend[np.isnan(day_summer_start_trend)] = -9999
  day_summer_end_trend[np.isnan(day_summer_end_trend)] = -9999

  ### create output directory
  data_info.createdir(outdir)

  #output into netCDF file
  fid = Dataset(outdir+'/'+outfilename,'w')
  # Define the dimensions
  nlat = fid.createDimension('lat', nlat) # Unlimited
  nlon = fid.createDimension('lon', nlon) # Unlimited
  nyear = fid.createDimension('year', ny) # Unlimited

  nc_var = fid.createVariable('lat', 'f8',('lat'),zlib=True)
  fid.variables['lat'][:] = lats
  fid.variables['lat'].standard_name='latitude'
  fid.variables['lat'].long_name='latitude'
  fid.variables['lat'].axis='Y'
  fid.variables['lat'].units='degrees_north'

  nc_var = fid.createVariable('lon', 'f8',('lon'),zlib=True)
  fid.variables['lon'][:] = lons
  fid.variables['lon'].standard_name='longitude'
  fid.variables['lon'].long_name='longitude'
  fid.variables['lon'].units='degrees_east'
  fid.variables['lon'].axis='X'

  nc_var = fid.createVariable('year', 'i4',('year'),zlib=True)
  fid.variables['year'][:] = range(startY,endY+1)
  fid.variables['year'].standard_name='Year'

  nc_var = fid.createVariable('day_min', 'i4',('lat','lon','year'),fill_value=-9999,zlib=True)
  fid.variables['day_min'][:] = day_min
  fid.variables['day_min'].standard_name='coldest day of the year'
  fid.variables['day_min'].units='Day number in a year'

  nc_var = fid.createVariable('day_max', 'i4',('lat','lon','year'),fill_value=-9999,zlib=True)
  fid.variables['day_max'][:] = day_max
  fid.variables['day_max'].standard_name='warmest day of the year'
  fid.variables['day_max'].units='Day number in a year'

  nc_var = fid.createVariable('data_min', 'f8',('lat','lon','year'),fill_value=-9999,zlib=True)
  fid.variables['data_min'][:] = data_min
  fid.variables['data_min'].standard_name='coldest data of the year'
  fid.variables['data_min'].units=''

  nc_var = fid.createVariable('data_max', 'f8',('lat','lon','year'),fill_value=-9999,zlib=True)
  fid.variables['data_max'][:] = data_max
  fid.variables['data_max'].standard_name='warmest data of the year'
  fid.variables['data_max'].units=''

  nc_var = fid.createVariable('day_summer_start', 'i4',('lat','lon','year'),fill_value=-9999,zlib=True)
  fid.variables['day_summer_start'][:] = day_summer_start
  fid.variables['day_summer_start'].standard_name='summer start day of the year'
  fid.variables['day_summer_start'].units='Day number in a year'

  nc_var = fid.createVariable('day_summer_end', 'i4',('lat','lon','year'),fill_value=-9999,zlib=True)
  fid.variables['day_summer_end'][:] = day_summer_end
  fid.variables['day_summer_end'].standard_name='summer start day of the year'
  fid.variables['day_summer_end'].units='Day number in a year'

  nc_var = fid.createVariable('day_summer_start_trend', 'f8',('lat','lon'),fill_value=-9999,zlib=True)
  fid.variables['day_summer_start_trend'][:] = day_summer_start_trend
  fid.variables['day_summer_start_trend'].standard_name='Summer End Day Trend per Year'
  fid.variables['day_summer_start_trend'].units='Day number per year'

  nc_var = fid.createVariable('day_summer_end_trend', 'f8',('lat','lon'),fill_value=-9999,zlib=True)
  fid.variables['day_summer_end_trend'][:] = day_summer_end_trend
  fid.variables['day_summer_end_trend'].standard_name='Summer Start Day Trend per Year'
  fid.variables['day_summer_end_trend'].units='Day number per Year'

  fid.cdm_data_type         = "Grid"
  fid.close()
  
  return

def climatologybox(lat_boxsize, lon_boxsize, lat_min, lat_max, lon_min, lon_max, lats, lons, variablename, clim_filename, outfilename):

  nstep = 1
  nd = 366

  ### get box info
  nx, ny, box_lat_start, box_lat_end, box_lon_start, box_lon_end = data_info.getboxinfo(lat_boxsize, lon_boxsize, lat_min, lat_max, lon_min, lon_max, lats, lons)

  # Read in climatology
  print(clim_filename)
  ncin = Dataset(clim_filename, 'r')
  sst_clim = ncin.variables['sst_climatology'][:]
  ncin.close()

  ### create box average of parameters
  sst_avg=np.empty((nd, nx, ny))

  for k in range(0, nd):
       for ii in range(0, nx):
          for jj in range(0, ny):
             atemp = sst_clim[int(box_lat_start[ii]):int(box_lat_end[ii]):nstep, int(box_lon_start[jj]):int(box_lon_end[jj]):nstep, k-1]
             atemp = atemp.flatten()
             b = (atemp < -100.0) | (atemp > 350.0)
             if np.count_nonzero(b) > 0:
               atemp[b]=np.nan
             sst_avg[k,ii,jj] = np.nanmean(atemp)

  fid = Dataset(outfilename,'w')
  # Define the dimensions
  nx = fid.createDimension('nx', nx) # Unlimited
  ny = fid.createDimension('ny', ny) # Unlimited
  day = fid.createDimension('daynumber', 366) # Unlimited

  nc_var = fid.createVariable('daynumber', 'i4',('daynumber'),zlib=True)
  fid.variables['daynumber'][:] = range(1, 367)
  fid.variables['daynumber'].standard_name='day number'
  fid.variables['daynumber'].units='day number'
  fid.variables['daynumber'].axis='T'
  fid.variables['daynumber'].comment="Starting on Jan 1 of earch year with day number 1 until the end of the year with either day number 365 or 366 depending on the leap year status."

  nc_var = fid.createVariable('sst_box_average', 'f8',('daynumber', 'nx', 'ny'),zlib=True)
  fid.variables['sst_box_average'][:] = sst_avg[:,:,:]
  fid.variables['sst_box_average'].standard_name='Data Average in the Box'
  fid.variables['sst_box_average'].units='kelvin'

  fid.close()

def climatologybox_chlor(lat_boxsize, lon_boxsize, lat_min, lat_max, lon_min, lon_max, lats, lons, variablename, clim_filename, outfilename):

  nstep = 1
  nd = 366

  ### get box info
  nx, ny, box_lat_start, box_lat_end, box_lon_start, box_lon_end = data_info.getboxinfo(lat_boxsize, lon_boxsize, lat_min, lat_max, lon_min, lon_max, lats, lons)

  # Read in climatology
  print(clim_filename)
  ncin = Dataset(clim_filename, 'r')
  sst_clim = ncin.variables['chlor_climatology'][:]
  ncin.close()

  ### create box average of parameters
  sst_avg=np.empty((nd, nx, ny))

  for k in range(0, nd):
       for ii in range(0, nx):
          for jj in range(0, ny):
             atemp = sst_clim[int(box_lat_start[ii]):int(box_lat_end[ii]):nstep, int(box_lon_start[jj]):int(box_lon_end[jj]):nstep, k-1]
             atemp = atemp.flatten()
             b = (atemp < -0.0)
             if np.count_nonzero(b) > 0:
               atemp[b]=np.nan
             sst_avg[k,ii,jj] = np.nanmean(atemp)

  fid = Dataset(outfilename,'w')
  # Define the dimensions
  nx = fid.createDimension('nx', nx) # Unlimited
  ny = fid.createDimension('ny', ny) # Unlimited
  day = fid.createDimension('daynumber', 366) # Unlimited

  nc_var = fid.createVariable('daynumber', 'i4',('daynumber'),zlib=True)
  fid.variables['daynumber'][:] = range(1, 367)
  fid.variables['daynumber'].standard_name='day number'
  fid.variables['daynumber'].units='day number'
  fid.variables['daynumber'].axis='T'
  fid.variables['daynumber'].comment="Starting on Jan 1 of earch year with day number 1 until the end of the year with either day number 365 or 366 depending on the leap year status."

  nc_var = fid.createVariable('chlor_box_average', 'f8',('daynumber', 'nx', 'ny'),zlib=True)
  fid.variables['chlor_box_average'][:] = sst_avg[:,:,:]
  fid.variables['chlor_box_average'].standard_name='Data Average in the Box'
  fid.variables['chlor_box_average'].units='kelvin'

  fid.close()
 
def indexProcessing(dataset_name, startY, endY, lat_boxsize, lon_boxsize, lat_min, lat_max, lon_min, lon_max, lats, lons, variablename, error_variable, clim_filename, dirname, output_dir):
  
  nstep = 1
  nd = 366

  ### get box info
  nx, ny, box_lat_start, box_lat_end, box_lon_start, box_lon_end = data_info.getboxinfo(lat_boxsize, lon_boxsize, lat_min, lat_max, lon_min, lon_max, lats, lons)

  ### Read in climatology
  ncin = Dataset(clim_filename, 'r')
  clim = ncin.variables['sst_box_average'][:]
  ncin.close()

  ### setup size for tempory matrix
  nlat0 = math.ceil((box_lat_end[0]-box_lat_start[0])/nstep)
  nlon0 = math.ceil((box_lon_end[0]-box_lon_start[0])/nstep)

  for i in range(startY, endY+1):
    print("Processing Year = " + str(i)+", please wait ...")
    k = 0

    sst_year=np.empty((nd, nlat0, nlon0, nx, ny))
    sst_avg=np.empty((nd, nx, ny))
    sst_std=np.empty((nd, nx, ny))
    sst_anomaly=np.empty((nd, nx, ny))
    error_year=np.empty((nd, nlat0, nlon0, nx, ny))
    mask_year=np.empty((nd, nlat0, nlon0, nx, ny))

    sst_year[:] = np.nan
    sst_avg[:] = np.nan
    sst_std[:] = np.nan
    sst_anomaly[:] = np.nan
    error_year[:] = np.nan
    mask_year[:] = np.nan

    atime=np.empty(nd)

    for j in range(1, nd+1):
     for filename in glob.glob(dirname+'/'+str(i)+'/'+"{0:0>3}".format(j)+'/*.nc'):
       ncfile = filename.rsplit( "/")[ -1 ]
       ncin = Dataset(filename, 'r')
       sst = ncin.variables[variablename][:]
       error = ncin.variables[error_variable][:]
       mask = ncin.variables['mask'][:]
       for ii in range(0, nx):
          for jj in range(0, ny):
             sst_year[k, 0:len(range(int(box_lat_start[ii]),int(box_lat_end[ii]),nstep)), 0:len(range(int(box_lon_start[jj]),int(box_lon_end[jj]),nstep)), ii,jj] = sst[0,int(box_lat_start[ii]):int(box_lat_end[ii]):nstep, int(box_lon_start[jj]):int(box_lon_end[jj]):nstep]
             error_year[k, 0:len(range(int(box_lat_start[ii]),int(box_lat_end[ii]),nstep)), 0:len(range(int(box_lon_start[jj]),int(box_lon_end[jj]),nstep)), ii,jj] = error[0,int(box_lat_start[ii]):int(box_lat_end[ii]):nstep, int(box_lon_start[jj]):int(box_lon_end[jj]):nstep]
             mask_year[k, 0:len(range(int(box_lat_start[ii]),int(box_lat_end[ii]),nstep)), 0:len(range(int(box_lon_start[jj]),int(box_lon_end[jj]),nstep)), ii,jj] = mask[0,int(box_lat_start[ii]):int(box_lat_end[ii]):nstep, int(box_lon_start[jj]):int(box_lon_end[jj]):nstep]
             btemp = mask_year[k, :,:, ii,jj]
             atemp = sst_year[k, :,:, ii,jj]
             atemp = atemp.flatten()
             btemp = btemp.flatten()
             b = (atemp < -250.0) | (atemp > 350.0) | (btemp > 1)
             if np.count_nonzero(b) > 0:
               atemp[b]=np.nan
             sst_avg[k,ii,jj] = np.nanmean(atemp)
             sst_std[k,ii,jj] = np.nanstd(atemp)
             if j == 366:
               sst_anomaly[k,ii,jj] = sst_avg[k,ii,jj] - clim[364,ii,jj]
             else:
               sst_anomaly[k,ii,jj] = sst_avg[k,ii,jj] - clim[j-1,ii,jj]
       ncin.close()
       atime[k] = j
       k = k + 1

    ### create output directory
    data_info.createdir(output_dir+'/'+str(i))

    for ii in range(0, nx):
     for jj in range(0, ny):
        outfilename = output_dir + '/'+str(i)+'/'+variablename+'_'+str(i)+'_'+dataset_name+'_'+str(lat_boxsize)+'_degree_box_'+str(jj+1)+'x'+str(ii+1)+'.nc'
        fid = Dataset(outfilename,'w')
        # Define the dimensions
        nlat = fid.createDimension('lat', math.ceil((box_lat_end[0]-box_lat_start[0])/nstep)) # Unlimited
        nlon = fid.createDimension('lon', math.ceil((box_lon_end[0]-box_lon_start[0])/nstep)) # Unlimited
        time = fid.createDimension('time', k) # Unlimited

        nc_var = fid.createVariable('lat', 'f8',('lat'),zlib=True)
        fid.variables['lat'][:] = lats[int(box_lat_start[ii]):int(box_lat_end[ii]):nstep]
        fid.variables['lat'].standard_name='latitude'
        fid.variables['lat'].long_name='latitude'
        fid.variables['lat'].axis='Y'
        fid.variables['lat'].units='degrees_north'
        fid.variables['lat'].comment='Uniform grid with centers from '+str(lats[int(box_lat_start[ii])])+' to '+str(lats[int(box_lat_end[ii])])+ ' degrees.'

        nc_var = fid.createVariable('lon', 'f8',('lon'),zlib=True)
        fid.variables['lon'][:] = lons[int(box_lon_start[jj]):int(box_lon_end[jj]):nstep]
        fid.variables['lon'].standard_name='longitude'
        fid.variables['lon'].long_name='longitude'
        fid.variables['lon'].units='degrees_east'
        fid.variables['lon'].axis='X'
        fid.variables['lon'].comment='Uniform grid with centers from '+str(lons[int(box_lon_start[jj])])+' to '+str(lons[int(box_lon_end[jj])])+' degrees.'

        nc_var = fid.createVariable('time', 'i4',('time'),zlib=True)
        fid.variables['time'][:] = atime[0:k]
        fid.variables['time'].standard_name='time'
        fid.variables['time'].long_name='day number'
        fid.variables['time'].units='day number'
        fid.variables['time'].axis='T'
        fid.variables['time'].comment="Starting on Jan 1 of earch year with day number 1 until the end of the year with either day number 365 or 366 depending on the leap year status."

        nc_var = fid.createVariable('sst_box_average', 'f8',('time'),zlib=True)
        fid.variables['sst_box_average'][:] = sst_avg[0:k,ii,jj]
        fid.variables['sst_box_average'].standard_name='SST Average in the Box'
        fid.variables['sst_box_average'].units='kelvin'

        nc_var = fid.createVariable('sst_box_anomaly', 'f8',('time'),zlib=True)
        fid.variables['sst_box_anomaly'][:] = sst_anomaly[0:k,ii,jj]
        fid.variables['sst_box_anomaly'].standard_name='SST Anomaly in the Box'
        fid.variables['sst_box_anomaly'].units='kelvin'

        nc_var = fid.createVariable('sst_box_std', 'f8',('time'),zlib=True)
        fid.variables['sst_box_std'][:] = sst_std[0:k,ii,jj]
        fid.variables['sst_box_std'].standard_name='SST Standard Deviation in the Box'
        fid.variables['sst_box_std'].units='kelvin'

        nc_var = fid.createVariable('analysed_sst', 'f8',('time', 'lat', 'lon'),zlib=True,fill_value=-32768)
        fid.variables['analysed_sst'][:] = sst_year[0:k,:,:,ii,jj]
        fid.variables['analysed_sst'].units='kelvin'
        #fid.variables['analysed_sst']._FillValue=-32768
        fid.variables['analysed_sst'].comment='SST defined at all grid points but no physical meaning is ascribed to values over land'
        fid.variables['analysed_sst'].source='https://podaac.jpl.nasa.gov/dataset/CMC0.2deg-CMC-L4-GLOB-v2.0'

        nc_var = fid.createVariable('analysis_error', 'f8',('time', 'lat', 'lon'),zlib=True,fill_value=-32768)
        fid.variables['analysis_error'][:] = error_year[0:k,:,:,ii,jj]
        fid.variables['analysis_error'].units='kelvin'
        fid.variables['analysis_error'].comment='Error defined at all grid points but no meaning is ascribed to values over land'
        fid.variables['analysis_error'].source='https://podaac.jpl.nasa.gov/dataset/CMC0.2deg-CMC-L4-GLOB-v2.0'

        nc_var = fid.createVariable('mask', 'i',('time', 'lat', 'lon'),zlib=True,fill_value=-1)
        fid.variables['mask'][:] = mask_year[0:k,:,:,ii,jj]
        fid.variables['mask'].valid_min = 0
        fid.variables['mask'].valid_max = 31
        fid.variables['mask'].flag_masks = '1b, 2b, 4b, 8b, 16b'
        fid.variables['mask'].flag_meanings = "water land optional_lake_surface sea_ice optional_river_surface"
        fid.variables['mask'].comment='Mask can be used to further filter the data'
        fid.variables['mask'].source='https://podaac.jpl.nasa.gov/dataset/CMC0.2deg-CMC-L4-GLOB-v2.0'

        fid.product_version     = "Version 1.0"
        #fid.spatial_resolution  = str(dlat*nstep)+" degree"
        fid.Conventions         = "CF-1.6"
        #fid.title               = "CMC 0.2 global "+str(dlat*nstep)+" deg daily sea surface temperature analysis using Multi-Resolution Variational Analysis (MRVA) method for interpolation"
        fid.institution         = "Jet Propulsion Laboratory"
        fid.summary             = "Applies the method of statistical interpolation to assimilate observations from in situ and satellite sources using, as the background, the analysis valid 24 hours prior assuming persistence of the anomalies."
        fid.source              = 'https://podaac.jpl.nasa.gov/dataset/CMC0.2deg-CMC-L4-GLOB-v2.0'
        fid.creator_name        = "Ed Armstrong, Yibo Jiang"
        fid.creator_email       = "Edward.M.Armstrong@jpl.nasa.gov; Yibo.Jiang@jpl.nasa.gov"
        fid.creator_url         = "http;//podaac.jpl.nasa.gov"
        fid.project             = "Phenology of North Atlantic Ocean"
        fid.acknowledgment      = "This project was supported in part by a grant from"
        fid.processing_level    = "L4"
        #fid.date_created        = strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
        fid.gds_version_id      = "2.0r5"
        fid.naming_authority    = "org.ghrsst"
        #fid.geospatial_lat_resolution    = str(dlat*nstep)+"f"
        fid.geospatial_lat_units         = "degrees_north"
        #fid.geospatial_lon_resolution    = str(dlat*nstep)+"f"
        fid.geospatial_lon_units         = "degrees_east"
        fid.westernmost_longitude        = str( lons[int(box_lon_start[jj])] )
        fid.easternmost_longitude        = str( lons[int(box_lon_end[jj])] )
        fid.southernmost_latitude        = str( lats[int(box_lat_start[ii])] )
        fid.northernmost_latitude        = str( lats[int(box_lat_end[ii])] )
        fid.cdm_data_type       = "Grid"

        fid.close()

