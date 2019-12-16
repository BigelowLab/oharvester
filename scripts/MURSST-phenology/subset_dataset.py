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
# a skeleton script to download a set of L3 and L4 file using OPeNDAP.
#
#
#   2011.04.06  mike chin, version 0
#   2011.04.25  mike chin, version 1
#   2012.08.29  mike chin, version 2
#   2013.10.30  mike chin, version 2 (bug fix)
#   2014.12.23  Y. Jiang, version 3 (Modified to be more general including L3 and L4 datasets)
#   2017.06.20  Y. Jiang, fix to work with MODIS datasets 
#   2017.09.26  Y. Jiang, fix python 2.6 and 3.0 version issues 

##################################
# user parameters to be editted: #
##################################

# Caution: This is a Python script, and Python takes indentation seriously.
# DO NOT CHANGE INDENTATION OF ANY LINE BELOW!

# Here is the web-address (URL) for the GHRSST L4 file OPeNDAP server.
# You probably do not need to change this URL.  In fact, you can
# use this URL to search for the filename of the product you want to download.

# Here you set the product names in 3 parts:
#   ncHead = main directory of the product.
#   ncBody = main file name of the product (file name minus "yyyymmdd").
#   ncTail = the part that indicates the compression method.
#
# You also specify the grid dimension of the product, which you need to
# know in advance, by setting the following 6 variables:
#   nlon = x-grid (longitudes) dimension in integer.
#   nlat = y-grid (latitudes) dimension in integer.
#   dint = grid interval in degrees in float (assuming it's the same for x & y).
#   lon0 = the smallest longitude of the grid, in degrees in float.
#   lat0 = the smallest latitude of the grid, in degrees in float.
#   order = how "time" ,"lon", "lat" dimensions are ordered in the L4 file,
#           in array of three integers where 0=time, 1=lon, 2=lat.
#
# The code below might look complex, but it's only setting these
# 3 names (strings) and 6 variables listed above.

# Example:
# % ./subset_dataset.py -s 20100101 -f 20100201 -b -140 -110 20 30 -x MUR-JPL-L4-GLOB-v4.1
# Subset the data from 1 Jan 2010 to 1 Feb 2010 in a box from -140 to -110 degrees longitude and 20 to 30 degrees latitude
# for shortName MUR-JPL-L4-GLOB-v4.1

import sys,os
import time
import datetime
from datetime import date, timedelta
from math import floor,ceil
from optparse import OptionParser
from netCDF4 import Dataset
import calendar
import ftplib
import numpy as np

from xml.dom import minidom

if sys.version_info >= (3,0):
  import subprocess
  import urllib.request
else:
  import commands
  import urllib

#####################
# Global Parameters #
#####################

itemsPerPage = 10
PODAAC_WEB = 'https://podaac.jpl.nasa.gov'

###############
# subroutines #
###############
def createdir(dirname):
  if not os.path.exists(dirname):
    os.makedirs(dirname)

def strmatch(a,b):
  return (a in b) and (b in a)

def ncname(body,year,yday):
  (day,month)=calday(yday,year)
  return('%04d%02d%02d%s'%(year,month,day,body))  # for GDS1.x and GDS2.0

def pathname(body,tail,y,d):
  return( '%04d/%03d/%s%s' %(y,d,ncname(body,y,d),tail) )

def yearday(day,month,year):
  months=[0,31,28,31,30,31,30,31,31,30,31,30,31]
  if isLeap(year):
    months[2]=29
  for m in range(month):
    day=day+months[m]
  return(day)

def span(i1,i2,i3=1):
   return range(i1,i2+int(i3/abs(i3)),i3)

def isLeap(year):
  flag = ( (year%4)==0) and ( not ( (year%100)==0 and (year%400)!=0 ))
  return(flag)

def calday(yday,year):
  months=[0,31,28,31,30,31,30,31,31,30,31,30,31]
  if isLeap(year):
    months[2]=29
  for m in span(1,12):
    months[m]=months[m]+months[m-1]
  for m in span(1,12):
    if (yday-1)/months[m]==0:
      month=m
      day=yday-months[m-1]
      return (day,month)
  import sys
  sys.exit('ERROR calday: yearday value out of range')

def cal2mjd(year,month,day):
  import math
  a = (14 - month) // 12
  y = year + 4800 - a
  m = month + (12 * a) - 3
  p = day + (((153 * m) + 2) // 5) + (365 * y)
  q = (y // 4) - (y // 100) + (y // 400) - 32045
  return int( math.floor(p + q - 2400000.5) )

def mjd2cal(mjd):
  y=(mjd-409)//366+1860
  while cal2mjd(y+1,1,1)<=mjd:
      y=y+1
  m=1
  while cal2mjd(y,m+1,1)<=mjd:
      m=m+1
  d=1
  while cal2mjd(y,m,d+1)<=mjd:
      d=d+1
  return(y,m,d)

def today():
  import datetime
  todays=datetime.date.today()
  return str(todays.year)+str(todays.month).zfill(2)+str(todays.day).zfill(2)

def yesterday():
  import datetime
  yesterdays=datetime.date.today() - timedelta(days=1)
  return str(yesterdays.year)+str(yesterdays.month).zfill(2)+str(yesterdays.day).zfill(2)

def getChildrenByTitle(node):
    for child in node.childNodes:
        if child.localName=='Title':
            yield child
# find index set based on the given "box":
def boundingindex(sort,nt,dmin,dint,length,boundary0,boundary1):
  if (sort == 'A'):
    inx0=max(int(floor((boundary0-dmin)/dint)),0)
    if (inx0 > 0):
      inx0 = inx0 - 1
    inx1=min(int(floor((boundary1-dmin)/dint)),length-1)
  else:
    inx1=nt-max(int(floor((boundary0-dmin)/dint)),0)
    inx0=nt-min(int(ceil((boundary1-dmin)/dint)),length-1)
    if (inx0 > 0):
      inx0 = inx0 - 1

  return [inx0,inx1]

def download_data(start_year, start_month, start_day, end_year, end_month, end_day, latmin, latmax, lonmin, lonmax, gridpoints, shortname, rootdir):

  # get command line options:

  year0=start_year; month0=str(start_month).zfill(2); day0=str(start_day).zfill(2);
  year1=end_year; month1=str(end_month).zfill(2); day1=str(end_day).zfill(2);

  timeStr = '&startTime='+year0+'-'+month0+'-'+day0+'&endTime='+year1+'-'+month1+'-'+day1

  ig = gridpoints
  box = [lonmin, lonmax, latmin, latmax]

  print ('Longitude range: %f to %f'%(box[0],box[1]))
  print ('Latitude range: %f to %f'%(box[2],box[3]))

  wsurl = PODAAC_WEB+'/ws/search/granule/?shortName='+shortname+timeStr+'&itemsPerPage=1&sortBy=timeAsc'

  print ('\nPlease wait while program searching for the granules ...\n')

  if sys.version_info >= (3,0):
    response = urllib.request.urlopen(wsurl)
  else:
    response = urllib.urlopen(wsurl)
  data = response.read()

  if (len(data.splitlines()) == 1):
    sys.exit('No granules found for dataset: '+shortname+'\nProgram will exit now !\n')

  numGranules = 0
  doc = minidom.parseString(data)
  for arrays in doc.getElementsByTagName('link'):
   names = arrays.getAttribute("title")
   if names == 'OPeNDAP URL':
      numGranules = numGranules + 1
      href = arrays.getAttribute("href")
      #if numGranules > 0:
      #  break

  if numGranules == 0 and len(data.splitlines()) < 30:
    sys.exit('No granules found for dataset: '+shortname+'\nProgram will exit now !\n')
  elif numGranules == 0 and len(data.splitlines()) > 30:
    sys.exit('No OpenDap access for dataset: '+shortname+'\nProgram will exit now !\n')

  samplefile = href.rsplit( ".", 1 )[ 0 ] + '.ddx'

  variable_list = []
  lon_order = 'X'
  lat_order = 'Y'
  nt = 1 #time dimension
  nd = 1 #depth dimension

  if sys.version_info >= (3,0):
    doc = minidom.parse(urllib.request.urlopen(samplefile))
  else:
    doc = minidom.parse(urllib.urlopen(samplefile))
  for arrays in doc.getElementsByTagName('Grid'):
   names = arrays.getAttribute("name")
   if names == 'lat' or names == 'latitude':
     for dimensions in arrays.getElementsByTagName("dimension"):
       size = dimensions.getAttribute("size")
       name = dimensions.getAttribute("name")
       ni = int(size)
     for attrs in arrays.getElementsByTagName("Attribute"):
       aname = attrs.getAttribute("name")
       if aname == 'axis':
         for nodes in attrs.getElementsByTagName("value"):
           for cn in nodes.childNodes:
              lat_order = cn.nodeValue
   elif names == 'lon' or names == 'longitude':
     for dimensions in arrays.getElementsByTagName("dimension"):
       size = dimensions.getAttribute("size")
       name = dimensions.getAttribute("name")
       nj = int(size)
     for attrs in arrays.getElementsByTagName("Attribute"):
       aname = attrs.getAttribute("name")
       if aname == 'axis':
         for nodes in attrs.getElementsByTagName("value"):
           for cn in nodes.childNodes:
              lon_order = cn.nodeValue
   else:
     variable_list.append(names)
  for arrays in doc.getElementsByTagName('Map'):
   names = arrays.getAttribute("name")
   if names == 'lat' or names == 'latitude':
     for dimensions in arrays.getElementsByTagName("dimension"):
       size = dimensions.getAttribute("size")
       name = dimensions.getAttribute("name")
       ni = int(size)
     for attrs in arrays.getElementsByTagName("Attribute"):
       aname = attrs.getAttribute("name")
       if aname == 'axis':
         for nodes in attrs.getElementsByTagName("value"):
           for cn in nodes.childNodes:
              lat_order = cn.nodeValue
   if names == 'lon' or names == 'longitude':
     for dimensions in arrays.getElementsByTagName("dimension"):
       size = dimensions.getAttribute("size")
       name = dimensions.getAttribute("name")
       nj = int(size)
     for attrs in arrays.getElementsByTagName("Attribute"):
       aname = attrs.getAttribute("name")
       if aname == 'axis':
         for nodes in attrs.getElementsByTagName("value"):
           for cn in nodes.childNodes:
              lon_order = cn.nodeValue
   if names == 'time':
     ntime = 1
     for dimensions in arrays.getElementsByTagName("dimension"):
       size = dimensions.getAttribute("size")
       name = dimensions.getAttribute("name")
       nt = int(size)
   if names == 'depth':
     ndepth = 1
     for dimensions in arrays.getElementsByTagName("dimension"):
       size = dimensions.getAttribute("size")
       name = dimensions.getAttribute("name")
       nd = int(size)

  try:
     ni # does a exist in the current namespace
  except NameError:
    sys.exit('Granule file format may not be in netcdf or no latitude or longitude info for dataset: '+shortname+'\n')

  #****************************************************************************
  #*** Default southernmost_latitude, northernmost_latitude *******************
  #*** and westernmost_longitude, easternmost_longitude     *******************
  #****************************************************************************
  lat_sort = 'A'
  lon_sort = 'A'
  lat0 = -90.0
  lat1 = 90.0
  lon0 = -180.0
  lon1 = 180.0
  if "AQUARIUS" in shortname:
    lon0 = 0.0
    lon1 = 360.0
  if "MODIS" in shortname:
    lat_sort = 'D'

  #****************************************************************************
  for arrays in doc.getElementsByTagName('Attribute'):
   names = arrays.getAttribute("name")
   if names == 'southernmost_latitude':
     for nodes in arrays.getElementsByTagName("value"):
         for cn in nodes.childNodes:
           lat0 = float(cn.nodeValue)
   if names == 'northernmost_latitude':
     for nodes in arrays.getElementsByTagName("value"):
         for cn in nodes.childNodes:
           lat1 = float(cn.nodeValue)
   if names == 'westernmost_longitude':
     for nodes in arrays.getElementsByTagName("value"):
         for cn in nodes.childNodes:
           lon0 = float(cn.nodeValue)
   if names == 'easternmost_longitude':
     for nodes in arrays.getElementsByTagName("value"):
         for cn in nodes.childNodes:
           lon1 = float(cn.nodeValue)

  try:
     lat0 # does a exist in the current namespace
  except NameError:
    sys.exit('No southernmost_latitude info for dataset: '+shortname+'\n')

  try:
     lat1 # does a exist in the current namespace
  except NameError:
    sys.exit('No northernmost_latitude info for dataset: '+shortname+'\n')

  try:
     lon0 # does a exist in the current namespace
  except NameError:
    sys.exit('No westernmost_longitude info for dataset: '+shortname+'\n')

  try:
     lon1 # does a exist in the current namespace
  except NameError:
    sys.exit('No easternmost_longitude info for dataset: '+shortname+'\n')

  nlon = nj
  nlat = ni

  dint_lon = (lon1-lon0)/float(nlon)
  dint_lat = (lat1-lat0)/float(nlat)

  [i0,i1]=boundingindex(lon_sort,nlon,lon0,dint_lon,nlon,box[0],box[1])
  [j0,j1]=boundingindex(lat_sort,nlat,lat0,dint_lat,nlat,box[2],box[3])
  if i0>i1 or j0>j1:
    sys.exit('No grid point in your domain box.')

  # modify the max grid indices, as necessary:
  if ig>1:
    i1=max(span(i0,i1,ig))
    j1=max(span(j0,j1,ig))

  #************************************************************************************
  if lon_order == 'X':
    try:
       ntime # does a exist in the current namespace
       order=[0,2,1]
       try:
         ndepth # does a exist in the current namespace
         order=[0,1,3,2]
       except NameError:
         order=[0,2,1]
    except NameError:
       order=[1,0]
  else:
    try:
       ntime # does a exist in the current namespace
       order=[0,1,2]
       try:
         ndepth # does a exist in the current namespace
         order=[0,1,2,3]
       except NameError:
         order=[0,1,2]
    except NameError:
       order=[0,1]
  #************************************************************************************
  # download size information:
  print (' ')
  print ('Longitude range: %f to %f'%(box[0],box[1]))
  print ('Latitude range: %f to %f'%(box[2],box[3]))
  print ('  every %d pixel(s) is obtained'%(ig))
  print (' ')
  print ('grid dimensions will be ( %d x %d )'%(len(span(i0,i1,ig)),len(span(j0,j1,ig))))
  print (' ')

  # Check if curl or wget commands exsit on your computer
  if sys.version_info >= (3,0):
    status_curl, result = subprocess.getstatusoutput('which curl')
    status_wget, result = subprocess.getstatusoutput('which wget')
  else:
    status_curl, result = commands.getstatusoutput("which curl")
    status_wget, result = commands.getstatusoutput("which wget")

  # form the index set for the command line:
  try:
     ntime # does a exist in the current namespace
     inx=[[0,1,nt-1],[i0,ig,i1],[j0,ig,j1]]
     try:
       ndepth # does a exist in the current namespace
       inx=[[0,1,nt-1],[0,1,nd-1],[i0,ig,i1],[j0,ig,j1]]
     except NameError:
       inx=[[0,1,nt-1],[i0,ig,i1],[j0,ig,j1]]
  except NameError:
     inx=[[i0,ig,i1],[j0,ig,j1]]
     try:
       ndepth # does a exist in the current namespace
       inx=[[0,1,nd-1],[i0,ig,i1],[j0,ig,j1]]
     except NameError:
       inx=[[i0,ig,i1],[j0,ig,j1]]

  index=''
  for i in order:
    index=index+'[%d:%d:%d]'%(inx[i][0],inx[i][1],inx[i][2])

  # main loop:
  start = time.time()
  bmore = 1
  bwrite = False
  while (bmore > 0):
   if (bmore == 1):
       urllink = PODAAC_WEB+'/ws/search/granule/?shortName='+shortname+timeStr+'&itemsPerPage=%d&sortBy=timeAsc'%itemsPerPage
   else:
       urllink = PODAAC_WEB+'/ws/search/granule/?shortName='+shortname+timeStr+'&itemsPerPage=%d&sortBy=timeAsc&startIndex=%d'%(itemsPerPage, (bmore-1)*itemsPerPage)
   bmore = bmore + 1
   if sys.version_info >= (3,0):
     response = urllib.request.urlopen(urllink)
   else:
     response = urllib.urlopen(urllink)
   data = response.read()
   doc = minidom.parseString(data)

   numGranules = 0
   for arrays in doc.getElementsByTagName('link'):
    names = arrays.getAttribute("title")
    if names == 'OPeNDAP URL':
      numGranules = numGranules + 1
      href = arrays.getAttribute("href")
      ncfile = href.rsplit( ".", 1 )[ 0 ]
      head, tail = os.path.split(ncfile)
      ncout = tail
      if ncout.endswith('.bz2') or ncout.endswith('.gz'):
        ncout = ncout.rsplit( ".", 1 )[ 0 ]
      ncout = ncout.rsplit( ".", 1 )[ 0 ]+'_subset.'+ncout.rsplit( ".", 1 )[ 1 ]
      day_of_year = datetime.date(int(ncout[0:4]), int(ncout[4:6]), int(ncout[6:8])).timetuple().tm_yday
      outdir = rootdir + '/' + ncout[0:4] + '/' + str(day_of_year).zfill(3)
      ncout = outdir + '/' + ncout
      createdir(outdir)
      cmd=ncfile+'.nc?'
      for item in variable_list:
        cmd=cmd+item+index+','
      cmd=cmd[0:(len(cmd)-1)]  # remove the extra "," at the end.

      if status_curl == 0:
        #cmd='curl -g -silent "'+cmd+'" -o '+ ncout
        cmd='curl -g  "'+cmd+'" -o '+ ncout
      elif status_wget == 0:
        cmd='wget -quite "'+cmd+'" -O '+ ncout
      else:
        sys.exit('\nThe script will need curl or wget on the system, please install them first before running the script !\nProgram will exit now !\n')

      os.system( cmd )
      print (ncout + ' download finished !')
      if bwrite == False:
        ncin = Dataset(ncout, 'r')
        lons = ncin.variables['lon'][:]
        lats = ncin.variables['lat'][:]
        ncin.close()

        outfilename = rootdir+"/"+shortname + '_info.nc'
        fid = Dataset(outfilename,'w')
        # Define the dimensions
        nlat = fid.createDimension('lat', lats.size) # Unlimited
        nlon = fid.createDimension('lon', lons.size) # Unlimited

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

        fid.cdm_data_type = "Grid"
        fid.close()
        bwrite = True
  
   if numGranules < itemsPerPage:
     bmore = 0

  end = time.time()
  print ('Time spend = ' + str(end - start) + ' seconds')


def download_data_seawifs(start_year, start_month, start_day, end_year, end_month, end_day, latmin, latmax, lonmin, lonmax, gridpoints, rootdir):

  GSFC_WEB = 'https://oceandata.sci.gsfc.nasa.gov:443/opendap/SeaWiFS/L3SMI/'

  # get command line options:

  year0=int(start_year); month0=str(start_month).zfill(2); day0=str(start_day).zfill(2);
  year1=int(end_year); month1=str(end_month).zfill(2); day1=str(end_day).zfill(2);

  box = [lonmin, lonmax, latmin, latmax]

  lat0 = -89.95834
  lat1 = 89.95834

  lon0 = -179.958
  lon1 = 179.958

  nlon = 4320
  nlat = 2160

  dint_lon = 0.083
  dint_lat = 0.0833

  i0=int(floor((lonmin-lon0)/dint_lon))
  i1=int(ceil((lonmax-lon0)/dint_lon))
  j0=int(floor((lat1-latmax)/dint_lat))
  j1=int(ceil((lat1-latmin)/dint_lat))

  print(i0)
  print(i1)
  print(j0)
  print(j1)
  #************************************************************************************
  # download size information:
  print (' ')
  print ('Longitude range: %f to %f'%(box[0],box[1]))
  print ('Latitude range: %f to %f'%(box[2],box[3]))
  print (' ')

  # Check if curl or wget commands exsit on your computer
  if sys.version_info >= (3,0):
    status_curl, result = subprocess.getstatusoutput('which curl')
    status_wget, result = subprocess.getstatusoutput('which wget')
  else:
    status_curl, result = commands.getstatusoutput("which curl")
    status_wget, result = commands.getstatusoutput("which wget")

  lat_index = '[%d:%d:%d]'%(j0,1,j1)
  lon_index = '[%d:%d:%d]'%(i0,1,i1)

  # main loop:
  start = time.time()

  day_of_year0 = datetime.date(int(start_year), int(start_month), int(start_day)).timetuple().tm_yday
  day_of_year1 = datetime.date(int(end_year), int(end_month), int(end_day)).timetuple().tm_yday

  for i in range(year0, year1+1):
    range_start = 1
    range_end = 367
    if(calendar.isleap(i)):
      range_end = 367
    else:
      range_end = 366
    if (i == year0):
      range_start = day_of_year0
    elif (i == year1):
      range_end = day_of_year1

    for j in range(range_start, range_end):
      filename = 'S'+str(i)+str(j).zfill(3)+'.L3m_DAY_CHL_chlor_a_9km.nc.nc'
      filename0 = 'S'+str(i)+str(j).zfill(3)+'.L3m_DAY_CHL_chlor_a_9km.nc'
      outdir = rootdir + '/' + str(i) + '/' + str(j).zfill(3)
      ncout = outdir + '/' + filename0
      createdir(outdir)

      cmd = GSFC_WEB+str(i)+'/'+str(j).zfill(3)+'/'+filename+'?chlor_a'+lat_index+lon_index+',lat'+lat_index+',lon'+lon_index+',palette[0:1:2][0:1:255]'

      if status_curl == 0:
        #cmd='curl -g -silent "'+cmd+'" -o '+ ncout
        cmd='curl -g  "'+cmd+'" -o '+ ncout
      elif status_wget == 0:
        cmd='wget -quite "'+cmd+'" -O '+ ncout
      else:
        sys.exit('\nThe script will need curl or wget on the system, please install them first before running the script !\nProgram will exit now !\n')

      os.system( cmd )
      print (ncout + ' download finished !')

      # remove file if file does not exist or damaged
      bsize = os.path.getsize(ncout)
      if(bsize < 1000):
        os.remove(ncout)
        os.rmdir(outdir)
  
  end = time.time()
  print ('Time spend = ' + str(end - start) + ' seconds')


def download_data_modis_bck(start_year, start_month, start_day, end_year, end_month, end_day, latmin, latmax, lonmin, lonmax, gridpoints, rootdir):

  GSFC_WEB = 'https://oceandata.sci.gsfc.nasa.gov:443/opendap/MODISA/L3SMI/'

  # get command line options:

  year0=int(start_year); month0=str(start_month).zfill(2); day0=str(start_day).zfill(2);
  year1=int(end_year); month1=str(end_month).zfill(2); day1=str(end_day).zfill(2);

  box = [lonmin, lonmax, latmin, latmax]

  lat0 = -89.95834
  lat1 = 89.95834

  lon0 = -179.958
  lon1 = 179.958

  nlon = 4320
  nlat = 2160

  dint_lon = 0.083
  dint_lat = 0.0833

  i0=int(floor((lonmin-lon0)/dint_lon))
  i1=int(ceil((lonmax-lon0)/dint_lon))
  j0=int(floor((lat1-latmax)/dint_lat))
  j1=int(ceil((lat1-latmin)/dint_lat))

  print(i0)
  print(i1)
  print(j0)
  print(j1)
  #************************************************************************************
  # download size information:
  print (' ')
  print ('Longitude range: %f to %f'%(box[0],box[1]))
  print ('Latitude range: %f to %f'%(box[2],box[3]))
  print (' ')

  # Check if curl or wget commands exsit on your computer
  if sys.version_info >= (3,0):
    status_curl, result = subprocess.getstatusoutput('which curl')
    status_wget, result = subprocess.getstatusoutput('which wget')
  else:
    status_curl, result = commands.getstatusoutput("which curl")
    status_wget, result = commands.getstatusoutput("which wget")

  lat_index = '[%d:%d:%d]'%(j0,1,j1)
  lon_index = '[%d:%d:%d]'%(i0,1,i1)

  # main loop:
  start = time.time()

  day_of_year0 = datetime.date(int(start_year), int(start_month), int(start_day)).timetuple().tm_yday
  day_of_year1 = datetime.date(int(end_year), int(end_month), int(end_day)).timetuple().tm_yday

  for i in range(year0, year1+1):
    range_start = 1
    range_end = 367
    if(calendar.isleap(i)):
      range_end = 367
    else:
      range_end = 366
    if (i == year0):
      range_start = day_of_year0
    elif (i == year1):
      range_end = day_of_year1

    for j in range(range_start, range_end):
      filename = 'A'+str(i)+str(j).zfill(3)+'.L3m_DAY_CHL_chlor_a_4km.nc.nc'
      filename0 = 'A'+str(i)+str(j).zfill(3)+'.L3m_DAY_CHL_chlor_a_4km.nc'
      outdir = rootdir + '/' + str(i) + '/' + str(j).zfill(3)
      ncout = outdir + '/' + filename0
      createdir(outdir)

      cmd = GSFC_WEB+str(i)+'/'+str(j).zfill(3)+'/'+filename+'?chlor_a'+lat_index+lon_index+',lat'+lat_index+',lon'+lon_index+',palette[0:1:2][0:1:255]'

      if status_curl == 0:
        #cmd='curl -g -silent "'+cmd+'" -o '+ ncout
        cmd='curl -g  "'+cmd+'" -o '+ ncout
      elif status_wget == 0:
        cmd='wget -quite "'+cmd+'" -O '+ ncout
      else:
        sys.exit('\nThe script will need curl or wget on the system, please install them first before running the script !\nProgram will exit now !\n')

      os.system( cmd )
      print (ncout + ' download finished !')

      # remove file if file does not exist or damaged
      bsize = os.path.getsize(ncout)
      if(bsize < 1000):
        os.remove(ncout)
        os.rmdir(outdir)
  
  end = time.time()
  print ('Time spend = ' + str(end - start) + ' seconds')


def download_data_modis(start_year, start_month, start_day, end_year, end_month, end_day, latmin, latmax, lonmin, lonmax, gridpoints, rootdir):

  #GSFC_WEB = 'https://oceandata.sci.gsfc.nasa.gov/MODIS-Aqua/Mapped/Daily/4km/chlor_a/'
  GSFC_WEB = "https://oceandata.sci.gsfc.nasa.gov/opendap/MODISA/L3SMI/"
  GSFC_WEB = "https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/"

  # get command line options:

  year0=int(start_year); month0=str(start_month).zfill(2); day0=str(start_day).zfill(2);
  year1=int(end_year); month1=str(end_month).zfill(2); day1=str(end_day).zfill(2);

  box = [lonmin, lonmax, latmin, latmax]

  lat0 = -89.95834
  lat1 = 89.95834

  lon0 = -179.958
  lon1 = 179.958

  nlon = 4320
  nlat = 2160

  dint_lon = 0.083
  dint_lat = 0.0833

  i0=int(floor((lonmin-lon0)/dint_lon))
  i1=int(ceil((lonmax-lon0)/dint_lon))
  j0=int(floor((lat1-latmax)/dint_lat))
  j1=int(ceil((lat1-latmin)/dint_lat))

  print(i0)
  print(i1)
  print(j0)
  print(j1)
  #************************************************************************************
  # download size information:
  print (' ')
  print ('Longitude range: %f to %f'%(box[0],box[1]))
  print ('Latitude range: %f to %f'%(box[2],box[3]))
  print (' ')

  # Check if curl or wget commands exsit on your computer
  if sys.version_info >= (3,0):
    status_curl, result = subprocess.getstatusoutput('which curl')
    status_wget, result = subprocess.getstatusoutput('which wget')
  else:
    status_curl, result = commands.getstatusoutput("which curl")
    status_wget, result = commands.getstatusoutput("which wget")

  lat_index = '[%d:%d:%d]'%(j0,1,j1)
  lon_index = '[%d:%d:%d]'%(i0,1,i1)

  # main loop:
  start = time.time()

  day_of_year0 = datetime.date(int(start_year), int(start_month), int(start_day)).timetuple().tm_yday
  day_of_year1 = datetime.date(int(end_year), int(end_month), int(end_day)).timetuple().tm_yday

  cmd = 'wget -q -O - "https://oceandata.sci.gsfc.nasa.gov/api/file_search?sensor=aqua&sdate='+start_year+'-'+str(int(start_month)).zfill(2)+'-'+str(int(start_day)).zfill(2)+'&edate='+end_year+'-'+str(int(end_month)).zfill(2)+'-'+str(int(end_day)).zfill(2)+'&dtype=L3m&search=A*.L3m_DAY_CHL_chlor_a_4km.nc&format=txt"'

  filelist = os.popen(cmd).read().split()

  for afile in filelist:
    outdir = rootdir + '/'+afile[1:5] + '/' + afile[5:8]
    ncout = outdir + '/'+afile
    createdir(outdir)

    #cmd = GSFC_WEB+str(i)+'/'+filename
    acmd = GSFC_WEB+afile

    if status_curl == 0:
        acmd='curl "'+acmd+'" -o '+ ncout
    elif status_wget == 0:
        acmd='wget "'+acmd+'" -O '+ ncout
    else:
        sys.exit('\nThe script will need curl or wget on the system, please install them first before running the script !\nProgram will exit now !\n')

    print(acmd)
    os.system( acmd )
    print (ncout + ' download finished !')

    # remove file if file does not exist or damaged
    bsize = os.path.getsize(ncout)
    if(bsize < 1000):
      os.remove(ncout)
      os.rmdir(outdir)
  
  end = time.time()
  print ('Time spend = ' + str(end - start) + ' seconds')



