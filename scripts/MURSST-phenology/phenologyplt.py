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
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.ticker import AutoMinorLocator
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
import datetime
import calendar
from pathlib import Path

import data_info

################################################################################################   

def plot_grid_box(startY, endY, lat_boxsize, lon_boxsize, lat_min, lat_max, lon_min, lon_max, lats, lons, variablename, error_variable, clim_filename, dirname):
  ### get box info
  nx, ny, box_lat_start, box_lat_end, box_lon_start, box_lon_end = data_info.getboxinfo(lat_boxsize, lon_boxsize, lat_min, lat_max, lon_min, lon_max, lats, lons)

  cmap=plt.get_cmap("jet")
  cmap.set_under("black")
  cmap.set_over("DarkViolet")
  clevs = np.linspace(270, 310, 21)

  parallels = np.arange(-90,90+5,5)
  meridians = np.arange(0,360+5,5)

  fig= plt.figure(figsize=(6,4), dpi=300)

  map = Basemap(projection='cyl', resolution='i', llcrnrlon=lon_min, llcrnrlat=lat_min, urcrnrlon=lon_max, urcrnrlat=lat_max)

  map.drawcoastlines(linewidth=0.25)
  map.drawcountries(linewidth=0.25)
  map.fillcontinents(color='0.8')
  map.drawparallels(parallels,labels=[1,0,0,0],fontsize=6)
  map.drawmeridians(meridians,labels=[1,0,0,1],fontsize=6)

  for i in range(1, ny+1):
   for j in range(1, nx+1):
       x, y = map([lon_min+(i-1)*lon_boxsize, lon_min+i*lon_boxsize, lon_min+i*lon_boxsize, lon_min+(i-1)*lon_boxsize, lon_min+(i-1)*lon_boxsize], [lat_min+j*lat_boxsize, lat_min+j*lat_boxsize,lat_min+(j-1)*lat_boxsize,lat_min+(j-1)*lat_boxsize,lat_min+j*lat_boxsize])

       map.plot(x,y, color='r', linewidth=2)

       plt.text((min(x)+max(x))*0.5, (min(y)+max(y))*0.5, str(j)+'x'+str(i),fontsize=14,fontweight='bold', ha='center',va='center',color='r')

  plt.title('Data Region Division ('+str(lat_boxsize)+'x'+str(lon_boxsize)+' degree)', fontsize=10, fontweight='bold')

  plt.show()

def plot_clim_box_contour(lat_boxsize, lon_boxsize, lat_min, lat_max, lon_min, lon_max, lats, lons, variablename, error_variable, clim_filename):
  ### get box info
  nx, ny, box_lat_start, box_lat_end, box_lon_start, box_lon_end = data_info.getboxinfo(lat_boxsize, lon_boxsize, lat_min, lat_max, lon_min, lon_max, lats, lons)

  #*** read climatology file
  ncin = Dataset(clim_filename, 'r')
  data = ncin.variables['sst_climatology'][:]
  lons = ncin.variables['lon'][:]
  lats = ncin.variables['lat'][:]
  ncin.close()

  parallels = np.arange(-90,90+5,5)
  meridians = np.arange(0,360+5,5)

  cmap=plt.get_cmap("jet")
  cmap.set_under("black")
  cmap.set_over("DarkViolet")
  clevs = np.linspace(270, 310, 21)

  fig= plt.figure(figsize=(6,4), dpi=300)

  map = Basemap(projection='cyl', resolution='i', llcrnrlon=lon_min, llcrnrlat=lat_min, urcrnrlon=lon_max, urcrnrlat=lat_max)

  map.drawcoastlines(linewidth=0.25)
  map.drawcountries(linewidth=0.25)
  map.fillcontinents(color='0.8')
  map.drawparallels(parallels,labels=[1,0,0,0],fontsize=6)
  map.drawmeridians(meridians,labels=[1,0,0,1],fontsize=6)

  x, y = map(*np.meshgrid(lons, lats))
  cs=map.contourf(x, y, data[:,:,0].squeeze(), clevs, cmap=cmap, extend='both')

  cb = map.colorbar(cs, 'right', size='2%', pad='0.5%')
  cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=6)
  cb.set_label('SST (K)', fontsize=7,fontweight="bold")
  cb.set_ticks(range(270,315,5))

  for i in range(1, ny+1):
   for j in range(1, nx+1):

       x, y = map([lon_min+(i-1)*lon_boxsize, lon_min+i*lon_boxsize, lon_min+i*lon_boxsize, lon_min+(i-1)*lon_boxsize, lon_min+(i-1)*lon_boxsize], [lat_min+j*lat_boxsize, lat_min+j*lat_boxsize,lat_min+(j-1)*lat_boxsize,lat_min+(j-1)*lat_boxsize,lat_min+j*lat_boxsize])

       map.plot(x,y, color='r', linewidth=1)

       x, y = map([lon_min+(i-1)*lon_boxsize, lon_min+i*lon_boxsize, lon_min+i*lon_boxsize, lon_min+(i-1)*lon_boxsize, lon_min+(i-1)*lon_boxsize], [lat_min+j*lat_boxsize, lat_min+j*lat_boxsize,lat_min+(j-1)*lat_boxsize,lat_min+(j-1)*lat_boxsize,lat_min+j*lat_boxsize])

       map.plot(x,y, color='r', linewidth=2)

       plt.text((min(x)+max(x))*0.5, (min(y)+max(y))*0.5, str(j)+'x'+str(i),fontsize=14,fontweight='bold', ha='center',va='center',color='r')

  plt.title('Data Region Division (5 degree)', fontsize=10, fontweight='bold')

  plt.show()

def plot_line_box_climatology(lat_boxsize, lon_boxsize, lat_min, lat_max, lon_min, lon_max, lats, lons, variablename, clim_filename, outfilename):

  ### get box info
  nx, ny, box_lat_start, box_lat_end, box_lon_start, box_lon_end = data_info.getboxinfo(lat_boxsize, lon_boxsize, lat_min, lat_max, lon_min, lon_max, lats, lons)

  # Read in climatology
  ncin = Dataset(outfilename, 'r')
  daynumber = ncin.variables['daynumber'][:]
  sst_box_average = ncin.variables['sst_box_average'][:]
  ncin.close()

  ### setup plot
  fig = plt.figure(figsize=(1.5 * ny, 1.2 * nx))
  fig.suptitle('Climatology in Each Grid Box ('+str(lat_boxsize)+'x'+str(lon_boxsize)+' degree)', fontsize=16, fontweight='bold')
  plt.subplots_adjust(bottom=0.05, left=.03, right=.99, top=.95, hspace=.02)
  
  cplot = 1
  for j in range(nx-1, -1,-1):
    for i in range(0, ny):
      ax=plt.subplot(nx, ny, cplot)
      ax.grid(True)
      plt.plot(daynumber, sst_box_average[:,j,i])
      plt.xlim(np.min(daynumber), np.max(daynumber))
      plt.ylim(np.nanmin(sst_box_average), np.nanmax(sst_box_average))
      if i == 0:
        plt.ylabel(variablename, fontsize='medium')
      if j == 0:
        plt.xlabel('Day Number', fontsize='medium')
      if i > 0:
        ax.set_yticklabels([])
      if j != 0:
        ax.set_xticklabels([])

      plt.text((np.min(daynumber)+np.max(daynumber))*0.5, (np.nanmin(sst_box_average)+np.nanmax(sst_box_average))*0.5, str(j+1)+'x'+str(i+1),
               fontsize=10,fontweight='bold', ha='center',va='center',color='r')
      cplot=cplot+1

  plt.show()

def plot_line_box_climatology_chlor(lat_boxsize, lon_boxsize, lat_min, lat_max, lon_min, lon_max, lats, lons, variablename, clim_filename, outfilename):

  ### get box info
  nx, ny, box_lat_start, box_lat_end, box_lon_start, box_lon_end = data_info.getboxinfo(lat_boxsize, lon_boxsize, lat_min, lat_max, lon_min, lon_max, lats, lons)

  # Read in climatology
  ncin = Dataset(outfilename, 'r')
  daynumber = ncin.variables['daynumber'][:]
  sst_box_average = ncin.variables['chlor_box_average'][:]
  ncin.close()

  ### setup plot
  fig = plt.figure(figsize=(1.5 * ny, 1.2 * nx))
  fig.suptitle('Climatology in Each Grid Box ('+str(lat_boxsize)+'x'+str(lon_boxsize)+' degree)', fontsize=16, fontweight='bold')
  plt.subplots_adjust(bottom=0.05, left=.03, right=.99, top=.95, hspace=.02)
  
  cplot = 1
  for j in range(nx-1, -1,-1):
    for i in range(0, ny):
      ax=plt.subplot(nx, ny, cplot)
      ax.grid(True)
      plt.plot(daynumber, sst_box_average[:,j,i])
      plt.xlim(np.min(daynumber), np.max(daynumber))
      plt.ylim(np.nanmin(sst_box_average), np.nanmax(sst_box_average))
      if i == 0:
        plt.ylabel(variablename, fontsize='medium')
      if j == 0:
        plt.xlabel('Day Number', fontsize='medium')
      if i > 0:
        ax.set_yticklabels([])
      if j != 0:
        ax.set_xticklabels([])

      plt.text((np.min(daynumber)+np.max(daynumber))*0.5, (np.nanmin(sst_box_average)+np.nanmax(sst_box_average))*0.5, str(j+1)+'x'+str(i+1),
               fontsize=10,fontweight='bold', ha='center',va='center',color='r')
      cplot=cplot+1

  plt.show()

def plot_clim_box_contour_chlor(lat_boxsize, lon_boxsize, lat_min, lat_max, lon_min, lon_max, lats, lons, variablename, error_variable, clim_filename):
  ### get box info
  nx, ny, box_lat_start, box_lat_end, box_lon_start, box_lon_end = data_info.getboxinfo(lat_boxsize, lon_boxsize, lat_min, lat_max, lon_min, lon_max, lats, lons)

  #*** read climatology file
  ncin = Dataset(clim_filename, 'r')
  data = ncin.variables['chlor_climatology'][:]
  #data = ncin.variables['chlor_num_days'][:]
  lons = ncin.variables['lon'][:]
  lats = ncin.variables['lat'][:]
  ncin.close()

  parallels = np.arange(-90,90+5,5)
  meridians = np.arange(0,360+5,5)

  cmap=plt.get_cmap("jet")
  cmap.set_under("black")
  cmap.set_over("DarkViolet")
  #clevs = np.linspace(0, 10, 21)
  clevs = np.linspace(0, 10, 21)

  fig= plt.figure(figsize=(6,4), dpi=300)

  map = Basemap(projection='cyl', resolution='i', llcrnrlon=lon_min, llcrnrlat=lat_min, urcrnrlon=lon_max, urcrnrlat=lat_max)

  map.drawcoastlines(linewidth=0.25)
  map.drawcountries(linewidth=0.25)
  map.fillcontinents(color='0.8')
  map.drawparallels(parallels,labels=[1,0,0,0],fontsize=6)
  map.drawmeridians(meridians,labels=[1,0,0,1],fontsize=6)

  x, y = map(*np.meshgrid(lons, lats))
  cs=map.contourf(x, y, data[:,:,0].squeeze(), clevs, cmap=cmap, extend='both')

  cb = map.colorbar(cs, 'right', size='2%', pad='0.5%')
  cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=6)
  cb.set_label('Chlorophyll (mg m^-3)', fontsize=7,fontweight="bold")
  cb.set_ticks(range(0,10,1))

  for i in range(1, ny+1):
   for j in range(1, nx+1):

       x, y = map([lon_min+(i-1)*lon_boxsize, lon_min+i*lon_boxsize, lon_min+i*lon_boxsize, lon_min+(i-1)*lon_boxsize, lon_min+(i-1)*lon_boxsize], [lat_min+j*lat_boxsize, lat_min+j*lat_boxsize,lat_min+(j-1)*lat_boxsize,lat_min+(j-1)*lat_boxsize,lat_min+j*lat_boxsize])

       map.plot(x,y, color='r', linewidth=1)

       x, y = map([lon_min+(i-1)*lon_boxsize, lon_min+i*lon_boxsize, lon_min+i*lon_boxsize, lon_min+(i-1)*lon_boxsize, lon_min+(i-1)*lon_boxsize], [lat_min+j*lat_boxsize, lat_min+j*lat_boxsize,lat_min+(j-1)*lat_boxsize,lat_min+(j-1)*lat_boxsize,lat_min+j*lat_boxsize])

       map.plot(x,y, color='r', linewidth=2)

       plt.text((min(x)+max(x))*0.5, (min(y)+max(y))*0.5, str(j)+'x'+str(i),fontsize=14,fontweight='bold', ha='center',va='center',color='r')

  plt.title('Data Region Division (5 degree)', fontsize=10, fontweight='bold')

  plt.show()

def plot_line_box_climatology(lat_boxsize, lon_boxsize, lat_min, lat_max, lon_min, lon_max, lats, lons, variablename, clim_filename, outfilename):

  ### get box info
  nx, ny, box_lat_start, box_lat_end, box_lon_start, box_lon_end = data_info.getboxinfo(lat_boxsize, lon_boxsize, lat_min, lat_max, lon_min, lon_max, lats, lons)

  # Read in climatology
  ncin = Dataset(outfilename, 'r')
  daynumber = ncin.variables['daynumber'][:]
  sst_box_average = ncin.variables['sst_box_average'][:]
  ncin.close()

  ### setup plot
  fig = plt.figure(figsize=(1.5 * ny, 1.2 * nx))
  fig.suptitle('Climatology in Each Grid Box ('+str(lat_boxsize)+'x'+str(lon_boxsize)+' degree)', fontsize=16, fontweight='bold')
  plt.subplots_adjust(bottom=0.05, left=.03, right=.99, top=.95, hspace=.02)
  
  cplot = 1
  for j in range(nx-1, -1,-1):
    for i in range(0, ny):
      ax=plt.subplot(nx, ny, cplot)
      ax.grid(True)
      plt.plot(daynumber, sst_box_average[:,j,i])
      plt.xlim(np.min(daynumber), np.max(daynumber))
      plt.ylim(np.nanmin(sst_box_average), np.nanmax(sst_box_average))
      if i == 0:
        plt.ylabel(variablename, fontsize='medium')
      if j == 0:
        plt.xlabel('Day Number', fontsize='medium')
      if i > 0:
        ax.set_yticklabels([])
      if j != 0:
        ax.set_xticklabels([])

      plt.text((np.min(daynumber)+np.max(daynumber))*0.5, (np.nanmin(sst_box_average)+np.nanmax(sst_box_average))*0.5, str(j+1)+'x'+str(i+1),
               fontsize=10,fontweight='bold', ha='center',va='center',color='r')
      cplot=cplot+1

  plt.show()

def plot_annual_maximum(lat_min, lat_max, lon_min, lon_max, lats, lons, variablename, metric_dir, metric_filename):

  # Read in climatology
  my_file = Path(metric_dir+"/"+metric_filename)
  if not my_file.is_file():
    print("\nPlease create the Annual Maximum first !\n")
    return None
  ncin = Dataset(metric_dir+"/"+metric_filename, 'r')
  data_max = ncin.variables['data_max'][:]
  day_max = ncin.variables['day_max'][:]
  year = ncin.variables['year'][:]
  lons = ncin.variables['lon'][:]
  lats = ncin.variables['lat'][:]
  ncin.close()

  ny = year.size

  for i in range(0,year.size):

    if(i%5 == 0):
      # setup plot
      fig = plt.figure(figsize=(3 * 3, 3 * 4))
      fig.suptitle(variablename + ' Annual Maximum and Its Day Number', fontsize=14, fontweight='bold')
      plt.subplots_adjust(left=0.06, right=0.92, top=0.94, bottom=0.05, hspace=.1, wspace=.4)
      nf = 1

    plt.subplot(5,2,nf)

    m = Basemap(projection='cyl', llcrnrlon=lon_min, llcrnrlat=lat_min, urcrnrlon=lon_max, urcrnrlat=lat_max)
    x, y = m(*np.meshgrid(lons, lats))
    clevs = np.linspace(np.nanmin(data_max[:,:,i].flatten()), np.nanmax(data_max[:,:,i].flatten()), 21)
    cs=m.contourf(x, y, data_max[:,:,i].squeeze(), clevs, cmap=plt.cm.RdBu_r)
    m.drawcoastlines()
    m.fillcontinents(color='#000000',lake_color='#99ffff')

    # draw parallels and meridians.
    # label parallels on right and top
    # meridians on bottom and left
    parallels = np.arange(lat_min,lat_max,10.)
    # labels = [left,right,top,bottom]
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
    meridians = np.arange(lon_min,lon_max,10.)
    meri = m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)

    cb = plt.colorbar(cs, orientation='vertical')
    cb.set_label('Maximum Value', fontsize=8)
    cb.ax.tick_params(labelsize=8)

    plt.text( min(lons)+0.8*(max(lons)-min(lons)), min(lats)+0.2*(max(lats)-min(lats)), str(year[i]),fontsize=8,fontweight='bold', ha='left',va='center',color='r', bbox=dict(boxstyle="round", ec=(0.8, 0.61, 0.11), fc=(0.8, 0.61, 0.11),))

    nf = nf + 1
    plt.subplot(5,2,nf)

    m = Basemap(projection='cyl', llcrnrlon=lon_min, llcrnrlat=lat_min, urcrnrlon=lon_max, urcrnrlat=lat_max)
    x, y = m(*np.meshgrid(lons, lats))
    clevs = np.linspace(np.nanmin(day_max[:,:,i].flatten()), np.nanmax(day_max[:,:,i].flatten()), 21)
    cs=m.contourf(x, y, day_max[:,:,i].squeeze(), clevs, cmap=plt.cm.RdBu_r)
    m.drawcoastlines()
    m.fillcontinents(color='#000000',lake_color='#99ffff')

    # draw parallels and meridians.
    # label parallels on right and top
    # meridians on bottom and left
    parallels = np.arange(lat_min,lat_max,10.)
    # labels = [left,right,top,bottom]
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
    meridians = np.arange(lon_min,lon_max,10.)
    meri = m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)

    cb = plt.colorbar(cs, orientation='vertical')
    cb.set_label('Day Numer', fontsize=8)
    cb.ax.tick_params(labelsize=8)

    plt.text( min(lons)+0.8*(max(lons)-min(lons)), min(lats)+0.2*(max(lats)-min(lats)), str(year[i]),fontsize=8,fontweight='bold', ha='left',va='center',color='r', bbox=dict(boxstyle="round", ec=(0.8, 0.61, 0.11), fc=(0.8, 0.61, 0.11),))

    nf = nf + 1

  plt.show()

def plot_annual_minimum(lat_min, lat_max, lon_min, lon_max, lats, lons, variablename, metric_dir, metric_filename):

  # Read in climatology
  my_file = Path(metric_dir+"/"+metric_filename)
  if not my_file.is_file():
    print("\nPlease create the Annual Minimum first !\n")
    return None
  ncin = Dataset(metric_dir+"/"+metric_filename, 'r')
  data = ncin.variables['data_min'][:]
  day = ncin.variables['day_min'][:]
  year = ncin.variables['year'][:]
  lons = ncin.variables['lon'][:]
  lats = ncin.variables['lat'][:]
  ncin.close()

  ny = year.size

  for i in range(0,year.size):

    if(i%5 == 0):
      # setup plot
      fig = plt.figure(figsize=(3 * 3, 3 * 4))
      fig.suptitle(variablename + ' Annual Minimum and Its Day Number', fontsize=14, fontweight='bold')
      plt.subplots_adjust(left=0.06, right=0.92, top=0.94, bottom=0.05, hspace=.1, wspace=.4)
      nf = 1

    plt.subplot(5,2,nf)

    m = Basemap(projection='cyl', llcrnrlon=lon_min, llcrnrlat=lat_min, urcrnrlon=lon_max, urcrnrlat=lat_max)
    x, y = m(*np.meshgrid(lons, lats))
    clevs = np.linspace(np.nanmin(data[:,:,i].flatten()), np.nanmax(data[:,:,i].flatten()), 21)
    cs=m.contourf(x, y, data[:,:,i].squeeze(), clevs, cmap=plt.cm.RdBu_r)
    m.drawcoastlines()
    m.fillcontinents(color='#000000',lake_color='#99ffff')

    # draw parallels and meridians.
    # label parallels on right and top
    # meridians on bottom and left
    parallels = np.arange(lat_min,lat_max,10.)
    # labels = [left,right,top,bottom]
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
    meridians = np.arange(lon_min,lon_max,10.)
    meri = m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)

    cb = plt.colorbar(cs, orientation='vertical')
    cb.set_label('Maximum Value', fontsize=8)
    cb.ax.tick_params(labelsize=8)

    plt.text( min(lons)+0.8*(max(lons)-min(lons)), min(lats)+0.2*(max(lats)-min(lats)), str(year[i]),fontsize=8,fontweight='bold', ha='left',va='center',color='r', bbox=dict(boxstyle="round", ec=(0.8, 0.61, 0.11), fc=(0.8, 0.61, 0.11),))

    nf = nf + 1
    plt.subplot(5,2,nf)

    m = Basemap(projection='cyl', llcrnrlon=lon_min, llcrnrlat=lat_min, urcrnrlon=lon_max, urcrnrlat=lat_max)
    x, y = m(*np.meshgrid(lons, lats))
    clevs = np.linspace(np.nanmin(day[:,:,i].flatten()), np.nanmax(day[:,:,i].flatten()), 21)
    cs=m.contourf(x, y, day[:,:,i].squeeze(), clevs, cmap=plt.cm.RdBu_r)
    m.drawcoastlines()
    m.fillcontinents(color='#000000',lake_color='#99ffff')

    # draw parallels and meridians.
    # label parallels on right and top
    # meridians on bottom and left
    parallels = np.arange(lat_min,lat_max,10.)
    # labels = [left,right,top,bottom]
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
    meridians = np.arange(lon_min,lon_max,10.)
    meri = m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)

    cb = plt.colorbar(cs, orientation='vertical')
    cb.set_label('Day Numer', fontsize=8)
    cb.ax.tick_params(labelsize=8)

    plt.text( min(lons)+0.8*(max(lons)-min(lons)), min(lats)+0.2*(max(lats)-min(lats)), str(year[i]),fontsize=8,fontweight='bold', ha='left',va='center',color='r', bbox=dict(boxstyle="round", ec=(0.8, 0.61, 0.11), fc=(0.8, 0.61, 0.11),))

    nf = nf + 1

  plt.show()

def plot_metric_spring_start(lat_min, lat_max, lon_min, lon_max, lats, lons, variablename, metricdir, metricfilename):

  # Read in climatology
  my_file = Path(metricdir+"/"+metricfilename)
  if not my_file.is_file():
    print("\nPlease create the Spring Start first !\n")
    return None
  ncin = Dataset(metricdir+"/"+metricfilename, 'r')
  day1 = ncin.variables['day_spring1'][:]
  day2 = ncin.variables['day_spring2'][:]
  year = ncin.variables['year'][:]
  lons = ncin.variables['lon'][:]
  lats = ncin.variables['lat'][:]
  ncin.close()

  ny = year.size

  npage = int(year.size/5.0)

  for i in range(0,year.size):

    if(i%5 == 0):
      # setup plot
      fig = plt.figure(figsize=(3 * 3, 3 * 4))
      fig.suptitle(variablename + ' Spring Start Day', fontsize=14, fontweight='bold')
      plt.subplots_adjust(left=0.06, right=0.92, top=0.94, bottom=0.05, hspace=.1, wspace=.4)
      nf = 1

    plt.subplot(5,2,nf)

    m = Basemap(projection='cyl', llcrnrlon=lon_min, llcrnrlat=lat_min, urcrnrlon=lon_max, urcrnrlat=lat_max)
    x, y = m(*np.meshgrid(lons, lats))
    #clevs = np.linspace(np.nanmin(day1[:,:,i].flatten()), np.nanmax(day1[:,:,i].flatten()), 21)
    clevs = np.linspace(0, 200, 101)
    cs=m.contourf(x, y, day1[:,:,i].squeeze(), clevs, cmap=plt.cm.RdBu_r)
    m.drawcoastlines()
    m.fillcontinents(color='#000000',lake_color='#99ffff')

    # draw parallels and meridians.
    # label parallels on right and top
    # meridians on bottom and left
    parallels = np.arange(lat_min,lat_max,10.)
    # labels = [left,right,top,bottom]
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
    meridians = np.arange(lon_min,lon_max,10.)
    meri = m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)

    cb = m.colorbar(cs, 'right', size='4%', pad='0.%')
    cb.set_label('Day Numer', fontsize=6)
    cb.ax.tick_params(labelsize=6)

    plt.text( min(lons)+0.8*(max(lons)-min(lons)), min(lats)+0.2*(max(lats)-min(lats)), str(year[i]),fontsize=8,fontweight='bold', ha='left',va='center',color='r', bbox=dict(boxstyle="round", ec=(0.8, 0.61, 0.11), fc=(0.8, 0.61, 0.11),))

    nf = nf + 1
    plt.subplot(5,2,nf)

    m = Basemap(projection='cyl', llcrnrlon=lon_min, llcrnrlat=lat_min, urcrnrlon=lon_max, urcrnrlat=lat_max)
    x, y = m(*np.meshgrid(lons, lats))
    #clevs = np.linspace(np.nanmin(day1[:,:,i].flatten()), np.nanmax(day1[:,:,i].flatten()), 21)
    clevs = np.linspace(0, 200, 101)
    cs=m.contourf(x, y, day2[:,:,i].squeeze(), clevs, cmap=plt.cm.RdBu_r)
    m.drawcoastlines()
    m.fillcontinents(color='#000000',lake_color='#99ffff')

    # draw parallels and meridians.
    # label parallels on right and top
    # meridians on bottom and left
    parallels = np.arange(lat_min,lat_max,10.)
    # labels = [left,right,top,bottom]
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
    meridians = np.arange(lon_min,lon_max,10.)
    meri = m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)

    cb = m.colorbar(cs, 'right', size='4%', pad='0.%')
    cb.set_label('Day Numer', fontsize=6)
    cb.ax.tick_params(labelsize=6)

    plt.text( min(lons)+0.8*(max(lons)-min(lons)), min(lats)+0.2*(max(lats)-min(lats)), str(year[i]),fontsize=8,fontweight='bold', ha='left',va='center',color='r', bbox=dict(boxstyle="round", ec=(0.8, 0.61, 0.11), fc=(0.8, 0.61, 0.11),))

    nf = nf + 1

    if nf == 10:
      plt.show()

def plot_metric_spring_start_trend(lat_min, lat_max, lon_min, lon_max, lats, lons, variablename, metricdir, metricfilename):

  # Read in metric file
  my_file = Path(metricdir+"/"+metricfilename)
  if not my_file.is_file():
    print("\nPlease create the Spring Start first !\n")
    return None
  ncin = Dataset(metricdir+"/"+metricfilename, 'r')
  day1 = ncin.variables['day_spring_trend1'][:]
  day2 = ncin.variables['day_spring_trend2'][:]
  year = ncin.variables['year'][:]
  lons = ncin.variables['lon'][:]
  lats = ncin.variables['lat'][:]
  ncin.close()

  ny = year.size

  # setup plot
  fig = plt.figure(figsize=(6, 6))
  fig.suptitle(variablename + ' Spring Start Day Trend', fontsize=14, fontweight='bold')
  plt.subplots_adjust(left=0.06, right=0.92, top=0.94, bottom=0.05, hspace=.1, wspace=.4)
  
  plt.subplot(2,1,1)

  m = Basemap(projection='cyl', resolution='h', llcrnrlon=lon_min, llcrnrlat=lat_min, urcrnrlon=lon_max, urcrnrlat=lat_max)
  x, y = m(*np.meshgrid(lons, lats))
  #clevs = np.linspace(np.nanmin(day1[:,:].flatten()), np.nanmax(day1[:,:].flatten()), 21)
  clevs = np.linspace(-1.0, 1.0, 21)
  cs=m.contourf(x, y, day1, clevs, cmap=plt.cm.RdBu_r)
  m.drawcoastlines()
  m.fillcontinents(color='#000000',lake_color='#99ffff')

  # draw parallels and meridians.
  # label parallels on right and top
  # meridians on bottom and left
  parallels = np.arange(lat_min,lat_max,10.)
  # labels = [left,right,top,bottom]
  m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
  meridians = np.arange(lon_min,lon_max,10.)
  meri = m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)

  cb = m.colorbar(cs, 'right', size='4%', pad='0.%')
  cb.set_label('Day Numer / Year', fontsize=8)
  cb.ax.tick_params(labelsize=8)

  plt.subplot(2,1,2)

  m = Basemap(projection='cyl', resolution='h', llcrnrlon=lon_min, llcrnrlat=lat_min, urcrnrlon=lon_max, urcrnrlat=lat_max)
  x, y = m(*np.meshgrid(lons, lats))
  #clevs = np.linspace(np.nanmin(day1[:,:].flatten()), np.nanmax(day1[:,:].flatten()), 21)
  clevs = np.linspace(-1.0, 1.0, 21)
  cs=m.contourf(x, y, day2, clevs, cmap=plt.cm.RdBu_r)
  m.drawcoastlines()
  m.fillcontinents(color='#000000',lake_color='#99ffff')

  # draw parallels and meridians.
  # label parallels on right and top
  # meridians on bottom and left
  parallels = np.arange(lat_min,lat_max,10.)
  # labels = [left,right,top,bottom]
  m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
  meridians = np.arange(lon_min,lon_max,10.)
  meri = m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)

  cb = m.colorbar(cs, 'right', size='4%', pad='0.%')
  cb.set_label('Day Numer / Year', fontsize=8)
  cb.ax.tick_params(labelsize=8)

  plt.show()

def plot_summer_start_end(lat_min, lat_max, lon_min, lon_max, lats, lons, variablename,  metricdir, metricfilename):

  # Read in metric file
  my_file = Path(metricdir+"/"+metricfilename)
  if not my_file.is_file():
    print("\nPlease create the Summer Start-End first !\n")
    return None
  ncin = Dataset(metricdir+"/"+metricfilename, 'r')
  day_start = ncin.variables['day_summer_start'][:]
  day_end = ncin.variables['day_summer_end'][:]
  year = ncin.variables['year'][:]
  lons = ncin.variables['lon'][:]
  lats = ncin.variables['lat'][:]
  ncin.close()

  ny = year.size
  npage = int(year.size/5.0)

  for i in range(0,year.size):

    if(i%5 == 0):
      # setup plot
      fig = plt.figure(figsize=(3 * 3, 3 * 4))
      fig.suptitle(variablename + ' Summer Start/End Day', fontsize=14, fontweight='bold')
      plt.subplots_adjust(left=0.06, right=0.92, top=0.94, bottom=0.05, hspace=.1, wspace=.4)
      nf = 1

    plt.subplot(5,2,nf)

    m = Basemap(projection='cyl', llcrnrlon=lon_min, llcrnrlat=lat_min, urcrnrlon=lon_max, urcrnrlat=lat_max)
    x, y = m(*np.meshgrid(lons, lats))
    atemp = day_start[:,:,i].flatten()
    atemp = atemp[atemp > 0] 
    #clevs = np.linspace(np.nanmin(day_start[:,:,i].flatten()), np.nanmax(day_start[:,:,i].flatten()), 21)
    if(atemp.size > 0):
      clevs = np.linspace(np.nanmin(atemp), np.nanmax(atemp), 21)
      cs=m.contourf(x, y, day_start[:,:,i].squeeze(), clevs, cmap=plt.cm.RdBu_r)
    m.drawcoastlines()
    m.fillcontinents(color='#000000',lake_color='#99ffff')

    # draw parallels and meridians.
    # label parallels on right and top
    # meridians on bottom and left
    parallels = np.arange(lat_min,lat_max,10.)
    # labels = [left,right,top,bottom]
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
    meridians = np.arange(lon_min,lon_max,10.)
    meri = m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)

    if(atemp.size > 0):
      cb = m.colorbar(cs, 'right', size='4%', pad='0.%')
      cb.set_label('Summber Start Day', fontsize=6)
      cb.ax.tick_params(labelsize=6)

    plt.text( min(lons)+0.8*(max(lons)-min(lons)), min(lats)+0.2*(max(lats)-min(lats)), str(year[i]),fontsize=8,fontweight='bold', ha='left',va='center',color='r', bbox=dict(boxstyle="round", ec=(0.8, 0.61, 0.11), fc=(0.8, 0.61, 0.11),))

    nf = nf + 1
    plt.subplot(5,2,nf)

    m = Basemap(projection='cyl', llcrnrlon=lon_min, llcrnrlat=lat_min, urcrnrlon=lon_max, urcrnrlat=lat_max)
    x, y = m(*np.meshgrid(lons, lats))
    atemp = day_end[:,:,i].flatten()
    atemp = atemp[atemp > 0] 
    #clevs = np.linspace(np.nanmin(day_end[:,:,i].flatten()), np.nanmax(day_end[:,:,i].flatten()), 21)
    if(atemp.size > 0):
      clevs = np.linspace(np.nanmin(atemp), np.nanmax(atemp), 21)
      cs=m.contourf(x, y, day_end[:,:,i].squeeze(), clevs, cmap=plt.cm.RdBu_r)
    m.drawcoastlines()
    m.fillcontinents(color='#000000',lake_color='#99ffff')

    # draw parallels and meridians.
    # label parallels on right and top
    # meridians on bottom and left
    parallels = np.arange(lat_min,lat_max,10.)
    # labels = [left,right,top,bottom]
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
    meridians = np.arange(lon_min,lon_max,10.)
    meri = m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)

    if(atemp.size > 0):
      cb = m.colorbar(cs, 'right', size='4%', pad='0.%')
      cb.set_label('Summber End Day', fontsize=6)
      cb.ax.tick_params(labelsize=8)

    plt.text( min(lons)+0.8*(max(lons)-min(lons)), min(lats)+0.2*(max(lats)-min(lats)), str(year[i]),fontsize=8,fontweight='bold', ha='left',va='center',color='r', bbox=dict(boxstyle="round", ec=(0.8, 0.61, 0.11), fc=(0.8, 0.61, 0.11),))

    nf = nf + 1

  #plt.show()

def plot_summer_start_end_trend(lat_min, lat_max, lon_min, lon_max, lats, lons, variablename,  metricdir, metricfilename):

  # Read in metric file
  my_file = Path(metricdir+"/"+metricfilename)
  if not my_file.is_file():
    print("\nPlease create the Summer Start-End first !\n")
    return None
  ncin = Dataset(metricdir+"/"+metricfilename, 'r')
  day_start_trend = ncin.variables['day_summer_start_trend'][:]
  day_end_trend = ncin.variables['day_summer_end_trend'][:]
  year = ncin.variables['year'][:]
  lons = ncin.variables['lon'][:]
  lats = ncin.variables['lat'][:]
  ncin.close()

  # setup plot
  fig = plt.figure(figsize=(6, 6))
  fig.suptitle(variablename + ' Summer Start/End Day Trend', fontsize=14, fontweight='bold')
  plt.subplots_adjust(left=0.06, right=0.92, top=0.94, bottom=0.05, hspace=.1, wspace=.4)
  
  plt.subplot(2,1,1)

  m = Basemap(projection='cyl', resolution='h', llcrnrlon=lon_min, llcrnrlat=lat_min, urcrnrlon=lon_max, urcrnrlat=lat_max)
  x, y = m(*np.meshgrid(lons, lats))
  #clevs = np.linspace(np.nanmin(day1[:,:].flatten()), np.nanmax(day1[:,:].flatten()), 21)
  clevs = np.linspace(-1.0, 0.0, 21)
  cs=m.contourf(x, y, day_start_trend, clevs, cmap=plt.cm.RdBu_r)
  m.drawcoastlines()
  m.fillcontinents(color='#000000',lake_color='#99ffff')

  # draw parallels and meridians.
  # label parallels on right and top
  # meridians on bottom and left
  parallels = np.arange(lat_min,lat_max,10.)
  # labels = [left,right,top,bottom]
  m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
  meridians = np.arange(lon_min,lon_max,10.)
  meri = m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)

  cb = m.colorbar(cs, 'right', size='4%', pad='0.%')
  cb.set_label('Summer Start Day / Year', fontsize=8)
  cb.ax.tick_params(labelsize=8)

  plt.subplot(2,1,2)

  m = Basemap(projection='cyl', resolution='h', llcrnrlon=lon_min, llcrnrlat=lat_min, urcrnrlon=lon_max, urcrnrlat=lat_max)
  x, y = m(*np.meshgrid(lons, lats))
  #clevs = np.linspace(np.nanmin(day1[:,:].flatten()), np.nanmax(day1[:,:].flatten()), 21)
  clevs = np.linspace(.0, 2.0, 21)
  cs=m.contourf(x, y, day_end_trend, clevs, cmap=plt.cm.RdBu_r)
  m.drawcoastlines()
  m.fillcontinents(color='#000000',lake_color='#99ffff')

  # draw parallels and meridians.
  # label parallels on right and top
  # meridians on bottom and left
  parallels = np.arange(lat_min,lat_max,10.)
  # labels = [left,right,top,bottom]
  m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
  meridians = np.arange(lon_min,lon_max,10.)
  meri = m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)

  cb = m.colorbar(cs, 'right', size='4%', pad='0.%')
  cb.set_label('Summer End Day / Year', fontsize=8)
  cb.ax.tick_params(labelsize=8)

  plt.show()

def plot_monthly_trend(lat_min, lat_max, lon_min, lon_max, lats, lons, variablename, clim_filename):

  # Read in climatology
  my_file = Path(clim_filename)
  if not my_file.is_file():
    print("\nPlease create the Climatology first !\n")
    return None
  ncin = Dataset(clim_filename, 'r')
  data_rate = ncin.variables['data_rate'][:]
  lons = ncin.variables['lon'][:]
  lats = ncin.variables['lat'][:]
  ncin.close()

  ### setup plot
  fig = plt.figure(figsize=(3 * 3, 3 * 4))
  fig.suptitle(variablename + ' Monthly Trend', fontsize=16, fontweight='bold')
  plt.subplots_adjust(bottom=0.05, left=.04, right=.99, top=.95, hspace=.02)
  
  for i in range(1,13):
    plt.subplot(4,3,i)

    m = Basemap(projection='cyl', llcrnrlon=lon_min, llcrnrlat=lat_min, urcrnrlon=lon_max, urcrnrlat=lat_max)
    x, y = m(*np.meshgrid(lons, lats))
    #clevs = np.linspace(data_rate[i-1,:,:].min(), data_rate[i-1,:,:].max(), 21)
    clevs = np.linspace(-0.2, 0.2, 21)
    cs=m.contourf(x, y, data_rate[i-1,:,:].squeeze(), clevs, cmap=plt.cm.RdBu_r)
    m.drawcoastlines()
    m.fillcontinents(color='#000000',lake_color='#99ffff')

    # draw parallels and meridians.
    # label parallels on right and top
    # meridians on bottom and left
    parallels = np.arange(lat_min,lat_max,10.)
    # labels = [left,right,top,bottom]
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=8)
    meridians = np.arange(lon_min,lon_max,10.)
    meri = m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=8)

    cb = plt.colorbar(cs, orientation='horizontal')
    cb.set_label('degree/decade', fontsize=8)
    cb.ax.tick_params(labelsize=8)

    if i == 1:
      plt.title('January', fontsize=12)
    elif i == 2:
      plt.title('February', fontsize=12)
    elif i == 3:
      plt.title('March', fontsize=12)
    elif i == 4:
      plt.title('April', fontsize=12)
    elif i == 5:
      plt.title('May', fontsize=12)
    elif i == 6:
      plt.title('June', fontsize=12)
    elif i == 7:
      plt.title('July', fontsize=12)
    elif i == 8:
      plt.title('August', fontsize=12)
    elif i == 9:
      plt.title('September', fontsize=12)
    elif i == 10:
      plt.title('October', fontsize=12)
    elif i == 11:
      plt.title('November', fontsize=12)
    else:
      plt.title('December', fontsize=12)

  #plt.tight_layout()

  plt.show()
 
def plot_line_box_parameter(pname, dataset_name, startY, endY, lat_boxsize, lon_boxsize, lat_min, lat_max, lon_min, lon_max, lats, lons, variablename, error_variable, clim_filename, dirname, output_dir):

  print("\nPlease wait, reading files ...\n")

  ### get box info
  nx, ny, box_lat_start, box_lat_end, box_lon_start, box_lon_end = data_info.getboxinfo(lat_boxsize, lon_boxsize, lat_min, lat_max, lon_min, lon_max, lats, lons)

  ### Read in data
  for ii in range(0, nx):
    ### setup plot
    fig = plt.figure(ii+1, figsize=(0.5* (endY-startY), 1.5 * ny))
    fig.suptitle(pname+' in Each Grid Box ('+str(lat_boxsize)+'x'+str(lon_boxsize)+' degree)', fontsize=16, fontweight='bold')
    plt.subplots_adjust(bottom=0.05, left=.06, right=.99, top=.95, hspace=.02)
    for jj in range(0, ny):
      ax=plt.subplot(ny, 1, jj+1)
      yearticks = np.arange(startY, endY+1)
      plt.xlim([startY, endY+1])
      plt.ylabel(pname, fontsize='medium')
      plt.xticks(yearticks, fontsize=8,rotation=90)
      plt.yticks(fontsize=8)
      if jj == ny-1:
        plt.xlabel('Year', fontsize='medium')
      if jj < ny-1:
        ax.set_xticklabels([])
        
      tpos = 0.0
      for i in range(startY, endY+1):
        boxfilename = output_dir + '/'+str(i)+'/'+variablename+'_'+str(i)+'_'+dataset_name+'_'+str(lat_boxsize)+'_degree_box_'+str(jj+1)+'x'+str(ii+1)+'.nc'
        ncin = Dataset(boxfilename, 'r')
        if pname == 'Average':
          data_pname = ncin.variables['sst_box_average'][:]
        elif pname == 'Anomaly':
          data_pname = ncin.variables['sst_box_anomaly'][:]
          plt.plot([startY, endY+1], [0,0], color='red', linewidth=1)
        elif pname == 'STD':
          data_pname = ncin.variables['sst_box_std'][:]
        atime = ncin.variables['time'][:]
        plt.plot(i+atime/float(np.max(atime)), data_pname)
        if np.nanmax(data_pname) == np.nanmax(data_pname):
          tpos = np.nanmax(data_pname)
        
      plt.text(endY, tpos, str(jj+1)+'x'+str(ii+1), fontsize=10,fontweight='bold', ha='center',va='center',color='r')

  plt.show()

def plot_line_box_parameter_merged(pname, dataset_name, startY, endY, lat_boxsize, lon_boxsize, lat_min, lat_max, lon_min, lon_max, lats, lons, variablename, error_variable, clim_filename, dirname, output_dir):

  print("\nPlease wait, reading files ...\n")

  ### get box info
  nx, ny, box_lat_start, box_lat_end, box_lon_start, box_lon_end = data_info.getboxinfo(lat_boxsize, lon_boxsize, lat_min, lat_max, lon_min, lon_max, lats, lons)

  ### Read in data
  for ii in range(0, nx):
    ### setup plot
    fig = plt.figure(ii+1, figsize=(0.5* (endY-startY), 1.5 * ny))
    fig.suptitle(pname+' in Each Grid Box ('+str(lat_boxsize)+'x'+str(lon_boxsize)+' degree)', fontsize=16, fontweight='bold')
    plt.subplots_adjust(bottom=0.05, left=.06, right=.99, top=.95, hspace=.02)
    for jj in range(0, ny):
      ax=plt.subplot(ny, 1, jj+1)
      yearticks = np.arange(startY, endY+1)
      plt.xlim([startY, endY+1])
      plt.ylabel(pname, fontsize='medium')
      plt.xticks(yearticks, fontsize=8,rotation=90)
      plt.yticks(fontsize=8)
      if jj == ny-1:
        plt.xlabel('Year', fontsize='medium')
      if jj < ny-1:
        ax.set_xticklabels([])
        
      tpos = 0.0
      for i in range(startY, endY+1):
        boxfilename = output_dir + '/'+str(i)+'/'+variablename+'_'+str(i)+'_'+dataset_name+'_'+str(lat_boxsize)+'_degree_box_'+str(jj+1)+'x'+str(ii+1)+'.nc'
        ncin = Dataset(boxfilename, 'r')
        if pname == 'Average':
          data_pname = ncin.variables['sst_box_average'][:]
        elif pname == 'Anomaly':
          data_pname = ncin.variables['sst_box_anomaly'][:]
          plt.plot([startY, endY+1], [0,0], color='red', linewidth=1)
        elif pname == 'STD':
          data_pname = ncin.variables['sst_box_std'][:]
        atime = ncin.variables['time'][:]
        plt.plot(i+atime/float(np.max(atime)), data_pname)
        if np.nanmax(data_pname) == np.nanmax(data_pname):
          tpos = np.nanmax(data_pname)
        
      plt.text(endY, tpos, str(jj+1)+'x'+str(ii+1), fontsize=10,fontweight='bold', ha='center',va='center',color='r')

  plt.show()

