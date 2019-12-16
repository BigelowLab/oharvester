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
#from mpl_toolkits.basemap import Basemap
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
from matplotlib.ticker import AutoMinorLocator
import datetime
import calendar
from sklearn import linear_model, datasets
from pathlib import Path

import data_info
import subset_dataset
import phenologyalg
import phenologyplt

################################################################################################   

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

def standalone_main():

  # get command line options:
  options=parseoptions()

  config_filename = options.config

  #*** get input file config file ***********************

  config_input = configparser.ConfigParser()
  config_input.read(config_filename)

  #*** Start Time ***************************************

  start = time.time()
  
  #*** Prepare Data *************************************
  sst_dataset_names = config_input.get('SST', 'Dataset_Name')
  sst_dataset_names = sst_dataset_names.split(" ")
  chlor_dataset_names = config_input.get('Chlorophyll', 'Dataset_Name')
  chlor_dataset_names = chlor_dataset_names.split(" ")

  sst_display_str = ""
  for i in range(0, len(sst_dataset_names)):
    sst_display_str = sst_display_str+str(i+1)+": "+sst_dataset_names[i] + " "
  #sst_display_str = sst_display_str+"-1: Quit\n"
  print("\nSea Surface Temperature Datasets: ")
  print(sst_display_str)

  chlor_display_str = ""
  for i in range(0, len(chlor_dataset_names)):
    chlor_display_str = chlor_display_str+str(len(sst_dataset_names)+i+1)+": "+chlor_dataset_names[i] + " "
  chlor_display_str = chlor_display_str+"\n\n-1: Quit\n"
  print("\nChlorophyll Datasets: ")
  print(chlor_display_str)

  userselection = 0

  while userselection == 0:
    aselection = input('Please make your selection: ')
    try:
      # Try to convert the user input to an integer
      userselection = int(aselection)
      # Catch the exception if the input was not a number
    except ValueError:
      userselection = 0

  if(userselection <= len(sst_dataset_names)):
    data_type = "SST"
  else:
    data_type = "Chlor"

  if(userselection != -1):
    if(data_type == "SST"):
      dataset_name = sst_dataset_names[userselection-1]
    elif(data_type == "Chlor"):
      dataset_name = chlor_dataset_names[userselection-1-len(sst_dataset_names)]
  else:
    sys.exit(1)

  if not config_input.has_section(dataset_name): 
    print("*** There is no config setup for "+dataset_name+", please try again later. ***")
    sys.exit(1)

  if(data_type == "SST"):
    #shortname = podaac_dataset_shorname(config_input.get('DEFAULT', 'Dataset_Name'))
    shortname = data_info.podaac_dataset_shorname(dataset_name)

  # Ask for user action
  userselection = 0
  while userselection != -1:

   print("\n*** 1. Subset and Download Dataset from PODAAC (this may take long time depending on the network speed)")
   print("*** 2. Climatology")
   print("*** 3. Climatology in the Box Calculation")
   print("*** 4. Parameter Calculations (Average, Anomaly, STD ...)")
   print("*** 5. Phenology Matrics: Spring Start Calculation")
   print("*** 6. Phenology Matrics: Summer Start/End Calculation")
   print("*** 7. Plotting")
   print("*** -1. Quit\n")
  
   aselection = input('Please make your selection: ')
   try:
      # Try to convert the user input to an integer
      userselection = int(aselection)
      # Catch the exception if the input was not a number
   except ValueError:
      userselection = -1

   if(userselection == 1 and data_type == "SST"):
    subset_dataset.download_data(config_input.get(dataset_name, 'Start_Year'),
                                 config_input.get(dataset_name, 'Start_Month'),
                                 config_input.get(dataset_name, 'Start_Day'),
                                 config_input.get(dataset_name, 'End_Year'),
                                 config_input.get(dataset_name, 'End_Month'),
                                 config_input.get(dataset_name, 'End_Day'),
                                 float(config_input.get('DEFAULT', 'Lat_Min')),
                                 float(config_input.get('DEFAULT', 'Lat_Max')),
                                 float(config_input.get('DEFAULT', 'Lon_Min')),
                                 float(config_input.get('DEFAULT', 'Lon_Max')),
                                 int(config_input.get(dataset_name, 'Grid_Skip')),
                                 shortname, config_input.get(dataset_name, 'Root_Name'))
    userselection = 0
   elif(userselection == 1 and data_type == "Chlor" and dataset_name == "SeaWiFS"):
    subset_dataset.download_data_seawifs(config_input.get(dataset_name, 'Start_Year'),
                                 config_input.get(dataset_name, 'Start_Month'),
                                 config_input.get(dataset_name, 'Start_Day'),
                                 config_input.get(dataset_name, 'End_Year'),
                                 config_input.get(dataset_name, 'End_Month'),
                                 config_input.get(dataset_name, 'End_Day'),
                                 float(config_input.get('DEFAULT', 'Lat_Min')),
                                 float(config_input.get('DEFAULT', 'Lat_Max')),
                                 float(config_input.get('DEFAULT', 'Lon_Min')),
                                 float(config_input.get('DEFAULT', 'Lon_Max')),
                                 int(config_input.get(dataset_name, 'Grid_Skip')),
                                 config_input.get(dataset_name, 'Root_Name'))
    userselection = 0
   elif(userselection == 1 and data_type == "Chlor" and dataset_name == "MODIS_Aqua_Chlor"):
    subset_dataset.download_data_modis(config_input.get(dataset_name, 'Start_Year'),
                                 config_input.get(dataset_name, 'Start_Month'),
                                 config_input.get(dataset_name, 'Start_Day'),
                                 config_input.get(dataset_name, 'End_Year'),
                                 config_input.get(dataset_name, 'End_Month'),
                                 config_input.get(dataset_name, 'End_Day'),
                                 float(config_input.get('DEFAULT', 'Lat_Min')),
                                 float(config_input.get('DEFAULT', 'Lat_Max')),
                                 float(config_input.get('DEFAULT', 'Lon_Min')),
                                 float(config_input.get('DEFAULT', 'Lon_Max')),
                                 int(config_input.get(dataset_name, 'Grid_Skip')),
                                 config_input.get(dataset_name, 'Root_Name'))
    userselection = 0
   elif(userselection == 2 and data_type == "SST"):
    lats, lons = data_info.latloninfo(config_input.get(dataset_name, 'Root_Name'), shortname)
    phenologyalg.climatology(lats, lons, int(config_input.get(dataset_name, 'Start_Year')), 
                int(config_input.get(dataset_name, 'End_Year')), 
                float(config_input.get('DEFAULT', 'Thresh_Spring1')), 
                float(config_input.get('DEFAULT', 'Thresh_Spring2')), 
                float(config_input.get('DEFAULT', 'Thresh_Offset_Summer')), 
                config_input.get(dataset_name, 'Parameter_Name'), 
		config_input.get(dataset_name, 'Root_Name'),
		config_input.get(dataset_name, 'Climatology_File'), 
		config_input.get(dataset_name, 'Root_Name'))
    userselection = 0
   elif(userselection == 2 and data_type == "Chlor"):
    lats, lons = data_info.latloninfo_seawifs(config_input.get(dataset_name, 'Root_Name'), int(config_input.get(dataset_name, 'Start_Year')))
    phenologyalg.climatology_chlor(lats, lons, 
                int(config_input.get(dataset_name, 'Start_Year')), 
                int(config_input.get(dataset_name, 'End_Year')), 
                float(config_input.get('DEFAULT', 'Lat_Min')),
                float(config_input.get('DEFAULT', 'Lat_Max')),
                float(config_input.get('DEFAULT', 'Lon_Min')),
                float(config_input.get('DEFAULT', 'Lon_Max')),
                #float(config_input.get('DEFAULT', 'Thresh_Spring')), 
                #float(config_input.get('DEFAULT', 'Thresh_Offset_Summer')), 
                config_input.get(dataset_name, 'Parameter_Name'), 
		config_input.get(dataset_name, 'Root_Name'),
		config_input.get(dataset_name, 'Climatology_File'), 
		config_input.get(dataset_name, 'Root_Name'))
    userselection = 0
   elif(userselection == 3 and data_type == "SST"):
    lats, lons = data_info.latloninfo(config_input.get(dataset_name, 'Root_Name'), shortname)
    phenologyalg.climatologybox(float(config_input.get('DEFAULT', 'Lat_Boxsize')),
                   float(config_input.get('DEFAULT', 'Lon_Boxsize')),
                   float(config_input.get('DEFAULT', 'Lat_Min')),
                   float(config_input.get('DEFAULT', 'Lat_Max')),
                   float(config_input.get('DEFAULT', 'Lon_Min')),
                   float(config_input.get('DEFAULT', 'Lon_Max')),
                   lats, lons, 
                   config_input.get(dataset_name, 'Parameter_Name'), 
		   config_input.get(dataset_name, 'Root_Name')+"/"+config_input.get(dataset_name, 'Climatology_File'), 
		   config_input.get(dataset_name, 'Root_Name')+"/"+config_input.get(dataset_name, 'Climatology_File_Box') )
    userselection = 0
   elif(userselection == 3 and data_type == "Chlor"):
    lats, lons = data_info.latloninfo_seawifs(config_input.get(dataset_name, 'Root_Name'), int(config_input.get(dataset_name, 'Start_Year')))
    phenologyalg.climatologybox_chlor(float(config_input.get('DEFAULT', 'Lat_Boxsize')),
                   float(config_input.get('DEFAULT', 'Lon_Boxsize')),
                   float(config_input.get('DEFAULT', 'Lat_Min')),
                   float(config_input.get('DEFAULT', 'Lat_Max')),
                   float(config_input.get('DEFAULT', 'Lon_Min')),
                   float(config_input.get('DEFAULT', 'Lon_Max')),
                   lats, lons, 
                   config_input.get(dataset_name, 'Parameter_Name'), 
		   config_input.get(dataset_name, 'Root_Name')+"/"+config_input.get(dataset_name, 'Climatology_File'), 
		   config_input.get(dataset_name, 'Root_Name')+"/"+config_input.get(dataset_name, 'Climatology_File_Box') )
    userselection = 0
   elif(userselection == 4):
    lats, lons = data_info.latloninfo(config_input.get(dataset_name, 'Root_Name'), shortname)
    phenologyalg.indexProcessing(dataset_name,
                    int(config_input.get(dataset_name, 'Start_Year')),
                    int(config_input.get(dataset_name, 'End_Year')),
                    float(config_input.get('DEFAULT', 'Lat_Boxsize')),
                    float(config_input.get('DEFAULT', 'Lon_Boxsize')),
                    float(config_input.get('DEFAULT', 'Lat_Min')),
                    float(config_input.get('DEFAULT', 'Lat_Max')),
                    float(config_input.get('DEFAULT', 'Lon_Min')),
                    float(config_input.get('DEFAULT', 'Lon_Max')),
                    lats, lons, 
                    config_input.get(dataset_name, 'Parameter_Name'), 
                    config_input.get(dataset_name, 'Parameter_Error_Name'), 
		    config_input.get(dataset_name, 'Root_Name')+"/"+config_input.get(dataset_name, 'Climatology_File_Box'), 
		    config_input.get(dataset_name, 'Root_Name'),
		    config_input.get(dataset_name, 'Root_Name') + "/" + config_input.get(dataset_name, 'Output_Dir') )
    userselection = 0
   elif(userselection == 5):
    lats, lons = data_info.latloninfo(config_input.get(dataset_name, 'Root_Name'), shortname)
    phenologyalg.metric_spring_start(lats, lons, int(config_input.get(dataset_name, 'Start_Year')), 
                int(config_input.get(dataset_name, 'End_Year')), 
                float(config_input.get('DEFAULT', 'Thresh_Spring1')), 
                float(config_input.get('DEFAULT', 'Thresh_Spring2')), 
                config_input.get(dataset_name, 'Parameter_Name'), 
		config_input.get(dataset_name, 'Root_Name')+"/"+config_input.get(dataset_name, 'Metric_Dir'), 
		config_input.get(dataset_name, 'Metric_Spring_File'), 
		config_input.get(dataset_name, 'Root_Name')+"/"+config_input.get(dataset_name, 'Climatology_File'), 
		config_input.get(dataset_name, 'Root_Name'))
    userselection = 0
   elif(userselection == 6):
    lats, lons = data_info.latloninfo(config_input.get(dataset_name, 'Root_Name'), shortname)
    phenologyalg.metric_summer_start_end(lats, lons, int(config_input.get(dataset_name, 'Start_Year')), 
                int(config_input.get(dataset_name, 'End_Year')), 
                float(config_input.get('DEFAULT', 'Thresh_Offset_Summer')), 
                config_input.get(dataset_name, 'Parameter_Name'), 
		config_input.get(dataset_name, 'Root_Name')+"/"+config_input.get(dataset_name, 'Metric_Dir'), 
		config_input.get(dataset_name, 'Metric_Summer_File'), 
		config_input.get(dataset_name, 'Root_Name')+"/"+config_input.get(dataset_name, 'Climatology_File'), 
		config_input.get(dataset_name, 'Root_Name'))
    userselection = 0
   elif(userselection == 7):
     aselection = 0 
     while aselection != -1:
       print("\n*** 1. Plot Grid Box")
       print("*** 2. Contour Plot Climatology")
       print("*** 3. Line Plot Climatology")
       print("*** 4. Plot Monthly Trend")
       print("*** 5. Plot Annual Maximum and its day number")
       print("*** 6. Plot Annual Minimum and its day number")
       print("*** 7. Plot Spring Start Day & Trend")
       print("*** 8. Plot Summer Start/End Day & Trend")
       print("*** 9. Plot Average")
       print("*** 10. Plot Anomaly")
       print("*** 11. Plot STD")
       print("*** 12. Plot Average (Datasets Comparision)")
       print("*** -1. Quit\n")

       aselection = input('Please make your plot selection: ')
       try:
           # Try to convert the user input to an integer
           aselection = int(aselection)
           # Catch the exception if the input was not a number
       except ValueError:
           aselection = -1

       if aselection == 1 and data_type == "SST":
         lats, lons = data_info.latloninfo(config_input.get(dataset_name, 'Root_Name'), shortname)
         phenologyplt.plot_grid_box(int(config_input.get(dataset_name, 'Start_Year')),
                    int(config_input.get(dataset_name, 'End_Year')),
                    float(config_input.get('DEFAULT', 'Lat_Boxsize')),
                    float(config_input.get('DEFAULT', 'Lon_Boxsize')),
                    float(config_input.get('DEFAULT', 'Lat_Min')),
                    float(config_input.get('DEFAULT', 'Lat_Max')),
                    float(config_input.get('DEFAULT', 'Lon_Min')),
                    float(config_input.get('DEFAULT', 'Lon_Max')),
                    lats, lons,
                    config_input.get(dataset_name, 'Parameter_Name'),
                    config_input.get(dataset_name, 'Parameter_Error_Name'),
                    config_input.get(dataset_name, 'Root_Name')+"/"+config_input.get(dataset_name, 'Climatology_File_Box'),
                    config_input.get(dataset_name, 'Root_Name'))
       elif aselection == 1 and data_type == "Chlor":
         lats, lons = data_info.latloninfo_seawifs(config_input.get(dataset_name, 'Root_Name'), int(config_input.get(dataset_name, 'Start_Year')))
         phenologyplt.plot_grid_box(int(config_input.get(dataset_name, 'Start_Year')),
                    int(config_input.get(dataset_name, 'End_Year')),
                    float(config_input.get('DEFAULT', 'Lat_Boxsize')),
                    float(config_input.get('DEFAULT', 'Lon_Boxsize')),
                    float(config_input.get('DEFAULT', 'Lat_Min')),
                    float(config_input.get('DEFAULT', 'Lat_Max')),
                    float(config_input.get('DEFAULT', 'Lon_Min')),
                    float(config_input.get('DEFAULT', 'Lon_Max')),
                    lats, lons,
                    config_input.get(dataset_name, 'Parameter_Name'),
                    config_input.get(dataset_name, 'Parameter_Error_Name'),
                    config_input.get(dataset_name, 'Root_Name')+"/"+config_input.get(dataset_name, 'Climatology_File_Box'),
                    config_input.get(dataset_name, 'Root_Name'))
       elif aselection == 2 and data_type == "SST":
         lats, lons = data_info.latloninfo(config_input.get(dataset_name, 'Root_Name'), shortname)
         phenologyplt.plot_clim_box_contour(float(config_input.get('DEFAULT', 'Lat_Boxsize')),
                    float(config_input.get('DEFAULT', 'Lon_Boxsize')),
                    float(config_input.get('DEFAULT', 'Lat_Min')),
                    float(config_input.get('DEFAULT', 'Lat_Max')),
                    float(config_input.get('DEFAULT', 'Lon_Min')),
                    float(config_input.get('DEFAULT', 'Lon_Max')),
                    lats, lons,
                    config_input.get(dataset_name, 'Parameter_Name'),
                    config_input.get(dataset_name, 'Parameter_Error_Name'),
                    config_input.get(dataset_name, 'Root_Name')+"/"+config_input.get(dataset_name, 'Climatology_File'))
       elif aselection == 2 and data_type == "Chlor":
         lats, lons = data_info.latloninfo_seawifs(config_input.get(dataset_name, 'Root_Name'), int(config_input.get(dataset_name, 'Start_Year')))
         phenologyplt.plot_clim_box_contour_chlor(float(config_input.get('DEFAULT', 'Lat_Boxsize')),
                    float(config_input.get('DEFAULT', 'Lon_Boxsize')),
                    float(config_input.get('DEFAULT', 'Lat_Min')),
                    float(config_input.get('DEFAULT', 'Lat_Max')),
                    float(config_input.get('DEFAULT', 'Lon_Min')),
                    float(config_input.get('DEFAULT', 'Lon_Max')),
                    lats, lons,
                    config_input.get(dataset_name, 'Parameter_Name'),
                    config_input.get(dataset_name, 'Parameter_Error_Name'),
                    config_input.get(dataset_name, 'Root_Name')+"/"+config_input.get(dataset_name, 'Climatology_File'))
       elif aselection == 3 and data_type == "SST":
         lats, lons = data_info.latloninfo(config_input.get(dataset_name, 'Root_Name'), shortname)
         phenologyplt.plot_line_box_climatology(float(config_input.get('DEFAULT', 'Lat_Boxsize')),
                    float(config_input.get('DEFAULT', 'Lon_Boxsize')),
                    float(config_input.get('DEFAULT', 'Lat_Min')),
                    float(config_input.get('DEFAULT', 'Lat_Max')),
                    float(config_input.get('DEFAULT', 'Lon_Min')),
                    float(config_input.get('DEFAULT', 'Lon_Max')),
                    lats, lons, 
                    config_input.get(dataset_name, 'Parameter_Name'), 
		    config_input.get(dataset_name, 'Root_Name')+"/"+config_input.get(dataset_name, 'Climatology_File'), 
		    config_input.get(dataset_name, 'Root_Name')+"/"+config_input.get(dataset_name, 'Climatology_File_Box') )
       elif aselection == 3 and data_type == "Chlor":
         lats, lons = data_info.latloninfo_seawifs(config_input.get(dataset_name, 'Root_Name'), int(config_input.get(dataset_name, 'Start_Year')))
         phenologyplt.plot_line_box_climatology_chlor(float(config_input.get('DEFAULT', 'Lat_Boxsize')),
                    float(config_input.get('DEFAULT', 'Lon_Boxsize')),
                    float(config_input.get('DEFAULT', 'Lat_Min')),
                    float(config_input.get('DEFAULT', 'Lat_Max')),
                    float(config_input.get('DEFAULT', 'Lon_Min')),
                    float(config_input.get('DEFAULT', 'Lon_Max')),
                    lats, lons, 
                    config_input.get(dataset_name, 'Parameter_Name'), 
		    config_input.get(dataset_name, 'Root_Name')+"/"+config_input.get(dataset_name, 'Climatology_File'), 
		    config_input.get(dataset_name, 'Root_Name')+"/"+config_input.get(dataset_name, 'Climatology_File_Box') )
       elif aselection == 4:
         lats, lons = data_info.latloninfo(config_input.get(dataset_name, 'Root_Name'), shortname)
         phenologyplt.plot_monthly_trend( float(config_input.get('DEFAULT', 'Lat_Min')),
                    float(config_input.get('DEFAULT', 'Lat_Max')),
                    float(config_input.get('DEFAULT', 'Lon_Min')),
                    float(config_input.get('DEFAULT', 'Lon_Max')),
                    lats, lons, 
                    config_input.get(dataset_name, 'Parameter_Name'), 
		    config_input.get(dataset_name, 'Root_Name')+"/"+config_input.get(dataset_name, 'Climatology_File') )
       elif aselection == 5:
         lats, lons = data_info.latloninfo(config_input.get(dataset_name,'Root_Name'), shortname)
         phenologyplt.plot_annual_maximum( float(config_input.get('DEFAULT', 'Lat_Min')),
                    float(config_input.get('DEFAULT', 'Lat_Max')),
                    float(config_input.get('DEFAULT', 'Lon_Min')),
                    float(config_input.get('DEFAULT', 'Lon_Max')),
                    lats, lons, 
                    config_input.get(dataset_name,'Parameter_Name'), 
		    config_input.get(dataset_name,'Root_Name')+"/"+config_input.get(dataset_name, 'Metric_Dir'),
		    config_input.get(dataset_name, 'Metric_Spring_File')) 
       elif aselection == 6:
         lats, lons = data_info.latloninfo(config_input.get(dataset_name,'Root_Name'), shortname)
         phenologyplt.plot_annual_minimum( float(config_input.get('DEFAULT', 'Lat_Min')),
                    float(config_input.get('DEFAULT', 'Lat_Max')),
                    float(config_input.get('DEFAULT', 'Lon_Min')),
                    float(config_input.get('DEFAULT', 'Lon_Max')),
                    lats, lons, 
                    config_input.get(dataset_name,'Parameter_Name'), 
		    config_input.get(dataset_name,'Root_Name')+"/"+config_input.get(dataset_name, 'Metric_Dir'),
		    config_input.get(dataset_name, 'Metric_Spring_File')) 
       elif aselection == 7:
         lats, lons = data_info.latloninfo(config_input.get(dataset_name,'Root_Name'), shortname)
         phenologyplt.plot_metric_spring_start( float(config_input.get('DEFAULT', 'Lat_Min')),
                    float(config_input.get('DEFAULT', 'Lat_Max')),
                    float(config_input.get('DEFAULT', 'Lon_Min')),
                    float(config_input.get('DEFAULT', 'Lon_Max')),
                    lats, lons, 
                    config_input.get(dataset_name,'Parameter_Name'), 
		    config_input.get(dataset_name,'Root_Name')+"/"+config_input.get(dataset_name, 'Metric_Dir'), 
		    config_input.get(dataset_name, 'Metric_Spring_File')) 
         phenologyplt.plot_metric_spring_start_trend( float(config_input.get('DEFAULT', 'Lat_Min')),
                    float(config_input.get('DEFAULT', 'Lat_Max')),
                    float(config_input.get('DEFAULT', 'Lon_Min')),
                    float(config_input.get('DEFAULT', 'Lon_Max')),
                    lats, lons, 
                    config_input.get(dataset_name,'Parameter_Name'), 
		    config_input.get(dataset_name,'Root_Name')+"/"+config_input.get(dataset_name, 'Metric_Dir'), 
		    config_input.get(dataset_name, 'Metric_Spring_File')) 
       elif aselection == 8:
         lats, lons = data_info.latloninfo(config_input.get(dataset_name,'Root_Name'), shortname)
         phenologyplt.plot_summer_start_end( float(config_input.get('DEFAULT', 'Lat_Min')),
                    float(config_input.get('DEFAULT', 'Lat_Max')),
                    float(config_input.get('DEFAULT', 'Lon_Min')),
                    float(config_input.get('DEFAULT', 'Lon_Max')),
                    lats, lons, 
                    config_input.get(dataset_name,'Parameter_Name'), 
		    config_input.get(dataset_name,'Root_Name')+"/"+config_input.get(dataset_name, 'Metric_Dir'), 
		    config_input.get(dataset_name, 'Metric_Summer_File'))
         phenologyplt.plot_summer_start_end_trend( float(config_input.get('DEFAULT', 'Lat_Min')),
                    float(config_input.get('DEFAULT', 'Lat_Max')),
                    float(config_input.get('DEFAULT', 'Lon_Min')),
                    float(config_input.get('DEFAULT', 'Lon_Max')),
                    lats, lons, 
                    config_input.get(dataset_name,'Parameter_Name'), 
		    config_input.get(dataset_name,'Root_Name')+"/"+config_input.get(dataset_name, 'Metric_Dir'), 
		    config_input.get(dataset_name, 'Metric_Summer_File'))
       elif aselection == 9:
         lats, lons = data_info.latloninfo(config_input.get(dataset_name,'Root_Name'), shortname)
         phenologyplt.plot_line_box_parameter('Average', dataset_name, int(config_input.get(dataset_name,'Start_Year')),
                    int(config_input.get(dataset_name,'End_Year')),
                    float(config_input.get('DEFAULT', 'Lat_Boxsize')),
                    float(config_input.get('DEFAULT', 'Lon_Boxsize')),
                    float(config_input.get('DEFAULT', 'Lat_Min')),
                    float(config_input.get('DEFAULT', 'Lat_Max')),
                    float(config_input.get('DEFAULT', 'Lon_Min')),
                    float(config_input.get('DEFAULT', 'Lon_Max')),
                    lats, lons, 
                    config_input.get(dataset_name,'Parameter_Name'), 
                    config_input.get(dataset_name,'Parameter_Error_Name'), 
		    config_input.get(dataset_name,'Root_Name')+"/"+config_input.get(dataset_name, 'Climatology_File_Box'), 
		    config_input.get(dataset_name,'Root_Name'),
		    config_input.get(dataset_name,'Root_Name') + "/" + config_input.get(dataset_name, 'Output_Dir') )
       elif aselection == 10:
         lats, lons = data_info.latloninfo(config_input.get(dataset_name,'Root_Name'), shortname)
         phenologyplt.plot_line_box_parameter('Anomaly', dataset_name, int(config_input.get(dataset_name,'Start_Year')),
                    int(config_input.get(dataset_name,'End_Year')),
                    float(config_input.get('DEFAULT', 'Lat_Boxsize')),
                    float(config_input.get('DEFAULT', 'Lon_Boxsize')),
                    float(config_input.get('DEFAULT', 'Lat_Min')),
                    float(config_input.get('DEFAULT', 'Lat_Max')),
                    float(config_input.get('DEFAULT', 'Lon_Min')),
                    float(config_input.get('DEFAULT', 'Lon_Max')),
                    lats, lons, 
                    config_input.get(dataset_name,'Parameter_Name'), 
                    config_input.get(dataset_name,'Parameter_Error_Name'), 
		    config_input.get(dataset_name,'Root_Name')+"/"+config_input.get(dataset_name, 'Climatology_File_Box'), 
		    config_input.get(dataset_name,'Root_Name'),
		    config_input.get(dataset_name,'Root_Name') + "/" + config_input.get(dataset_name, 'Output_Dir') )
       elif aselection == 11:
         lats, lons = data_info.latloninfo(config_input.get(dataset_name,'Root_Name'), shortname)
         phenologyplt.plot_line_box_parameter('STD', dataset_name, int(config_input.get(dataset_name,'Start_Year')),
                    int(config_input.get(dataset_name,'End_Year')),
                    float(config_input.get('DEFAULT', 'Lat_Boxsize')),
                    float(config_input.get('DEFAULT', 'Lon_Boxsize')),
                    float(config_input.get('DEFAULT', 'Lat_Min')),
                    float(config_input.get('DEFAULT', 'Lat_Max')),
                    float(config_input.get('DEFAULT', 'Lon_Min')),
                    float(config_input.get('DEFAULT', 'Lon_Max')),
                    lats, lons, 
                    config_input.get(dataset_name,'Parameter_Name'), 
                    config_input.get(dataset_name,'Parameter_Error_Name'), 
		    config_input.get(dataset_name,'Root_Name')+"/"+config_input.get(dataset_name, 'Climatology_File_Box'), 
		    config_input.get(dataset_name,'Root_Name'),
		    config_input.get(dataset_name,'Root_Name') + "/" + config_input.get(dataset_name, 'Output_Dir') )
       elif aselection == 12:
         stemp = config_input.get('MERGE_DATASET', 'Dataset_Name')
         dataset_names = stemp.split()
         stemp = config_input.get('MERGE_DATASET', 'Root_Name')
         root_names = stemp.split()

         for i in range(0, len(dataset_names)):
           if not os.path.exists(root_names[i]):
             print("\n>>> Please create dataset '"+dataset_names[i]+"' output first ! <<<\n")
             sys.exit(1)

         lats, lons = data_info.latloninfo(config_input.get(dataset_name,'Root_Name'), shortname)
         phenologyplt.plot_line_box_parameter_merged('STD', dataset_name, int(config_input.get(dataset_name,'Start_Year')),
                    int(config_input.get(dataset_name,'End_Year')),
                    float(config_input.get('DEFAULT', 'Lat_Boxsize')),
                    float(config_input.get('DEFAULT', 'Lon_Boxsize')),
                    float(config_input.get('DEFAULT', 'Lat_Min')),
                    float(config_input.get('DEFAULT', 'Lat_Max')),
                    float(config_input.get('DEFAULT', 'Lon_Min')),
                    float(config_input.get('DEFAULT', 'Lon_Max')),
                    lats, lons, 
                    config_input.get(dataset_name,'Parameter_Name'), 
                    config_input.get(dataset_name,'Parameter_Error_Name'), 
		    config_input.get(dataset_name,'Root_Name')+"/"+config_input.get(dataset_name, 'Climatology_File_Box'), 
		    config_input.get(dataset_name,'Root_Name'),
		    config_input.get(dataset_name,'Root_Name') + "/" + config_input.get(dataset_name, 'Output_Dir') )
       elif aselection == -1:
         sys.exit(1)
   else:
    sys.exit(1)


  #*** End Time *****************************************
  end = time.time()

  #*** Total Time ***************************************
  print (' Time spend = ' + str(end - start) + ' seconds')

################################################################################################   

if __name__ == "__main__":
        standalone_main()

