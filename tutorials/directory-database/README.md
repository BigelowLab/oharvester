---
title: "README"
author: "btupper"
date: "12/30/2019"
output: 
  html_document:
    keep_md: yes
---



Harvesting data means storing data, and storing data presents both challenges and opportunities. This tutorial is really a biased sales pitch - to consider how to store downloaded data and how to access it. We present a number of ideas (none are original, but all are time-and-practice tested.)

### Biased sales pitch (time and practice tested)

  + Directories as databases

    Leverage directories as if they are databases.  Often datasets are served using a directory-like heirarchy.  These are often serving data at a granularity that makes for easy searching and storage. [OBPG](https://oceandata.sci.gsfc.nasa.gov/opendap/) for provides a very good example where data is stored by platform/product/year/monthday (ala `/MODISA/L3SMI/2019/0101`).  It is worth trying to mimic the source data oragnizational paradigm in your own dataset.  Often you may need to alter the paradigm slightly if you are downloading only a specific region - so something like platform/product/region/year/monthday (ala `/MODISA/L3SMI/gulfofmaine/2019/0101`) might suite.
    
  + Retain original filenames

    It is tempting to rename files as we download and store them. We have found it is better to retain the original filename, which makes it simpler to retrace to the original source.
    
  + Use a simple database file to track files

    For each directory-database it is useful to have a system for keeping track of what you have in hand - something that you can access and filter programmatically. For example, if you have a `/MODISA/L3SMI/gulfofmaine` with multiple years of files, create a CSV file (or similar) that keeps track of the files included.  
    
  + Use software functions to read and filter the database

    Once you have the database, use that to determine what is in your file collection rathaer than repeatedly searching the directories.  You can select certain files based upon criteria stored in the database.
    
  + Build software functions to read your files using the database as input

    Wrap your file reading function with one that takes a filtered database to read one or more file. 
    
### An example 

We present a dummy example where we create a mock dataset (from [real data](https://dplyr.tidyverse.org/reference/storms.html)) that spans two years, organized by month.  The data stored in each file is tiny for this exampe, but it could represent a much larger set of data files.

