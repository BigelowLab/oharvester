---
title: "README"
author: "btupper"
date: "12/30/2019"
output: 
  html_document:
    keep_md: yes
---



Harvesting data means storing data, and storing data presents both challenges and opportunities. This tutorial is really a biased sales pitch - to urge you to consider early in your process how you will store downloaded data and how you will access it. We present a number of ideas (none are original, but all are time-and-practice tested.)

### Biased sales pitch (time and practice tested)

  + **Treat directories as databases**

    Leverage directories/subdirectories as if they are databases.  Often datasets are served using a directory-like heirarchy; replicating that structure as appropriate can simplify the downloading process.
    
  + **Choose your granularity**
  
    These data sources are often serving data at a granularity that makes for easy searching and storage. [OBPG](https://oceandata.sci.gsfc.nasa.gov/opendap/) provides a very good example where data is stored by `platform/product/year/monthday` (ala `/MODISA/L3SMI/2019/0101`).  It is worth trying to mimic the source data oragnizational paradigm in your own dataset.  Often you may need to alter the paradigm slightly if you are downloading only a specific region - so something like `platform/product/region/year/monthday` (ala `/MODISA/L3SMI/gulfofmaine/2019/0101`) might suite.
    
  + **Retain original filenames if you can**

    It is tempting to rename files as we download and store them. We have found it is better to retain the original filename, which makes it simpler to retrace to the original source.
    
  + **Use a simple companion database file to track files**

    For each directory-database it is useful to have a system for keeping track of what you have in hand - something that you can access and filter programmatically. For example, if you have a `/MODISA/L3SMI/gulfofmaine` with multiple years of files, create a CSV file (or similar) that keeps track of the files included. We have used [feather](https://CRAN.R-project.org/package=feather) format for very large databases.
    
  + **Design a smart database, but not too smart**
  
    The database should have enough info about your data that you can easily select/filter based upon rules you find helpful.  On the other hand, the database should not know too much as to be unweildly.  The worst thing to do is allow your database to know where the data resides in an absolute way.  It should only be permitted enough info to construct relative data locations.  We use a general rule that only info encoded in filenames goes into the database (hence use the original source data filename where possible.)
    
  + **Use software functions to read and filter the database**

    Once you have the database, use that to determine what is in your file collection rather than repeatedly searching the subdirectories.  You can select certain files based upon criteria stored in the database.
    
  + **Build software functions to read your files using the database as input**

    Wrap your file reading function with one that takes a database argument and a path to read one or more data files. 
    
### An example 

We present a [dummy example](https://github.com/BigelowLab/oharvester/blob/master/tutorials/directory-database/database.R) where we create a subset of the [storm database](https://dplyr.tidyverse.org/reference/storms.html) that spans two years, organized by month.  We store a single day's data in each file which makes for a small example, but it could, and usually does, represent a much larger set of big data files. Here is an abstract mock-up.

```
<root_path>/
  database.txt
  year1/
    month1/
      datafile-1-1-1
      datafile-1-1-2
    month7/
      datafile-1-7-1
      datafile-1-7-2
  year2/
    month5/
      datafile-2-5-1
      datafile-2-5-2
    month10/
      datafile-2-10-1
      datafile-2-10-7  
```

#### Root path, reading and writing

Defining the root path to the dataset is one of the maddening parts of this technique. The question is where to store this?  Hardwire it into the code (expedient, but could be painful if you ever move your database)?  Store it in an environmental variable or in .Rprofile (then it is tied to individual accounts)?  Force the user to define a global variable in each new R session (portable but annoying)?  There is no easy answer to this, so for now we'll avoid this issue by defining a function the returns the root path.  In the example below, we assume that you have set your R working directory to the root dataset path.



```r
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

# Define the dataset root path
# @param ... file path segments to jon with the root path
# @param root the root path - with a default value
# @return a file or directory path
storm_path <- function(..., root = "data"){
  file.path(root, ...)
}

# try it!
storm_path()
```

```
## [1] "data"
```

```r
storm_path("2014", "07", "foobar.txt")
```

```
## [1] "data/2014/07/foobar.txt"
```

Next we need to make wrappers around reading and writing the database file.


```r
# Read a database file
# @param path character, the dataset root path
# @param add_date logical, if TRUE add a date column to the database
# @return table of the database
read_database <- function(path = storm_path(), add_date = TRUE){
  filename <- storm_path("database.csv", root = path)
  stopifnot(file.exists(filename))
  x <- suppressMessages(readr::read_csv(filename,
                                        col_types ='ccc' ))
  if (add_date) x <- x %>%
    dplyr::mutate(date = as.Date(sprintf("%s-%s-%s", year, month, day)))
  x
}

# Write a database file
# @param x table of the database (typically tibble)
# @param path character, the dataset path
# @return the input table
write_database <- function(x, path = storm_root()){
  filename <- file.path(path, "database.csv")
  readr::write_csv(x, filename)
}

# try it
path <- storm_path()
db <- read_database(path)
db
```

```
## # A tibble: 93 x 4
##    year  month day   date      
##    <chr> <chr> <chr> <date>    
##  1 2014  07    01    2014-07-01
##  2 2014  07    02    2014-07-02
##  3 2014  07    03    2014-07-03
##  4 2014  07    04    2014-07-04
##  5 2014  07    05    2014-07-05
##  6 2014  07    21    2014-07-21
##  7 2014  07    22    2014-07-22
##  8 2014  07    23    2014-07-23
##  9 2014  08    23    2014-08-23
## 10 2014  08    24    2014-08-24
## # … with 83 more rows
```


### Reading data

Now we can leverage the tools in the [dplyr](https://dplyr.tidyverse.org/) to filter our full database.  Let's find the data files that fall between `2014-10-10` and `2015-06-30`.  


```r
db2 <- db %>% 
    dplyr::filter(dplyr::between(date, as.Date("2014-10-10"),as.Date("2015-06-30")))
cat("nrow(db2) = ", nrow(db2), "\n")
```

```
## nrow(db2) =  18
```

Next, make a function to accept this 18-record database and the root path to read
the files.


```r
# Read one or more data files - these are tables, so we can optionally 
# row bind them.
#
# @param x a database of file descriptors
# @param path character, the dataset root path
# @param do_bind logical, if TRUE row bind the results
# @return list of table  if do_bind is TRUE then return a table, otherwise
#         return a losyt of tables
read_storms <- function(x, path = storm_path(), do_bind = TRUE){
  opath <- file.path(path, x$year, x$month)
  ofile <- sprintf("%s-%s-%s-storms.csv",x$year, x$month, x$day)
  files <- file.path(opath, ofile)
  stopifnot(all(file.exists(files)))
  x <- lapply(files,
              function(f){
                suppressMessages(readr::read_csv(f, col_types = 'cccc?????????'))
              })
  if (do_bind) x <- dplyr::bind_rows(x)
  x
}

x <- read_storms(db2, path)
x
```

```
## # A tibble: 78 x 13
##    name  year  month day    hour   lat  long status category  wind pressure
##    <chr> <chr> <chr> <chr> <dbl> <dbl> <dbl> <chr>     <dbl> <dbl>    <dbl>
##  1 Fay   2014  10    11        6  26   -65   tropi…        0    55      994
##  2 Fay   2014  10    11       12  27.2 -65.3 tropi…        0    60      990
##  3 Fay   2014  10    11       18  28.7 -65.4 tropi…        0    60      988
##  4 Fay   2014  10    12        0  30.2 -65.2 tropi…        0    60      987
##  5 Fay   2014  10    12        6  31.7 -64.9 hurri…        1    65      985
##  6 Fay   2014  10    12        8  32.3 -64.7 hurri…        1    70      984
##  7 Fay   2014  10    12       12  33.1 -63.9 hurri…        1    70      983
##  8 Fay   2014  10    12       18  33.6 -61.9 tropi…        0    60      986
##  9 Gonz… 2014  10    12        0  16.4 -55.9 tropi…       -1    30     1010
## 10 Gonz… 2014  10    12        6  16.4 -56.9 tropi…       -1    30     1008
## # … with 68 more rows, and 2 more variables: ts_diameter <dbl>,
## #   hu_diameter <dbl>
```

### Changing the database

As files are added/removed to the dataset collection you should add/remove entries in your database.  These you can acheive using the [dplyr](https://dplyr.tidyverse.org/) tools and your read/write_database functions.  Once in a while it is useful to simply rebuild a database from scratch.  For that you'll want a function.


```r
# Build a database and possibly save
#
# @param path character, the path to the dataset
# @param save_db logical, if TRUE save the database possibly overwriting an existing one
# @return database table
build_database <- function(path = storm_path(), save_db = FALSE){
  ff <- list.files(path, 
                   pattern = "^.*-storms\\.csv$",
                   recursive = TRUE)
  if (length(ff) == 0) stop("no files found")
  ss <- strsplit(basename(ff), "-", fixed = TRUE)
  x <- dplyr::tibble(
    year = sapply(ss, "[[", 1),
    month = sapply(ss, "[[", 2),
    day = sapply(ss, "[[", 3))
  if (save_db) x <- write_database(x, path=path)
  x
}

# read and verify that it matches what we had originally (without date)
db3 <- build_database(path, save_db = FALSE)
identical(dim(db %>% dplyr::select(-date)), dim(db3))
```

```
## [1] TRUE
```
