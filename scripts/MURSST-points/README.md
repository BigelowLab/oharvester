MURsst-points
================

# MUR Sea Surface Temperature at points

This tool will harvest sea surface temperature at points (know locations
and dates) from
[MURsst](https://podaac.jpl.nasa.gov/dataset/MUR-JPL-L4-GLOB-v4.1?ids=&values=&search=MUR%20v4.1)
and (if you want) write the results to files.

#### Setup

Have the file `MURSST-points.R` in your working directory

#### Install packages (optional)

If you don’t have the [ncdf4](https://CRAN.R-project.org/package=ncdf4),
[dplyr]() and [readr]() packages installed then please install them.

    install.packages(c('ncdf4', 'dplyr', 'readr'))

#### Source

Source the file: (be sure it is in your working directory)

``` r
source('MURsst-points.R')
```

#### Prep point data

Gather you point data into a data frame
([tibbles](https://tibble.tidyverse.org/) are handy). You can do this by
reading values in from a file, such as CSV, or making one on the fly. In
either case, the data frame must have columns for longitude, latitude
and date (in `YYYY-mm-dd` string format or as the R `Date` class). Below
we make on on the fly.

``` r
points <- tribble(
  ~lon, ~lat, ~date,
  -50,  38,   "2017-01-01",
  -50,  38,   "2017-02-01",
  -50,  38,   "2017-03-01",  
  -50,  38,   "2017-04-01",
  -14,  -34,  "2017-01-01",
  -14,  -34,  "2017-02-01",
  -14,  -34,  "2017-03-01",  
  -14,  -34,  "2017-04-01",  
  -80,  -11,  "2017-01-01",
  -80,  -11,  "2017-02-01",
  -80,  -11,  "2017-03-01",  
  -80,  -11,  "2017-04-01"
)
```

#### Extract

Extract the data, save the data to a location of your choice, and plot.

``` r
mur <- MURpoints(points, outfile = "mur_points.csv")
mur
```

    ## # A tibble: 12 x 8
    ##      lon   lat date       status analysed_sst analysis_error  mask
    ##    <dbl> <dbl> <date>     <lgl>         <dbl>          <dbl> <int>
    ##  1   -50    38 2017-01-01 TRUE           293.           0.39     1
    ##  2   -14   -34 2017-01-01 TRUE           292.           0.38     1
    ##  3   -80   -11 2017-01-01 TRUE           297.           0.4      1
    ##  4   -50    38 2017-02-01 TRUE           291.           0.38     1
    ##  5   -14   -34 2017-02-01 TRUE           295.           0.38     1
    ##  6   -80   -11 2017-02-01 TRUE           301.           0.39     1
    ##  7   -50    38 2017-03-01 TRUE           291.           0.39     1
    ##  8   -14   -34 2017-03-01 TRUE           295.           0.38     1
    ##  9   -80   -11 2017-03-01 TRUE           300.           0.37     1
    ## 10   -50    38 2017-04-01 TRUE           291.           0.39     1
    ## 11   -14   -34 2017-04-01 TRUE           295.           0.4      1
    ## 12   -80   -11 2017-04-01 TRUE           300.           0.38     1
    ## # … with 1 more variable: sea_ice_fraction <dbl>

This will get you the daily sea surface temperature data at the
specified places and times. The data will be written to the file you
name (or it won’t if you don’t provide a name.)

*Note: if you don’t want to write files, but rather only load the data
into R, then omit the `outfile` argument. See the script for other
available arguments.*

------------------------------------------------------------------------

Developer notes:

-   This is quick-and-dirty

-   There is limited error checking, so use at your own risk
