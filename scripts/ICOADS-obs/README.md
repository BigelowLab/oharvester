# ICOADS observational data

This tool will harvest tables of sobservations entered into [ICOADS](https://icoads.noaa.gov/) and optionally write the results to files.

Keep in mind that you have options to specify arguments such as `dates`, `varnames` ("SST", "AT", "WH", etc), as well as others. Please consult the script file for details.

1) Have the file `ICOADS-obs.R` in your working directory
2) Be sure to install the following packages (if you haven't already) 

   + [ncdf4](https://CRAN.R-project.org/package=ncdf4)
   
   + [dplyr](https://CRAN.R-project.org/package=dplyr) 
   
   + [readr](https://CRAN.R-project.org/package=readr)

```
install.packages(c('ncdf4', 'dplyr', 'readr'))
```

3) Source the file:

```
source('ICOADS-obs.R')
```

4) Extract the data, optionally save the data to a location of your choice. Note that if you opt to save files they will be saved by month as gzipped CSV.

```
xx <- ICOADS_get_obs(dates = c("1818-01-01", "1918-01-01", "2018-01-01"),
                     varnames = c("SST", "AT"),
                     outpath = "./my_icoads_obs")
names(xx)
# [1] "1818-01" "1918-01" "2018-01"

sapply(xx, nrow)
# 1818-01 1918-01 2018-01 
#     871   30283 3661770 

xx[['1818-0']]
# # A tibble: 871 x 5
#    time                     lat   lon   SST    AT
#    <chr>                  <dbl> <dbl> <dbl> <dbl>
#  1 1818-01-01T00:00:00Z  36.7   283.   23.1  13.3
#  2 1818-01-01T00:00:00Z  35.1   352.   NA    NA  
#  3 1818-01-01T00:00:00Z   8.87  331    NA    NA  
#  4 1818-01-01T00:00:00Z  -8.23   81.2  NA    NA  
#  5 1818-01-01T00:59:59Z -39.0    73.8  NA    NA  
#  6 1818-01-01T03:59:59Z  -5.63  116.   NA    27.8
#  7 1818-01-01T05:00:00Z  -0.570 106.   NA    25.6
#  8 1818-01-01T05:00:00Z  -6.57  105.   NA    28.3
#  9 1818-01-01T05:00:00Z -38.9    74.5  NA    NA  
# 10 1818-01-01T06:00:00Z  17.3    83.2  NA    26.9
# # … with 861 more rows
```

*Note: if you don't want to write files, but rather only load the data into R, then omit the `outfile` argument.  See the script for other available arguments.*

---

Developer notes:

 - This is quick-and-dirty

 - There is limited error checking, so use at your own risk
 
 - This is abusive of the concept of the THREDDS system as we have hardwired the change from ICOADS v3.0.0 to v3.0.1 by date.
