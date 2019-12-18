# OISST Optimally Interpolated Sea Surface Temperature gridded time series

This tool will harvest a gridded sea surface temperature time series from [OISST](https://www.ncdc.noaa.gov/oissthttps://www.ncdc.noaa.gov/oisst) and (if you want) write the results to files. 

Keep in mind that you have options to select such as `type` ("AVHRR" or "AVHRR-AMSR"), `daterange`, `varname` ("sst", "anom", "ice", etc), as well as others. Please consult the script file for details.

1) Have the file `OISST-grid.R` in your working directory
2) If you don't have the [ncdf4](https://CRAN.R-project.org/package=ncdf4) and [raster](https://CRAN.R-project.org/package=raster) packages already installed, install them:

```
install.packages(c('ncdf4', 'raster'))
```

3) Source the file:

```
source('OISST-grid.R')
```

4) Extract the data, optionally save the data to a location of your choice, and plot.

```
sst <- OISST_get_grid(bb =  c(-72, -63, 39, 46),
                      daterange = c('2018-01-01', '2018-01-04'),
                      outpath = "./my_oi_sst")
plot(sst)
```
![](OISST-grid.png)

This will get you the daily sea surface temperature 'sst' data within the specified latitude and longitude bounding box, and between the specified start and end dates. The data will be written to the directory "./my_oi_sst" (or whatever you name it).

*Note: if you don't want to write files, but rather only load the data into R, then omit the `outfile` argument.  See the script for other available arguments.*

---

Developer notes:

 - This is quick-and-dirty
 
 - There is limited error checking, so use at your own risk
