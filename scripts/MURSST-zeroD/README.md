# MUR Sea Surface Temperature time series

This tool will harvest a sea surface temperature time series from MURsst 
and (if you want) write it to a file 

1) Have the file `MURsst-zeroD.R` in your working directory
2) If you don't have the `rerddap` package already installed, install it:

```install.packages('rerddap')```

3) Source the file:

```source('MURsst-zeroD.R')```

4) Extract the data at the location of your choice:

```sst <- MURtimeseries(latitude = 43.862131, longitude = -69.573120, starttime = '2018-01-01', 
endtime = '2018-01-31', outfile = 'bigelow.csv')```

This will get you the daily sea surface temperature data at the specified latitude and longitude, 
and between the specified start and end dates. 
The data will be written to the file "bigelow.csv" (or whatever you name it).

*Note: if you don't want to write a file, but rather load the data into R, then omit the `outfile` argument*


