library("rerddap")

MURtimeseries <- function(latitude = 43.862131, 
                      longitude = -69.573120,
                      starttime = '2018-01-01',
                      endtime = '2018-01-31',
                      outfile = NA)
{
  latitude <- c(latitude,latitude)
  longitude <- c(longitude,longitude)
  timerange <- c(starttime,endtime)
  murSST <- griddap(info('jplMURSST41'),
                    latitude=latitude,
                    longitude=longitude,
                    time=timerange,
                    fields='analysed_sst')
  if (!is.na(outfile)) {
    write.csv(x=murSST$data,file=outfile)
  }
  murSST
}

