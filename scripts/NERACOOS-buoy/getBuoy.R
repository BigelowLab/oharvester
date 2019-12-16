library(thredds)
library(ncdf4)

buoylookup <- function(buoyname='B01',baseURL='http://www.neracoos.org/thredds/dodsC/UMO/DSG/SOS')
{
  switch(buoyname[1],
         'A01'=file.path(baseURL,'A01/CTD20m/HistoricRealtime'),
         'B01'=file.path(baseURL,'B01/CTD20m/HistoricRealtime')
  )
}

getbuoy <- function(buoyname='A01',daterange=c('2018-01-01','2018-01-10'))
{
  baseurl <- buoylookup(buoyname)
  NC <- ncdf4::nc_open(baseurl)
  dates <- NC$dim$time$vals*86400 + as.POSIXct('1858-11-17 00:00:00',format='%Y-%m-%d %H:%M:%S',tz='UTC')
  daterange <- as.POSIXct(daterange,tz='UTC')
  idx <- findInterval(dates,daterange)
  len <- sapply(NC$var, '[[','size')
  len <- len[len==length(idx)]
  idxstart <- which(idx==1)[1]
  idxlength <- sum(idx==1)
  XX <- sapply(names(len),
               function(nm){
                 ncdf4::ncvar_get(NC,varid=nm, start = idxstart, count = idxlength)
               }, simplify=F )
  
}