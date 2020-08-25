library(ncdf4)
library(raster)

#https://www.ncei.noaa.gov/thredds/dodsC/OisstBase/NetCDF/AVHRR/201801/avhrr-only-v2.20180101.nc

# https://www.ncei.noaa.gov/thredds/dodsC/OisstBase/NetCDF/AVHRR-AMSR/201109/amsr-avhrr-v2.20110901.nc

#' Convert bounding box [0,360] longitudes to [-180, 180]
#'
#' Bounding boxes are 4 element vectors of [left, right, bottom, top]
#'
#' @export
#' @param x numeric bounding box vector, no check is done for being withing 0,360 range
#' @return numeric bounding box vector
to180BB <- function(x) {x[1:2] <- to180(x[1:2]) ; x}

#' Convert [-180,180] bounding box longitudes to [0,360]
#'
#' Bounding boxes are 4 element vectors of [left, right, bottom, top]
#'
#' @export
#' @param x numeric bounding box vector, no check is done for being withing 0,360 range
#' @return numeric bounding box vector
to360BB <- function(x) {x[1:2] <- to360(x[1:2]) ; x}

#' Convert [0,360] longitudes to [-180, 180]
#'
#' @export
#' @param x numeric vector, no check is done for being withing [0, 360] range
#' @return numeric vector
to180 <- function(x) {ix <- x > 180 ; x[ix] <- x[ix]-360; x}

#' Convert [-180,180] longitudes to [0, 360]
#'
#' @export
#' @param x numeric vector, no check is done for being within [0,3 60] range
#' @return numeric vector
to360 <- function(x) {ix <- x < 0 ; x[ix] <- x[ix]+ 360; x}

#' Construct opendap urls for provided dates
#'
#' @param dates Date-class, one or more dates
#' @param type character, either "AVHRR" (default) or "AVHRR+AMSR" (not available for V2.1)
#' @param version character, the OISST version, defaults to "V2.1"
#' @param release character, the OISST version release, defaults to "v02r01"
#' @param base_url character, the base URL
#' @return character URLs
OISST_nc_url <- function(dates,
  type = c("AVHRR", "AVHRR-AMSR")[1],
  version = "V2.1",
  release = "v02r01",
  base_url = "https://www.ncei.noaa.gov/thredds/dodsC/OisstBase/NetCDF"){

  # AVHRR-AMSR/201109/amsr-avhrr-v2.20110901.nc
  # AVHRR/201801/avhrr-only-v2.20180101.nc
  # https://www.ncei.noaa.gov/thredds/dodsC/OisstBase/NetCDF/V2.1/AVHRR/202008/oisst-avhrr-v02r01.20200801.nc
  ym <- format(dates, "%Y%m")
  if (version[1] == 'V2.1'){
    typ <- switch(type[1],
                  'AVHRR' = sprintf("oisst-avhrr-%s.", release),
                  'AVHRR-AMSR' = stop("AVHRR-AMSR not available for V2.1"),
                  stop("type not known:", type[1]))
  } else {
    typ <- switch(type[1],
                  'AVHRR' = "avhrr-only-v2.",
                  'AVHRR-AMSR' = "amsr-avhrr-v2.",
                  stop("type not known:", type[1]))
  }
  ym <- format(dates, "%Y%m")
  d <- paste0(typ, format(dates, "%Y%m%d"), ".nc")
  if (version[1] == 'V2.1'){
    uri <- file.path(base_url, version[1], 'AVHRR', ym, d)
  } else {
    uri <- file.path(base_url, type[1], ym, d)
  }
  uri
}
#' Retrieve OISST navigation values (start, count, lons, lats)
#'
#' @param x ncdfd4 object
#' @param bb numeric, 4 element requested bounding box [west, east, south, north]
#' @param res numeric, 2 element resolution [res_x,res_y]
#' @param varname character the name of the variable.
#'        Options include "ice", "err", "anom", and "sst" (default)
#' @return list with
#' \itemize{
#'   \item{bb the requested bounding box}
#'   \item{res the resolution}
#'   \item{start vector of start values}
#'   \item{count vector of count values}
#'   \item{ext vector of extent (for raster::raster)}
#'   \item{crs character, proj string for raster::raster}
#'   \item{varname character}
#' }
OISST_nc_nav <- function(x,
                    bb = c(-69.6, -69.5, 43.8, 43.9),
                    res = c(0.25, 0.25),
                    varname = "sst"){
  stopifnot(inherits(x, 'ncdf4'))
  if (!(varname[1] %in% names(x$var))) stop("varname not known:", varname[1])
  if (length(res) == 1) res <- c(res[1],res[1])
  ix <- sapply(bb[1:2],
               function(xbb) which.min(abs(x$dim$lon$vals-xbb)))
  we <- x$dim$lon$vals[ix]
  iy <- sapply(bb[3:4],
               function(ybb) which.min(abs(x$dim$lat$vals-ybb)))
  sn <- x$dim$lat$vals[iy]

  list(bb = bb,
       res = res,
       start = c(ix[1], iy[1], 1, 1),
       count = c(ix[2] - ix[1] + 1, iy[2] - iy[1] + 1, -1, -1),
       ext = c(we + (res[1]/2 * c(-1,1)), sn + (res[2]/2 * c(-1,1)) ),
       crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",
       varname = varname)
}

#' Fetch a single OISST raster
#'
#' @param x ncdf4 object or URL
#' @param nav list as returned by \code{OISST_nc_nav}
#' @param outpath characater or NA to prevent file writing
#' @param overwrite logical, see \code{raster}
#' @param fmt character either 'raster' or 'GTiff' (default)
#' @return a raster layer
OISST_fetch <- function(x,
                      nav,
                      outpath = NA,
                      overwrite = TRUE,
                      fmt = "GTiff"){
    if (!inherits(x, "ncdf4")) x <- ncdf4::nc_open(x)
    m <- ncdf4::ncvar_get(x,
                          varid = nav$varname,
                          start = nav$start,
                          count = nav$count)
    ncdf4::nc_close(x)
    bb <- to180BB(nav$ext)
    r <- raster::raster(t(m),
                   crs = nav$crs,
                   xmn = bb[1], xmx = bb[2],
                   ymn = bb[3], ymx = bb[4])
    r <- raster::flip(r,"y")
    if (!is.na(outpath)){
      if (!dir.exists(outpath)) {
        stopifnot(dir.create(outpath, recursive = TRUE, showWarnings = FALSE))
      }
      ext <- c("raster" = ".grd" , "GTiff" = ".tif" )[fmt]
      ofile <- gsub(".nc", ext, basename(x$filename), fixed = TRUE)
      raster::writeRaster(r, file.path(outpath, ofile),
                          overwrite = TRUE,
                          format = fmt)
    }
  r
}

#' Retrieve one or more raster images
#'
#' Note that spatial boundaries are requests.  Actual bounds returned are to the
#' nearest whole pixel.
#'
#' @seealso https://www.ncdc.noaa.gov/oisst
#' @param latitude numeric, two element vector of [south, north] boundaries, positive north
#' @param bb numeric, 4 element vector of [west, east, south, north] boundaries,
#'        where west and south are negative.
#' @param type character, either "AVHRR" (default) or "AVHRR+AMSR"
#' @param varname character the name of the variable.
#'        Options include "ice", "err", "anom", and "sst" (default)
#' @param daterange 2 element character in YYYY-mm-dd format or Date day, inclusive start and end dates
#' @param outpath character or NA, optional output path to save rasters in GeoTiff format
#' @param ... other arguments passed to OISST_fetch.
#' @return raster stack
OISST_get_grid <- function(bb = c(-72, -63, 39, 46),
                    type = c("AVHRR", "AVHRR-AMSR")[1],
                    daterange = c('2018-01-01', '2018-01-04'),
                    varname = 'sst',
                    outpath = NA,
                    ...){
  if (!inherits(daterange, "Date")) daterange <- as.Date(daterange)
  dates <- seq.Date(from = daterange[1], to = daterange[2], by = 'day')
  oi <- as.list(OISST_nc_url(dates))
  oi[[1]] <- ncdf4::nc_open(oi[[1]])
  nav <- OISST_nc_nav(oi[[1]],
                    bb = to360BB(bb),
                    varname = varname)

  S <- raster::stack(lapply(oi, OISST_fetch, nav, outpath = outpath, ...))
  S <- raster::setZ(S, dates)
  names(S) <- format(dates, "%b-%d-%Y")
  S
}
