library(ncdf4)
library(raster)

#https://podaac-opendap.jpl.nasa.gov/opendap/allData/ghrsst/data/GDS2/L4/GLOB/JPL/MUR/v4.1/2018/001/20180101090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc

#' Construct opendap urls for provided dates
#'
#' @param dates Date-class, one or more dates
#' @param base_url character, the base URL
#' @param time character, we assume the time is always 9am
#' @return character URLs
MUR_nc_url <- function(dates,
  base_url = file.path("https://podaac-opendap.jpl.nasa.gov/opendap",
                       "allData/ghrsst/data/GDS2/L4/GLOB/JPL/MUR/v4.1"),
  time = "090000"){

  # 2018/001/20180101090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc
  y <- format(dates, "%Y")
  j <- format(dates, "%j")
  d <- sprintf("%s%s%s",
               format(dates, "%Y%m%d"),
               time[1],
               "-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc")

  file.path(base_url,y,j,d)

  }

#' Retrieve MUR navigation values (start, count, lons, lats)
#'
#' @param x ncdfd4 object
#' @param bb numeric, 4 element requested bounding box [west, east, south, north]
#' @param res numeric, 2 element resolution [res_x,res_y]
#' @param varname character the name of the variable
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
mur_nc_nav <- function(x,
                    bb = c(-69.6, -69.5, 43.8, 43.9),
                    res = c(0.01, 0.01),
                    varname = "analysed_sst"){
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
       start = c(ix[1], iy[1], 1),
       count = c(ix[2] - ix[1] + 1, iy[2] - iy[1] + 1, -1),
       ext = c(we + (res[1]/2 * c(-1,1)), sn + (res[2]/2 * c(-1,1)) ),
       crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",
       varname = varname)
}

#' Fetch a single MUR raster
#'
#' @param x ncdf4 object or URL
#' @param nav list as returned by \code{MUR_nc_nav}
#' @param outpath characater or NA to prevent file writing
#' @param overwrite logical, see \code{raster}
#' @param fmt character either 'raster' or 'GTiff' (default)
#' @return a raster layer
MUR_fetch <- function(x,
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
    r <- raster::raster(t(m),
                   crs = nav$crs,
                   xmn = nav$ext[1], xmx = nav$ext[2],
                   ymn = nav$ext[3], ymx = nav$ext[4])
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
#' @seealso https://podaac.jpl.nasa.gov/dataset/MUR-JPL-L4-GLOB-v4.1
#' @param latitude numeric, two element vector of [south, north] boundaries, positive north
#' @param longitude numeric, two element vector of [west, east] boundaries, positive east
#' @param starttime character in YYYY-mm-dd format or Date day, inclusive start date
#' @param endtime character in YYYY-mm-dd format or Date day, inclusive end date
#' @param outpath character or NA, optional output path to save rasters in GeoTiff format
#' @param ... other arguments passed to MUR_fetch
#' @return raster stack
MURgrid <- function(latitude =  c(39, 46),
                    longitude = c(-72, -63),
                    startdate = '2018-01-01',
                    enddate = '2018-01-04',
                    varname = 'analysed_sst',
                    outpath = NA,
                    ...){

  if (!inherits(startdate, "Date")) startdate <- as.Date(startdate[1])
  if (!inherits(enddate, "Date")) enddate <- as.Date(enddate[1])
  dates <- seq.Date(from = startdate[1], to = enddate[1], by = 'day')
  mur <- as.list(MUR_nc_url(dates))
  mur[[1]] <- ncdf4::nc_open(mur[[1]])
  nav <- mur_nc_nav(mur[[1]],
                    bb = c( sort(longitude[1:2]) , sort(latitude[1:2]) ),
                    varname = varname)

  S <- raster::stack(lapply(mur, MUR_fetch, nav, outpath = outpath, ...))
  S <- raster::setZ(S, dates)
  names(S) <- format(dates, "%b-%d-%Y")
  S
}
