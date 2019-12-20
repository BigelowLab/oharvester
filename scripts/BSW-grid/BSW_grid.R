library(ncdf4)
library(raster)

# https://www.ncdc.noaa.gov/data-access/marineocean-data/blended-global/blended-sea-winds


#' Convert bounding box [0,360] longitudes to [-180, 180]
#'
#' Bounding boxes are 4 element vectors of [left, right, bottom, top]
#'
#' @export
#' @param x numeric bounding box vector, no check is done for being withing 0,360 range
#' @return numeric bounding box vector
to180BB <- function(x) {
  x[1:2] <- to180(x[1:2])
  if (identical(x[1], 180)) x[1] <- -180
  x}

#' Convert [-180,180] bounding box longitudes to [0,360]
#'
#' Bounding boxes are 4 element vectors of [left, right, bottom, top]
#'
#' @export
#' @param x numeric bounding box vector, no check is done for being withing 0,360 range
#' @return numeric bounding box vector
to360BB <- function(x) {
  x[1:2] <- to360(x[1:2])
  if (identical(x[1], 360)) x[1] <- 0
  x
}

#' Convert [0,360] longitudes to [-180, 180]
#'
#' @export
#' @param x numeric vector, no check is done for being withing [0, 360] range
#' @return numeric vector
to180 <- function(x) {((x + 180) %% 360) - 180}

#' Convert [-180,180] longitudes to [0, 360]
#'
#' @export
#' @param x numeric vector, no check is done for being within [0,3 60] range
#' @return numeric vector
to360 <- function(x) {x %% 360}



#' Compose a BSW url
#'
#' @param period character one of "6h", "daily" (default), "monthly", "climatology"
#' @param type character, one of "stress", "wind" (default)
#' @param base_url character, the base url
BSW_nc_url <- function(period = c("6h", "daily", "monthly", "climatology")[2],
                       type = c("stress", "wind")[2],
                       base_url = "https://www.ncei.noaa.gov/thredds/dodsC/uv"){

  if (tolower(type[1]) == "stress"){
    fname <- switch(tolower(period[1]),
                    "6h" = "6h_strs_agg/Aggregation_of_6h_Ocean_Wind_Stress_best.ncd",
                    "daily" = "daily_strs_agg/Aggregation_of_Daily_Ocean_Wind_Stress_best.ncd",
                    "monthly" = "monthly_strs_agg/Aggregation_of_Monthly_Ocean_Wind_Stress_best.ncd",
                    stop("not available for stress:", period[1]))
  } else {
    fname <- switch(tolower(period[1]),
                    "6h" = "6h_agg_rt/Preliminary_Aggregation_of_6h_Ocean_Wind_best.ncd",
                    "daily" = "daily_agg_rt/Preliminary_Aggregation_of_Daily_Ocean_Wind_best.ncd",
                    "monthly" = "monthly_agg/Aggregation_of_Monthly_Ocean_Wind_best.ncd",
                    "climatology" = "clm_agg/Latest_Climatology_of_Ocean_Winds_best.ncd")
  }
  file.path(base_url, fname)
}

#' Compute the navigation info for BSW
#'
#' @param x ncdf4 object
#' @param bb numeric, 4 element vector of [west, east, south, north] boundaries,
#'        where west and south are negative.
#' @param dates Date-class, vector of ordered sequence (no gaps) Dates
#' @param type character, the data type
#' @param period character, period of the data
#' @param res numeric, 2 element resolution [res_x,res_y]
#' @return list of navigation elements
BSW_nc_nav <- function(x,
  bb = c(-72, -63, 39, 46),
  dates = as.Date(c("2018-01-01","2018-01-02","2018-01-03","2018-01-04")),
  type = 'wind',
  period = 'daily',
  res = c(1,1)/4
  ){

  period <- period[1]
  if (period == 'climatology'){
    time <- as.POSIXct(format(seq(from = as.Date("2018-01-01"),
                       to = as.Date("2018-12-01"),
                       by = 'month'),
                   "%Y-%m-%d 00:00:00"), tz = "UTC")
    idx = seq_len(12)
  } else if (type == 'monthly') {
    time <- (x$dim$time$vals * 3600) +
      as.POSIXct("1987-07-15 00:00:00",
                 format = "%Y-%m-%d %H:%M:%S",
                 tz = 'UTC')
    idx <- range(findInterval(dates, as.Date(time)))
    idx <- seq(from = idx[1], to = idx[2], by = 1)
    time <- time[idx]
  } else if (period == "daily") {
    time <- (x$dim$time$vals * 3600) +
      as.POSIXct("2011-10-01 09:00:00",
                 format = "%Y-%m-%d %H:%M:%S",
                 tz = 'UTC')
    idx <- range(findInterval(dates, as.Date(time)))
    idx <- seq(from = idx[1], to = idx[2], by = 1)
    time <- time[idx]
  } else if (period == "6h"){
    time <- (x$dim$time$vals * 3600) +
      as.POSIXct("2011-10-01 09:00:00",
                 format = "%Y-%m-%d %H:%M:%S",
                 tz = 'UTC')
    idx <- range(findInterval(dates, as.Date(time)))
    idx <- seq(from = idx[1], to = idx[2], by = 1)
    time <- time[idx]
  } else {
    stop("period not known:", period)
  }


  if (length(res) == 1) res <- c(res,res)
  r2 <- res/2
  if (!inherits(x, 'ncdf4')) x <- ncdf4::nc_open(x)
  lat <- x$dim$lat$vals
  lon <- x$dim$lon$vals
  bb2 <- to360BB(bb + c(-r2[1], r2[1], -r2[2], r2[2]))

  iW <- which.min(abs(lon - bb2[1]))
  iE <- which.min(abs(lon - bb2[2]))
  iS <- which.min(abs(lat - bb2[3]))
  iN <- which.min(abs(lat - bb2[4]))

  list(bb = bb,
       res = res,
       start = c(iW, iS, 1, idx[1]),
       count = c(iE - iW + 1, iN - iS + 1, -1, length(idx)),
       ext = c(lon[iW] - r2[1], lon[iE] + r2[1],
               lat[iS] - r2[2], lat[iN] + r2[2]),
       crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",
       type = type[1],
       period = period,
       time = format_time(time, period))
}

#' Convert POSIXct time to a period-sepcifc format
#'
#' @param x POSIXct, one or more timestamps
#' @param period character, the period to convert time to
#' @return character time formated as requested
format_time <- function(x, period = 'daily'){
  switch(period[1],
         'climatology' = month.abb,
         '6h'          = format(x, "%Y-%m-%dT%H"),
         "daily"       = format(x, "%Y-%m-%d"),
         "monthly"     = format(x, "%Y-%m")
    )
}

#' Fetch a stack of BSW rasters
#'
#' @param x ncdf4 object or a character url
#' @param nav list of navigation elements
#' @return a list of raster stacks (wind and sheer each have 3 components)
BSW_fetch <- function(x, nav){

  if (!inherits(x, "ncdf4")) x <- ncdf4::nc_open(x)
  vars <- switch(nav$type,
    'wind' = c("u", "v", "w"),
    'stress' = c("taux", "tauy", "tau"))
  # not every wind resource has "w"
  vars <- vars[vars %in% names(x$var)]
  mm <- sapply(vars,
    function(v){
      ncdf4::ncvar_get(x, varid = v,
                       start = nav$start,
                       count = nav$count)
    }, simplify = FALSE)

  d <- dim(mm[[1]])
  if (length(d) < 3) n <- 1 else n <- d[3]

  if (n == 1){
    ss <- sapply(names(mm),
      function(v){
        raster::raster(t(mm[[v]]),
                       crs = nav$crs,
                       xmn = nav$ext[1],
                       xmx = nav$ext[2],
                       ymn = nav$ext[3],
                       ymx = nav$ext[4])
      }, simplify = FALSE)
    ss <- suppressWarnings(raster::rotate(raster::flip(raster::stack(ss), "y")))
    ss <- list(ss)
  } else {
    ss <- lapply(seq_len(n),
      function(i){
        ss <- sapply(names(mm),
          function(v){
            raster::raster(t(mm[[v]][,,i]),
                           crs = nav$crs,
                           xmn = nav$ext[1],
                           xmx = nav$ext[2],
                           ymn = nav$ext[3],
                           ymx = nav$ext[4])
          }, simplify = FALSE)
        suppressWarnings(raster::rotate(raster::flip(raster::stack(ss),"y")))
      })
  }
  names(ss) <- nav$time
  ss
}
#' Retrieve gridded Blended Sea Winds raster files
#'
#' @param period character one of "6h", "daily" (default), "monthly", "climatology"
#' @param type character, one of "stress", "wind" (default)
#' @param bb numeric, 4 element vector of [west, east, south, north] boundaries with
#'        longitude in the range [-180,180]
#' @param daterange 2 element character in YYYY-mm-dd format or Date day, inclusive
#'        start and end dates
#' @param outpath character or NA, optional output path to save rasters
#' @param overwrite logical, see \code{raster}
#' @param fmt character either 'raster' (default) or 'GTiff' Note that wind and
#'        stress are multiparted (u, x and possible w, or taux, tauy and possibly
#'        tau). Each period is saved a separate file (one per day, on per 6h, etc),
#'        and each file is a brick (single file with multiple layers).  Geotiff
#'        does not preserve layer names - so if you need to know the names of the
#'        layers use the 'raster' option.
#' @return raster stack
BSW_get_grid <- function(period = c("6h", "daily", "monthly", "climatology")[2],
                         type = "wind",
                         bb = c(-88, -48, 24, 52),
                         daterange = c('2018-01-01', '2018-01-04'),
                         outpath = NA,
                         overwrite = TRUE,
                         fmt = c("raster", "GTiff")[1]){

  if (!inherits(daterange, "Date")) daterange <- as.Date(daterange)
  dates <- switch(tolower(period[1]),
      "6h"    = seq.Date(from = daterange[1], to = daterange[2], by = 'day'),
      "daily" = seq.Date(from = daterange[1], to = daterange[2], by = 'day'),
      "month" = seq.Date(from = daterange[1], to = daterange[2], by = 'month'))

  uri <- BSW_nc_url(type = type[1], period = period[1])
  x <- ncdf4::nc_open(uri)
  nav <- BSW_nc_nav(x,
                    bb = bb,
                    dates = dates,
                    period = period,
                    type = type)
  SS <- BSW_fetch(x, nav)
  ncdf4::nc_close(x)

  if (!is.na(outpath)){
    if (!dir.exists(outpath[1])){
      ok <- dir.create(outpath[1], recursive = TRUE, showWarnings = FALSE)
      if (!ok) stop("unable to create outpath:", outpath[1])
    }
    ext <- switch(tolower(fmt[1]),
                  'gtiff' = "tif",
                  'raster' = "grd")
    ofile <- sprintf("%s_%s.%s", type[1], names(SS), ext)
    for (i in seq_along(SS)){
      SS[[i]] <- writeRaster(SS[[i]], file.path(outpath, ofile[i]),
                             overwrite = overwrite, fmt = fmt[1])
    }
  }
  SS
}


