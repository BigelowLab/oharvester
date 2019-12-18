library(ncdf4)
library(dplyr)
library(readr)

#' Construct URLs for ICOADS observations
#'
#' @param x Date-class one or more dates
#' @param base_url charcater, the base URL to the resources
#' @return vector of URLs
ICOADS_nc_url <- function(x,
                          base_url = "https://data.nodc.noaa.gov/thredds/dodsC/icoads"){

  if (!inherits(x, "Date")) x <- as.Date(x)
  # 2010s/2010s/ICOADS_R3.0.1_2018-01.nc
  Ys  <- paste0(substring(format(x, "%Y"), 1,3), "0s")
  Ym <- format(x, "%Y-%m")
  vdates <- c("R3.0.0" = as.Date("1662-10-15"), "R3.0.1" = as.Date("2015-01-01"))
  iv <- findInterval(x,vdates)
  d <- sprintf("ICOADS_%s_%s.nc", names(vdates)[iv], Ym)
  file.path(base_url, Ys, Ys, d)
}


#' Fetch the data from one opendap resource
#'
#' @param x character or ncdf4 object, if character then it should be a URL
#' @param varnames character one or more variable names to add tot he core variables
#' @param corenames character one or more core variables to retrieve (time and location)
#' @param origin POSIXct the origin of the timestamp
#' @param outpath character or NA, if character write this table to the specified output path
#' @return table (tibble) of the requested variables
ICOADS_fetch <- function(x,
  varnames = c("SST", "WP", "WH","AT", "RH", "SLP"),
  corenames = c("time", "lat", "lon"),
  origin = as.POSIXct("1662-10-15 12:00:00",
                      format = "%Y-%m-%d %H:%M:%S",
                      tz = 'UTC'),
  outpath = NA){

  if (!inherits(x, 'ncdf4')) x <- ncdf4::nc_open(x)
  vn <- names(x$var)
  ix <- varnames %in% vn
  if (!all(ix)){
    stop("one or more varnames not known:", paste(varnames[!ix], sep = " "))
  }
  ix <- corenames %in% vn
  if (!all(ix)){
    stop("one or more corenames not known:", paste(corenames[!ix], sep = " "))
  }

  xx <- sapply(c(corenames, varnames),
               function(vn){
                 ncdf4::ncvar_get(x, vn)
               }, simplify = FALSE)
  if ("time" %in% names(xx)){
    xx$time <- as.POSIXct(xx$time*86400, origin = origin, tz = 'UTC')
  }
  filename <- x$filename
  ncdf4::nc_close(x)
  x <- dplyr::as_tibble(xx)
  if (!is.na(outpath)){
    if (!dir.exists(outpath)){
      ok <- dir.create(outpath, recursive = TRUE, showWarnings = FALSE)
      if (!ok) stop("unable to create output directory:", outpath)
    }
    ofile <- gsub(".nc", ".csv.gz", basename(filename), fixed = TRUE)
    x <- readr::write_csv(x, file.path(outpath, ofile))
  }
  x
}

#' Retrieve ICOADS observation date in month-sized tables
#'
#' @param dates charcater or Date-class, one or more dates ("YYYY-mm-dd" format)
#'   which indicate the months to retrieve.  Note that they need not be sequential nor
#'   contiguous.  Also note that if you request multiple dates within one month that
#'   only on table is returned - thus requesting three dates like this
#'   c("2018-01-01", "2018-01-17", "2018-02-04") will yield just two tables.
#' @param varnames character one of more variables to retrieve
#' @param outpath NA or character, an optional output path.  Leave as NA to skip
#'   writing to files to disk
#' @param ... further arguments for ICOADS_fetch
#' @return list of tables (tibbles) named by "YYYY-mm"
ICOADS_get_obs <- function(dates = c("1818-01-01", "1918-01-01", "2018-01-01"),
                           varnames = c("SST", "AT"),
                           outpath = NA,
                           ...){

  if (!inherits(dates, 'Date')) dates <- as.Date(dates)
  if (!is.na(outpath)) {
    if (!dir.exists(outpath)){
      ok <- dir.create(outpath, recursive = TRUE, showWarnings = FALSE)
      if (!ok) stop("unable to create output directory:", outpath)
    }
  }
  icoads_url <- unique(ICOADS_nc_url(dates))
  icoads <- as.list(icoads_url)
  icoads[[1]] <- ncdf4::nc_open(icoads[[1]][1])
  vn <- names(icoads[[1]]$var)
  ix <- varnames %in% vn
  if (!all(ix)){
    stop("one or more varnames not known:", paste(varnames[!ix], collapse = " "))
  }
  xx <- lapply(icoads, ICOADS_fetch, varnames = varnames, outpath = outpath, ...)
  names(xx) <- substring(basename(icoads_url),15, 21)
  xx
}
