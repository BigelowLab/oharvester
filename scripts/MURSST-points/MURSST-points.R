library(ncdf4)
library(dplyr)
library(readr)

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


#' Find the index location for a given value to a vector of values. In the case of ties
#' the first closest match is returned.
#' 
#' @param x one or more values
#' @param v vector of values
#' @return index of the closest value in v to each element in x  
closest_index <- function(x, v){
  
  ix <- sapply(x,
               function(x, vec = NULL){
                 which.min(abs(vec-x))
               }, v = v)
  return(ix)
}


#' Fetch MUR data for a set of points on the same date
#' 
#' @param tbl data frame of lon, lat and date
#' @param key data frame, convenience from tidyverse
#' @param vars the name(s) of the variables to retrieve
#' @return the input data frame with added columns 
#' \itemize{
#' \item{status logical, TRUE if the connection was successful}
#' \item{var1 ... depends upon the variables requested}
#' \item{var2 ... etc}
#' }
fetch_MUR <- function(tbl, key, 
                      vars = c("analysed_sst", 
                               "analysis_error", 
                               "mask",
                               "sea_ice_fraction")){
  
  N <- nrow(tbl)
  mur_url <- MUR_nc_url(key$date[1])
  mur <- try(ncdf4::nc_open(mur_url))
  if (!inherits(mur, "try-error")){
      ilon <- closest_index(tbl$lon, mur$dim$lon$vals)
      ilat <- closest_index(tbl$lat, mur$dim$lat$vals)
      r <- dplyr::tibble(status = rep(TRUE, N))
      for (v in vars) r <- r %>% 
        dplyr::mutate(!!v := sapply(seq_len(N),
                                    function(i){
                                      ncdf4::ncvar_get(mur, varid = v,
                                                       start = c(ilon[i], ilat[i], 1),
                                                       count = c(1,1,1))
                                    }))
    
  } else {
    # no resource found
    r <- dplyr::tibble(status = rep(FALSE, N))
    for (v in vars) r <- r %>% 
        dplyr::mutate(!!v := rep(NA_real_, N))
  }
  
  tbl %>%
    dplyr::bind_cols(r)
}

#' Retrieve one or more points of MUR data
#'
#' @seealso https://podaac.jpl.nasa.gov/dataset/MUR-JPL-L4-GLOB-v4.1
#' @param x data frame with lon, lat and date columns
#' @param vars the name(s) of the variables to retrieve
#' \itemize{
#' \item{"analysed_sst" in kelvins}
#' \item{"analysis_error" in kelvins, estimated error standard deviation of analysed_sst}
#' \item{"mask" unitless, 1=open-sea, 2=land, 5=open-lake, 9=open-sea with ice in the grid, 13=open-lake with ice in the grid}
#' \item{"sea_ice_fraction", 0-1 but NA possible}
#' }
#' @param outfile character or NA, optional output path to the results as CSV
#' @return the input data frame with requested MUR data appended as a columns
MURpoints <- function(x,
                      vars = c("analysed_sst", 
                               "analysis_error", 
                               "mask",
                               "sea_ice_fraction"),
                      outfile = NA){

  if (!inherits(x$date, "Date")) x$date <- as.Date(x$date)
  
  x <- x %>%
    dplyr::group_by(date) %>%
    dplyr::group_map(fetch_MUR, .keep = TRUE) %>%
    dplyr::bind_rows()
  
  if (!is.na(outfile)) x <- readr::write_csv(x, outfile)
  
  x
}
