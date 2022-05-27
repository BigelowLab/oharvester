suppressPackageStartupMessages({
  required <- c("dplyr","readr", "rnoaa", "ggplot2")
  installed <- installed.packages() |>
    dplyr::as_tibble() 
  needed <- !(required %in% installed$Package)
  if (any(needed)){
    install.packages(required[needed])
  }
  ok <- sapply(required, library, character.only = TRUE)
})



#' Retrieve data for a given buoy and dataset type for one or more years
#' 
#' @seealso https://github.com/BigelowLab/oharvester/tree/master/tutorials/NDBC-buoydata
#' @param dataset character, see \code{?rnoaa::buoy}
#' @param buoyid character, see \code{?rnoaa::buoy}
#' @param years numeric, one or more years to retrieve.  Note the example here
#'   shows noncontiguous years
#' @param outfile character or NA, in not NA, then write the output to this file as CSV.
#'   Use extension ".csv.gz" to save the CSV gzipped.
#' @return a data frame (tibble) of results
get_buoy_data <- function(dataset = 'stdmet',
                          buoyid = "CFWM1",
                          years = c(2012:2015, 2017:2022),
                          outfile = NA){
  
  x <- lapply(years,
    function(year, dataset = "", buoyid = ""){
      rnoaa::buoy(dataset[1], buoyid[1], year = year)$data
    }, dataset = dataset, buoyid = buoyid) |>
    dplyr::bind_rows() |>
    dplyr::mutate(time = as.POSIXct(time, format = "%Y-%m-%dT%H:%M:%SZ"))
  
  if (!is.na(outfile)) x <- readr::write_csv(x, outfile)
}