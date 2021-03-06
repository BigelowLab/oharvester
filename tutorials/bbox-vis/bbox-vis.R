library(maps)
library(sf)
library(dplyr)
library(raster)

#' Convert bounding box [0,360] longitudes to [-180, 180]
#'
#' Bounding boxes are 4 element vectors of [west, east, south, north]
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
#' Bounding boxes are 4 element vectors of [west, east, south, north]
#'
#' @export
#' @param x numeric bounding box vector, no check is done for being withing 0,360 range
#' @return numeric bounding box vector
to360BB <- function(x) {
  x[1:2] <- to360(x[1:2])
  if (identical(x[1], 360)) x[1] <- 0   # western edge
  if (identical(x[2], 0)) x[2] <- 360   # eastern edge
  
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


#' Convert a bounding box to a sf of type GEOMETRY
#' 
#' @param bb numeric, 4 element bouding box or a list with
#'        multiple bb vectors
#' @param crs character, coordinate reference
#' @return sf of geometry type POLYGON
bb_to_polygon <- function(bb = c(-170,-10,-60,60),
                          crs = "+init=epsg:4326"){
  if (!is.list(bb)) bb <- list(bb1 = bb)
  pp <- lapply(seq_len(length(bb)),
               function(i){
                 x <- bb[[i]]
                 sf::st_polygon(x = list(cbind(x[c(1,2,2,1,1)], 
                                               x[c(3,3,4,4,3)]))) %>%
                   sf::st_sfc(crs = crs) %>%
                   sf::st_sf() %>%
                   dplyr::mutate(ID = i)
               })
  do.call(rbind, pp)
}  


#' Split a bounding box into two at \code{at}
#' 
#' @param bb numeric, 4 element bouding box
#' @param at numeric, longitude to split around
#' @return list of one or more bounding box vectors
bb_split <- function(bb = c(-170,50,-60,60),
                     at = 0){
  if (bb[1] < at && bb[2] > at){
    x <- list(
      bb1 = c(bb[1], 0, bb[3:4]),
      bb2 = c(0, bb[2:4])
    )
  } else {
    x <- list(
      bb1 = bb
    )
  }
  x
}

#' Test if a blunding box straddles a longitude
#' 
#' @param bb numeric, 4 element bouding box
#' @param at numeric, longitude to split around
#' @return logical
bb_straddles <- function(bb = c(-170,50,-60,60),
                         at = 0){
  bb[1] < at && bb[2] > at
}

#' Draw a world map with optional raster and polygon overlays
#'
#' @param database character, either 'world' [-180,180] or 'world2' [0,360]
#' @param base_col character, color of base map vectors
#' @param R raster or NA, if raster then add to plot
#' @param bb list or numeric, bounding box or a list of bounding boxes
#' @param main character, optional title
#' @param lwd numeric, optipnal line weight for polygon edges
#' @param ofile charcater, optional PNG filename
draw_map <- function(database = 'world', 
                     base_col = 'gray50',
                     R = NA, 
                     bb = NULL, 
                     main = NA,
                     lwd = 2,
                     ofile = NA){
  
  if (!is.na(ofile)) png(ofile, width = 400, height = 200)
  xlim <- switch(database[1],
                 'world' = c(-180, 180),
                 "world2" = c(0, 360))
  map(database = database,
      xlim = xlim,
      col = base_col,
      mar = c(0.1,0.1,0.2,0.1))
  box(col = 'black')
  if (!is.null(R)){
    plot(R, add = TRUE, legend = FALSE)
    map(database = database,col = base_col)
  }
  if (!is.na(main)) title(main)
  abline(v = mean(xlim), lty = 'dotted')
  if (!is.null(bb)){
    pp <- bb_to_polygon(bb)
    plot(st_geometry(pp), 
         col = sf.colors(n = nrow(pp), alpha = 0.5, categorical = TRUE), 
         alpha = 0.7,
         border = "black", 
         lwd = 3,
         add = TRUE)
  }
  if (!is.na(ofile)) dev.off()
}