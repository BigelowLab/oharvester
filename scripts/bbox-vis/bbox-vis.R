library(maps)

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

draw_basemap <- function(xlim = c(-180, 180), col = 'gray50',
                         bb = NULL){

  map(xlim = xlim,
      col = col)
  box(col = 'black')
  abline(v = mean(xlim), lty = 'dotted')
  if (!is.null(bb)){

  }


}

