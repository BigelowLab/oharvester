library(rwfs)
library(sf)

#' Fetch a table of data from Eye on Water
#'
#' Be advised that GDAL warns that the source data has an unknown CRS.  We
#' simply assign the
#'
#' @param uri character, the URL of the WFS service
#' @param layer_name character, the name of the layer to retrieve
#' @param crs character, proj string for CRS
#' @param path character or NA, if not NA then try to save the table to this location
#'        as <path>/layer_name.driver  Note that we have best luck when we specify
#'        a full path specification (/home/myname/data) rather than short hand
#'        ("./data" or "~/data").
#' @param driver character, one of the drivers available to \code{sf::writer_sf()}
#' @return simple feature POINT table
fetch_eow <- function(
  uri = file.path("https://geo-service.maris.nl/project/eyeonwater",
                  "wfs?request=getcapabilities&service=wfs"),
  layer_name = "eow_all",
  crs = "+proj=longlat +datum=WGS84",
  path = NA,
  driver = c("GeoJSON", "CSV")[1]){


  fileName <- tempfile()
  ok <- download.file(uri, fileName)
  request <- rwfs::GMLFile$new(fileName)
  client <- rwfs::WFSCachingClient$new(request)
  layers <- suppressWarnings(client$listLayers())
  if (!any(grepl(layer_name[1], layers$name, fixed = TRUE))){
    stop("layer name not found in available layers:", layer_name[1])
  }
  layer  <- suppressWarnings(client$getLayer(layer_name[1], crs = crs[1]))
  unlink(fileName)
  if (!is.na(path)){
    if (!dir.exists(path[1])) stop("path not found:", path[1])
    ofile <- paste0(layer_name[1],".", tolower(driver))
    sf::write_sf(layer, dsn = file.path(path, ofile), driver = driver)
  }
  layer
}

#' Plot the EOW observation locations with FU values scaled by color
#'
#' @param x sf table of points
#' @return leaflet map object
plot_eow <- function(x = fetch_eow()){
  stopifnot(require(leaflet))
  pal <- leaflet::colorNumeric(palette = "YlOrBr",
                               domain = c(0,25),
                               reverse = FALSE)
  leaflet::leaflet(layer) %>%
    addTiles() %>%
    addCircleMarkers(clusterOptions = markerClusterOptions(),
                     color = ~pal(fu_value)) %>%
    addLegend("bottomleft",
              pal = pal,
              values = ~fu_value,
              title = "FU Values",
              opacity = 1 )
}
