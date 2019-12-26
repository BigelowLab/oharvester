suppressPackageStartupMessages(library(ncdf4))
suppressPackageStartupMessages(library(raster))
mur_url <- file.path("https://podaac-opendap.jpl.nasa.gov/opendap/allData",
                     "ghrsst/data/GDS2/L4/GLOB/JPL/MUR/v4.1",
                     "2018/001/20180101090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc")
res <- c(0.01, 0.01)
varname <- "analysed_sst"
crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
bb <-  c(-72, -63, 39, 46)

x <- ncdf4::nc_open(mur_url)
lon <- x$dim$lon$vals
head(lon)
tail(lon)
lat <- x$dim$lat$vals
head(lat)
tail(lat)

half <- res/2                                            # half-res
pad_bb <- bb + c(-half[1], half[1], -half[2], half[2])   # pad bb
west_index  <- which.min(abs(lon - pad_bb[1]))           # closest lon to bb[1]
east_index  <- which.min(abs(lon - pad_bb[2]))           # closest lon to bb[2]
south_index <- which.min(abs(lat - pad_bb[3]))           # closest lat to bb[3]
north_index <- which.min(abs(lat - pad_bb[4]))           # closest lat to bb[4]

start <- c(west_index, south_index, 1)      # time layer starts at index 1
count <- c(
  east_index - west_index + 1,              # number of cells west-to-east
  north_index - south_index + 1,            # number of cells south-to-north
  -1                                        # -1 means get entire slice of time
)
start
count
m <- ncdf4::ncvar_get(x, varid = varname, start = start, count = count)
dim(m)
m[1:5, 1:5]
