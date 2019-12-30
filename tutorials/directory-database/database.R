# Create a dummy dataset using dplyr's storms data
# @param the output datapath
# @return database of datasets as a tibble
generate_dummy_database <- function(path = "data"){
  stopifnot(require(dplyr))
  stopifnot(require(readr))
  if (!dir.exists(path)) stopifnot(dir.create(path,
                                              recursive = TRUE,
                                              showWarnings = FALSE))

  # create output subdirectory
  # @param .x tibble of data for year/month grouping
  # @param .y tibble of data with grouping cols (year and month)
  # @param path the root directory
  do_subdir <- function(.x, .y, path = "."){
    subdir <- file.path(path, .y$year[1], .y$month[1])
    if (!dir.exists(subdir)) stopifnot(dir.create(subdir,
                                                  recursive = TRUE,
                                                  showWarnings = FALSE))
  } # do_subdir

  # write a file to the appropriate path
  # @param x table of data
  # @param path the root directory
  do_day <- function(x, path = "."){
    opath <- file.path(path, x$year[1], x$month[1])
    ofile <- sprintf("%s-%s-%s-storms.csv",x$year[1], x$month[1], x$day[1])
    readr::write_csv(x, file.path(opath, ofile))
  }

  x <- dplyr::storms %>%
    dplyr::mutate(year  = sprintf("%0.4i", year),
                  month = sprintf("%0.2i", month),
                  day   = sprintf("%0.2i", day)) %>%
    dplyr::filter(year %in% c("2014", "2015")) %>%
    dplyr::group_by(year, month) %>%
    dplyr::group_walk(do_subdir, path = path) %>%
    dplyr::group_by(year, month, day) %>%
    dplyr::group_split() %>%
    lapply(do_day, path = path) %>%
    dplyr::bind_rows() %>%
    dplyr::select(year, month, day) %>%
    dplyr::distinct() %>%
    write_csv(file.path(path, "database.csv"))

  x
}
