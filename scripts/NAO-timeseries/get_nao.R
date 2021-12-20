library(stringr)
library(readr)
library(dplyr)
library(tidyr)

### https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/nao.shtml
nao <- readLines("https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/norm.nao.monthly.b5001.current.ascii.table")


wide_form <- nao  %>%
  stringr::str_trim(side = "both") %>%
  stringr::str_squish() %>%
  stringr::str_replace_all(stringr::fixed(" "), stringr::fixed(",")) %>%
  paste(collapse = "\n") %>%
  readr::read_csv(col_names = c("Year", month.abb),
                  col_types = strrep("d", 13))

long_form <- wide_form %>%
  tidyr::pivot_longer(all_of(month.abb), 
                      names_to = "Month",
                      values_to = "AMO")
