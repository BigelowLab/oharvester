library(stringr)
library(readr)
library(dplyr)
library(tidyr)

amo <- readLines("https://psl.noaa.gov/data/correlation/amon.us.long.data") 

dims <- amo[1] %>%
  stringr::str_trim(side = "both") %>%
  stringr::str_squish() %>%
  stringr::str_split(pattern = " ") %>%
  `[[`(1) %>%
  as.numeric()

n <- dims[2]-dims[1] + 2

wide_form <- amo[2:n] %>%
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
  
