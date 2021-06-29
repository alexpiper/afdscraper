
<!-- README.md is generated from README.Rmd. Please edit that file -->

# afdscraper

<!-- badges: start -->
<!-- badges: end -->

afdscraper is an r package for scraping data from the [Australian Faunal
Directory](https://biodiversity.org.au/afd/home)

## Installation

You can install the development version of afdscraper from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("alexpiper/afdscraper")
```

# Working examples

afdscraper can be used for a number of different applications

## Checking presence

## Retrieving all data

One of the main uses of afdscraper is to retrieve CSV data for a large
amount of taxa, overcoming the current 1000 taxon limit on the website.

``` r
# Load the afdscraper package
library(afdscraper)

# Also load tidyverse packages for their handy filtering and summarising functiosn
library(tidyverse)

afd_insecta <- fetch_afd_checklist("Insecta", retry_attempt=3, retry_wait=5, quiet=FALSE)

head(afd_insecta)
```
