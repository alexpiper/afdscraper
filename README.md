
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

## Retrieving checklist data for large taxonomic groups

One of the main uses of afdscraper is to retrieve CSV data for a large
amount of taxa, overcoming the current 1000 taxon limit on the website.

``` r
# Load the afdscraper package
library(afdscraper)

# Also load tidyverse packages for their handy filtering and summarising functiosn
library(tidyverse)

# Download all checklist data for Insecta
afd_insecta <- fetch_afd_checklist("Insecta", retry_attempt=3, retry_wait=5, quiet=FALSE)

head(afd_insecta)
```

## Checking presence of a taxon in Australia

The presence of a species on AFD can be checked. This is useful for
automated cross referencing of species detections to identify potential
new species introductions.

``` r
check_afd_presence("Bactrocera tryoni")
```

## Get the taxonomic rank of a query taxon

The taxonomic rank of a query taxon name can be obtained.

``` r
get_afd_rank("Scaptodrosophila")
```

## Retrieve the immediate children of a query taxon

The taxonomic children can be retrieved for taxonomic ranks above
species.

``` r
get_afd_children("Tephritidae")
```

## Retrieve distribution information

Distribution data can be retrieved for taxa, either by states (type =
“states”), or by [Australia’s
Bioregions](https://www.environment.gov.au/land/nrs/science/ibra) (type
= “ibra”). Note that occurance data is lacking for many of the taxa
listed on AFD.

``` r
get_afd_occurance("Bactrocera tryoni", type = "states") 
```

### Multiple taxa

All of the above functions can also handle multiple taxa provided to the
function or read in from a file.

``` r
taxa <- c("Bactrocera tryoni", "Scaptodrosophila", "Nitidulidae")
check_afd_presence(taxa)

# Read in a list from a file
taxa <- readLines("taxon_list.txt")
check_afd_presence(taxa)
```
