
#' Fetch data from the Australian Faunal Directory
#'
#' @param taxa (Character) A query taxonomic name, or vector of taxonomic names
#' @param retry_attempt (Integer) The number of query attempts in case of query failure due to poor internet connection.
#' @param retry_wait (Integer) How long to wait between query attempts.
#' @param quiet (Logical) Option to print warnings to console
#']
#' @return A dataframe containing checklist data for a taxon, or a list of data frames if multiple taxon names were provided.
#' @export
#' @importFrom purrr map
#' @importFrom purrr map_dfr
#' @importFrom dplyr pull
#' @importFrom methods is
#'
#'
#' @examples
fetch_afd_checklist <- function(taxa, retry_attempt=3, retry_wait=5, quiet=FALSE){
  if(!is.character(taxa)){
    stop("taxa must be a character vector containing family, genus, or species names")
  }
  dupvec <- split(seq_along(taxa), taxa)
  out <- vector("list", length(taxa))
  names(out)<- taxa
  for(i in 1:length(dupvec)){
    index <- dupvec[[i]]
    # First check if the taxonomic name is valid
    taxon <- check_afd_name(names(dupvec)[i], quiet=quiet)

    # Check the query is valid
    query <- check_afd_query(taxon)

    need_child <- taxon[query == FALSE]
    if(length(need_child) > 0){
      to_download <- list()
    } else {
      to_download <- taxon
    }

    while(length(need_child) > 0){
      child <- need_child %>%
        purrr::map(~{
          .x %>%
            get_afd_children %>%
            dplyr::pull(nameKey)
        }) %>%
        unlist()

      query_child <- child %>%
        purrr::map(check_afd_query, check_name = FALSE)%>%
        unlist()

      # Add successfully resolved to to_download
      to_download <- c(to_download, child[!query_child==FALSE] )

      # Repeat on unresolved
      need_child <- child[query_child==FALSE]
    }

    # Retrieve CSV files
    out[[index]] <- to_download %>%
      purrr::map_dfr(fetch_afd_csv, retry_attempt = retry_attempt, retry_wait=retry_wait, quiet=quiet, check_name=FALSE)
  }
  if(methods::is(out, "list") & !length(out) > 1){
    out <- out[[1]]
  }
  return(out)
}


#' Get the taxonomic rank of a query on AFD
#'
#' @param taxa (Character) A query taxonomic name, or vector of taxonomic names
#' @param quiet (Logical) Option to print warnings to console
#'
#' @return A character vector of the same length as `taxa` containing the taxonomic rank/s of the query taxa.
#' @export
#' @importFrom rvest read_html
#' @importFrom rvest html_nodes
#' @importFrom rvest html_text
#' @importFrom purrr pluck
#' @importFrom stringr str_remove
#'
#' @examples
get_afd_rank <- function(taxa, quiet=FALSE){
  if(!is.character(taxa)){
    stop("taxa must be a character vector containing family, genus, or species names")
  }
  dupvec <- split(seq_along(taxa), taxa)
  out <- vector("character", length(taxa))
  names(out)<- taxa
  for(i in 1:length(dupvec)){
    index <- dupvec[[i]]
    taxname <- check_afd_name(names(dupvec)[i], quiet=quiet)
    if(!is.null(taxname)){
      url <- paste0("https://biodiversity.org.au/afd/taxa/", taxname)
      out[index] <- rvest::read_html(url) %>%
        rvest::html_nodes('h1') %>%
        rvest::html_text() %>%
        purrr::pluck(2) %>%
        stringr::str_remove(" .*$")
    } else{
      out[index] <- NULL
     }
   }
  return(out)
}

#' Retrieve the immediate children of a query taxa
#'
#' @param taxa (Character) A query taxonomic name, or vector of taxonomic names
#' @param quiet (Logical) Option to print warnings to console
#' @return A dataframe containing taxonomic children data for a taxon, or a list of data frames if multiple taxon names were provided.
#' @export
#' @importFrom jsonlite fromJSON
#'
#' @examples
get_afd_children <- function(taxa, quiet=FALSE){
  if(!is.character(taxa)){
    stop("taxa must be a character vector containing family, genus, or species names")
  }
  dupvec <- split(seq_along(taxa), taxa)
  out <- vector("list", length(taxa))
  names(out)<- taxa
  for(i in 1:length(dupvec)){
    index <- dupvec[[i]]
    taxname <- check_afd_name(names(dupvec)[i], quiet=quiet)
    if(!is.null(taxname)){
      url <- paste0("https://biodiversity.org.au/afd/taxa/", taxname, "/checklist-subtaxa.json")
      raw_json <- jsonlite::fromJSON(url)
      if(length(raw_json) > 0){
        out[[index]] <- raw_json$metadata
      } else{
        out[index] <- list(NULL)
      }
    } else{
      out[index] <- list(NULL)
    }
  }
  if(is(out, "list") & !length(out) > 1){
    out <- out[[1]]
  }
  return(out)
}


#' Get occurance data from the Australian Faunal Directory
#'
#' @param taxa (Character) A query taxonomic name, or vector of taxonomic names
#' @param type The type of occurance data to return. Options are "states" for just the state occurance data,
#' or "ibra" for Australia's Bioregions occurance data (see https://www.environment.gov.au/land/nrs/science/ibra)
#' @param quiet (Logical) Option to print warnings to console
#'
#' @return A character vector containing taxonomic children data for a taxon, or a list of character vectors if multiple taxon names were provided.
#' @export
#' @importFrom rvest read_html
#' @importFrom rvest html_nodes
#' @importFrom rvest html_text
#' @importFrom stringr str_to_lower
#' @importFrom stringr str_trim
#' @importFrom stringr str_remove_all
#' @importFrom stringr str_split
#' @importFrom purrr pluck
#'
#' @examples
get_afd_occurance <- function(taxa, type="states", quiet=FALSE){
  if(!is.character(taxa)){
    stop("taxa must be a character vector containing family, genus, or species names")
  }
  dupvec <- split(seq_along(taxa), taxa)
  out <- vector("list", length(taxa))
  names(out)<- taxa
  for(i in 1:length(dupvec)){
    index <- dupvec[[i]]
    taxname <- check_afd_name(names(dupvec)[i], quiet=quiet)
    if(is.null(taxname) | length(taxname) < 1){
      out[index] <- list(NULL)
    }
    url <- paste0("https://biodiversity.org.au/afd/taxa/", taxname)

    query_data <- rvest::read_html(url) %>%
      rvest::html_nodes('h4, p')

    #get the index of the p div to pluck
    if(stringr::str_to_lower(type) == "states"){
      pdiv <- which(query_data %>% rvest::html_text() =="States") + 1
      if(length(pdiv)== 0){
        # If its cosmopolitan its put in a different div for some reason
        pdiv <- which(query_data %>% rvest::html_text() =="Extra Distribution Information") + 1
      }
    } else if(stringr::str_to_lower(type) == "ibra"){
      pdiv <- which(query_data %>% rvest::html_text() =="IBRA") + 1
    } else {
      stop("type must be either 'states' or 'ibra'")
    }
    if (!length(pdiv) > 0){
      if(!quiet){warning("No type data available for: ", taxname)}
      out[index] <- list(NULL)
    }
    # pull out the pdiv line and clean
    dist_data <- query_data %>%
      purrr::pluck(pdiv) %>%
      rvest::html_text()%>%
      stringr::str_trim(side="both") %>%
      stringr::str_remove_all("\n +")%>%
      stringr::str_remove_all("  +")%>%
      stringr::str_remove_all("\t")

    out[[index]] <- stringr::str_split(dist_data, ",", n=Inf) %>%
      unlist()%>%
      stringr::str_trim(side="both")
  }
  if(is(out, "list") & !length(out) > 1){
    out <- unlist(out)
    names(out) <- NULL
  }
  return(out)
}


#' Check if name is valid on the Australian Faunal Directory
#'
#' @param taxa (Character) A query taxonomic name, or vector of taxonomic names
#' @param quiet (Logical) Option to print warnings to console
#'
#' @return (character) A character vector of the same length as `taxa` containing the valid taxonomic name/s on AFD, or NULL if the name was not found.
#' @importFrom rvest read_html
#' @importFrom rvest html_nodes
#' @importFrom rvest html_text
#' @importFrom rvest html_table
#' @importFrom rvest html_attr
#' @importFrom purrr pluck
#' @importFrom purrr map
#' @importFrom dplyr mutate_all
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom dplyr slice_min
#' @importFrom dplyr mutate
#' @importFrom stringr str_trim
#' @importFrom stringr str_detect
#' @importFrom stringr str_remove
#' @importFrom stringr str_remove_all
#' @importFrom stringr str_to_lower
#' @importFrom tibble rownames_to_column
#'
#' @examples
check_afd_name <- function(taxa, quiet=FALSE){
  if(!is.character(taxa)){
    stop("taxa must be a character vector containing family, genus, or species names")
  }
  dupvec <- split(seq_along(taxa), taxa)
  out <- vector("character", length(taxa))
  names(out)<- taxa
  for(i in 1:length(dupvec)){
    index <- dupvec[[i]]
    taxon <- names(dupvec)[i]
    if(!is.null(attr(taxon, "valid_name"))){
      return(taxon)
    }
    taxon <- taxon %>% stringr::str_replace_all("_", " ")
    if(stringr::str_detect(taxon, " ")){
      taxon <- taxon %>% stringr::str_replace_all(" ", "%25")
    }
    url <- paste0("https://biodiversity.org.au/afd/search/names?search=true&keyword=", taxon)

    # Check page type (single result or multiple choice)
    table_header <- rvest::read_html(url) %>%
      rvest::html_nodes('th') %>%
      rvest::html_text()

    if(length(table_header)==0){
      # Check why it could have failed
      failed_reason <- rvest::read_html(url) %>%
        rvest::html_nodes('h3') %>%
        rvest::html_text() %>%
        stringr::str_replace_all("%", " ")
      if(!quiet){warning(failed_reason)}
      return(NULL)
    }

    if(table_header[1] == "Published"){
      # taxon is valid and unique
      out[index] <- taxon
    } else if (table_header[1] == "Rank"){
      # taxon is not unique
      table_data <- rvest::read_html(url) %>%
        rvest::html_table() %>%
        purrr::pluck(1) %>%
        dplyr::mutate_all(~{
          .x %>%
            stringr::str_trim(side="both") %>%
            stringr::str_remove_all("\n +")
        })

      # Get the taxon links
      rows <-  rvest::read_html(url)  %>%
        rvest::html_nodes("table") %>%
        rvest::html_nodes("tr") %>%
        purrr::map(function(x){
          x %>%
            rvest::html_nodes( "a") %>%
            rvest::html_attr("href") %>%
            purrr::pluck(1)
        }) %>%
        unlist()

      # Get the index of the valid taxon at the highest rank
      possible_ranks <- c("Order", "Suborder", "Family","Subfamily","Genus","Subgenus", "Species", "Subspecies")
      name_check <- taxon %>%
        stringr::str_replace("%[0-9]+", "_") %>%
        stringr::str_split(pattern="_", n=Inf) %>%
        unlist()

      if(any(table_data$Status == "Valid Name")){
        # At least one returned name is valid
        valid_name <- TRUE
      } else if(!any(table_data$Status == "Valid Name") & any(stringr::str_detect(table_data$Status, "Generic Combination"))){
        # Returned name is generic combination (synonym)
        valid_name <- FALSE
      } else{
        return(NULL)
      }

      index_to_keep <- table_data %>%
        tibble::rownames_to_column("index") %>%
        dplyr::mutate(rank_n = match(Rank, possible_ranks)) %>%
        dplyr::filter(dplyr::case_when(
          valid_name ~ Status == "Valid Name",
          !valid_name ~ stringr::str_detect(table_data$Status, "Generic Combination")
        ))%>%
        dplyr::mutate(name_check = purrr::map(Name, function(x){
          if(valid_name){
            all(stringr::str_detect(stringr::str_to_lower(x), paste0("\\b", stringr::str_to_lower(name_check), "\\b")))
          } else {
            any(stringr::str_detect(stringr::str_to_lower(x), paste0("\\b", stringr::str_to_lower(name_check), "\\b")))
          }
        })) %>%
        tidyr::unnest(name_check) %>%
        dplyr::filter(name_check) %>%
        dplyr::slice_min(rank_n) %>%
        dplyr::pull(index) %>%
        as.numeric()

      out[index] <- stringr::str_remove(rows[index_to_keep], "../taxa/")
    }
  }
  attr(out, "valid_name") <- TRUE
  return(out)
}

#' Check if taxonomic name is present in the Australian Faunal Directory
#'
#' @param taxa (Character) A query taxonomic name, or vector of taxonomic names
#' @param quiet (Logical) Option to print warnings to console
#'
#' @return (character) A valid taxon name
#' @export
#'
check_afd_presence <- function(taxa, quiet=FALSE) {
  if(!is.character(taxa)){
    stop("taxa must be a character vector containing family, genus, or species names")
  }
  dupvec <- split(seq_along(taxa), taxa)

  out <- vector("logical", length(taxa))
  names(out)<- taxa

  for(i in 1:length(dupvec)){
    index <- dupvec[[i]]
    to_search <- names(dupvec[i])
    out[index] <- !is.null(check_afd_name(to_search, quiet=quiet))
  }
  return(out)
}



#' Internal function to check that query is not too large for AFD server
#'
#' @param taxon (Character) A query taxon
#' @param check_name (Logical) Whether to check the name is valid first using check_afd_name
#'
#' @return a valid URL query, or FALSE if the query is too large
#' @importFrom rvest read_html
#' @importFrom rvest html_nodes
#' @importFrom rvest html_text
#' @importFrom purrr pluck
#' @importFrom stringr str_detect
#'
#' @examples
check_afd_query <- function(taxon, check_name = TRUE){
  # Check name is compatible
  if(is.null(attr(taxon, "valid_name")) & isTRUE(check_name)){
    taxon <- check_afd_name(taxon)
  }
  # Check if query is of appropriate size (<1000 records)
  url <- paste0("https://biodiversity.org.au/afd/taxa/", taxon, "/names/csv/")
  html_header <- rvest::read_html(url) %>%
    rvest::html_nodes('h1') %>%
    rvest::html_text() %>%
    purrr::pluck(2)
  if(stringr::str_detect(html_header, "limit exceeded")){
    return(FALSE)
  } else if(!stringr::str_detect(html_header, "limit exceeded")) {
    return(url)
  }
}

#' Internal function to download CSV from AFD
#'
#' @param taxon (Character) A query taxon
#' @param retry_attempt (Integer) The number of query attempts in case of query failure due to poor internet connection.
#' @param retry_wait (Integer) How long to wait between query attempts.
#' @param quiet (Logical) Whether the function should print progress to console.
#' @param check_name (Logical) Whether to check the name is valid first using check_afd_name
#'
#' @return
#' @importFrom readr read_csv
#' @importFrom janitor clean_names
#' @import readr
#'
#' @examples
fetch_afd_csv <- function(taxon, retry_attempt=3, retry_wait=5, quiet=FALSE, check_name=TRUE){
  # Check name is compatible
  if(check_name){
    taxon <- check_afd_name(taxon)
  }
  # check if query is of appropriate size (<1000 records)
  url <- check_afd_query(taxon, check_name = check_name)
  if(is.character(url)){
    csv <- NULL
    attempt <- 1
    while(is.null(csv) && attempt <= (retry_attempt+1)) {
      csv <- tryCatch({
        readr::read_csv(
          paste0(url,taxon,".csv" ),
          col_types = readr::cols(
            CAVS_CODE = readr::col_number(),
            CAAB_CODE = readr::col_number(),
            NAMES_VARIOUS = readr::col_character(),
            SCIENTIFIC_NAME = readr::col_character(),
            FAMILY = readr::col_character(),
            GENUS = readr::col_character(),
            SUBGENUS = readr::col_character(),
            SPECIES = readr::col_character(),
            SUBSPECIES = readr::col_character(),
            NAME_TYPE = readr::col_character(),
            NAME_SUBTYPE = readr::col_character(),
            RANK = readr::col_character(),
            QUALIFICATION = readr::col_character(),
            AUTHOR = readr::col_character(),
            YEAR = readr::col_number(),
            ORIG_COMBINATION = readr::col_character(),
            NAME_GUID = readr::col_character(),
            NAME_LAST_UPDATE = readr::col_datetime(),
            TAXON_GUID = readr::col_character(),
            TAXON_LAST_UPDATE = readr::col_datetime(),
            PARENT_TAXON_GUID = readr::col_character(),
            CONCEPT_GUID = readr::col_character(),
            PUB_PUB_AUTHOR = readr::col_character(),
            PUB_PUB_YEAR = readr::col_character(),
            PUB_PUB_TITLE = readr::col_character(),
            PUB_PUB_PAGES = readr::col_character(),
            PUB_PUB_PARENT_BOOK_TITLE = readr::col_character(),
            PUB_PUB_PARENT_JOURNAL_TITLE = readr::col_character(),
            PUB_PUB_PARENT_ARTICLE_TITLE = readr::col_character(),
            PUB_PUB_PUBLICATION_DATE = readr::col_character(),
            PUB_PUB_PUBLISHER = readr::col_character(),
            PUB_PUB_FORMATTED = readr::col_character(),
            PUB_PUB_QUALIFICATION = readr::col_character(),
            PUB_PUB_TYPE = readr::col_character(),
            PUBLICATION_GUID = readr::col_character(),
            PUBLICATION_LAST_UPDATE = readr::col_datetime(),
            PARENT_PUBLICATION_GUID = readr::col_character()
          ))
      }, error = function(e){
        if (!quiet) {cat(paste("Failed attempt ", attempt,"\n"))}
        Sys.sleep(retry_wait)
        NULL
      })
      if(!is.null(csv)){
        out <- csv %>%
          janitor::clean_names()
      }else {
        attempt <- attempt + 1
      }
    }
    if(is.null(csv)) {
      warning("Failed to download ", taxon, " after ", attempt-1, " attempts")
      out <- NULL
    }
  } else if (isFALSE(url)){
    message("Query limit of of 1000 taxa exceeded")
    out <- NULL
  }
  return(out)
}

