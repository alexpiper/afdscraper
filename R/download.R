
#' Fetch data from the Australian Faunal Directory
#'
#' @param taxon (Character) A query taxon
#' @param retry_attempt (Integer) The number of query attempts in case of query failure due to poor internet connection.
#' @param retry_wait (Integer) How long to wait between query attempts.
#' @param quiet (Logical) Whether the function should print progress to console.
#'
#' @return A dataframe containing
#' @export
#' @importFrom purrr map
#' @importFrom purrr map_dfr
#' @importFrom dplyr pull
#'
#'
#' @examples
fetch_afd_checklist <- function(taxon, retry_attempt=3, retry_wait=5, quiet=FALSE){

  # First check if the query is valid
  taxon <- check_afd_name(taxon)

  # Try a query, see if it fails
  query <- check_afd_query(taxon)

  need_child <- taxon[query==FALSE]
  if(length(need_child) > 0){
    to_download <- list()
  } else {
    to_download <- query
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
  out <- to_download %>%
    purrr::map_dfr(fetch_afd_csv, retry_attempt=retry_attempt, retry_wait=retry_wait, quiet=quiet, check_name=FALSE)
  return(out)
}


#' Get the taxonomic rank of a query on AFD
#'
#' @param taxon (Character) A query taxon
#'
#' @return
#' @export
#' @importFrom rvest read_html
#' @importFrom rvest html_nodes
#' @importFrom rvest html_text
#' @importFrom purrr pluck
#' @importFrom stringr str_remove
#'
#' @examples
get_afd_rank <- function(taxon){
  url <- paste0("https://biodiversity.org.au/afd/taxa/", taxon)
  rank <- rvest::read_html(url) %>%
    rvest::html_nodes('h1') %>%
    rvest::html_text() %>%
    purrr::pluck(2) %>%
    stringr::str_remove(" .*$")
  return(rank)
}

#' Retrieve the immediate children of a query taxa
#'
#' @param taxon (Character) A query taxon
#'
#' @return a dataframe containing child data
#' @export
#' @importFrom jsonlite fromJSON
#'
#' @examples
get_afd_children <- function(taxon){
  url <- paste0("https://biodiversity.org.au/afd/taxa/", taxon, "/checklist-subtaxa.json")
  raw_json <- jsonlite::fromJSON(url)
  children <- raw_json$metadata
  return(children)
}


#' Get occurance data from the Australian Faunal Directory
#'
#' @param taxon (Character) A query taxon
#' @param occurance The type of occurance data to return. Options are "states" for just the state occurance data,
#' or "ibra" for Australia's Bioregions occurance data (see https://www.environment.gov.au/land/nrs/science/ibra)
#'
#' @return
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
get_afd_occurance <- function(taxon, occurance="states"){
  taxon <- check_afd_name(taxon)
  if(is.null(taxon) | length(taxon) < 1){
    return(NA)
  }
  url <- paste0("https://biodiversity.org.au/afd/taxa/", taxon)

  query_data <- rvest::read_html(url) %>%
    rvest::html_nodes('h4, p')

  #get the index of the p div to pluck
  if(stringr::str_to_lower(occurance) == "states"){
    index <- which(query_data %>% rvest::html_text() =="States") + 1
    if(length(index)== 0){
      # If its cosmopolitan its put in a different div for some reason
      index <- which(query_data %>% rvest::html_text() =="Extra Distribution Information") + 1
    }
  } else if(stringr::str_to_lower(occurance) == "ibra"){
    index <- which(query_data %>% rvest::html_text() =="IBRA") + 1
  } else {
    message("occurance must be either 'states' or 'ibra'")
  }
  if (!length(index) > 0){
    message("No occurance data available for: ", taxon)
    return(NA)
  }
  # pull out the index line and clean
  dist_data <- query_data %>%
    purrr::pluck(index) %>%
    rvest::html_text()%>%
    stringr::str_trim(side="both") %>%
    stringr::str_remove_all("\n +")%>%
    stringr::str_remove_all("  +")%>%
    stringr::str_remove_all("\t")

  out <- stringr::str_split(dist_data, ",", n=Inf) %>%
    unlist()%>%
    stringr::str_trim(side="both")
  return(out)
}


#' Check if name is valid on the Australian Faunal Directory
#'
#' @param taxon (Character) A query taxon
#'
#' @return (character) A valid taxon name
#' @export
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
check_afd_name <- function(taxon){
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
    warning(failed_reason)
    return(NULL)
  }

  if(table_header[1] == "Published"){
    # taxon is valid and unique
    out <- taxon
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

    out <- stringr::str_remove(rows[index_to_keep], "../taxa/")
  }
  attr(out, "valid_name") <- TRUE
  return(out)
}


#' Check query is not too large for the AFD server
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
#' @importFrom readr cols
#'
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
            CAVS_CODE = col_number(),
            CAAB_CODE = col_number(),
            NAMES_VARIOUS = col_character(),
            SCIENTIFIC_NAME = col_character(),
            FAMILY = col_character(),
            GENUS = col_character(),
            SUBGENUS = col_character(),
            SPECIES = col_character(),
            SUBSPECIES = col_character(),
            NAME_TYPE = col_character(),
            NAME_SUBTYPE = col_character(),
            RANK = col_character(),
            QUALIFICATION = col_character(),
            AUTHOR = col_character(),
            YEAR = col_number(),
            ORIG_COMBINATION = col_character(),
            NAME_GUID = col_character(),
            NAME_LAST_UPDATE = col_datetime(),
            TAXON_GUID = col_character(),
            TAXON_LAST_UPDATE = col_datetime(),
            PARENT_TAXON_GUID = col_character(),
            CONCEPT_GUID = col_character(),
            PUB_PUB_AUTHOR = col_character(),
            PUB_PUB_YEAR = col_character(),
            PUB_PUB_TITLE = col_character(),
            PUB_PUB_PAGES = col_character(),
            PUB_PUB_PARENT_BOOK_TITLE = col_character(),
            PUB_PUB_PARENT_JOURNAL_TITLE = col_character(),
            PUB_PUB_PARENT_ARTICLE_TITLE = col_character(),
            PUB_PUB_PUBLICATION_DATE = col_character(),
            PUB_PUB_PUBLISHER = col_character(),
            PUB_PUB_FORMATTED = col_character(),
            PUB_PUB_QUALIFICATION = col_character(),
            PUB_PUB_TYPE = col_character(),
            PUBLICATION_GUID = col_character(),
            PUBLICATION_LAST_UPDATE = col_datetime(),
            PARENT_PUBLICATION_GUID = col_character()
          ))
      }, error = function(e){
        if (!quiet) {cat(paste("Failed attempt ", attempt,"\n"))}
        Sys.sleep(retry_wait)
        NULL
      })
      if(!is.null(csv)){
        out <- csv

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

