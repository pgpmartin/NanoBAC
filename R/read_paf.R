#' Read a .paf file such as produced by minimap/minimap2.
#' @description
#'     This function is only slighly modified from https://raw.githubusercontent.com/thackl/thacklr/master/R/read.R
#'     By default, this function corrects the start values to take into account that the paf file is 0-based while R/bioconductor is generally 1-based
#' @inheritParams readr::read_tsv
#' @importFrom readr read_tsv
#' @importFrom magrittr %>%
#' @importFrom methods is
#' @param max_tags maximum number of optional fields to include
#' @param fix0based Logical. Are the ranges in the paf file 0-based and should they be converted to 1-based (default is TRUE)?
#' @export
#' @return tibble
#' @examples
#' ## Example data set:
#'     Path2paf <- system.file("extdata", "pafFileExample.paf", package = "NanoBAC")
#' ## Import the data (ignore the parsing failures):
#'    mypaf <- read_paf(Path2paf)
#'

read_paf <- function (file, max_tags = 20, fix0based = TRUE){

  ## Fix the time zone issue if present:
  tt <- tryCatch(readr::default_locale(), error=function(e) e)

  if (methods::is(tt, "simpleError")) {
    options(readr.default_locale=readr::locale(tz=Sys.timezone()))
  }

  ## Define column names and types
  col_names <- c("query_name", "query_length", "query_start",
                 "query_end", "strand", "target_name", "target_length",
                 "target_start", "target_end", "map_match", "map_length",
                 "map_quality")
  col_types <- "ciiicciiiiin"

  if(max_tags > 0){
    col_names <- c(col_names, paste0("tag_", seq_len(max_tags)))
    col_types <- paste0(col_types, paste(rep("?", max_tags), collapse=""))
  }

  res <- readr::read_tsv(file, col_names = col_names, col_types = col_types) %>%
    tidy_paf_tags

  # Fix paf file (0-based) in order to respect the 1-based system used in bioconductor
  if (fix0based) {
    res$query_start <- res$query_start + 1
    res$target_start <- res$target_start + 1
  }

  return(res)

}

#' Tidy the flag columns imported from a paf file
#'
#' @importFrom tidyselect starts_with
#' @importFrom dplyr select mutate_at bind_cols
#' @importFrom magrittr %>%
#' @importFrom rlang warn inform
#' @importFrom tibble has_name
#' @importFrom stringr str_split str_glue
#' @importFrom stats na.omit
#' @param Somedata A \code{tibble} as imported during the [NanoBAC::read_paf()] function
#' @details
#'     This function is not intended to be used directly but is used in the [NanoBAC::read_paf()] function
#' @return tibble

tidy_paf_tags <- function(Somedata){
  tag_df <- tibble::tibble(.rows=nrow(Somedata))
  tag_types <- c()
  seen_empty_tag_col <- FALSE

  for (x in dplyr::select(Somedata, tidyselect::starts_with("tag_"))){
    tag_mx <- stringr::str_split(x, ":", 3, simplify = TRUE)
    tag_mx_nr <- stats::na.omit(unique(tag_mx[,1:2]))
    if(nrow(tag_mx_nr) == 0){
      seen_empty_tag_col <- TRUE
      break; # empty col -> seen all tags
    }
    tags <- tag_mx_nr[,1]
    tag_type <- tag_mx_nr[,2]
    names(tag_type) <- tags
    # add to global tag_type vec
    tag_types <- c(tag_types, tag_type)
    tag_types <- tag_types[unique(names(tag_types))]
    # sort tag values into tidy tag tibble
    for (tag in tags){
      if(!tibble::has_name(tag_df, tag)){ # init tag
        tag_df[[tag]] <- NA
      }
      tag_idx <- tag_mx[,1] %in% tag
      tag_df[[tag]][tag_idx] <- tag_mx[tag_idx,3]
    }
  }

  tag_df <- tag_df %>%
    dplyr::mutate_at(names(tag_types)[tag_types == "i"], as.integer) %>%
    dplyr::mutate_at(names(tag_types)[tag_types == "f"], as.numeric)

  if(!seen_empty_tag_col)
    rlang::warn("Found tags in max_tags column, you should increase max_tags to ensure all tags for all entries got included")

  rlang::inform(
    stringr::str_glue("Read and tidied up a .paf file with {n_tags} optional tag fields:\n{s_tags}",
                      s_tags = toString(names(tag_types)),
                      n_tags = length(tag_types)))

  dplyr::bind_cols(dplyr::select(Somedata, -tidyselect::starts_with("tag_")), tag_df)
}

