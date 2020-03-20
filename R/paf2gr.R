#' Convert an paf file imported with \code{\link[NanoBAC]{read_paf}} to a \code{GRanges}
#'
#' @description
#'     the ranges in the resulting \code{GRanges} object correspond either to the reads or to the target depending on the value of the focus argument
#' @importFrom S4Vectors DataFrame mcols
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb Seqinfo
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr pull distinct select
#' @param paf A tibble with a paf file imported with the read_paf function
#' @param focus Character string. Either "target" or "reads" indicating which ranges will be represented in the GRanges
#' @param includetags Logical. Should the tags in the paf file be present in mcols from the GRanges (default is TRUE)?
#' @export
#' @return a \code{GRanges} object
#' @seealso \code{\link[NanoBAC]{read_paf}}
#' @examples
#' ## Example data set:
#'     Path2paf <- system.file("extdata", "pafFileExample.paf", package = "NanoBAC")
#' ## Import the data (ignore the parsing failures):
#'    mypaf <- read_paf(Path2paf)
#' ## Convert the dataset to a GRanges focusing on the reads
#'    mypafgr <- paf2gr(mypaf)
#'

paf2gr <- function(paf,
                   focus = c("target", "reads"),
                   includetags = TRUE) {

  stdColnames <- 	c("query_name", "query_length", "query_start",
                    "query_end", "strand", "target_name",
                    "target_length", "target_start", "target_end",
                    "map_match", "map_length", "map_quality")
  stopifnot(ncol(paf) >=12)
  stopifnot(colnames(paf)[1:12] == stdColnames)

  focus <- match.arg(focus, several.ok = FALSE)

  if (focus == "target") {
    #Metadata:
    if (includetags && ncol(paf)>12) {
      metacols <- S4Vectors::DataFrame(paf[,c(1, 10:ncol(paf))])
      colnames(metacols) <- c("ReadName", colnames(metacols)[-1])
    } else {
      metacols <- S4Vectors::DataFrame(paf[,c(1:4, 10:12)])
      colnames(metacols) <- c("ReadName", "ReadLength",
                              "query_start", "query_end",
                              "map_match", "map_length",
                              "map_quality")
    }

    #Seqinfo
    unikChr <- paf %>%
      dplyr::select(.data$target_name, .data$target_length) %>%
      dplyr::distinct()
    myseqinfo <- GenomeInfoDb::Seqinfo(
      seqnames = unikChr %>% dplyr::pull(.data$target_name),
      seqlengths = unikChr %>% dplyr::pull(.data$target_length)
      )

    #GRanges:
    resgr <- GRanges(seqnames = paf %>% dplyr::pull(.data$target_name),
                     ranges = IRanges(
                       start = paf %>% dplyr::pull(.data$target_start),
                       end = paf %>% dplyr::pull(.data$target_end)),
                     strand = paf %>% dplyr::pull(.data$strand),
                     seqinfo = myseqinfo)
    S4Vectors::mcols(resgr) <- metacols
  } else {
    #Metadata:
    if (includetags && ncol(paf)>12) {
      metacols <- S4Vectors::DataFrame(paf[,c(6:ncol(paf))])
    } else {
      metacols <- S4Vectors::DataFrame(paf[,6:12])
    }

    #Seqinfo
    unikReads <- paf %>%
      dplyr::select(.data$query_name, .data$query_length) %>%
      dplyr::distinct()
    myseqinfo <- GenomeInfoDb::Seqinfo(
      seqnames = unikReads %>% dplyr::pull(.data$query_name),
      seqlengths = unikReads %>% dplyr::pull(.data$query_length)
      )

    #GRanges:
    resgr <- GRanges(seqnames = paf %>% dplyr::pull(.data$query_name),
                     ranges = IRanges(
                       start = paf %>% dplyr::pull(.data$query_start),
                       end = paf %>% dplyr::pull(.data$query_end)),
                     strand = paf %>% dplyr::pull(.data$strand),
                     seqinfo = myseqinfo)
    S4Vectors::mcols(resgr) <- metacols
  }
  return(resgr)

}

