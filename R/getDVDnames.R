#' Identify DVD reads based on alignment of the vector on the reads
#'
#' @description
#'   DVD reads contain DNA-vector-DNA
#'   DVD reads have a single vector alignment covering at least PercentVecLength\% (default 95\%) of the vector length.
#    They also have 2 pieces of insert (non-vector DNA) on each side that are at least MinDNASides bp long (defaults to 10kb)
#'
#' @param alignGR a \code{GRanges} object containing the alignment of the vector on the reads
#' @param vectorLength integer Length of the vector (in bp)
#' @param PercentVecLength numeric in [0.5,1] indicating the \% of the vector length that must be aligned on the read to be considered a DVD read
#' @param MinDNASides minimum length of non-vector DNA on each side of the central vector alignment
#'
#' @importFrom GenomeInfoDb seqinfo seqnames
#' @importFrom methods as
#' @importFrom IRanges overlapsAny
#' @importFrom tibble tibble as_tibble
#' @importFrom GenomicRanges resize width
#' @importFrom magrittr %>%
#' @importFrom dplyr pull count mutate group_by summarize filter left_join
#' @importFrom rlang .data
#'
#' @return a character vector with the names of DVD reads
#' @export
#'
#' @examples
#' ## Create a GRanges. only Read1 and Read2 are DVD reads
#' ## vector is 10kb
#' ## alignment of vector on Read2 covers 9501 bp (>95\% of vector length)
#'   rgr <- GenomicRanges::GRanges(c("Read1:12e3-22e3",
#'                                   "Read2:20e3-29.6e3",
#'                                   "Read3:1-2000", "Read3:98001-1e5"),
#'                                 seqlengths = c("Read1"=5e4, "Read2"=6e4,
#'                                                "Read3"=1e5, "Read4"=5e4),
#'                                  QueryRange.width = c(1e4, 9501, 2000, 2000))
#' ## Names of VDV reads (using 95\% of vector length as threshold):
#'   getDVDnames(rgr, 10000, 0.95)
#' ## With 98\% of vector length only Read1 will be selected as a DVD read
#'   getDVDnames(rgr, 10000, 0.98)
#' ## If requiring at least 15kb on each side of the vector then only Read2 is DVD
#'   getDVDnames(rgr, 10000, 0.95, 15e3)

getDVDnames <- function(alignGR,
                        vectorLength,
                        PercentVecLength = 0.95,
                        MinDNASides = 10e3L) {

  stopifnot(is.numeric(vectorLength))

  if (PercentVecLength > 1 | PercentVecLength < 0.5) {
    stop("PercentVecLength should be in [0.5,1]")
  }

  #Get a GRanges of all reads
  readGR <- as(GenomeInfoDb::seqinfo(alignGR), "GRanges")
  readNames <- tibble::tibble(seqnames = as.character(GenomeInfoDb::seqnames(readGR)))

  # Is there at least one alignment
  hasOVL <- (readNames %>%
               dplyr::pull(.data$seqnames)) %in%
    unique(as.character(GenomeInfoDb::seqnames(alignGR)))

  # Does the read have a unique alignment?
  isUniqueAlignment <- dplyr::left_join(readNames,
                                        tibble::as_tibble(alignGR) %>%
                                          dplyr::count(.data$seqnames) %>%
                                          dplyr::mutate(isUnique = (.data$n==1)) %>%
                                          dplyr::mutate(seqnames = as.character(.data$seqnames)),
                                        by = "seqnames") %>%
    dplyr::pull(.data$isUnique)

  # Is the alignment length long enough (>95\% of vector length)?
  isAlnLengthOK <- dplyr::left_join(readNames,
                                    tibble::as_tibble(alignGR) %>%
                                      dplyr::group_by(.data$seqnames) %>%
                                      dplyr::summarize(maxAlnLength = max(.data$QueryRange.width)) %>%
                                      dplyr::mutate(isAlnLongEnough = (.data$maxAlnLength >= ceiling(PercentVecLength*vectorLength))) %>%
                                      dplyr::mutate(seqnames = as.character(.data$seqnames)),
                                    by = "seqnames") %>%
    dplyr::pull(.data$isAlnLongEnough)

  # Is there no alignment in the first MinDNASides bp?
  isStartOK <- !IRanges::overlapsAny(GenomicRanges::resize(readGR,
                                                           width = MinDNASides,
                                                           fix="start"),
                                     alignGR,
                                     ignore.strand = TRUE)
  # Is there no alignment in the last MinDNASides bp
  isEndOK <- !IRanges::overlapsAny(GenomicRanges::resize(readGR,
                                                         width = MinDNASides,
                                                         fix="end"),
                                   alignGR,
                                   ignore.strand = TRUE)
  ## Return names of reads with a single alignment in [MinDNASides:end-MinDNASides] that covers >PercentVecLength\% of the vector
  return(readNames %>%
           dplyr::filter(hasOVL &
                           isUniqueAlignment &
                           isAlnLengthOK &
                           isStartOK &
                           isEndOK) %>%
           dplyr::pull(.data$seqnames))
}
