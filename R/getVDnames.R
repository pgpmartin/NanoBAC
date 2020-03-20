#' Identify VD reads based on alignment of the vector on the reads
#' @description
#'   VD reads contain DNA-vector or vector-DNA
#'   These reads contain a single region with an alignment to the vector
#'   The vector alignment must cover at least tolerance\% of the first (or last) EndWindow bp at the read extremity
#'
#' @param alignGR a \code{GRanges} object containing the alignment of the vector on the reads
#' @param tolerance numeric in [0.5,1] indicating the pecentage of the read extremity that must be covered by vector alignment to be considered a VD read
#' @param EndWindow integer. Number of bp at the extremity of the read to look for vector alignment
#'
#' @importFrom GenomeInfoDb seqinfo seqnames
#' @importFrom methods as
#' @importFrom tibble tibble as_tibble
#' @importFrom GenomicRanges resize width findOverlaps pintersect
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom magrittr %>%
#' @importFrom dplyr pull count mutate summarize filter left_join
#' @importFrom rlang .data
#'
#' @return a character vector with the names of VD reads
#' @export
#'
#' @examples
#' ## Create a GRanges. only Read1 and Read2 are VD reads
#' ## vector is 10kb
#' ## alignment of vector on Read2 covers 410bp in the first 500bp
#' ## (i.e. >80\% but <90\% of the first 500bp)
#'   rgr <- GenomicRanges::GRanges(c("Read1:48e3-5e4",
#'                                   "Read2:90-500",
#'                                   "Read3:1-2000", "Read3:98001-1e5"),
#'                                 seqlengths = c("Read1"=5e4, "Read2"=6e4,
#'                                                "Read3"=1e5, "Read4"=5e4))
#' ## Names of VD reads (using default 80\% of the 500bp at the extremity of the read):
#'   getVDnames(rgr, 0.8, 500)
#' ## With 90\% of vector length only Read1 will be selected as a VD read
#'   getVDnames(rgr, 0.9, 500)

getVDnames <- function(alignGR,
                       tolerance = 0.8,
                       EndWindow = 500) {

  if (tolerance > 1 | tolerance < 0.5) {
    stop("tolerance should be in [0.5,1]")
  }


  #Get a GRanges of all reads
  readGR <- as(GenomeInfoDb::seqinfo(alignGR), "GRanges")
  readNames <- tibble::tibble(seqnames = as.character(GenomeInfoDb::seqnames(readGR)))

  # Is there at least one alignment?
  hasOVL <- (readNames %>% dplyr::pull(seqnames)) %in%
    unique(as.character(GenomeInfoDb::seqnames(alignGR)))

  # Does the read have a unique alignment?
  # We cold improve this by imposing that all alignements be located in a defined window (e.g. vector size +/- 5-10%)
  isUniqueAlignment <- dplyr::left_join(readNames,
                                        tibble::as_tibble(alignGR) %>%
                                          dplyr::count(.data$seqnames) %>%
                                          dplyr::mutate(isUnique = (.data$n==1)) %>%
                                          dplyr::mutate(seqnames = as.character(.data$seqnames)),
                                        by = "seqnames") %>%
    dplyr::pull(.data$isUnique)

  # Does the alignment cover at least e.g. 400bp within the first or last 500 bp (i.e. 80% of the 500bp at read end is covered by vector sequence) ?
  #Read start: coverage of alignment(s) in EndWindow
  ReadStart <- GenomicRanges::resize(readGR, width = EndWindow, fix = "start")
  fovStart <- GenomicRanges::findOverlaps(ReadStart, alignGR)
  interStart <- GenomicRanges::pintersect(ReadStart[S4Vectors::queryHits(fovStart)],
                                          alignGR[S4Vectors::subjectHits(fovStart)])
  CovStart <- sum(GenomicRanges::coverage(GenomicRanges::reduce(interStart))) #Number of bases aligned at the beginning of the read
  #Read end: coverage of alignment(s) in EndWindow
  ReadEnd <- GenomicRanges::resize(readGR, width = EndWindow, fix = "end")
  fovEnd <- GenomicRanges::findOverlaps(ReadEnd, alignGR)
  interEnd <- GenomicRanges::pintersect(ReadEnd[S4Vectors::queryHits(fovEnd)],
                                        alignGR[S4Vectors::subjectHits(fovEnd)])
  CovEnd <- sum(GenomicRanges::coverage(GenomicRanges::reduce(interEnd))) #Number of bases aligned at the end of the read

  stopifnot(identical(readNames %>% dplyr::pull(.data$seqnames),
                      names(CovStart))) # should be guaranteed because seqinfo is passed to coverage function

  isAlnLengthOK <- (CovStart >= ceiling(tolerance * EndWindow)) |
    (CovEnd >= ceiling(tolerance * EndWindow))

  return(readNames %>%
           dplyr::filter(hasOVL & isUniqueAlignment & isAlnLengthOK) %>%
           dplyr::pull(.data$seqnames))

}
