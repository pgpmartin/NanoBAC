#' Filter a blast result table for alignment overlapping other alignments
#'
#' @description
#'     Filter a blast result table imported with \code{\link[NanoBAC]{readBlast}}
#'     to remove alignments with large overlaps (e.g. > 50\%) with other alignments
#'
#' @importFrom GenomicRanges findOverlaps pintersect width setdiff
#' @importFrom S4Vectors queryHits subjectHits
#'
#' @param aln A tibble obtained with the \code{\link[NanoBAC]{readBlast}} function (containing strand info)
#' @param rl A ReadLength table, i.e. a data frame with 2 columns: ReadName and ReadLength
#' @param threshold Number in ]0,1]. Remove all alignments that overlap with other alignments on > threshold \% of their length
#'
#' @export
#' @return a character vector with the row numbers corresponding to alignments to keep
#'     (i.e. that don't overlap on >threshold% of their length with other alignments)
#'
#' @seealso \link{readBlast}
#' @references
#'   The code for this function is based on a suggestion by Michael Lawrence: \url{https://support.bioconductor.org/p/72656/}
#'
#' @examples
#' ## Example datasets
#'   Path2OVL <- system.file("extdata", "blastExampleOverlap.tab", package = "NanoBAC")
#'   alignmt <- readBlast(Path2OVL)
#' ## Read lengths
#'   rl <- data.frame(ReadName = c("BAC43R75552", "BAC43R80634"),
#'                    ReadLength = c(11e4, 1e4))
#' ## Filter the data to keep alignments not overlapping on more than 50% of their length
#'   filtaln <- alignmt[SelectSingularBlastALN(alignmt, rl),]
#'

SelectSingularBlastALN <- function(aln,
                                   rl,
                                   threshold = 0.5) {
  ## Test threshold argument
    stopifnot(threshold > 0 && threshold <= 1)

  ## Convert to GRanges
  alngr <- NanoBAC::blaST2GR(aln, readlength = rl)

  ## Identify the overlaps between the alignments
  hits <- GenomicRanges::findOverlaps(alngr, drop.self = TRUE)
  # Obtain the overlaps between the alignments
  overlaps <- GenomicRanges::pintersect(alngr[S4Vectors::queryHits(hits)],
                                        alngr[S4Vectors::subjectHits(hits)])
  # Calculate the % of overlap
  percentOverlap <- GenomicRanges::width(overlaps) /
    GenomicRanges::width(alngr[S4Vectors::subjectHits(hits)])

  # Return the index of the alignments that do not have any overlap of >threshold% of their length with other alignments
  AlnToKeep <-
    setdiff(1:length(alngr),
            unique(S4Vectors::subjectHits(hits[percentOverlap > threshold])))

  return(AlnToKeep)
}

