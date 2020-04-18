#' Identify VDV reads based on alignment of the vector on the reads
#'
#' @description
#'   VDV reads contain vector-DNA-vector sequences.
#'   Only reads longer than minReadLength (>=2kb) can be annotated as VDV reads.
#'   VDV reads contain an alignment of the vector in their first 1kb and in their last 1kb
#'   but do not contain any alignment to the vector in the center region (i.e. excluding 1.15x the length of the vector on each side)
#'
#' @param alignGR a \code{GRanges} object containing the alignment of the vector on the reads
#' @param vectorLength integer Length of the vector (in bp)
#' @param minReadLength integer. Minimum read length to be a VDV read (defaults to 2kb)
#' @param SideWidth integer. Length on each side of the read that must align to some extent to the vector
#'
#' @importFrom GenomeInfoDb seqinfo
#' @importFrom methods as
#' @importFrom IRanges overlapsAny
#' @importFrom GenomicRanges resize narrow width
#'
#' @return a character vector of names of the VDV reads
#' @export
#'
#' @examples
#' ## Create a GRanges. Only Read1 and Read2 are VDV reads
#'   rgr <- GenomicRanges::GRanges(c("Read1:1-2000", "Read1:98001-1e5",
#'                                   "Read2:100-1800", "Read2:99e3-1e5",
#'                                   "Read3:1e4-1.4e4"),
#'                                 seqlengths = c("Read1"=1e5, "Read2"=1e5,
#'                                                "Read3"=2e4, "Read4"=5e4))
#' ## Names of VDV reads:
#'   getVDVnames(rgr, 4000)

getVDVnames <- function(alignGR,
                        vectorLength,
                        minReadLength = 2e3,
                        SideWidth = 1e3) {
  vectorLength <- as.integer(vectorLength)
  ReadsGR <- as(GenomeInfoDb::seqinfo(alignGR), "GRanges")
    # Only reads longer than 2kb can be annotated as VDV reads
  minReadLength <- as.integer(minReadLength)
  minReadLength <- max(c(minReadLength, 2000))
  isLengthOK <- width(ReadsGR) > minReadLength
  # The first 1kb of the read must align with the vector
  isStartOK <- IRanges::overlapsAny(
    GenomicRanges::resize(ReadsGR, width = SideWidth, fix = "start"),
    alignGR)
  # The last 1kb of the read must align with the vector
  isEndOK <- IRanges::overlapsAny(
    GenomicRanges::resize(ReadsGR, width = SideWidth, fix = "end"),
    alignGR)
  # The center of the read must not align with the vector
  isCenterOK <- !IRanges::overlapsAny(
    GenomicRanges::narrow(ReadsGR,
                          start = pmin(ceiling(GenomicRanges::width(ReadsGR)/2),
                                       floor(1.15*vectorLength)),
                          end = pmax(floor(width(ReadsGR)/2)+1,
                                     GenomicRanges::width(ReadsGR) -
                                       floor(1.15*vectorLength))),
                             alignGR)
  vdvreads <- names(ReadsGR[isLengthOK & isStartOK & isEndOK & isCenterOK])
  return(vdvreads)
}
