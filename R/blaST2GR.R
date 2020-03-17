#' Convert an imported blast result table to a GRanges
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb Seqinfo
#' @importFrom magrittr %>%
#' @importFrom dplyr pull
#' @importFrom rlang .data
#' @param blastresults A tibble obtained with the readBlast function (containing strand info)
#' @param readlength A tibble with 2 columns named ReadName and ReadLength.
#'                   This table is expected to contain all reads, not only those having an alignment to the vector
#' @export
#' @return GRanges
#' @examples
#' ## Import some blast result data:
#'     Path2Blast <- system.file("extdata", "blastExample.tab", package = "NanoBAC")
#'     myblastResult <- NanoBAC::readBlast(Path2Blast)
#' ## Import a 2-column table with ReadName and ReadLength:
#'     Path2ReadLength <- system.file("extdata", "ReadLengthExample.tsv", package = "NanoBAC")
#'     ReadLengthTable <- read.table(Path2ReadLength,
#'                                   sep = "\t", header = FALSE,
#'                                   stringsAsFactors = FALSE)
#'     colnames(ReadLengthTable) <- c("ReadName", "ReadLength")
#' ## Convert the blast result table to a GRanges
#'     myBlastGR <- blaST2GR(myblastResult, ReadLengthTable)


blaST2GR <- function(blastresults, readlength) {
# Note that readLength should contain all reads (without filtering compared to the reads used to align the vector)
        resgr <- GenomicRanges::GRanges(seqnames = as.character(blastresults$SubjectACC),
                         ranges = IRanges(start = ifelse(blastresults$Strand == "-",
                                                         blastresults$SubjectEnd,
                                                         blastresults$SubjectStart),
                                         end = ifelse(blastresults$Strand == "-",
                                                      blastresults$SubjectStart,
                                                      blastresults$SubjectEnd)),
                         strand = blastresults$Strand,
                         "PercentID" = blastresults$PercentID, "AlnLength" = blastresults$AlnLength,
                         "NumMismatch" = blastresults$NumMismatch, "NumGapOpen" = blastresults$NumGapOpen,
                         "QueryRange" = IRanges(start = blastresults$QueryStart, end = blastresults$QueryEnd),
                         "evalue" = blastresults$evalue, "bitscore" = blastresults$bitscore,
                         seqinfo = GenomeInfoDb::Seqinfo(seqnames = readlength %>% dplyr::pull(.data$ReadName),
                                                         seqlengths = readlength %>% dplyr::pull(.data$ReadLength),
                                                         isCircular= rep(FALSE, nrow(readlength))))
        return(resgr)
        }

