#' Select DVD reads from a set of Nanopore reads and split the reads in 2 parts both containing the vector (DV and VD)
#'
#' The function does the following:
#' \itemize{
#'   \item{Selects DVD reads}{This is done using \code{\link{FilterBACreads}}}
#'   \item{Split the read sequence in DV and VD}{Based on vector alignment, split the read sequence in DV and VD}
#'   \item{Filter based on size}{Keep only the split reads with a DNA fragment that is at least \code{MinDNAlength}bp long}
#'   \item{reverse complement reads on minus strand}{strand is determined based on vector alignment}
#'   \item{Return the split reads}{The split reads are returned as a DNAString object}
#'   }
#' By default, if alignemnt to the host genome is provided in the \code{ReadClass} object (column \code{HostAlign}),
#' then the selected DVD reads are selected to not show any significant alignment to the host genome
#'
#' @param ReadClass Either a tibble obtained with the \code{\link{AnnotateBACreads}} function
#'                  or a path to an rds file containing such a file
#' @param blastvec Either a table imported with \code{\link{readBlast}} or a path to a blast file
#'                 obtained by aligning th evectors on the reads and using \code{-outfmt 6}
#' @param FastaFile Either a \code{DNAStringSet} object containing the full read sequences
#'                  or a path to a fasta file containing these sequences
#' @param WithGeneA Logical. Should the VDV reads align with GeneA? Default is NULL, i.e. no filtering on GeneA alignment
#' @param WithGeneB Logical. Should the VDV reads align with GeneB? Default is NULL, i.e. no filtering on GeneB alignment
#' @param MinDNAlength Integer. Minimum length of the DNA fragment to keep the reads in the results
#'
#' @return A list with:
#'     \itemize{
#'         \item{ReadDefinition}{ a \code{DNAStringSet} with the split reads}
#'         \item{ReadSequence}{ a \code{GRanges} object with the definition of the DV/VD reads}
#'     }
#'     Note that reads with alignment on the opposite strand of the vector ("-" strand)
#'     are automatically reverse complemented
#'        If no reads are selected, the function returns NULL and a warning.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter left_join mutate transmute select rename
#' @importFrom rlang .data
#' @importFrom Rsamtools FaFile indexFa scanFa
#' @importFrom methods is
#' @importFrom Biostrings readDNAStringSet
#' @importFrom GenomicRanges GRanges mcols width
#' @importFrom BSgenome getSeq
#' @importClassesFrom Biostrings DNAString DNAStringSet
#'
#' @author Pascal GP Martin
#'
#' @export
#'
#' @examples
#' ## For simplicity (and to limit file size) we only keep the data for 5 pre-selected DVD reads
#' ## Path to file (.rds) created with the AnnotateBACreads function
#' pathRC <- system.file("extdata", "BAC02_ReadClass.rds", package = "NanoBAC")
#' RC <- readRDS(pathRC)
#' selectedReads <- c("BAC02R5572", "BAC02R21438", "BAC02R1152",
#'                    "BAC02R20794", "BAC02R6278" )
#' RC <- RC[RC$ReadName %in% selectedReads,]
#' ## Path to a fasta file containing the sequence of the 5 DVD reads
#' pathFasta <- system.file("extdata", "BAC02_5DVDreads.fa", package = "NanoBAC")
#' ## Path to the file containing the result from the Blast alignment of the vector on the reads
#' pathBlast <- system.file("extdata", "BAC02_BlastVector.res", package = "NanoBAC")
#' ## Select DVD reads and split the reads
#' myDVDreads <- splitDVDreads(ReadClass = RC,
#'                             blastvec = pathBlast,
#'                             FastaFile = pathFasta,
#'                             WithGeneA = TRUE,
#'                             WithGeneB = TRUE,
#'                             MinDNAlength = 35000)
#' ## Read sequences:
#' myDVDreads$ReadSequence
#' ## Read definitions:
#' myDVDreads$ReadDefinition

splitDVDreads <- function(
                          ReadClass = NULL,
                          blastvec = NULL,
                          FastaFile = NULL,
                          WithGeneA = NULL,
                          WithGeneB = NULL,
                          MinDNAlength = 1e4L
                          ) {

#---------------------
# Required Arguments
#---------------------

## ReadClass argument
  if (is.null(ReadClass)) {
    stop(strwrap(prefix = " ", initial = "",
                 "Please provide a ReadClass object.
                 either a tibble obtained using the AnnotateBACreadsfunction
                 or a path to an RDS file containing such a tibble"))
  }


## Blast results
  if (is.null(blastvec)) {
    stop(strwrap(prefix = " ", initial = "",
                 "Please provide a blastvec object
                 (result of blast alignment of the vector on the reads)"))
  }


## Fastafile
  if (is.null(FastaFile)) {
    stop(strwrap(prefix = " ", initial = "",
                 "Please provide a FastaFile object.
                  either a DNAString objet with read sequences
                  or a path to a fasta file"))
  }

  if (!is(FastaFile, "DNAStringSet")) {
    if (!file.exists(FastaFile)) {
      stop(FastaFile, " not found")
    }
  }

#---------------------
# Import ReadClass and check other arguments
#---------------------

## import ReadClass
  RC <- TestArg(ReadClass, readRDS, "tbl_df")

  expectedColNames <- c("ReadName", "ReadLength", "Strand",
                        "DNAlength", "LongestDNA", "ShortestDNA", "ReadType")

  if (ncol(RC)<7 || !identical(colnames(RC)[1:7], expectedColNames)) {
    stop("ReadClass does not have the expected first 7 columns")
  }

## WithGeneA
  if (is.na(WithGeneA) || is.null(WithGeneA)) {WithGeneA <- NULL}
  if (!is.null(WithGeneA) && !is.logical(WithGeneA)) {
    stop("WithGeneA should be TRUE/FALSE (or NULL)")
  }
  if (!is.null(WithGeneA) && !("AlignedToGeneA" %in% colnames(RC))) {
    warning(strwrap(prefix = " ", initial = "",
                    "No AlignedToGeneA column found in ReadClass table.
                    No filtering done on GeneA alignement"))
  }

## WithGeneB
  if (is.na(WithGeneB) || is.null(WithGeneB)) {WithGeneB <- NULL}
  if (!is.null(WithGeneB) && !is.logical(WithGeneB)) {
    stop("WithGeneB should be TRUE/FALSE (or NULL)")
  }
  if (!is.null(WithGeneB) && !("AlignedToGeneB" %in% colnames(RC))) {
    warning(strwrap(prefix = " ", initial = "",
                    "No AlignedToGeneB column found in ReadClass table.
                    No filtering done on GeneB alignement"))
  }


#---------------------
# Select DVD reads
#---------------------
  if ("AlignedToGeneA" %in% colnames(RC)) {
    GnAfilter = WithGeneA
  } else {
    GnAfilter = NULL
  }

  if ("AlignedToGeneB" %in% colnames(RC)) {
    GnBfilter = WithGeneB
  } else {
    GnBfilter = NULL
  }

  if ("HostAlign" %in% colnames(RC)) {
    Hostfilter = FALSE #do not select reads that align to the host genome
  } else {
    Hostfilter = NULL
  }

#Select DVD reads
  allDVD <- FilterBACreads(
                           RC,
                           readtype = "DVD",
                           alnGeneA = GnAfilter,
                           alnGeneB = GnBfilter,
                           isHostAlign = Hostfilter
                           )



#---------------------
# If DVD reads are present split their sequences in 2 reads containing the vector sequence
#---------------------

  if (is.null(allDVD)) {

    warning("No DVD reads in the dataset")
    res <- NULL

  } else {
    #keep only reads with a DNA fragment longer than MinDNAlength
    allDVD <- dplyr::filter(allDVD, .data$LongestDNA >= MinDNAlength)

    if (nrow(allDVD) == 0) {
      warning("No DVD reads selected after filtering for DNA size")
      res <- NULL

    } else {

    ## import blastvec, filter and keep only alignments for DVD reads
    vecaln <- checkBlastVar(blastvec)
    vecaln <- dplyr::filter(vecaln, .data$SubjectACC %in% allDVD$ReadName)
    vecaln <- vecaln[SelectSingularBlastALN(vecaln,
                                            allDVD[,c("ReadName", "ReadLength")],
                                            threshold = 0.5),]

    ## Create a GRanges with the regions to extract
    #Define the reads to extract
    readDef <- vecaln %>%
      dplyr::left_join(allDVD %>%
                         dplyr::select(.data$ReadName, .data$ReadLength) %>%
                         dplyr::rename("SubjectACC" = "ReadName"),
                       by = "SubjectACC") %>%
      dplyr::mutate(vecStart = pmin(.data$SubjectStart, .data$SubjectEnd),
                    vecEnd = pmax(.data$SubjectStart, .data$SubjectEnd)) %>%
      dplyr::transmute(Read1 = paste0(.data$SubjectACC, ":",
                                   1, "-", .data$vecEnd,
                                   ":", .data$Strand),
                       Read2 = paste0(.data$SubjectACC, ":",
                                      .data$vecStart, "-", .data$ReadLength,
                                      ":", .data$Strand),
                       vecLength = .data$vecEnd - .data$vecStart + 1,
                       ReadName = .data$SubjectACC)
    #Build a GRanges
    readGR <- GenomicRanges::GRanges(c(readDef$Read1,
                                       readDef$Read2))
    names(readGR) <- paste(rep(readDef$ReadName, 2),
                           rep(c("R1", "R2"), each = nrow(readDef)),
                           sep="_")
    GenomicRanges::mcols(readGR)$DNAlength <- GenomicRanges::width(readGR) -
      rep(readDef$vecLength, 2)
    readGR <- GenomicRanges::sort(readGR)

    # Remove reads with a length of DNA (vector excluded) shorter than MinDNAlength
    readGR <- readGR[GenomicRanges::mcols(readGR)$DNAlength >= MinDNAlength]

    ## Import and split sequences
    if (!is(FastaFile, "DNAStringSet")) {

      if (!file.exists(paste0(FastaFile, ".fai"))) {
        Rsamtools::indexFa(FastaFile)
      }

      hemiDVDreads <- Rsamtools::scanFa(Rsamtools::FaFile(FastaFile),
                                        readGR)
      names(hemiDVDreads) <- names(readGR)
      isMinusStrand <- GenomicRanges::strand(readGR)=="-"
      hemiDVDreads[isMinusStrand] <-
        Biostrings::reverseComplement(hemiDVDreads[isMinusStrand])

    } else {

      hemiDVDreads <- BSgenome::getSeq(FastaFile, readGR)
      names(hemiDVDreads) <- names(readGR)

    }

    res <- list(ReadDefinition = readGR,
                ReadSequence = hemiDVDreads)

    }
  }

return(res)

}

