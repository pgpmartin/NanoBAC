#' Find the vector-insert junctions at the beginning and end of a VDV read
#'
#' @param ReadName character string. Name of the read
#' @param ReadDNA A DNAString or DNAStringSet with the read sequence
#' @param ReadVecAlign Table with the Blast results from aligning th evector on the read
#' @param RestrictionSite Character string in the for "G^AATTC" indicating the sequence and the cut site
#' @param VectorSequence A DNAString or DNAStringSet of length 1 with the vector sequence starting and ending
#'                       with the full sequence of the restriction site used for cloning
#' @param UnalignedVectorLength Integer. If more than \code{UnalignedVectorLength}bp of expected vector sequence
#'                                       is not correctly aligned at the vector-insert junction,
#'                                       the function will return a message
#' @param SideSeqSearch Integer. If the expected restriction site is not found at a vector-insert junction
#'                               then the algorithm will try to search for the vector sequence of length
#'                               \code{SideSeqSearch} bp that is adjacent to the restriction site
#'                               (Default to 10 bp)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter top_n mutate pull
#' @importFrom rlang .data
#' @importFrom Biostrings subseq matchPattern reverseComplement
#' @importFrom methods is as
#' @importFrom BiocGenerics start end width
#' @importClassesFrom Biostrings DNAString DNAStringSet
#'
#' @export
#'
#' @return a named character vector with the read name, the read strand and the coordinates of the first and last bases of the insert
#'
#' @examples
#' # Some dummy sequences (use paste0 for clarity / comparison of sequences):
#' ## Vector sequence start and ends with the HindIII site used for cloning ("A^AGCTT")
#' vector <- Biostrings::DNAString(paste0(
#'             "AAGCTTTATTAAGACACCCGGTATGCTTCAGGATCGTTCGGACTAA",
#'             "ACCGTAACTGCGATATTTTAGGCGTGTTACAAGCTT"))
#' read <- Biostrings::DNAString(paste0(
#'             "ACCGTAACTGCGATATTTTAGGCGTGTTACAAGCTT",
#'             "GCTAGATCGCGCGATATGTG",
#'             "AAGCTTTATTAAGACACCCGGTATGCTTCAGGATCGTTCGGACTAA"))
#' noisyread <- Biostrings::DNAString(paste0(
#'             "ACCGTAACTGCGTTTTTTTAGGCGTGTTACAAGCTT",
#'             "GCTAAATCGCGCGCTATGTG",
#'             "GGGCTTTATTAAGACACCCGGTATGCTTTCAGGATCGTTCGGACTAA"))
#' # Import the blastn results for these sequences (alignment of vector on the reads):
#' readaln <- readBlast(system.file("extdata",
#'                                  "juncEx_vec_read.res",
#'                                  package = "NanoBAC"))
#' noisyreadaln <- readBlast(system.file("extdata",
#'                                       "juncEx_vec_noisyread.res",
#'                                       package = "NanoBAC"))
#' # Get the coordinates of the insert sequence:
#' FindVDVjunctions("read", read, readaln, "A^AGCTT", vector)
#' # With the noisy read, the restriction site is not found but an adjacent sequence is:
#' FindVDVjunctions("noisyread", noisyread, noisyreadaln, "A^AGCTT", vector)

FindVDVjunctions <- function(ReadName = NULL,
                             ReadDNA,
                             ReadVecAlign,
                             RestrictionSite = "G^AATTC",
                             VectorSequence = NULL,
                             UnalignedVectorLength = 1000L,
                             SideSeqSearch = 10L) {


  #---------------------
  #Test arguments and import data if necessary
  #---------------------

  ## ReadName
  if (is.null(ReadName)) {
    stop("Please provide a read name in the ReadName argument (no default)")
  } else {
    rnn <- as.character(ReadName)
  }

  ## ReadDNA
  if (!is(ReadDNA, "DNAStringSet") && !is(ReadDNA, "DNAString")) {
    stop("ReadDNA should be a DNAString or a DNAStringSet")
  }

  ## Make ReadDNA a DNAStringSet of length 1 with names(ReadDNA)==ReadName
  if (is(ReadDNA, "DNAString")) {
    ReadDNA <- as(ReadDNA, "DNAStringSet")
    names(ReadDNA) <- rnn
  } else {
    if (!(rnn %in% names(ReadDNA))) {
        stop("Read name not found in ReadDNA")
    } else {
        ReadDNA <- ReadDNA[rnn]
    }
  }

  ## Filter read/vector alignments. Throw a message if there are >2 alignments
  ## If this message appears for many reads, then consider filtering Blast results with SelectSingularBlastALN function
  if (!is.data.frame(ReadVecAlign)) {
    stop("ReadVecAlign is not a data frame")
  }

  if (!(rnn %in% (ReadVecAlign %>%
                  dplyr::pull(.data$SubjectACC)))) {
    stop("Read name not found in ReadVecAlign")
  } else {
    ReadVecAlign <- ReadVecAlign %>%
      dplyr::filter(.data$SubjectACC == rnn)
    if (nrow(ReadVecAlign) > 2) {
        message(rnn, ": found ",
                nrow(ReadVecAlign),
                " alignments. Should ReadVecAlign be filtered?")
    }
  }

  ## Vector sequence and restriction site
  ### The vector sequence is expected to start by the restriction site
  if (is.null(VectorSequence)) {
    stop(strwrap(prefix = " ", initial = "",
                 "Please provide the vector sequence, starting with
                 the restriction site used for cloning"))
  }

  if (is.character(VectorSequence)) {
    VectorSequence <- as(VectorSequence, "DNAString")
  }
  if (is(VectorSequence, "DNAStringSet")) {
    VectorSequence <- VectorSequence[[1]]
  }

  ### Get the vector size
  VectorSize <- length(VectorSequence)


  ### Check that vector sequence starts and end with restriction site

  CutSiteLocation <- which(strsplit(RestrictionSite, "")[[1]] == "^")
  if (!length(CutSiteLocation) == 1) {
    stop(strwrap(prefix = " ", initial = "",
                 "Please provide RestrictionSite in the following format
                 (ex for EcoRI): G^AATTC"))
  }
  RestrictionSiteLength <- nchar(RestrictionSite) - 1
  BasesBeforeCut <- CutSiteLocation - 1
  BasesAfterCut <- RestrictionSiteLength - BasesBeforeCut
  RestrictionSiteSequence <- as(gsub("\\^", "", RestrictionSite), "DNAString")

  if (!identical(toString(RestrictionSiteSequence),
                 toString(Biostrings::subseq(VectorSequence,
                       1,
                       length(RestrictionSiteSequence))))) {
    stop("Vector sequence does not start with restriction site")
  }
  if (!identical(toString(RestrictionSiteSequence),
                 toString(Biostrings::subseq(VectorSequence,
                                             VectorSize-RestrictionSiteLength+1,
                                             VectorSize)))) {
    stop("Vector sequence does not end with restriction site")
  }


  # Get the strand of the read
    readStrand <- unique(ReadVecAlign$Strand)
    if (length(readStrand)!=1) {
        stop(rnn, " has alignment on both strands. Read cannot be analyzed")
    }


  #---------------------
  #Initial search for potential vector-insert borders at the beginning and end of the reads
  #---------------------

  ## Detect if the aligments are at the beginning or at the end of the read.
  ## Add this info to the alignment table ReadVecAlign
    ReadVecAlign <- ReadVecAlign %>%
      dplyr::mutate(ReadSide = ifelse(pmax(.data$SubjectStart,
                                           .data$SubjectEnd) <
                                        length(ReadDNA[[1]])/2,
                                      "Begin", "End"))


  ## Identify the theoretical location of the beginning/end of the insert sequence (based on alignment positions only):
    ### For reads with alignment on the "+" strand:
      if (readStrand == "+") {
### 1) find the largest plasmid coordinate that aligns at the beginning of the read
###    then estimate where the vector-insert border should be if the alignment went all the way to the end of the vector (InsertBorder)
###    and obtain the distance between the end of the vector and where the alignment stops (DistToVectorEnd)
          beginAlign <- ReadVecAlign %>%
                          dplyr::filter(.data$ReadSide == "Begin") %>%
                          dplyr::top_n(1, .data$QueryEnd) %>%
                          dplyr::mutate(InsertBorder = .data$SubjectEnd + 1 +
                                          VectorSize - .data$QueryEnd,
                                        DistToVectorEnd = (VectorSize -
                                                             .data$QueryEnd))
### 2) find the smallest plasmid coordinate that aligns at the end of the read
###    then estimate where the insert-vector border should be if the alignment started at the first base of the vector (InsertBorder)
###    and obtain the distance between the start of the vector (1) and where the alignment starts (DistToVectorEnd)
          endAlign <- ReadVecAlign %>%
                          dplyr::filter(.data$ReadSide == "End") %>%
                          dplyr::top_n(1, -.data$QueryStart) %>%
                          dplyr::mutate(InsertBorder = .data$SubjectStart - 1 -
                                          .data$QueryStart + 1,
                                        DistToVectorEnd = .data$QueryStart - 1)
     } else {
     ## For reads with alignment on the "-" strand:
### 1) find the smallest vector coordinate that aligns at the beginning of the read
          beginAlign <- ReadVecAlign %>%
                          dplyr::filter(.data$ReadSide == "Begin") %>%
                          dplyr::top_n(1,-.data$QueryStart) %>%
                          dplyr::mutate(InsertBorder = .data$SubjectStart + 1 +
                                          .data$QueryStart - 1,
                                        DistToVectorEnd = .data$QueryStart - 1)
### 2) find the largest vector coordinate that aligns at the end of the read
          endAlign <- ReadVecAlign %>%
                          filter(.data$ReadSide == "End") %>%
                          top_n(1, .data$QueryEnd) %>%
                          mutate(InsertBorder = .data$SubjectEnd - 1 -
                                   VectorSize + .data$QueryEnd,
                                 DistToVectorEnd = (VectorSize -
                                                      .data$QueryEnd))
    }


  ## If the alignment does not cover at all the first (resp. last) UnalignedVectorLength bp of the vector sequence, throw a message
    if (dplyr::pull(beginAlign, .data$DistToVectorEnd) >= UnalignedVectorLength) {
        message(rnn,
                ": at the beginning of the read, more than ",
                UnalignedVectorLength,
                "bp of expected vector sequence are not aligned")
    }
    if (dplyr::pull(endAlign, .data$DistToVectorEnd) >= UnalignedVectorLength) {
      message(rnn,
              ": at the end of the read, more than ",
              UnalignedVectorLength,
              "bp of expected vector sequence are not aligned")
    }

# Search for the restriction site at the predicted vector-insert junction (i.e. in a window around InsertBorder) allowing only 1 mismatch/indel
# If found, the location of the restriction site determines the location of the vector-insert junction
# If the restriction sit is not found or found multiple times, we search for the vector sequence (currently 10bp) adjacent to the restriction site, allowing 20% mismatch/indel and using a wider search window.
# If not found then JuncAtBegin wil be NA

  #----------------------------------------
  # Search the vector-insert junction at the beginning of the read:
  #----------------------------------------
    PointOfSearch <- beginAlign %>%
      dplyr::pull(.data$InsertBorder) #This is the position of the first base following the restriction site

    # We search for the restriction site in a window covering +/-15bp or +/- 5% of DistToVectorEnd if larger
    # i.e. the farther we are from vector end, the wider the window:
    WidthOfwindow <- max(c(15,
                           round(0.05*(beginAlign %>%
                                         dplyr::pull(.data$DistToVectorEnd)),
                                 0)))
    # Sequence to search
    SearchStart <- PointOfSearch - RestrictionSiteLength - WidthOfwindow
    SearchEnd <- PointOfSearch + WidthOfwindow - 1
    SeqToSearch <- Biostrings::subseq(ReadDNA[[rnn]],
                                      SearchStart,
                                      SearchEnd) # +/- WidthOfwindow bases on each side of the restriction motif
    # Search for the restriction site
    findRestrSite <- Biostrings::matchPattern(RestrictionSiteSequence,
                                              SeqToSearch,
                                              max.mismatch = 1,
                                              with.indels= TRUE)
    #Define JuncAtBegin which is the first base of the insert
    if (length(findRestrSite) != 1) {
        JuncAtBegin <- NA
        if (length(findRestrSite) == 0) {
            message(rnn,
                    strwrap(prefix = " ", initial = "",
                            ": restriction site not found near the expected
                            position at the beginning of the read"))
        } else {
            message(rnn,
                    strwrap(prefix = " ", initial = "",
                            ": multiple restriction sites found near the
                            expected position at the beginning of the read"))
        }
    } else {
      JuncAtBegin <-  PointOfSearch - RestrictionSiteLength -
                          WidthOfwindow + end(findRestrSite)
    }


    # If the restriction site is not found:
    #   then try to search for the vector sequence located just upstream of the restriction site
    if (is.na(JuncAtBegin)) {
        SideSeqLength = as.integer(SideSeqSearch)
        message(rnn, ": Searching for a ", SideSeqLength,
                strwrap(prefix = " ", initial = "",
                        "bp sequence adjacent to the restriction site near the
                        predicted vector-insert junction
                        at the beginning of the read"))
        MisMatchTolerance <- ceiling(0.2 * SideSeqLength)
    ## Increase the WidthOfWindow by the max of 2*SideSeqLength or +50%
        WidthOfwindow <- max(c(ceiling(1.5 * WidthOfwindow),
                               ceiling(WidthOfwindow + 2 * SideSeqLength)))
    ## Get the corresponding larger sequence on the read that is supposed to contain the junction between vector and insert
        SeqToSearch <- subseq(ReadDNA[[rnn]],
                              PointOfSearch -
                                RestrictionSiteLength -
                                WidthOfwindow,
                              PointOfSearch + WidthOfwindow - 1)
    ## Get the sequence located upstream of the restriction site
        if (readStrand == "+") {
            MotifToSearch <- Biostrings::subseq(VectorSequence,
                                                VectorSize -
                                                  RestrictionSiteLength -
                                                  SideSeqLength + 1,
                                                VectorSize -
                                                  RestrictionSiteLength)
#~             DownSeqToSearch <- subseq(VectorSequence, RestrictionSiteLength+1, RestrictionSiteLength + SideSeqLength)
        } else {
            MotifToSearch <- Biostrings::reverseComplement(
              Biostrings::subseq(VectorSequence,
                                 RestrictionSiteLength + 1,
                                 RestrictionSiteLength + SideSeqLength))
#~             DownSeqToSearch <- reverseComplement(subseq(VectorSequence, VectorSize - RestrictionSiteLength - SideSeqLength + 1, VectorSize - RestrictionSiteLength))
        }
    ## Search for the motif around the predicted vector-insert junction
        findAdjacentSequence <- Biostrings::matchPattern(
          MotifToSearch,
          SeqToSearch,
          max.mismatch = MisMatchTolerance,
          with.indels= TRUE)
    ## If the motif is found only once, then set the vector-insert junction based on this
        if (length(findAdjacentSequence) != 1) {
            JuncAtBegin <- NA
            if (length(findAdjacentSequence) == 0) {
                message(rnn,
                        strwrap(prefix = " ", initial = "",
                                ": sequence adjacent to the restriction site
                                not found near the expected position
                                at the beginning of the read"))
            } else {
                message(rnn,
                        strwrap(prefix = " ", initial = "",
                                ": sequence adjacent to the restriction sites
                                found multiple times near the expected position
                                at the beginning of the read"))
            }
        } else {
            message(rnn,
                    strwrap(prefix = " ", initial = "",
                            ": sequence found at the beginning of the read!
                            Predicted vector-insert junction available"))
            JuncAtBegin <-  PointOfSearch -
              RestrictionSiteLength -
              WidthOfwindow +
              end(findAdjacentSequence) +
              RestrictionSiteLength  #Predicted first base of the insert
        }
    }


  #----------------------------------------
  # Search the vector-insert junction at the end of the read:
  #----------------------------------------
  ##PointOfSearch is the position of the first base preceding the restriction site
    PointOfSearch <- endAlign %>% dplyr::pull(.data$InsertBorder)
  ## Search for the restriction site in a window covering +/- 5% of DistToVectorEnd
  ## the farther we are from vector end, the wider the window
    WidthOfwindow <- max(c(15,
                           round(0.05 * (endAlign %>%
                                           dplyr::pull(.data$DistToVectorEnd)),
                                 0)))
  ## Sequence to search
  ## defined as +/- WidthOfwindow bases on each side of the expected position of the restriction site
    SeqToSearch <- Biostrings::subseq(ReadDNA[[rnn]],
                                      PointOfSearch - WidthOfwindow + 1,
                                      PointOfSearch + RestrictionSiteLength +
                                        WidthOfwindow)
  ## Search for the restriction site
    findRestrSite <- Biostrings::matchPattern(RestrictionSiteSequence,
                                              SeqToSearch,
                                              max.mismatch = 1,
                                              with.indels = TRUE)
  # Define JuncAtBegin which is the last base of the insert
    if (length(findRestrSite) != 1) {
      JuncAtEnd <- NA
      if (length(findRestrSite) == 0) {
        message(rnn,
                strwrap(prefix = " ", initial = "",
                        ": restriction site not found near the
                        expected position at the end of the read"))
      } else {
        message(rnn,
                strwrap(prefix = " ", initial = "",
                        ": multiple restriction sites found near the
                        expected position at the end of the read"))
      }
    } else {
      JuncAtEnd <-  PointOfSearch - WidthOfwindow + start(findRestrSite) - 1
    }

  # If search for the restriction site fails
  # then try to search for the vector sequence located just downstream of the restriction site
    if (is.na(JuncAtEnd)) {
        SideSeqLength = as.integer(SideSeqSearch)
        message(rnn,
                ": Searching for a ",
                SideSeqLength,
                strwrap(prefix = " ", initial = "",
                        "bp sequence adjacent to the restriction site
                        near the predicted vector-insert junction
                        at the end of the read"))
        MisMatchTolerance <- ceiling(0.2 * SideSeqLength)
    ## Increase the WidthOfWindow by the max of 2*SideSeqLength or +50%
        WidthOfwindow <- max(c(ceiling(1.5 * WidthOfwindow),
                               ceiling(WidthOfwindow + 2 * SideSeqLength)))
    ## Get the corresponding larger sequence on the read that is supposed to contain the vector-insert junction
        SeqToSearch <- Biostrings::subseq(ReadDNA[[rnn]],
                                          PointOfSearch -
                                            WidthOfwindow + 1,
                                          PointOfSearch +
                                            RestrictionSiteLength +
                                            WidthOfwindow)
    ## Get the sequence located upstream of the restriction site
        if (readStrand == "+") {
            MotifToSearch <- Biostrings::subseq(VectorSequence,
                                                RestrictionSiteLength + 1,
                                                RestrictionSiteLength +
                                                  SideSeqLength)
        } else {
            MotifToSearch <- Biostrings::reverseComplement(
              Biostrings::subseq(VectorSequence,
                                 VectorSize - RestrictionSiteLength -
                                   SideSeqLength + 1,
                                 VectorSize - RestrictionSiteLength))
        }
    ## Search for the motif around the predicted vector-insert junction
        findAdjacentSequence <- Biostrings::matchPattern(
          MotifToSearch,
          SeqToSearch,
          max.mismatch = MisMatchTolerance,
          with.indels= TRUE)
    ## If the motif is found only once, then set the vector-insert junction based on this
    ## JuncAtEnd is the position of the last base of the insert
        if (length(findAdjacentSequence) != 1) {
            JuncAtEnd <- NA
            if (length(findAdjacentSequence) == 0) {
                message(rnn,
                        strwrap(prefix = " ", initial = "",
                                ": sequence adjacent to the restriction site
                                not found near the expected position
                                at the end of the read"))
            } else {
                message(rnn,
                        strwrap(prefix = " ", initial = "",
                                ": sequence adjacent to the restriction sites
                                found multiple times near the expected position
                                at the end of the read"))
            }
        } else {
            message(rnn,
                    strwrap(prefix = " ", initial = "",
                            ": sequence found at the end of the read!
                            Predicted vector-insert junction available"))
            JuncAtEnd <-  PointOfSearch -
              WidthOfwindow +
              start(findAdjacentSequence) - 1 -
              RestrictionSiteLength
        }
    }

# For each read, report the base at which the junction occurs with the vector and where it ends
  c("ReadName" = rnn,
    "Strand" = readStrand,
    "InsertStart" = JuncAtBegin,
    "InsertEnd" = JuncAtEnd)

}
