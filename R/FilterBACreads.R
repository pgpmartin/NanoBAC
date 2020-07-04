#' Select specific reads from a table of BAC read annotations obtained by the function \code{AnnotateBACreads}
#'
#' @param ReadClass Either a tibble obtained with the \code{\link{AnnotateBACreads}} function or
#'                  a path to an RDS file containing such a tibble
#' @param readtype Character vector. Type(s) of reads to select.
#'                                   Values should be in c("VDV", "DVD", "VD", "D", "V", "Chimeric", NA).
#'                                   (Defaults to "VDV")
#' @param MinReadLength Integer. Minimum length of the selected reads. Reads below this length are dropped.
#'                               (defaults to NULL, i.e. no filtering on read length)
#' @param alnGeneA Logical. Should the selected read align to GeneA ? (defaults to NULL, i.e no filtering on GeneA alignment)
#' @param alnGeneB Logical. Should the selected read align to GeneB ? (defaults to NULL, i.e no filtering on GeneB alignment)
#' @param isHostAlign Logical. Should the selected reads align to the host genome (TRUE) or not (FALSE)?
#'                             (defaults to NULL, i.e no filtering on host genome alignment alignment)
#' @param VDVInsertLengthTarget Integer. What is the estimated length of the DNA insert (without the vector) in the VDV reads to be selected?
#'                                     (Defaults to NULL, no filtering on VDV read length)
#' @param VDVInsertLengthTolerance Numeric. Number in [0,1[. What tolerance (in \%) is accepted around VDVInsertLengthTarget?
#'
#' @importFrom dplyr filter pull
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export
#'
#' @return tibble with selected reads.
#'
#' @details The function uses the column "LongestDNA" to filter for the size of VDV reads
#' @examples
#' ## Path to file (.rds) created with the AnnotateBACreads function
#' pathRC <- system.file("extdata", "BAC02_ReadClass.rds", package = "NanoBAC")
#' ## Extract the annotation of all VDV reads longer than 111 kb
#' FilterBACreads(pathRC, "VDV", 111e3)
#' ## Extract the annotation of the VDV reads which have an insert that is 104 kb long (+/- 0.1%)
#' FilterBACreads(pathRC, "VDV",
#'                VDVInsertLengthTarget = 104e3,
#'                VDVInsertLengthTolerance = 0.001)
#' ## Extract the annotation of the D reads that map to the E. coli genome
#' FilterBACreads(pathRC, "D", isHostAlign = TRUE)
#' ## Extract the annotation of the VD and the VDV reads that align to both GeneA and GeneB
#' FilterBACreads(pathRC, c("VD", "VDV"), alnGeneA = TRUE, alnGeneB = TRUE)
#' ## There are no VDV reads that align to the E. coli host genome
#' \dontrun{
#' FilterBACreads(pathRC, "VDV", isHostAlign = TRUE)
#' }


FilterBACreads <- function(ReadClass = NULL,
                           readtype = "VDV",
                           MinReadLength = NULL,
                           alnGeneA = NULL,
                           alnGeneB = NULL,
                           isHostAlign = NULL,
                           VDVInsertLengthTarget = NULL,
                           VDVInsertLengthTolerance = 0.05) {


#-------------------
# Test on arguments
#-------------------
  ## ReadClass argument
  if (is.null(ReadClass)) {
    stop(strwrap(prefix = " ", initial = "",
                 "Please provide a ReadClass object.
                 either a tibble obtained using the AnnotateBACreadsfunction
                 or a path to an RDS file containing such a tibble"))
  }

  RC <- TestArg(ReadClass, readRDS, "tbl_df")

  expectedColNames <- c("ReadName", "ReadLength", "Strand",
                        "DNAlength", "LongestDNA", "ShortestDNA", "ReadType")

  if (ncol(RC)<7 || !identical(colnames(RC)[1:7], expectedColNames)) {
    stop("ReadClass does not have the expected first 7 columns")
  }

  ## readtype
  readtype <- match.arg(readtype,
                        choices=c("VDV", "DVD", "VD",
                                  "D", "V", "Chimeric", NA),
                        several.ok = TRUE)

  ## alnGeneA
  if (!is.null(alnGeneA)) {
    if (!("AlignedToGeneA" %in% colnames(RC))) {
      warning(strwrap(prefix = " ", initial = "",
                   "ReadClass does not contain a columns corresponding to
                   the alignment to GeneA. Use another ReadClass object
                   or set alnGeneA to NULL.
                   No filtering is done on alignment to GeneA."))
      alnGeneA <- NULL
    } else {
      alnGeneA <- as.logical(alnGeneA)
    }
  }

  ## alnGeneB
  if (!is.null(alnGeneB)) {
    if (!("AlignedToGeneB" %in% colnames(RC))) {
      warning(strwrap(prefix = " ", initial = "",
                   "ReadClass does not contain a columns corresponding to
                   the alignment to GeneB. Use another ReadClass object
                   or set alnGeneB to NULL.
                   No filtering is done on alignment to GeneB."))
      alnGeneB <- NULL
    } else {
      alnGeneB <- as.logical(alnGeneB)
    }
  }

  ## MinReadLength
  if (!is.null(MinReadLength)) {
    MinReadLength <- as.numeric(MinReadLength)
    if (!is.numeric(MinReadLength)) {
      stop("Provide a valid read length for MinReadLength")
    }
  }

  ## VDVInsertLengthTarget and VDVInsertLengthTolerance
  if (!is.null(VDVInsertLengthTolerance)) {
    VDVInsertLengthTolerance <- as.numeric(VDVInsertLengthTolerance)
  }

  if (!is.null(VDVInsertLengthTarget)) {
    if (!is.numeric(VDVInsertLengthTolerance) ||
        VDVInsertLengthTolerance < 0 ||
        VDVInsertLengthTolerance >= 1) {
      stop(strwrap(prefix = " ", initial = "",
                   "VDVInsertLengthTolerance should be a number
                   in [0,1[ in order to filter on VDV read length"))
    } else {
      VDVInsertLengthTarget <- as.integer(VDVInsertLengthTarget)
      if (!is.numeric(VDVInsertLengthTarget)) {
        stop("Provide a valid read length for VDVInsertLengthTarget")
      }
    }
  }

  ## Host alignment
  if (!is.null(isHostAlign)) {
    isHostAlign <- as.logical(isHostAlign)
    if (is.na(isHostAlign)) {
        warning(strwrap(prefix = " ", initial = "",
                        "The value entered for isHostAlign is not TRUE or FALSE.
                        No filtering is done on HostAlign"))
        isHostAlign <- NULL
    }
  }

  if (!is.null(isHostAlign)) {
    if (!("HostAlign" %in% colnames(RC))) {
      warning(strwrap(prefix = " ", initial = "",
                      "ReadClass object does not contain a HostAlign column.
                      Either provide another ReadClass object or
                      set isHostAlign to NULL.
                      No filtering is done on HostAlign."))
      isHostAlign <- NULL
    }
  }


  ###-----------------
  ### Filter for read type
  ###-----------------

  RC <- RC %>%
      dplyr::filter(.data$ReadType %in% readtype)

  ###-----------------
  ### Filter for min read length
  ###-----------------
  if (!is.null(MinReadLength)) {
    RC <- RC %>%
      dplyr::filter(.data$ReadLength >= MinReadLength)
  }

  ###-----------------
  ### Filter for GeneA alignment
  ###-----------------
  if (!is.null(alnGeneA)) {
    RC <- RC %>%
      dplyr::filter(.data$AlignedToGeneA == alnGeneA)
  }

  ###-----------------
  ### Filter for GeneB alignment
  ###-----------------
  if (!is.null(alnGeneB)) {
    RC <- RC %>%
      dplyr::filter(.data$AlignedToGeneB == alnGeneB)
  }

  ###-----------------
  ### Filter for VDV read length
  ###-----------------

  if (("VDV"  %in% readtype) & !is.null(VDVInsertLengthTarget)) {
    InsertSize <- VDVInsertLengthTarget
    tolerance <- VDVInsertLengthTolerance
    vdvReadsOK <- RC %>%
            dplyr::filter(.data$ReadType == "VDV") %>%
            dplyr::filter(.data$LongestDNA  >= (1 - tolerance) * InsertSize,
                          .data$LongestDNA <= (1 + tolerance) * InsertSize) %>%
            dplyr::pull(.data$ReadName)
    RC <- RC %>%
      dplyr::filter(.data$ReadType != "VDV" |
                      .data$ReadName %in% vdvReadsOK)
  }



  ###-----------------
  ### Filter for reads with or without host alignment
  ###-----------------

  if (!is.null(isHostAlign)) {
    RC <- RC %>%
      dplyr::filter(.data$HostAlign == isHostAlign)
  }

  ###-----------------
  ### Return result
  ###-----------------
  if (nrow(RC) == 0) {
      warning("Zero reads selected")
      return(NULL)
  } else {
      message(nrow(RC), " reads selected")
      return(RC)
  }

}
