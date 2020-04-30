#' Obtain the consensus sequence from an MSF file
#'
#' @param msfALN Character string. Path to an MSF file obtained with a multiple aligment program
#' @param threshold Numeric in ]0,1[. (Default is 0.5)
#' @param removegaps Logical. Should the gaps be removed in the output sequence? (Default is TRUE)
#' @param consname Character string. Name of the consensus sequence. (Default is "cons")
#'
#' @return A DNAStringSet of length 1 with the consensus sequence
#'
#' @importFrom seqinr read.alignment
#' @importFrom Biostrings DNAMultipleAlignment consensusMatrix DNAStringSet
#' @importClassesFrom Biostrings DNAStringSet
#'
#' @export
#'
#' @examples
#' ## example msf file (obtained with kalign2 using the varseq100.fa file)
#' exMSFfile <- system.file("extdata", "varseq100.msf", package = "NanoBAC")
#' ## get the consensus sequence
#' consensusFromMSF(exMSFfile)
#'
consensusFromMSF <- function(msfALN = NULL,
                             threshold = 0.5,
                             removegaps = TRUE,
                             consname = "cons") {

  #---------------------
  #Test arguments and import data if necessary
  #---------------------

  if (is.null(msfALN)) {
    stop("Arguments msfALN is missing")
  }
  if (!file.exists(msfALN)) {
    stop("File ", msfALN, " not found")
  }

  threshold <- as.numeric(threshold)
  if (threshold<0) {
    stop("threshold cannot be negative")
  }
  if (threshold>1) {
    stop("threshold cannot be > 1")
  }

  #---------------------
  # import the multiple alignment MSF file
  #---------------------
  msf <- seqinr::read.alignment(msfALN,
                                format = "msf",
                                forceToLower = FALSE)
  seqs <- msf$seq
  names(seqs) <- msf$nam
  msf <- Biostrings::DNAMultipleAlignment(seqs, use.names=TRUE)
  # Also possible to use the msa::msaConvertfunction but would add a dependency

  #---------------------
  # get the consensus String
  #---------------------

  cons <- Biostrings::consensusString(msf,
                                      ambiguityMap="N",
                                      threshold = threshold)

  #---------------------
  # Convert to a DNAStringSet and remove gaps if required
  #---------------------
  if (removegaps) {
    cons <- Biostrings::DNAStringSet(gsub("-", "", cons))
  } else {
    cons <- Biostrings::DNAStringSet(cons)
  }
  names(cons) <- consname

  # Return the consensus sequence
  return(cons)
}
