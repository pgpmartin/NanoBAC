#' Obtain the consensus sequence from an MSF file
#'
#' @param msfALN Character string. Path to an MSF file obtained with a multiple aligment program
#' @param gapthreshold Numeric in ]0,1[. (Default is 0.5)
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
#'
consensusFromMSF <- function(msfALN = NULL,
                             gapthreshold = 0.5,
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

  gapthreshold <- as.numeric(gapthreshold)
  if (gapthreshold<0) {
    stop("gapthreshold cannot be negative")
  }
  if (gapthreshold>1) {
    stop("gapthreshold cannot be > 1")
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
  # get the consensus matrix
  #---------------------

  mmat <- Biostrings::consensusMatrix(msf)

  #---------------------
  # Obtain the consensus
  #---------------------

  cons <- consmat2seq(mmat,
                      ambiguityLetter = "N",
                      threshold = gapthreshold)

  #---------------------
  # Convert to a DNAStringSet
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
