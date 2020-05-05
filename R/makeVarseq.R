#' Introduce variations in a DNA sequence
#'
#' @param dnaseq Either a character string, a DNAString or a character vector with individual characters.
#'               The DNA sequence in which variations are introduced
#' @param lettrs character vector
#' @param subst numeric vector of length 2 with values in [0,1].
#'              Percentage of substitutions (upper and lower bounds).
#'              Default is 3-6\% substition.
#' @param del numeric vector of length 2 with values in [0,1].
#'            Percentage of deletions (upper and lower bounds).
#'            Default is 1-5\% deletion.
#' @param ins numeric vector of length 2 with values in [0,1].
#'            Percentage of insertions (upper and lower bounds).
#'            Default is 1-5\% insertions.
#' @param returnString Logical. Should the function return a single character string?
#'
#' @return either a vector of individual characters
#'         or a single character string if \code{returnString==TRUE}
#'
#' @importFrom stats runif
#' @importFrom methods is
#' @importFrom R.utils insert
#' @importFrom stringi stri_rand_strings
#'
#' @export
#'
#' @examples
#' set.seed(12345) # for reproducibility
#' ## Introduce ~10% variation in a sequence
#' makeVarseq("ATGCATGCATGCATGCATGCATGC",
#'            subst = c(0.08, 0.12),
#'            del = c(0.08, 0.12),
#'            ins = c(0.08, 0.12))

makeVarseq <- function(dnaseq,
                       lettrs = c("A", "T", "G", "C"),
                       subst = c(0.01, 0.04),
                       del = c(0.01, 0.02),
                       ins = c(0.0075, 0.0125),
                       returnString = TRUE) {

  #convert dnaseq to a list of individual chararacters
  if (is(dnaseq, "DNAString")) (dnaseq <- toString(dnaseq))
  if (is.character(dnaseq) && length(dnaseq)==1 && nchar(dnaseq)>1) {
    dnaseq <- unlist(strsplit(dnaseq, ""), use.names = FALSE)
  }
  if (!is.character(dnaseq)) {
    stop("Could not convert dnaseq to a character vector")
  }
  if (!all(nchar(dnaseq)==1)) {
    stop("dnaseq could not be converted to a vector of individual characters")
  }

  #Check the substition and indel rates
  ## substitution rate
  if (!is.numeric(subst)) {
    stop("subst is not numeric")
  }
  if (length(subst)!=2) {
    stop("subst is not of length 2")
  }
  if (!all(subst<=1 & subst>=0)) {
    stop("subst contains value not in [0,1]")
  }
  subst <- c(min(subst), max(subst))

  ## deletion rate
  if (!is.numeric(del)) {
    stop("del is not numeric")
  }
  if (length(del)!=2) {
    stop("del is not of length 2")
  }
  if (!all(del<=1 & del>=0)) {
    stop("del contains value not in [0,1]")
  }
  del <- c(min(del), max(del))

  ## insertion rate
  if (!is.numeric(ins)) {
    stop("ins is not numeric")
  }
  if (length(ins)!=2) {
    stop("ins is not of length 2")
  }
  if (!all(ins<=1 & ins>=0)) {
    stop("ins contains value not in [0,1]")
  }
  ins <- c(min(ins), max(ins))

  #insertion/deletion size probabilities
  indelSizeRange <- 1:15
  indelSizeProbs <- c(0.67, 0.22, 0.07, 0.02, 0.007, 0.002, 0.001,
                      4e-4, 2.6e-4, 7.6e-5, 4.9e-5, 3.3e-5, 2.7e-5,
                      2.17e-5, 1.6e-5)
  ## This is close to what we observed in our Nanopore data

  #number of characters
  lseq <- length(dnaseq)

  #Initialize the object
  dnares <- dnaseq

  # Substitute bases
  ## The number of susbtitutions is increased to account for the fact that some substitutions will be silent
  numsub <- round(runif(1, subst[1], subst[2]) *
                    length(lettrs) /
                    (length(lettrs)-1) * lseq, 0) # number of substitutions
  subpos <- sample.int(lseq, numsub) #position of the substitutions
  dnares[subpos] <- sample(lettrs, numsub, replace = TRUE)

  # Remove bases (deletions)
  numdel <- round(runif(n = 1,
                        min = del[1],
                        max = del[2]) * lseq, 0) # number of deletions
  delpos <- sample.int(lseq, numdel) # position of deletions
  dellen <- sample(indelSizeRange,
                   numdel,
                   prob = indelSizeProbs,
                   replace = TRUE) # length of the deletions
  delpos2 <- unique(as.integer(
    rep(delpos, times = dellen) +
      sequence(dellen) - 1)) #final positions to delete
  dnares <- dnares[-delpos2]

  #Add bases (insertions)
  numins <- round(runif(n = 1,
                        min = ins[1],
                        max = ins[2]) * lseq, 0) # number of insertions
  inspos <- sample.int(length(dnares), numins) #position of the insertions
  inslen <- sample(indelSizeRange,
                   numins,
                   prob = indelSizeProbs,
                   replace = TRUE) #length of the insertions
  inseq <- stringi::stri_rand_strings(n =numins,
                                      length = inslen,
                                      pattern = paste0("[",
                                                      paste0(lettrs,
                                                             collapse=""),
                                                      "]",
                                                      collapse = ""))
  dnares <- R.utils::insert(x = dnares,
                            at = inspos,
                            values = inseq,
                            useNames = FALSE)

  ## Return the result
  if (returnString) {
    return(paste0(dnares, collapse = ""))
  } else {
    dnares <- unlist(strsplit(dnares, ""), use.names = FALSE)
    return(dnares)
  }
}
