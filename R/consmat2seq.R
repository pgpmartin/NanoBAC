#' Get the consensus sequence from a consensus matrix
#'
#' This function is a simplified version of the \code{\link[Biostrings:letterFrequency]{consensusString}}
#' function that:
#' \enumerate{
#'   \item only allows bases A, T, G, C and possibly a single ambiguity letter (default to N)
#'   \item does not use the IUPAC ambiguity codes in the input and output
#'   \item prioritizes the insertion of a gap when the \% of sequence with a gap is >\code{threshold}
#'   \item uses only the non gap sequences to evaluate the \% of each letter at a given position
#' }
#' Note that the function will not work if one of the base A, T, G or C or
#' if the gap "-" is completely absent from the alignment.
#' Surprisingly, the \code{\link[Biostrings:letterFrequency]{consensusString}} function does not give identical
#' results when used on a \code{DNAalignment} object or a frequency matrix given
#' by \code{consensusMatrix} with \code{prob=TRUE}. Using the default values,
#' this function will give the same result as
#' \code{consensusString(DNAalignment-object, ambiguityMap="N", threshold=0.5)}.
#'
#' @param x a consensus matrix (i.e. a matrix of integers) obtained generally
#'          using the \code{\link[Biostrings:letterFrequency]{consensusMatrix}} function.
#' @param ambiguityLetter Letter used when there is an ambiguity. (Default is "N")
#' @param threshold \% above which a base (or a gap) is selected as the consensus
#'
#' @return A character string
#'
#' @importFrom S4Vectors isSingleNumber
#'
#' @keywords internal
#'
#' @examples
#' # Create a function to compare the consmat2seq and the consensusString functions:
#' ccons <- function(intmat) {
#'     print(Biostrings::consensusString(intmat / colSums(intmat),
#'                                       "N", threshold = 0.5))
#'     print(NanoBAC:::consmat2seq(intmat, "N", threshold = 0.5))
#'     }
#' # Create a simple matrix:
#' consmat <- matrix(c(10L, 0L, 0L, 0L, 0L, 0L,
#'                     0L, 10L, 0L, 0L, 0L, 0L,
#'                     0L, 0L, 10L, 0L, 0L, 0L,
#'                     0L, 0L, 0L, 10L, 0L, 0L),
#'                     nrow = 6,
#'                     dimnames = list(c("A", "T", "G", "C", "N", "-")))
#' ccons(consmat) # same result
#' consmat[,3] <- c(4L, 0L, 0L, 0L, 0L, 6L)
#' ccons(consmat) # same
#' consmat[,3] <- c(4L, 1L, 0L, 0L, 0L, 5L)
#' ccons(consmat) # same
#' consmat[,3] <- c(5L, 0L, 0L, 0L, 0L, 5L)
#' ccons(consmat) # different. (favor the gap)
#' consmat[,3] <- c(4L, 2L, 0L, 0L, 0L, 4L)
#' ccons(consmat) #different (only consider letters)
#' consmat[,3] <- c(4L, 2L, 1L, 1L, 0L, 2L)
#' ccons(consmat) #different (only consider letters)
#' consmat[,3] <- c(4L, 2L, 1L, 1L, 1L, 1L)
#' ccons(consmat) #same
#'
consmat2seq <- function(x,
                        ambiguityLetter = "N",
                        threshold = 0.5) {

  ## x must be an integer matrix (e.g. given by consensusMatrix)
  if (!is.integer(x)) {
    stop("x should be a matrix of integers")
  }

  ## All sums by column must be the same
  if (any(colSums(x, na.rm =TRUE) != sum(x[,1], na.rm = TRUE))) {
    stop("All columns of x do not have the same sum")
  }

  ## threshold must be a single number in ]0,1]
  if (!S4Vectors::isSingleNumber(threshold) ||
      threshold <= 0 ||
      threshold > 1) {
    stop("'threshold' must be a numeric in ]0, 1]")
  }

  # Remove rows with only 0 values, except if they are in A, T, G, C, N or -
  keepRow <- rowSums(x, na.rm = TRUE) > 0 |
    rownames(x) %in% c("A", "T", "G", "C", "-")
  x <- x[keepRow, , drop = FALSE]

  # Stop if we don't have 5 or 6 rows with the expected letters
  if (ambiguityLetter %in% rownames(x)) {
    isN <- TRUE
    stopifnot(nrow(x) == 6)
    stopifnot(all(c("A", "T", "G", "C", ambiguityLetter, "-") %in% rownames(x)))
  } else {
    isN <- FALSE
    stopifnot(nrow(x) == 5)
    stopifnot(all(c("A", "T", "G", "C", "-") %in% rownames(x)))
  }

  all_letters <- rownames(x)
  true_letters <- all_letters[all_letters != "-"]

  # Convert x to probabilities
  ## A, T, G, C and N have proportions of reads with this letter, without counting the gaps
  ## and "-" has total proportion of gap
  numAlign <- sum(x[,1], na.rm = TRUE)
  numLetters <- colSums(x[true_letters,], na.rm = TRUE)
  numLetters[numLetters == 0] <- 1
  xtot <- rbind(x[true_letters,] / rep(numLetters, each = 4 + isN),
                x["-",] / numAlign)
  rownames(xtot) <- c(true_letters, "-")

  #Function to assign the letter based on the probabilities
  consensusLetter <- function(col, hasN) {
    if (hasN) (colnum=6) else (colnum=5)
    if (col[colnum] >= threshold) {
      "-"
    } else {
      i <- which(col[-colnum] >= threshold)
      if (length(i) == 1)
        true_letters[i]
      else ambiguityLetter
    }
  }

  paste(apply(xtot, 2, consensusLetter, hasN=isN), collapse = "")
}
