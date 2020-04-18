#' Check if an object was imported using \code{readBlast} or import an object using \code{readBlast}
#'
#' @param bvar either a path to a Blast result table (gerenated with outfmt=6) or
#'             such a table imported with the readBlast function
#' @param bform character string. Currently can only be "outfmt6"
#'
#' @return a Blast result table imported with the readBlast function
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' ## Create a basic data frame with the right columns
#' df <- data.frame("QueryACC" = "MyVector",
#'                  "SubjectACC" = "Read1",
#'                  "PercentID" = "97.5",
#'                  "AlnLength" = 2530,
#'                  "NumMismatch" = 20,
#'                  "NumGapOpen" = 30,
#'                  "QueryStart" = 1,
#'                  "QueryEnd" = 2500,
#'                  "SubjectStart" = 102490,
#'                  "SubjectEnd" = 100000,
#'                  "evalue" = 0.0,
#'                  "bitscore" = 3000,
#'                  "Strand")
#' ## checkBlastVar returns the object itself
#' mytab <- checkBlastVar(df)
#' ## If we mofidy a column name
#' colnames(df)[1] <- "newcolname"
#' ## then checkBlastVar throws an error
#' \dontrun{
#' mytab <- checkBlastVar(df)
#' }
#' }
checkBlastVar <- function(bvar, bform="outfmt6") {
  if (is.character(bvar)) {
    if (!file.exists(bvar)) {
      stop("File not found: ", bvar)
    } else {
      bform <- match.arg(bform, several.ok = FALSE)
      res <- readBlast(bvar, bform)
    }
  } else {
    if (is.null(bvar)) {
      res <- NULL
    } else {
      stdcoln <- c("QueryACC", "SubjectACC",
                   "PercentID", "AlnLength",
                   "NumMismatch", "NumGapOpen",
                   "QueryStart", "QueryEnd",
                   "SubjectStart", "SubjectEnd",
                   "evalue", "bitscore",
                   "Strand")
      if (!is.data.frame(bvar) || !all(colnames(bvar) == stdcoln)) {
        stop(bvar, " is not correctly formatted. Has it been imported with the readBlast function?")
      } else {
        res <- bvar
      }
    }
  }
  return(res)
}
