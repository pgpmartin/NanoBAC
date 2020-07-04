#' Test if the argument of a function is of a given class or import the file corresponding to the argument and test its class
#'
#' @description
#'   If the \code{avar} object is of the class \code{expectedClass}, then the function returns the \code{avar} object as is
#'   If \code{avar} is NULL, then the function returns NULL
#'   If \code{avar} is a character string corresponding to an existing file, then the file is imported using \code{importFUN}
#'   Else, the function stops and returns an error
#'
#' @param avar Either an object or a path to a file to import with \code{importFUN}
#' @param importFUN Function to import \code{avar} if it is a path. The function must return an object of class \code{expectedClass}
#' @param expectedClass Character vector. expected class of the object avar
#' @param ... further arguments passed to \code{importFUN}
#'
#' @return either NULL if \code{avar} is NULL or an object of class \code{expectedClass}
#'
#' @importFrom methods is
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # With an object in R
#' TestArg(1:10, expectedClass = "integer") #returns 1:10
#' TestArg(1:10, expectedClass = "character") #returns an error
#'
#' #With a path to a file
#' Path2Blast <- system.file("extdata", "BAC02_Blast18S.res", package = "NanoBAC")
#' TestArg(Path2Blast, readBlast, "tbl_df") #returns the imported file
#' TestArg(Path2Blast, readBlast, "integer") #returns an error
#' }

TestArg <- function(avar, importFUN, expectedClass, ...) {
  if (is.character(avar) && length(avar)==1) {
    if (!file.exists(avar)) {
      stop("File not found: ", avar)
    } else {
      res <- importFUN(avar, ...)
      if (!is(res, expectedClass)) {
        stop("Imported object is not of the expected class: ", expectedClass)
      }
    }
  } else {
    if (is.null(avar)) {
      res <- NULL
    } else {
      if (!is(avar, expectedClass)) {
        stop("Variable is not of the expected class: ", expectedClass)
      } else {
        res <- avar
      }
    }
  }
  return(res)
}


