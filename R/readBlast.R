#' Read a blast result file
#'   (currently only for files created using the blast parameter outfmt=6)
#'
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom utils read.table
#' @importFrom rlang .data
#' @param blastpath character string. Path to a Blast result file
#' @param blastformat format of the Blast result file (currently only the default outfmt=6 is supported)
#' @export
#' @return a tibble
#' @examples
#' ## Example data set:
#'     Path2Blast <- system.file("extdata", "BAC02_Blast18S.res", package = "NanoBAC")
#' ## Import the data:
#'    myblastResult <- readBlast(Path2Blast)
#'

readBlast <- function(blastpath, blastformat="outfmt6") {

    stopifnot(file.exists(blastpath))
    blastformat <- match.arg(blastformat, several.ok = FALSE)
    stopifnot(blastformat == "outfmt6")

# Import table
    res <- utils::read.table(blastpath,
                        sep="\t", dec=".",
                        header = FALSE,
                        stringsAsFactors = FALSE)

# Convert to tibble
    res <- tibble::as_tibble(res)

# Check columns and name them
    stopifnot(ncol(res) == 12)
    colnames(res) <- c("QueryACC", "SubjectACC",
                       "PercentID", "AlnLength",
                       "NumMismatch", "NumGapOpen",
                       "QueryStart", "QueryEnd",
                       "SubjectStart", "SubjectEnd",
                       "evalue", "bitscore")

# Add strand info
    res <- res %>%
        dplyr::mutate(Strand = ifelse((.data$SubjectStart - .data$SubjectEnd) > 0,
                                      "-", "+"))

# Return result
    return(res)

}

