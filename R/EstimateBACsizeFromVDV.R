#' Estimate the size of a BAC using the length of all VDV reads
#'
#' Cluster the sizes of VDV read for different number of clusters
#' (use \code{LongestDNA} rather than \code{ReadLength}, see details).
#' For each clustering, obtain the median size for each cluster and take the largest one.
#' The estimated BAC size is the median of these larger medians
#'
#' @param vdvLength integer vector of lengths of VDV reads
#' @param nclust integer vector of numbers of clusters to use (default to 2:10)
#' @param method character string indicating which method to use ("jenks" or "pam"). Default is "jenks" (faster)
#'
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom BAMMtools getJenksBreaks
#' @importFrom cluster pam
#' @importFrom stats median
#'
#' @export
#' @return tibble
#' @details In order to use the result of this function with the \code{\link{FilterBACreads}} function
#'           you should input the "LongestDNA" values instead of the ReadLength values for VDV reads.
#' @examples
#' ## Generate some random insert sizes
#' set.seed(12345)
#' InsertLengths <- c(sample(102e3:106e3, 30),
#'                    sample(43e3:45e3, 10),
#'                    sample(23e3:25e3, 10))
#' ## Estimate the BAC size (smaller lengths are considered likely recombinants)
#' estimateBACsizeFromVDV(InsertLengths, method = "jenks")$BACsize
#'


estimateBACsizeFromVDV <- function(vdvLength,
                                    nclust = 2:10,
                                    method = c("jenks", "pam")) {

#Arguments
nclust <- as.integer(nclust)
method <- match.arg(method, several.ok = FALSE)

stopifnot(max(nclust) < length(vdvLength))

# Get the median of the clusters for nclust values
estBACsizes <- rep(NA, length(nclust))
clsters <- matrix(nrow = length(vdvLength), ncol = length(nclust))
clsmed <- list()

if (method == "jenks") {

    for (i in 1:length(nclust)) {
        bks <- BAMMtools::getJenksBreaks(vdvLength, k = nclust[i] + 1, subset = NULL)
        clsters[,i] <- findInterval(vdvLength, bks, rightmost.closed = TRUE)
        clsmed[[i]] <- sapply(split(vdvLength, clsters[,i]), median, na.rm=TRUE)
        numobs <- lengths(split(clsters[,i], clsters[,i]))
        estBACsizes[i] <- max(clsmed[[i]])
    }

} else {

    for (i in 1:length(nclust)) {
        clsters[,i] <- cluster::pam(vdvLength, k = nclust[i], cluster.only = TRUE)
        clsmed[[i]] <- sapply(split(vdvLength, clsters[,i]), median, na.rm = TRUE)
        estBACsizes[i] <- max(clsmed[[i]])
    }

}


colnames(clsters) <- paste0("k", nclust)
names(clsmed) <- paste0("k", nclust)
names(estBACsizes) <- paste0("k", nclust)

res <- list("clusters" = clsters,
            "cluster_medians" = clsmed,
            "maxMedian" = estBACsizes,
            "BACsize" = median(estBACsizes))

## TODO: Verify that the largest estimated BAC size does not correspond to a class with very few members (like 1...)

return(res)
}



