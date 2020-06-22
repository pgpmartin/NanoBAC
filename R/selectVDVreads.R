#' Select VDV reads from a set of Nanopore reads
#'
#' The function does the following:
#' \itemize{
#'   \item{Selects VDV reads}{This is done using \code{\link{FilterBACreads}}}
#'   \item{Estimate BAC size}{The size of the insert sequence is estimated using \code{\link{estimateBACsizeFromVDV}}}
#'   \item{Select VDV reads of the right size}{VDV reads within +/- \code{SizeTolerance}\% of the estimated size are selected}
#'   \item{Plot VDV read size}{The plot illustrates the size of all VDV reads and the selection process}
#'   }
#' By default, if alignment to the host genome is provided in the \code{ReadClass} object (column \code{HostAlign}),
#' then the selected VDV reads do not align to the host genome.
#' The function allows the user to exclude reads from the analysis.
#'
#' @param ReadClass Either a tibble obtained with the \code{\link{AnnotateBACreads}} function
#'                  or a path to an rds file containing such a file
#' @param ignoredReads (optional) vector of read names that should be ignored.
#' @param SizeTolerance A single number in [0,1[.
#'                      Reads with a size corresponding to the estimated size +/- \code{SizeTolerance}\% are selected
#' @param WithGeneA Logical. Should the VDV reads align with GeneA? Default is NULL, i.e. no filtering on GeneA alignment
#' @param WithGeneB Logical. Should the VDV reads align with GeneB? Default is NULL, i.e. no filtering on GeneB alignment
#' @param MaxClusters Integer. Maximum number of clusters to make with VDV reads to estimate the size of the DNA insert.
#'                             Default to 10 which in our hands works in most situations.
#' @param makePlot Logical. Should a plot be produced (and saved in the result)
#' @param plotVar Character string. Variable to plot: either \code{'ReadLength'} or \code{InsertLength}.
#'
#' @return A list with the following objects:
#' \itemize{
#'   \item{VDVreads}{A tibble with info on VDV reads: name, length, selection by the procedure, cluster}
#'   \item{InsertSizeEstimate}{A list with result from the \code{\link{estimateBACsizeFromVDV}} function}
#'   }
#' if makePlot is TRUE, then the list also contains:
#' \itemize{
#'   \item{VDVlengthPlot}{The plot (ggplot2 object)}
#'   \item{PlotType}{The variable used to make the plot, i.e. the plotVar argument that was used to create the object}
#'   \item{xcoord}{The x coordinates of the points (obtained using geom_jitter) in order to reproduce the plot easily}
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr bind_cols bind_rows select filter mutate pull
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_jitter layer_data geom_point labs
#'                     scale_y_continuous expand_scale theme_bw theme
#'                     element_blank element_text
#' @importFrom scales viridis_pal hue_pal
#'
#' @author Pascal GP Martin
#'
#' @export
#'
#' @examples
#' ## Path to file (.rds) created with the AnnotateBACreads function
#' pathRC <- system.file("extdata", "BAC02_ReadClass.rds", package = "NanoBAC")
#' ## Import the data
#' annotatedReads <- readRDS(pathRC)
#' ## Select VDV reads that contain alignment to GeneA and GeneB
#' myVDVreads <- selectVDVreads(ReadClass = annotatedReads,
#'                              WithGeneA = TRUE,
#'                              WithGeneB = TRUE)
#' ## Same but making less clusters and producing a plot on insert length rather than read length
#' myVDVreads <- selectVDVreads(ReadClass = annotatedReads,
#'                              WithGeneA = TRUE,
#'                              WithGeneB = TRUE,
#'                              MaxClusters = 8,
#'                              plotVar = "InsertLength")


selectVDVreads <- function(ReadClass = NULL,
                           SizeTolerance = 0.05,
                           WithGeneA = NULL,
                           WithGeneB = NULL,
                           ignoredReads = NULL,
                           MaxClusters = 10L,
                           makePlot = TRUE,
                           plotVar = c("ReadLength", "InsertLength")) {

#---------------------
# Arguments
#---------------------

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

## Remove ignored reads
  if (!is.null(ignoredReads)) {
    if (!all(ignoredReads %in% RC$ReadName)) {
       stop("Some ignored reads were not found in the dataset")
    }
    RC <- RC %>% dplyr::filter(!(.data$ReadName %in% ignoredReads))
    message(length(ignoredReads), " reads were removed from the dataset before analysis (ignoredReads).")
  }

## SizeTolerance
  SizeTolerance <- as.numeric(SizeTolerance)
  if (length(SizeTolerance)!=1 ||
      SizeTolerance < 0 ||
      SizeTolerance >= 1) {
    stop("SizeTolerance should be a number in [0,1[")
  }

## WithGeneA
  if (is.na(WithGeneA)) {WithGeneA <- NULL}
  if (!is.null(WithGeneA) && !is.logical(WithGeneA)) {
    stop("WithGeneA should be TRUE/FALSE (or NULL)")
  }
  if (!is.null(WithGeneA) && !("AlignedToGeneA" %in% colnames(RC))) {
    warning(strwrap(prefix = " ", initial = "",
                    "No AlignedToGeneA column found in ReadClass table.
                    No filtering done on GeneA alignement"))
  }

## WithGeneB
  if (is.na(WithGeneB)) {WithGeneB <- NULL}
  if (!is.null(WithGeneB) && !is.logical(WithGeneB)) {
    stop("WithGeneB should be TRUE/FALSE (or NULL)")
  }
  if (!is.null(WithGeneB) && !("AlignedToGeneB" %in% colnames(RC))) {
    warning(strwrap(prefix = " ", initial = "",
                    "No AlignedToGeneB column found in ReadClass table.
                    No filtering done on GeneB alignement"))
  }


## makePlot
  if (!is.finite(makePlot) || is.null(makePlot)) {
    stop("makePlot should be TRUE or FALSE")
  } else {
    makePlot <- as.logical(makePlot)
  }

## plotVar
  plotVar <- match.arg(plotVar, several.ok = FALSE)

#---------------------
# Select VDV reads
#---------------------
  if ("AlignedToGeneA" %in% colnames(RC)) {
    GnAfilter = WithGeneA
  } else {
    GnAfilter = NULL
  }

  if ("AlignedToGeneB" %in% colnames(RC)) {
    GnBfilter = WithGeneB
  } else {
    GnBfilter = NULL
  }

  if ("HostAlign" %in% colnames(RC)) {
    Hostfilter = FALSE #do not select reads that align to the host genome
  } else {
    Hostfilter = NULL
  }

  allVDV <- FilterBACreads(RC,
                           readtype = "VDV",
                           alnGeneA = GnAfilter,
                           alnGeneB = GnBfilter,
                           isHostAlign = Hostfilter)

  if (nrow(allVDV) == 0) {
    stop("No VDV reads in the dataset")
  }

#---------------------
# Estimate the length of the insert DNA
#---------------------
  ## Define the number of clusters to test:
  if (nrow(allVDV) < 3) {
    ## Only 1 cluster if 1 or 2 reads
    warning(strwrap(prefix = " ", initial = "",
                    "Less than 3 VDV reads found.
                     No clustering is done (k=1)"))
    clustNUM <- 1

  } else {
    ## k = 2 to number of VDV reads - 1 (with a max of MaxClusters)
    clustNUM = 2:min(c(MaxClusters, nrow(allVDV)-1))

  }

  BACSizeEst <- estimateBACsizeFromVDV(allVDV %>%
                                         dplyr::pull(.data$LongestDNA),
                                       nclust = clustNUM,
                                       method = "jenks")


#---------------------
# Filter VDV reads to keep those within +/- SizeTolerance% of the estimated insert size
#---------------------
  selVDV <- FilterBACreads(allVDV,
                           readtype = "VDV",
                           VDVInsertLengthTarget = BACSizeEst$BACsize,
                           VDVInsertLengthTolerance = SizeTolerance)


#---------------------
# Prepare a data frame for plotting
#---------------------

  ## Identify the clustering with the smallest k that is closest to the estimated BAC insert size
  whichk <- which.min(abs(BACSizeEst$maxMedian - BACSizeEst$BACsize))

  vdvdf <- dplyr::bind_cols(
    allVDV %>%
      dplyr::select(.data$ReadName,
                    .data$ReadLength,
                    .data$LongestDNA) %>%
      dplyr::mutate(Selected = .data$ReadName %in%
                      (selVDV %>% dplyr::pull(.data$ReadName))),
    data.frame("Cluster" = BACSizeEst$clusters[,whichk]))


#---------------------
# Plotting
#---------------------

  if (makePlot) {
  ## Create a dataset with the right variable (ReadLength or LongestDNA) as DNAlength variable
    vdvdf2 <- vdvdf
      vdvdf2 <- vdvdf2 %>%
        dplyr::rename("DNAlength" = ifelse(plotVar == "ReadLength",
                                           "ReadLength",
                                           "LongestDNA"))

  ## Create a dummy plot in order to get the x coordinates from geom_jitter:
    pdum <- ggplot2::ggplot(vdvdf2,
                            ggplot2::aes(x = "BAC",
                                         y = .data$DNAlength/1000)) +
            ggplot2::geom_jitter(height = 0)
    xPoints <- ggplot2::layer_data(pdum, i=1)$x

    vdvdf2 <- dplyr::bind_cols(vdvdf2,
                              "xcoord" = xPoints)

    NumberOfReads <- nrow(vdvdf)
    NumberOfSelectedReads <- sum(vdvdf$Selected)
    NumberOfClusters <- max(vdvdf$Cluster)

    if (max(clustNUM) > 12) {
      clustcols <- scales::viridis_pal(direction = -1)(NumberOfClusters)
    } else {
#      clustcols <- rev(scales::dichromat_pal("Categorical.12")(NumberOfClusters))
      clustcols <- scales::hue_pal()(NumberOfClusters)
    }

    newdf <- bind_rows(
      vdvdf2 %>%
        dplyr::mutate(PlotType = paste0("All reads (n=",
                                        NumberOfReads,
                                        ")"),
                      pointcolor = "#000000"),
      vdvdf2 %>%
        dplyr::mutate(PlotType = paste0("Clusters (k=",
                                        NumberOfClusters,
                                        ")"),
                      pointcolor = clustcols[vdvdf$Cluster]),
      vdvdf2 %>%
        dplyr::filter(.data$Selected == TRUE) %>%
        dplyr::mutate(PlotType = paste0("Selected Reads (n=",
                                        NumberOfSelectedReads,
                                        ")"),
                      pointcolor = "#0000FF")
    )

pp <- newdf %>%
  ggplot2::ggplot(ggplot2::aes(x = .data$xcoord,
                               y = .data$DNAlength/1000,
                               color = .data$pointcolor)) +
    ggplot2::geom_point() +
    ggplot2::scale_colour_identity() +
    ggplot2::labs(x=NULL,
                  y = paste0(ifelse(plotVar == "ReadLength",
                                    "Read", "Insert"),
                             " Length (kb)")) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.line.x = ggplot2::element_blank(),
                   legend.position = "none",
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank()
                   ) +
    ggplot2::facet_grid(cols = ggplot2::vars(.data$PlotType))

print(pp)
  }

  if (makePlot) {
    fres <- list("VDVreads" = vdvdf,
                "InsertSizeEstimate" = BACSizeEst,
                "VDVlengthPlot" = pp,
                "PlotType" = plotVar,
                "xcoord" = xPoints)
  } else {
    fres <- list("VDVreads" = vdvdf,
                "InsertSizeEstimate" = BACSizeEst)
  }

return(fres)

}
