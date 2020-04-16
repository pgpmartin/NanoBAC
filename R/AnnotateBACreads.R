#' Check a Blast result variable
#'
#' @param bvar either a path to a Blast result table (gerenated with outfmt=6) or
#'             such a table imported with the readBlast function
#' @param bform character string. Currently can only be "outfmt6"
#'
#' @return a Blast result table imported with the readBlast function
#'
#' @keyword internal
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
      if (!is.data.frame(bvar) || colnames(bvar) != stdcoln) {
        stop(bvar, " is not correctly formatted. Has it been imported with the readBlast function?")
      } else {
        res <- bvar
      }
    }
  }
  return(res)
}


#' Test the argument of a function
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
#'
#' @return either NULL if \code{avar} is NULL or an object of class \code{expectedClass}
#'
#' @keyword internal
#'
#' @examples
#' \donttest{
#' # With an object in R
#' TestArg(1:10, expectedClass = "integer") #returns 1:10
#' \dontrun{
#' TestArg(1:10, expectedClass = "character") #returns an error
#' }
#'
#' #With a path to a file
#' Path2Blast <- system.file("extdata", "blastExample.tab", package = "NanoBAC")
#' TestArg(Path2Blast, readBlast, "tbl_df") #returns the imported file
#' \dontrun{
#' TestArg(Path2Blast, readBlast, "integer") #returns an error
#' }
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




#' Annotate BAC reads
#'
#' @description
#'     The function annotates BAC reads based on:
#'     \itemize{
#'       \item{"readfq"}{Read quality}
#'       \item{"blastvec"}{Blast alignment of the vector on the reads}
#'       \item{"blastGeneA"}{Optional Blast alignment of gene A on the reads}
#'       \item{"blastGeneB"}{Optional Blast alignment of gene B on the reads}
#'       \item{"pafHost"}{Optional minimap2 alignment of the reads on the host genome}
#'     }
#'
#' @param blastvec Alignment of the vector on the reads. A blast result table imported with readBlast or a path to such a table.
#' @param blastGeneA Alignment of gene A (e.g. 18S) on the reads. A blast result table imported with readBlast or a path to such a table.
#' @param blastGeneB Alignment of gene B (e.g. 28S) on the reads. A blast result table imported with readBlast or a path to such a table.
#' @param pafHost Alignment of the reads on host genome (e.g. E Coli). A table obtained with the read_paf function or a path to a paf file (e.g. obtained with minimap2)
#' @param minHost_mapQ Integer. Minimum mapQ value to consider the alignment to the host genome as significant
#' @param readLength data.frame with 2 colums: ReadName and ReadLength. Provide directly the data frame or the path to the tab-delimited file containing the data (without header).
#' @param vectorSequence Either a DNAStringSet with 1 element or the path to a fasta file
#' @param minaln Integer. Minimum alignment length on GeneA and GeneB
#' @param MinDVDsides. Integer. Minimum length on each side of the vector sequence to be considered a DVD read
#'
#' @importFrom Biostrings readDNAStringSet
#' @importFrom dplyr pull select distinct add_count filter rename left_join group_by top_n
#' @importFrom rlang .data
#' @importFrom tibble as_tibble tibble
#' @importFrom magrittr %>%
#' @importFrom GenomicRanges seqnames setdiff width coverage reduce
#' @importFrom GenomeInfoDb seqinfo
#' @importFrom methods as
#' @importFrom S4Vectors split
#'
#' @return A tibble
#'
#' @export
#'
#' @examples
#'
AnnotateBACreads <- function(
                             blastvec,
                             blastGeneA = NULL,
                             blastGeneB = NULL,
                             pafHost = NULL,
                             minHost_mapQ = 10L,
                             readLength = NULL,
                             vectorSequence = NULL,
                             minaln = 1L,
                             MinDVDsides = 10e3L
                             ) {

  #---------------------
  #Test arguments and import data if necessary
  #---------------------

  ## Blast results
    vecaln <- checkBlastVar(blastvec)
    gnAaln <- checkBlastVar(blastGeneA)
    gnBaln <- checkBlastVar(blastGeneB)

    ## pafHost (minimap2 result as paf file)
    hostaln <- TestArg(pafHost,
                       importFUN = read_paf,
                       expectedClass = "tbl_df")
    ExpectedPafColnames <- c("query_name", "query_length",
                             "query_start", "query_end",
                             "strand", "target_name",
                             "target_length", "target_start",
                             "target_end", "map_match",
                             "map_length", "map_quality")
    if (!is.null(hostaln)) {
      if (ncol(hostaln) < 12 ||
          colnames(hostaln)[1:12] != ExpectedPafColnames) {
        stop("pafHost does not have the expected first 12 columns")
      }
    }

    ## readLength
    RL <- TestArg(readLength,
                  importFUN = read.table,
                  expectedClass = "data.frame",
                  sep="\t",
                  header = FALSE,
                  stringsAsFactors = FALSE,
                  col.names = c("ReadName", "ReadLength"))
    if (ncol(RL)<2 || colnames(RL)[1:2]!=c("ReadName","ReadLength")) {
      stop("readLength does not have the ReadName and ReadLength columns")
    }

    ### Verify that all reads in the alignment files are also in readLength
    isOK_blastvec <- all(unique(vecaln$SubjectACC) %in%
                           RL$ReadName)
    if (!isOK_blastvec) {
      stop("Some reads in blastvec are not in readLength")
    }

    if (is.null(gnAaln)) {
      isOK_blastGeneA <- TRUE
    } else {
      isOK_blastGeneA <- all(unique(gnAaln$SubjectACC) %in%
                               RL$ReadName)
    }
    if (!isOK_blastGeneA) {
      stop("Some reads in blastGeneA are not in readLength")
    }

    if (is.null(gnBaln)) {
      isOK_blastGeneB <- TRUE
    } else {
      isOK_blastGeneB <- all(unique(gnBaln$SubjectACC) %in%
                               RL$ReadName)
    }
    if (!isOK_blastGeneB) {
      stop("Some reads in blastGeneB are not in readLength")
    }

    if (is.null(hostaln)) {
      isOK_pafHost <- TRUE
    } else {
      isOK_pafHost <- all((hostaln %>%
                            dplyr::pull(.data$query_name) %>%
                            unique()) %in%
                          RL$ReadName)
    }
    if (!isOK_pafHost) {
      stop("Some reads in pafHost are not in readLength")
    }


    ## vectorSequence
    backboneSeq <- TestArg(vectorSequence,
                           importFUN =
                             function(x){Biostrings::readDNAStringSet(x)[[1]]},
                           expectedClass = "DNAString")
    backboneLength <- length(backboneSeq)
    message("Backbone/Vector size is ", backboneLength, " bp")


    #---------------------
    #Filter blastvec
    # Remove all alignments that overlap with another alignment on more than 50% of their length
    #---------------------

    vecalngr <- blaST2GR(vecaln, readlength = RL)
    alnToKeep <- SelectSingularBlastALN(vecaln, RL, threshold = 0.5)

    vecaln <- vecaln[alnToKeep,]
    vecalngr <- vecalngr[alnToKeep]


    #---------------------
    #Identify chimeric reads
    # These reads have a vector alignment on both strands
    #---------------------

    ChimericReads <- vecaln %>%
      dplyr::select(.data$SubjectACC, .data$Strand) %>%
      dplyr::distinct() %>%
      dplyr::add_count(.data$SubjectACC) %>%
      dplyr::filter(n==2) %>%
      dplyr::pull(.data$SubjectACC) %>%
      as.character() %>%
      unique()

    if (length(ChimericReads)!=0) {
      message("Found reads with vector alignments on both strands: ",
              paste(ChimericReads, collapse = ", "))
    }

    # Remove Chimeric reads
    vecaln <- vecaln %>%
      dplyr::filter(!(.data$SubjectACC %in%
                        ChimericReads))
    vecalngr <- vecalngr[!(GenomicRanges::seqnames(vecalngr) %in%
                             ChimericReads)]

    #---------------------
    # Get the strand of (non chimeric) reads that have an alignment with the vector
    # For the reads that don't have an alignment with the vector we don't define the strand
    # It could be defined afterwards based on alignment of e.g. GeneA or GeneB
    #---------------------

    ReadStrand <- tibble::as_tibble(vecalngr) %>%
      dplyr::select(.data$seqnames, .data$strand) %>%
      dplyr::distinct() %>%
      dplyr::rename("ReadName" = "seqnames",
                    "Strand" = "strand")

    #---------------------
    # Get the length of DNA fragments not aligned with the vector
    # And the length of the longest such non-vector DNA fragment
    #---------------------
    ## For each read, get the regions not aligned with the vector
    NonVectorDNA <- GenomicRanges::setdiff(
      methods::as(GenomeInfoDb::seqinfo(vecalngr),
                  "GRanges"),
      vecalngr,
      ignore.strand = TRUE)

    ## Get the length of these regions as an IntegerList
    NonVectorDNA_length <- S4Vectors::split(GenomicRanges::width(NonVectorDNA),
                                            GenomicRanges::seqnames(NonVectorDNA))

    ## Get the longest DNA region
    NonVectorDNA_Longest <- max(NonVectorDNA_length)

    ## Get the shortest DNA region
    NonVectorDNA_Shortest <- min(NonVectorDNA_length)


    #---------------------
    # D reads (DNA only)
    # Redas with NO alignment to the vector
    #---------------------

    #These reads have no alignments to the vector and are not ChimericReads
    Dreads <- setdiff(RL$ReadName,
                      vecaln %>%
                        dplyr::pull(.data$SubjectACC) %>%
                        as.character %>%
                        unique)

    #---------------------
    # V reads (Vector only)
    # These reads align with the vector on >30% of their length
    #---------------------

    ## Coverage of the read by vector alignment
    ReadVecCov <- GenomicRanges::coverage(GenomicRanges::reduce(vecalngr))
    ## Select reads with >30% of vector alignment
    Vreads <- names(ReadVecCov)[mean(ReadVecCov) > 0.3]
    ## Remove potential chimeric reads
    Vreads <- setdiff(Vreads, ChimericReads)

    #---------------------
    # VDV reads (Vector-DNA-Vector reads)
    #   vector        DNA          vector
    #  --------|=================|-------
    # VDV reads have a vector alignment within their first 1kb and within their last 1kb
    #     they have NO vector alignment within [+10kb and -10kb] (10kb replaced by backboneLength*1.15)
    #---------------------

    VDVreads <- suppressWarnings(
      getVDVnames(alignGR = vecalngr,
                  vectorLength = backboneLength,
                  minReadLength = 2e3,
                  SideWidth = 1e3))
    VDVreads <- setdiff(VDVreads, ChimericReads)
    VDVreads <- setdiff(VDVreads, Vreads)

    #---------------------
    # DVD reads (DNA-vector-DNA)
    #       DNA    -   vector   -  DNA
    #    ==========|------------|===========
    # DVD reads have a single vector alignment covering at least 95% (default) of the vector length
    #     and have 2 pieces of insert/DNA on each side that are at least MinDVDsides long
    #---------------------
    DVDreads <- suppressWarnings(
      getDVDnames(alignGR = vecalngr,
                  vectorLength = backboneLength,
                  PercentVecLength = 0.95,
                  MinDNASides = MinDVDsides))
    # check that there is no overlap between VDV and DVD reads
    stopifnot(length(intersect(DVDreads, VDVreads)) == 0)
    # Remove potential overlap with ChimericReads:
    DVDreads <- setdiff(DVDreads, ChimericReads)
    # Remove potential overlap with Vreads:
    DVDreads <- setdiff(DVDreads, Vreads)

    #---------------------
    # VD reads (vector-DNA)
    #      vector   DNA
    #     --------|==========
    # or     DNA       vector
    #     ===========|---------
    # VD reads have a single vector alignment that covers tolerance% of their first (or last) EndWindow bp
    #---------------------

    VDreads <- suppressWarnings(
      getVDnames(alignGR = vecalngr,
                 tolerance = 0.8,
                 EndWindow = 500)
    )

    VDreads <- setdiff(VDreads, ChimericReads) #Should not be necessary because Chimeric Reads have 2 alignments
    VDreads <- setdiff(VDreads, Vreads)
    VDreads <- setdiff(VDreads, DVDreads) #there shouldn't  be any intersection if MinDNASides and EndWindow are compatible
    VDreads <- setdiff(VDreads, VDVreads) #there shouldn't  be any intersection since VDV have 2 alignments and VD only one


    #-----------------
    # Format the output
    #-----------------

    rn <- RL %>% dplyr::pull(.data$ReadName)


    res <- dplyr::left_join(as_tibble(RL),
                            ReadStrand,
                            by = "ReadName") %>%
      dplyr::left_join( . ,
                 tibble::tibble("ReadName" = names(NonVectorDNA_length),
                                "DNAlength" = as.list(NonVectorDNA_length),
                                "LongestDNA" = NonVectorDNA_Longest,
                                "ShortestDNA" = NonVectorDNA_Shortest),
                 by = "ReadName") %>%
      dplyr::left_join( . ,
                 tibble::tibble(ReadName = rn,
                        ReadType = case_when(
                          rn %in% ChimericReads ~ "Chimeric",
                          rn %in% Dreads ~ "D",
                          rn %in% Vreads ~ "V",
                          rn %in% DVDreads ~ "DVD",
                          rn %in% VDVreads ~ "VDV",
                          rn %in% VDreads ~ "VD"
                        )),
                 by = "ReadName")

    #---------------------
    # Add (optional) info about alignments to GeneA and GeneB
    #---------------------
    # If the data is available:
    if (!is.null(gnAaln)) {
      # Get the names of reads with an alignment (longer than minaln) to GeneA:
      hasGeneA <- gnAaln %>%
        dplyr::filter(.data$AlnLength >= minaln) %>%
        dplyr::pull(.data$SubjectACC) %>% unique()

      # Add the corresponding column to the results:
      res <- res %>%
          dplyr::mutate(AlignedToGeneA = (.data$ReadName %in% hasGeneA))
    }


    if (!is.null(gnBaln)) {
      hasGeneB <- gnBaln %>%
        dplyr::filter(.data$AlnLength >= minaln) %>%
        dplyr::pull(.data$SubjectACC) %>% unique()

      res <- res %>%
        dplyr::mutate(AlignedToGeneB = (.data$ReadName %in% hasGeneB))
    }

    #---------------------
    # Add (optional) info about alignment to the host genome
    #---------------------
    if (!is.null(hostaln)) {
        #convert paf to GRanges
        pafgr <- paf2gr(hostaln, focus = "read", includetags=TRUE)
        # Keep only primary alignments
        pafgr <- pafgr[pafgr$tp == "P"]
        # If necessary filter on mapping quality
        # alignments with mapQ<10 are mostly secondary alignments, short alignments (e.g. 75% under 750bp) and alignments with <50% matches
        if (!is.na(minHost_mapQ) && length(minHost_mapQ)==1) {
          pafgr <- pafgr[pafgr$map_quality >= minHost_mapQ]
        }
        # Get the total length of the reads that is aligned to the host genome
        ALNreadlength <- sum(GenomicRanges::coverage(GenomicRanges::reduce(pafgr)))
        # Create a tibble
        ALNreadlength <- tibble::tibble(ReadName = names(ALNreadlength),
                                        HostAlign = TRUE,
                                        HostAlignedLength = ALNreadlength)
        # Join to the ReadLength table
        Raln <- dplyr::left_join(tibble::as_tibble(RL),
                                 ALNreadlength,
                                 by = "ReadName") %>%
          dplyr::mutate(HostAlign = ifelse(is.na(.data$HostAlign), FALSE, TRUE))
        # Calculate the percentage of read length that is aligned to the E coli genome
        Raln <- Raln %>%
          dplyr::mutate(HostAlignedPercent = 100 *
                          .data$HostAlignedLength /
                          .data$ReadLength)
        # Keep only relevant columns
        Raln <- Raln %>%
          dplyr::select(.data$ReadName, .data$HostAlign,
                        .data$HostAlignedLength, .data$HostAlignedPercent)
        # Get the best alignment for each read
        # Top quality, then top length, then top # of matches
        bestAlign <- hostaln %>%
          dplyr::group_by(.data$query_name) %>%
          dplyr::top_n(1, .data$map_quality) %>%
          dplyr::top_n(1, .data$map_length) %>%
          dplyr::top_n(1, .data$map_match) %>%
          dplyr::ungroup() %>%
          dplyr::distinct(.data$query_name, .keep_all = TRUE) %>%
          dplyr::mutate(BestHostAln_Readlength = .data$query_end - .data$query_start + 1) %>%
          dplyr::select(.data$query_name, .data$map_quality, .data$map_match,
                        .data$map_length, .data$BestHostAln_Readlength) %>%
          dplyr::rename("ReadName" = "query_name",
                        "BestHostAln_mapQ" =  "map_quality",
                        "BestHostAln_match" = "map_match",
                        "BestHostAln_length" = "map_length")

        # Merge this table with Raln
        Raln <- dplyr::left_join(Raln, bestAlign, by="ReadName")

        # merge with the result table
        res <- dplyr::left_join(res, Raln, by = "ReadName")
      }

  return(res)
}
