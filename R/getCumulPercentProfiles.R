#' @title For gene sets, obtain the cumulative percentage of genes with single point annotations around their borders.
#'
#' @description For different sets of genes defined in \code{genesets},
#'              calculate the cumulative percentage of genes that have at least
#'              one single-point annotation when the distance from the gene borders increases.
#'
#' @param extMat A list with 4 matrices named upS, upAS, dnS and dnAS.
#'               Typically obtained with the \code{\link{extendPointPresence}} function
#' @param genesets a list of gene sets (character vectors) for which cumulative percentages should be obtained
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr rename mutate select
#' @importFrom reshape2 melt
#' @importFrom rlang .data
#'
#' @export
#'
#' @return a data frame with columns:
#'  \itemize{
#'    \item \code{Percent}: cumulative percentage when the distance from the gene borders increases
#'    \item \code{GeneSet}: Name of the gene sets. Taken from \code{names(genesets)} if not NULL
#'    \item \code{Position}: Integer giving the genomic position (Negative values for upstream, positive for downstream)
#'    \item \code{Side}: Character string indicating "Upstream" or "Downstream"
#'    \item \code{Strand}: Character string indicating "sense" or "antisense"
#'  }
#'
#' @seealso \code{\link{extendPointPresence}}
#'
#' @examples
#' ## Get the TES of all genes:
#'   tes <- getTES(Genegr)
#' ## Get their (presence/absence) coverage around (+/-50bp) genes (with 3 bins/gene):
#'   tescov <- annotationCoverageAroundFeatures(annot = tes,
#'                                              features = Genegr,
#'                                              sidedist = 50,
#'                                              usePercent = TRUE,
#'                                              nbins = 3)
#' ## For each gene, find the closest TES and extend its presence to subsequent positions:
#'   extTEScov <- extendPointPresence(tescov, sidedist=50)
#' ## Identify genes with a TES at less than 6bp from their TES, on any strand
#'   require(GenomicRanges)
#'   Dist2TES <- mcols(distanceToNearest(tes, ignore.strand = TRUE))$distance
#'   CloseTES <- names(Genegr)[Dist2TES<6]
#'   NotCloseTES <- names(Genegr)[Dist2TES>=6] #Also get the complementary set
#' ## Calculate the cumulative percentage as the distance increases:
#'   CP <- getCumulPercentProfiles(extTEScov,
#'                                 list("All" = names(Genegr),
#'                                      "CloseTES" = CloseTES,
#'                                      "FarTES" = NotCloseTES))
#'
#' @author Pascal GP Martin


getCumulPercentProfiles <- function(extMat,
                                 genesets = NULL) {

# Check extMat
  if (!is.list(extMat) || !identical(names(extMat), c("upS", "upAS", "dnS", "dnAS"))) {
    stop("extMat should be a list with elements upS, upAS, dnS and dnAS")
  }

  if (!all(sapply(extMat, is.matrix))) {
    stop("extMat should contains 4 matrices")
  }

  isNAextMat <- sapply(extMat, function(x) all(is.na(x)))
  if (all(isNAextMat)) {
    stop("extMat only has NA values")
  }

  ncolextMat <- sapply(extMat[!isNAextMat], ncol)
  if (stats::var(ncolextMat)!=0) {
    stop("All matrices in extMat should have the same number of columns")
  }
  NCOLS <- ncolextMat[1]

# function to calculate the percentages per column)
  percFUN <- function(x) {100 * colMeans(x, na.rm = TRUE)}
##TODO: Add the upper and lower bounds of the CI for the percentage

# if genesets is null, we take all genes
  if (is.null(genesets)) {

    res <- list()
    res[[1]] <- list()
    for (i in 1:4) {
      if (isNAextMat[i]) {
        res[[1]][[i]] <- NA
      } else {
        res[[1]][[i]] <- percFUN(extMat[[i]])
      }
    }
    names(res) <- "all"

    #otherwise we use genesets
  } else {

    #if genesets is a vector, convert to a list of length 1
    if (!is.list(genesets)) {
      stopifnot(is.vector(genesets))
      genesets <- list(genesets)
    }

    #Get rownames of extMat
    rn <- lapply(extMat, rownames)

    #Calculate the percentages:
    res <- list()
    for (i in 1:length(genesets)) {
      ## Define genesInMat
      if (length(genesets)==1 && is.na(genesets[[i]])) {
        genesInMat <- rn
        warning("Gene Set #", i, " is NA. Taking all available genes for this gene set")
      } else {
        genesInMat <- lapply(rn, intersect, genesets[[i]][!is.na(genesets[[i]])])
      }
      ## Get the number of selected genes for each matrix
      ng <- sapply(genesInMat, length)
      ## For each matrix, calculate the colMeans
      res[[i]] <- list()
      for (j in 1:4) {
        if (ng[j] == 0) {
          res[[i]][[j]] <- rep(0, NCOLS)
          warning("Gene Set #", i, " has zero intersection with element ", names(extMat)[j], " of extMat")
        } else {
          res[[i]][[j]] <- percFUN(extMat[[j]][genesInMat[[j]],])
        }
      }
      ## Name the elements
      names(res[[i]]) <- names(extMat)
    }
    names(res) <- if (!is.null(names(genesets))) (names(genesets)) else (paste0("GeneSet", 1:length(genesets)))
  }

  #Get the number of groups
  ngrp <- ifelse(is.null(genesets), 1, length(genesets))

  #Get the distance on each side of the borders
  sidedist <- NCOLS

  #create the genopos vector
  genopos <- c((-1):(-sidedist),
               (-1):(-sidedist),
               1:sidedist,
               1:sidedist)


  #Format the results
  res <- reshape2::melt(res) %>%
    dplyr::mutate(Position = rep(genopos, ngrp),
                  Side = factor(ifelse(substring(.data$L2,1,2)=="up",
                                       "Upstream",
                                       "Downstream"),
                                levels=c("Upstream", "Downstream"),
                                ordered = TRUE),
                  Strand = factor(ifelse(substring(.data$L2,3,4)=="S",
                                         "sense",
                                         "antisense"),
                                  levels=c("sense", "antisense"),
                                  ordered = TRUE) ) %>%
    dplyr::select(-.data$L2) %>%
    dplyr::rename("Percent" = "value",
                  "GeneSet" = "L1")


  #return the result
  return(res)

}

