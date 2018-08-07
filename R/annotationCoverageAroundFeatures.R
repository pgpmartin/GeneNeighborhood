#' Extract stranded annotation profiles at and around a set of features/intervals.
#'
#' @param annot A \code{\link{GRanges}} of (stranded) genomic annotations
#' @param features An optional \code{GRanges} or a character vector with a subset of annot names to extract profiles at. Defaults to NULL for all annot features.
#' @param sidedist Integer. Distance on each side of features at which to extract the profiles
#' @param usePercent Logical. Should annotation coverage be expressed as presence/absence (0/1)? (instead of counting the number of features aligning at each position)
#' @param BinFeatures Logical. Should \code{features} be binned? Default to TRUE and will be forced to TRUE if \code{features} have different width.
#' @param ... additional arguments passed to \code{\link{BinFeatureProfiles}}
#'
#' @importFrom GenomicRanges strand coverage reduce width trim promoters
#' @importFrom S4Vectors Rle
#'
#' @export
#'
#'@section DETAILS:
#' The \code{usePercent} option is intended to give a percentage of bases covered by annotations when performing average profiles from several features. \cr
#' When calculating a coverage with usePercent=TRUE any value > 1 is replaced by 1. \cr
#' This is equivalent to calculating \code{coverage(reduce(annot))} instead of \code{coverage(annot)}.
#'
#' @return
#' A list with the following elements:
#' \itemize{
#'   \item \code{Feature_Sense}: Coverage of annot (possibly binned) at features, on the sense strand (i.e. the strand of the feature)
#'   \item \code{Feature_Antisense}: Coverage of annot (possibly binned) at features, on the antisense strand (i.e. the opposite strand of features)
#'   \item \code{UpstreamBorder_Sense}: Coverage of annot around (+/- sidedist) the upstream border of features (Typically the TSS), on the sense strand
#'   \item \code{UpstreamBorder_Antisense}: Coverage of annot around (+/- sidedist) the upstream border of features, on the antisense strand
#'   \item \code{DownstreamBorder_Sense}: Coverage of annot around (+/- sidedist) the downstream border of features (Typically the TES), on the sense strand
#'   \item \code{DownstreamBorder_Antisense}: Coverage of annot around (+/- sidedist) the downstream border of features, on the antisense strand
#' }
#'
#' @seealso \code{\link{BinFeatureProfiles}}
#'
#' @examples
#' ## Extract the profiles around (+/-50bp) all genes. We bin the genes in 3 bins only.
#' AllProf <- annotationCoverageAroundFeatures(Genegr, sidedist = 50, usePercent = TRUE, nbins=3)
#'
#' @author Pascal GP Martin

annotationCoverageAroundFeatures <- function(annot,
                                             features = NULL,
                                             sidedist = 2000L,
                                             usePercent = FALSE,
                                             BinFeatures = TRUE,
                                             ...) {

  #tests
  if (!is(annot, "GRanges")) {
    stop("annot should be a GRanges")
  }

  #Remove annotations with undefined strands:
  undefinedStrand <- as.logical(GenomicRanges::strand(annot)=="*")
  if (any(undefinedStrand)) {
    annot <- annot[!undefinedStrand]
    message(sum(undefinedStrand), "annotations with undefined strand are removed")
  }

  #Get annotation coverages for each strand
if (usePercent) {
  covr_plus <- GenomicRanges::coverage(GenomicRanges::reduce(annot[GenomicRanges::strand(annot)=="+"]))
  covr_minus <- GenomicRanges::coverage(GenomicRanges::reduce(annot[GenomicRanges::strand(annot)=="-"]))
} else {
  covr_plus <- GenomicRanges::coverage(annot[GenomicRanges::strand(annot)=="+"])
  covr_minus <- GenomicRanges::coverage(annot[GenomicRanges::strand(annot)=="-"])
}


  #Define windows of interests (woi) based on features
  if (is.null(features)) {
    features <- 1:length(annot)
  }

  if (is(features, "GRanges")) {
    undefStrandFeature <- as.logical(GenomicRanges::strand(features)=="*")
    if (any(undefStrandFeature)) {
      woi <- features[!undefStrandFeature]
      message(sum(undefStrandFeature), " features with undefined strand are removed")
    } else {
      woi <- features
    }
  } else {
    woi <- annot[features]
    undefStrandwoi <- as.logical(GenomicRanges::strand(woi)=="*")
    if (any(undefStrandwoi)) {
      woi <- woi[!undefStrandwoi]
      message(sum(undefStrandwoi), " features with undefined strand are removed")
    }
  }

  #Add names to woi if necessary
  if (is.null(names(woi))) {
    names(woi) <- paste("woi", 1:length(woi), sep="_")
  }

  #Define ranges around woi borders
  ## Ranges around upstream border
  upwoi <- suppressWarnings(GenomicRanges::promoters(woi,
                                                     upstream = sidedist,
                                                     downstream = sidedist+1))
  ## Ranges around downstream border
  downwoi <- suppressWarnings(getTESregion(woi,
                                           upstream = sidedist,
                                           downstream = sidedist+1,
                                           limitToTSS = FALSE,
                                           trim = FALSE))

  #Remove ranges that exceed chromosome borders
  outOfChrom <- (GenomicRanges::width(upwoi) > GenomicRanges::width(GenomicRanges::trim(upwoi))) |
    (GenomicRanges::width(downwoi) > GenomicRanges::width(GenomicRanges::trim(downwoi)))

  if (any(outOfChrom)) {
    message(sum(outOfChrom), " windows exceeding chromosome borders are removed")
    woi <- woi[!outOfChrom]
    upwoi <- upwoi[!outOfChrom]
    downwoi <- downwoi[!outOfChrom]
  }


  #Extract profiles at woi and around woi borders
  prof_woi_Plus <- profcomp(covr_plus, woi)
  prof_woi_Minus <- profcomp(covr_minus, woi)

  prof_upwoi_Plus <- profcomp(covr_plus, upwoi)
  prof_upwoi_Minus <- profcomp(covr_minus, upwoi)

  prof_dnwoi_Plus <- profcomp(covr_plus, downwoi)
  prof_dnwoi_Minus <- profcomp(covr_minus, downwoi)


  #Combine Plus/Minus strands as sense/antisense strands
  stwoi <- as.character(GenomicRanges::strand(woi))
  nmwoi <- names(woi)

  ## WOI
  prof_woi_sense <- append(prof_woi_Plus[stwoi=="+"],
                           prof_woi_Minus[stwoi=="-"])
  prof_woi_sense <- prof_woi_sense[nmwoi]

  prof_woi_antisense <- append(prof_woi_Minus[stwoi=="+"],
                               prof_woi_Plus[stwoi=="-"])
  prof_woi_antisense <- prof_woi_antisense[nmwoi]

  ## Upstream border of WOI
  prof_upwoi_sense <- append(prof_upwoi_Plus[stwoi=="+"],
                             prof_upwoi_Minus[stwoi=="-"])
  prof_upwoi_sense <- prof_upwoi_sense[nmwoi]

  prof_upwoi_antisense <- append(prof_upwoi_Minus[stwoi=="+"],
                                 prof_upwoi_Plus[stwoi=="-"])
  prof_upwoi_antisense <- prof_upwoi_antisense[nmwoi]

  ## Downstream border of WOI
  prof_dnwoi_sense <- append(prof_dnwoi_Plus[stwoi=="+"],
                             prof_dnwoi_Minus[stwoi=="-"])
  prof_dnwoi_sense <- prof_dnwoi_sense[nmwoi]

  prof_dnwoi_antisense <- append(prof_dnwoi_Minus[stwoi=="+"],
                                 prof_dnwoi_Plus[stwoi=="-"])
  prof_dnwoi_antisense <- prof_dnwoi_antisense[nmwoi]


  # if woi have different width or if BinFeatures=TRUE, then use BinFeatureProfiles on feature bodies:
  featureSize <- GenomicRanges::width(woi)
  isVariableSize <- var(featureSize)!=0
  if (isVariableSize) {
    BinFeatures <- TRUE
  }

  if (BinFeatures) {
    message("Features are binned using BinFeatureProfiles")
    ##Define the kind of binning performed (nbins or binwidth) and if isNBins, the number of bins used
    supArgs <- list(...)
    if (length(supArgs)==0) {
      isNBins = TRUE
      minFeatureSize = 100L
    } else {
      ArgMatch <- pmatch(names(supArgs), c("nbins", "binwidth"))
      if (any(ArgMatch==1)) {
        isNBins = TRUE
        names(supArgs)[ArgMatch==1] <- "nbins"
        minFeatureSize = supArgs$nbins
      } else {
        if (any(ArgMatch==2)) {
          isNBins = FALSE
        } else {
          isNBins = TRUE
          minFeatureSize = 100L
        }
      }
    }


    ## Get the sense strand profile
      ### Note that if usePercent is TRUE AND features is a subset of annot, then the sense strand profile is always == 1
    if (usePercent && !is(features, "GRanges") && isNBins) {
      isOKFeatureSize <- featureSize>=minFeatureSize
      nn <- names(prof_woi_sense)[isOKFeatureSize]
      prof_woi_sense <- as(rep(list(S4Vectors::Rle(rep(1, minFeatureSize))),
                               sum(isOKFeatureSize)),
                           "CompressedRleList")
      names(prof_woi_sense) <- nn
    } else {
      prof_woi_sense <- BinFeatureProfiles(prof_woi_sense, ...)
    }

    ## Get the antisense strand profile
    prof_woi_antisense <- BinFeatureProfiles(prof_woi_antisense, ...)

    #Features of width<NBins are removed so we filter and order the upstream and downstream profiles accordingly
    nmwoi <- names(prof_woi_sense)
    prof_upwoi_sense <- prof_upwoi_sense[nmwoi]
    prof_upwoi_antisense <- prof_upwoi_antisense[nmwoi]
    prof_dnwoi_sense <- prof_dnwoi_sense[nmwoi]
    prof_dnwoi_antisense <- prof_dnwoi_antisense[nmwoi]

  }

  res <- list("Feature_Sense" = prof_woi_sense,
              "Feature_Antisense" = prof_woi_antisense,
              "UpstreamBorder_Sense" = prof_upwoi_sense,
              "UpstreamBorder_Antisense" = prof_upwoi_antisense,
              "DownstreamBorder_Sense" = prof_dnwoi_sense,
              "DownstreamBorder_Antisense" = prof_dnwoi_antisense)

  return(res)

}


