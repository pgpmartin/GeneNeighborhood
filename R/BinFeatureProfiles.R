#' Bin coverage profiles given as a matrix or RleList of profiles.
#'
#' @param FeatureProfiles A \code{\link{matrix}} (with features in row and genomic positions in column) or \code{\link{RleList}}, typically obtained with the \code{\link[GeneNeighborhood]{profcomp}} function.
#' @param nbins Integer. Number of bins to define along each feature. Default value is 100 bins.
#' @param binwidth Integer. Size of the bins to define on the features. Default value is NULL and 100 bins are used instead.
#' @param aggregFUN A character string (or a function name), such as "median", "mean", etc., with the name of a function used to aggregate the data within the bins. Defaults to "mean".
#'                  A user-defined aggregFUN should take an \code{\link{Rle}} as an argument and return a single value.
#' @param asMatrix Logical (default to FALSE). If TRUE, the binned profiles are returned as a matrix (only possible if featrures have the same number of bins).
#' @param ... Further arguments passed to aggregFUN.
#'
#' @importFrom IRanges Views IRanges tile viewMins viewMaxs viewMeans viewSums viewApply mean
#' @importFrom methods as
#' @importFrom S4Vectors elementNROWS width Rle
#' @importFrom stats var
#'
#' @export
#'
#' @return An \code{RleList} or a \code{matrix} (if \code{asMatrix} is TRUE) of binned profiles.
#'
#'@section DETAILS:
#'  The function removes features of length < nbins (or binwidth) and returns an message with the number of removed features.
#'  The remaining features are sliced into \code{nbins} (possibly different bin sizes for different features) or into bins of size \code{binwidth} (possibly different number of bins for different features).
#'  Then the data in each bins is aggregated using the aggregFUN function
#'
#' @seealso profcomp AnnotationCoverageAroundFeatures
#'
#' @examples
#' ## Obtain coverage for all genes:
#'   covr <- GenomicRanges::coverage(Genegr)
#' ## Select a (random) set of (100) genes
#'   set.seed(123)
#'   randGenes <- sample(names(Genegr), 100)
#' ## Windows of interests covering the gene bodies + 50bp on each side
#'   woi <- Genegr[randGenes] +50
#' ## Coverage on these windows of interest:
#'   profs <- profcomp(covr, woi)
#' ## Make 10 bins:
#'   BinFeatureProfiles(profs, nbins = 10)
#' ## Return a matrix instead:
#'   BinFeatureProfiles(profs[1:3], nbins=10, asMatrix = TRUE)
#' ## Make bins of size 20bp:
#'   BinFeatureProfiles(profs, binwidth=20)
#' ##use your own function to aggregate the data within the bins
#'   mymean <- function(x) {mean(as.numeric(x), trim=0.1)}
#'   BinFeatureProfiles(profs[1:3], aggregFUN=mymean, nbins=10, asMatrix = TRUE)
#'
#' @author Pascal GP Martin

BinFeatureProfiles <- function(FeatureProfiles,
                               nbins = 100L,
                               binwidth = NULL,
                               aggregFUN = "mean",
                               asMatrix = FALSE,
                               ...)
{
# Some controls on arguments
  if (!missing(nbins) && !missing(binwidth)) {
    message("Both nbins and binwidth were provided. Using nbins only")
  }

  if (missing(nbins) && missing(binwidth)) {
    message("Working with 100 bins (default). Note that different features may have different bin sizes\n")
  }

  if (!missing(nbins) && (length(nbins) != 1 || any(is.na(nbins)))) {
    stop("nbins should be a single integer value (or should not be provided to use binwidth instead)")
  }

  if (!missing(nbins) && (!is.numeric(nbins) || nbins!=as.integer(nbins))) {
    stop("nbins must be an integer (or should not be provided to use binwidth instead)")
  }

  if (!missing(binwidth) && (length(binwidth) != 1 || any(is.na(binwidth)))) {
    stop("binwidth should be a single integer value (or should not be provided to use nbins instead)")
  }

    if (!missing(binwidth) && (!is.numeric(binwidth) || binwidth!=as.integer(binwidth))) {
    stop("binwidth must be an integer (or should not be provided to use nbins instead)")
    }

  if (!is.matrix(FeatureProfiles) && !is(FeatureProfiles,"RleList")) {
    stop("FeatureProfiles should be a matrix or an RleList")
  }

  if (!is.function(aggregFUN) && !is.character(aggregFUN)) {
    stop("aggregFUN must be a function or a character string with a function name")
  } else {
    FUN <- match.fun(aggregFUN)
  }

  if (length(FUN(S4Vectors::Rle(1:5), ...))!=1 || !is.numeric(FUN(S4Vectors::Rle(1:5), ...))) {
    stop("aggregFUN should return a single numeric value")
  }

# If FeatureProfile is a matrix, convert to an RleList
  if (is.matrix(FeatureProfiles)){
    FeatureProfiles <- methods::as(apply(FeatureProfiles, 1, Rle), "RleList")
  }

# Define the type of binning: FixedNumBin=TRUE (nbins is given) or FALSE (binwith is given)
  FixedNumBin <- TRUE
  if (!missing(binwidth) && missing(nbins)){
    FixedNumBin <- FALSE
  }

# Remove Features not compatible with nbins or binwith
  FeatureLengths = S4Vectors::elementNROWS(FeatureProfiles)

  if (FixedNumBin){
    FeatureIsOK <- FeatureLengths >= nbins
    if (sum(!FeatureIsOK) > 0) {
      message("Removing ", sum(!FeatureIsOK), " features of size lower than ", nbins, "bp")
    }
  } else {
    FeatureIsOK <- FeatureLengths >= binwidth
    if (sum(!FeatureIsOK) > 0) {
      message("Removing ", sum(!FeatureIsOK), " features of size lower than ", binwidth, "bp")
    }
  }

  FeatureProfiles <- FeatureProfiles[FeatureIsOK]
  FeatureLengths <- FeatureLengths[FeatureIsOK]

# Create the tiling
  if (FixedNumBin){

    tileFeatures = IRanges::Views(FeatureProfiles,
                                  IRanges::tile(IRanges::IRanges(start = 1,
                                                                 end = FeatureLengths),
                                                n = nbins))
    avgtilewidth <- IRanges::mean(S4Vectors::width(tileFeatures))
    message("Bin size is: ",
            round(mean(avgtilewidth), 2),
            " +/- ",
            round(sd(avgtilewidth), 2),
            "bp (mean +/- sd)")

  } else {

    tileFeatures=Views(FeatureProfiles,
                       tile(IRanges(start = 1,
                                    end = FeatureLengths),
                            width = binwidth))
    message("Average number of bins is: ",
            round(mean(S4Vectors::elementNROWS(S4Vectors::width(tileFeatures))), 2),
            " bins")
  }


# Perform the binning:

  ## If aggregFUN is "min", "max", "sum" or "mean" (default), use the view* functions (faster than viewApply)
   ### Define if FUN is one of "min", "max", "sum" or "mean" (default) base function:
      isBaseFUN <- c("min" = identical(FUN, base::min),
                     "max" = identical(FUN, base::max),
                     "sum" = identical(FUN, base::sum),
                     "mean" = identical(FUN, base::mean))
    ### Capture the extra arguments for aggregFUN
      supArgs <- list(...)
    ### Evaluate if supArgs is compatible with view* functions
      isOKargs <- (length(supArgs)==0) || (length(supArgs)==1 & names(supArgs)=="na.rm")

  if (any(isBaseFUN) && isOKargs) {
    RemoveNA <- if (length(supArgs)==0) (FALSE) else (supArgs$na.rm)

    res <- base::switch(names(isBaseFUN)[which(isBaseFUN)],
                        min = methods::as(IRanges::viewMins(tileFeatures, na.rm = RemoveNA), "CompressedRleList"),
                        max = methods::as(IRanges::viewMaxs(tileFeatures, na.rm = RemoveNA), "CompressedRleList"),
                        sum = methods::as(IRanges::viewSums(tileFeatures, na.rm = RemoveNA), "CompressedRleList"),
                        mean = methods::as(IRanges::viewMeans(tileFeatures, na.rm = RemoveNA), "CompressedRleList"))
  } else {
  ## Otherwise, use aggregFUN -> FUN
    res <- methods::as(IRanges::viewApply(tileFeatures, FUN, ..., simplify = TRUE), "CompressedRleList")
  }


# Convert to a matrix if asMatrix=TRUE and if all binned features have the same length
  if (asMatrix) {
    if(stats::var(S4Vectors::elementNROWS(res))!=0) {
      message("Binned features have different lengths. Cannot return a matrix. Returning an RleList")
    } else {
      res <- matrix(as.numeric(unlist(res, use.names = FALSE)),
                    nrow = length(res),
                    byrow = TRUE,
                    dimnames = list(names(res), NULL))
    }
  }

  return(res)
}

