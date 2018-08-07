#' @title Assemble metagene profiles from the results of the annotationCoverageAroundFeatures function
#'
#' @description Assemble metagene/metafeature profiles from the results of the annotationCoverageAroundFeatures function.   \cr
#'              A metagene profile is defined as:
#' \enumerate{
#'   \item A region upstream of the feature ending with the upstream border of the feature (typically the TSS)
#'   \item A region covering the body of the feature (possibly binned to normalize for different features' sizes)
#'   \item A region downstream of the feature starting with the downstream border of the feature (typically the TES)
#'   }
#'
#' @param annotcovr a list with 6 \code{\link[IRanges]{RleList}}s named:
#' \itemize{
#'   \item \code{Feature_Sense}: Coverage (possibly binned) on the body of features, on the sense strand (i.e. the strand of the feature)
#'   \item \code{Feature_Antisense}: Coverage (possibly binned) on the body of features, on the antisense strand (i.e. the opposite strand of features)
#'   \item \code{UpstreamBorder_Sense}: Coverage centered on the upstream border of features (Typically the TSS), on the sense strand
#'   \item \code{UpstreamBorder_Antisense}: Coverage centered on the the upstream border of features, on the antisense strand
#'   \item \code{DownstreamBorder_Sense}: Coverage centered on the downstream border of features (Typically the TES), on the sense strand
#'   \item \code{DownstreamBorder_Antisense}: Coverage centered on the the downstream border of features, on the antisense strand
#' }
#' as produced by the \code{\link{annotationCoverageAroundFeatures}} function
#' @param sideDist Integer (optional). Distance to keep on each side of the features. Defaults to the maximum distance availabe in \code{annotcovr}
#'
#' @return a list of 2 RleLists corresponding to:
#'  \itemize{
#'   \item \code{Metaprofiles_Sense}: RleList of metagene profiles for the sense strand
#'   \item \code{Metaprofiles_Antisense}: RleList of metagene profiles for the antisense strand
#'  }
#'
#' @export
#'
#' @seealso \code{\link{annotationCoverageAroundFeatures}}
#'
#' @examples
#' ## Extract the profiles around (+/-50bp) the first 10 genes. We bin the genes in 3 bins only.
#'   top10Prof <- annotationCoverageAroundFeatures(Genegr, features=1:10,
#'                                                 sidedist = 50,
#'                                                 usePercent = TRUE,
#'                                                 nbins=3)
#' ## Assemble the metagene profiles
#'   metaProf <- assembleMetagene(top10Prof)

assembleMetagene <- function(annotcovr, sideDist=NULL) {

 #Test arguments
if (!all(c(is.list(annotcovr),
           length(annotcovr) == 6,
           names(annotcovr) %in% c("Feature_Sense", "Feature_Antisense",
                                   "UpstreamBorder_Sense", "UpstreamBorder_Antisense",
                                   "DownstreamBorder_Sense", "DownstreamBorder_Antisense")))) {
  stop("annotcovr should be a list of RleLists as produced by the annotationCoverageAroundFeatures function")
}

if (!is.null(sideDist) && !is.numeric(sideDist)) {
  stop("sideDist should be an integer")
}

if (!is.null(sideDist) && is.finite(sideDist) && sideDist < 0) {
  sideDist = -sideDist
}

 # Convert RleLists to matrices
  annotcovr <- lapply(annotcovr, RleList2matrix)

 # Adjust sideDist
  sideDistMat <- (ncol(annotcovr$UpstreamBorder_Sense)-1)/2
  sideDistArg <- sideDist
  sideDist <- min(c(sideDistArg, sideDistMat), na.rm=TRUE)

 # Sense strand
  res_sense <- cbind(annotcovr$UpstreamBorder_Sense[,(sideDistMat+1-sideDist):(sideDistMat+1)],
                     annotcovr$Feature_Sense,
                     annotcovr$DownstreamBorder_Sense[,(sideDistMat+1):(sideDistMat+1+sideDist)])
  names(res_sense) <- rownames(annotcovr$Feature_Sense)
  res_sense <- matrix2RleList(res_sense)

 # Antisense strand
  res_antisense <- cbind(annotcovr$UpstreamBorder_Antisense[,(sideDistMat+1-sideDist):(sideDistMat+1)],
                         annotcovr$Feature_Antisense,
                         annotcovr$DownstreamBorder_Antisense[,(sideDistMat+1):(sideDistMat+1+sideDist)])
  names(res_antisense) <- rownames(annotcovr$Feature_Antisense)
  res_antisense <- matrix2RleList(res_antisense)

  return(list("Metaprofiles_Sense" = res_sense,
              "Metaprofiles_Antisense" = res_antisense))
  }
