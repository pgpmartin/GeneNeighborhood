#' @title Extend the presence value (1) of the closest single-point annotation to subsequent positions.
#'
#' @description Single-point annotations are annotations that cover only 1 base in the genome.
#'              For example, a gene border such as the TSS or TES.
#'              To study the neighborhood of genes/features in terms of single-point annotation,
#'              we first extract annotation profiles using \code{\link{annotationCoverageAroundFeatures}} (with \code{usePercent=TRUE}).
#'              Then we use this function \code{extendPointPresence} which, for each gene/feature,
#'              defines the location of the closest single-point annotation at each border and for each strand
#'              and extends its presence to all subsequent positions.
#'
#' @param annotprofiles An RleList similar to the one produced by the \code{\link{annotationCoverageAroundFeatures}} function.
#' @param sidedist an integer distance (in bp) around borders
#'
#' @export
#'
#' @return a list of 4 matrices with (extended) presence profiles at:
#'  \itemize{
#'  \item \code{upS}: the upstream border on the same strand ("sense")
#'  \item \code{upAS}: the upstream border on the opposite strand ("antisense")
#'  \item \code{dnS}: the downstream border on the same strand ("sense")
#'  \item \code{dnAS}: the downstream border on the opposite strand ("antisense")
#'  }
#'
#' @examples
#' ## Get the TSS of all genes:
#'   tss <- GenomicRanges::promoters(Genegr, upstream = 0, downstream = 1)
#' ## Get their (presence/absence) coverage around (+/-50bp) genes (with 3 bins/gene):
#'   tsscov <- annotationCoverageAroundFeatures(annot = tss,
#'                                              features = Genegr,
#'                                              sidedist = 50,
#'                                              usePercent = TRUE,
#'                                              nbins = 3)
#' ## For each gene, find the closest TSS and extend its presence to subsequent positions:
#'   extTSScov <- extendPointPresence(tsscov, sidedist=50)
#'
#' @author Pascal GP Martin

extendPointPresence <- function(annotprofiles,
                                sidedist = 2000L) {

# Check arguments
  if (!all(c("UpstreamBorder_Sense", "UpstreamBorder_Antisense",
             "DownstreamBorder_Sense", "DownstreamBorder_Antisense") %in% names(annotprofiles))) {
    stop("annotprofiles should have names as produced by the annotationCoverageAroundFeatures function")
  }

# Extract relevant matrices
  SideMatrices <- list(
    "upS" = RleList2matrix(annotprofiles[["UpstreamBorder_Sense"]])[,sidedist:1],
    "upAS" = RleList2matrix(annotprofiles[["UpstreamBorder_Antisense"]])[,sidedist:1],
    "dnS" = RleList2matrix(annotprofiles[["DownstreamBorder_Sense"]])[,(sidedist+2):(2*sidedist+1)],
    "dnAS" = RleList2matrix(annotprofiles[["DownstreamBorder_Antisense"]])[,(sidedist+2):(2*sidedist+1)])

# Check rownames
  stopifnot(all(sapply(lapply(SideMatrices, rownames),
                       identical,
                       rownames(SideMatrices[[1]]))))

# Check that usePercent has been used (i.e. that matrices are binary)
  if (!all(sapply(SideMatrices, isBinaryMat))) {
    stop("The data in annotprofiles is not binary.\nYou should run annotationCoverageAroundFeatures with usePercent=TRUE")
  }

# Number of features
  nfeat <- nrow(SideMatrices$upS)

# Extract profiles with non zero values on the range
  SideMatrices <- lapply(SideMatrices, function(mat) {mat[rowSums(mat)!=0,]})
  NumberProfiles <- sapply(SideMatrices, nrow)

# identify the first 1 in each row
  isPres <- lapply(SideMatrices, max.col, ties= "first") #Find the first 1 in each row
  ZeroTime <- lapply(isPres, `-`, 1) #Number of zeros
  OneTime <- lapply(isPres, function(x) {sidedist-x+1}) #number of ones

# function to fill the rows with 0s and 1s
  fillRow <- function(zero, one) {
    rep(c(0,1), times=c(zero, one))
  }

# Obtain the "extended" profiles
  resMat <- list()
  for (i in 1:4) {
    resMat[[i]] <- t(mapply(fillRow, ZeroTime[[i]], OneTime[[i]]))
    rownames(resMat[[i]]) <- rownames(SideMatrices[[i]])
  }

# If there is no profile, we replace by NA
  resMat[NumberProfiles==0] <- NA

# Set names
  names(resMat) <- names(SideMatrices)

  return(resMat)

}

