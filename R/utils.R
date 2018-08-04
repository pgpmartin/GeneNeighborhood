#' @title Extract windows of interest from a coverage.
#'
#' @description Extract profiles from a genome-wide coverage on a set of windows of interest. Reverses the profiles for windows that are on the '-' (minus) strand.
#'
#' @param covr RleList. A genome-wide coverage (typically obtained with the \code{\link[GenomicRanges]{coverage}} function).
#' @param woi A \code{\link{GRanges}}. Windows of interest, over which to extract the profiles.
#'
#' @importFrom GenomicRanges strand
#' @importFrom S4Vectors revElements
#' @importFrom GenomeInfoDb seqlevels
#'
#' @export
#'
#' @return An RleList with stranded profiles at woi
#'
#' @section DETAILS:
#' Note that names of covr and seqnames of gr must match
#'
#' @examples
#' ## Obtain coverage for all genes:
#'   covr <- GenomicRanges::coverage(Genegr)
#' ## Select a (random) set of (100) genes
#'   set.seed(123)
#'   randGenes <- sample(names(Genegr), 100)
#' ## Windows of interests cover the gene bodies + 50bp on each side
#'   woi <- Genegr[randGenes] + 50
#' ## Coverage on these windows of interest:
#'   profcomp(covr, woi)
#'
#'
#' @author Pascal GP Martin
#'
profcomp <- function(covr=NULL, woi=NULL) {

## Check objects class
  if (any(c(is.null(covr), is.null(woi)))) {
    stop("Both covr and woi arguments should be provided to the profcomp function")
  }
  if (!is(covr, "RleList")) {
    stop("covr should be an RleList")
  }
  if (!is(woi, "GRanges")) {
    stop("woi should be a GRanges")
  }

## Check that seqlevels(gr) and names(covr) match
  if (is.null(names(covr)) & length(covr)>1) {
    message()
  }
  if (length(intersect(GenomeInfoDb::seqlevels(woi),
                       names(covr))) == 0) {
    stop("Names of covr and seqlevels of woi don't match")
  }

## Check for undefined strand
  isStar <- as.logical(GenomicRanges::strand(woi)=="*")
  if (any(isStar)) {
    message("'*' ranges were treated as '+' ranges")
  }

## Extract the profiles
  prof = covr[woi]

## Reverse the profiles for features on the minus strand
  prof <- revElements(prof, as.logical(GenomicRanges::strand(woi)=='-'))

## Add names to the profiles:
  if (is.null(names(woi))) {
    message("woi has no names to propagate to the extracted profiles")
  } else {
    names(prof)=names(woi)
  }

  return(prof)
}


#' @title Extract the TES (=end of the gene) from a \code{\link{GRanges}}

#' @description Extract the TES (=end of the gene) from a \code{\link{GRanges}}
#'
#' @param gr GRanges.
#'
#' @importFrom GenomicRanges strand promoters
#'
#' @export
#'
#' @return A GRanges of the TES positions
#'
#' @section DETAILS:
#' The TES is defined as start(gr) when strand(gr)=="+" and as end(gr) when strand(gr)=="-"
#' The function considers that strand=="*" is equivalent to strand=="+"
#'
#' @examples
#' ## Get the TSS of genes in GeneGR:
#'   GenomicRanges::promoters(Genegr, upstream = 0, downstream = 1)
#' ## Get their TES:
#'   getTES(Genegr)
#'
#' @author Pascal GP Martin

getTES <- function(gr) {

## Check for undefined strand
  isStar <- as.logical(GenomicRanges::strand(gr)=="*")
  if (any(isStar)) {
    message("'*' ranges were treated as '+' ranges")
  }

## Reverse gr strand
  revgr <- gr
  GenomicRanges::strand(revgr) <- ifelse(as.logical(GenomicRanges::strand(gr)=="-"), "+", "-")

## Get TES coordinates
  tes <- GenomicRanges::promoters(revgr, upstream=0, downstream=1)

## Reverse back the strand
  GenomicRanges::strand(tes) <- GenomicRanges::strand(gr)

  return(tes)
}

#' @title Convert an \code{RleList} to a \code{matrix}
#'
#' @description Convert an \code{RleList} to a \code{matrix}. The elements of the \code{RleList} should all have the same length.
#'
#' @param rlelist An \code{\link{RleList}}
#'
#' @return a \code{\link{matrix}} with one row for each element of the \code{RleList}
#'
#' @importFrom IRanges elementNROWS
#' @importFrom stats var
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' myrlelist <- IRanges::RleList(
#'                  x = S4Vectors::Rle(sample.int(100, 5)),
#'                  y = S4Vectors::Rle(sample.int(100, 5)),
#'                  z = S4Vectors::Rle(sample.int(100, 5)))
#' RleList2matrix(myrlelist)
#'
#' @author Pascal GP Martin


RleList2matrix <- function(rlelist = NULL)
{
## Check argument
  if (is.null(rlelist)) {
    stop("Please provide an RleList to be converted")
  }

  if (!is(rlelist, "RleList")) {
    stop("rlelist should be an RleList")
  }

  if (stats::var(IRanges::elementNROWS(rlelist)) != 0) {
    stop("The elements of rlelist should all have the same length")
  }
##Convert to matrix
  matrix(as.numeric(unlist(rlelist, use.names=F)),
         nrow = length(rlelist),
         byrow = TRUE,
         dimnames = list(names(rlelist),
                         NULL))
}


#' @title Convert an \code{matrix} to an \code{RleList}
#'
#' @description Convert an \code{matrix} to an \code{RleList}. if mat is not a matrix, a conversion will be attempted.
#'
#'
#' @param mat A \code{\link{matrix}} to convert as an \code{\link{RleList}}
#'
#' @return an \code{\link{RleList}}
#'
#' @importFrom methods as
#' @importFrom S4Vectors Rle
#'
#' @export
#'
#' @examples
#' matrix2RleList(matrix(c(rep(1:3, 3), rep(4:6, 3)), nrow=3))
#'
#' @author Pascal GP Martin

matrix2RleList <- function(mat=NULL)
{
## Check argument
if (is.null(mat)) {
  stop("Please provide a matrix to be converted")
}

## If necessary convert mat
  if (!is.matrix(mat)) {
    mat <- as.matrix(mat)
    message("mat was converted to a matrix")
  }

## Convert to RleList
  as(apply(mat, 1 , Rle), "RleList")
}

