#' @title Extract windows of interest from a coverage.
#'
#' @description Extract profiles from a genome-wide coverage on a set
#'     of windows of interest. Reverses the profiles for windows that
#'     are on the '-' (minus) strand.
#'
#' @param covr RleList. A genome-wide coverage (typically obtained with
#'             the \code{\link[GenomicRanges]{coverage}} function).
#' @param woi A \code{\link{GRanges}}. Windows of interest,
#'             over which to extract the profiles.
#'
#' @importFrom GenomicRanges strand
#' @importFrom S4Vectors revElements
#' @importFrom GenomeInfoDb seqlevels
#'
#' @export
#'
#' @return An RleList with stranded profiles at woi
#'
#' @details
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
    stop("Both covr and woi arguments should be provided
         to the profcomp function")
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
#' @param gr A \code{\link{GRanges}}.
#'
#' @importFrom GenomicRanges strand promoters
#'
#' @export
#'
#' @return A GRanges of the TES positions
#'
#' @details
#' The TES is defined as start(gr) when strand(gr)=="+" and
#' as end(gr) when strand(gr)=="-"
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
GenomicRanges::strand(revgr) <- ifelse(
                                    as.logical(GenomicRanges::strand(gr)=="-"),
                                    "+",
                                    "-")

## Get TES coordinates
tes <- GenomicRanges::promoters(revgr, upstream=0, downstream=1)

## Reverse back the strand
GenomicRanges::strand(tes) <- GenomicRanges::strand(gr)

return(tes)
}



#' Extract a region around the TES (=end of the gene) from
#'     a \code{\link{GRanges}}
#'
#' @param gr A \code{\link{GRanges}}.
#' @param upstream An integer defining the distance upstream of the
#'                 TES to extract.
#' @param downstream An integer defining the distance downstream of
#'                   the TES to extract.
#' @param limitToTSS A logical (default to FALSE) indicating if ranges should
#'                   be trimmed to don't extend upstream of the TSS
#' @param trim A logical (default to FALSE) indicating if out-of-bound ranges
#'             should be trimmed (see ?`trim,GenomicRanges-method`)
#'
#' @importFrom GenomicRanges strand promoters trim start end
#'
#' @export
#'
#' @return A GRanges of regions around the TES positions.
#'
#' @details
#' The TES is defined as start(gr) when strand(gr)=="+" and
#' as end(gr) when strand(gr)=="-"
#' The function considers that strand=="*" is equivalent to strand=="+"
#' To extract 500bp on each side of the TES,
#' use \code{upstream=500} and \code{downstream=501}
#'
#'@seealso \code{\link[GenomicRanges:intra-range-methods]{promoters}},
#'         \code{\link[GenomicRanges:intra-range-methods]{flank}} and
#'         \code{\link[GenomicRanges:intra-range-methods]{trim}}
#'         in \code{\link[GenomicRanges:intra-range-methods]{intra-range-methods}}
#'
#' @examples
#' ## Get a 50bp region on each side of the TSS of genes in GeneGR:
#'   suppressWarnings(
#'    GenomicRanges::promoters(Genegr,
#'                               upstream = 50,
#'                               downstream = 51)
#'                   )
#' ## Same around the TES:
#'   suppressWarnings(
#'     getTESregion(Genegr, upstream = 50, downstream = 51)
#'     )
#' ## Do not extend regions beyond the TSS or Chr1 borders
#'   getTESregion(Genegr,
#'                upstream = 50,
#'                downstream = 51,
#'                limitToTSS = TRUE,
#'                trim = TRUE)
#'
#' @author Pascal GP Martin

getTESregion <- function(gr,
                         upstream=500L,
                         downstream=501L,
                         limitToTSS = FALSE,
                         trim = FALSE) {

## Check for undefined strand
isStar <- as.logical(GenomicRanges::strand(gr)=="*")
if (any(isStar)) {
    message("'*' ranges were treated as '+' ranges")
}

## Reverse gr strand
revgr <- gr
GenomicRanges::strand(revgr) <- ifelse(
                                    as.logical(GenomicRanges::strand(gr)=="-"),
                                    "+",
                                    "-")

## Get TES region coordinates
if (trim) {
    tesRegion <- suppressWarnings(
        GenomicRanges::promoters(revgr,
                                 upstream = downstream-1,
                                 downstream = upstream+1))
} else {
    tesRegion <- GenomicRanges::promoters(revgr,
                                          upstream = downstream-1,
                                          downstream = upstream+1)
}

## Reverse back the strand
GenomicRanges::strand(tesRegion) <- GenomicRanges::strand(gr)

## Correct the borders if limitToTSS is TRUE
if (limitToTSS) {
    GenomicRanges::start(tesRegion) <- ifelse(
        as.logical(GenomicRanges::strand(gr)=="+"),
        pmax(GenomicRanges::start(gr), GenomicRanges::start(tesRegion)),
        GenomicRanges::start(tesRegion))
    GenomicRanges::end(tesRegion) <- ifelse(
        as.logical(GenomicRanges::strand(gr)=="+"),
        GenomicRanges::end(tesRegion),
        pmin(GenomicRanges::end(gr), GenomicRanges::end(tesRegion)))
}

#Limit the region to chromosome borders:
if (trim) {
    tesRegion=GenomicRanges::trim(tesRegion)
}

return(tesRegion)
}


#' @title Convert an \code{RleList} to a \code{matrix}
#'
#' @description Convert an \code{RleList} to a \code{matrix}. The elements of
#'              the \code{RleList} should all have the same length.
#'
#' @param rlelist An \code{\link{RleList}}
#'
#' @return a \code{\link{matrix}} with one row for each element
#'         of the \code{RleList}
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
matrix(as.numeric(unlist(rlelist, use.names = FALSE)),
       nrow = length(rlelist),
       byrow = TRUE,
       dimnames = list(names(rlelist),
                       NULL))
}


#' @title Convert an \code{matrix} to an \code{RleList}
#'
#' @description Convert an \code{matrix} to an \code{RleList}.
#'              if mat is not a matrix, a conversion will be attempted.
#'
#'
#' @param mat A \code{\link{matrix}} to convert as an
#'            \code{\link[IRanges]{RleList}}
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


#' @title Test if a matrix is a binary matrix
#'
#' @description Test if a matrix is a binary matrix
#'              (i.e. contains only 0s and 1s, and possibly NA).
#'              Adapted from https://stackoverflow.com/a/23276062
#'
#' @param mat A matrix. If mat is a data frame it will
#'            be converted by data.matrix with a warning.
#'
#' @return A logical. TRUE if the matrix is binary, FALSE otherwise
#' @export
#'
#' @examples
#' set.seed(123)
#' x <- matrix(rnorm(100), ncol= 10) # not a binary matrix
#' y <- matrix(sample(0:1, 100, rep=TRUE), ncol=10) # a binary matrix
#' z <- y ; z[x>1.5] <- NA #a binary matrix with some NA
#' isBinaryMat(x) #FALSE
#' isBinaryMat(y) #TRUE
#' isBinaryMat(z) #TRUE

isBinaryMat <- function(mat) {
    stopifnot(is.matrix(mat) | is.data.frame(mat))
    if (is.data.frame(mat)) {
        mat <- data.matrix(mat)
        warning("mat is a data.frame.
                A conversion was attempted using data.matrix")
    }
    return(identical(as.numeric(mat), as.numeric(as.logical(mat))))
}

