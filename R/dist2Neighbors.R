#' @title Extract distances to the upstream and downstream neighbors
#'
#' @description Extract distances to the upstream and downstream neighbors.
#'              This is essentially a wrapper around \code{getDistSide} with
#'              the possibility to use several gene sets and with
#'              extra checks on the gene set(s)
#'
#' @param GeneNeighborhood A \code{tibble} obtained with the
#'        \code{\link{getGeneNeighborhood}} function.
#'        Can also be any data.frame with the relevant information.
#' @param geneset Either a character vector of identifiers defining the genes
#'                of interest or a (named) list of such character vectors.
#' @param genesetName Character string. Name of the gene set (only used if
#'                    \code{geneset} is a character vector).
#'                    If \code{geneset} is a list, the names of its elements
#'                    are used as names of the gene sets.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select rename slice transmute mutate pull
#'             group_by summarise bind_rows one_of
#' @importFrom reshape2 melt
#' @importFrom rlang .data !!
#' @importFrom tibble as_tibble
#'
#' @export
#'
#' @return A tibble with the following columns:
#' \itemize{
#'   \item GeneName. ID of the focus gene.
#'   \item Neighbor. ID of the gene neighbor
#'   \item Orientation. Orientation of the neighbor (same or opposite strand).
#'   \item Side. Upstream or Downstream.
#'   \item Distance in bp
#'   \item GeneSet. Name of the gene set
#' }
#'
#' @seealso \code{\link{getGeneNeighborhood}}
#'          \code{\link{getDistSide}}
#'
#' @examples
#' ## Obtain gene neighborhood information:
#'   GeneNeighbors <- getGeneNeighborhood(Genegr)
#' ## Get a (random) set of (100) genes:
#'   set.seed(123)
#'   randGenes <- sample(names(Genegr), 100)
#' ## Extract distances to the upstream/downstream neighbors for this set:
#'   distForRandGenes <- dist2Neighbors(GeneNeighborhood = GeneNeighbors,
#'                                      geneset = randGenes,
#'                                      genesetName = "RandomGenes")
#' ## Or for this set of genes and for all genes (the 'gene universe'):
#'  myGeneSets <- list("RandomGenes" = randGenes,
#'                     "AllGenes" = GeneNeighbors$GeneName)
#'  distForGeneSets <- dist2Neighbors(GeneNeighbors,
#'                                    myGeneSets)


dist2Neighbors <- function(GeneNeighborhood = NULL,
                           geneset = NULL,
                           genesetName = "GeneSet") {

#-----------
# Check arguments
#-----------

## Check GeneNeighborhood
##-------
if (is.null(GeneNeighborhood)) {
    stop("GeneNeighborhood dataset should be provided")
}

### Check that required columns are present:
ColumnsNames <- c("GeneName",
                  "Upstream", "UpstreamClass", "UpstreamDistance",
                  "Downstream", "DownstreamClass", "DownstreamDistance")

if (!all(ColumnsNames %in% colnames(GeneNeighborhood))) {
    stop("The GeneNeighborhood object must contain the following columns:\n",
         paste(ColumnsNames, collapse = ", ")
    )
}

## Check geneset
##-------
### Take all genes if geneset is NULL
if (is.null(geneset)) {
    geneset <- as.character(GeneNeighborhood[,"GeneName"])
}

### Convert factor to character vector (with a message)
if (is.factor(geneset)) {
    geneset <- as.character(geneset)
    message("Converting geneset from factor to character")
}

### check that geneset is a character vector or a list
if (!is.character(geneset) & !is.list(geneset)) {
    stop("geneset should be a character vector or a list of character vectors")
}

### If geneset is a character vector:
if (is.character(geneset)) {

    #### Remove NA values if any
    if (any(is.na(geneset))) {
        geneset <- geneset[!is.na(geneset)]
        message("NA values in geneset are removed")
        if (length(geneset)==0) {
            stop("Length of geneset is 0 after removing NAs")
        }
    }

    #### Check that IDs in geneset match those in GeneNeighborhood
    ids <- intersect(geneset,
                     GeneNeighborhood %>% dplyr::pull(.data$GeneName))
    if (length(ids)==0) {
        stop("IDs in geneset do not match those in GeneNeighborhood")
    }

    #### Keep only the IDs that match those in GeneNeighborhood
    inilength <- length(geneset)
    message("Initial number of genes in geneset (NA removed): ",
            inilength)
    if (length(ids)!=length(geneset)) {
        geneset <- ids
        message(length(geneset)-inilength,
                " genes were removed from geneset because they are ",
                "not in GeneNeighborhood")
    }
}

### If geneset is a list
if (is.list(geneset)) {

    #### Convert factors to character vectors (with a message)
    gsisfactor <- sapply(geneset, is.factor)
    if (any(gsisfactor)) {
        geneset[gsisfactor] <- lapply(geneset[gsisfactor], as.character)
        message(sum(gsisfactor),
                " genesets are converted from factor to character")
    }

    #### Check that all elements are character vectors
    if (!all(sapply(geneset, is.character))) {
        stop("Gene sets are not all character vectors")
    }

    #### Remove NA values if any
    nas <- lapply(geneset, is.na)
    containsNAs <- sapply(nas, sum)!=0
    if (any(containsNAs)) {
        geneset <- lapply(geneset, function(x) {x[!is.na(x)]})
        message("genesets ",
                paste(which(containsNAs), collapse = ", "),
                " contain NA values that have been removed")
        newLengths <- sapply(geneset, length)
        if (any(newLengths==0)) {
            stop("genesets",
                 paste(which(newLengths==0), collapse = ", "),
                 " have length zero after removing NAs")
        }
    }

    #### Check that IDs in gene sets match those in GeneNeighborhood
    ids <- lapply(geneset,
                  intersect,
                  y = GeneNeighborhood %>% dplyr::pull(.data$GeneName))
    idsLength <- sapply(ids, length)
    if (any(idsLength==0)) {
        stop("IDs in gene sets ",
             paste(which(idsLength==0), collapse = ", "),
             "do not match those in GeneNeighborhood")
    }

    #### Make names if they are absent
    if (length(names(geneset))==0) {
        names(geneset) <- paste0("GeneSet", 1:length(geneset))
    }
    if (length(unique(names(geneset)))!=length(geneset)) {
        message("Names of gene sets are not unique. They will be replaced")
        names(geneset) <- paste0("GeneSet", 1:length(geneset))
    }

    #### Keep only the IDs that match those in GeneNeighborhood
    inilength <- sapply(geneset, length)
    message("Initial number of genes in gene sets (NA removed): \n",
            paste(paste0(names(geneset),
                         ": ",
                         inilength,
                         " genes"),
                  collapse="\n"))

    IDsNotInData <- inilength - idsLength
    if (any(IDsNotInData > 0)) {
        geneset[IDsNotInData > 0] <- ids[IDsNotInData > 0]
        message("Some genes were removed from the gene sets ",
                "because they were not in GeneNeighborhood")
    }

}

#-----------
#Extract distances to the upstream and downstream genes
#-----------

## For a single gene set:
if (is.character(geneset)) {
    res <- dplyr::bind_rows(
       getDistSide(GNN = GeneNeighborhood,
                   glist = geneset,
                   Side = "Upstream"),
       getDistSide(GNN = GeneNeighborhood,
                   glist = geneset,
                   Side = "Downstream")
            ) %>%
        dplyr::mutate("GeneSet" = genesetName)
}

## For multiple gene sets:
if (is.list(geneset)) {
res <- list()
    for (i in 1:length(geneset)) {
        message("\nGeneset: ", names(geneset)[i])

        res[[i]] <- dplyr::bind_rows(
            getDistSide(GNN = GeneNeighborhood,
                        glist = geneset[[i]],
                        Side = "Upstream"),
            getDistSide(GNN = GeneNeighborhood,
                        glist = geneset[[i]],
                        Side = "Downstream")
                )
        names(res)[i] <- names(geneset)[i]
    }

    res <- reshape2::melt(res,
                          measure.vars = "Distance",
                          value.name = "Distance") %>%
            dplyr::select(-.data$variable) %>%
            dplyr::rename("GeneSet" = "L1") %>%
            tibble::as_tibble()
}

#Adjust the Side and Orientation columns (for stat tables and plotting)
res <- res %>%
    dplyr::mutate(Side = factor(.data$Side,
                                levels=c("Upstream",
                                         "Downstream"),
                                ordered = TRUE),
                  Orientation = factor(.data$Orientation,
                                       levels = c("SameStrand",
                                                  "OppositeStrand"),
                                       ordered = TRUE))

return(res)

}
