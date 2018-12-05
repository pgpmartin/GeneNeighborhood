#' @title Filter a gene set and extract upstream or downstream distances
#'        do their neighbors
#' @description Extract distances to the upstream or downstream neighbors
#'
#' @param GNN A \code{tibble} obtained with the
#'            \code{\link{getGeneNeighborhood}} function.
#'            Can also be any data.frame with the relevant information.
#' @param glist A character vector of gene IDs (must match GNN$GeneName)
#' @param Side One of 'Upstream' or 'Downstream'
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter select rename slice mutate transmute pull
#'             group_by summarise bind_rows one_of
#' @importFrom rlang .data !!
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
#' }
#'
#' @seealso \code{\link{getGeneNeighborhood}}
#'          \code{\link{dist2Neighbors}}
#'
#' @examples
#' ## Obtain gene neighborhood information:
#'   GeneNeighbors <- getGeneNeighborhood(Genegr)
#' ## Get a (random) set of (10) genes:
#'   set.seed(123)
#'   randGenes <- sample(names(Genegr), 10)
#' ## Extract distances to the upstream neighbors for this set:
#' getDistSide(GeneNeighbors, randGenes, "Upstream")
#' ## And for the downstream genes:
#' getDistSide(GeneNeighbors, randGenes, "Downstream")

getDistSide <- function(GNN,
                        glist,
                        Side = c("Upstream", "Downstream")) {

    Side <- match.arg(Side)
    className <- paste0(Side, "Class")
    DistName <- paste0(Side, "Distance")

    selectedColumns <- c("GeneName",
                         Side,
                         className,
                         DistName)

    GNNdata <- GNN %>%
        dplyr::select(dplyr::one_of(selectedColumns)) %>%
        dplyr::rename("SideGene" = !!Side,
                      "SideClass" = !!className,
                      "Distance" = !!DistName)

    #Remove gene that have no upstream/downstream data (i.e. class=="other")
    hasSide <- GNNdata %>%
        dplyr::slice(match(glist, .data$GeneName)) %>%
        dplyr::transmute(hasSide = !is.na(.data$SideGene) &
                             !is.na(.data$Distance)) %>%
        dplyr::pull(hasSide)

    if (any(!hasSide)) {
        message(sum(!hasSide), paste0(" genes with undefined ",
                                      tolower(Side),
                                      " gene are removed from Geneset"))
    }

    glistSel <- glist[hasSide]

    ##Remove genes with overlapping upstream/downstream gene
    isOVL <- GNNdata %>%
        dplyr::slice(match(glistSel, .data$GeneName)) %>%
        dplyr::transmute(isOVL = grepl("overlap",
                                       tolower(.data$SideClass)) &
                             .data$Distance == 0) %>%
        dplyr::pull(isOVL)

    if (any(isOVL)) {
        message(sum(isOVL), paste0(" genes with an overlapping ",
                                   tolower(Side),
                                   " gene are removed from Geneset"))
        glistSel <- glistSel[!isOVL]
    }

    ##Count the remaining genes
    ngSel <- length(glistSel) #Number of genes with Upstream or Downstream info
    cat("Final number of genes with",
        tolower(Side),
        "gene data:",
        ngSel,
        "genes\n")

    ##Extract distances
    res <- GNNdata %>%
        dplyr::filter(.data$GeneName %in% glistSel) %>%
        dplyr::rename("Neighbor" = "SideGene", "Orientation" = "SideClass") %>%
        dplyr::mutate("Side" = Side)

    return(res)
}
