#' @title Stacked barplot of gene orientations for upstream and
#'     downstream genes.
#'
#' @description Plots stacked bars of the percentage of each orientation for
#'     the upstream and downstream genes.
#'
#' @param orientationStats A data frame produced by the
#'     \code{\link{analyzeNeighborsOrientation}} function
#'
#' @importFrom ggplot2 ggplot aes_string geom_bar position_stack
#'     scale_fill_manual scale_y_continuous
#'     theme_bw theme element_text labs rel
#'
#' @export
#'
#' @return A \code{\link{ggplot}} object
#'
#' @details
#' The function only plots the following orientations, if present:
#' \itemize{
#'   \item SameStrand. The upstream/downstream gene is on the same strand
#'       as the focus gene
#'   \item OppositeStrand. The upstream/downstream gene is on
#'       the opposite strand of the focus gene
#'   \item OppositeOverlap. The upstream/downstream gene is on
#'       the opposite strand of the focus gene and overlaps with the focus gene
#'   \item SameOverlap. The upstream/downstream gene is on
#'       the same strand of the focus gene and overlaps with the focus gene
#'   \item other. Any other, more complex situation (e.g. multiple overlaps,
#'       focus gene contained within another gene, etc.)
#' }
#'
#' @seealso \code{\link{ggplot2}}, \code{\link{ggplot}}
#'
#' @examples
#' ## Obtain gene neighborhood information:
#'   GeneNeighbors <- getGeneNeighborhood(Genegr)
#' ## Define a (random) set of (100) genes:
#'   set.seed(123)
#'   randGenes <- sample(names(Genegr), 100)
#' ## Analyze the orientation of their neighbors:
#'   NOS <- analyzeNeighborsOrientation(randGenes, GeneNeighborhood = GeneNeighbors)
#' ## Plot the results:
#'   plotNeighborsOrientation(NOS)
#
#' @author Pascal GP Martin

plotNeighborsOrientation <- function(orientationStats) {

    statres <- orientationStats

    if (any(tolower(statres$Orientation)=="other")) {
    statres$Orientation <- factor(statres$Orientation,
                                  levels = c("SameStrand", "SameOverlap",
                                             "OppositeStrand", "OppositeOverlap",
                                             "other"),
                                  ordered = TRUE)
  } else {
    statres$Orientation <- factor(statres$Orientation,
                                  levels = c("SameStrand", "SameOverlap",
                                             "OppositeStrand", "OppositeOverlap"),
                                  ordered = TRUE)
  }

  statres$Side <- factor(statres$Side,
                         levels = c("Upstream", "Downstream"),
                         ordered = TRUE)

  pp <- ggplot2::ggplot(statres, ggplot2::aes_string(x = "Side",
                                                     y = "Percentage",
                                                     fill = "Orientation")) +
            ggplot2::geom_bar(stat = "identity",
                              position = ggplot2::position_stack(reverse = TRUE)) +
            ggplot2::scale_fill_manual(name = "Neighbor\norientation",
                                       values = c("SameStrand" = "#417dbe",
                                                  "SameOverlap" = "#8db1d8",
                                                  "OppositeStrand" = "#55aa55",
                                                  "other" = "#b3b3b3",
                                                  "OppositeOverlap" = "#99cc99"),
                                       label = c("SameStrand" = "Same Strand",
                                                 "SameOverlap" = "Same Overlap",
                                                 "OppositeStrand" = "Opposite Strand",
                                                 "other" = "other",
                                                 "OppositeOverlap" = "Opposite Overlap")) +
            ggplot2::scale_y_continuous(expand = c(0.01, 0), limits = c(0, 100)) +
            ggplot2::theme_bw() +
            ggplot2::theme(legend.position = "right",
                           plot.title = ggplot2::element_text(hjust = 0.5,
                                                              size = ggplot2::rel(1.5),
                                                              face = "bold"),
                           axis.title.y = ggplot2::element_text(size = ggplot2::rel(1.5)),
                           axis.text.y = ggplot2::element_text(size = ggplot2::rel(1.5)),
                           axis.text.x = ggplot2::element_text(size = ggplot2::rel(1.5)),
                           legend.title = ggplot2::element_text(size = ggplot2::rel(1.4),
                                                                face = "bold"),
                           legend.text = ggplot2::element_text(size = ggplot2::rel(1.1)),
                           legend.title.align = 0.5) +
            ggplot2::labs(x = "")

    return(pp)
}
