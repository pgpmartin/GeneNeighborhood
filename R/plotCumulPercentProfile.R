#' @title Plot the cumulative proportions obtained from \code{\link{getCumulPercentProfiles}}
#'
#' @description Plot the cumulative proportions obtained from \code{\link{getCumulPercentProfiles}}
#'
#' @param cumulprofs A data frame with cumulative percentage profiles for different gene sets. As produced by the \code{\link{getCumulPercentProfiles}} function.
#' @param gscolors Optional named character vector with colors per gene sets
#' @param freeY Logical (defaults to FALSE) indicating if y-axis scales for panels/facets should be different
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr rename mutate
#' @importFrom rlang .data
#' @import ggplot2
#'
#' @export
#'
#' @return A ggplot object
#'
#' @seealso \code{\link{getCumulPercentProfiles}}
#'
#'#' @examples
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
#' ## Plot
#'   plotCumulPercentProfile(CP)
#'
#' @author Pascal GP Martin

plotCumulPercentProfile <- function(cumulprofs,
                                    gscolors=NULL,
                                    freeY = FALSE) {

# Check arguments
  if (!all(c("Percent", "GeneSet", "Position", "Side", "Strand") %in% colnames(cumulprofs))) {
    stop("cumulprofs should have the following columns: 'Percent', 'GeneSet', 'Position', 'Side', 'Strand'\n
         See the ")
  }

  if (!missing(gscolors) && !is.null(gscolors) &&
      (is.null(names(gscolors)) ||
       !all(levels(factor(cumulprofs$GeneSet)) %in% names(gscolors)))) {
    stop("names of gscolors should match the levels of the GeneSet column")
  }


# Fix Strand and Side columns
  cumulprofs <- cumulprofs %>%
    dplyr::mutate(Strand = factor(tolower(.data$Strand),
                                  levels=c("sense", "antisense"),
                                  ordered = TRUE),
                  Side = factor(.data$Side,
                                levels = c("Upstream", "Downstream"),
                                ordered = TRUE))


  #Define the labels for x coordinates
  maxxpos <- max(cumulprofs$Position)
  if (maxxpos %% 1000 == 0) {
    labx <- c(paste0("-", maxxpos/1000, "Kb"),
              paste0("-", maxxpos/2000, "Kb"),
              "Border",
              paste0("+", maxxpos/2000, "Kb"),
              paste0("-", maxxpos/1000, "Kb"))
    atx <- c(-maxxpos,
             -maxxpos/2,
             0,
             maxxpos/2,
             maxxpos)
  } else {
    labx <- c(paste0("-", maxxpos, "bp"),
              paste0("-", maxxpos/2, "bp"),
              "Border",
              paste0("+", maxxpos/2, "bp"),
              paste0("+", maxxpos, "bp"))
    atx <- c(-maxxpos,
             -maxxpos/2,
             0,
             maxxpos/2,
             maxxpos)
  }



  #Plot
  pp <- cumulprofs %>%
            ggplot2::ggplot(ggplot2::aes_string(x = "Position",
                                                y = "Percent",
                                                group = "GeneSet")) +
            ggplot2::xlab("Maximum distance from gene border") +
            ggplot2::ylab("Cumulative percentage of genes") +
            ggplot2::scale_x_continuous(breaks = atx,
                                        labels = labx) +
            ggplot2::geom_vline(xintercept = c(-0.1,0.1),
                                color = "darkgrey",
                                linetype = "longdash") +
            ggplot2::geom_line(ggplot2::aes_string(x = "Position",
                                                   y = "Percent",
                                                   color = "GeneSet")) +
            ggplot2::theme_bw() +
            ggplot2::theme(text = ggplot2::element_text(size=16, color="black"),
                           axis.text.x = ggplot2::element_text(size=14, angle=90, hjust=1, vjust=0.5),
                           axis.text.y = ggplot2::element_text(size=14),
                           panel.grid.minor = ggplot2::element_blank(),
                           panel.grid.major.x = ggplot2::element_blank(),
                           plot.title = ggplot2::element_text(size=16, face="bold", hjust=0.5),
                           strip.background = ggplot2::element_blank(),
                           strip.text = ggplot2::element_text(size=16, color="black")) +
            ggplot2::facet_grid(Strand~Side,
                                scales = ifelse(freeY, "free", "free_x"),
                                space = "free_x")


  if (!missing(gscolors) && !is.null(gscolors)) {
    pp <- pp +
            ggplot2::scale_colour_manual(values=gscolors,
                                         name = "Gene Set")
  } else {
    pp <- pp +
            ggplot2::scale_colour_discrete(name = "Gene Set")
  }

  return(pp)

}
