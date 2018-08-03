#' @title plot the distribution of intergenic distances
#'
#' @description Plots the densities of intergenic distances by side (upstream/downstream) and orientation (Same or Opposite Strand) for different sets of genes.
#'
#' @param distdf A data frame or tibble containing distances for different sets of gene (ideally defined in a "GeneSet" column)
#' @param groupcolumn A character string with the name of the column containing the description of the gene sets
#' @param genesetcols An optional named character vector with colors for the different gene sets (names should correspond to levels of the groupcolumn)
#' @param newlabs An optional named character vector with new labels for the different gene sets (names should correspond to levels of the groupcolumn)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr rename mutate
#' @importFrom rlang .data !! enquo
#' @importFrom ggplot2 ggplot aes geom_vline scale_x_continuous scale_y_discrete labs theme_bw theme facet_grid scale_fill_hue scale_fill_manual element_text element_blank
#' @importFrom ggridges geom_density_ridges
#'
#' @export
#'
#' @return A ggplot object
#'
#' @examples
#' ## Obtain gene neighborhood information:
#'   GeneNeighbors <- getGeneNeighborhood(Genegr)
#'
#' ## Get a (random) set of (100) genes:
#'   set.seed(123)
#'   randGenes <- sample(names(Genegr), 100)
#'
#' ## Extract intergenic distances and stats for these genes:
#'   randDist <- analyzeNeighborsDistance(GeneList = randGenes,
#'                                        GeneNeighborhood = GeneNeighbors)
#'
#' ## Extract the intergenic distances for all non-overlapping genes:
#'   alldist <- analyzeNeighborsDistance(GeneList = names(Genegr),
#'                                       GeneNeighborhood = GeneNeighbors,
#'                                       DistriTest = FALSE)
#'
#' ## Select a set of genes with a short upstream distances:
#'   ### get the upstream distances:
#'     updist <- alldist$distances$Distance[alldist$distances$Side=="Upstream"]
#'     names(updist) <- alldist$distances$GeneName[alldist$distances$Side=="Upstream"]
#'   ### Define sampling probabilities inversely proportional to the distance:
#'     probs <- (max(updist) - updist) / sum(max(updist) - updist)
#'   ### Sample 100 genes with these probabilities:
#'     set.seed(1234)
#'     lessRandGenes <- sample(names(updist), 100, prob=probs)
#'
#' ## Extract intergenic distances and stats for this new set of genes:
#'   lessRandDist <- analyzeNeighborsDistance(GeneList = lessRandGenes,
#'                                            GeneNeighborhood = GeneNeighbors)
#'
#' ## Assemble the data for all genes, random genes and "less random" genes in a single data frame:
#'   mydist <- rbind.data.frame(alldist$distances,
#'                              randDist$distances,
#'                              lessRandDist$distances)
#'   mydist$GeneSet <- rep(c("All Genes", "Random Genes", "Less Random genes"),
#'                         times = c(nrow(alldist$distances),
#'                                   nrow(randDist$distances),
#'                                   nrow(lessRandDist$distances)))
#' ## Finally, plot the distribution:
#'   plotDistanceDensities(mydist)
#'
#' @author Pascal GP Martin

plotDistanceDensities <- function(distdf, groupcolumn = "GeneSet", genesetcols=NULL, newlabs=NULL) {

  . <- Distance <- GeneSet <- NULL #avoids check errors due to non standard eval

  ## Rename the column defining gene sets if necessary
  if (groupcolumn != "GeneSet") {
    groupcolumn <- rlang::enquo(groupcolumn)
    distdf <- distdf %>% dplyr::rename("GeneSet" = !!groupcolumn)
  }

  ## Convert as factor if necessary
  if (!is.factor(distdf$GeneSet)) {
    distdf <- distdf %>% dplyr::mutate(GeneSet = factor(.data$GeneSet))
  }


  ## Define the limits of the x axis
    maxlogX <- ceiling(log10(max(distdf$Distance+1)))
    ### We keep the largest log interval that contains at least 1% of the data:
        maxX <- distdf %>%
               dplyr::mutate(grp = cut(.data$Distance + 1,
                                       breaks = c(0, 10^(1:maxlogX)))) %>%
               dplyr::count(.data$grp) %>%
               dplyr::mutate(prop = 100*.data$n / sum(.data$n)) %>%
               dplyr::slice(max(which(.data$prop > 1))) %>%
               dplyr::pull(.data$grp) %>%
               as.character() %>%
               gsub(pattern="^\\(.+,|]$", "", .) %>%
               as.numeric()
    ### Define limit, breaks and labels for the X axis:
  limX <- c(1, maxX)
  breakX <- c(0, 10^(1:log10(maxX)))
  labsX <- c("1bp", ifelse(breakX[-1]<1e3,
                           breakX[-1],
                           ifelse(breakX[-1]<1e6,
                                  paste0(breakX[-1]/1e3, "Kb"),
                                  paste0(breakX[-1]/1e6, "Mb"))
                           ))



  ## Prepare the plot
    pp <- distdf %>%
    dplyr::mutate(Orientation = factor(gsub("Strand$", " Strand", distdf$Orientation),
                                       levels = c("Same Strand", "Opposite Strand"),
                                       ordered = TRUE)) %>%
    ggplot2::ggplot(ggplot2::aes(x=Distance+1, y=GeneSet, fill=GeneSet)) +
      ggplot2::geom_vline(xintercept=breakX+1,
                          linetype="dashed",
                          color="lightgrey") +
      ggridges::geom_density_ridges(scale=0.9) +
      ggplot2::scale_x_continuous(limits=limX,
                                  breaks=breakX+1,
                                  labels=labsX,
                                  trans="log10") +
      ggplot2::labs(x="Distance to the neighbor gene (bp)",
                    y="Gene sets") +
      ggplot2::theme_bw() +
      ggplot2::theme(text = ggplot2::element_text(size = 16),
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
                     axis.text.x = ggplot2::element_text(size = 16, vjust = 1, angle = 45, hjust = 1),
                     axis.text.y = ggplot2::element_text(size = 16, vjust = 0.5),
                     strip.background = ggplot2::element_blank(),
                     strip.text = ggplot2::element_text(size = 16, color = "black")) +
      ggplot2::facet_grid(Orientation~Side)

  ## Set colors
  if (is.null(genesetcols)) {
    pp <- pp + ggplot2::scale_fill_hue(guide = FALSE)
  } else {
    if (is.null(names(genesetcols)) | !all(names(genesetcols) %in% levels(distdf$GeneSet))) {
      message("genesetcols should be a named vector with names corresponding to gene sets\n")
      pp <- pp + ggplot2::scale_fill_hue(guide = FALSE)
    } else {
      pp <- pp + ggplot2::scale_fill_manual(values = genesetcols,
                                            guide = FALSE)
    }
  }

## Set y axis, change the labels with newlabs if necessary
if (is.null(newlabs)) {
    pp <- pp + ggplot2::scale_y_discrete(limits = rev(levels(distdf$GeneSet)),
                                         breaks = levels(distdf$GeneSet),
                                         expand = c(0.02, 0))
} else {
    if (is.null(names(newlabs)) | !all(names(newlabs) %in% levels(distdf$GeneSet))) {
        message("newlabs should be a named vector with names corresponding to gene sets\n")
        pp <- pp + ggplot2::scale_y_discrete(limits = rev(levels(distdf$GeneSet)),
                                             breaks = levels(distdf$GeneSet),
                                             expand = c(0.02, 0))
    } else {
    pp <- pp + ggplot2::scale_y_discrete(limits = rev(levels(distdf$GeneSet)),
                                         breaks = levels(distdf$GeneSet),
                                         labels = newlabs,
                                         expand = c(0.02, 0))
    }
}

## Return plot
return(pp)

}
