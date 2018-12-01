#' @title plot the distribution of intergenic distances
#'
#' @description Plots the distribution of intergenic distances by side
#'    (upstream/downstream) and orientation (Same or Opposite Strand)
#'    for different sets of genes, using density or violin plot.
#'
#' @param distdf A data frame or tibble containing distances for different sets
#'               of genes (ideally defined in a "GeneSet" column)
#' @param groupcolumn A character string with the name of the column containing
#'                    the description of the gene sets
#' @param type A character string in \code{c("ridge", "violin", "jitterbox")}
#'             indicating the type of plot/geom to use
#' @param genesetcols An optional named character vector with colors for the
#'                    different gene sets
#'                    (names should correspond to levels of the groupcolumn)
#' @param newlabs An optional named character vector with new labels for the
#'                different gene sets
#'                (names should correspond to levels of the groupcolumn)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr rename mutate
#' @importFrom rlang .data !! enquo
#' @importFrom ggplot2 ggplot aes geom_vline
#'    scale_x_continuous scale_y_discrete labs theme_bw theme facet_grid
#'    scale_fill_hue scale_fill_manual element_text element_blank
#'    rel expand_scale
#' @importFrom ggridges geom_density_ridges
#' @importFrom ggstance geom_violinh geom_boxploth
#'
#' @export
#'
#' @return A ggplot object
#'
#' @seealso \code{\link[ggplot2]{ggplot2}},
#'          \code{\link[ggplot2]{ggplot}},
#'          \code{\link[ggstance]{geom_violinh}},
#'          \code{\link[ggstance]{geom_boxploth}},
#'          \code{\link[ggridges]{geom_density_ridges}}
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
#'                                        GeneNeighborhood = GeneNeighbors,
#'                                        nboot=1e3, CItype="perc",
#'                                        ncores = 2)
#'
#' ## Extract the intergenic distances for all non-overlapping genes:
#'   alldist <- analyzeNeighborsDistance(GeneList = names(Genegr),
#'                                       GeneNeighborhood = GeneNeighbors,
#'                                       DistriTest = FALSE,
#'                                       nboot=1e3, CItype="perc",
#'                                       ncores = 2)
#'
#' ## Select a set of genes with a short upstream distances:
#'   ### get the upstream distances:
#'     updist <- alldist$distances$Distance[
#'                   alldist$distances$Side=="Upstream"]
#'     names(updist) <- alldist$distances$GeneName[
#'                          alldist$distances$Side=="Upstream"]
#'   ### Define sampling probabilities inversely proportional to the distance:
#'     probs <- (max(updist) - updist) / sum(max(updist) - updist)
#'   ### Sample 100 genes with these probabilities:
#'     set.seed(1234)
#'     lessRandGenes <- sample(names(updist), 100, prob=probs)
#'
#' ## Extract intergenic distances and stats for this new set of genes:
#'   lessRandDist <- analyzeNeighborsDistance(GeneList = lessRandGenes,
#'                                            GeneNeighborhood = GeneNeighbors,
#'                                            nboot=1e3, CItype="perc",
#'                                            ncores = 2)
#'
#' ## Assemble the data for all genes, random genes and "less random" genes
#' ## in a single data frame:
#'   mydist <- rbind.data.frame(alldist$distances,
#'                              randDist$distances,
#'                              lessRandDist$distances)
#'   mydist$GeneSet <- rep(c("All Genes", "Random Genes", "Less Random genes"),
#'                         times = c(nrow(alldist$distances),
#'                                   nrow(randDist$distances),
#'                                   nrow(lessRandDist$distances)))
#' ## Finally, plot the distribution of distances (density plot):
#'   plotDistanceDistrib(mydist)
#' ## Using violin plots and specific colors and labels
#'   plotDistanceDistrib(mydist,
#'                       type = "violin",
#'                       genesetcols = c("All Genes" = "grey",
#'                                       "Random Genes" = "lightblue",
#'                                       "Less Random genes" = "pink"),
#'                       newlabs = c("All Genes" = "ALL",
#'                                   "Random Genes" = "Random",
#'                                   "Less Random genes" = "Close Upstream"))
#'  ## Adding jitter points and a boxplot under the density plot
#'   plotDistanceDistrib(mydist, type = "jitterbox")
#'
#' @author Pascal GP Martin

plotDistanceDistrib <- function(distdf,
                                  groupcolumn = "GeneSet",
                                  type = c("ridge", "violin", "jitterbox"),
                                  genesetcols=NULL, newlabs=NULL) {

. <- Distance <- GeneSet <- NULL #avoids check errors due to non std eval

type <- match.arg(type)

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
    limX <- c(1, maxX+1)
    lineXintercept <- c(0, 10^(1:maxlogX))
    breakX <- c(0, 10^(1:log10(maxX)))
    labsX <- c("1bp", ifelse(breakX[-1]<1e3,
                             breakX[-1],
                             ifelse(breakX[-1]<1e6,
                                    paste0(breakX[-1]/1e3, "Kb"),
                                    paste0(breakX[-1]/1e6, "Mb"))
                             ))


## Prepare the dataset
  distdf <- distdf %>%
    dplyr::mutate(Orientation = factor(gsub("Strand$",
                                            " Strand",
                                            distdf$Orientation),
                                       levels = c("Same Strand",
                                                  "Opposite Strand"),
                                       ordered = TRUE))

## Define geom
  if (type=="ridge") {
    geomFUN <- ggridges::geom_density_ridges(scale=.9, alpha = .8)
  }
  if (type=="violin") {
    geomFUN <- ggstance::geom_violinh(draw_quantiles = .5, alpha = .8)
  }
  if (type=="jitterbox") {
      geomFUN <- ggridges::geom_density_ridges(scale=.7,
                                    position = position_nudge(x = 0,
                                                              y = .1),
                                    alpha=.8)
  }

## Prepare the plot
pp <- ggplot2::ggplot(distdf,
                      ggplot2::aes(x = Distance + 1,
                                   y = GeneSet,
                                   fill = GeneSet)) +
      ggplot2::scale_x_continuous(limits = limX,
                                  breaks = breakX + 1,
                                  labels = labsX,
                                  trans = "log10") +
      ggplot2::geom_vline(xintercept = lineXintercept + 1,
                          linetype = "dashed",
                          color = "lightgrey") +
      geomFUN +
      ggplot2::labs(x="Distance to the neighbor gene (bp)",
                    y="Gene sets") +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(
                         size = ggplot2::rel(1.2),
                         face = "bold",
                         hjust = 0.5),
                     axis.text.x = ggplot2::element_text(
                         size = ggplot2::rel(1.2),
                         vjust = 1,
                         angle = 45,
                         hjust = 1),
                     axis.text.y = ggplot2::element_text(
                         size = ggplot2::rel(1.2),
                         vjust = 0.5),
                     strip.background = ggplot2::element_blank(),
                     strip.text = ggplot2::element_text(
                         size = ggplot2::rel(1.2),
                         color = "black"))

if (type == "jitterbox") {
    pp <- pp +
        ggplot2::geom_jitter(height = .1,
                             size = .5,
                             alpha = 0.6,
                             show.legend = FALSE,
                             col="#999999") +
        ggstance::geom_boxploth(width=0.1,
                                fatten = 1.5,
                                outlier.shape = NA,
                                show.legend = FALSE,
                                alpha = 0.3)
}

## Facet by Side and Orientation
pp <- pp + ggplot2::facet_grid(Orientation~Side)

## Set colors
if (is.null(genesetcols)) {
    pp <- pp + ggplot2::scale_fill_hue(guide = FALSE)
} else {
    if (is.null(names(genesetcols)) |
        !all(names(genesetcols) %in% levels(distdf$GeneSet))) {
        message("genesetcols should be a named vector with names
                corresponding to gene sets\n")
        pp <- pp + ggplot2::scale_fill_hue(guide = FALSE)
    } else {
        pp <- pp + ggplot2::scale_fill_manual(values = genesetcols,
                                              guide = FALSE)
    }
}

## Set y axis, change the labels with newlabs if necessary
if (type=="ridge") {
    expY <- c(0.2, 1)
}
if (type == "violin")
{
  expY <- c(0.6, 0.6)
}

if (type == "jitterbox")
{
    expY <- c(0.2, 1)
}

if (is.null(newlabs)) {
    pp <- pp + ggplot2::scale_y_discrete(limits = rev(levels(distdf$GeneSet)),
                                         breaks = levels(distdf$GeneSet),
                                         expand =  ggplot2::expand_scale(
                                           add = expY))
} else {
    if (is.null(names(newlabs)) |
        !all(names(newlabs) %in% levels(distdf$GeneSet))) {
        message("newlabs should be a named vector with names
                corresponding to gene sets\n")
        pp <- pp +
            ggplot2::scale_y_discrete(limits = rev(levels(distdf$GeneSet)),
                                      breaks = levels(distdf$GeneSet),
                                      expand = ggplot2::expand_scale(
                                        add = expY))
    } else {
    pp <- pp +
        ggplot2::scale_y_discrete(limits = rev(levels(distdf$GeneSet)),
                                  breaks = levels(distdf$GeneSet),
                                  labels = newlabs,
                                  expand = ggplot2::expand_scale(
                                    add = expY))
    }
}

## Return plot
return(pp)

}
