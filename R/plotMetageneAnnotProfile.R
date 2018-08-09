#' @title Plot metagene annotation profiles for different groups of genes.
#'
#' @description Plot metagene annotation profiles for different groups of genes.
#'
#' @param avgprofs A data frame with average profiles for different groups of genes (or "gene sets"). The object should contain the following columns:
#' \itemize{
#'   \item \code{Position}. Integer vector of genomic positions. The following rules should apply:
#'     \itemize{
#'       \item{The region upstream of the focus features has negative positions.}
#'       \item{The region downstream of the focus features has positive coordinates. The highest number in the \code{Position} column should correspond to the end of this region}
#'       \item{Positions of the TSS (beginning of focus features) and TES (end of focus features) should both be 0}
#'     }
#'   \item \code{Xcoord}. Numeric in [0,5]. Coordinates to use on the x axis. The following rules should apply:
#'     \itemize{
#'       \item{The \code{Xcoord} of the upstream region cover [0,2[}
#'       \item{The \code{Xcoord} of the TSS (beginning of focus feature) is 2}
#'       \item{The \code{Xcoord} of the body of the focus feature cover ]2,3[}
#'       \item{The \code{Xcoord} of the TES (end of focus feature) is 3}
#'       \item{The \code{Xcoord} of the downstream region cover ]3,5]}
#'     }
#'   \item \code{Strand}. Character vector of factor (with values "sense" or "antisense") for strand information.
#'   \item \code{Profile}. Numeric. The actual signal that is plotted.
#'   \item \code{Upper}. Numeric. Upper bound of the confidence interval for \code{AnnotationCoverage}.
#'   \item \code{Lower} Numeric. Lower bound of the confidence interval for \code{AnnotationCoverage}.
#'   \item \code{GeneSet} Character vector or factor defining the different gene sets.
#' }
#' @param gscolors Optional named character vector with colors per gene sets. Names of \code{gscolors} should match the levels of the \code{GeneSet} column.
#' @param usePercent Logical (defaults to FALSE) indicating if annotation coverages were expressed as presence/absence (and thus the \code{AnnotationCoverage} column is a percentage)
#' @param freeY Logical (defaults to FALSE) indicating if y-axis scales for sense and antisense strands should be different
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#' @import ggplot2
#'
#' @export
#'
#' @return A ggplot object
#'
#' @examples
#' ## Extract profiles around (+/-50bp) all genes.
#'   Prof <- annotationCoverageAroundFeatures(Genegr,
#'                                            sidedist = 50,
#'                                            usePercent = TRUE,
#'                                            nbins=3)
#' ## Assemble the profiles:
#'   Prof <- assembleProfiles(Prof)
#' ## Get the average profile for all genes on the sense and antisense strands
#'   avgALLsense = getAvgProfileWithCI(Prof$Profiles_Sense,
#'                                     pos = c(-50:0, 1:3, 0:50))
#'   avgALLantisense = getAvgProfileWithCI(Prof$Profiles_Antisense,
#'                                         pos = c(-50:0, 1:3, 0:50))
#' ## Select genes with a close neighbor on the same strand
#'   selGenes <- names(Genegr)[S4Vectors::mcols(
#'                   GenomicRanges::distanceToNearest(Genegr))$distance <= 5]
#' ## Get their average profiles
#'   avgSELsense = getAvgProfileWithCI(Prof$Profiles_Sense,
#'                                     selFeatures = selGenes,
#'                                     pos = c(-50:0, 1:3, 0:50))
#'   avgSELantisense = getAvgProfileWithCI(Prof$Profiles_Antisense,
#'                                         selFeatures = selGenes,
#'                                         pos = c(-50:0, 1:3, 0:50))
#' ## Define the xcoordinates:
#'   xcoord = c(seq(0, 2, length.out = 51),
#'              seq(2, 3, length.out = 5)[2:4],
#'              seq(3, 5, length.out = 51))
#' ## Assemble the data for plotting:
#'   avgprof <- rbind(avgALLsense, avgALLantisense,
#'                    avgSELsense, avgSELantisense)
#'   avgprof$Strand = rep(rep(c("sense", "antisense"),
#'                            each = 105),
#'                        times =2)
#'   avgprof$Xcoord = rep(xcoord, 4)
#'   avgprof$GeneSet = rep(c("All genes", "Close tandem neighbor"),
#'                         each = 210)
#' ## Finally plot the results:
#'   plotMetageneAnnotProfile(avgprof, usePercent = TRUE)
#' ## Focus on a region closer to the genes:
#'   plotMetageneAnnotProfile(avgprof, usePercent = TRUE) +
#'     ggplot2::coord_cartesian(xlim=c(1,4))
#'
#' @author Pascal GP Martin

plotMetageneAnnotProfile <- function(avgprofs,
                                     gscolors = NULL,
                                     usePercent = FALSE,
                                     freeY = FALSE) {

# Check arguments
  if (!all(c("Position", "Xcoord", "Strand",
             "Profile", "Upper", "Lower",
              "GeneSet") %in% colnames(avgprofs))) {
    stop("avgprofs should have the following columns: 'Position', 'Xcoord' 'Strand',
          'Profile', 'Upper', 'Lower', 'GeneSet'")
  }

  if (!missing(gscolors) && !is.null(gscolors) &&
        (is.null(names(gscolors)) ||
           !all(levels(factor(avgprofs$GeneSet)) %in% names(gscolors)))) {
    stop("names of gscolors should match the levels of the GeneSet column")
  }

  if (max(avgprofs$Xcoord)!=5 || min(avgprofs$Xcoord)!=0) {
    stop("Xcoord should be in [0,5]")
  }

# Fix strand column
  avgprofs <- avgprofs %>%
                dplyr::mutate(Strand = factor(tolower(.data$Strand),
                                              levels=c("sense", "antisense"),
                                              ordered = TRUE))


# Transform coverage if UsePercent
  if (usePercent) {
    if (max(avgprofs[avgprofs$Strand == "sense", "Profile"]) == 1) {
      avgprofs <- avgprofs %>%
                    dplyr::mutate(Profile = 100 * .data$Profile,
                                  Upper = 100 * .data$Upper,
                                  Lower = 100 * .data$Lower)
    }
  }


# Guess the labels for x coordinates
  maxxpos <- max(avgprofs[avgprofs$Strand=="sense","Position"])
  if (maxxpos %% 1000 == 0) {
    labx <- c(paste0("-", maxxpos/1000, "Kb"),
              paste0("-", maxxpos/2000, "Kb"),
              "TSS", "50%", "TES",
              paste0("+", maxxpos/2000, "Kb"),
              paste0("+", maxxpos/1000, "Kb"))
  } else {
    labx <- c(paste0("-", maxxpos, "bp"),
              paste0("-", maxxpos/2, "bp"),
              "TSS", "50%", "TES",
              paste0("+", maxxpos/2, "bp"),
              paste0("+", maxxpos, "bp"))
  }


#Plot
  pp <- avgprofs %>%
    ggplot2::ggplot(ggplot2::aes_string(x = "Xcoord",
                                 y = "Profile",
                                 group = "GeneSet")) +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = "Lower",
                                        ymax = "Upper",
                                        fill = "GeneSet"),
                           linetype=0,
                           alpha=0.2) +
      ggplot2::xlab("Genomic coordinate") +
      ggplot2::scale_x_continuous(breaks = c(0, 1, 2, 2.5, 3, 4, 5),
                                  labels = labx) +
      ggplot2::geom_vline(xintercept = c(2,3),
                          colour = "darkgrey",
                          linetype = "longdash") +
      ggplot2::geom_line(ggplot2::aes_string(x = "Xcoord",
                                             y = "Profile",
                                             color = "GeneSet"),
                                             lwd = 1) +
      ggplot2::theme_bw() +
      ggplot2::theme(text = ggplot2::element_text(size=16, color="black"),
                     axis.text.x = ggplot2::element_text(size=14),
                     axis.text.y = ggplot2::element_text(size=14),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.grid.major.x = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(size=16, face="bold", hjust=0.5),
                     strip.background = ggplot2::element_blank(),
                     strip.text = ggplot2::element_text(size=16, color="black")) +
      ggplot2::facet_grid(Strand~.,
                          scales = ifelse(freeY, "free_y", "fixed"))


  if (!missing(gscolors) && !is.null(gscolors)) {
    pp <- pp +
      ggplot2::scale_colour_manual(values = gscolors,
                                   name = "Gene Set") +
      ggplot2::scale_fill_manual(values = gscolors,
                                 name = "Gene Set")
  } else {
    pp <- pp +
      ggplot2::scale_colour_discrete(name = "Gene Set") +
      ggplot2::scale_fill_discrete(name = "Gene Set")
  }



  if (usePercent) {
    pp <- pp + ggplot2::ylab("% bases covered by annotations")
  } else {
    pp <- pp + ggplot2::ylab("Annotation coverage")
  }

  return(pp)

}
