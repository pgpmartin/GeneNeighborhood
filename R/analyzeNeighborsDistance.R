#' @title Analyze intergenic distances
#'
#' @description Extracts the upstream OR downstream intergenic distances and compares their distribution to a reference universe.
#'
#' @param GeneNeighborhood A \code{tibble} obtained with the \code{\link{getGeneNeighborhood}} function. Can also be any data.frame with the relevant information.
#'
#' @param glist A character vector of genes of interest
#' @param guniv A character vector of genes in the universe.
#' @param distest A logical (default to TRUE) indicating if genes in glist should be compared to the Universe
#' @param Side A character string indicating if the distances to the "Upstream" or "Downstream" genes should be analyzed
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom dplyr select rename slice transmute pull group_by summarise bind_rows
#' @importFrom coin pvalue independence_test
#' @importFrom stats quantile median sd ks.test wilcox.test pnorm
#'
#' @export
#'
#' @return A list containing:
#' \itemize{
#'   \item distances. All the intergenic distances for glist genes
#'   \item stats. The statistics on the distances
#' }
#'
#' @section DETAILS:
#' Note that the function removes overlapping up/downstream genes but not adjacent genes (for which distance = 0).
#' When  following \code{distest} is \code{TRUE}, the following tests are performed:
#' \itemize{
#'   \item Kolmogorov-Smirnov test. See \code{\link{ks.test}}. Although not adapted to integer values that generate ties, this test is generally consistent with the other two.
#'   \item Mann-Whithney U test. See \code{\link{wilcox.test}}.
#'   \item Independence test. See \code{\link{independence_test}} in the \code{coin} package
#' }
#'
#' @examples
#' ## Obtain gene neighborhood information:
#'   GeneNeighbors <- getGeneNeighborhood(Genegr)
#' ## Get a (random) set of (100) genes:
#'   set.seed(123)
#'   randGenes <- sample(names(Genegr), 100)
#' ## Extract their upstream distances and compare to non selected genes:
#'   statDistanceSide(GeneNeighborhood = GeneNeighbors,
#'                    glist = randGenes,
#'                    Side = "Upstream")
#'
#' @author Pascal GP Martin
#'

statDistanceSide <- function(GeneNeighborhood = NULL,
                             glist = NULL,
                             guniv = NULL,
                             distest = TRUE,
                             Side=c("Upstream", "Downstream")) {

#-----------
# Check arguments
#-----------

if (is.null(GeneNeighborhood)) {
  stop("GeneNeighborhood dataset should be provided")
}

if (Side == "Upstream") {
  if (!all(c("GeneName", "Upstream", "UpstreamClass", "UpstreamDistance") %in% colnames(GeneNeighborhood))) {
    stop("The GeneNeighborhood object must contain columns 'GeneName', 'Upstream', 'UpstreamClass' and 'UpstreamDistance'")
  }
} else {
  if (Side == "Downstream") {
    if (!all(c("GeneName", "Downstream", "DownstreamClass", "DownstreamDistance") %in% colnames(GeneNeighborhood))) {
      stop("The GeneNeighborhood object must contain columns 'GeneName', 'Downstream', 'DownstreamClass' and 'DownstreamDistance'")
    }
  } else {
    stop("Side must be either 'Upstream' or 'Downstream'")
  }
}


if (is.null(glist)) {
  guniv = NULL
  distest = FALSE
  glist <- as.character(GeneNeighborhood$GeneName)
} else {
  if (is.null(guniv) && distest == TRUE) {
    guniv <- as.character(GeneNeighborhood$GeneName)
  }
  if (length(intersect(glist, GeneNeighborhood$GeneName)) == 0) {
    stop("No intersection between the gene set and the GeneName column in GeneNeighborhood")
  }
}


  #Select data based on Side
  if (Side == "Upstream") {
        GNNdata <- GeneNeighborhood %>%
      dplyr::select(.data$GeneName, .data$Upstream, .data$UpstreamClass, .data$UpstreamDistance) %>%
      dplyr::rename("SideGene" = "Upstream", "SideClass" = "UpstreamClass", "Distance" = "UpstreamDistance")
  }

    if (Side == "Downstream") {
      GNNdata <- GeneNeighborhood %>%
        dplyr::select(.data$GeneName, .data$Downstream, .data$DownstreamClass, .data$DownstreamDistance) %>%
        dplyr::rename("SideGene" = "Downstream", "SideClass" = "DownstreamClass", "Distance" = "DownstreamDistance")
    }


  #Start the analysis
  cat("\nAnalysis of Distances to", tolower(Side), "gene:\n")
  cat("========================================\n")

  #Remove gene that have no upstream/downstream data (i.e. class=="other")
  hasSide <- GNNdata %>%
    dplyr::slice(match(glist, .data$GeneName)) %>%
    dplyr::transmute(hasSide = !is.na(.data$SideGene) & !is.na(.data$Distance)) %>%
    dplyr::pull(hasSide)

  if (any(!hasSide)) {
    cat(sum(!hasSide), paste0("genes with undefined ", tolower(Side), " gene are removed\n"))
  }

  glistSel <- glist[hasSide]

  ##Remove genes with overlapping upstream/downstream gene
  isOVL <- GNNdata %>%
    dplyr::slice(match(glistSel, .data$GeneName)) %>%
    dplyr::transmute(isOVL = grepl("overlap", tolower(.data$SideClass)) & .data$Distance == 0) %>%
    dplyr::pull(isOVL)

  if (any(isOVL)) {
    cat(sum(isOVL), paste0("genes with an overlapping ", tolower(Side), " gene are removed from GeneList\n"))
    glistSel <- glistSel[!isOVL]
  }


  ##Count the remaining genes
  ngSel <- length(glistSel) #Number of genes with UPSTREAM/DOWNSTREAM info
  cat("GeneList for", tolower(Side), "gene analysis has", ngSel, "genes\n")

  ##Extract distances
  distSel <- GNNdata %>%
    dplyr::filter(.data$GeneName %in% glistSel)

  #Summarize the distances
  ## 95% of the values are contained with [Lower95:Upper95]
  distsum <- distSel %>%
    dplyr::group_by(.data$SideClass) %>%
    dplyr::summarise(n = sum(!is.na(.data$Distance)),
                     Min = min(.data$Distance, na.rm=TRUE),
                     Lower95 = quantile(.data$Distance, 0.025, na.rm=TRUE),
                     Q1 = quantile(.data$Distance, 0.25, na.rm=TRUE),
                     Median = median(.data$Distance, na.rm=TRUE),
                     Mean = mean(.data$Distance, na.rm=TRUE),
                     Q3 = quantile(.data$Distance, 0.75, na.rm=TRUE),
                     Upper95 = quantile(.data$Distance, 0.975, na.rm=TRUE),
                     Max = max(.data$Distance, na.rm=TRUE),
                     SD = sd(.data$Distance, na.rm=TRUE),
                     SEM = sd(.data$Distance, na.rm=TRUE)/sqrt(sum(!is.na(.data$Distance)))
    ) %>%
    as.data.frame(stringsAsFactors = FALSE)

  #Compare to universe
  if (distest) {

    #Remove from universe genes without a defined upstream/downstream gene
    univhasSide <- GNNdata %>%
      dplyr::slice(match(guniv, .data$GeneName)) %>%
      dplyr::transmute(hasSide = !is.na(.data$SideGene) & !is.na(.data$Distance)) %>%
      dplyr::pull(hasSide)

    if (any(!univhasSide)) {
      cat(sum(!univhasSide), paste0("genes from GeneUniverse with undefined ", tolower(Side), " gene are removed\n"))
    }

    gunivSel <- guniv[univhasSide]

    #Remove genes with overlapping upstream/downstream gene
    univIsOVL <- GNNdata %>%
      dplyr::slice(match(gunivSel, .data$GeneName)) %>%
      dplyr::transmute(isOVL = grepl("overlap", tolower(.data$SideClass)) & .data$Distance == 0) %>%
      dplyr::pull(isOVL)

    if (any(univIsOVL)) {
      cat(sum(univIsOVL), paste0("genes from GeneUniverse with an overlapping ", tolower(Side), " gene are removed\n"))
      gunivSel <- gunivSel[!univIsOVL]
    }


    #Count the remaining genes
    ngunivSel <- length(gunivSel)
    cat("Universe for", tolower(Side), "gene analysis has", ngunivSel, "genes\n")

    #Get distances for Universe
    distSelUniv <- GNNdata %>%
                       dplyr::filter(.data$GeneName %in% gunivSel)

    #Summarize the distances for the Universe
    ## 95% of the values are contained within [Lower95:Upper95]
    distsumUniv <- distSelUniv %>%
      dplyr::group_by(.data$SideClass) %>%
      dplyr::summarise(n = sum(!is.na(.data$Distance)),
                       Min = min(.data$Distance, na.rm=TRUE),
                       Lower95 = quantile(.data$Distance, 0.025, na.rm=TRUE),
                       Q1 = quantile(.data$Distance, 0.25, na.rm=TRUE),
                       Median = median(.data$Distance, na.rm=TRUE),
                       Mean = mean(.data$Distance, na.rm=TRUE),
                       Q3 = quantile(.data$Distance, 0.75, na.rm=TRUE),
                       Upper95 = quantile(.data$Distance, 0.975, na.rm=TRUE),
                       Max = max(.data$Distance, na.rm=TRUE),
                       SD = sd(.data$Distance, na.rm=TRUE),
                       SEM = sd(.data$Distance, na.rm=TRUE)/sqrt(sum(!is.na(.data$Distance)))
      ) %>%
      as.data.frame(stringsAsFactors = FALSE)

    #Merge the summary datasets
    nr <- nrow(distsum)
    distsum <- dplyr::bind_rows(distsum, distsumUniv)
    distsum$GeneGroup <- factor(rep(c("GeneList", "GeneUniverse"),
                                    times=c(nr, nrow(distsumUniv))),
                                levels=c("GeneList", "GeneUniverse"), ordered=T)

    #Add statistics to the table
    distsum$KS.pvalue <- rep(NA, nrow(distsum))
    distsum$Wilcox.pvalue <- rep(NA, nrow(distsum))
    distsum$Independ.pvalue <- rep(NA, nrow(distsum))

    for (i in 1:nr) {
      oriclass <- as.character(distsum$SideClass[i])
      testSet <- distSel %>%
        dplyr::filter(.data$SideClass == oriclass)

      univSet <- distSelUniv %>%
        dplyr::filter(.data$SideClass==oriclass)

      isInTestSet <- factor(ifelse(
        (univSet %>% dplyr::pull(.data$GeneName)) %in% (testSet %>% dplyr::pull(.data$GeneName)),
        "GeneList",
        "Universe"),
        levels=c("GeneList", "Universe"),
        ordered = TRUE)

      ##Kolmogorov-Smirnov test
      ### May not be the best choice because genomic distances are integers (not continuous) which generates ties
      ### But overal KS p-values appear well correlated to those of other tests
      suppressWarnings(
        distsum$KS.pvalue[i] <- ks.test(testSet %>% dplyr::pull(.data$Distance),
                                        univSet[isInTestSet=="Universe",] %>% dplyr::pull(.data$Distance),
                                        exact=F)$p.value
      )

      ##Mann-Whitney U test
      distsum$Wilcox.pvalue[i] <- wilcox.test(testSet %>% dplyr::pull(.data$Distance),
                                              univSet[isInTestSet=="Universe",] %>% dplyr::pull(.data$Distance),
                                              paired = FALSE)$p.value

      ##Test of independence (coin package, default asymptotic p-value)
      distsum$Independ.pvalue[i] <- coin::pvalue(
                                      coin::independence_test(
                                          univSet %>% dplyr::pull(.data$Distance) ~ isInTestSet)
                                      )

    }
}

return(list("distances" = distSel,
            "stats" = distsum))
}



#' @title Analyze the distances to the upstream and downstream neighbors, for different categories of orientations (SameStrand / OppositeStrand)
#'
#' @description Essentially a wrapper that does some prefiltering on gene list,
#'   applies \code{\link{statDistanceSide}} to both Upstream and Downstream distances and formats the results
#'
#' @param GeneList A character vector of "focus" genes to analyze
#' @param GeneNeighborhood A \code{\link{tibble}} obtained with the \code{\link{getGeneNeighborhood}} function or any data frame
#'  with the following columns:
#' \itemize{
#'   \item \code{GeneName}: Name of the focus gene
#'   \item \code{Upstream}: Name of the gene located upstream of the focus gene
#'   \item \code{Downstream}: Name of the gene located downstream of the focus gene
#'   \item \code{UpstreamClass}. Class of the Upstream gene (e.g. SameStrand / OppositeStrand). The class name must contain the string 'overlap' in case of overlap
#'   \item \code{DownstreamClass}. Class of the Downstream gene
#'   \item \code{UpstreamDistance}. Distance between the focus gene and the Upstream gene
#'   \item \code{DownstreamDistance}. Distance between the focus gene and the Downstream gene
#' }
#' @param DistriTest A logical (default is TRUE) indicating if the distribution of distances should be compared to the reference universe
#' @param GeneUniverse An optional character vector of genes in the universe. By default all genes in \code{GeneNeighborhood} are considered in the Universe.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select bind_rows mutate rename
#' @importFrom tibble as_tibble
#' @importFrom rlang .data
#'
#' @export
#'
#' @return A list with:
#' \itemize{
#'   \item distances. Upstream and downstream intergenic distances for non overlapping genes
#'   \item stats. Statistics on the upstream/downstream intergenic distances for non overlapping genes
#' }
#'
#' @seealso \code{\link{statDistanceSide}}
#'
#' @section DETAILS:
#' The function only returns the distances for the GeneList, not for the GeneUniverse.
#' To get the distances for the GeneUniverse, use the GeneUniverse as the GeneList and DistriTest=FALSE (see examples)
#'
#' Overlapping genes are removed.
#' The function should thus return results for the following orientations:
#' \itemize{
#'   \item \code{SameStrand}: The upstream/downstream gene is on the same strand as the focus gene
#'   \item \code{OppositeStrand}: The upstream/downstream gene is on the opposite strand of the focus gene
#' }
#'
#' @examples
#' ## Obtain gene neighborhood information:
#'   GeneNeighbors <- getGeneNeighborhood(Genegr)
#' ## Get a (random) set of (100) genes:
#'   set.seed(123)
#'   randGenes <- sample(names(Genegr), 100)
#' ## Extract all distances and compare them to reference genes:
#'   NDS <- analyzeNeighborsDistance(GeneList = randGenes,
#'                                   GeneNeighborhood = GeneNeighbors)
#'   NDS$stats
#' ## Extract the distances for all non-overlapping genes:
#'   alldist <- analyzeNeighborsDistance(GeneList = names(Genegr),
#'                                       GeneNeighborhood = GeneNeighbors,
#'                                       DistriTest = FALSE)
#' ## Some stats on these distances:
#'   alldist$stats
#'
#' @author Pascal GP Martin
#'

analyzeNeighborsDistance <- function(GeneList,
                                     GeneNeighborhood=NULL,
                                     DistriTest = TRUE,
                                     GeneUniverse = NULL) {


  ##----------
  ## Check arguments
  ##----------

  if (is.null(GeneNeighborhood)) {
    stop("GeneNeighborhood dataset should be provided")
  }

  if (!all(c("GeneName", "Upstream", "UpstreamClass", "UpstreamDistance",
             "Downstream", "DownstreamClass", "DownstreamDistance") %in% colnames(GeneNeighborhood))) {
      stop("The GeneNeighborhood object must contain columns 'Upstream', 'Downstream', 'UpstreamClass', 'DownstreamClass', 'UpstreamDistance' and 'DownstreamDistance'")
  }

  if (length(intersect(as.character(GeneList), as.character(GeneNeighborhood$GeneName))) == 0) {
      stop("No intersection between GeneList and the GeneName column in GeneNeighborhood")
  }


  ##----------
  ## Extract the relevant info from GeneNeighborhood
  ##----------
  GNN <- GeneNeighborhood[,c("GeneName", "Upstream", "Downstream",  "UpstreamClass",
                             "DownstreamClass", "UpstreamDistance", "DownstreamDistance")]
  if (any(apply(GNN, 2, is.factor))) {
    GNN <- tibble::as_tibble(apply(GNN, 2, as.character))
  }

  ##----------
  ## Prefilter GeneList
  ##----------
  cat("\nPrefiltering of GeneList and GeneUniverse:\n")
  cat("==========================================\n")

  #Remove NAs from GeneList
  isNAGene <- is.na(GeneList)
  if (any(isNAGene)) {
    cat(sum(isNAGene), "NA values removed from GeneList")
    GeneList <- GeneList[!isNAGene]
  }

  ## Remove genes that are not in the GeneNeighborhood table
  isAnalyzed <- GeneList %in% GNN$GeneName
  if (any(!isAnalyzed)) {
    cat(sum(!isAnalyzed), "genes from GeneList are not in the GeneNeighborhood table and are removed\n")
    GeneList <- GeneList[isAnalyzed]
  }

  ## Total number of genes in GeneList
  ng <- length(GeneList)
  cat("Total number of genes in GeneList:", ng, "\n")


  ##----------
  ## Get info for the universe
  ##----------

  if (DistriTest) {
    if (is.null(GeneUniverse)) {
      GeneUniverse <- GNN$GeneName
    } else {

      #Check for NA values in GeneUniverse
      isNAuniv <- is.na(GeneUniverse)
      if (any(isNAuniv)) {
        cat(sum(isNAuniv), " NA values removed from GeneUniverse")
        GeneUniverse <- GeneUniverse[!isNAuniv]
      }

      #Check for genes in GeneUniverse that are not in GeneNeighborhood table
      UnivHasNeighbor <- GeneUniverse %in% GNN$GeneName
      if (any(!UnivHasNeighbor)) {
        cat(sum(!UnivHasNeighbor), "genes from GeneUniverse are not in the GeneNeighborhood table and are removed\n")
        GeneUniverse <- GeneUniverse[UnivHasNeighbor]
      }

      #Merge GeneUniverse and GeneList
      GeneUniverse <- union(GeneUniverse, GeneList)
    }
    cat("Total number of genes in GeneUniverse (including GeneList) is", length(GeneUniverse), "\n")
  } else {GeneUniverse = NULL}


  ##----------
  ## For upstream genes
  ##----------
  upres <- statDistanceSide(GeneNeighborhood = GNN,
                            glist = GeneList,
                            distest = DistriTest,
                            guniv = GeneUniverse,
                            Side = "Upstream")


  ##----------
  ## For downstream genes
  ##----------
  dnres <- statDistanceSide(GeneNeighborhood = GNN,
                            glist = GeneList,
                            distest = DistriTest,
                            guniv = GeneUniverse,
                            Side = "Downstream")


  ##----------
  #Format results
  ##----------
  res <- list()
  #distances
  res$distances <- dplyr::bind_rows(upres$distances, dnres$distances) %>%
    dplyr::mutate(Side = factor(rep(c("Upstream", "Downstream"),
                                    times=c(nrow(upres$distances), nrow(dnres$distances))),
                                levels=c("Upstream", "Downstream"),
                                ordered=T)) %>%
    dplyr::rename("Neighbor" = "SideGene",
                  "Orientation" = "SideClass") %>%
    dplyr::select(.data$GeneName, .data$Neighbor, .data$Side, .data$Orientation, .data$Distance)

  #statistics
  res$stats <- dplyr::bind_rows(upres$stats, dnres$stats) %>%
    dplyr::mutate(Side = factor(rep(c("Upstream", "Downstream"),
                                    times=c(nrow(upres$stats), nrow(dnres$stats))),
                                levels=c("Upstream", "Downstream"),
                                ordered=T)) %>%
    dplyr::rename("Orientation" = "SideClass")

  if (DistriTest) {
    res$stats <- res$stats %>%
      dplyr::select(.data$GeneGroup, .data$Side, .data$Orientation, .data$n:.data$SEM, .data$KS.pvalue:.data$Independ.pvalue)
  } else {
    res$stats <- res$stats %>%
      dplyr::select(.data$Side, .data$Orientation, .data$n:.data$SEM)
  }
  #Return res
  return(res)
}

