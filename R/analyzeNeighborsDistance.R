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
#' @param confLevel Confidence level for the intervals on the mean and median
#' @param nboot Number of bootstrap replicates used to estimate confidence intervals of the mean and median (default to 1e4).
#' @param CItype Type of bootstrap confidence interval ("perc" for classical percentile or "bca" for bias-corrected and accelerated intervals). Passed to \code{\link[boot]{boot.ci}} function.
#' @param ncores Number of processes for parallel computation (passed to \code{\link[boot]{boot}} function).
#'
#' @importFrom magrittr %>% %<>%
#' @importFrom rlang .data !!
#' @importFrom dplyr select rename slice transmute pull group_by summarise bind_rows one_of
#' @importFrom tidyr nest unnest
#' @importFrom purrr map
#' @importFrom coin pvalue independence_test
#' @importFrom stats quantile median sd ks.test wilcox.test pnorm
#' @importFrom boot boot boot.ci
#' @importFrom parallel detectCores
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
#'   \item Kolmogorov-Smirnov test. See \code{\link[stats]{ks.test}}. Although not adapted to integer values, it should give concervative p-values for large enough gene sets.
#'   \item Mann-Whithney U test. See \code{\link[stats]{wilcox.test}}.
#'   \item Independence test. See \code{\link[coin]{independence_test}} in the \code{coin} package
#' }
#'
#' @seealso
#' \code{\link[stats]{ks.test}}, \code{\link[stats]{wilcox.test}},
#' \code{\link[coin]{independence_test}},
#' \code{\link[boot]{boot}}, \code{\link[boot]{boot.ci}}
#'
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
#'                    Side = "Upstream",
#'                    nboot = 1e3,
#'                    ncores = 3)
#'
#' @author Pascal GP Martin
#'

statDistanceSide <- function(GeneNeighborhood = NULL,
                             glist = NULL,
                             guniv = NULL,
                             distest = TRUE,
                             Side=c("Upstream", "Downstream"),
                             confLevel = 0.95,
                             nboot = 1e4,
                             CItype=c("bca", "perc"),
                             ncores = NULL) {

#-----------
# Check arguments
#-----------


if (is.null(GeneNeighborhood)) {
  stop("GeneNeighborhood dataset should be provided")
}

## Check that required columns are present:
className <- paste0(Side, "Class")
DistName <- paste0(Side, "Distance")

selectedColumns <- c("GeneName",
                     Side,
                     className,
                     DistName)

if (!all(selectedColumns %in% colnames(GeneNeighborhood))) {
    stop("The GeneNeighborhood object must contain columns: ",
         paste(selectedColumns, collapse=" "))
}

## Check gene Universe:
if (is.null(glist)) {
  guniv = NULL
  distest = FALSE
  glist <- as.character(GeneNeighborhood$GeneName)
} else {
  if (is.null(guniv) && distest == TRUE) {
    guniv <- as.character(GeneNeighborhood$GeneName)
  }
  if (length(intersect(glist, GeneNeighborhood$GeneName)) == 0) {
    stop("No intersection between the gene set and
         the GeneName column in GeneNeighborhood")
  }
}

## SetCItype
CItype <- match.arg(CItype)

## check confLevel
if (confLevel > 100 || confLevel <=0) {
    stop("conflevel should not be >100 or <=0")
}

if (confLevel >1) {
    confLevel = confLevel / 100
}

## Parameters for parallel computing of bootstrap samples
parallelType <- "multicore"
if (is.null(ncores)) {
    availCores <- parallel::detectCores() -1
    ncores <- ifelse(!is.na(availCores) && availCores > 0,
                     availCores,
                     1)
}

if (.Platform$OS.type=="windows") {
    ncores <- 1
}

if (ncores == 1) {parallelType <- "no"}

  #Select data based on Side
GNNdata <- GeneNeighborhood %>%
            dplyr::select(dplyr::one_of(selectedColumns)) %>%
            dplyr::rename("SideGene" = !!Side,
                          "SideClass" = !!className,
                          "Distance" = !!DistName)


  #Start the analysis
cat("\nAnalysis of Distances to", tolower(Side), "gene:\n")
cat("========================================\n")

  #Remove gene that have no upstream/downstream data (i.e. class=="other")
hasSide <- GNNdata %>%
    dplyr::slice(match(glist, .data$GeneName)) %>%
    dplyr::transmute(hasSide = !is.na(.data$SideGene) &
                         !is.na(.data$Distance)) %>%
    dplyr::pull(hasSide)

if (any(!hasSide)) {
    cat(sum(!hasSide), paste0("genes with undefined ",
                              tolower(Side),
                              " gene are removed\n"))
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
    cat(sum(isOVL), paste0("genes with an overlapping ",
                           tolower(Side),
                           " gene are removed from GeneList\n"))
    glistSel <- glistSel[!isOVL]
}


  ##Count the remaining genes
ngSel <- length(glistSel) #Number of genes with UPSTREAM/DOWNSTREAM info
cat("GeneList for",
    tolower(Side),
    "gene analysis has",
    ngSel,
    "genes\n")

  ##Extract distances
distSel <- GNNdata %>%
    dplyr::filter(.data$GeneName %in% glistSel)

  ## Get statistics on the distances
# Code for bootstrap with bca CI is adapted from: https://www.painblogr.org/2017-10-18-purrring-through-bootstraps.html

## split/nest the data by SideClass
distByClass <- distSel %>%
    dplyr::group_by(.data$SideClass) %>%
    tidyr::nest()


##Functions on which we want boostrap confidence intervals
mystatfun <- function(genodist, indices) {
    c(mean(genodist[indices], na.rm = TRUE),
      median(genodist[indices], na.rm = TRUE))
    }


##Add the boostrap samples to the data structure:
distByClass %<>%
        dplyr::mutate(bootSamples =
                          purrr::map(.x = .data$data,
                                     ~boot::boot(data = .x$Distance,
                                                 statistic = mystatfun,
                                                 R = nboot,
                                                 stype = "i",
                                                 parallel = parallelType,
                                                 ncpus = ncores)))

## Get confidence interval for the mean
# When nboot is too small, we can get the error "estimated adjustment 'a' is NA"
# Increasing nboot should fix the issue
distByClass %<>%
    dplyr::mutate(bootCI_mean = purrr::map(.x = .data$bootSamples,
                                           ~boot::boot.ci(.x,
                                                          conf = confLevel,
                                                          type = CItype,
                                                          index = 1)))

## Get confidence interval for the median
distByClass %<>%
    dplyr::mutate(bootCI_median = purrr::map(.x = .data$bootSamples,
                                           ~boot::boot.ci(.x,
                                                          conf = confLevel,
                                                          type = CItype,
                                                          index = 2)))

## Get all statistics in a table
citype <- ifelse(CItype=="perc", "percent", "bca")

distsum <- distByClass %>%
    dplyr::mutate(n = purrr::map(.x = .data$data,
                                 ~sum(!is.na(.x$Distance))),
                  Min = purrr::map(.x = .data$data,
                                   ~min(.x$Distance, na.rm = TRUE)),
                  Q1 = purrr::map(.x = .data$data,
                                  ~quantile(.x$Distance, 0.25, na.rm = TRUE)),
                  Median = purrr::map(.x = .data$bootCI_median,
                                      ~.x$t0),
                  Median_lowerCI = purrr::map(.x = .data$bootCI_median,
                                              ~.x[[citype]][[4]]),
                  Median_upperCI = purrr::map(.x = .data$bootCI_median,
                                              ~.x[[citype]][[5]]),
                  Mean = purrr::map(.x = .data$bootCI_mean,
                                    ~.x$t0),
                  Mean_lowerCI = purrr::map(.x = .data$bootCI_mean,
                                              ~.x[[citype]][[4]]),
                  Mean_upperCI = purrr::map(.x = .data$bootCI_mean,
                                              ~.x[[citype]][[5]]),
                  Q3 = purrr::map(.x = .data$data,
                                  ~quantile(.x$Distance, 0.75, na.rm = TRUE)),
                  Max = purrr::map(.x = .data$data,
                                   ~max(.x$Distance, na.rm = TRUE)),
                  SD = purrr::map(.x = .data$data,
                                  ~sd(.x$Distance, na.rm = TRUE)),
                  SEM = purrr::map(.x = .data$data,
                                   ~sd(.x$Distance, na.rm = TRUE)/
                                       sqrt(sum(!is.na(.x$Distance))))
    ) %>%
    dplyr::select(-.data$data, -.data$bootSamples,
                  -.data$bootCI_median, -.data$bootCI_mean) %>%
    tidyr::unnest(cols = c(.data$n, .data$Min, .data$Q1, .data$Median,
                           .data$Median_lowerCI, .data$Median_upperCI,
                           .data$Mean, .data$Mean_lowerCI, .data$Mean_upperCI,
                           .data$Q3, .data$Max, .data$SD, .data$SEM))

  #Compare to universe
  if (distest) {

    #Remove from universe genes without a defined upstream/downstream gene
    univhasSide <- GNNdata %>%
      dplyr::slice(match(guniv, .data$GeneName)) %>%
      dplyr::transmute(hasSide = !is.na(.data$SideGene) &
                           !is.na(.data$Distance)) %>%
      dplyr::pull(hasSide)

    if (any(!univhasSide)) {
      cat(sum(!univhasSide), paste0("genes from GeneUniverse with undefined ",
                                    tolower(Side),
                                    " gene are removed\n"))
    }

    gunivSel <- guniv[univhasSide]

    #Remove genes with overlapping upstream/downstream gene
    univIsOVL <- GNNdata %>%
      dplyr::slice(match(gunivSel, .data$GeneName)) %>%
      dplyr::transmute(isOVL = grepl("overlap",
                                     tolower(.data$SideClass)) &
                           .data$Distance == 0) %>%
      dplyr::pull(isOVL)

    if (any(univIsOVL)) {
      cat(sum(univIsOVL),
          paste0("genes from GeneUniverse with an overlapping ",
                 tolower(Side),
                 " gene are removed\n"))
      gunivSel <- gunivSel[!univIsOVL]
    }


    #Count the remaining genes
    ngunivSel <- length(gunivSel)
    cat("Universe for",
        tolower(Side),
        "gene analysis has",
        ngunivSel,
        "genes\n")

    #Get distances for Universe
    distSelUniv <- GNNdata %>%
                       dplyr::filter(.data$GeneName %in% gunivSel)


    #Summarize the distances for the Universe
    distsumUniv <- distSelUniv %>%
      dplyr::group_by(.data$SideClass) %>%
      dplyr::summarise(n = sum(!is.na(.data$Distance)),
                       Min = min(.data$Distance, na.rm=TRUE),
                       Q1 = quantile(.data$Distance, 0.25, na.rm=TRUE),
                       Median = median(.data$Distance, na.rm=TRUE),
                       Median_lowerCI = NA,
                       Median_upperCI = NA,
                       Mean = mean(.data$Distance, na.rm=TRUE),
                       Mean_lowerCI = NA,
                       Mean_upperCI = NA,
                       Q3 = quantile(.data$Distance, 0.75, na.rm=TRUE),
                       Max = max(.data$Distance, na.rm=TRUE),
                       SD = sd(.data$Distance, na.rm=TRUE),
                       SEM = sd(.data$Distance, na.rm=TRUE)/
                           sqrt(sum(!is.na(.data$Distance)))
      )

    #Merge the summary datasets
    nr <- nrow(distsum)
    distsum <- dplyr::bind_rows(distsum, distsumUniv)
    distsum$GeneGroup <- factor(rep(c("GeneList", "GeneUniverse"),
                                    times=c(nr, nrow(distsumUniv))),
                                levels=c("GeneList", "GeneUniverse"),
                                ordered = TRUE)

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

      isInTestSet <- (univSet %>% dplyr::pull(.data$GeneName)) %in%
                         (testSet %>% dplyr::pull(.data$GeneName))

      ##Kolmogorov-Smirnov test
      ### May not be the best choice because genomic distances are integers (not continuous) which generates ties
      ### But the KS test is known to produce conservative p-values for discrete distributions when sample size is large enough
      ### Here we use the continuous version but it should be better to use the version for discrete variables from the dgof package:
      ### dgof::ks.test(testSet, ecdf(universe))
      suppressWarnings(
        distsum$KS.pvalue[i] <- ks.test(testSet %>%
                                            dplyr::pull(.data$Distance),
                                        univSet[!isInTestSet,] %>%
                                            dplyr::pull(.data$Distance),
                                        exact = FALSE)$p.value
      )

      ##Mann-Whitney U test
      distsum$Wilcox.pvalue[i] <- wilcox.test(testSet %>%
                                                  dplyr::pull(.data$Distance),
                                              univSet[!isInTestSet,] %>%
                                                  dplyr::pull(.data$Distance),
                                              paired = FALSE)$p.value

      ##Test of independence (coin package, default asymptotic p-value)
      distsum$Independ.pvalue[i] <- coin::pvalue(
                                      coin::independence_test(
                                          univSet %>%
                                              dplyr::pull(.data$Distance) ~
                                              isInTestSet)
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
#' @param confLevel Confidence level for the intervals on the mean and median
#' @param nboot Number of bootstrap replicates used to estimate confidence intervals of the mean and median
#' @param CItype Type of bootstrap confidence interval ("perc" for classical percentile or "bca" for bias-corrected and accelerated intervals). Passed to \code{\link[boot]{boot.ci}} function
#' @param ncores Number of processes for parallel computation (passed to \code{\link[boot]{boot}} function).
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
#'                                   GeneNeighborhood = GeneNeighbors,
#'                                   nboot = 1e3,
#'                                   ncores = 2)
#'   NDS$stats
#' ## Extract the distances for all non-overlapping genes:
#'   alldist <- analyzeNeighborsDistance(GeneList = names(Genegr),
#'                                       GeneNeighborhood = GeneNeighbors,
#'                                       DistriTest = FALSE,
#'                                       nboot = 1e3,
#'                                       ncores = 2)
#' ## Some stats on these distances:
#'   alldist$stats
#'
#' @author Pascal GP Martin
#'

analyzeNeighborsDistance <- function(GeneList,
                                     GeneNeighborhood=NULL,
                                     DistriTest = TRUE,
                                     GeneUniverse = NULL,
                                     confLevel = 0.95,
                                     nboot = 1e4,
                                     CItype = "bca",
                                     ncores = NULL) {


  ##----------
  ## Check arguments
  ##----------

    if (is.null(GeneNeighborhood)) {
      stop("GeneNeighborhood dataset should be provided")
    }

    if (!all(c("GeneName",
               "Upstream", "UpstreamClass", "UpstreamDistance",
               "Downstream", "DownstreamClass", "DownstreamDistance") %in%
             colnames(GeneNeighborhood))) {
      stop("The GeneNeighborhood object must contain columns
           'Upstream', 'Downstream',
           'UpstreamClass', 'DownstreamClass',
           'UpstreamDistance' and 'DownstreamDistance'")
  }

  if (length(intersect(as.character(GeneList),
                       as.character(GeneNeighborhood$GeneName))) == 0) {
      stop("No intersection between GeneList and
           the GeneName column in GeneNeighborhood")
  }

    if (confLevel > 100 || confLevel <=0) {
        stop("conflevel should not be >100 or <=0")
    }

    if (confLevel >1) {
        confLevel = confLevel / 100
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

cat("\n",
    paste0(round(100*confLevel, 1),
           "% bootstrap CI of the mean and median are obtained from ",
           nboot,
           " replicates"),
    "\n")

  ##----------
  ## For upstream genes
  ##----------
  upres <- statDistanceSide(GeneNeighborhood = GNN,
                            glist = GeneList,
                            distest = DistriTest,
                            guniv = GeneUniverse,
                            Side = "Upstream",
                            confLevel = confLevel,
                            nboot = nboot,
                            CItype = CItype,
                            ncores = ncores)


  ##----------
  ## For downstream genes
  ##----------
  dnres <- statDistanceSide(GeneNeighborhood = GNN,
                            glist = GeneList,
                            distest = DistriTest,
                            guniv = GeneUniverse,
                            Side = "Downstream",
                            confLevel = confLevel,
                            nboot = nboot,
                            CItype = CItype,
                            ncores = ncores)


  ##----------
  #Format results
  ##----------
  res <- list()
  #distances
  res$distances <- dplyr::bind_rows(upres$distances,
                                    dnres$distances) %>%
    dplyr::mutate(Side = factor(rep(c("Upstream", "Downstream"),
                                    times=c(nrow(upres$distances),
                                            nrow(dnres$distances))),
                                levels=c("Upstream", "Downstream"),
                                ordered = TRUE)) %>%
    dplyr::rename("Neighbor" = "SideGene",
                  "Orientation" = "SideClass") %>%
    dplyr::select(.data$GeneName,
                  .data$Neighbor,
                  .data$Side,
                  .data$Orientation,
                  .data$Distance)

  #statistics
  res$stats <- dplyr::bind_rows("Upstream" = upres$stats,
                                "Downstream" = dnres$stats,
                                .id="Side") %>%
    dplyr::mutate(Side = factor(.data$Side,
                                levels = c("Upstream", "Downstream"),
                                ordered = TRUE)) %>%
    dplyr::rename("Orientation" = "SideClass")

  if (DistriTest) {
    res$stats <- res$stats %>%
      dplyr::select(.data$GeneGroup,
                    .data$Side,
                    .data$Orientation,
                    .data$n:.data$SEM,
                    .data$KS.pvalue:.data$Independ.pvalue)
  } else {
    res$stats <- res$stats %>%
      dplyr::select(.data$Side,
                    .data$Orientation,
                    .data$n:.data$SEM)
  }
  #Return res
  return(res)
}

