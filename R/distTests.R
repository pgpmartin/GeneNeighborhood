#' Compare 2 sets of distances with a Kolmogorov-Smirnov test
#' (only returns the p-value)
#'
#' @param .x A data frame with distances for the gene set of interest.
#' @param .y A data frame with distances for the control/reference gene set.
#'
#' @section DETAILS:
#' Both \code{.x} and \code{.y} should contain at least the following columns:
#' \itemize{
#'     \item Distance. Distance in bp.
#'     \item GeneName. Name of the focus gene.
#' }
#' Note that the Kolmogorov-Smirnov test is not adapted to integer values such
#' as intergenic distances in bp. However, it gives conservative p-values for
#' large enough gene sets.
#'
#' @return The p-value of the Kolmorovov-Smirnov test.
#'
#' @seealso \code{\link[stats]{ks.test}}
#'
#' @examples
#' ## Create some random distance data
#'   set.seed(123)
#'   mydistData <- data.frame(GeneName = stringi::stri_rand_strings(400, 5),
#'                            Distance = sample(1:5000, 400, replace=TRUE),
#'                            GeneSet = sample(c("TestSet", "RefSet"),
#'                                             400, replace= TRUE))
#' ## Compute the p-value of the KS test comparing TestSet to RefSet
#'   ksfun(.x = mydistData[mydistData$GeneSet == "TestSet",],
#'         .y = mydistData[mydistData$GeneSet == "RefSet",])
#' ## If .x is empty, the function returns NA
#'   ksfun(.x = mydistData[mydistData$GeneSet == "SomeRandomName",],
#'         .y = mydistData[mydistData$GeneSet == "RefSet",])
#' ## If .y is empty, the function returns an error
#' \dontrun{
#'   ksfun(.x = mydistData[mydistData$GeneSet == "TestSet",],
#'         .y = mydistData[mydistData$GeneSet == "SomeRandomName",])
#'         }

ksfun <- function(.x, .y) {
    if (nrow(.x)>0) {
        suppressWarnings(
            ks.test(.x %>% dplyr::pull(.data$Distance),
                    .y %>%
                        dplyr::filter(!(.data$GeneName %in% .x$GeneName)) %>%
                        dplyr::pull(.data$Distance),
                    exact = FALSE)$p.value
        )
    } else {NA}
}


#' Compare 2 sets of distances with a Mann-Whitney U test
#' (only returns the p-value)
#'
#' @param .x A data frame with distances for the gene set of interest.
#' @param .y A data frame with distances for the control/reference gene set.
#'
#' @section DETAILS:
#' Both \code{.x} and \code{.y} should contain at least the following columns:
#' \itemize{
#'     \item Distance. Distance in bp.
#'     \item GeneName. Name of the focus gene.
#' }
#' The Mann-Whithney U test is also called the Wilcoxon rank sum test.
#'
#' @return The p-value of the Mann-Whitney U test.
#'
#' @seealso \code{\link[stats]{wilcox.test}}
#'
#' @examples
#' ## Create some random distance data
#'   set.seed(123)
#'   mydistData <- data.frame(GeneName = stringi::stri_rand_strings(400, 5),
#'                            Distance = sample(1:5000, 400, replace=TRUE),
#'                            GeneSet = sample(c("TestSet", "RefSet"),
#'                                             400, replace= TRUE))
#' ## Compute the p-value of the Mann-Whitney U-test comparing TestSet to RefSet
#'   Utest(.x = mydistData[mydistData$GeneSet == "TestSet",],
#'         .y = mydistData[mydistData$GeneSet == "RefSet",])
#' ## If .x is empty, the function returns NA
#'   Utest(.x = mydistData[mydistData$GeneSet == "SomeRandomName",],
#'         .y = mydistData[mydistData$GeneSet == "RefSet",])
#' ## If .y is empty, the function returns an error
#' \dontrun{
#'   Utest(.x = mydistData[mydistData$GeneSet == "TestSet",],
#'         .y = mydistData[mydistData$GeneSet == "SomeRandomName",])
#'         }

Utest <- function(.x, .y) {
    if (nrow(.x)>0) {
        wilcox.test(.x %>% dplyr::pull(.data$Distance),
                    .y %>%
                        dplyr::filter(!(.data$GeneName %in% .x$GeneName)) %>%
                        dplyr::pull(.data$Distance),
                    paired = FALSE)$p.value
    } else {NA}
}

## TODO: compare to coin::wilcox_test which handles ties better
## See: https://stackoverflow.com/questions/23450221/coinwilcox-test-versus-wilcox-test-in-r
## See: https://stats.stackexchange.com/questions/31417/what-is-the-difference-between-wilcox-test-and-coinwilcox-test-in-r
## See: http://r.789695.n4.nabble.com/wilcox-test-function-in-coin-package-td4668314.html


#' Test the independence of 2 sets of distances with the coin package
#' (only returns the p-value)
#'
#' @param .x A data frame with distances for the gene set of interest.
#' @param .y A data frame with distances for the control/reference gene set.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select mutate
#' @importFrom rlang .data
#' @importFrom coin pvalue independence_test

#' @section DETAILS:
#' Both \code{.x} and \code{.y} should contain at least the following columns:
#' \itemize{
#'     \item Distance. Distance in bp.
#'     \item GeneName. Name of the focus gene.
#'     \item GeneSet. Name of the Gene Set
#' }
#'
#' @return The p-value of the independence test from the coin package.
#'
#' @seealso \code{\link[coin]{independence_test}}
#'
#' @examples
#' ## Create some random distance data
#'   set.seed(123)
#'   mydistData <- data.frame(GeneName = stringi::stri_rand_strings(400, 5),
#'                            Distance = sample(1:5000, 400, replace=TRUE),
#'                            GeneSet = sample(c("TestSet", "RefSet"),
#'                                             400, replace= TRUE))
#' ## Compute the p-value of the independence test
#'   coinIndep(.x = mydistData[mydistData$GeneSet == "TestSet",],
#'             .y = mydistData[mydistData$GeneSet == "RefSet",])
#' ## If .x is empty, the function returns NA
#'   coinIndep(.x = mydistData[mydistData$GeneSet == "SomeRandomName",],
#'             .y = mydistData[mydistData$GeneSet == "RefSet",])
#' ## If .y is empty, the function returns an error
#' \dontrun{
#'   coinIndep(.x = mydistData[mydistData$GeneSet == "TestSet",],
#'             .y = mydistData[mydistData$GeneSet == "SomeRandomName",])
#'         }

coinIndep <- function(.x, .y) {
    stopifnot(nrow(.y)>0)
    if (nrow(.x)>0) {
        df <- dplyr::union(.x %>% dplyr::select(-.data$GeneSet),
                           .y %>% dplyr::select(-.data$GeneSet)) %>%
            dplyr::mutate(isInTestSet = .data$GeneName %in% .x$GeneName)
        coin::pvalue(coin::independence_test(df$Distance ~ df$isInTestSet))
    } else {NA}
}


#' Resampling test on the median (only returns the p-value)
#'
#' @param .x A data frame with distances for the gene set of interest.
#' @param .y A data frame with distances for the control/reference gene set.
#' @param R Number of samples to draw from \code{.y}
#'
#' @section DETAILS:
#' Both \code{.x} and \code{.y} should contain at least the following columns:
#' \itemize{
#'     \item Distance. Distance in bp.
#'     \item GeneName. Name of the focus gene.
#'     \item GeneSet. Name of the Gene Set
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select pull union
#' @importFrom rlang .data
#' @importFrom matrixStats colMedians
#'
#' @return A p-value representing the percentage of times that the median of
#'         \code{R} random samples drawn from the union of \code{.x} and \code{.y}
#'         is below the observed median of \code{.x}.
#'
#' @seealso
#' \code{\link{replicate}}
#' \code{\link[matrixStats]{colMedians}}
#'
#' @examples
#' ## Create some random distance data
#'   set.seed(123)
#'   mydistData <- data.frame(GeneName = stringi::stri_rand_strings(400, 5),
#'                            Distance = sample(1:5000, 400, replace=TRUE),
#'                            GeneSet = sample(c("TestSet", "RefSet"),
#'                                             400, replace= TRUE))
#' ## Evaluate (by resampling) the probability to get a median lower that that of TestSet
#' ## We use only 1000 random samples here for speed purposes but you should use >1e4
#'   resampMed(.x = mydistData[mydistData$GeneSet == "TestSet",],
#'             .y = mydistData[mydistData$GeneSet == "RefSet",],
#'             R = 1e3)
#' ## If .x is empty, the function returns NA
#'   resampMed(.x = mydistData[mydistData$GeneSet == "SomeRandomName",],
#'             .y = mydistData[mydistData$GeneSet == "RefSet",],
#'             R = 1e3)
#' ## If .y is empty, the function returns an error
#' \dontrun{
#'   resampMed(.x = mydistData[mydistData$GeneSet == "TestSet",],
#'             .y = mydistData[mydistData$GeneSet == "SomeRandomName",],
#'             R = 1e3)
#'         }

resampMed <- function(.x, .y, R = 1e4) {
    stopifnot(nrow(.y)>0)
    if (nrow(.x)>0) {
        refmed <- median(.x %>% dplyr::pull(.data$Distance))
        univdist <- dplyr::union(.y %>% dplyr::select(-.data$GeneSet),
                                 .x %>% dplyr::select(-.data$GeneSet)) %>%
                        dplyr::pull(.data$Distance)
        ngenes <- nrow(.x)
        mean(matrixStats::colMedians(
             replicate(R,
                       sample(univdist, ngenes))) <= refmed)
    } else {NA}
}

#' Apply a test to different genesets using the same refereence/control set
#'
#' @param tbdist A data frame with the following columns:
#'     \itemize{
#'         \item Distance. Distance in bp.
#'         \item GeneName. Name of the focus gene.
#'         \item GeneSet. Name of the Gene Set
#'     }
#' @param univ Character string giving the name of the control set
#'             (should be in the GeneSet column).
#' @param FUN A function returning a single pvalue and that takes at least 2 arguments:
#'            the data frame .x and .y, of the same form as \code{tbdist},
#'            containing the distances for the test set and for the control/reference set.
#'            and returns a single p-value.
#' @param ... further arguments passed to FUN
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter pull union
#' @importFrom rlang .data
#' @importFrom tibble tibble
#'
#' @return A tibble with 2 columns: GeneSet and pvalue
#'
#' @seealso
#' \code{\link{replicate}}
#' \code{\link[matrixStats]{colMedians}}
#'
#' @examples
#' ## Create some random distance data
#'   set.seed(123)
#'   mydistData <- data.frame(GeneName = stringi::stri_rand_strings(600, 5),
#'                            Distance = sample(1:5000, 600, replace=TRUE),
#'                            GeneSet = sample(c("TestSet1", "TestSet2", "RefSet"),
#'                                             600, replace= TRUE))
#' ## Apply Kolmogorov-Smirnov test to the different test sets
#'   pvalByGeneSet(tbdist = mydistData,
#'                 univ = "RefSet",
#'                 FUN = ksfun)

pvalByGeneSet <- function(tbdist, univ = Universe, FUN, ...) {
    ud <- tbdist %>%
        dplyr::filter(.data$GeneSet == univ)
    td <- tbdist %>%
        dplyr::filter(.data$GeneSet != univ)
    pval <- unlist(
        by(td,
           td$GeneSet,
           FUN,
           .y = ud,
           ...,
           simplify = FALSE))
    pval <- tibble::tibble(GeneSet = names(pval),
                           pvalue = pval)
    return(pval)
}


#' @title Statistical tests for intergenic distance data
#'
#' @description Statistical tests for intergenic distance data
#'
#' @param GeneSetDistances A \code{tibble} with intergenic distances for the
#'                         different gene sets, as generated by the
#'                         \code{\link{dist2Neighbors}} function.
#' @param Universe Character string indicating which set should be considered
#'                  as the universe (or control set)
#' @param MedianResample Logical. Should the resample test of the median
#'                       be performed (defaults to TRUE)
#' @param R integer giving the number of resampling to perform for the
#'           resampling test (Default to 1e4)
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select rename slice transmute mutate pull
#'             group_by summarise bind_rows one_of
#' @importFrom reshape2 melt
#' @importFrom rlang .data !!
#' @importFrom tibble as_tibble
#' @importFrom matrixStats colMedians
#'
#' @export
#'
#' @return A tibble with the following columns:
#' \itemize{
#'   \item GeneSet. Name of the gene set.
#'   \item Orientation. Orientation of the neighbor (same or opposite strand).
#'   \item Side. Upstream or Downstream.
#'   \item KS.pvalue. p-value of Kolmogorov-Smirnof test
#'   \item Wilcox.pvalue. p-value of Wilcoxon rank sum test (or Mann-Whitney U test).
#'   \item Independ.pvalue. p-value of the independance test
#' }
#' if \code{MedianResample} is \code{TRUE} the tibble will also contain this additional column:
#' \itemize{
#'   \item Resample.pvalue. p-value from the resampling test.
#' }
#'
#' @section DETAILS:
#' The following tests are possible:
#' \itemize{
#'   \item Kolmogorov-Smirnov test. See \code{\link[stats]{ks.test}}.
#'         Although not adapted to integer values, it gives conservative
#'         p-values for large enough gene sets.
#'   \item Wilcoxon rank sum test (or Mann-Whithney U test).
#'         See \code{\link[stats]{wilcox.test}}.
#'   \item Independence test. See \code{\link[coin]{independence_test}}
#'         in the \code{coin} package
#'   \item resample. A test based on random resampling of the universe distances.
#' }
#'
#' @seealso \code{\link[stats]{ks.test}},
#'          \code{\link[stats]{wilcox.test}},
#'          \code{\link[coin]{independence_test}}
#'
#' @examples
#' #' ## Obtain gene neighborhood information:
#'   GeneNeighbors <- getGeneNeighborhood(Genegr)
#' ## Get a (random) set of (100) genes:
#'   set.seed(123)
#'   randGenes <- sample(names(Genegr), 100)
#' ## Create a set enriched for close upstream genes:
#'   GenePool <- GeneNeighbors[!is.na(GeneNeighbors$UpstreamDistance),]
#'   Proba <- (max(GenePool$UpstreamDistance)-GenePool$UpstreamDistance) /
#'              sum(max(GenePool$UpstreamDistance)-GenePool$UpstreamDistance)
#'   Proba <- (1/(GenePool$UpstreamDistance+1)) / sum(1/(GenePool$UpstreamDistance+1))
#'   CloseUpstream <- sample(GenePool$GeneName, size = 100, prob = Proba)
#' ## Extract distances for this set of genes and for all genes :
#'   myGeneSets <- list("RandomGenes" = randGenes,
#'                      "CloseUpstream" = CloseUpstream,
#'                      "AllGenes" = GeneNeighbors$GeneName)
#'   distForGeneSets <- dist2Neighbors(GeneNeighbors,
#'                                     myGeneSets)
#' ## Compare distances for genesets to a control set (here "AllGene")
#' ## (using only 1K permutations fo speed purposes here, prefer using > 1e4)
#'   distTests(distForGeneSets,
#'             Universe="AllGenes",
#'             MedianResample = TRUE,
#'             R=1e3)


distTests <- function(GeneSetDistances,
                      Universe = NULL,
                      MedianResample = TRUE,
                      R = 1e4) {

#-----------
# Check arguments
#-----------

## GeneSetDistances
if (!is.data.frame(GeneSetDistances)) {
    stop("GeneSetDistances should be a data frame")
}

if (!all(c("GeneName", "Neighbor", "Orientation",
           "Side", "Distance", "GeneSet") %in% colnames(GeneSetDistances))) {
    stop("GeneSet should have colums 'GeneName', 'Neighbor', 'Orientation', ",
         "'Side', 'Distance', and 'GeneSet'. See dist2Neighbors function.")
}


if (!is.character(Universe) || length(Universe)!=1) {
    stop("Universe should be a single character string")
}

if (!(Universe %in%
        unique(GeneSetDistances %>% dplyr::pull(.data$GeneSet)))) {
    stop("Universe was not found in GeneSet column")
}


gs <- GeneSetDistances %>%
    dplyr::select(-.data$Neighbor) %>%
    dplyr::group_by(.data$Side, .data$Orientation) %>%
    tidyr::nest(.key = "DistSet")


#-----------
# Add p-values to the table
#-----------

gs %<>%
    dplyr::mutate(KS = purrr::map(.x = .data$DistSet,
                                  ~ pvalByGeneSet(.x,
                                                 univ = Universe,
                                                 FUN = ksfun)),
                  Wilcox = purrr::map(.x = .data$DistSet,
                                      ~ pvalByGeneSet(.x,
                                                     univ = Universe,
                                                     FUN = Utest)),
                  Indep = purrr::map(.x = .data$DistSet,
                                     ~ pvalByGeneSet(.x,
                                                    univ = Universe,
                                                    FUN = coinIndep)))

if (MedianResample) {
    gs %<>%
        dplyr::mutate(Median_resample =
                          purrr::map(.x = .data$DistSet,
                                     ~ pvalByGeneSet(.x,
                                                     univ = Universe,
                                                     FUN = resampMed,
                                                     R = R)))
}

# gs <- gs %>% dplyr::select(-.data$DistSet)
gs <- gs %>%
      dplyr::select(-.data$DistSet) %>%
        tidyr::unnest() %>%
        dplyr::select(-.data$GeneSet1, -.data$GeneSet2) %>%
        dplyr::arrange(.data$GeneSet, .data$Side, .data$Orientation) %>%
        dplyr::select(.data$GeneSet, dplyr::everything()) %>%
        dplyr::rename("KS.pvalue" = "pvalue",
                      "Wilcox.pvalue" = "pvalue1",
                      "Indep.pvalue" = "pvalue2")

if (MedianResample) {
    gs <- gs %>%
            dplyr::select(-.data$GeneSet3) %>%
            dplyr::rename("Median_resample" = "pvalue3")
}

return(gs)

## See https://robertamezquita.github.io/post/2017-05-23-using-map-with-generic-functions-like-t-test/
## https://stackoverflow.com/questions/35558766/purrr-map-a-t-test-onto-a-split-df
## https://community.rstudio.com/t/applying-dunn-test-using-purrr-map/15155

# sets<- GeneSetDistances %>%
#          dplyr::filter(.data$GeneSet!=Universe)
#
# univ <- GeneSetDistances %>%
#     dplyr::filter(.data$GeneSet==Universe)

# GeneSetDistances %>%
#     dplyr::group_by(.data$Side, .data$Orientation) %>%
#     dplyr::mutate(KS.pvalue = ks.test( .data %>%
#                                            dplyr::filter(GeneSet!=Universe) %>%
#                                            dplyr::pull(Distance),
#                                        .data %>%
#                                            dplyr::filter(GeneSet==Universe) %>%
#                                            dplyr::pull(Distance),
#                                        exact = FALSE)$p.value)
##OK problem here is that we mix all gene sets together.
## We could probably apply this with map
## See https://jennybc.github.io/purrr-tutorial/index.html

# replicate(R, sample.int(nrow(x)))
# https://stats.stackexchange.com/questions/247004/statistical-significance-of-difference-between-distances
# https://stats.stackexchange.com/questions/24300/how-to-resample-in-r-without-repeating-permutations
# https://websites.pmc.ucsc.edu/~mclapham/Rtips/resampling.htm
# http://danielnee.com/2015/01/random-permutation-tests/

}
