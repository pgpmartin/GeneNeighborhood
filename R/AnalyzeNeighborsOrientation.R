#' @title Analyze the orientation of gene neighbors for a set of pre-defined genes
#'
#' @description Analyze the orientation of gene neighbors (upstream/downstream genes) for a set of pre-defined focus genes.
#'  Gives the number and percentage of genes in each orientation.
#'  Evaluate the enrichment of the different orientations when using a reference gene universe.
#'
#' @param GeneList A character vector of focus genes to analyze
#' @param GeneNeighborhood A tibble obtained with the GetGeneNeighborhood function. Can also be any data.frame with the relevant information.
#' @param colNames A named character vector of length 5 containing the names of the columns from \code{GeneNeighborhood} containing the following information (in this order):
#' \itemize{
#'   \item \code{GeneName}. Name of the focus gene
#'   \item \code{Upstream}. Name of the gene located upstream of the focus gene
#'   \item \code{Downstream}. Name of the gene located downstream of the focus gene
#'   \item \code{UpstreamClass}. Class of the Upstream gene
#'   \item \code{DownstreamClass}. Class of the Downstream gene
#' }
#' @param keepOther A logical (default to \code{TRUE}) indicating if the category "other" should be analyzed or not
#' @param EnrichTest A logical (default to \code{TRUE}) indicating if Enrichment test should be performed
#' @param GeneUniverse An optional character vector of genes in the universe. By default all genes in GeneNeighborhood are considered in the Universe.
#'
#' @importFrom dplyr slice mutate select pull filter count right_join bind_rows
#' @importFrom magrittr %>%
#' @importFrom stats fisher.test
#'
#' @export
#'
#' @return A \code{tibble} with orientation and distance for upstream/downstream gene.
#'  If \code{EnrichTest} is \code{TRUE}, the \code{tibble} also contains data for the universe and a p-value from a Fisher exact test evaluating the enrichment for each orientation in \code{GeneList}.
#'
#' @section DETAILS:
#' The function analyzes the following orientations:
#' \itemize{
#'   \item \code{SameStrand}. The upstream/downstream gene is on the same strand as the focus gene
#'   \item \code{OppositeStrand}. The upstream/downstream gene is on the opposite strand of the focus gene
#'   \item \code{OppositeOverlap}. The upstream/downstream gene is on the opposite strand of the focus gene and overlaps with the focus gene
#'   \item \code{SameOverlap}. The upstream/downstream gene is on the same strand of the focus gene and overlaps with the focus gene
#'   \item \code{other}. (only if \code{keepOther} is \code{TRUE}) Any other, more complex situation (e.g. multiple overlaps, focus gene contained within another gene, etc.)
#' }
#'
#' @seealso \code{\link{fisher.test}}
#'
#' @author Pascal GP Martin
#'
#' @examples
#' ## Obtain gene neighborhood information:
#'   GeneNeighbors <- GetGeneNeighborhood(Genegr)
#' ## get a random set of 100 genes:
#'   set.seed(123)
#'   randGenes <- names(Genegr)[sample.int(676, 100)]
#' ## Analyze the orientation of their neighbors
#'   AnalyzeNeighborsOrientation(randGenes,
#'                               GeneNeighborhood = GeneNeighbors)
#' ## Select a less random set of genes:
#'   isOpposite <- grepl("Opposite",GeneNeighbors$UpstreamClass)
#'   probs <- ifelse(isOpposite, 0.6/sum(isOpposite), 0.4/sum(!isOpposite))
#'   set.seed(123)
#'   lessrandGenes <- GeneNeighbors$GeneName[sample.int(676, 100, prob=probs)]
#'   AnalyzeNeighborsOrientation(lessrandGenes,
#'                               GeneNeighborhood = GeneNeighbors)
#'

AnalyzeNeighborsOrientation <- function(GeneList,
                                        GeneNeighborhood = GeneNeighbors,
                                        colNames = c("GeneName"="GeneName",
                                                     "Upstream"="Upstream", "Downstream"="Downstream",
                                                     "UpstreamClass"="UpstreamClass", "DownstreamClass"="DownstreamClass"),
                                        keepOther = TRUE,
                                        EnrichTest = TRUE,
                                        GeneUniverse = NULL) {

##----------
## verify colNames
##----------
if (length(colNames)!=5) {
  stop("colNames should be of length 5")
  }
if (is.null(names(colNames))) {
  warning("colNames is not a named vector, adding the following names:\n  c('GeneName', 'Upstream', 'Downstream', 'UpstreamClass', 'DownstreamClass')",
          "\nPlease make sure they are provided in this order")
  names(colNames) <- c("GeneName", "Upstream", "Downstream", "UpstreamClass", "DownstreamClass")
  }
if (!(all(c("GeneName", "Upstream", "Downstream", "UpstreamClass", "DownstreamClass") %in% names(colNames)))) {
  stop("names of colNames should be c('GeneName', 'Upstream', 'Downstream', 'UpstreamClass', 'DownstreamClass')")
  }

##----------
## Extract the relevant info from GeneNeighborhood
##----------
GNN <- GeneNeighborhood[,colNames[c("GeneName", "Upstream", "Downstream",  "UpstreamClass", "DownstreamClass")]]
names(GNN) <- c("GeneName", "Upstream", "Downstream", "UpstreamClass", "DownstreamClass")
if (!is.character(GNN$GeneName)) {
    GNN$GeneName <- as.character(GNN$GeneName)
}

##----------
## Prefilter GeneList
##----------

#Remove NAs from GeneList
isNAGene <- is.na(GeneList)
if (any(isNAGene)) {
  warning(sum(isNAGene), " genes are NA. They are removed from GeneList")
  GeneList <- GeneList[!isNAGene]
}

## Remove genes that are not in the GeneNeighborhood table
isAnalyzed <- GeneList %in% GNN$GeneName
if (any(!isAnalyzed)) {
  warning(sum(!isAnalyzed), " genes from GeneList are not in the GeneNeighborhood table and are removed")
  GeneList <- GeneList[isAnalyzed]
}

## Total number of genes in GeneList
ng <- length(GeneList)
cat("Total number of genes in GeneList:", ng, "\n")


##----------
## Get info for the universe
##----------

if (EnrichTest) {
if (is.null(GeneUniverse)) {
  GeneUniverse <- GNN$GeneName
} else {
  GeneUniverse <- union(GeneUniverse, GeneList)
  UnivHasNeighbor <- GeneUniverse %in% GNN$GeneName
  if (any(!UnivHasNeighbor)) {
  warning(sum(!UnivHasNeighbor), "genes removed from universe because they have no neighborhood info\n")
  GeneUniverse <- GeneUniverse[UnivHasNeighbor]
  }
}
cat("Length of Gene Universe is", length(GeneUniverse), "\n")
}

##----------
## For upstream genes
##----------
cat("\nAnalysis of upstream gene orientation:\n")
cat("======================================\n")

##Remove genes that do not have upstream info
hasUP <- GNN %>%
    dplyr::slice(match(GeneList, GeneName)) %>%
    dplyr::mutate(hasUP = !is.na(UpstreamClass)) %>%
    dplyr::pull(hasUP)

if (any(!hasUP)) {
  cat(sum(!hasUP), "genes with no info for their upstream gene are removed\n")
}
gselUP <- GeneList[hasUP]
nUP <- length(gselUP) #Number of genes with UP info

if (!keepOther) {
    isOther <- GNN %>%
                   dplyr::slice(match(gselUP, GeneName)) %>%
                   dplyr::mutate(isOther = tolower(UpstreamClass) == "other") %>%
                   dplyr::pull(isOther)

    if (any(isOther)) {
        cat(sum(isOther), "genes with 'other' info for their upstream gene are removed\n")
        gselUP <- gselUP[!isOther]
        nUP <- length(gselUP)
    }
}

cat("Gene set for upstream gene analysis has", nUP, "genes\n")

##Counts and Percentages
up <- GNN %>%
            dplyr::filter(GeneName %in% gselUP) %>%
            dplyr::count(UpstreamClass) %>% as.data.frame

up$Percentage <- 100* up$n / nUP

##Hypergeomtric tests
if (EnrichTest) {
#Remove the genes without upstream data from universe
univhasUP <- GNN %>%
    dplyr::slice(match(GeneUniverse,GeneName)) %>%
    dplyr::mutate(univhasUP = !is.na(UpstreamClass)) %>%
    dplyr::pull(univhasUP)

gunivUP <- GeneUniverse[univhasUP]
nunivUP <- length(gunivUP)

if (any(!univhasUP)) {
  cat(sum(!univhasUP), "genes from universe have missing data for upstream gene\n")
}


if (!keepOther) {
    isOther <- GNN %>%
                   dplyr::slice(match(gunivUP, GeneName)) %>%
                   dplyr::mutate(isOther = tolower(UpstreamClass) == "other") %>%
                   dplyr::pull(isOther)

    if (any(isOther)) {
        cat(sum(isOther), "genes from universe with 'other' info for their upstream gene are removed\n")
        gunivUP <- gunivUP[!isOther]
        nunivUP <- length(gunivUP)
    }
}

cat("Universe for upstream gene analysis has", nunivUP, "genes\n")


#Get counts and percentages for Universe
upUniv <- GNN %>%
                dplyr::filter(GeneName %in% gunivUP) %>%
                dplyr::count(UpstreamClass) %>% as.data.frame
upUniv$Percentage <- 100* upUniv$n / nunivUP

up <- dplyr::right_join(up, upUniv, by="UpstreamClass", suffix=c("", "_Universe"))

up[is.na(up)] <- 0 #if some categories are not represented in GeneList, we give them 0 counts = 0 %

up$p.value <- apply(up[,-1],1, function(x){
                                           tablo <- matrix(c(x["n"],
                                                             x["n_Universe"]-x["n"],
                                                             nUP-x["n"],
                                                             nunivUP-nUP-x["n_Universe"]+x["n"]),
                                                           nrow=2, byrow=T)
                                           fisher.test(tablo,
                                                       alternative="greater")$p.value
                   })
}



##----------
## For downstream genes
##----------
cat("\nAnalysis of downstream gene orientation:\n")
cat("========================================\n")

##Remove genes that do not have downstream info
hasDN <- GNN %>%
            dplyr::slice(match(GeneList, GeneName)) %>%
            dplyr::mutate(hasDN = !is.na(DownstreamClass)) %>%
            dplyr::pull(hasDN)

if (any(!hasDN)) {
  cat(sum(!hasDN), "genes with no info for their downstream gene are removed\n")
}
gselDN <- GeneList[hasDN]
nDN <- length(gselDN) #Number of genes with DN info

if (!keepOther) {
    isOther <- GNN %>%
                   dplyr::slice(match(gselDN, GeneName)) %>%
                   dplyr::mutate(isOther = tolower(DownstreamClass)=="other") %>%
                   dplyr::pull(isOther)

    if (any(isOther)) {
        cat(sum(isOther), "genes with 'other' info for their downstream gene are removed\n")
        gselDN <- gselDN[!isOther]
        nDN <- length(gselDN)
    }
}

cat("Gene set for downstream gene analysis has", nDN, "genes\n")


##Counts and Percentages
dn <- GNN %>%
            dplyr::filter(GeneName %in% gselDN) %>%
            dplyr::count(DownstreamClass) %>% as.data.frame

dn$Percentage <- 100* dn$n / sum(dn$n)


##Hypergeomtric tests
if (EnrichTest) {
#Remove from universe genes without upstream data
univhasDN <- GNN %>%
                 dplyr::slice(match(GeneUniverse,GeneName)) %>%
                 dplyr::mutate(univhasDN = !is.na(DownstreamClass)) %>%
                 dplyr::pull(univhasDN)

gunivDN <- GeneUniverse[univhasDN]
nunivDN <- length(gunivDN)

if (any(!univhasDN)) {
  cat(sum(!univhasDN), "genes from universe have missing data for downstream gene\n")
}

if (!keepOther) {
    isOther <- GNN %>%
                   dplyr::slice(match(gunivDN, GeneName)) %>%
                   dplyr::mutate(isOther = tolower(DownstreamClass) == "other") %>%
                   dplyr::pull(isOther)

    if (any(isOther)) {
        cat(sum(isOther), "genes from universe with 'other' info for their downstream gene are removed\n")
        gunivDN <- gunivDN[!isOther]
        nunivDN <- length(gunivDN)
    }
}

cat("Universe for downstream gene analysis has", nunivDN, "genes\n\n")


#Get counts and percentages for Universe
dnUniv <- GNN %>%
                dplyr::filter(GeneName %in% gunivDN) %>%
                dplyr::count(DownstreamClass) %>% as.data.frame
dnUniv$Percentage <- 100* dnUniv$n / nunivDN

dn <- dplyr::right_join(dn, dnUniv, by="DownstreamClass", suffix=c("", "_Universe"))

dn[is.na(dn)] <- 0 #if some categories are not represented in GeneList, we give them 0 counts = 0 %

dn$p.value <- apply(dn[,-1],1, function(x){
                                           tablo <- matrix(c(x["n"],
                                                             x["n_Universe"]-x["n"],
                                                             nDN-x["n"],
                                                             nunivDN-nDN-x["n_Universe"]+x["n"]),
                                                           nrow=2, byrow=T)
                                           fisher.test(tablo,
                                                       alternative="greater")$p.value
                   })
}


res <- dplyr::bind_rows(up, dn) %>%
           dplyr::mutate(Orientation = ifelse(is.na(UpstreamClass),
                                              DownstreamClass,
                                              UpstreamClass),
                         Side = factor(c(rep("Upstream", nrow(up)),
                                         rep("Downstream", nrow(dn))),
                                       levels=c("Upstream", "Downstream"),
                                       ordered = TRUE))

if (EnrichTest) {
  res <- dplyr::select(res, Side, Orientation, n:p.value)
} else {
  res <- dplyr::select(res, Side, Orientation, n:Percentage)
}

return(res)
}
