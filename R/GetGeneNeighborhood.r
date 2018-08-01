#' @title Extract info on upstream and downstream features/genes
#'
#' @description Get orientation and distance info on the upstream and downstream features/genes from a GRanges
#'
#' @param GeneGRanges A \code{GRanges} of feature/gene annotations
#'
#' @importFrom GenomicRanges sort seqnames strand start end precede follow distance countOverlaps findOverlaps ranges
#' @importFrom IRanges poverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom tibble tibble
#' @importFrom dplyr left_join recode
#'
#' @export
#'
#' @return A \code{\link{tibble}} with the following columns:
#' \itemize{
#'   \item \code{GeneName}: Name of the focus gene
#'   \item \code{Chr}: Seqnames of the focus gene
#'   \item \code{Strand}: Strand of the focus gene
#'   \item \code{Precede}: Name of the gene that the focus gene precedes (ignoring strand info)
#'   \item \code{Follow}: Name of the gene that the focus gene follows (ignoring strand info)
#'   \item \code{Upstream}: Name of the gene located upstream of the focus gene
#'   \item \code{Downstream}: Name of the gene located downstream of the focus gene
#'   \item \code{StrandUpstream}: Strand of the gene located upstream of the focus gene
#'   \item \code{StrandDownstream}: Strand of the gene located downstream of the focus gene
#'   \item \code{UpstreamOrientation}: Orientation of the upstream gene relative to the focus gene (S=same strand, O=Opposite strand or other)
#'   \item \code{DownstreamOrientation}: Orientation of the downstream gene relative to the focus gene (S, O or other)
#'   \item \code{UpstreamDistance}: Distance between the focus gene and its upstream neighbor
#'   \item \code{DownstreamDistance}: Distance between the focus gene and its downstream neighbor
#'   \item \code{CountOvlAnyStrand}: Number of genes with which the focus gene overlaps
#'   \item \code{ovlType}: Type of overlaps (see DETAILS)
#'   \item \code{ovlGene}: Name of the gene with which the focus gene overlaps (when it overlaps with a single gene)
#'   \item \code{NeighborClass}: Class of gene neighborhood (see DETAILS)
#'   \item \code{UpstreamClass}: Class of upstream Neighborhood (takes into account the presence of an overlap)
#'   \item \code{DownstreamClass}: Class of downstream Neighborhood (takes into account the presence of an overlap)
#'   \item \code{GenePair}: Type of gene pair formed between the focus gene and the gene that follows (see DETAILS)
#'   \item \code{Distance2Pair}: Intergenic distance between the focus gene and the gene that follows
#' }
#'
#' @section DETAILS:
#' Up/DownstreamOrientation columns do not precise if the up/downstream gene overlaps with the focus gene.
#'
#' A distance of 0 in Up/DownstreamDistance columns can indicate that the genes overlap or that they are adjacent.
#'
#' The type of overlaps in the ovlType column are coded as follow:
#' \itemize{
#'   \item c. The focus gene overlaps with more than 1 gene ("complex" overlap pattern)
#'   \item I. The focus gene shares the same borders (start AND end) with a gene on the same strand
#'   \item i. The focus gene shares the same borders (start AND end) with a gene on the opposite strand
#'   \item H. The focus gene "hosts" (=contains entirely) a gene on the same strand
#'   \item h. The focus gene "hosts" (=contains entirely) a gene on the opposite strand
#'   \item G. The focus gene is the "guest" of (=is entirely contained within) a gene that is on the same strand
#'   \item g. The focus gene is the "guest" of (=is entirely contained within) a gene that is on the opposite strand
#'   \item U. The 5' end of the focus gene (TSS) overlaps with a gene that is on the same strand
#'   \item u. The 5' end of the focus gene (TSS) overlaps with a gene that is on the opposite strand
#'   \item D. The 3' end of the focus gene (TES) overlaps with a gene that is on the same strand
#'   \item d. The 3' end of the focus gene (TES) overlaps with a gene that is on the opposite strand
#' }
#'
#'  The Neighborclass column is defined as follow:
#' \itemize{
#'   \item{In the absence of overlap:}{ A 2-letter code using S (Same Strand) and O (Opposite Strand) to indicate the orientation of the uptstream (first letter) and downstream genes (second letter)}
#'   \item{When ovlType is u/U or d/D:}{ ovlType is combined with the information on the downstream (resp. upstream) gene to form the code (e.g. uS or uO)}
#'   \item{In all other cases:}{ NeighborClass is equal to ovlType}
#'}
#'
#'  The codes in the GenePair column are defined as follow:
#' \itemize{
#'   \item{H2H:}{ Head-to-Head orientation: <---   --->}
#'   \item{T2H:}{ Tail-to-Head orientation: --->   --->}
#'   \item{T2T:}{ Tail-to-Tail orientation: --->   <---}
#' }
#'
#' @author Pascal GP Martin
#'
#' @examples \dontrun{
#' library(TxDb.Athaliana.BioMart.plantsmart25)
#' GeneNeighbors <- GetGeneNeighborhood(GenomicFeatures::genes(TxDb.Athaliana.BioMart.plantsmart25))
#' # There are 155 genes (0.461%) that overlap with >1 gene
#' table(GeneNeighbors$NeighborClass)
#' #
#' #   c    g    G    h    H    i    I   Od   OD   OO   OS   Sd   SD   SO   SS   uO   UO   uS   US
#' # 155  257  137  190   92    4   64 1448   52 6203 7461 1415   65 6237 9358  152   37  196   65
#' table(table(GeneNeighbors$GenePair))
#' #
#' #  H2H   T2H   T2T
#' # 7762 17173  7754
#' }
#'
#'
#' set.seed(123)
#' #The package includes a \code{GRanges} named Genegr with 676 random genes on a single chromosome:
#'   Genegr
#'
#' #Extract info on the neighbors of these genes:
#'   GeneNeighbors <- GetGeneNeighborhood(Genegr)
#'
#' #Classes of neighborhoods:
#'   table(GeneNeighbors$NeighborClass)
#' #Classes of gene pairs:
#'   table(GeneNeighbors$GenePair)
#'
#' #UpstreamOrientation and DownstreamOrientation do not take into account potential overlaps of the upstream/downstream gene:
#'   table(GeneNeighbors$NeighborClass[GeneNeighbors$UpstreamOrientation=="O"])
#' #This information is present in the Up/DownstreamClass columns:
#'   table(GeneNeighbors$NeighborClass[GeneNeighbors$DownstreamClass=="OppositeOverlap"])
#'

GetGeneNeighborhood <- function(GeneGRanges) {

if (!is(GeneGRanges, "GRanges")) {
  stop("GeneGRanges should be a GRanges object")
}

  if (is.null(names(GeneGRanges))) {
  stop("No names for the GRanges provided")
}

if (!(length(unique(names(GeneGRanges))) == length(GeneGRanges))) {
  stop("Gene names of GeneGRanges are not unique")
}

#Sort the input GeneGranges (for reproducibility)
GeneGRanges <- GenomicRanges::sort(GeneGRanges)

##----------------------
## Extract raw info on upstream/downstream genes and on number of overlaps
##----------------------
#Precede/Follow info
res <- tibble::tibble(GeneName = as.character(names(GeneGRanges)),
                      Chr = as.character(GenomicRanges::seqnames(GeneGRanges)),
                      Strand = as.character(GenomicRanges::strand(GeneGRanges)),
                      Precede = as.character(names(GeneGRanges)[GenomicRanges::precede(GeneGRanges, select="first", ignore.strand=TRUE)]),
                      Follow = as.character(names(GeneGRanges)[GenomicRanges::follow(GeneGRanges, select="last", ignore.strand=TRUE)]))
#Note that if a gene has 2 neighbors at the exact same distance, only one will be selected by precede/follow
#For these genes, different results could be obtained with different sorting of the GeneGRanges

##Strand information
strgn <- as.character(GenomicRanges::strand(GeneGRanges))
names(strgn) <- names(GeneGRanges)


#Upstream/Downstream info based on the strand of the focus gene
res$Upstream <- res$Follow
res$Upstream[strgn == "*"] <- NA
res$Upstream[strgn == "-"] <- res$Precede[strgn == "-"]

res$Downstream <- res$Precede
res$Downstream[strgn == "*"] <- NA
res$Downstream[strgn == "-"] <- res$Follow[strgn == "-"]


#Strand of the upstream/downstream gene
res$StrandUpstream <- strgn[res$Upstream]
res$StrandDownstream <- strgn[res$Downstream]

#Define if the upstream/downstream gene is on the same (S) or opposite (O) orientation as the focus gene
res$UpstreamOrientation <- ifelse(res$StrandUpstream == res$Strand, "S", "O")
res$DownstreamOrientation <- ifelse(res$StrandDownstream == res$Strand, "S", "O")

#Distances to upstream and downstream genes
res$UpstreamDistance <- rep(NA, nrow(res))
res$UpstreamDistance[!is.na(res$Upstream)] <- GenomicRanges::distance(GeneGRanges[!is.na(res$Upstream)],
                                                                      GeneGRanges[res$Upstream[!is.na(res$Upstream)]],
                                                                      ignore.strand = TRUE)

res$DownstreamDistance <- rep(NA, nrow(res))
res$DownstreamDistance[!is.na(res$Downstream)] <- GenomicRanges::distance(GeneGRanges[!is.na(res$Downstream)],
                                                                          GeneGRanges[res$Downstream[!is.na(res$Downstream)]],
                                                                          ignore.strand = TRUE)

#Number of overlapping genes on any strand:
res$CountOvlAnyStrand <- GenomicRanges::countOverlaps(GeneGRanges, ignore.strand = TRUE)-1


##----------------------
## Annotate the genes that overlap with a single gene
##----------------------
## Get the index of the genes that have a single overlapping gene:
WhichSingleOVL <- which(!is.na(res$CountOvlAnyStrand) & res$CountOvlAnyStrand == 1)
## Find all overlaps between all genes:
AllFOV <- GenomicRanges::findOverlaps(GeneGRanges, ignore.strand = TRUE)
## Extract the hits corresponding to single overlaps:
SingleOVL <- AllFOV[S4Vectors::queryHits(AllFOV) %in% WhichSingleOVL &
                      S4Vectors::queryHits(AllFOV) != S4Vectors::subjectHits(AllFOV)]
stopifnot(identical(S4Vectors::queryHits(SingleOVL),
                    WhichSingleOVL))

## Get the corresponding GRanges:
GenesWithSingleOVL <- GeneGRanges[queryHits(SingleOVL)]
SingleOVLappingGenes <- GeneGRanges[subjectHits(SingleOVL)]

isSameBorders <- IRanges::poverlaps(GenomicRanges::ranges(GenesWithSingleOVL),
                                    GenomicRanges::ranges(SingleOVLappingGenes),
                                    type="equal")
isGuest <- IRanges::poverlaps(GenomicRanges::ranges(GenesWithSingleOVL),
                              GenomicRanges::ranges(SingleOVLappingGenes),
                              type="within") & !isSameBorders
isHost <- IRanges::poverlaps(GenomicRanges::ranges(SingleOVLappingGenes),
                             ranges(GenesWithSingleOVL),
                             type="within") & !isSameBorders
isUpstream <- ifelse(as.logical(GenomicRanges::strand(GenesWithSingleOVL) == "+"),
                     GenomicRanges::start(SingleOVLappingGenes) < GenomicRanges::start(GenesWithSingleOVL) &
                         GenomicRanges::end(SingleOVLappingGenes) < GenomicRanges::end(GenesWithSingleOVL),
                     GenomicRanges::end(SingleOVLappingGenes) > GenomicRanges::end(GenesWithSingleOVL) &
                         GenomicRanges::start(SingleOVLappingGenes) > GenomicRanges::start(GenesWithSingleOVL))
isDownstream <- ifelse(as.logical(GenomicRanges::strand(GenesWithSingleOVL) == "+"),
                       GenomicRanges::start(SingleOVLappingGenes) > GenomicRanges::start(GenesWithSingleOVL) &
                           GenomicRanges::end(SingleOVLappingGenes) > GenomicRanges::end(GenesWithSingleOVL),
                       GenomicRanges::end(SingleOVLappingGenes) < GenomicRanges::end(GenesWithSingleOVL) &
                           GenomicRanges::start(SingleOVLappingGenes) < GenomicRanges::start(GenesWithSingleOVL))
isSameStrand <- as.logical(GenomicRanges::strand(GenesWithSingleOVL) == GenomicRanges::strand(SingleOVLappingGenes))


#Build result table
tov <- tibble::tibble(GeneName = names(GenesWithSingleOVL),
                      ovlType = NA,
                      ovlGene = names(SingleOVLappingGenes))

tov$ovlType[isSameBorders] <- "i"
tov$ovlType[isHost] <- "h"
tov$ovlType[isGuest] <- "g"
tov$ovlType[isUpstream] <- "u"
tov$ovlType[isDownstream] <- "d"
tov$ovlType[isSameStrand] <- toupper(tov$ovlType[isSameStrand])


##----------------------
## Combine the information in the NeighborClass column
##----------------------
# Add the info on the gene with a single overlap:
res <- dplyr::left_join(res, tov, by="GeneName")

#Define the NeighborhoodClass column (combining upstream/downstream genomic context)
res$NeighborClass <- paste0(res$UpstreamOrientation,res$DownstreamOrientation)
res$NeighborClass[res$CountOvlAnyStrand>=2] <- "c"

isuSuO <- !is.na(res$ovlType) & tolower(res$ovlType)=="u"
res$NeighborClass[isuSuO] <- paste0(res$ovlType[isuSuO], res$DownstreamOrientation[isuSuO])

isSdOd <- !is.na(res$ovlType) & tolower(res$ovlType)=="d"
res$NeighborClass[isSdOd] <- paste0(res$UpstreamOrientation[isSdOd], res$ovlType[isSdOd])

isIHG <- !is.na(res$ovlType) & tolower(res$ovlType) %in% c("i", "h", "g")
res$NeighborClass[isIHG] <- res$ovlType[isIHG]

res$NeighborClass[is.na(res$UpstreamOrientation) | is.na(res$DownstreamOrientation)] <- NA

##----------------------
## Correct the Upstream/Downstream columns to take into account the genes that overlap with other genes
## A/ For genes that have the same border of another gene ("i/I") or that host another gene ("h/H"), we keep the upstream/downstream data as is
## B/ For genes that have an overlapping gene at their 5' or 3' end only, we redefine the upstream/downstream gene accordingly and set the distance to 0
## C/ For genes that are guest of another gene ("g/G") or that interact with >1 gene ("c"), we set the upstream and downstream gene as NA and the orientation as "other"
##----------------------

#--------
# B/ Genes with overlapping gene at the 5' (upstream) or 3' end (downstream)
# Correct the upstream columns to take into account the genes with a single overlap
#--------
## Upstream gene:
isOVLUpstream <- !(is.na(res$NeighborClass)) & (toupper(res$NeighborClass) %in% c("UO", "US"))
res$Upstream[isOVLUpstream] <- res$ovlGene[isOVLUpstream]
res$UpstreamDistance[isOVLUpstream] <- 0
res$StrandUpstream[isOVLUpstream] <- strgn[res$ovlGene[isOVLUpstream]]
res$UpstreamOrientation <- ifelse(res$StrandUpstream==res$Strand, "S", "O")
#~ message("The UpstreamOrientation column does not distinguish whether the upstream gene is overlapping or not")

## Downstream gene:
isOVLDownstream <- !(is.na(res$NeighborClass)) & (toupper(res$NeighborClass) %in% c("OD", "SD"))
res$Downstream[isOVLDownstream] <- res$ovlGene[isOVLDownstream]
res$DownstreamDistance[isOVLDownstream] <- 0
res$StrandDownstream[isOVLDownstream] <- strgn[res$ovlGene[isOVLDownstream]]
res$DownstreamOrientation <- ifelse(res$StrandDownstream==res$Strand, "S", "O")
#~ message("The DownstreamOrientation column does not distinguish whether the upstream gene is overlapping or not")

#--------
# C/ Genes with >1 overlapping genes and "guest" genes ("g/G")
#--------
## Genes with multiple overlaps:
isMultiOVL <- !is.na(res$CountOvlAnyStrand) & res$CountOvlAnyStrand > 1
percMultiOVL <- 100*mean(isMultiOVL)
cat(paste0("There are ",
           sum(isMultiOVL),
           " genes (",
           ifelse(percMultiOVL>0.1,
                  signif(percMultiOVL, 3),
                  signif(percMultiOVL, 1)),
           "%) that overlap with >1 gene\n"))
if (percMultiOVL>=10) {
    message("More than 10% of the genes overlap with multiple genes")
}

## "Guest" genes:
isGuestGene <- !is.na(res$NeighborClass) & toupper(res$NeighborClass)=="G"

#Upstream info:
res$Upstream[isMultiOVL | isGuestGene] <- NA
res$UpstreamDistance[isMultiOVL | isGuestGene] <- NA
res$StrandUpstream[isMultiOVL | isGuestGene] <- NA
res$UpstreamOrientation[isMultiOVL | isGuestGene] <- "other"
#Downstream info:
res$Downstream[isMultiOVL | isGuestGene] <- NA
res$DownstreamDistance[isMultiOVL | isGuestGene] <- NA
res$StrandDownstream[isMultiOVL | isGuestGene] <- NA
res$DownstreamOrientation[isMultiOVL | isGuestGene] <- "other"
#Genes with multiple overlaps or that are within another gene have NA values for both upstream and downstream genes

##----------------------
## Create UpstreamClass/DownstreamClass columns in which we keep the info on overlapping upstream/downstream genes
##----------------------
#Create simplified columns for classes of upstream/downstream orientations
res$UpstreamClass <- dplyr::recode(res$UpstreamOrientation,
                                   "S" = "SameStrand",
                                   "O" = "OppositeStrand")
res$UpstreamClass[!is.na(res$NeighborClass) & res$NeighborClass %in% c("uS", "uO")] <- "OppositeOverlap"
res$UpstreamClass[!is.na(res$NeighborClass) & res$NeighborClass %in% c("US", "UO")] <- "SameOverlap"

res$DownstreamClass <- dplyr::recode(res$DownstreamOrientation,
                                     "S" = "SameStrand",
                                     "O" = "OppositeStrand")
res$DownstreamClass[!is.na(res$NeighborClass) & res$NeighborClass %in% c("Sd", "Od")] <- "OppositeOverlap"
res$DownstreamClass[!is.na(res$NeighborClass) & res$NeighborClass %in% c("SD", "OD")] <- "SameOverlap"


##----------------------
## Add information on gene pairs
##----------------------

##For each gene, indicate the kind of gene pair that it forms with the following gene
### GenePair definitions:
#  T2H:   ---->  ---->
#  H2H:   <----  ---->
#  T2T:   ---->  <----

res$GenePair <- ifelse(!is.na(res$NeighborClass) & nchar(res$NeighborClass)==2 & res$Strand=="+",
                       ifelse(res$DownstreamOrientation=="O", "T2T", "T2H"),
                       ifelse(!is.na(res$NeighborClass) & nchar(res$NeighborClass)==2 & res$Strand=="-",
                              ifelse(res$UpstreamOrientation=="O", "H2H", "T2H"),
                              NA))

## And the distance within the pair:
res$Distance2Pair <- ifelse(!is.na(res$GenePair),
                            ifelse(res$Strand=="+",
                                   res$DownstreamDistance,
                                   res$UpstreamDistance),
                            NA)


# Return results
return(res)
}
