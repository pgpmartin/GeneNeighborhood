% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getDistSide.R
\name{getDistSide}
\alias{getDistSide}
\title{Filter a gene set and extract upstream or downstream distances
       do their neighbors}
\usage{
getDistSide(GNN, glist, Side = c("Upstream", "Downstream"))
}
\arguments{
\item{GNN}{A \code{tibble} obtained with the
\code{\link{getGeneNeighborhood}} function.
Can also be any data.frame with the relevant information.}

\item{glist}{A character vector of gene IDs (must match GNN$GeneName)}

\item{Side}{One of 'Upstream' or 'Downstream'}
}
\value{
A tibble with the following columns:
\itemize{
  \item GeneName. ID of the focus gene.
  \item Neighbor. ID of the gene neighbor
  \item Orientation. Orientation of the neighbor (same or opposite strand).
  \item Side. Upstream or Downstream.
  \item Distance in bp
}
}
\description{
Extract distances to the upstream or downstream neighbors
}
\section{DETAILS}{

Note that the function removes overlapping upstream/downstream genes
but not adjacent genes (for which distance = 0).
}

\examples{
## Obtain gene neighborhood information:
  GeneNeighbors <- getGeneNeighborhood(Genegr)
## Get a (random) set of (10) genes:
  set.seed(123)
  randGenes <- sample(names(Genegr), 10)
## Extract distances to the upstream neighbors for this set:
getDistSide(GeneNeighbors, randGenes, "Upstream")
## And for the downstream genes:
getDistSide(GeneNeighbors, randGenes, "Downstream")
}
\seealso{
\code{\link{getGeneNeighborhood}}
         \code{\link{dist2Neighbors}}
}
