
<!-- README.md is generated from README.Rmd. Please edit that file -->
GeneNeighborhood
================

The goal of GeneNeighborhood is to extract and analyze the orientation and proximity of upstream/downstream genes.

Installation
------------

You can install the released version of GeneNeighborhood from [CRAN](https://CRAN.R-project.org) with (NOT THERE YET!):

``` r
install.packages("GeneNeighborhood")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("pgpmartin/GeneNeighborhood")
```

Example
-------

Let's start by creating a GRanges with random features/genes:

``` r
set.seed(123)
gr <- GenomicRanges::GRanges(seqnames = "Chr1",
                             ranges = IRanges::IRanges(start = sample.int(10000, 676, 
                                                                          replace = TRUE),
                                                       width = sample.int(10, 676, 
                                                                          replace = TRUE),
                                                       names = paste0(rep(LETTERS[1:26], 
                                                                          times = 26),
                                                                      rep(LETTERS[1:26], 
                                                                          each = 26))),
                                                       strand = sample(c("+", "-"), 
                                                                       size = 676, 
                                                                       replace = TRUE))
```

For each feature/gene, we want to extract information about their upstream/downstream neighbors. We can do this using:

``` r
GeneNeighbors <- GeneNeighborhood::GetGeneNeighborhood(gr)
#> There are 96 genes (14.2%) that overlap with >1 gene
#> More than 10% of the genes overlap with multiple genes
```
