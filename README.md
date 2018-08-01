
<!-- README.md is generated from README.Rmd. Please edit that file -->
GeneNeighborhood
================

The goal of GeneNeighborhood is to extract and analyze the orientation and proximity of upstream/downstream genes.

Installation
------------

(NOT THERE YET!) You can install the released version of GeneNeighborhood from [CRAN](https://CRAN.R-project.org) with:

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

Load the library:

``` r
library(GeneNeighborhood)
```

To run examples, the package includes a *GRanges* named **Genegr** that contains 676 genes with random coordinates on a single chromosome *Chr1*:

``` r
Genegr
#> GRanges object with 676 ranges and 0 metadata columns:
#>      seqnames       ranges strand
#>         <Rle>    <IRanges>  <Rle>
#>   AA     Chr1 [2876, 2884]      +
#>   AB     Chr1 [7884, 7886]      +
#>   AC     Chr1 [4090, 4094]      +
#>   AD     Chr1 [8831, 8835]      +
#>   AE     Chr1 [9405, 9412]      +
#>   ..      ...          ...    ...
#>   ZV     Chr1 [5214, 5220]      +
#>   ZW     Chr1 [ 863,  868]      +
#>   ZX     Chr1 [2831, 2836]      +
#>   ZY     Chr1 [4205, 4205]      -
#>   ZZ     Chr1 [5873, 5873]      -
#>   -------
#>   seqinfo: 1 sequence from mock genome
```

For each feature/gene, we extract information (orientation and distance, potential overlaps) about their upstream/downstream neighbors with:

``` r
GeneNeighbors <- GetGeneNeighborhood(Genegr)
#> There are 96 genes (14.2%) that overlap with >1 gene
#> More than 10% of the genes overlap with multiple genes
```

We define a random set of 100 genes:

``` r
set.seed(123)
randGenes <- names(Genegr)[sample.int(676, 100)]
```

We extract statistics about the orientation of their neighbors using:

``` r
## Neighbors Orientation Statistics:
NOS <- AnalyzeNeighborsOrientation(randGenes, 
                                   GeneNeighborhood = GeneNeighbors)
#> Total number of genes in GeneList: 100 
#> Length of Gene Universe is 676 
#> 
#> Analysis of upstream gene orientation:
#> ======================================
#> Gene set for upstream gene analysis has 100 genes
#> 2 genes from universe have missing data for upstream gene
#> Universe for upstream gene analysis has 674 genes
#> 
#> Analysis of downstream gene orientation:
#> ========================================
#> Gene set for downstream gene analysis has 100 genes
#> 1 genes from universe have missing data for downstream gene
#> Universe for downstream gene analysis has 675 genes
```

By default all genes are used as a universe and an enrichment test is performed.
By default, the function also analyzes the *"other"* orientation which may be hard to interpret. We can remove this orientation using:

``` r
NOS <- AnalyzeNeighborsOrientation(randGenes, 
                                   GeneNeighborhood = GeneNeighbors,
                                   keepOther = FALSE)
#> Total number of genes in GeneList: 100 
#> Length of Gene Universe is 676 
#> 
#> Analysis of upstream gene orientation:
#> ======================================
#> 21 genes with 'other' info for their upstream gene are removed
#> Gene set for upstream gene analysis has 79 genes
#> 2 genes from universe have missing data for upstream gene
#> 157 genes from universe with 'other' info for their upstream gene are removed
#> Universe for upstream gene analysis has 517 genes
#> 
#> Analysis of downstream gene orientation:
#> ========================================
#> 21 genes with 'other' info for their downstream gene are removed
#> Gene set for downstream gene analysis has 79 genes
#> 1 genes from universe have missing data for downstream gene
#> 157 genes from universe with 'other' info for their downstream gene are removed
#> Universe for downstream gene analysis has 518 genes
```

We obtain the following table:

| Side       | Orientation     |    n|  Percentage|  n\_Universe|  Percentage\_Universe|  p.value|
|:-----------|:----------------|----:|-----------:|------------:|---------------------:|--------:|
| Upstream   | OppositeOverlap |    4|        5.06|           35|                  6.77|    0.810|
| Upstream   | OppositeStrand  |   38|       48.10|          236|                 45.65|    0.360|
| Upstream   | SameOverlap     |    4|        5.06|           31|                  6.00|    0.730|
| Upstream   | SameStrand      |   33|       41.77|          215|                 41.59|    0.530|
| Downstream | OppositeOverlap |    5|        6.33|           16|                  3.09|    0.081|
| Downstream | OppositeStrand  |   40|       50.63|          248|                 47.88|    0.340|
| Downstream | SameOverlap     |    4|        5.06|           31|                  5.98|    0.720|
| Downstream | SameStrand      |   30|       37.97|          223|                 43.05|    0.870|
