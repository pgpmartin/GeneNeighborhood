
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

Load the library:

``` r
library(GeneNeighborhood)
```

To run examples, the package includes a *GRanges* named gr that contains 676 genes with random coordinates on a single chromosome *Chr1*:

``` r
gr
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
GeneNeighbors <- GetGeneNeighborhood(gr)
#> There are 96 genes (14.2%) that overlap with >1 gene
#> More than 10% of the genes overlap with multiple genes
```
