
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

### Obtain data on the gene neighbors

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
GeneNeighbors <- getGeneNeighborhood(Genegr)
#> There are 96 genes (14.2%) that overlap with >1 gene
#> More than 10% of the genes overlap with multiple genes
```

### Analyze the orientation of the genes' neighbors

We define a random set of 100 genes:

``` r
set.seed(123)
randGenes <- names(Genegr)[sample.int(676, 100)]
```

We extract statistics about the orientation of their neighbors using:

``` r
## Neighbors Orientation Statistics:
NOS <- analyzeNeighborsOrientation(randGenes, 
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
NOS <- analyzeNeighborsOrientation(randGenes, 
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

We can plot the corresponding percentages using:

``` r
plotNeighborsOrientation(NOS)
```

<img src="man/figures/README-plotOrientationNOS-1.png" width="100%" />

### Analyze the proximity of the genes' neighbors

We can analyze specifically the distances to the upstream genes with:

``` r
randUpstreamDist <- statDistanceSide(GeneNeighborhood = GeneNeighbors,
                                     glist = randGenes,
                                     Side = "Upstream")
#> 
#> Analysis of Distances to upstream gene:
#> ========================================
#> 21 genes with undefined upstream gene are removed
#> 8 genes with an overlapping upstream gene are removed from GeneList
#> GeneList for upstream gene analysis has 71 genes
#> 159 genes from GeneUniverse with undefined upstream gene are removed
#> 66 genes from GeneUniverse with an overlapping upstream gene are removed
#> Universe for upstream gene analysis has 451 genes
```

Which gives the following (simplified) table:

| GeneGroup    | SideClass      |    n|  Median|   Mean|     SD|  KS.pvalue|  Wilcox.pvalue|  Independ.pvalue|
|:-------------|:---------------|----:|-------:|------:|------:|----------:|--------------:|----------------:|
| GeneList     | OppositeStrand |   38|     9.5|  14.47|  13.81|       0.65|           0.63|             0.52|
| GeneList     | SameStrand     |   33|    10.0|  13.00|   9.80|       0.91|           0.98|             0.51|
| GeneUniverse | OppositeStrand |  236|    10.0|  13.30|  12.43|         NA|             NA|               NA|
| GeneUniverse | SameStrand     |  215|    11.0|  14.47|  13.80|         NA|             NA|               NA|

Or we directly analyze the distances to both upstream and downstream genes:

``` r
randDist <- analyzeNeighborsDistance(GeneList = randGenes,
                                     GeneNeighborhood = GeneNeighbors)
#> 
#> Prefiltering of GeneList and GeneUniverse:
#> ==========================================
#> Total number of genes in GeneList: 100 
#> Total number of genes in GeneUniverse (including GeneList) is 676 
#> 
#> Analysis of Distances to upstream gene:
#> ========================================
#> 21 genes with undefined upstream gene are removed
#> 8 genes with an overlapping upstream gene are removed from GeneList
#> GeneList for upstream gene analysis has 71 genes
#> 159 genes from GeneUniverse with undefined upstream gene are removed
#> 66 genes from GeneUniverse with an overlapping upstream gene are removed
#> Universe for upstream gene analysis has 451 genes
#> 
#> Analysis of Distances to downstream gene:
#> ========================================
#> 21 genes with undefined downstream gene are removed
#> 9 genes with an overlapping downstream gene are removed from GeneList
#> GeneList for downstream gene analysis has 70 genes
#> 158 genes from GeneUniverse with undefined downstream gene are removed
#> 47 genes from GeneUniverse with an overlapping downstream gene are removed
#> Universe for downstream gene analysis has 471 genes
```

Which gives the following (simplified) table:

<table>
<colgroup>
<col width="12%" />
<col width="10%" />
<col width="14%" />
<col width="4%" />
<col width="7%" />
<col width="6%" />
<col width="6%" />
<col width="9%" />
<col width="13%" />
<col width="15%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">GeneGroup</th>
<th align="left">Side</th>
<th align="left">Orientation</th>
<th align="right">n</th>
<th align="right">Median</th>
<th align="right">Mean</th>
<th align="right">SD</th>
<th align="right">KS.pvalue</th>
<th align="right">Wilcox.pvalue</th>
<th align="right">Independ.pvalue</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">GeneList</td>
<td align="left">Upstream</td>
<td align="left">OppositeStrand</td>
<td align="right">38</td>
<td align="right">9.5</td>
<td align="right">14.47</td>
<td align="right">13.81</td>
<td align="right">0.65</td>
<td align="right">0.63</td>
<td align="right">0.520</td>
</tr>
<tr class="even">
<td align="left">GeneList</td>
<td align="left">Upstream</td>
<td align="left">SameStrand</td>
<td align="right">33</td>
<td align="right">10.0</td>
<td align="right">13.00</td>
<td align="right">9.80</td>
<td align="right">0.91</td>
<td align="right">0.98</td>
<td align="right">0.510</td>
</tr>
<tr class="odd">
<td align="left">GeneUniverse</td>
<td align="left">Upstream</td>
<td align="left">OppositeStrand</td>
<td align="right">236</td>
<td align="right">10.0</td>
<td align="right">13.30</td>
<td align="right">12.43</td>
<td align="right">NA</td>
<td align="right">NA</td>
<td align="right">NA</td>
</tr>
<tr class="even">
<td align="left">GeneUniverse</td>
<td align="left">Upstream</td>
<td align="left">SameStrand</td>
<td align="right">215</td>
<td align="right">11.0</td>
<td align="right">14.47</td>
<td align="right">13.80</td>
<td align="right">NA</td>
<td align="right">NA</td>
<td align="right">NA</td>
</tr>
<tr class="odd">
<td align="left">GeneList</td>
<td align="left">Downstream</td>
<td align="left">OppositeStrand</td>
<td align="right">40</td>
<td align="right">7.0</td>
<td align="right">15.20</td>
<td align="right">18.69</td>
<td align="right">0.62</td>
<td align="right">0.50</td>
<td align="right">0.810</td>
</tr>
<tr class="even">
<td align="left">GeneList</td>
<td align="left">Downstream</td>
<td align="left">SameStrand</td>
<td align="right">30</td>
<td align="right">15.5</td>
<td align="right">18.70</td>
<td align="right">15.54</td>
<td align="right">0.29</td>
<td align="right">0.10</td>
<td align="right">0.098</td>
</tr>
<tr class="odd">
<td align="left">GeneUniverse</td>
<td align="left">Downstream</td>
<td align="left">OppositeStrand</td>
<td align="right">248</td>
<td align="right">10.0</td>
<td align="right">14.65</td>
<td align="right">15.66</td>
<td align="right">NA</td>
<td align="right">NA</td>
<td align="right">NA</td>
</tr>
<tr class="even">
<td align="left">GeneUniverse</td>
<td align="left">Downstream</td>
<td align="left">SameStrand</td>
<td align="right">223</td>
<td align="right">12.0</td>
<td align="right">14.85</td>
<td align="right">13.69</td>
<td align="right">NA</td>
<td align="right">NA</td>
<td align="right">NA</td>
</tr>
</tbody>
</table>
