# Name regions by the cell types they are enriched for

The regions found by
[`lisaClust()`](https://github.com/ellispatrick/lisaClust/reference/lisaClust.md)
are numbered arbitrarily. `nameRegions()` gives each region the name of
the marker cell type it is most enriched for, as measured by
`regionMap(type = "table")`. Regions most enriched for the same marker
get the same name, so several regions can be merged into one domain.

## Usage

``` r
nameRegions(
  cells,
  markers,
  region = "region",
  cellType = "cellType",
  regionName = region
)
```

## Arguments

- cells:

  A SingleCellExperiment, SpatialExperiment or data frame with a region
  column.

- markers:

  A character vector of cell types. Its names, when given, are used as
  the names of the regions; otherwise the cell types are.

- region:

  The column storing the regions.

- cellType:

  The column storing the cell types.

- regionName:

  The output column for the named regions; by default the regions are
  replaced.

## Value

`cells` with the named regions in the column `regionName`.

## Examples

``` r
set.seed(51773)
x <- round(c(
  runif(200), runif(200) + 1, runif(200) + 2, runif(200) + 3,
  runif(200) + 3, runif(200) + 2, runif(200) + 1, runif(200)
), 4) * 100
y <- round(c(
  runif(200), runif(200) + 1, runif(200) + 2, runif(200) + 3,
  runif(200), runif(200) + 1, runif(200) + 2, runif(200) + 3
), 4) * 100
cellType <- factor(paste("c", rep(rep(c(1:2), rep(200, 2)), 4), sep = ""))
imageID <- rep(c("s1", "s2"), c(800, 800))
cells <- data.frame(x, y, cellType, imageID)
cells <- lisaClust(cells, k = 2)
#> Generating local indicators of spatial association.

cells <- nameRegions(cells, c(first = "c1", second = "c2"), regionName = "domain")
table(cells$domain, cells$cellType)
#>         
#>           c1  c2
#>   first  800   0
#>   second   0 800
```
