# Plot heatmap of cell type enrichment for lisaClust regions

Plot heatmap of cell type enrichment for lisaClust regions

## Usage

``` r
regionMap(
  cells,
  type = "bubble",
  cellType = "cellType",
  region = "region",
  limit = c(0.33, 3),
  ...
)
```

## Arguments

- cells:

  SingleCellExperiment, SpatialExperiment or data.frame

- type:

  Make a "bubble" or "heatmap" plot, or return the relative frequencies
  as a "table".

- cellType:

  The column storing the cell types

- region:

  The column storing the regions

- limit:

  limits to the lower and upper relative frequencies

- ...:

  Any arguments to be passed to the pheatmap package

## Value

A bubble plot or heatmap, or with `type = "table"` a matrix of the
relative frequencies with one row per cell type and one column per
region: how much more often the cell type is found in the region than if
cell types were spread evenly over the regions.

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

regionMap(cells)

regionMap(cells, type = "table")
#>     
#>      region_1 region_2
#>   c1        0        2
#>   c2        2        0
```
