# Use k-means clustering to cluster local indicators of spatial association. For other clustering use lisa.

Use k-means clustering to cluster local indicators of spatial
association. For other clustering use lisa.

## Usage

``` r
lisaClust(
  cells,
  k = 2,
  r = NULL,
  imageID = "imageID",
  cellType = "cellType",
  spatialCoords = c("x", "y"),
  regionName = "region",
  cores = 1,
  window = "convex",
  window.length = NULL,
  whichParallel = "imimageID",
  sigma = NULL,
  lisaFunc = "K",
  minLambda = 0.05,
  BPPARAM = NULL,
  Rs = r
)
```

## Arguments

- cells:

  A SingleCellExperiment, SpatialExperiment or data frame that contains
  at least the variables x and y, giving the coordinates of each cell,
  imageID and cellType.

- k:

  The number of regions to cluster.

- r:

  A vector of the radii that the measures of association should be
  calculated.

- imageID:

  The column which contains image identifiers.

- cellType:

  The column which contains the cell types.

- spatialCoords:

  The columns which contain the x and y spatial coordinates.

- regionName:

  The output column for the lisaClust regions.

- cores:

  Number of cores to use for parallel processing, or a BiocParallel
  MulticoreParam or SerialParam object.

- window:

  Should the window around the regions be 'square', 'convex' or
  'concave'.

- window.length:

  A tuning parameter for controlling the level of concavity when
  estimating concave windows.

- whichParallel:

  Should the function use parallization on the imageID or the cellType.

- sigma:

  A numeric variable used for scaling when filting inhomogeneous
  L-curves.

- lisaFunc:

  Either "K" or "L" curve.

- minLambda:

  Minimum value for density for scaling when fitting inhomogeneous
  L-curves.

- BPPARAM:

  {DEPRECATED} A BiocParalell MulticoreParam or SerialParam object.

- Rs:

  A vector of the radii that the measures of association should be
  calculated.

## Value

A matrix of LISA curves

## Examples

``` r
## Generate toy data: two images, two cell types that sit in separate bands
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

# Cluster the cells into two regions
cells <- lisaClust(cells, k = 2)
#> Generating local indicators of spatial association.
table(cells$region, cells$cellType)
#>           
#>             c1  c2
#>   region_1   1 800
#>   region_2 799   0
```
