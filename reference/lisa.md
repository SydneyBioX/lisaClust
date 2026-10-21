# Generate local indicators of spatial association

Generate local indicators of spatial association

## Usage

``` r
lisa(
  cells,
  r = NULL,
  imageID = "imageID",
  cellType = "cellType",
  spatialCoords = c("x", "y"),
  cores = 1,
  window = "convex",
  window.length = NULL,
  whichParallel = "imageID",
  sigma = NULL,
  lisaFunc = "K",
  minLambda = 0.05,
  BPPARAM = NULL,
  Rs = r,
  cellWeights = NULL
)
```

## Arguments

- cells:

  A SingleCellExperiment, SpatialExperiment or data frame that contains
  at least the variables x and y, giving the coordinates of each cell,
  imageID and cellType.

- r:

  A vector of the radii that the measures of association should be
  calculated.

- imageID:

  The column which contains image identifiers.

- cellType:

  The column which contains the cell types.

- spatialCoords:

  The columns which contain the x and y spatial coordinates.

- cores:

  Number of cores to use for parallel processing, or a BiocParallel
  param object (such as MulticoreParam, SnowParam or SerialParam).

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

  {DEPRECATED} A BiocParallel param object.

- Rs:

  {DEPRECATED} A vector of the radii that the measures of association
  should be calculated.

- cellWeights:

  Optional cell weights in \\0, 1\\, in the order of \`cells\`: a
  numeric vector with one weight per cell for its own cell type (for
  example a quality score), or a numeric matrix (or sparse \`Matrix\`)
  of soft cell type membership with one row per cell and one named
  column per label (for example a classifier's confidence). The curves
  are then one per label: each neighbour counts towards a label in
  proportion to its weight, and a label's density is its summed weight
  per unit area. Labels may overlap, for example a group and its
  subtypes. Weights of 1 give the curves of the cell types.

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

# Generate LISA curves
lisaCurves <- lisa(cells, r = c(10, 20, 50))
#> Generating local indicators of spatial association.

# Cluster the LISA curves
kM <- kmeans(lisaCurves, 2)
```
