# Calculate the inhomogenous local K function.

Calculate the inhomogenous local K function.

## Usage

``` r
inhomLocalK(
  data,
  Rs = c(20, 50, 100, 200),
  sigma = 10000,
  window = "convex",
  window.length = NULL,
  minLambda = 0.05,
  lisaFunc = "K"
)
```

## Arguments

- data:

  The data.

- Rs:

  A vector of the radii that the measures of association should be
  calculated.

- sigma:

  A numeric variable used for scaling when filting inhomogeneous
  L-curves.

- window:

  Should the window around the regions be 'square', 'convex' or
  'concave'.

- window.length:

  A tuning parameter for controlling the level of concavity.

- minLambda:

  Minimum value for density for scaling when fitting inhomogeneous
  L-curves.

- lisaFunc:

  Either "K" or "L" curve.

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
cells$cellID <- paste0("cell_", seq_len(nrow(cells)))

# The curves of the first image
inhom <- inhomLocalK(cells[cells$imageID == "s1", ], Rs = c(10, 20, 50))
```
