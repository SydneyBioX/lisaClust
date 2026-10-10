# Introduction to lisaClust

Abstract

Which parts of a tumour are tumour cells, which are immune infiltrate
and which are stroma? Which local neighbourhoods of cells are more
common in one group of patients than another? lisaClust describes, for
every cell, which cell types surround it, and groups cells with similar
surroundings into regions. With a few regions these are tissue domains;
with many they are cellular niches. This vignette finds domains in
triple-negative breast cancers and niches in colorectal cancers, and
compares them between groups of patients.

## Introduction

For each cell, lisaClust (Patrick et al. 2023) counts the cells of every
type within a few radii of it and compares each count with what would be
expected if that type were spread evenly over the image. These local
indicators of spatial association (LISA) describe the cell’s
surroundings: a tumour cell deep in a tumour nest has many tumour
neighbours and few immune neighbours, while a tumour cell at the
invasive margin has both. Clustering the cells by their LISA groups
cells with similar surroundings, whatever their own type, into
**regions**.

The number of regions sets the scale of the question. With a few
regions, they are broad tissue **domains**, such as tumour, stroma and
immune infiltrate. With 20 or more, they are cellular **niches**:
recurring local combinations of cell types, such as T cells mixed with
macrophages at a tumour border.

lisaClust needs only the type and position of every cell, so it applies
to any technology that gives segmented, typed cells: spatial proteomics
such as imaging mass cytometry, MIBI and CODEX, and spatial
transcriptomics such as Xenium, CosMx, MERSCOPE or high-definition
spatial transcriptomics (Patrick et al. 2023). It is not designed for
spot-based data such as Visium, where a spot holds several cells. It
works alongside other SydneyBioX packages: *spicyR* tests whether pairs
of cell types co-localise differently between groups of patients, and
*Statial* tests co-localisation within tissue compartments. Here we use
spicyR’s
[`colTest()`](https://sydneybiox.github.io/spicyR/reference/colTest.html)
to compare regions between groups of patients, and
*tidySpatialExperiment* to summarise the cells of each image with dplyr
verbs.

## Installation

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("lisaClust", "spicyR", "tidySpatialExperiment", "SpatialDatasets"))
```

``` r

library(lisaClust)
library(spicyR)
library(tidySpatialExperiment)
library(SpatialDatasets)
library(ggplot2)
theme_set(theme_classic())
```

## Tissue domains in triple-negative breast cancer

### The data

Keren et al. (2018) imaged tumours from 40 patients with triple-negative
breast cancer by MIBI-TOF, one image per patient, and classified each
tumour by how its immune and tumour cells are arranged:
*compartmentalised*, where immune cells and tumour cells occupy separate
areas; *mixed*, where they are interspersed; and *cold*, with few immune
cells. The data are a `SpatialExperiment` in the SpatialDatasets
package, with the coordinates in pixels.

``` r

kerenSPE <- spe_Keren_2018()
kerenSPE |> distinct(imageID, tumour_type) |> count(tumour_type)
#> # A tibble: 3 × 2
#>   tumour_type           n
#>   <fct>             <int>
#> 1 cold                  6
#> 2 compartmentalised    15
#> 3 mixed                19
sort(table(kerenSPE$cellType), decreasing = TRUE)
#> 
#> Keratin_Tumour    Macrophages     CD8_T_cell     CD4_T_cell         B_cell 
#>          99487          20616          15698          12438           9115 
#>    Mesenchymal   Other_Immune     DC_or_Mono       dn_T_CD3         Tumour 
#>           8170           6891           5049           3848           3167 
#>    Mono_or_Neu    Neutrophils    Endothelial   Unidentified          Tregs 
#>           3110           3018           2086           1725           1341 
#>             DC             NK 
#>           1245            674
```

### Finding domains

[`lisaClust()`](https://github.com/ellispatrick/lisaClust/reference/lisaClust.md)
computes the LISA of every cell and clusters them with k-means into `k`
regions, stored in a new `colData` column (`region` by default). It
looks for columns called `imageID` and `cellType`, and takes the
coordinates of a `SpatialExperiment` from
[`spatialCoords()`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment-methods.html).
Here we ask for four domains, describing each cell’s surroundings within
20, 50 and 100 pixels.

``` r

set.seed(51773)
kerenSPE <- lisaClust(kerenSPE, k = 4, r = c(20, 50, 100))
table(kerenSPE$region)
#> 
#> region_1 region_2 region_3 region_4 
#>    62148   104863    23516     7151
```

k-means starts from random centres, so set a seed to make the regions
reproducible.

### What is in each domain?

[`regionMap()`](https://github.com/ellispatrick/lisaClust/reference/regionMap.md)
shows how much more often each cell type is found in each region than if
cell types were spread evenly over the regions: above one, the type is
enriched in that region.

``` r

regionMap(kerenSPE, type = "bubble")
```

![](lisaClust_files/figure-html/keren-region-map-1.png)

The regions are numbered arbitrarily.
[`nameRegions()`](https://github.com/ellispatrick/lisaClust/reference/nameRegions.md)
names each one by the marker cell type it is most enriched for, here
among tumour cells, CD8 T cells, B cells and mesenchymal cells. Regions
most enriched for the same marker get the same name, so they are merged
into one domain.

``` r

markers <- c(tumour = "Keratin_Tumour", `T cell` = "CD8_T_cell", `B cell` = "B_cell", stroma = "Mesenchymal")
kerenSPE <- nameRegions(kerenSPE, markers, regionName = "domain")
round(100 * prop.table(table(kerenSPE$domain)), 1)
#> 
#> B cell stroma T cell tumour 
#>    3.6   31.4   11.9   53.0
```

The cell types in each domain:

``` r

regionMap(kerenSPE, region = "domain", type = "bubble")
```

![](lisaClust_files/figure-html/keren-domain-map-1.png)

### Domains in two tumours

[`hatchingPlot()`](https://github.com/ellispatrick/lisaClust/reference/hatchingPlot.md)
draws the cells coloured by type, with each region outlined and filled
with its own hatching, so types and regions can be seen together. Here
are the compartmentalised tumour and the mixed tumour with the most
cells in the T cell domain. With tidySpatialExperiment, `group_by()` and
`summarise()` work on the cells of a `SpatialExperiment` and return a
table.

``` r

tcellShare <- kerenSPE |>
  group_by(imageID, tumour_type) |>
  summarise(share = mean(domain == "T cell"), .groups = "drop")
examples <- tcellShare |>
  filter(tumour_type != "cold") |>
  group_by(tumour_type) |>
  slice_max(share)
examples
#> # A tibble: 2 × 3
#> # Groups:   tumour_type [2]
#>   imageID tumour_type       share
#>   <chr>   <fct>             <dbl>
#> 1 16      compartmentalised 0.403
#> 2 29      mixed             0.170
```

``` r

hatchingPlot(kerenSPE, useImages = examples$imageID, region = "domain")
```

![](lisaClust_files/figure-html/keren-hatching-1.png)

### Do the domains differ between kinds of tumour?

Each patient has one image, so the share of each domain in an image is a
value per patient.
[`regionBoxPlot()`](https://github.com/ellispatrick/lisaClust/reference/regionBoxPlot.md)
draws these shares by group.

``` r

regionBoxPlot(kerenSPE, condition = "tumour_type", region = "domain")
```

![](lisaClust_files/figure-html/keren-domain-share-1.png)

spicyR’s
[`colTest()`](https://sydneybiox.github.io/spicyR/reference/colTest.html)
computes the share of each domain in each image and compares the groups
with a Wilcoxon test; here we compare compartmentalised with mixed
tumours.

``` r

domainTests <- kerenSPE |>
  filter(tumour_type != "cold") |>
  colTest(condition = "tumour_type", feature = "domain", type = "wilcox")
domainTests
#>        tval.W    pval adjPval cluster
#> tumour      6 3.2e-08 1.3e-07  tumour
#> T cell    250 9.1e-05 1.8e-04  T cell
#> B cell    240 7.5e-04 1.0e-03  B cell
#> stroma    230 3.1e-03 3.1e-03  stroma
medians <- with(tcellShare, tapply(share, tumour_type, median))
round(100 * medians, 1)
#>              cold compartmentalised             mixed 
#>               1.3              15.2               4.5
pT <- domainTests["T cell", "pval"]
```

The T cell domain makes up a median of 15.2% of the cells in a
compartmentalised tumour and 4.5% in a mixed one (Wilcoxon p = 9.1e-05).
In compartmentalised tumours the immune cells form their own areas,
which become a domain; in mixed tumours they sit among the tumour cells
and fall into the tumour domain.

Keren et al. (2018) classified the tumours from how their immune and
tumour cells are arranged, so this difference is a check that the
domains capture that arrangement, not an independent finding.

## Cellular niches in colorectal cancer

### The data

Schürch et al. (2020) imaged colorectal cancers from 35 patients by
CODEX, four tissue microarray cores per patient. The patients had one of
two patterns of immune response: a Crohn’s-like reaction (CLR, with
lymphoid follicles; 17 patients) or diffuse inflammatory infiltration
(DII; 18 patients). We leave out the objects the authors labelled as
debris (`dirt`).

``` r

schurchSPE <- spe_Schurch_2020() |>
  filter(cellType != "dirt") |>
  mutate(condition = ifelse(group == 1, "CLR", "DII"))
schurchSPE |> distinct(patients, condition) |> count(condition)
#> # A tibble: 2 × 2
#>   condition     n
#>   <chr>     <int>
#> 1 CLR          17
#> 2 DII          18
length(unique(schurchSPE$cellType))
#> [1] 28
```

### Finding niches

With 20 regions, the regions are finer combinations of cell types. We
store them in a column called `niche`.

``` r

set.seed(51773)
schurchSPE <- lisaClust(schurchSPE, k = 20, r = c(20, 50, 100), regionName = "niche")
```

``` r

regionMap(schurchSPE, region = "niche", type = "bubble") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
```

![](lisaClust_files/figure-html/schurch-region-map-1.png)

### Comparing niches between groups of patients

The question is whether a niche makes up more of the tissue in one group
of patients than in the other. The four cores of a patient come from the
same tumour, so they are not independent: we pool the cells of each
patient’s cores, so that each patient contributes one share of each
niche, and compare the two groups of patients. Setting
`imageID = "patients"` in
[`colTest()`](https://sydneybiox.github.io/spicyR/reference/colTest.html)
does this pooling. With 20 niches tested,
[`colTest()`](https://sydneybiox.github.io/spicyR/reference/colTest.html)
also adjusts the p-values for multiple testing (`adjPval`).

``` r

patientTests <- colTest(schurchSPE, condition = "condition", feature = "niche",
                        imageID = "patients", type = "wilcox")
head(patientTests, 5)
#>           tval.W   pval adjPval   cluster
#> region_13    250 0.0016   0.019 region_13
#> region_7     240 0.0019   0.019  region_7
#> region_11     68 0.0043   0.029 region_11
#> region_4      74 0.0083   0.042  region_4
#> region_5     220 0.0200   0.080  region_5
```

For comparison, the same test treating each core as an independent
sample:

``` r

coreTests <- colTest(schurchSPE, condition = "condition", feature = "niche",
                     imageID = "imageID", type = "wilcox")
c(patients = sum(patientTests$adjPval < 0.05), cores = sum(coreTests$adjPval < 0.05))
#> patients    cores 
#>        4        4
```

With the 35 patients as the units, 4 niches differ between CLR and DII
patients at a 5% false discovery rate (smallest adjusted p-value 0.019).
Treating each of the 140 cores as an independent sample instead, 4
niches differ, with a smallest adjusted p-value of 0.00088.

Cores or images from the same patient are alike, so treating them as
independent overstates the evidence and can turn up differences that the
patients do not support. Compare niches with patients as the units, and
adjust for the number of niches tested.

`regionMap(type = "table")` gives the enrichment of each cell type in
each niche as a matrix, so we can list the cell types the niche with the
smallest p-value is most enriched for, and its median share in each
group.

``` r

top <- patientTests$cluster[1]
enrichment <- regionMap(schurchSPE, region = "niche", type = "table")
topTypes <- sort(enrichment[, top], decreasing = TRUE)
round(topTypes[1:5], 1)
#>            b_cells        cd4_t_cells        cd3_t_cells         cd11c_d_cs 
#>                8.3                2.1                2.1                1.4 
#> cd4_t_cells_cd45ro 
#>                1.4

topMedians <- schurchSPE |>
  group_by(patients, condition) |>
  summarise(share = mean(niche == top), .groups = "drop") |>
  with(tapply(share, condition, median))
round(100 * topMedians, 1)
#> CLR DII 
#> 6.4 3.4
```

The niche with the smallest p-value, region_13, is most enriched for
b_cells, cd4_t_cells, cd3_t_cells. It makes up a median of 6.4% of the
cells in CLR patients and 3.4% in DII patients, as expected for the
lymphoid follicles that define a Crohn’s-like reaction:

``` r

regionBoxPlot(schurchSPE, condition = "condition", region = "niche", imageID = "patients", regions = top)
```

![](lisaClust_files/figure-html/schurch-top-niche-1.png)

## Using your own clustering

[`lisa()`](https://github.com/ellispatrick/lisaClust/reference/lisa.md)
returns the LISA themselves, one row per cell and one column per radius
and cell type, so you can cluster them any way you like.

``` r

curves <- lisa(kerenSPE[, kerenSPE$imageID %in% examples$imageID], r = c(20, 50, 100))
dim(curves)
#> [1] 13031    51
curves[1:3, 1:4]
#>        20_Keratin_Tumour 20_CD8_T_cell 20_dn_T_CD3 20_CD4_T_cell
#> cell_1         3.7228902    -0.8831332  -0.2041394    -0.8778542
#> cell_2        -0.7323238     5.5545686  -0.2041394    -0.8778542
#> cell_3        -0.7323238     5.4766730  -0.2041394    -0.8778542

# for example, k-means with several random starts
set.seed(51773)
km <- kmeans(curves, centers = 4, nstart = 10, iter.max = 50)
table(km$cluster)
#> 
#>    1    2    3    4 
#> 5466 1362 6043  160
```

## Choosing the settings

**The number of regions.** A few regions describe the large compartments
of a tissue; many regions separate finer neighbourhoods but split rarer
combinations across several regions and make each one harder to
interpret. Look at
[`regionMap()`](https://github.com/ellispatrick/lisaClust/reference/regionMap.md)
for a few values of `k` and choose the coarsest one that separates the
structures you care about; generic rules for the number of clusters,
such as the elbow or silhouette, can merge tissue structures that are
known to be distinct (Patrick et al. 2023).

**The radii.** Choose them from the scales of interaction you expect, in
the units of the coordinates: a radius of two or three cell diameters
describes the cell’s immediate neighbours, and larger radii describe its
wider surroundings. Several radii (`r`) are used together; run time
grows with the number and size of the radii, so a few well-chosen ones
are better than many.

**Uneven cell density.** By default the expected count assumes cells are
spread evenly over the image. With `sigma`, each neighbour is weighted
by the inverse of the local density of all cells, estimated with a
Gaussian kernel of bandwidth `sigma`, so that a dense area does not by
itself look enriched for every cell type.

## How it works

For cell *i*, cell type *j* and radius *r*, let *n_(ij)(r)* be the
number of cells of type *j* within *r* of cell *i*, counting cell *i*
itself when it is of type *j*. If the cells of type *j* were spread
evenly over the image window *W*, with density
$`\lambda_j = N_j / |W|`$, the expected count would be

``` math
E_{ij}(r) = \lambda_j \, \pi r^2 e_i(r),
```

where $`e_i(r)`$ is the share of the disc of radius *r* around cell *i*
that lies inside the window, which corrects for cells near the edge. The
LISA is the standardised difference

``` math
\frac{n_{ij}(r) - E_{ij}(r)}{\sqrt{E_{ij}(r)}},
```

a local version of Ripley’s K function for the pair (`lisaFunc = "L"`
gives $`\sqrt{n_{ij}(r)} - \sqrt{E_{ij}(r)}`$ instead). The window is
the convex hull of the image’s cells by default (`window`). The
neighbour search and the edge correction run in C++, so the LISA of
hundreds of thousands of cells take seconds.

## Reporting results

A methods sentence might read: “We used lisaClust (version 1.21.7) to
compute, for each cell, local indicators of spatial association with
every cell type at radii of 20, 50 and 100 pixels, and clustered the
cells into 20 niches by k-means. The share of each niche in each patient
was compared between groups with a Wilcoxon test (spicyR’s colTest),
with patients as the units, and p-values were adjusted across niches by
the Benjamini–Hochberg method.” Show
[`regionMap()`](https://github.com/ellispatrick/lisaClust/reference/regionMap.md)
and a
[`hatchingPlot()`](https://github.com/ellispatrick/lisaClust/reference/hatchingPlot.md)
of example images with the results.

``` r

citation("lisaClust")
#> To cite package 'lisaClust' in publications use:
#> 
#>   Patrick E, Canete NP, Iyengar SS, Harman AN, Sutherland GT, Yang P
#>   (2023). "Spatial analysis for highly multiplexed imaging data to
#>   identify tissue microenvironments." _Cytometry Part A_, *103*,
#>   593-599. doi:10.1002/cyto.a.24729
#>   <https://doi.org/10.1002/cyto.a.24729>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {Spatial analysis for highly multiplexed imaging data to identify tissue microenvironments},
#>     author = {Ellis Patrick and Nicolas P. Canete and Sourish S. Iyengar and Andrew N. Harman and Greg T. Sutherland and Pengyi Yang},
#>     journal = {Cytometry Part A},
#>     volume = {103},
#>     pages = {593--599},
#>     year = {2023},
#>     doi = {10.1002/cyto.a.24729},
#>   }
```

## Session information

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] SpatialDatasets_1.10.0          ExperimentHub_3.2.2            
#>  [3] AnnotationHub_4.2.2             BiocFileCache_3.2.0            
#>  [5] dbplyr_2.6.0                    tidySpatialExperiment_1.8.0    
#>  [7] ggplot2_4.0.3                   tidyr_1.3.2                    
#>  [9] dplyr_1.2.1                     tidySingleCellExperiment_1.22.0
#> [11] ttservice_0.5.3                 SpatialExperiment_1.22.0       
#> [13] SingleCellExperiment_1.34.0     SummarizedExperiment_1.42.0    
#> [15] Biobase_2.72.0                  GenomicRanges_1.64.0           
#> [17] Seqinfo_1.2.0                   IRanges_2.46.0                 
#> [19] S4Vectors_0.50.3                BiocGenerics_0.58.1            
#> [21] generics_0.1.4                  MatrixGenerics_1.24.0          
#> [23] matrixStats_1.5.0               spicyR_1.24.0                  
#> [25] lisaClust_1.21.7                BiocStyle_2.40.0               
#> 
#> loaded via a namespace (and not attached):
#>   [1] splines_4.6.1               later_1.4.8                
#>   [3] filelock_1.0.3              tibble_3.3.1               
#>   [5] polyclip_1.10-7             httr2_1.3.0                
#>   [7] lifecycle_1.0.5             Rdpack_2.6.6               
#>   [9] rstatix_1.1.0               lattice_0.22-9             
#>  [11] MASS_7.3-65                 MultiAssayExperiment_1.38.0
#>  [13] backports_1.5.1             magrittr_2.0.5             
#>  [15] plotly_4.12.1               sass_0.4.10                
#>  [17] rmarkdown_2.32              jquerylib_0.1.4            
#>  [19] yaml_2.3.12                 httpuv_1.6.17              
#>  [21] otel_0.2.0                  doRNG_1.8.6.3              
#>  [23] ClassifyR_3.16.0            dcanr_1.28.0               
#>  [25] spatstat.sparse_3.2-0       DBI_1.3.0                  
#>  [27] minqa_1.2.8                 RColorBrewer_1.1-3         
#>  [29] abind_1.4-8                 purrr_1.2.2                
#>  [31] rappdirs_0.3.4              tweenr_2.0.3               
#>  [33] spatstat.utils_3.2-5        pheatmap_1.0.13            
#>  [35] goftest_1.2-3               spatstat.random_3.5-2      
#>  [37] pkgdown_2.2.1               codetools_0.2-20           
#>  [39] DelayedArray_0.38.2         ggforce_0.5.0              
#>  [41] tidyselect_1.2.1            farver_2.1.2               
#>  [43] lme4_2.0-6                  viridis_0.6.5              
#>  [45] spatstat.explore_3.8-3      jsonlite_2.0.0             
#>  [47] ellipsis_0.3.3              Formula_1.2-6              
#>  [49] survival_3.8-6              iterators_1.0.14           
#>  [51] systemfonts_1.3.2           foreach_1.5.2              
#>  [53] tools_4.6.1                 ggnewscale_0.5.2           
#>  [55] ragg_1.5.2                  Rcpp_1.1.2                 
#>  [57] glue_1.8.1                  gridExtra_2.3.1            
#>  [59] SparseArray_1.12.3          xfun_0.61                  
#>  [61] mgcv_1.9-4                  ggthemes_7.0.0             
#>  [63] scam_1.2-22                 withr_3.0.3                
#>  [65] numDeriv_2016.8-1.1         BiocManager_1.30.27        
#>  [67] fastmap_1.2.0               ggh4x_0.3.1                
#>  [69] boot_1.3-32                 fansi_1.0.7                
#>  [71] digest_0.6.39               R6_2.6.1                   
#>  [73] mime_0.13                   textshaping_1.0.5          
#>  [75] tensor_1.5.1                spatstat.data_3.1-9        
#>  [77] RSQLite_3.53.3              utf8_1.2.6                 
#>  [79] data.table_1.18.6.1         httr_1.4.9                 
#>  [81] htmlwidgets_1.6.4           S4Arrays_1.12.1            
#>  [83] pkgconfig_2.0.3             gtable_0.3.6               
#>  [85] blob_1.3.0                  S7_0.2.2                   
#>  [87] XVector_0.52.0              htmltools_0.5.9            
#>  [89] carData_3.0-6               bookdown_0.48              
#>  [91] fftwtools_0.9-11            scales_1.4.0               
#>  [93] ggupset_0.4.1               png_0.1-9                  
#>  [95] spatstat.univar_3.2-0       reformulas_0.4.4           
#>  [97] knitr_1.52                  reshape2_1.4.5             
#>  [99] rjson_0.2.23                curl_8.0.0                 
#> [101] nlme_3.1-169                nloptr_2.2.1               
#> [103] bdsmatrix_1.3-7             cachem_1.1.0               
#> [105] stringr_1.6.0               BiocVersion_3.23.1         
#> [107] parallel_4.6.1              concaveman_1.2.0           
#> [109] AnnotationDbi_1.74.0        desc_1.4.3                 
#> [111] pillar_1.11.1               grid_4.6.1                 
#> [113] vctrs_0.7.3                 coxme_2.2-22               
#> [115] promises_1.5.0              ggpubr_1.0.0               
#> [117] car_3.1-5                   xtable_1.8-8               
#> [119] evaluate_1.0.5              magick_2.9.1               
#> [121] cli_3.6.6                   compiler_4.6.1             
#> [123] crayon_1.5.3                rlang_1.3.0                
#> [125] rngtools_1.5.2              ggsignif_0.6.4             
#> [127] labeling_0.4.3              plyr_1.8.9                 
#> [129] fs_2.1.0                    stringi_1.8.9              
#> [131] viridisLite_0.4.3           deldir_2.0-4               
#> [133] BiocParallel_1.46.0         Biostrings_2.80.2          
#> [135] lmerTest_3.2-1              spatstat.geom_3.8-3        
#> [137] Matrix_1.7-5                bit64_4.8.6                
#> [139] KEGGREST_1.52.2             shiny_1.14.0               
#> [141] rbibutils_2.4.1             tidygate_1.0.19            
#> [143] memoise_2.0.1               igraph_2.3.4               
#> [145] broom_1.0.13                bslib_0.12.0               
#> [147] bit_4.6.0
```

## References

Keren, L, M Bosse, D Marquez, et al. 2018. “A Structured Tumor-Immune
Microenvironment in Triple Negative Breast Cancer Revealed by
Multiplexed Ion Beam Imaging.” *Cell* 174 (6): 1373–1387.e19.

Patrick, Ellis, Nicolas P Canete, Sourish S Iyengar, Andrew N Harman,
Greg T Sutherland, and Pengyi Yang. 2023. “Spatial Analysis for Highly
Multiplexed Imaging Data to Identify Tissue Microenvironments.”
*Cytometry Part A* 103: 593–99. <https://doi.org/10.1002/cyto.a.24729>.

Schürch, Christian M et al. 2020. “Coordinated Cellular Neighborhoods
Orchestrate Antitumoral Immunity at the Colorectal Cancer Invasive
Front.” *Cell* 182 (5): 1341–59.
<https://doi.org/10.1016/j.cell.2020.07.005>.
