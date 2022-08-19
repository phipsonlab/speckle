
# speckle

<!-- badges: start -->
<!-- badges: end -->

The speckle package currently contains functions to analyse differences in cell 
type proportions in single cell RNA-seq data. As our research into specialised 
analyses of single cell data continues we anticipate that the package will be 
updated with new functions.

The propeller, propeller.ttest and propeller.anova functions perform 
statistical tests for differences in cell type composition in single cell data. 
In order to test for differences in cell type proportions between multiple 
experimental conditions at least one of the groups must have some form of 
biological replication (i.e. at least two samples). For a two group scenario, 
the absolute minimum sample size is thus 
three. Since there are many technical aspects which can affect cell type 
proportion estimates, having biological replication is essential for a 
meaningful analysis.

The propeller function takes a SingleCellExperiment or Seurat object as input,
extracts the relevant cell information, and tests whether the cell type 
proportions are statistically significantly different between experimental
conditions/groups. The user can also explicitly pass the cluster, sample and 
experimental group information to propeller. P-values and false discovery rates 
are outputted for each cell type. 

The propeller.ttest() and propeller.anova() are more general functions that can 
account for additional covariates in the analysis.

## Installation

If you would like to view the speckle vignette, you can install the released 
version of speckle from [github](https://github.com/phipsonlab/speckle) using the 
following commands:

``` r
# devtools/remotes won't install Suggested packages from Bioconductor
BiocManager::install(c("CellBench", "BiocStyle", "scater"))

remotes::install_github("phipsonlab/speckle", build_vignettes = TRUE, 
dependencies = "Suggest")
```

In order to view the vignette for speckle use the following command:

``` r
browseVignettes("speckle")
```

If you don't care to view the glorious vignette you can also install speckle as 
follows:

``` r
library(devtools)
devtools::install_github("phipsonlab/speckle")
```

## propeller example

This is a basic example which shows you how to test for differences in cell 
type proportions in a two group experimental set up.

``` r
library(speckle)
library(limma)
library(ggplot2)

# Get some example data which has two groups, three cell types and two 
# biological replicates in each group
cell_info <- speckle_example_data()
head(cell_info)

# Run propeller testing for cell type proportion differences between the two 
# groups
propeller(clusters = cell_info$clusters, sample = cell_info$samples, 
group = cell_info$group)

# Plot cell type proportions
plotCellTypeProps(clusters=cell_info$clusters, sample=cell_info$samples)
```

Please note that this basic implementation is for when you are only modelling
group information. When you have additional covariates that you would like to 
account for, please use the propeller.ttest() and propeller.anova() functions
directly. Please read the vignette for examples on how to model a continuous 
variable, account for additional covariates and include a random effect in the 
analysis. 


