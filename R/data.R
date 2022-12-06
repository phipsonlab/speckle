#' Cell type proportions from single cell PBMC data
#'
#' This dataset is from a paper published in PNAS that looked at differences in 
#' immune functioning between young and old, male and female samples: \
#' Huang Z.  et al.  (2021) Effects of sex and aging on the immune cell 
#' landscape as assessed by single-cell transcriptomic analysis. Proc. Natl. 
#' Acad. Sci. USA, 118, e2023216118. 
#' 
#' The cell type proportions were extracted from the Supplementary Material 
#' "Dataset_S02".
#' The sample information was extracted from the sample names (Y=young, M=male,
#' O=old, F=female). The total number of cells profiled across all samples
#' is 174684, but the number of cells per sample is unknown.
#'
#' @format ## `pbmc_props`
#' A list object with the following components:
#' \describe{
#'   \item{proportions }{A data frame of cell type proportions, where the rows are 
#'   cell types and the columns are the samples. There are 24 rows and 20 
#'   columns}
#'   \item{sample_info }{A data frame with age and sex information for each 
#'   sample}
#'   \item{total_cells }{Numeric, the total number of cells profiled across
#'   all samples in the single cell experiment}
#' }
#' @source <https://www.pnas.org/doi/10.1073/pnas.2023216118>,
#' <https://www.pnas.org/doi/suppl/10.1073/pnas.2023216118/suppl_file/pnas.2023216118.sd02.xlsx>
"pbmc_props"