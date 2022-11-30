#' Get PBMC proportions data
#'
#' @return a list object with cell type proportions, sample information and
#' total number of cells
#' 
#' @export
#'
#' @examples
#' 
#' pbmc_list <- load_pbmc_props()
#' names(pbmc_list)
#' pbmc_props <- pbmc_list$proportions
#' pbmc_sample_info <- pbmc_list$sample_info
#' tot_cells <- pbmc_list$total_cells
#' 
#' pbmc_props
#' 
#' barplot(as.matrix(pbmc_props),col=ggplotColors(nrow(pbmc_props)),
#' ylab="Cell type proportions", xlab="Samples")
#' 
load_pbmc_props <- function(){
    # Get the cell type proportions from internal data
    props <- pbmc_props
    # Get the sample annotation from the column names
    age <- rep(NA, 20)
    sex <- rep(NA, 20)
    age[grep("Y", colnames(props))] <- "young"
    age[grep("O", colnames(props))] <- "old"
    sex[grep("F", colnames(props))] <- "female"
    sex[grep("M", colnames(props))] <- "male"
    sample_info <- data.frame(age, sex)
    tot_cells <- 174684
    return(list(proportions=props, sample_info=sample_info, 
                total_cells=tot_cells))
}