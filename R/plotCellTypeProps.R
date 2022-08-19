#' Plot cell type proportions for each sample
#'
#' This is a plotting function that shows the cell type composition for each
#' sample as a stacked barplot. The \code{plotCellTypeProps} returns a
#' \code{ggplot2} object enabling the user to make style changes as required.
#'
#' @param x object of class \code{SingleCellExperiment} or \code{Seurat}
#' @param clusters a factor specifying the cluster or cell type for every cell.
#' For \code{SingleCellExperiment} objects this should correspond to a column
#' called \code{clusters} in the \code{colData} assay. For \code{Seurat}
#' objects this will be extracted by a call to \code{Idents(x)}.
#' @param sample a factor specifying the biological replicate for each cell.
#' For \code{SingleCellExperiment} objects this should correspond to a column
#' called \code{sample} in the \code{colData} assay and for \code{Seurat}
#' objects this should correspond to \code{x$sample}.
#'
#' @return a ggplot2 object
#'
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 theme
#' @export
#'
#' @author Belinda Phipson
#'
#' @examples
#'
#' library(speckle)
#' library(ggplot2)
#' library(limma)
#'
#' # Generate some fake data from a multinomial distribution
#' # Group A, 4 samples, 1000 cells in each sample
#' countsA <- rmultinom(4, size=1000, prob=c(0.1,0.3,0.6))
#' colnames(countsA) <- paste("s",1:4,sep="")
#'
#' # Group B, 3 samples, 800 cells in each sample
#'
#' countsB <- rmultinom(3, size=800, prob=c(0.2,0.05,0.75))
#' colnames(countsB) <- paste("s",5:7,sep="")
#' rownames(countsA) <- rownames(countsB) <- paste("c",0:2,sep="")
#'
#' allcounts <- cbind(countsA, countsB)
#' sample <- c(rep(colnames(allcounts),allcounts[1,]),
#'           rep(colnames(allcounts),allcounts[2,]),
#'           rep(colnames(allcounts),allcounts[3,]))
#' clust <- rep(rownames(allcounts),rowSums(allcounts))
#'
#' plotCellTypeProps(clusters=clust, sample=sample)
#'
plotCellTypeProps <- function(x=NULL, clusters=NULL, sample=NULL)
{
    if(is.null(x) & is.null(sample) & is.null(clusters))
        stop("Please provide either a SingleCellExperiment object or Seurat 
                object with required annotation metadata, or explicitly provide
                clusters and sample information")

    if((is.null(clusters) | is.null(sample)) & !is.null(x)){
        # Extract cluster, sample and group info from SCE object
        if(is(x,"SingleCellExperiment"))
            y <- .extractSCE(x)

        # Extract cluster, sample and group info from Seurat object
        if(is(x,"Seurat"))
            y <- .extractSeurat(x)
            clusters <- y$clusters
            sample <- y$sample
    }

    prop.list <- getTransformedProps(clusters, sample)

    Proportions <- as.vector(t(prop.list$Proportions))
    Samples <- rep(colnames(prop.list$Proportions), nrow(prop.list$Proportions))
    Clusters <- rep(rownames(prop.list$Proportions),
                    each=ncol(prop.list$Proportions))

    plotdf <- data.frame(Samples=Samples, Clusters=Clusters,
                            Proportions=Proportions)

    ggplot(plotdf,aes(x=Samples,y=Proportions,fill=Clusters)) +
        geom_bar(stat="identity") +
        theme(axis.text.x = element_text(size=12),
                axis.text.y = element_text(size=12),
                axis.title = element_text(size=14),
                legend.text = element_text(size=12),
                legend.title = element_text(size=14))

}
