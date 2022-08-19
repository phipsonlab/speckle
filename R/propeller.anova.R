#' Performs F-tests for transformed cell type proportions
#'
#' This function is called by \code{propeller} and performs F-tests between
#' multiple experimental groups or conditions (> 2) on transformed cell type 
#' proportions.
#'
#' In order to run this function, the user needs to run the
#' \code{getTransformedProps} function first. The output from
#' \code{getTransformedProps} is used as input. The \code{propeller.anova}
#' function expects that the design matrix is not in the intercept format.
#' This \code{coef} vector will identify the columns in the design matrix that
#' correspond to the groups being tested.
#' Note that additional confounding covariates can be accounted for as extra
#' columns in the design matrix, but need to come after the group-specific
#' columns.
#'
#' The \code{propeller.anova} function uses the \code{limma} functions
#' \code{lmFit} and \code{eBayes} to extract F statistics and p-values.
#' This has the additional advantage that empirical Bayes shrinkage of the
#' variances are performed.
#'
#' @param prop.list a list object containing two matrices:
#' \code{TransformedProps} and \code{Proportions}
#' @param design a design matrix with rows corresponding to samples and columns
#' to coefficients to be estimated
#' @param coef a vector specifying which the columns of the design matrix
#' corresponding to the groups to test
#' @param robust logical, should robust variance estimation be used. Defaults to
#' TRUE.
#' @param trend logical, should a trend between means and variances be accounted
#' for. Defaults to FALSE.
#' @param sort logical, should the output be sorted by P-value.
#'
#' @return produces a dataframe of results
#'
#' @importFrom stats p.adjust
#' @importFrom limma lmFit
#' @importFrom limma eBayes
#' @export propeller.anova
#'
#' @author Belinda Phipson
#'
#'@seealso \code{\link{propeller}}, \code{\link{getTransformedProps}},
#' \code{\link{lmFit}}, \code{\link{eBayes}}
#'
#' @examples
#'   library(speckle)
#'   library(ggplot2)
#'   library(limma)
#'
#'   # Make up some data
#'
#'   # True cell type proportions for 4 samples
#'   p_s1 <- c(0.5,0.3,0.2)
#'   p_s2 <- c(0.6,0.3,0.1)
#'   p_s3 <- c(0.3,0.4,0.3)
#'   p_s4 <- c(0.4,0.3,0.3)
#'   p_s5 <- c(0.8,0.1,0.1)
#'   p_s6 <- c(0.75,0.2,0.05)
#'
#'   # Total numbers of cells per sample
#'   numcells <- c(1000,1500,900,1200,1000,800)
#'
#'   # Generate cell-level vector for sample info
#'   biorep <- rep(c("s1","s2","s3","s4","s5","s6"),numcells)
#'   length(biorep)
#'
#'   # Numbers of cells for each of 3 clusters per sample
#'   n_s1 <- p_s1*numcells[1]
#'   n_s2 <- p_s2*numcells[2]
#'   n_s3 <- p_s3*numcells[3]
#'   n_s4 <- p_s4*numcells[4]
#'   n_s5 <- p_s5*numcells[5]
#'   n_s6 <- p_s6*numcells[6]
#'
#'   cl_s1 <- rep(c("c0","c1","c2"),n_s1)
#'   cl_s2 <- rep(c("c0","c1","c2"),n_s2)
#'   cl_s3 <- rep(c("c0","c1","c2"),n_s3)
#'   cl_s4 <- rep(c("c0","c1","c2"),n_s4)
#'   cl_s5 <- rep(c("c0","c1","c2"),n_s5)
#'   cl_s6 <- rep(c("c0","c1","c2"),n_s6)
#'
#'   # Generate cell-level vector for cluster info
#'   clust <- c(cl_s1,cl_s2,cl_s3,cl_s4,cl_s5,cl_s6)
#'   length(clust)
#'
#'   prop.list <- getTransformedProps(clusters = clust, sample = biorep)
#'
#'   # Assume s1 and s2 belong to group A, s3 and s4 belong to group B, s5 and
#'   # s6 belong to group C
#'   grp <- rep(c("A","B","C"), each=2)
#'
#'   # Make sure design matrix does not have an intercept term
#'   design <- model.matrix(~0+grp)
#'   design
#'
#'   propeller.anova(prop.list, design=design, coef=c(1,2,3), robust=TRUE,
#'   trend=FALSE, sort=TRUE)
#'
propeller.anova <- function(prop.list=prop.list, design=design, coef = coef,
                            robust=robust, trend=trend, sort=sort)
{
    prop.trans <- prop.list$TransformedProps
    prop <- prop.list$Proportions
    
    # Add check for fewer than 3 cell types
    # Robust eBayes doesn't work with fewer than 3 cell types
    if(nrow(prop.trans)<=2){
        message("Robust eBayes doesn't work with fewer than 3 cell types.
                    Setting robust to FALSE")
        robust <- FALSE
    }

    # get cell type mean proportions ignoring other variables
    # this assumes that the design matrix is not in Intercept format
    fit.prop <- lmFit(prop, design[,coef])

    # Change design matrix to intercept format
    design[,1] <- 1
    colnames(design)[1] <- "Int"

    # Fit linear model taking into account all confounding variables
    fit <- lmFit(prop.trans,design)

    # Get F statistics corresponding to group information only
    # You have to remove the intercept term for this to work
    fit <- eBayes(fit[,coef[-1]], robust=robust, trend=trend)

    # Extract F p-value
    p.value <- fit$F.p.value
    # and perform FDR adjustment
    fdr <- p.adjust(fit$F.p.value, method="BH")

    out <- data.frame(PropMean=fit.prop$coefficients, Fstatistic= fit$F,
                    P.Value=p.value, FDR=fdr)
    if(sort){
        o <- order(out$P.Value)
        out[o,]
    }
    else out
}
