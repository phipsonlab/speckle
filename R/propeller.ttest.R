#' Performs t-tests of transformed cell type proportions
#'
#' This function is called by \code{propeller} and performs t-tests between two
#' experimental groups or conditions on the transformed cell type proportions.
#'
#' In order to run this function, the user needs to run the
#' \code{getTransformedProps} function first. The output from
#' \code{getTransformedProps} is used as input. The \code{propeller.ttest}
#' function expects that the design matrix is not in the intercept format
#' and a contrast vector needs to be supplied. This contrast vector will
#' identify the two groups to be tested. Note that additional confounding
#' covariates can be accounted for as extra columns in the design matrix after
#' the group-specific columns.
#'
#' The \code{propeller.ttest} function uses the \code{limma} functions
#' \code{lmFit}, \code{contrasts.fit} and \code{eBayes} which has the additional
#' advantage that empirical Bayes shrinkage of the variances are performed.
#'
#' @param prop.list a list object containing two matrices:
#' \code{TransformedProps} and \code{Proportions}
#' @param design a design matrix with rows corresponding to samples and columns
#' to coefficients to be estimated
#' @param contrasts a vector specifying which columns of the design matrix
#' correspond to the two groups to test
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
#' @importFrom limma contrasts.fit
#' @export propeller.ttest
#'
#' @author Belinda Phipson
#'
#' @seealso \code{\link{propeller}}, \code{\link{getTransformedProps}},
#' \code{\link{lmFit}}, \code{\link{contrasts.fit}}, \code{\link{eBayes}}
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
#'
#'   # Total numbers of cells per sample
#'   numcells <- c(1000,1500,900,1200)
#'
#'   # Generate cell-level vector for sample info
#'   biorep <- rep(c("s1","s2","s3","s4"),numcells)
#'   length(biorep)
#'
#'   # Numbers of cells for each of 3 clusters per sample
#'   n_s1 <- p_s1*numcells[1]
#'   n_s2 <- p_s2*numcells[2]
#'   n_s3 <- p_s3*numcells[3]
#'   n_s4 <- p_s4*numcells[4]
#'
#'   cl_s1 <- rep(c("c0","c1","c2"),n_s1)
#'   cl_s2 <- rep(c("c0","c1","c2"),n_s2)
#'   cl_s3 <- rep(c("c0","c1","c2"),n_s3)
#'   cl_s4 <- rep(c("c0","c1","c2"),n_s4)
#'
#'   # Generate cell-level vector for cluster info
#'   clust <- c(cl_s1,cl_s2,cl_s3,cl_s4)
#'   length(clust)
#'
#'   prop.list <- getTransformedProps(clusters = clust, sample = biorep)
#'
#'   # Assume s1 and s2 belong to group 1 and s3 and s4 belong to group 2
#'   grp <- rep(c("A","B"), each=2)
#'
#'   design <- model.matrix(~0+grp)
#'   design
#'
#'   # Compare Grp A to B
#'   contrasts <- c(1,-1)
#'
#'   propeller.ttest(prop.list, design=design, contrasts=contrasts, robust=TRUE,
#'   trend=FALSE, sort=TRUE)
#'
#'   # Pretend additional sex variable exists and we want to control for it
#'   # in the linear model
#'   sex <- rep(c(0,1),2)
#'   des.sex <- model.matrix(~0+grp+sex)
#'   des.sex
#'
#'   # Compare Grp A to B
#'   cont.sex <- c(1,-1,0)
#'
#'   propeller.ttest(prop.list, design=des.sex, contrasts=cont.sex, robust=TRUE,
#'   trend=FALSE, sort=TRUE)
#'
#'
propeller.ttest <- function(prop.list=prop.list, design=design,
                            contrasts=contrasts, robust=robust, trend=trend,
                            sort=sort)
{
    prop.trans <- prop.list$TransformedProps
    prop <- prop.list$Proportions
    
    # Add check for fewer than 3 cell types
    # Robust eBayes doesn't work with fewer than 3 cell types
    if(nrow(prop.trans)<=2){
        message("Setting robust to FALSE for eBayes for less than 3 cell types")
        robust <- FALSE
    }
    
    fit <- lmFit(prop.trans, design)
    fit.cont <- contrasts.fit(fit, contrasts=contrasts)
    fit.cont <- eBayes(fit.cont, robust=robust, trend=trend)

    # Get mean cell type proportions and relative risk for output
    # If no confounding variable included in design matrix
    if(length(contrasts)==2){
        fit.prop <- lmFit(prop, design)
        z <- apply(fit.prop$coefficients, 1, function(x) x^contrasts)
        RR <- apply(z, 2, prod)
    }
    # If confounding variables included in design matrix exclude them
    else{
        new.des <- design[,contrasts!=0]
        fit.prop <- lmFit(prop,new.des)
        new.cont <- contrasts[contrasts!=0]
        z <- apply(fit.prop$coefficients, 1, function(x) x^new.cont)
        RR <- apply(z, 2, prod)
    }

    fdr <- p.adjust(fit.cont$p.value[,1], method="BH")

    out <- data.frame(PropMean=fit.prop$coefficients, PropRatio=RR,
                    Tstatistic=fit.cont$t[,1], P.Value=fit.cont$p.value[,1],
                    FDR=fdr)
    if(sort){
        o <- order(out$P.Value)
        out[o,]
    }
    else out
}
