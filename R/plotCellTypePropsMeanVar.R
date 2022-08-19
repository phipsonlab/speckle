#' Plot cell type proportions versus variances
#' 
#' This function returns a plot of the log10(proportion) versus log10(variance)
#' given a matrix of cell type counts. The rows are the clusters/cell types and
#' the columns are the samples.
#' 
#' The expected variance under a binomial distribution is shown in the solid 
#' line, and the points represent the observed variance for each cell type/row 
#' in the counts table. The blue line shows the empirical Bayes variance 
#' that is used in \code{propeller}.
#' 
#' The mean and variance for each cell type is calculated across all samples.
#' 
#' @param x a matrix or table of counts
#'
#' @return a base R plot
#' 
#' @importFrom limma lmFit
#' @importFrom limma eBayes
#' @importFrom graphics legend lines par title
#' @importFrom stats lowess
#' 
#' @export
#' 
#' @author Belinda Phipson
#'
#' @examples
#' library(limma)
#' # Generate some data
#' # Total number of samples
#' nsamp <- 10
#' # True cell type proportions
#' p <- c(0.05, 0.15, 0.35, 0.45)
#' 
#' # Parameters for beta distribution
#' a <- 40
#' b <- a*(1-p)/p
#' 
#' # Sample total cell counts per sample from negative binomial distribution
#' numcells <- rnbinom(nsamp,size=20,mu=5000)
#' true.p <- matrix(c(rbeta(nsamp,a,b[1]),rbeta(nsamp,a,b[2]),
#'           rbeta(nsamp,a,b[3]),rbeta(nsamp,a,b[4])),byrow=TRUE, ncol=nsamp)
#' 
#' counts <- matrix(NA,ncol=nsamp, nrow=nrow(true.p))
#' rownames(counts) <- paste("c",0:(nrow(true.p)-1), sep="")
#' for(j in 1:length(p)){
#'     counts[j,] <- rbinom(nsamp, size=numcells, prob=true.p[j,])
#' }
#' 
#' plotCellTypePropsMeanVar(counts)
#' 
plotCellTypePropsMeanVar <- function(x){
    x <- as.matrix(x)
    params <- estimateBetaParamsFromCounts(x)
    tot.cells <- colSums(x)
    props <- t(t(x)/tot.cells)

    varp <- params$pi*(1-params$pi)/params$n
    var.props <- apply(props,1,var)
    fitp <- lmFit(props)
    fitp <- eBayes(fitp, robust=TRUE)

    ylimits.min <- min(log10(var.props),log10(varp), log10(fitp$s2.post))
    ylimits.max <- max(log10(var.props),log10(varp), log10(fitp$s2.post))

    # par(mfrow=c(1,1))
    plot(log10(rowMeans(props)), log10(var.props), pch=16, cex=2,
            xlab="log10(proportion)", ylab="log10(variance)",
            ylim=c(ylimits.min,ylimits.max), cex.lab=1.5, cex.axis=1.5)
    lines(lowess(log10(rowMeans(props)),log10(varp)))

    lines(lowess(log10(rowMeans(props)), log10(fitp$s2.post)), lwd=2, col=4)
    legend("bottomright", 
            legend=c("Empirical Bayes variance","Binomial variance"),
            col=c(4,1), lty=1, lwd=2)

    title("Mean-Var: cell type proportions", cex.main=1.5)
}
