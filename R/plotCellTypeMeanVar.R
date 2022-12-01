#' Plot cell type counts means versus variances
#' 
#' This function returns a plot of the log10(mean) versus log10(variance) of 
#' the cell type counts. The function takes a matrix of cell type counts as 
#' input. The rows are the clusters/cell types and the columns are the samples.
#' 
#' The expected variance under a binomial distribution is shown in the solid 
#' line, and the points represent the observed variance for each cell type/row 
#' in the counts table. The expected variance under different model assumptions
#' are shown in the different coloured lines.
#' 
#' The mean and variance for each cell type is calculated across all samples.
#'
#' @param x a matrix or table of counts
#'
#' @return a base R plot
#' 
#' @importFrom edgeR estimateDisp
#' @importFrom edgeR DGEList
#' @importFrom graphics legend lines par title
#' @importFrom stats lowess
#' 
#' @export
#' 
#' @author Belinda Phipson
#' 
#' @examples
#' library(edgeR)
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
#' plotCellTypeMeanVar(counts)
#' 
plotCellTypeMeanVar <- function(x){
    # Make sure input is a matrix
    x <- as.matrix(x)
    # Normalise the cell type counts 
    nc <- normCounts(x)

    # Beta binomial variance
    params <- estimateBetaParamsFromCounts(x)

    #Observed variance
    means <- rowMeans(nc)
    means.mat <- matrix(rep(means,ncol(nc)),ncol=ncol(nc))
    vars <- rowSums((nc-means.mat)^2)/(ncol(nc)-1)

    # Binomial variance
    ebv <- params$n*params$pi*(1-params$pi)

    # Negative binomial variance
    y <- estimateDisp(DGEList(x))
    tagvars <- means + means^2 * y$tagwise.dispersion

    ylimits.min <- min(log10(vars),log10(ebv), log10(params$var),
                        log10(tagvars), log10(params$n*params$pi))
    ylimits.max <- max(log10(vars),log10(ebv), log10(params$var),
                        log10(tagvars), log10(params$n*params$pi))

    par(mar=c(5,5,2,2))
    #par(mfrow=c(1,1))
    plot(log10(means),log10(vars), pch=16, cex=1.5, cex.lab=1.5, cex.axis=1.5,
            ylim=c(ylimits.min,ylimits.max),
            xlab = "log10(mean)", ylab = "log10(variance)")
    lines(lowess(log10(means), log10(ebv)), col=1, lwd=2)
    lines(lowess(log10(means), log10(params$var)), lwd=2, col=4)
    lines(lowess(log10(means),log10(tagvars)),lwd=2,col="purple", lty=2)
    lines(lowess(log10(means),log10(params$n*params$pi)),lwd=2,col="red", lty=2)
    legend("bottomright", legend=c("Beta-binomial","Negative binomial", 
                                    "Binomial", "Poisson"),
            col=c(4,"purple",1,2), lty=c(1,2,1,2), lwd=2)
    title("Mean-Var: cell type counts", cex.main=1.5)

}
