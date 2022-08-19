#' Estimate parameters of a Beta distribution from counts
#' 
#' This function estimates the two parameters of the Beta distribution, alpha
#' and beta for each cell type. The input is a matrix of cell type counts, 
#' where the rows are the cell types/clusters and the columns are the samples.
#' 
#' This function is called from the plotting function \code{plotCellTypeMeanVar}
#' in order to estimate the variance for the Beta-Binomial distribution for 
#' each cell type.
#' 
#' @param x a matrix of counts
#'
#' @return outputs a list object with the following components
#' \item{n }{Normalised library size}
#' \item{alpha }{a vector of alpha parameters for the Beta distribution for 
#' each cell type}
#' \item{beta }{vector of beta parameters for the Beta distribution for 
#' each cell type}
#' \item{pi }{Estimate of the true proportion for each cell type}
#' \item{dispersion }{Dispersion estimates for each cell type}
#' \item{var }{Variance estimates for each cell type}
#' 
#' @export
#' 
#' @author Belinda Phipson
#'
#' @examples
#' data <- speckle_example_data()
#' x <- table(data$clusters, data$samples)
#' estimateBetaParamsFromCounts(x)
#' 
#' 
estimateBetaParamsFromCounts <- function(x){
    # Make sure input is a matrix
    counts <- as.matrix(x)
    # Normalise the counts so that the total number of counts 
    # per sample is equal
    nc <- normCounts(counts)
    # Get cell type means
    m1 <- rowMeans(nc)
    # Get variance estimate for each cell type
    m2 <- rowSums(nc^2)/ncol(nc)
    n <- mean(colSums(nc))
    alpha <- (n*m1-m2)/(n*(m2/m1-m1-1)+m1)
    beta <- ((n-m1)*(n-m2/m1))/(n*(m2/m1-m1-1)+m1)
    disp <- 1/(alpha+beta)
    pi <- alpha/(alpha+beta)
    var <- n*pi*(1-pi)*(n*disp+1)/(1+disp)
    output <- list(n=n, alpha=alpha, beta=beta, pi=pi, dispersion=disp, var=var)
    output
}
