#' Normalise a counts matrix to the median library size
#'
#' This function takes a \code{DGEList} object or matrix of counts and
#' normalises the counts to the median library size. This puts the normalised
#' counts on a similar scale to the original counts.
#'
#' If the input is a DGEList object, the normalisation factors in
#' \code{norm.factors} are taken into account in the normalisation. The prior
#' counts are added proportionally to the library size
#'
#' @param x a \code{DGEList} object or matrix of counts.
#' @param log logical, indicates whether the output should be on the log2 scale
#' or counts scale. Default is FALSE.
#' @param prior.count The prior count to add if the data is log2 normalised.
#' Default is a small count of 0.5.
#' @param lib.size a vector of library sizes to be used during the normalisation
#' step. Default is NULL and will be computed from the counts matrix.
#'
#' @return a matrix of normalised counts
#' 
#' @export normCounts
#' @importFrom stats median
#' @author Belinda Phipson
#'
#' @examples
#' # Simulate some data from a negative binomial distribution with mean equal
#' # to 100 and dispersion set to 1. Simulate 1000 genes and 6 samples.
#' y <- matrix(rnbinom(6000, mu = 100, size = 1), ncol = 6)
#'
#' # Normalise the counts
#' norm.y <- normCounts(y)
#'
#' # Return log2 normalised counts
#' lnorm.y <- normCounts(y, log=TRUE)
#'
#' # Return log2 normalised counts with prior.count = 2
#' lnorm.y2 <- normCounts(y, log=TRUE, prior.count=2)
#'
#' par(mfrow=c(1,2))
#' boxplot(norm.y, main="Normalised counts")
#' boxplot(lnorm.y, main="Log2-normalised counts")
#'
normCounts <-function(x, log=FALSE, prior.count=0.5, lib.size=NULL)
    # Function to normalise to median library size instead of counts per million
    # Input is DGEList object or matrix
    # Belinda Phipson
    # 30 November 2015
{
    if(is(x, "DGEList")){
        lib.size <- x$samples$lib.size*x$samples$norm.factors
        counts <- x$counts
    }
    else{
        counts <- as.matrix(x)
        if(is.null(lib.size)){
            lib.size <- colSums(counts)
        }
        else{
            if(length(lib.size)==ncol(x))
                lib.size <- as.vector(lib.size)
            else{
                message("Vector of library sizes does not match dimensions 
                            of input data. Calculating library sizes 
                            from the counts matrix.")
                lib.size <- colSums(counts)
            }
        }
    }

    M <- median(lib.size)
    if(log){
        prior.count.scaled <- lib.size/mean(lib.size)*prior.count
        lib.size <- lib.size + 2*prior.count.scaled
        log2(t((t(counts)+prior.count.scaled)/lib.size*M))
    }
    else t(t(counts)/lib.size*M)
}
