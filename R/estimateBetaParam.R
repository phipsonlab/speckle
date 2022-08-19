#' Estimate the parameters of a Beta distribution
#'
#' This function estimates the two parameters of the Beta distribution, alpha
#' and beta, given a vector of proportions. It uses the method of moments to
#' do this.
#'
#' @param x a vector of proportions.
#'
#' @return a list object with the estimate of alpha in \code{a} and beta in
#' \code{b}.
#'
#' @export estimateBetaParam
#' @importFrom stats var
#'
#' @author Belinda Phipson
#'
#' @examples
#' # Generate proportions from a beta distribution
#' props <- rbeta(1000, shape1=2, shape2=10)
#' estimateBetaParam(props)
#'
estimateBetaParam <- function(x){
    # solve for the hyperparameters of the beta distribution given a vector
    # of proportions
    mu <- mean(x)
    V <- var(x)
    a <- ((1-mu)/V - 1/mu)*mu^2
    b <- ((1-mu)/V - 1/mu)*mu*(1-mu)
    list(a=a,b=b)
}
