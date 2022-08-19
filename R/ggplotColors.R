#' Output a vector of colours based on the ggplot colour scheme
#'
#' This function takes as input the number of colours the user would like, and
#' outputs a vector of colours in the ggplot colour scheme.
#'
#' @param g the number of colours to be generated.
#'
#' @return a vector with the names of the colours.
#' @export ggplotColors
#' 
#' @importFrom grDevices hcl
#'
#' @author Belinda Phipson
#'
#' @examples
#' # Generate a palette of 6 colours
#' cols <- ggplotColors(6)
#' cols
#'
#' # Generate some count data
#' y <- matrix(rnbinom(600, mu=100, size=1), ncol=6)
#'
#' par(mfrow=c(1,1))
#' boxplot(y, col=cols)
#'
ggplotColors <- function(g){

    d <- 360/g

    h <- cumsum(c(15, rep(d,g - 1)))

    hcl(h = h, c = 100, l = 65)

}
