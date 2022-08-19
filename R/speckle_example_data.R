#' Generate example data
#'
#' @return a dataframe containing cell-level information for sample, group and
#' clusters
#'
#' @export
#'
#' @examples
#'
#' speckle_example_data()
#'
speckle_example_data <- function(){
    # Make up some data with two groups, two biological replicates in each
    # group and three cell types

    # True cell type proportions for 4 samples
    p_s1 <- c(0.5,0.3,0.2)
    p_s2 <- c(0.6,0.3,0.1)
    p_s3 <- c(0.3,0.4,0.3)
    p_s4 <- c(0.4,0.3,0.3)

    # Total numbers of cells per sample
    numcells <- c(1000,1500,900,1200)

    # Generate cell-level vector for sample info
    biorep <- rep(c("s1","s2","s3","s4"),numcells)
    length(biorep)

    # Numbers of cells for each of the 3 clusters per sample
    n_s1 <- p_s1*numcells[1]
    n_s2 <- p_s2*numcells[2]
    n_s3 <- p_s3*numcells[3]
    n_s4 <- p_s4*numcells[4]

    # Assign cluster labels for 4 samples
    cl_s1 <- rep(c("c0","c1","c2"),n_s1)
    cl_s2 <- rep(c("c0","c1","c2"),n_s2)
    cl_s3 <- rep(c("c0","c1","c2"),n_s3)
    cl_s4 <- rep(c("c0","c1","c2"),n_s4)

    # Generate cell-level vector for cluster info
    clust <- c(cl_s1,cl_s2,cl_s3,cl_s4)
    length(clust)

    # Assume s1 and s2 belong to group 1 and s3 and s4 belong to group 2
    grp <- rep(c("grp1","grp2"),c(sum(numcells[c(1,2)]),sum(numcells[c(3,4)])))

    data.frame(clusters=clust, samples=biorep, group=grp)
}
