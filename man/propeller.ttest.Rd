% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/propeller.ttest.R
\name{propeller.ttest}
\alias{propeller.ttest}
\title{Performs t-tests of transformed cell type proportions}
\usage{
propeller.ttest(
    prop.list = prop.list,
    design = design,
    contrasts = contrasts,
    robust = robust,
    trend = trend,
    sort = sort
)
}
\arguments{
\item{prop.list}{a list object containing two matrices:
\code{TransformedProps} and \code{Proportions}}

\item{design}{a design matrix with rows corresponding to samples and columns
to coefficients to be estimated}

\item{contrasts}{a vector specifying which columns of the design matrix
correspond to the two groups to test}

\item{robust}{logical, should robust variance estimation be used. Defaults to
TRUE.}

\item{trend}{logical, should a trend between means and variances be accounted
for. Defaults to FALSE.}

\item{sort}{logical, should the output be sorted by P-value.}
}
\value{
produces a dataframe of results
}
\description{
This function is called by \code{propeller} and performs t-tests between two
experimental groups or conditions on the transformed cell type proportions.
}
\details{
In order to run this function, the user needs to run the
\code{getTransformedProps} function first. The output from
\code{getTransformedProps} is used as input. The \code{propeller.ttest}
function expects that the design matrix is not in the intercept format
and a contrast vector needs to be supplied. This contrast vector will
identify the two groups to be tested. Note that additional confounding
covariates can be accounted for as extra columns in the design matrix after
the group-specific columns.

The \code{propeller.ttest} function uses the \code{limma} functions
\code{lmFit}, \code{contrasts.fit} and \code{eBayes} which has the additional
advantage that empirical Bayes shrinkage of the variances are performed.
}
\examples{
    library(speckle)
    library(ggplot2)
    library(limma)

    # Make up some data

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

    # Numbers of cells for each of 3 clusters per sample
    n_s1 <- p_s1*numcells[1]
    n_s2 <- p_s2*numcells[2]
    n_s3 <- p_s3*numcells[3]
    n_s4 <- p_s4*numcells[4]

    cl_s1 <- rep(c("c0","c1","c2"),n_s1)
    cl_s2 <- rep(c("c0","c1","c2"),n_s2)
    cl_s3 <- rep(c("c0","c1","c2"),n_s3)
    cl_s4 <- rep(c("c0","c1","c2"),n_s4)

    # Generate cell-level vector for cluster info
    clust <- c(cl_s1,cl_s2,cl_s3,cl_s4)
    length(clust)

    prop.list <- getTransformedProps(clusters = clust, sample = biorep)

    # Assume s1 and s2 belong to group 1 and s3 and s4 belong to group 2
    grp <- rep(c("A","B"), each=2)

    design <- model.matrix(~0+grp)
    design

    # Compare Grp A to B
    contrasts <- c(1,-1)

    propeller.ttest(prop.list, design=design, contrasts=contrasts, robust=TRUE,
    trend=FALSE, sort=TRUE)

    # Pretend additional sex variable exists and we want to control for it
    # in the linear model
    sex <- rep(c(0,1),2)
    des.sex <- model.matrix(~0+grp+sex)
    des.sex

    # Compare Grp A to B
    cont.sex <- c(1,-1,0)

    propeller.ttest(prop.list, design=des.sex, contrasts=cont.sex, robust=TRUE,
    trend=FALSE, sort=TRUE)


}
\seealso{
\code{\link{propeller}}, \code{\link{getTransformedProps}},
\code{\link{lmFit}}, \code{\link{contrasts.fit}}, \code{\link{eBayes}}
}
\author{
Belinda Phipson
}
