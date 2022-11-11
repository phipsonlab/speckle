
# Generate some data
nsamp <- 10
# True cell type proportions
p <- c(0.05, 0.15, 0.35, 0.45)

# Parameters for beta distribution
a <- 40
b <- a*(1-p)/p

set.seed(986)
# Sample total cell counts per sample from negative binomial distribution
numcells <- rnbinom(nsamp,size=20,mu=5000)
true.p <- matrix(c(rbeta(nsamp,a,b[1]),rbeta(nsamp,a,b[2]),
rbeta(nsamp,a,b[3]),rbeta(nsamp,a,b[4])),byrow=TRUE, ncol=nsamp)
counts <- matrix(NA,ncol=nsamp, nrow=nrow(true.p))
rownames(counts) <- paste("c",0:(nrow(true.p)-1), sep="")
for(j in 1:length(p)){
  counts[j,] <- rbinom(nsamp, size=numcells, prob=true.p[j,])
}

test_that("`plotCellTypeMeanVar()` plots as expected", {
  #local_edition(3)
  vdiffr::expect_doppelganger(
    title = "plotCellTypeMeanVar",
    fig = plotCellTypeMeanVar(counts)
  )
})
