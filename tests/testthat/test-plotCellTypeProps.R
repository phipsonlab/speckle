

# Generate some fake data from a multinomial distribution
# Group A, 4 samples, 1000 cells in each sample
set.seed(8596)
countsA <- rmultinom(4, size=1000, prob=c(0.1,0.3,0.6))
colnames(countsA) <- paste("s",1:4,sep="")

# Group B, 3 samples, 800 cells in each sample
set.seed(1658)
countsB <- rmultinom(3, size=800, prob=c(0.2,0.05,0.75))
colnames(countsB) <- paste("s",5:7,sep="")
rownames(countsA) <- rownames(countsB) <- paste("c",0:2,sep="")
allcounts <- cbind(countsA, countsB)
sample <- c(rep(colnames(allcounts),allcounts[1,]),
           rep(colnames(allcounts),allcounts[2,]),
           rep(colnames(allcounts),allcounts[3,]))
clust <- rep(rownames(allcounts),rowSums(allcounts))


test_that("`plotCellTypeProps()` plots as expected", {
  #local_edition(3)
  vdiffr::expect_doppelganger(
    title = "plotCellTypeProps",
    fig = plotCellTypeProps(clusters=clust, sample=sample)
  )
})
