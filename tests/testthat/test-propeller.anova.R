# True cell type proportions for 4 samples
p_s1 <- c(0.5,0.3,0.2)
p_s2 <- c(0.6,0.3,0.1)
p_s3 <- c(0.3,0.4,0.3)
p_s4 <- c(0.4,0.3,0.3)
p_s5 <- c(0.8,0.1,0.1)
p_s6 <- c(0.75,0.2,0.05)

# Total numbers of cells per sample
numcells <- c(1000,1500,900,1200,1000,800)

# Generate cell-level vector for sample info
biorep <- rep(c("s1","s2","s3","s4","s5","s6"),numcells)

# Numbers of cells for each of 3 clusters per sample
n_s1 <- p_s1*numcells[1]
n_s2 <- p_s2*numcells[2]
n_s3 <- p_s3*numcells[3]
n_s4 <- p_s4*numcells[4]
n_s5 <- p_s5*numcells[5]
n_s6 <- p_s6*numcells[6]

cl_s1 <- rep(c("c0","c1","c2"),n_s1)
cl_s2 <- rep(c("c0","c1","c2"),n_s2)
cl_s3 <- rep(c("c0","c1","c2"),n_s3)
cl_s4 <- rep(c("c0","c1","c2"),n_s4)
cl_s5 <- rep(c("c0","c1","c2"),n_s5)
cl_s6 <- rep(c("c0","c1","c2"),n_s6)

# Generate cell-level vector for cluster info
clust <- c(cl_s1,cl_s2,cl_s3,cl_s4,cl_s5,cl_s6)

prop.list <- getTransformedProps(clusters = clust, sample = biorep)

# Assume s1 and s2 belong to group A, s3 and s4 belong to group B, s5 and
#s6 belong to group C
grp <- rep(c("A","B","C"), each=2)

# Make sure design matrix does not have an intercept term
design <- model.matrix(~0+grp)

res = propeller.anova(prop.list, design=design, coef=c(1,2,3), robust=TRUE,
trend=FALSE, sort=TRUE)

test_that("`propeller.anova` returns as expected", {

  r_to_json <- function(x) {
    path <- tempfile(fileext = ".json")
    jsonlite::write_json(x, path)
    path
  }
  expect_equal(nrow(res), 3)
  expect_equal(ncol(res), 6)
  expect_snapshot_file(r_to_json(res), "propeller_anova_table.json")
})

