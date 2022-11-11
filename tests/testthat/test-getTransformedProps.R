# True cell type proportions for 4 samples
p_s1 <- c(0.5,0.3,0.2)
p_s2 <- c(0.6,0.3,0.1)
p_s3 <- c(0.3,0.4,0.3)
p_s4 <- c(0.4,0.3,0.3)

# Total numbers of cells per sample
numcells <- c(1000,1500,900,1200)
# Generate cell-level vector for sample info
biorep <- rep(c("s1","s2","s3","s4"),numcells)
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

res = as.data.frame(getTransformedProps(clusters = clust, sample = biorep))


test_that("`getTransformedProps` returns as expected", {

  r_to_json <- function(x) {
    path <- tempfile(fileext = ".json")
    jsonlite::write_json(x, path)
    path
  }

  expect_snapshot_file(r_to_json(res), "transformed_props.json")
})

