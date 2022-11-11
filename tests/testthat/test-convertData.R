
# Make up some data with two groups, two biological replicates in each
# group and three cell types
# True cell type proportions for 4 samples
props <- matrix(c(0.5,0.3,0.2,0.6,0.3,0.1,0.3,0.4,0.3,0.4,0.3,0.3),
                ncol=4, nrow=3, byrow=FALSE)
rownames(props) <- c("C0","C1","C2")
colnames(props) <- paste("S",c(1,2,3,4),sep="")
# Total numbers of cells per sample
numcells <- c(1000,1500,900,1200)

# Get data into list object to use as input to propeller.ttest
propslist <- convertDataToList(props, data.type="proportions",
                               transform="asin",
                               scale.fac=numcells)

test_that("`convertDataToList` returns as expected", {
 # local_edition(3)

  r_to_json <- function(x) {
    path <- tempfile(fileext = ".json")
    jsonlite::write_json(x, path)
    path
  }

  expect_snapshot_file(r_to_json(propslist), "converted_data.json")
})

