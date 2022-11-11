
set.seed(3485)
y <- matrix(rnbinom(6000, mu = 100, size = 1), ncol = 6)
lnorm.y2 <- normCounts(y, log=TRUE, prior.count=2)

test_that("`normCounts` returns as expected", {
  r_to_json <- function(x) {
    path <- tempfile(fileext = ".json")
    jsonlite::write_json(x, path)
    path
  }

  expect_snapshot_file(r_to_json(lnorm.y2), "lnorm_y2.json")
})
