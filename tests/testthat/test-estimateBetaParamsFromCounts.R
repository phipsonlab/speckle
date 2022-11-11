
data <- speckle_example_data()
x <- table(data$clusters, data$samples)
res = estimateBetaParamsFromCounts(x)

test_that("estimateBetaParamsFromCounts works", {

  expect_equal(length(res), 6)

  expect_equal(res$n, 1100)

  expect_equal(length(res$alpha),ncol(data))
  expect_equal(length(res$beta),ncol(data))
  expect_equal(length(res$pi),ncol(data))
  expect_equal(length(res$dispersion),ncol(data))
  expect_equal(length(res$var),ncol(data))

})
