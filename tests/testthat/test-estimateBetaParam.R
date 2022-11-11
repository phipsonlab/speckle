

set.seed(6879)
props <- rbeta(1000, shape1=2, shape2=10)
res = estimateBetaParam(props)

test_that("estimateBetaParam works", {

  expect_equal(length(res), 2)

  expect_equal(round(res$a,2), 2.03)

  expect_equal(round(res$b,2),9.68)
})
