
res<- speckle_example_data()

test_that("`speckle_example_data` returns as expected", {
  r_to_json <- function(x) {
    path <- tempfile(fileext = ".json")
    jsonlite::write_json(x, path)
    path
  }

  expect_snapshot_file(r_to_json(res), "speckle_example_data.json")
})
