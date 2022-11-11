test_that("`ggplotColors()` returns as expected", {
  expect_equal( ggplotColors(6),
    c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3")
  )
})

test_that("`ggplotColors()` returns expected number of colors", {
  expect_equal(length(ggplotColors(10)),10)
})
