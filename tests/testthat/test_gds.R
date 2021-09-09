library(ggmanh)

test_that("Check that the default gds file has been saved.", {
  expect_vector(default_gds_path())
})
