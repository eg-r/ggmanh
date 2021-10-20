library(ggmanh)

test_that("Check that good arguments are provided", {
  df <- data.frame("A" = 1:3, "B" = c('a', 'b','c'), "C" = c(0.05, 0.005, 0.0005), "D" = c("A", "A", "B"))
  expect_error(manhattan_data_preprocess(df))
  expect_error(manhattan_data_preprocess(df, pval.colname = "C", chr.colname = "B", pos.colname = "B"))
  expect_error(manhattan_data_preprocess(
    df, pval.colname = "C", chr.colname = "B", pos.colname = "A",
    highlight.colname = "D", highlight.col = c("A" = "grey")
    )
  )

  tmpdf <- df
  tmpdf[2,"B"] <- NA
  expect_warning(manhattan_data_preprocess(tmpdf, pval.colname = "C", pos.colname = "A", chr.colname = "B"))

  tmpdf <- df
  tmpdf[3,"C"] <- NA
  expect_warning(manhattan_data_preprocess(tmpdf, pval.colname = "C", pos.colname = "A", chr.colname = "B"))

  tmpdf <- df
  tmpdf[1,"A"] <- NA
  expect_warning(manhattan_data_preprocess(tmpdf, pval.colname = "C", pos.colname = "A", chr.colname = "B"))

  tmpdf <- df
  tmpdf[1,"A"] <- NA; tmpdf[2,"B"] <- NA; tmpdf[3,"C"] <- NA
  expect_warning(expect_error(manhattan_data_preprocess(tmpdf, pval.colname = "C", pos.colname = "A", chr.colname = "B")))
})

test_that("Check that preprocess works as intended", {
  df <- data.frame("pos" = c(1,3,4,5,2,4), "chr" = c('a','a','a','b','b','c'), "pval" = c(0.05,0.05,0.0005,0.005,0.000005,0.0005))
  mpdat1 <- manhattan_data_preprocess(
    x = df, pval.colname = "pval", chr.colname = "chr", pos.colname = "pos",
    chr.order = c('a','b','c'), preserve.position = FALSE, thin = FALSE
  )
  mpdat2 <- manhattan_data_preprocess(
    x = df, pval.colname = "pval", chr.colname = "chr", pos.colname = "pos",
    chr.order = c('a','b','c'), preserve.position = TRUE, thin = FALSE
  )
  lg <- 0.15 / 26 * 3 # gap between chromosomes - hard coded in manhattan_preprocess function
  expect_equal(mpdat1$data$new_pos, c(0, 1/2, 1, 1 + lg, 2 + lg, 2 + lg*2 + 1/2))
  expect_equal(mpdat2$data$new_pos, c(0, 2/3 * 1.5, 1.5, 1.5 + lg, 2.5 + lg, 2.5 + lg*2))
  expect_equal(mpdat1$data$pval, mpdat2$data$pval)
  expect_equal(mpdat1$data$pval, c(0.05,0.05,0.0005,0.000005,0.005,0.0005))
})

test_that("Test that the function generates a warning when too many labels are supplied.", {
  df <- data.frame(
    pval = runif(100),
    chr = sample(1:12,100,replace=TRUE),
    pos = sample(100),
    label = as.character(sample(100))
  )
  expect_warning(p <- manhattan_plot(x = df, label.colname = "label"))
})

test_that("Test that the thinPoint function subsets correctly.", {
  dat <- data.frame(
    A1 = c(1:20, 20, 20),
    A2 = c(rep(1, 12), rep(1,5), rep(20, 3), 20, 20) ,
    B = rep(c("a", "b", "c", "d"), times = c(5, 7, 8, 2))
  )
  expect_equal(nrow(thinPoints(dat, value = "A1", n = 6, nbins = 2)), 12)
  expect_equal(nrow(thinPoints(dat, value = "A2", n = 6, nbins = 2)), 11)
  expect_equal(nrow(thinPoints(dat, value = "A1", n = 3, nbins = 2, groupBy = "B")), 19)

  tmp <- thinPoints(dat, value = "A2", n = 3, nbins = 2, groupBy = "B")
  expect_equal(nrow(tmp), 14)
  expect_true(with(tmp, (sum(B == "a") == 3) && (sum(B == "b") == 3) && (sum(B == "c") == 6) && (sum(B == "d") == 2)))
})
