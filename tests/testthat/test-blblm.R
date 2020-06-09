test_that("Our blblm package outputs what we expect it to", {
  m = 3
  B = 100
  fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100, parallel = TRUE)

  lmFit <- lm(mpg ~ wt * hp, data = mtcars)

  expect_s3_class(fit, "blblm")

  #expect_equivalent(coef(fit), coef(lmFit))
  expect_identical(dim(coef(lmFit)), dim(coef(fit)))
})
