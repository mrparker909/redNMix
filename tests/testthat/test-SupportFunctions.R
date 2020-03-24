context("Testing Support Functions")

test_that("Support Functions", {
  set.seed(1234)
  nit = gen_Nmix_closed(5,5,50,0.5)$nit

  expect_equal(as.numeric(mcd(nit,1)), 25)
})
