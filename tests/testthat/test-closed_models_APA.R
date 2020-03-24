context("Testing fit_red_Nmix_closed")

test_that("testing APA against NAPA models", {

  set.seed(4321)
  Y <- gen_Nmix_closed(num_sites = 2, num_times = 3, lambda = 3, pdet = 0.5)

  time1 <- system.time(
    mod1 <- fit_red_Nmix_closed(nit = Y$nit,
                                red = 2,
                                K = 7,
                                starts=c(log(25),0),
                                method="DFP",
                                APA=FALSE, maxSteps=100, tolerance = 10^-5)
  )

  time2 <- system.time(
    mod2 <- fit_red_Nmix_closed(nit = Y$nit,
                                red = 2,
                                K = 7,
                                starts=c(log(25),0),
                                method="DFP",
                                APA=TRUE, precBits = 128,
                                maxSteps=100, tolerance = 10^-5)
  )

  expect_equal(as.numeric(mod1$x), as.numeric(mod2$x), tolerance=10^-2)

})
