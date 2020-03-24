context("Testing fit_red_Nmix_open")

test_that("testing APA against NAPA models", {

  set.seed(4321)
  Y <- gen_Nmix_open(num_sites = 4, num_times = 3, lambda = 3, pdet = 0.85, omega = 0.75, gamma = 1)

  suppressWarnings({
    time1 <- system.time(
      mod1 <- fit_red_Nmix_open(nit = Y$nit,
                                  red = 1,
                                  K = 6,
                                  starts=c(log(5),0,0,0),
                                  method="DFP",
                                  APA=FALSE, maxSteps=5, tolerance = 10^-5)
    )
  })

  suppressWarnings({
    time2 <- system.time(
      mod2 <- fit_red_Nmix_open(nit = Y$nit,
                                red = 1,
                                K = 6,
                                starts=c(log(5),0,0,0),
                                method="DFP",
                                APA=TRUE, precBits = 53,
                                maxSteps=5, tolerance = 10^-5)
    )
  })

  expect_equal(as.numeric(mod1$x), as.numeric(mod2$x), tolerance=10^-4)

})
