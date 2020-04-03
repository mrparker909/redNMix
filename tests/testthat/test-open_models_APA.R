context("Testing fit_red_Nmix_open APA")

test_that("testing APA against NAPA models", {

  set.seed(4321)
  Y <- gen_Nmix_open(num_sites = 4, num_times = 3, lambda = 3, pdet = 0.85, omega = 0.75, gamma = 1)

  suppressWarnings({
    time1 <- system.time(
      mod1 <- fit_red_Nmix_open(nit = Y$nit,
                                  red = 2,
                                  K = 5,
                                  starts=c(log(5),0,0,0),
                                  method="DFP",
                                  APA=FALSE, maxSteps=2, tolerance = 10^-5)
    )
  })

  suppressWarnings({
    time2 <- system.time(
      mod2 <- fit_red_Nmix_open(nit = Y$nit,
                                red = 2,
                                K = 5,
                                starts=c(log(5),0,0,0),
                                method="DFP",
                                APA=TRUE, precBits = 53,
                                maxSteps=2, tolerance = 10^-5)
    )
  })

  expect_equal(as.numeric(mod1$x), as.numeric(mod2$x), tolerance=10^-4)


  expect_error({
    time3 <- system.time(
      mod3 <- fit_red_Nmix_open(nit = Y$nit,
                                red = 2,
                                K = 5,
                                starts=c(log(25),0,0,0,0),
                                method="DFP",
                                lambda_site_covariates = list(c(0,0,1,1)),
                                fixed_omega = 0.75,
                                APA=TRUE, precBits = 53,
                                maxSteps=2, keepValues=T,
                                tolerance = 10^-3)
    )
  }, NA)

  expect_error({
    time3 <- system.time(
      mod3 <- fit_red_Nmix_open(nit = Y$nit,
                                red = 1,
                                K = 6,
                                starts=c(log(25),0,0,0,0),
                                method="DFP",
                                lambda_site_covariates = list(c(0,0,1,1)),
                                fixed_gamma = 0,
                                APA=T,
                                maxSteps=5, keepValues=T,
                                tolerance = 10^-3)
    )
  }, NA)

  expect_error({
    time6 <- system.time(
      mod6 <- fit_red_Nmix_open(nit = Y$nit,
                                red = 1,
                                K = 6,
                                starts=c(log(25),0,0,0,0),
                                method="DFP",
                                gamma_time_covariates = list(c(0,1,2)),
                                fixed_gamma = 0.25,
                                APA=T,
                                maxSteps=5, keepValues=T,
                                tolerance = 10^-3)
    )
  }, NA)
})
