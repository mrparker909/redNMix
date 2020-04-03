context("Testing fit_red_Nmix_open")

test_that("testing against unmarked", {

  set.seed(4321)
  Y <- gen_Nmix_open(num_sites = 4, num_times = 3, lambda = 3, pdet = 0.5, omega = 0.75, gamma = 1)

  if(require(unmarked)) {
    uFrameC <- unmarkedFramePCO(y = Y$nit, numPrimary = 3)
    time0 <- system.time(
      mod0 <- pcountOpen(lambdaformula = ~1,
                         gammaformula  = ~1,
                         omegaformula  = ~1,
                         pformula      = ~1,
                         data          = uFrameC,
                         K = 6,
                         control=list(maxit=200),
                         se = FALSE,
                         starts = c(log(25),0,0,0), method="BFGS")
    )


    suppressWarnings({
      time1 <- system.time(
        mod1 <- fit_red_Nmix_open(nit = Y$nit,
                                  red = 1,
                                  K = 6,
                                  starts=c(log(25),0,0,0),
                                  method="DFP",
                                  APA=FALSE,
                                  maxSteps=200,
                                  tolerance = 10^-3)
      )
    })

    expect_equal(as.numeric(mod1$f), as.numeric(mod0@opt$value), tolerance=10^-3)
    expect_equal(exp(as.numeric(mod1$x[1])), exp(as.numeric(mod0@opt$par[1])), tolerance=10^-2)
    expect_equal(exp(as.numeric(mod1$x[2])), exp(as.numeric(mod0@opt$par[2])), tolerance=10^-2)
    expect_equal(plogis(as.numeric(mod1$x[3])), plogis(as.numeric(mod0@opt$par[3])), tolerance=10^-2)
    expect_equal(plogis(as.numeric(mod1$x[4])), plogis(as.numeric(mod0@opt$par[4])), tolerance=10^-3)
  }

  expect_error({
    time2 <- system.time(
      mod2 <- fit_red_Nmix_open(nit = Y$nit,
                                red = 1,
                                K = 6,
                                starts=c(log(25),0,0,0,0),
                                method="DFP",
                                lambda_site_covariates = list(c(0,0,1,1)),
                                fixed_omega = 0.75,
                                APA=FALSE,
                                maxSteps=5, keepValues=T,
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
                                APA=FALSE,
                                maxSteps=5, keepValues=T,
                                tolerance = 10^-3)
    )
  }, NA)

  expect_error({
    time4 <- system.time(
      mod4 <- fit_red_Nmix_open(nit = Y$nit,
                                red = 1,
                                K = 6,
                                starts=c(log(25),0,0,0,0),
                                method="DFP",
                                pdet_site_covariates = list(c(0,0,1,1)),
                                fixed_gamma = 1,
                                APA=FALSE,
                                maxSteps=5, keepValues=T,
                                tolerance = 10^-3)
    )
  }, NA)


  expect_error({
    time5 <- system.time(
      mod5 <- fit_red_Nmix_open(nit = Y$nit,
                                red = 1,
                                K = 6,
                                starts=c(log(25),0,0,0,0),
                                method="DFP",
                                omega_site_covariates = list(c(0,0,1,1)),
                                fixed_gamma = 0.25,
                                APA=FALSE,
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
                                APA=FALSE,
                                maxSteps=5, keepValues=T,
                                tolerance = 10^-3)
    )
  }, NA)
})
