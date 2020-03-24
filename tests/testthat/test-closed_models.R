context("Testing fit_red_Nmix_closed")

test_that("testing against unmarked", {

  set.seed(4321)
  Y <- gen_Nmix_closed(num_sites = 2, num_times = 3, lambda = 3, pdet = 0.5)

  if(require(unmarked)) {
    uFrameC <- unmarkedFramePCount(y = Y$nit)
    time0 <- system.time(
      mod0 <- pcount(~1 ~1,
                     data = uFrameC,
                     K = 7,
                     se = FALSE,
                     starts = c(log(25),0), method="BFGS")
    )

    time1 <- system.time(
      mod1 <- fit_red_Nmix_closed(nit = Y$nit,
                                  red = 1,
                                  K = 7,
                                  starts=c(log(25),0),
                                  method="DFP",
                                  APA=FALSE, maxSteps=100, tolerance = 10^-5)
    )

    mod2 <- fit_red_Nmix_closed(nit = Y$nit,
                                red = 1,
                                K = 7,
                                starts=c(log(25),0),
                                method="DFP",
                                APA=TRUE, maxSteps=100, tolerance = 10^-5, outFile = "./test_outputFile.csv")

    expect_equal(file.exists("./test_outputFile.csv"), T)


    expect_equal(as.numeric(mod1$f), as.numeric(mod0@opt$value))
    expect_equal(as.numeric(mod1$x), as.numeric(mod0@opt$par), tolerance=10^-3)
  }
})

test_that("testing closed PARALLEL", {

  set.seed(4321)
  Y <- gen_Nmix_closed(num_sites = 2, num_times = 3, lambda = 1, pdet = 0.5)


  time1 <- system.time(
    mod1 <- fit_red_Nmix_closed(nit = Y$nit,
                                red = 1,
                                K = 3,
                                starts=c(log(25),0),
                                method="DFP",
                                APA=TRUE, precBits=128, tolerance = 10^-5)
  )

  START_PARALLEL(2)
  time2 <- system.time({
    mod2 <- fit_red_Nmix_closed(nit = Y$nit,
                                red = 1,
                                K = 3,
                                starts=c(log(25),0),
                                method="DFP",
                                APA=TRUE, precBits=128, tolerance = 10^-5,
                                PARALLELIZE = TRUE)
  })
  END_PARALLEL()

  expect_equal(as.numeric(mod1$f), as.numeric(mod2$f), tolerance=10^-3)

})



test_that("testing against unmarked: covariates", {

  # site lambda covariate
  set.seed(43215)
  Y1 <- gen_Nmix_closed(num_sites = 1, num_times = 3, lambda = 1, pdet = 0.75)
  Y2 <- gen_Nmix_closed(num_sites = 1, num_times = 3, lambda = 3, pdet = 0.75)
  Y = rbind(Y1$nit,Y2$nit)

  X1 = c(0,1)
  X = data.frame(X1=X1)

  if(require(unmarked)) {
    uFrameC <- unmarkedFramePCount(y = Y, siteCovs = X)
    time0 <- system.time(
      mod0 <- pcount(~1 ~X1,
                     data = uFrameC,
                     K = 7,
                     se = FALSE,
                     starts = c(log(25),0,0), method="BFGS")
    )

    time1 <- system.time(
      mod1 <- fit_red_Nmix_closed(nit = Y,
                                  lambda_site_covariates = list(X1=X1),
                                  red = 1,
                                  K = 7,
                                  starts=c(log(25),0,0),
                                  method="DFP",
                                  APA=FALSE, maxSteps=100, tolerance = 10^-5)
    )

    mod2 <- fit_red_Nmix_closed(nit = Y,
                                lambda_site_covariates = list(X1=X1),
                                red = 1,
                                K = 7,
                                starts=c(log(25),0,0),
                                method="DFP",
                                APA=TRUE, maxSteps=100, tolerance = 10^-5)


    expect_equal(as.numeric(mod1$f), as.numeric(mod0@opt$value))
    expect_equal(as.numeric(mod1$x), as.numeric(mod0@opt$par), tolerance=10^-3)
    expect_equal(as.numeric(mod1$f), as.numeric(mod2$f), tolerance=10^-3)

  }
})

