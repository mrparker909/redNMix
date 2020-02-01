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

    # time2 <- system.time(
    #   mod2 <- fit_red_Nmix_closed(nit = Y$nit,
    #                               red = 1,
    #                               K = 7,
    #                               starts=c(log(25),0),
    #                               method="DFP",
    #                               APA=TRUE, precBits=128, tolerance = 10^-5)
    # )
    #
    # time3 <- system.time({
    #   START_PARALLEL(2)
    #   mod3 <- fit_red_Nmix_closed(nit = Y$nit,
    #                               red = 1,
    #                               K = 7,
    #                               starts=c(log(25),0),
    #                               method="DFP",
    #                               APA=TRUE, precBits=128, PARALLELIZE = TRUE)
    #   END_PARALLEL()}
    # )

    expect_equal(as.numeric(mod1$f), as.numeric(mod0@opt$value))
    expect_equal(as.numeric(mod1$x), as.numeric(mod0@opt$par), tolerance=10^-3)
  }
})
