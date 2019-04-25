context("Testing fit_red_Nmix_closed")

test_that("testing against unmarked", {

  set.seed(4321)
  Y <- gen_Nmix_closed(num_sites = 2, num_times = 3, lambda = 3, pdet = 0.5)

  if(require(unmarked)) {
    uFrameC <- unmarkedFramePCount(y = Y$nit)
    time0 <- system.time(
      mod0 <- pcount(~1 ~1,
                     data = uFrameC,
                     K = 6,
                     se = FALSE,
                     starts = c(log(25),0), method="BFGS")
    )

    time1 <- system.time(
      mod1 <- fit_red_Nmix_closed(nit = Y$nit,
                                  red = 1,
                                  K = 6,
                                  starts=c(log(25),0),
                                  method="DFP",
                                  APA=FALSE, maxSteps=100, tolerance = 10^-3)
    )
    # # ~300 seconds
    # time2 <- system.time(
    #   mod2 <- fit_red_Nmix_closed(nit = Y$nit,
    #                               red = 1,
    #                               K = 6,
    #                               starts=c(log(25),0),
    #                               method="DFP",
    #                               APA=TRUE, precBits=128)
    # )
    expect_equal(as.numeric(mod1$f), as.numeric(mod0@opt$value))
    expect_equal(as.numeric(mod1$x), as.numeric(mod0@opt$par), tolerance=10^-3)
  }
})
