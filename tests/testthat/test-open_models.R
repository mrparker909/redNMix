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
})
