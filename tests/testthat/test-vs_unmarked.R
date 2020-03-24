context("Testing models against unmarked")

test_that("testing closed against unmarked", {
  if(!require(unmarked)) { testthat::skip() }

  set.seed(4321)
  Y <- gen_Nmix_closed(num_sites = 4, num_times = 3, lambda = 3, pdet = 0.5)

  uFrameC <- unmarkedFramePCount(y = Y$nit,
                              siteCovs = data.frame(site1=as.factor(c(0,0,1,1))))
  # NO COVARIATES
  mod1 <- pcount(formula = ~1 ~1,
                     data = uFrameC,
                     K = 6,
                     se = FALSE,
                     method="BFGS")

  mod2 <- redNMix::fit_red_Nmix_closed(nit = Y$nit,
                                       lambda_site_covariates = NULL,
                                       pdet_site_covariates   = NULL,
                                       red = 1,
                                       K = 6,
                                       starts=c(1,0),
                                       method="DFP",
                                       tolerance = 10^-6)
  expect_equal(mod1@opt$value, mod2$f, 10^-3)
  # l_s COVARIATES
  mod1 <- pcount(formula = ~1 ~1+site1,
                 data = uFrameC,
                 K = 6,
                 se = FALSE,
                 method="BFGS")

  mod2 <- redNMix::fit_red_Nmix_closed(nit = Y$nit,
                                       lambda_site_covariates = list(site1=c(0,0,1,1)),
                                       pdet_site_covariates   = NULL,
                                       red = 1,
                                       K = 6,
                                       starts=c(1,0),
                                       method="DFP",
                                       tolerance = 10^-6)
  expect_equal(mod1@opt$value, mod2$f, 10^-3)
  # p_s COVARIATES
  mod1 <- pcount(formula = ~1+site1 ~1,
                 data = uFrameC,
                 K = 6,
                 se = FALSE,
                 method="BFGS")

  mod2 <- redNMix::fit_red_Nmix_closed(nit = Y$nit,
                                       lambda_site_covariates = NULL,
                                       pdet_site_covariates   = list(site1=c(0,0,1,1)),
                                       red = 1,
                                       K = 6,
                                       starts=c(1,0),
                                       method="DFP",
                                       tolerance = 10^-6)
  expect_equal(mod1@opt$value, mod2$f, 10^-3)
  # p_s and l_s COVARIATES
  mod1 <- pcount(formula = ~1+site1 ~1+site1,
                 data = uFrameC,
                 K = 6,
                 se = FALSE,
                 method="BFGS")

  mod2 <- redNMix::fit_red_Nmix_closed(nit = Y$nit,
                                       lambda_site_covariates = list(site1=c(0,0,1,1)),
                                       pdet_site_covariates   = list(site1=c(0,0,1,1)),
                                       red = 1,
                                       K = 6,
                                       starts=c(1,1,0,0),
                                       method="DFP",
                                       tolerance = 10^-6)
  expect_equal(mod1@opt$value, mod2$f, 10^-3)

})


test_that("testing open against unmarked", {
  if(!require(unmarked)) { testthat::skip() }

  set.seed(4321)
  Y <- gen_Nmix_closed(num_sites = 4, num_times = 3, lambda = 3, pdet = 0.5)

  uFrameC <- unmarkedFramePCO(y = Y$nit,
                              siteCovs = data.frame(site1=as.factor(c(0,0,1,1))),
                              obsCovs = NULL,
                              yearlySiteCovs = data.frame(time1=as.factor(rep(c(0,1,0),2))),
                              numPrimary = 3)
  # NO COVARIATES
  mod1 <- pcountOpen(lambdaformula = ~1, gammaformula = ~1, omegaformula = ~1, pformula = ~1,
                     data = uFrameC,
                     K = 6,
                     se = FALSE,
                     method="BFGS")

  mod2 <- redNMix::fit_red_Nmix_open(nit = Y$nit,
                                     lambda_site_covariates = NULL,
                                     gamma_site_covariates  = NULL,
                                     gamma_time_covariates  = NULL,#list(a=rep(c(1,2,3,4,5), each=2)),
                                     omega_site_covariates = NULL,
                                     omega_time_covariates = NULL,#list(time1=c(0,1,2)),
                                     pdet_time_covariates   = NULL,
                                     red = 1,
                                     K = 5,
                                     starts=NULL,
                                     method="BFGS",
                                     tolerance = 10^-6)
  expect_equal(mod1@opt$value, mod2$value, 10^-3)
  # l_s COVARIATES
  mod1 <- pcountOpen(lambdaformula = ~1+site1, gammaformula = ~1, omegaformula = ~1, pformula = ~1,
                     data = uFrameC,
                     K = 6,
                     se = FALSE,
                     method="BFGS")

  mod2 <- redNMix::fit_red_Nmix_open(nit = Y$nit,
                                     lambda_site_covariates = list(site1=c(0,0,1,1)),
                                     gamma_site_covariates  = NULL,
                                     gamma_time_covariates  = NULL,#list(a=rep(c(1,2,3,4,5), each=2)),
                                     omega_site_covariates = NULL,
                                     omega_time_covariates = NULL,#list(time1=c(0,1,2)),
                                     pdet_time_covariates   = NULL,
                                     red = 1,
                                     K = 6,
                                     starts=NULL,
                                     method="DFP",
                                     tolerance = 10^-6)
  expect_equal(mod1@opt$value, mod2$f, 10^-3)
  # o_t COVARIATES
  mod1 <- pcountOpen(lambdaformula = ~1, gammaformula = ~1, omegaformula = ~1+time1, pformula = ~1,
                     data = uFrameC,
                     K = 6,
                     se = FALSE,
                     method="BFGS")

  mod2 <- redNMix::fit_red_Nmix_open(nit = Y$nit,
                                     lambda_site_covariates = NULL,
                                     gamma_site_covariates  = NULL,
                                     gamma_time_covariates  = NULL,#list(a=rep(c(1,2,3,4,5), each=2)),
                                     omega_site_covariates = NULL,
                                     omega_time_covariates = list(time1=c(0,1,0)),#list(time1=c(0,1,2)),
                                     pdet_time_covariates   = NULL,
                                     red = 1,
                                     K = 6,
                                     starts=NULL,
                                     method="BFGS",
                                     tolerance = 10^-6)
  expect_equal(mod1@opt$value, mod2$value, 10^-3)
  # p_t COVARIATES
  mod1 <- pcountOpen(lambdaformula = ~1, gammaformula = ~1, omegaformula = ~1, pformula = ~1+time1,
                     data = uFrameC,
                     K = 6,
                     se = FALSE,
                     method="BFGS")

  mod2 <- redNMix::fit_red_Nmix_open(nit = Y$nit,
                                     lambda_site_covariates = NULL,
                                     gamma_site_covariates  = NULL,
                                     gamma_time_covariates  = NULL,#list(a=rep(c(1,2,3,4,5), each=2)),
                                     omega_site_covariates = NULL,
                                     omega_time_covariates = NULL,#list(time1=c(0,1,2)),
                                     pdet_time_covariates   = list(time1=c(0,1,0)),
                                     red = 1,
                                     K = 6,
                                     starts=NULL,
                                     method="BFGS",
                                     tolerance = 10^-6)
  expect_equal(mod1@opt$value, mod2$value, 10^-3)

})
