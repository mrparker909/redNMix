context("Testing likelihoods")

test_that("testing closed: APA vs NAPA", {
  set.seed(4321)
  Y <- gen_Nmix_closed(num_sites = 2, num_times = 3, lambda = 3, pdet = 0.5)


  ll_NAPA <- redNMix::red_Like_closed(
    par=c(log(25),0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = NULL,
    K = matrix(50, nrow=3),
    red = matrix(1, nrow=3),
    APA = FALSE
  )

  ll_APA <- redNMix::red_Like_closed(
    par=c(log(25),0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = NULL,
    K = matrix(50, nrow=3),
    red = matrix(1, nrow=3),
    APA = TRUE
  )

    expect_equal(as.numeric(ll_NAPA), as.numeric(ll_APA))
})



test_that("testing open: APA vs NAPA", {

  set.seed(4321)
  Y <- gen_Nmix_closed(num_sites = 2, num_times = 3, lambda = 3, pdet = 0.5)

  ll_NAPA <- redNMix::red_Like_open(
    par=c(log(25),log(1),0,0),
    nit = Y$nit,
    l_s_c = NULL, g_s_c = NULL, g_t_c = NULL, o_s_c = NULL, o_t_c = NULL, p_s_c = NULL, p_t_c = NULL,
    K = matrix(5, nrow=3),
    red = matrix(1, nrow=2, ncol=3),
    APA = FALSE
  )

  ll_APA <- redNMix::red_Like_open(
    par=c(log(25),log(1),0,0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = NULL, g_s_c = NULL, g_t_c = NULL, o_s_c = NULL, o_t_c = NULL, p_t_c = NULL,
    K = matrix(5, nrow=3),
    red = matrix(1, nrow=2, ncol=3),
    APA = TRUE,
    precBits = 53
  )

  expect_equal(as.numeric(ll_NAPA), as.numeric(ll_APA), 10^-3)

  # testing p_t
  ll_NAPA <- redNMix::red_Like_open(
    par=c(1,1,0,0,0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = NULL, g_s_c = NULL, g_t_c = NULL,
    o_s_c = NULL, o_t_c = NULL, p_t_c = list(c(0,1,1)),
    K = matrix(5, nrow=2),
    red = matrix(1, nrow=2, ncol=10),
    APA = FALSE,
    precBits = 53
  )
  ll_APA <- redNMix::red_Like_open(
    par=c(1,1,0,0,0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = NULL, g_s_c = NULL, g_t_c = NULL,
    o_s_c = NULL, o_t_c = NULL, p_t_c = list(c(0,1,1)),
    K = matrix(5, nrow=2),
    red = matrix(1, nrow=2, ncol=10),
    APA = TRUE,
    precBits = 53
  )
  expect_equal(as.numeric(ll_NAPA), as.numeric(ll_APA), 10^-3)

  # testing o_t
  ll_NAPA <- redNMix::red_Like_open(
    par=c(1,1,0,0,0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = NULL, g_s_c = NULL, g_t_c = NULL,
    o_s_c = NULL, o_t_c = list(c(0,1,1)), p_t_c = NULL,
    K = matrix(5, nrow=2),
    red = matrix(1, nrow=2, ncol=10),
    APA = FALSE,
    precBits = 53
  )
  ll_APA <- redNMix::red_Like_open(
    par=c(1,1,0,0,0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = NULL, g_s_c = NULL, g_t_c = NULL,
    o_s_c = NULL, o_t_c = list(c(0,1,1)), p_t_c = NULL,
    K = matrix(5, nrow=2),
    red = matrix(1, nrow=2, ncol=10),
    APA = TRUE,
    precBits = 53
  )
  expect_equal(as.numeric(ll_NAPA), as.numeric(ll_APA), 10^-3)


  # testing g_t
  ll_NAPA <- redNMix::red_Like_open(
    par=c(1,1,0,0,0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = NULL, g_s_c = NULL, g_t_c = list(c(0,1,1)),
    o_s_c = NULL, o_t_c = NULL, p_t_c = NULL,
    K = matrix(5, nrow=2),
    red = matrix(1, nrow=2, ncol=10),
    APA = FALSE,
    precBits = 53
  )
  ll_APA <- redNMix::red_Like_open(
    par=c(1,1,0,0,0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = NULL, g_s_c = NULL, g_t_c = list(c(0,1,1)),
    o_s_c = NULL, o_t_c = NULL, p_t_c = NULL,
    K = matrix(5, nrow=2),
    red = matrix(1, nrow=2, ncol=10),
    APA = TRUE,
    precBits = 53
  )
  expect_equal(as.numeric(ll_NAPA), as.numeric(ll_APA), 10^-3)

})





# ### DELETE WHEN FIXED
# set.seed(4321)
# Y <- gen_Nmix_closed(num_sites = 4, num_times = 4, lambda = 3, pdet = 0.5)
#
# ll_NAPA <- redNMix::red_Like_open(
#   par=c(1,1,0,0,0),
#   nit = Y$nit,
#   l_s_c = NULL, p_s_c = NULL, g_s_c = NULL, g_t_c = NULL,
#   o_s_c = NULL, o_t_c = NULL, p_t_c = list(c(0,1,1)),
#   K = matrix(5, nrow=4),
#   red = matrix(1, nrow=4, ncol=4),
#   APA = FALSE,
#   precBits = 53
# )
#
# ll_NAPA <- red_Like_open2(
#   par=c(1,1,0,0,0),
#   nit = Y$nit,
#   l_s_c = NULL, p_s_c = NULL, g_s_c = NULL, g_t_c = NULL,
#   o_s_c = NULL, o_t_c = NULL, p_t_c = list(c(0,1,1)),
#   K = matrix(5, nrow=2),
#   red = matrix(1, nrow=2, ncol=10),
#   APA = FALSE,
#   precBits = 53
# )
#
# uFrameC <- unmarkedFramePCO(y = Y$nit,
#                             siteCovs = data.frame(site1=as.factor(c(0,0,1,1))),
#                             obsCovs = NULL,
#                             yearlySiteCovs = data.frame(time1=as.factor(rep(c(0,0,1,1),2))),
#                             numPrimary = 4)
# # NO COVARIATES
# mod1 <- pcountOpen(lambdaformula = ~1, gammaformula = ~1, omegaformula = ~1, pformula = ~1,
#                    data = uFrameC,
#                    K = 6,
#                    se = FALSE,
#                    method="BFGS")
#
# mod2 <- fit_red_Nmix_open2(nit = Y$nit,
#                                    lambda_site_covariates = NULL,#list(site1=c(0,0,1,1)),
#                                    gamma_site_covariates  = NULL,#list(site1=c(0,0,1,1)),
#                                    gamma_time_covariates  = NULL,#list(a=rep(c(1,2,3,4,5), each=2)),
#                                    omega_site_covariates = NULL,
#                                    omega_time_covariates = NULL,#list(time1=c(0,1,2)),
#                                    pdet_time_covariates   = NULL,#list(time1=c(0,0,1,1)),
#                                    red = 1,
#                                    K = 6,
#                                    starts=NULL,
#                                    method="BFGS",
#                                    tolerance = 10^-6)
#
#
# # mod3 <- redNMix::fit_red_Nmix_open(nit = Y$nit,
# #                                    lambda_site_covariates = NULL,#list(site1=c(0,0,1,1)),
# #                                    gamma_site_covariates  = NULL,#list(site1=c(0,0,1,1)),
# #                                    gamma_time_covariates  = NULL,#list(a=rep(c(1,2,3,4,5), each=2)),
# #                                    omega_site_covariates = NULL,
# #                                    omega_time_covariates = NULL,#list(time1=c(0,1,2)),
# #                                    pdet_time_covariates   = NULL,#list(time1=c(0,0,1,1)),
# #                                    red = 1,
# #                                    K = 6,
# #                                    starts=NULL,
# #                                    method="DFP",
# #                                    APA = TRUE,
# #                                    precBits = 64,
# #                                    tolerance = 10^-4)
# mod1@opt$par
# mod1@opt$value
# mod2$par
# mod2$value
