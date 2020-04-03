context("Testing likelihoods: open")

test_that("testing open: APA vs NAPA", {

  set.seed(4321)
  Y <- gen_Nmix_closed(num_sites = 2, num_times = 3, lambda = 3, pdet = 0.5)
  red = matrix(1, nrow=2, ncol=3)
  K   = matrix(5, nrow=2, ncol=3)

  ll_NAPA <- redNMix::red_Like_open(
    par=c(log(25),log(1),0,0),
    nit = Y$nit,
    l_s_c = NULL, g_s_c = NULL, g_t_c = NULL, o_s_c = NULL, o_t_c = NULL, p_s_c = NULL, p_t_c = NULL,
    K = K,
    red = red,
    APA = FALSE
  )

  ll_APA <- redNMix::red_Like_open(
    par=c(log(25),log(1),0,0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = NULL, g_s_c = NULL, g_t_c = NULL, o_s_c = NULL, o_t_c = NULL, p_t_c = NULL,
    K = K,
    red = red,
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
    K = K,
    red = red,
    APA = FALSE,
    precBits = 53
  )
  ll_APA <- redNMix::red_Like_open(
    par=c(1,1,0,0,0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = NULL, g_s_c = NULL, g_t_c = NULL,
    o_s_c = NULL, o_t_c = NULL, p_t_c = list(c(0,1,1)),
    K = K,
    red = red,
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
    K = K,
    red = red,
    APA = FALSE,
    precBits = 53
  )
  ll_APA <- redNMix::red_Like_open(
    par=c(1,1,0,0,0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = NULL, g_s_c = NULL, g_t_c = NULL,
    o_s_c = NULL, o_t_c = list(c(0,1,1)), p_t_c = NULL,
    K = K,
    red = red,
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
    K = K,
    red = red,
    APA = FALSE,
    precBits = 53
  )
  ll_APA <- redNMix::red_Like_open(
    par=c(1,1,0,0,0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = NULL, g_s_c = NULL, g_t_c = list(c(0,1,1)),
    o_s_c = NULL, o_t_c = NULL, p_t_c = NULL,
    K = K,
    red = red,
    APA = TRUE,
    precBits = 53
  )
  expect_equal(as.numeric(ll_NAPA), as.numeric(ll_APA), 10^-3)



  # testing l_s_c
  ll_NAPA <- redNMix::red_Like_open(
    par=c(1,1,0,0,0),
    nit = Y$nit,
    l_s_c = list(c(0,1)), p_s_c = NULL, g_s_c = NULL, g_t_c = NULL,
    o_s_c = NULL, o_t_c = NULL, p_t_c = NULL,
    K = K,
    red = red,
    APA = FALSE,
    precBits = 53
  )
  ll_APA <- redNMix::red_Like_open(
    par=c(1,1,0,0,0),
    nit = Y$nit,
    l_s_c = list(c(0,1)), p_s_c = NULL, g_s_c = NULL, g_t_c = NULL,
    o_s_c = NULL, o_t_c = NULL, p_t_c = NULL,
    K = K,
    red = red,
    APA = TRUE,
    precBits = 53
  )
  expect_equal(as.numeric(ll_NAPA), as.numeric(ll_APA), 10^-3)


  # testing p_s_c
  ll_NAPA <- redNMix::red_Like_open(
    par=c(1,1,0,0,0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = list(c(0,1)), g_s_c = NULL, g_t_c = NULL,
    o_s_c = NULL, o_t_c = NULL, p_t_c = NULL,
    K = K,
    red = red,
    APA = FALSE,
    precBits = 53
  )
  ll_APA <- redNMix::red_Like_open(
    par=c(1,1,0,0,0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = list(c(0,1)), g_s_c = NULL, g_t_c = NULL,
    o_s_c = NULL, o_t_c = NULL, p_t_c = NULL,
    K = K,
    red = red,
    APA = TRUE,
    precBits = 53
  )
  expect_equal(as.numeric(ll_NAPA), as.numeric(ll_APA), 10^-3)


  # testing g_s_c
  ll_NAPA <- redNMix::red_Like_open(
    par=c(1,1,0,0,0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = NULL, g_s_c = list(c(0,1)), g_t_c = NULL,
    o_s_c = NULL, o_t_c = NULL, p_t_c = NULL,
    K = K,
    red = red,
    APA = FALSE,
    precBits = 53
  )
  ll_APA <- redNMix::red_Like_open(
    par=c(1,1,0,0,0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = NULL, g_s_c = list(c(0,1)), g_t_c = NULL,
    o_s_c = NULL, o_t_c = NULL, p_t_c = NULL,
    K = K,
    red = red,
    APA = TRUE,
    precBits = 53
  )
  expect_equal(as.numeric(ll_NAPA), as.numeric(ll_APA), 10^-3)



  # testing p_s_c and p_t_c
  ll_NAPA <- redNMix::red_Like_open(
    par=c(1,1,0,0,0,0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = list(c(0,1)), g_s_c = NULL, g_t_c = NULL,
    o_s_c = NULL, o_t_c = NULL, p_t_c = list(c(0,1,1)),
    K = K,
    red = red,
    APA = FALSE,
    precBits = 53
  )
  ll_APA <- redNMix::red_Like_open(
    par=c(1,1,0,0,0,0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = list(c(0,1)), g_s_c = NULL, g_t_c = NULL,
    o_s_c = NULL, o_t_c = NULL, p_t_c = list(c(0,1,1)),
    K = K,
    red = red,
    APA = TRUE,
    precBits = 53
  )
  expect_equal(as.numeric(ll_NAPA), as.numeric(ll_APA), 10^-3)



  # testing ALL covariates
  ll_NAPA <- redNMix::red_Like_open(
    par=c(1,1,0,0,0,0,0,0,0,0,0),
    nit = Y$nit,
    l_s_c = list(c(0,1)), p_s_c = list(c(0,1)), g_s_c = list(c(0,1)), g_t_c = list(c(0,1,1)),
    o_s_c = list(c(0,1)), o_t_c = list(c(0,1,1)), p_t_c = list(c(0,1,1)),
    K = K,
    red = red,
    APA = FALSE,
    precBits = 53
  )
  ll_APA <- redNMix::red_Like_open(
    par=c(1,1,0,0,0,0,0,0,0,0,0),
    nit = Y$nit,
    l_s_c = list(c(0,1)), p_s_c = list(c(0,1)), g_s_c = list(c(0,1)), g_t_c = list(c(0,1,1)),
    o_s_c = list(c(0,1)), o_t_c = list(c(0,1,1)), p_t_c = list(c(0,1,1)),
    K = K,
    red = red,
    APA = TRUE,
    precBits = 53
  )
  expect_equal(as.numeric(ll_NAPA), as.numeric(ll_APA), 10^-3)

  # testing ALL covariates, K, red
  red = matrix(c(1,2), nrow=2, ncol=3)
  K   = matrix(c(4,8), nrow=2, ncol=3)

  ll_NAPA <- redNMix::red_Like_open(
    par=c(1,1,0,0,0,0,0,0,0,0,0),
    nit = Y$nit,
    l_s_c = list(c(0,1)), p_s_c = list(c(0,1)), g_s_c = list(c(0,1)), g_t_c = list(c(0,1,1)),
    o_s_c = list(c(0,1)), o_t_c = list(c(0,1,1)), p_t_c = list(c(0,1,1)),
    K = K,
    red = red,
    APA = FALSE,
    precBits = 53
  )
  ll_APA <- redNMix::red_Like_open(
    par=c(1,1,0,0,0,0,0,0,0,0,0),
    nit = Y$nit,
    l_s_c = list(c(0,1)), p_s_c = list(c(0,1)), g_s_c = list(c(0,1)), g_t_c = list(c(0,1,1)),
    o_s_c = list(c(0,1)), o_t_c = list(c(0,1,1)), p_t_c = list(c(0,1,1)),
    K = K,
    red = red,
    APA = TRUE,
    precBits = 53
  )
  expect_equal(as.numeric(ll_NAPA), as.numeric(ll_APA), 10^-3)

})




