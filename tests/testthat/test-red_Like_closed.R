context("Testing likelihoods: closed")

test_that("testing closed: APA vs NAPA", {
  set.seed(4321)
  Y <- gen_Nmix_closed(num_sites = 2, num_times = 3, lambda = 3, pdet = 0.5)
  red = matrix(1, nrow=2, ncol=3)
  K   = matrix(50, nrow=2,ncol=3)

  ll_NAPA <- redNMix::red_Like_closed(
    par=c(log(25),0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = NULL,
    K = K,
    red = red,
    APA = FALSE
  )

  ll_APA <- redNMix::red_Like_closed(
    par=c(log(25),0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = NULL,
    K = K,
    red = red,
    APA = TRUE
  )

  expect_equal(as.numeric(ll_NAPA), as.numeric(ll_APA))


  ll_NAPA <- redNMix::red_Like_closed(
    par=c(log(25),0.5,0),
    nit = Y$nit,
    l_s_c = list(c(0,1)), p_s_c = NULL,
    K = K,
    red = red,
    APA = FALSE
  )

  ll_APA <- redNMix::red_Like_closed(
    par=c(log(25),0.5,0),
    nit = Y$nit,
    l_s_c = list(x=c(0,1)), p_s_c = NULL,
    K = K,
    red = red,
    APA = TRUE
  )

  expect_equal(as.numeric(ll_NAPA), as.numeric(ll_APA))


  ll_NAPA <- redNMix::red_Like_closed(
    par=c(log(25),0.5,0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = list(c(0,1)),
    K = K,
    red = red,
    APA = FALSE
  )

  ll_APA <- redNMix::red_Like_closed(
    par=c(log(25),0.5,0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = list(c(0,1)),
    K = K,
    red = red,
    APA = TRUE
  )

  expect_equal(as.numeric(ll_NAPA), as.numeric(ll_APA))


  ll_NAPA <- redNMix::red_Like_closed(
    par=c(log(25),0.5,0,0),
    nit = Y$nit,
    l_s_c = list(c(0,1)), p_s_c = list(c(0,1)),
    K = K,
    red = red,
    APA = FALSE
  )

  ll_APA <- redNMix::red_Like_closed(
    par=c(log(25),0.5,0,0),
    nit = Y$nit,
    l_s_c = list(c(0,1)), p_s_c = list(c(0,1)),
    K = K,
    red = red,
    APA = TRUE
  )

  expect_equal(as.numeric(ll_NAPA), as.numeric(ll_APA))


  ## testing reduction factor
  red = matrix(1, nrow=2, ncol=3)
  red[2,3] <- 5

  ll_NAPA <- redNMix::red_Like_closed(
    par=c(log(25),0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = NULL,
    K = K,
    red = red,
    APA = FALSE
  )

  ll_APA <- redNMix::red_Like_closed(
    par=c(log(25),0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = NULL,
    K = K,
    red = red,
    APA = TRUE
  )

  expect_equal(as.numeric(ll_NAPA), as.numeric(ll_APA))


  ll_NAPA <- redNMix::red_Like_closed(
    par=c(log(25),0.5,0),
    nit = Y$nit,
    l_s_c = list(c(0,1)), p_s_c = NULL,
    K = K,
    red = red,
    APA = FALSE
  )

  ll_APA <- redNMix::red_Like_closed(
    par=c(log(25),0.5,0),
    nit = Y$nit,
    l_s_c = list(x=c(0,1)), p_s_c = NULL,
    K = K,
    red = red,
    APA = TRUE
  )

  expect_equal(as.numeric(ll_NAPA), as.numeric(ll_APA))


  ll_NAPA <- redNMix::red_Like_closed(
    par=c(log(25),0.5,0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = list(c(0,1)),
    K = K,
    red = red,
    APA = FALSE
  )

  ll_APA <- redNMix::red_Like_closed(
    par=c(log(25),0.5,0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = list(c(0,1)),
    K = K,
    red = red,
    APA = TRUE
  )

  expect_equal(as.numeric(ll_NAPA), as.numeric(ll_APA))


  ll_NAPA <- redNMix::red_Like_closed(
    par=c(log(25),0.5,0,0),
    nit = Y$nit,
    l_s_c = list(c(0,1)), p_s_c = list(c(0,1)),
    K = K,
    red = red,
    APA = FALSE
  )

  ll_APA <- redNMix::red_Like_closed(
    par=c(log(25),0.5,0,0),
    nit = Y$nit,
    l_s_c = list(c(0,1)), p_s_c = list(c(0,1)),
    K = K,
    red = red,
    APA = TRUE
  )

  expect_equal(as.numeric(ll_NAPA), as.numeric(ll_APA))


  ## testing K
  red = matrix(1, nrow=2, ncol=3)
  K   = matrix(50, nrow=2,ncol=3)
  K[2,3] <- 60

  ll_NAPA <- redNMix::red_Like_closed(
    par=c(log(25),0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = NULL,
    K = K,
    red = red,
    APA = FALSE
  )

  ll_APA <- redNMix::red_Like_closed(
    par=c(log(25),0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = NULL,
    K = K,
    red = red,
    APA = TRUE
  )

  expect_equal(as.numeric(ll_NAPA), as.numeric(ll_APA))


  ll_NAPA <- redNMix::red_Like_closed(
    par=c(log(25),0.5,0),
    nit = Y$nit,
    l_s_c = list(c(0,1)), p_s_c = NULL,
    K = K,
    red = red,
    APA = FALSE
  )

  ll_APA <- redNMix::red_Like_closed(
    par=c(log(25),0.5,0),
    nit = Y$nit,
    l_s_c = list(x=c(0,1)), p_s_c = NULL,
    K = K,
    red = red,
    APA = TRUE
  )

  expect_equal(as.numeric(ll_NAPA), as.numeric(ll_APA))


  ll_NAPA <- redNMix::red_Like_closed(
    par=c(log(25),0.5,0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = list(c(0,1)),
    K = K,
    red = red,
    APA = FALSE
  )

  ll_APA <- redNMix::red_Like_closed(
    par=c(log(25),0.5,0),
    nit = Y$nit,
    l_s_c = NULL, p_s_c = list(c(0,1)),
    K = K,
    red = red,
    APA = TRUE
  )

  expect_equal(as.numeric(ll_NAPA), as.numeric(ll_APA))


  ll_NAPA <- redNMix::red_Like_closed(
    par=c(log(25),0.5,0,0),
    nit = Y$nit,
    l_s_c = list(c(0,1)), p_s_c = list(c(0,1)),
    K = K,
    red = red,
    APA = FALSE
  )

  ll_APA <- redNMix::red_Like_closed(
    par=c(log(25),0.5,0,0),
    nit = Y$nit,
    l_s_c = list(c(0,1)), p_s_c = list(c(0,1)),
    K = K,
    red = red,
    APA = TRUE
  )

  expect_equal(as.numeric(ll_NAPA), as.numeric(ll_APA))
})

