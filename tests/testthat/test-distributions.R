context("Testing APA and NAPA Distributions")

test_that("drbinom_ and drbinomAPA_", {
  expect_equal(drbinom(x=0:5, size=5, prob=0.5, red=1, log=FALSE),dbinom(x=0:5, size=5, prob=0.5, log=FALSE))
  expect_equal(drbinom(x=5, size=0:10, prob=0.5, red=1, log=FALSE),dbinom(x=5, size=0:10, prob=0.5, log=FALSE))
  expect_equal(drbinom(x=1:5, size=5, prob=0.5, red=1, log=FALSE),as.numeric(drbinomAPA(x=1:5, size=5, prob=0.5, red=1, log=FALSE, precBits=53)))

  # NAPA vs base R
  expect_equal(drbinom(x=reduction(0:19,1), size=19, prob=0.25, red=1, log=FALSE),dbinom(x=0:19, size=19, prob=0.25, log=FALSE))
  expect_equal(sum(drbinom(x=0:10, size=19, prob=0.25, red=2, log=FALSE)),sum(dbinom(x=0:19, size=19, prob=0.25, log=FALSE)))
  expect_equal(drbinom(x=0, size=19, prob=0.25, red=2, log=FALSE),sum(dbinom(x=0:0, size=19, prob=0.25, log=FALSE)))
  expect_equal(drbinom(x=1, size=19, prob=0.25, red=2, log=FALSE),sum(dbinom(x=1:2, size=19, prob=0.25, log=FALSE)))
  expect_equal(drbinom(x=2, size=19, prob=0.25, red=2, log=FALSE),sum(dbinom(x=3:4, size=19, prob=0.25, log=FALSE)))
  expect_equal(drbinom(x=3, size=19, prob=0.25, red=2, log=FALSE),sum(dbinom(x=5:6, size=19, prob=0.25, log=FALSE)))
  expect_equal(drbinom(x=4, size=19, prob=0.25, red=2, log=FALSE),sum(dbinom(x=7:8, size=19, prob=0.25, log=FALSE)))
  expect_equal(drbinom(x=5, size=19, prob=0.25, red=2, log=FALSE),sum(dbinom(x=9:10, size=19, prob=0.25, log=FALSE)))
  expect_equal(drbinom(x=6, size=19, prob=0.25, red=2, log=FALSE),sum(dbinom(x=11:12, size=19, prob=0.25, log=FALSE)))
  expect_equal(drbinom(x=7, size=19, prob=0.25, red=2, log=FALSE),sum(dbinom(x=13:14, size=19, prob=0.25, log=FALSE)))
  expect_equal(drbinom(x=8, size=19, prob=0.25, red=2, log=FALSE),sum(dbinom(x=15:16, size=19, prob=0.25, log=FALSE)))
  expect_equal(drbinom(x=9, size=19, prob=0.25, red=2, log=FALSE),sum(dbinom(x=17:18, size=19, prob=0.25, log=FALSE)))
  expect_equal(drbinom(x=10,size=19, prob=0.25, red=2, log=FALSE),sum(dbinom(x=19, size=19, prob=0.25, log=FALSE)))

  expect_equal(sum(drbinom(x=0:3, size=10, prob=0.25, red=3, log=FALSE)),sum(dbinom(x=0:10, size=10, prob=0.75, log=FALSE)))
  expect_equal(drbinom(x=0, size=10, prob=0.75, red=3, log=FALSE),sum(dbinom(x=0:1, size=10, prob=0.75, log=FALSE)))
  expect_equal(drbinom(x=1, size=10, prob=0.75, red=3, log=FALSE),sum(dbinom(x=2:4, size=10, prob=0.75, log=FALSE)))
  expect_equal(drbinom(x=2, size=10, prob=0.75, red=3, log=FALSE),sum(dbinom(x=5:7, size=10, prob=0.75, log=FALSE)))
  expect_equal(drbinom(x=3, size=10, prob=0.75, red=3, log=FALSE),sum(dbinom(x=8:10, size=10, prob=0.75, log=FALSE)))

  expect_equal(sum(drbinom(x=0:3, size=10, prob=0.5, red=4, log=FALSE)),sum(dbinom(x=0:10, size=10, prob=0.5, log=FALSE)))
  expect_equal(drbinom(x=0, size=10, prob=0.5, red=4, log=FALSE),sum(dbinom(x=0:1, size=10, prob=0.5, log=FALSE)))
  expect_equal(drbinom(x=1, size=10, prob=0.5, red=4, log=FALSE),sum(dbinom(x=2:5, size=10, prob=0.5, log=FALSE)))
  expect_equal(drbinom(x=2, size=10, prob=0.5, red=4, log=FALSE),sum(dbinom(x=6:9, size=10, prob=0.5, log=FALSE)))
  expect_equal(drbinom(x=3, size=10, prob=0.5, red=4, log=FALSE),sum(dbinom(x=10, size=10, prob=0.5, log=FALSE)))

  expect_equal(sum(drbinom(x=0:2, size=10, prob=0.5, red=5, log=FALSE)),sum(dbinom(x=0:10, size=10, prob=0.5, log=FALSE)))
  expect_equal(drbinom(x=0, size=10, prob=0.5, red=5, log=FALSE),sum(dbinom(x=0:2, size=10, prob=0.5, log=FALSE)))
  expect_equal(drbinom(x=1, size=10, prob=0.5, red=5, log=FALSE),sum(dbinom(x=3:7, size=10, prob=0.5, log=FALSE)))
  expect_equal(drbinom(x=2, size=10, prob=0.5, red=5, log=FALSE),sum(dbinom(x=8:10, size=10, prob=0.5, log=FALSE)))

  # APA vs NAPA
  expect_equal(sum(drbinom(x=0:10, size=19, prob=0.25, red=2, log=FALSE)),as.numeric(sum(drbinomAPA(x=0:10, size=19, prob=0.25, red=2, log=FALSE, precBits=53))))
  expect_equal(drbinom(x=0, size=19, prob=0.25, red=2, log=FALSE),as.numeric(drbinomAPA(x=0, size=19, prob=0.25, red=2, log=FALSE, precBits=53)))
  expect_equal(drbinom(x=1, size=19, prob=0.25, red=2, log=FALSE),as.numeric(drbinomAPA(x=1, size=19, prob=0.25, red=2, log=FALSE, precBits=53)))
  expect_equal(drbinom(x=2, size=19, prob=0.25, red=2, log=FALSE),as.numeric(drbinomAPA(x=2, size=19, prob=0.25, red=2, log=FALSE, precBits=53)))
  expect_equal(drbinom(x=3, size=19, prob=0.25, red=2, log=FALSE),as.numeric(drbinomAPA(x=3, size=19, prob=0.25, red=2, log=FALSE, precBits=53)))
  expect_equal(drbinom(x=4, size=19, prob=0.25, red=2, log=FALSE),as.numeric(drbinomAPA(x=4, size=19, prob=0.25, red=2, log=FALSE, precBits=53)))
  expect_equal(drbinom(x=5, size=19, prob=0.25, red=2, log=FALSE),as.numeric(drbinomAPA(x=5, size=19, prob=0.25, red=2, log=FALSE, precBits=53)))
  expect_equal(drbinom(x=6, size=19, prob=0.25, red=2, log=FALSE),as.numeric(drbinomAPA(x=6, size=19, prob=0.25, red=2, log=FALSE, precBits=53)))
  expect_equal(drbinom(x=7, size=19, prob=0.25, red=2, log=FALSE),as.numeric(drbinomAPA(x=7, size=19, prob=0.25, red=2, log=FALSE, precBits=53)))
  expect_equal(drbinom(x=8, size=19, prob=0.25, red=2, log=FALSE),as.numeric(drbinomAPA(x=8, size=19, prob=0.25, red=2, log=FALSE, precBits=53)))
  expect_equal(drbinom(x=9, size=19, prob=0.25, red=2, log=FALSE),as.numeric(drbinomAPA(x=9, size=19, prob=0.25, red=2, log=FALSE, precBits=53)))
  expect_equal(drbinom(x=10,size=19, prob=0.25, red=2, log=FALSE),as.numeric(drbinomAPA(x=10, size=19, prob=0.25, red=2, log=FALSE, precBits=53)))

  expect_equal(sum(drbinom(x=0:3, size=10, prob=0.25, red=3, log=FALSE)),as.numeric(sum(drbinomAPA(x=0:3, size=10, prob=0.25, red=3, log=FALSE, precBits=53))))
  expect_equal(drbinom(x=0, size=10, prob=0.75, red=3, log=FALSE),as.numeric(drbinomAPA(x=0, size=10, prob=0.75, red=3, log=FALSE, precBits=53)))
  expect_equal(drbinom(x=1, size=10, prob=0.75, red=3, log=FALSE),as.numeric(drbinomAPA(x=1, size=10, prob=0.75, red=3, log=FALSE, precBits=53)))
  expect_equal(drbinom(x=2, size=10, prob=0.75, red=3, log=FALSE),as.numeric(drbinomAPA(x=2, size=10, prob=0.75, red=3, log=FALSE, precBits=53)))
  expect_equal(drbinom(x=3, size=10, prob=0.75, red=3, log=FALSE),as.numeric(drbinomAPA(x=3, size=10, prob=0.75, red=3, log=FALSE, precBits=53)))

  expect_equal(sum(drbinom(x=0:3, size=10, prob=0.5, red=4, log=FALSE)),as.numeric(sum(drbinomAPA(x=0:3, size=10, prob=0.5, red=4, log=FALSE, precBits=53))))
  expect_equal(drbinom(x=0, size=10, prob=0.5, red=4, log=FALSE),as.numeric(drbinomAPA(x=0, size=10, prob=0.5, red=4, log=FALSE, precBits=53)))
  expect_equal(drbinom(x=1, size=10, prob=0.5, red=4, log=FALSE),as.numeric(drbinomAPA(x=1, size=10, prob=0.5, red=4, log=FALSE, precBits=53)))
  expect_equal(drbinom(x=2, size=10, prob=0.5, red=4, log=FALSE),as.numeric(drbinomAPA(x=2, size=10, prob=0.5, red=4, log=FALSE, precBits=53)))
  expect_equal(drbinom(x=3, size=10, prob=0.5, red=4, log=FALSE),as.numeric(drbinomAPA(x=3, size=10, prob=0.5, red=4, log=FALSE, precBits=53)))

})


test_that("drpois and drpoisAPA", {
  expect_equal(drpois(x=0:15, lambda=5, red=1, log=FALSE), dpois(x=0:15, lambda=5, log=FALSE))
  expect_equal(drpois(x=5, lambda=0:10, red=1, log=FALSE), dpois(x=5, lambda=0:10, log=FALSE))
  expect_equal(drpois(x=0:15, lambda=5, red=1, log=FALSE), as.numeric(drpoisAPA(x=0:15, lambda=5, red=1, log=FALSE)))

  # NAPA vs base R
  expect_equal(drpois(x=0, lambda=10, red=2, log=FALSE),sum(dpois(x=0, lambda=10, log=FALSE)))
  expect_equal(drpois(x=1, lambda=10, red=2, log=FALSE),sum(dpois(x=1:2, lambda=10, log=FALSE)))
  expect_equal(drpois(x=2, lambda=10, red=2, log=FALSE),sum(dpois(x=3:4, lambda=10, log=FALSE)))
  expect_equal(drpois(x=3, lambda=10, red=2, log=FALSE),sum(dpois(x=5:6, lambda=10, log=FALSE)))
  expect_equal(drpois(x=4, lambda=10, red=2, log=FALSE),sum(dpois(x=7:8, lambda=10, log=FALSE)))
  expect_equal(drpois(x=5, lambda=10, red=2, log=FALSE),sum(dpois(x=9:10, lambda=10, log=FALSE)))

  expect_equal(drpois(x=0, lambda=10, red=3, log=FALSE),sum(dpois(x=0:1, lambda=10, log=FALSE)))
  expect_equal(drpois(x=1, lambda=10, red=3, log=FALSE),sum(dpois(x=2:4, lambda=10, log=FALSE)))
  expect_equal(drpois(x=2, lambda=10, red=3, log=FALSE),sum(dpois(x=5:7, lambda=10, log=FALSE)))
  expect_equal(drpois(x=3, lambda=10, red=3, log=FALSE),sum(dpois(x=8:10, lambda=10, log=FALSE)))
  expect_equal(drpois(x=4, lambda=10, red=3, log=FALSE),sum(dpois(x=11:13, lambda=10, log=FALSE)))
  expect_equal(drpois(x=5, lambda=10, red=3, log=FALSE),sum(dpois(x=14:16, lambda=10, log=FALSE)))

  expect_equal(drpois(x=0, lambda=10, red=4, log=FALSE),sum(dpois(x=0:1, lambda=10, log=FALSE)))
  expect_equal(drpois(x=1, lambda=10, red=4, log=FALSE),sum(dpois(x=2:5, lambda=10, log=FALSE)))
  expect_equal(drpois(x=2, lambda=10, red=4, log=FALSE),sum(dpois(x=6:9, lambda=10, log=FALSE)))
  expect_equal(drpois(x=3, lambda=10, red=4, log=FALSE),sum(dpois(x=10:13, lambda=10, log=FALSE)))
  expect_equal(drpois(x=4, lambda=10, red=4, log=FALSE),sum(dpois(x=14:17, lambda=10, log=FALSE)))
  expect_equal(drpois(x=5, lambda=10, red=4, log=FALSE),sum(dpois(x=18:21, lambda=10, log=FALSE)))

  # APA vs NAPA
  expect_equal(drpois(x=0, lambda=10, red=2, log=FALSE),as.numeric(drpoisAPA(x=0, lambda=10, red=2, log=FALSE)))
  expect_equal(drpois(x=1, lambda=10, red=2, log=FALSE),as.numeric(drpoisAPA(x=1, lambda=10, red=2, log=FALSE)))
  expect_equal(drpois(x=2, lambda=10, red=2, log=FALSE),as.numeric(drpoisAPA(x=2, lambda=10, red=2, log=FALSE)))
  expect_equal(drpois(x=3, lambda=10, red=2, log=FALSE),as.numeric(drpoisAPA(x=3, lambda=10, red=2, log=FALSE)))
  expect_equal(drpois(x=4, lambda=10, red=2, log=FALSE),as.numeric(drpoisAPA(x=4, lambda=10, red=2, log=FALSE)))
  expect_equal(drpois(x=5, lambda=10, red=2, log=FALSE),as.numeric(drpoisAPA(x=5, lambda=10, red=2, log=FALSE)))

  expect_equal(drpois(x=0, lambda=10, red=3, log=FALSE),as.numeric(drpoisAPA(x=0, lambda=10, red=3, log=FALSE)))
  expect_equal(drpois(x=1, lambda=10, red=3, log=FALSE),as.numeric(drpoisAPA(x=1, lambda=10, red=3, log=FALSE)))
  expect_equal(drpois(x=2, lambda=10, red=3, log=FALSE),as.numeric(drpoisAPA(x=2, lambda=10, red=3, log=FALSE)))
  expect_equal(drpois(x=3, lambda=10, red=3, log=FALSE),as.numeric(drpoisAPA(x=3, lambda=10, red=3, log=FALSE)))
  expect_equal(drpois(x=4, lambda=10, red=3, log=FALSE),as.numeric(drpoisAPA(x=4, lambda=10, red=3, log=FALSE)))
  expect_equal(drpois(x=5, lambda=10, red=3, log=FALSE),as.numeric(drpoisAPA(x=5, lambda=10, red=3, log=FALSE)))

  expect_equal(drpois(x=0, lambda=10, red=4, log=FALSE),as.numeric(drpoisAPA(x=0, lambda=10, red=4, log=FALSE)))
  expect_equal(drpois(x=1, lambda=10, red=4, log=FALSE),as.numeric(drpoisAPA(x=1, lambda=10, red=4, log=FALSE)))
  expect_equal(drpois(x=2, lambda=10, red=4, log=FALSE),as.numeric(drpoisAPA(x=2, lambda=10, red=4, log=FALSE)))
  expect_equal(drpois(x=3, lambda=10, red=4, log=FALSE),as.numeric(drpoisAPA(x=3, lambda=10, red=4, log=FALSE)))
  expect_equal(drpois(x=4, lambda=10, red=4, log=FALSE),as.numeric(drpoisAPA(x=4, lambda=10, red=4, log=FALSE)))
  expect_equal(drpois(x=5, lambda=10, red=4, log=FALSE),as.numeric(drpoisAPA(x=5, lambda=10, red=4, log=FALSE)))

  expect_equal(drpois(x=0, lambda=10, red=5, log=FALSE),as.numeric(drpoisAPA(x=0, lambda=10, red=5, log=FALSE)))
  expect_equal(drpois(x=1, lambda=10, red=5, log=FALSE),as.numeric(drpoisAPA(x=1, lambda=10, red=5, log=FALSE)))
  expect_equal(drpois(x=2, lambda=10, red=5, log=FALSE),as.numeric(drpoisAPA(x=2, lambda=10, red=5, log=FALSE)))
  expect_equal(drpois(x=3, lambda=10, red=5, log=FALSE),as.numeric(drpoisAPA(x=3, lambda=10, red=5, log=FALSE)))
  expect_equal(drpois(x=4, lambda=10, red=5, log=FALSE),as.numeric(drpoisAPA(x=4, lambda=10, red=5, log=FALSE)))
  expect_equal(drpois(x=5, lambda=10, red=5, log=FALSE),as.numeric(drpoisAPA(x=5, lambda=10, red=5, log=FALSE)))

})
