---
title: "README"
output:
  html_document:
    keep_md: yes
---



# redNMix

_redNMix_ is an R package for fitting N-mixture models. The package has several unique features, allowing for:

1. grouped/coarse counts data  
2. parallel processing  
3. applications to large populations

The package can be installed from github:


```r
remotes::install_github("mrparker909/redNMix")
```


```r
library(redNMix)
```

# Generate Data

There are helper functions for quickly generating N-mixtures data. These are particularily useful for testing and for simulation studies.

## Closed Population

We generate a closed population with 3 sampling sites and 4 sampling occasions, with site abundance parameter lambda=5, and a probability of detection of 0.80.


```r
pop1 = redNMix::gen_Nmix_closed(num_sites = 3,
                                num_times = 4,
                                lambda    = 5,
                                pdet      = 0.8)
```

_pop1_ contains two matrices, the first is the true population size for each sampling:


```r
pop1$Ni
```

```
##      [,1] [,2] [,3] [,4]
## [1,]    6    6    6    6
## [2,]    8    8    8    8
## [3,]    6    6    6    6
```

while the second matrix are the observed counts:


```r
pop1$nit
```

```
##      [,1] [,2] [,3] [,4]
## [1,]    4    5    2    4
## [2,]    7    6    8    8
## [3,]    6    4    6    5
```

## Open Population

We generate an open population with 3 sampling sites and 4 sampling occasions, with initial abundance parameter lambda=5, recruitment parameter gamma=1, survival probability omega=0.25, and a probability of detection of 0.80.


```r
pop2 = redNMix::gen_Nmix_open(num_sites = 3,
                              num_times = 4,
                              lambda    = 5,
                              gamma     = 1,
                              omega     = 0.25,
                              pdet      = 0.8)
```

_pop2_ contains two matrices, the first is the true population size for each sampling:


```r
pop2$Ni
```

```
##      [,1] [,2] [,3] [,4]
## [1,]    5    0    1    0
## [2,]    4    6    3    1
## [3,]    4    2    0    1
```

while the second matrix are the observed counts:


```r
pop2$nit
```

```
##      [,1] [,2] [,3] [,4]
## [1,]    4    0    1    0
## [2,]    4    3    2    1
## [3,]    2    2    0    0
```

# Model Fitting

There are many options for model fitting, see the documentation in R for more examples (eg: `?redNMix::fit_red_Nmix_open`). 

## Closed Population Example


```r
out <- redNMix::fit_red_Nmix_closed(nit=pop1$nit, red=1, K=20)
```


```r
out
```

```
## $x
##           [,1]
## B_l_0 1.920851
## B_p_0 1.345819
## 
## $f
## [1] 21.07789
## 
## $grad
## [1] -2.807343e-07 -1.839293e-07
## 
## $inv_Hessian
##             [,1]       [,2]
## [1,]  0.06163461 -0.0613004
## [2,] -0.06130040  0.3705894
## 
## $steps
## [1] 32
## 
## $converged
## [1] TRUE
```

The optimization algorithm converged after 32 iterations. The value of the likelihood function at its maximum is 21.0778874. The parameter estimates are shown under `out$x`, and need to be transformed back to their usual units; lambda was log transformed, so we need to exponentiate it; pdet was logit transformed, so we need to inverse logit it.

The lambda estimate is:


```r
exp(out$x[1])
```

```
## [1] 6.826768
```

The probability of detection estimate is:


```r
plogis(out$x[2])
```

```
## [1] 0.7934453
```

The estimates are reasonably close to the true parameter values of 5 and 0.8.

# How to Cite

Parker, M.R.P. (2020). redNMix: An R package for N-mixtures models. R package version 1.0.1. https://github.com/mrparker909/redNMix



