# Sample size for a Two-Sample Logrank Test
Dan Chaltiel

This section addresses methods to compute the total sample size (for a
1:1 randomization) for a two-sample logrank test.

## Sample size: Package `npsurvSS`

![](https://img.shields.io/badge/East-Untested-blue.svg)

![](https://img.shields.io/badge/nQuery-Untested-blue.svg)

The package `npsurvSS` allows to compute the sample size for a logrank
test.

During arms creation, the `size` parameter is irrelevant and any
value\>0 will yield the same results. Survival and dropout are modelized
using a Weibull distribution with default shape of 1 (exponential
distribution). See `?npsurvSS::create_arm` for more information.

``` r
library(npsurvSS)

arm0 = create_arm(size=1, accr_time=6, follow_time=12, 
                  surv_scale=0.05, loss_scale=0)
arm1 = create_arm(size=1, accr_time=6, follow_time=12, 
                  surv_scale=0.03, loss_scale=0)
size_two_arm(arm0, arm1, alpha=0.025, sides=1, power=0.8)
```

           n0        n1         n        d0        d1         d 
    137.66566 137.66566 275.33133  72.39288  49.76761 122.16049 

### Grid search

The following example considers a difference between a control arm with
a 3-year PFS of 50% and an experimental arm with a 3-year PFS of 70%.

It uses `expand_grid()` to calculate the sample size for an abitrary
number of scenarios of accrual/follow-up time, alpha, and power.

``` r
library(tidyverse)
pfs3y_ctl = 0.5
pfs3y_exp = 0.7
rslt = expand_grid(accrual=3:5, followup=2:3,
                   alpha=c(0.1, 0.2),
                   power=c(0.8, 0.9)) %>%
  rowwise() %>% 
  mutate(tbl = {
    arm0 = create_arm(size=1, accr_time=accrual, follow_time=followup, 
                      surv_scale=-log(pfs3y_ctl)/3, loss_scale=0.005)
    arm1 = create_arm(size=1, accr_time=accrual, follow_time=followup, 
                      surv_scale=-log(pfs3y_exp)/3, loss_scale=0.005)
    as_tibble_row(size_two_arm(arm0, arm1, alpha=alpha, sides=2, power=power))
  }) %>% 
  unpack(tbl)

rslt
```

    # A tibble: 24 x 10
       accrual followup alpha power    n0    n1     n    d0    d1     d
         <int>    <int> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
     1       3        2   0.1   0.8  65.1  65.1 130.   35.3  21.8  57.0
     2       3        2   0.1   0.9  90.4  90.4 181.   48.9  30.2  79.1
     3       3        2   0.2   0.8  47.4  47.4  94.7  25.6  15.8  41.5
     4       3        2   0.2   0.9  69.2  69.2 138.   37.5  23.1  60.6
     5       3        3   0.1   0.8  54.3  54.3 109.   34.4  22.1  56.5
     6       3        3   0.1   0.9  75.5  75.5 151.   47.8  30.7  78.5
     7       3        3   0.2   0.8  39.5  39.5  79.0  25.0  16.1  41.1
     8       3        3   0.2   0.9  57.8  57.8 116.   36.6  23.5  60.1
     9       4        2   0.1   0.8  59.6  59.6 119.   34.8  22.0  56.7
    10       4        2   0.1   0.9  82.7  82.7 165.   48.3  30.5  78.8
    # i 14 more rows

## Power calculation: package `Hmisc`

The power can be found using F. Harrell’s `Hmisc::cpower()` function.

However, the interface if not very easy to use, so I use the following
wrapper.

``` r
#' Power of Cox/log-rank Two-Sample Test 
#'
#' @description
#' Wrapper around [Hmisc::cpower()]
#'
#'
#' @param tref time of reference
#' @param surv_ctl,surv_exp survival at time `tref` in the control/experimental arm
#' @param n_ctl,n_exp sample size in the control/experimental arm
#' @param dropout proportion of non-compliants
#' @param alpha type I error probability (2-sided)
#' @param verbose whether to print informations
#' @param ... passed on to `Hmisc::cpower()`
#' 
#' @author Dan Chaltiel
logrank_power = function(tref, surv_ctl, surv_exp, 
                         n_ctl, n_exp, accrual, tmin,
                         dropout=0, alpha=0.05, verbose=TRUE, ...){
  mortality_c = 1-surv_ctl
  mortality_i = 1-surv_exp
  dropout=100*dropout
  r = 100*(1-mortality_i/mortality_c)
  Hmisc::cpower(tref, nc=n_ctl, ni=n_exp, mc=mortality_c, r=r, 
                accrual=accrual, tmin=tmin, noncomp.c=dropout, 
                noncomp.i=dropout, alpha=alpha, pr=verbose, ...)
}
```

In our previous example with scale 0.03 versus 0.05, we found a total
sample size of 276, so here we can find our 80% power back:

``` r
tref=10 #arbitraty, works for any value >0 
surv_ctl = exp(-0.03*tref) 
surv_exp = exp(-0.05*tref)

logrank_power(tref, surv_ctl, surv_exp, 
              n_ctl=276/2, n_exp=276/2, accrual=6, tmin=12,
              dropout=0, alpha=0.05, verbose=FALSE)
```

        Power 
    0.7931923 
