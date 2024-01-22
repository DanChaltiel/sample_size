# Sample size for a Two-Sample Logrank Test
Dan Chaltiel

This section addresses methods to compute the total sample size (for a
1:1 randomization) for a two-sample logrank test.

## Sample size: Package `npsurvSS`

![](https://img.shields.io/badge/East-Validated-green.svg)

![](https://img.shields.io/badge/nQuery-Untested-blue.svg)

WIP: nQuery gives estimates lower than the ones calculated here.

The package `npsurvSS` allows to compute the sample size for a logrank
test.

Here is a simple wrapper around `npsurvSS::size_two_arm()`:

``` r
samplesize_surv = function(surv_ctl, surv_exp, surv_time, 
                           accr_time, follow_time, 
                           ratio=1, loss_scale=0,
                           alpha=0.05, sides=2, power=0.8){ 
  arm0 = npsurvSS::create_arm(size=1, 
                              accr_time=accr_time, follow_time=follow_time, 
                              surv_scale=-log(surv_exp)/surv_time, 
                              loss_scale=loss_scale)
  arm1 = npsurvSS::create_arm(size=ratio, 
                              accr_time=accr_time, follow_time=follow_time, 
                              surv_scale=-log(surv_ctl)/surv_time, 
                              loss_scale=loss_scale)
  npsurvSS::size_two_arm(arm0, arm1, alpha=alpha, sides=sides, power=power) |> 
    setNames(c("N_exp", "N_ctl", "N_tot", "E_exp", "E_ctl", "E_tot")) |> 
    ceiling()
}

#3y-survival of 61% versus 49%
samplesize_surv(surv_ctl=0.49, surv_exp=0.61, surv_time=3, 
                accr_time=3, follow_time=1.5, 
                alpha=0.05, sides=2, power=0.9) 
```

    N_exp N_ctl N_tot E_exp E_ctl E_tot 
      356   356   712   137   178   315 

``` r
#Same design using East v6.4: 
c(N_exp=355, N_ctl=354, N_tot=709, E_exp=137, E_ctl=176, E_tot=313)
```

    N_exp N_ctl N_tot E_exp E_ctl E_tot 
      355   354   709   137   176   313 

### Grid search

The following example considers a difference between a control arm with
a 3-year PFS of 50% and an experimental arm with a 3-year PFS of 70%.
The dropout was estimated exponential, with 1.5% patients missing after
3 years.

It uses `expand_grid()` to calculate the sample size for an abitrary
number of scenarios of accrual/follow-up time, alpha, and power.

``` r
library(tidyverse)
loss_scale = -log(0.985)/3

rslt = expand_grid(accrual=3:5, followup=2:3,
                   alpha=c(0.1, 0.2),
                   power=c(0.8, 0.9)) %>%
  rowwise() %>% 
  mutate(tbl = {
    samplesize_surv(surv_ctl=0.5, surv_exp=0.7, surv_time=3, loss_scale=loss_scale,
                    accr_time=accrual, follow_time=followup,
                    alpha=alpha, sides=2, power=power) |>
      as_tibble_row()
  }) %>% 
  unpack(tbl)
rslt
```

    # A tibble: 24 x 10
       accrual followup alpha power N_exp N_ctl N_tot E_exp E_ctl E_tot
         <int>    <int> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
     1       3        2   0.1   0.8    66    66   131    22    36    58
     2       3        2   0.1   0.9    91    91   181    31    49    80
     3       3        2   0.2   0.8    48    48    95    16    26    42
     4       3        2   0.2   0.9    70    70   139    24    38    61
     5       3        3   0.1   0.8    55    55   109    23    35    57
     6       3        3   0.1   0.9    76    76   151    31    48    79
     7       3        3   0.2   0.8    40    40    79    17    26    42
     8       3        3   0.2   0.9    58    58   116    24    37    61
     9       4        2   0.1   0.8    60    60   120    22    35    57
    10       4        2   0.1   0.9    83    83   166    31    49    79
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

In our previous example with 3y-survival of 61% versus 49%, we found a
total sample size of 709, so here we can find our 90% power back:

``` r
logrank_power(tref=3, surv_ctl=0.49, surv_exp=0.61, 
              n_ctl=709/2, n_exp=709/2, accrual=3, tmin=1.5,
              dropout=0, alpha=0.05, verbose=FALSE)
```

        Power 
    0.8957575 
