# A’Hern design for a 1-sample Test
Dan Chaltiel

This section addresses methods to compute the sample size for a
single-stage one-sample A’Hern design.

## Custom function

![](https://img.shields.io/badge/Validation-TODO-blue.svg)

The following function is an implementation of the [A’Hern
design](https://stat.ethz.ch/education/semesters/as2012/bio/ahernSampleSize.pdf).

``` r
#' Compute exact single-stage sample size using A'Hern's design
#'
#' This function computes the optimal sample size `n` and decision threshold `r`
#' for a single-stage phase II design based on A'Hern (2001), using the exact
#' binomial distribution. 
#'
#' @param p0 Numeric. The unacceptable response rate under the null hypothesis (H0).
#' @param p1 Numeric. The desirable minimal response rate under the alternative hypothesis (H1).
#' @param alpha Numeric. One-sided type I error rate. 
#' @param power Numeric. Desired statistical power (1 - type II error). 
#' @param n_min,n_max Integer. Minimum & maximum sample size to consider. 
#' @param tolerance Numeric. Acceptable numerical slack for comparison of alpha and power. Default is 0.001.
#'
#' @references
#' A’Hern RP. Sample size tables for exact single-stage phase II designs. Stat Med. 2001;20(6):859–866.
#' https://stat.ethz.ch/education/semesters/as2012/bio/ahernSampleSize.pdf
sample_size_ahern = function(p0, p1, alpha=0.05, power=0.8, 
                             n_min=1, n_max=200, 
                             tolerance=0.001) {
  stopifnot(p0<p1)
  rtn = tidyr::expand_grid(n = n_min:n_max,
                    r = 0:n_max) |>
    dplyr::filter(r <= n) |>
    dplyr::mutate(
      alpha_real = 1 - pbinom(r - 1, size = n, prob = p0),
      power_real = 1 - pbinom(r - 1, size = n, prob = p1),
    ) |>
    dplyr::filter(alpha_real <= alpha+tolerance, 
                  power_real >= power-tolerance) |>
    dplyr::slice_min(r, by=n) |>
    dplyr::slice_min(n) |>
    as.data.frame()
  if(nrow(rtn)==0) {
    warning("Increase `n_max` ", 
            "p0=", p0, ", p1=",p1, ", alpha=",alpha, ", power=",power)
  }
  rtn
}
```

### One-line example

The function outputs a dataframe with input probabilities and alpha, the
actual attained power, and the total sample size. Each arm will
therefore hold `n_total/2` patients.

``` r
sample_size_ahern(p0=0.35, p1=0.55, alpha=0.05, power=0.8)
```

       n  r alpha_real power_real
    1 41 20 0.04806186  0.8309021

``` r
sample_size_ahern(p0=0.15, p1=0.3, alpha=0.1, power=0.9)
```

       n  r alpha_real power_real
    1 53 12 0.09066874  0.9094409

### Grid search

This example uses `expand_grid()` to calculate the sample size for an
abitrary number of scenarios.

``` r
library(tidyverse)
rslt = expand_grid(p0=c(0.5,0.6),
                   p1=c(0.7,0.8),
                   power=c(0.8, 0.9)
) %>%
  filter(p0<p1) %>% 
  rowwise() %>% 
  mutate(tbl = list(sample_size_ahern(p0=p0, p1=p1, power=power))) %>% 
  unnest(tbl)
rslt
```

    # A tibble: 8 × 7
         p0    p1 power     n     r alpha_real power_real
      <dbl> <dbl> <dbl> <int> <int>      <dbl>      <dbl>
    1   0.5   0.7   0.8    37    24     0.0494      0.807
    2   0.5   0.7   0.9    53    33     0.0492      0.914
    3   0.5   0.8   0.8    18    13     0.0481      0.867
    4   0.5   0.8   0.9    23    16     0.0466      0.928
    5   0.6   0.7   0.8   143    96     0.0477      0.800
    6   0.6   0.7   0.9   194   128     0.0509      0.902
    7   0.6   0.8   0.8    33    25     0.0444      0.800
    8   0.6   0.8   0.9    45    33     0.0446      0.901

### Tolerance

As you can see, the real type I error rate can be slightly higher than
5%. Use the `tolerance` argument if you need to be strict:

``` r
sample_size_ahern(p0=0.6, p1=0.7, alpha=0.05, power=0.9)
```

        n   r alpha_real power_real
    1 194 128 0.05089715  0.9021587

``` r
sample_size_ahern(p0=0.6, p1=0.7, alpha=0.05, power=0.9, tolerance=0)
```

        n   r alpha_real power_real
    1 197 130 0.04915422  0.9031079
