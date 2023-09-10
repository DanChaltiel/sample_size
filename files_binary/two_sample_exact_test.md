# Sample size for a Two-Sample Exact Test
Dan Chaltiel

This section addresses methods to compute the total sample size (for a
1:1 randomization) for a two-sample exact test.

## Custom function

![](https://img.shields.io/badge/East-Same-green.svg)

![](https://img.shields.io/badge/nQuery-Not%20tested-blue.svg)

The following function is a wrapper around
[`Exact::power.exact.test()`](https://www.rdocumentation.org/packages/Exact/versions/3.2/topics/power.exact.test)
and it has been validated against East on several cases.

``` r
#' Compute sample size for a Two-Sample Exact Test
#'
#' @description
#' Wrapper around [Exact::power.exact.test] to compute total 
#' sample size (for a 1:1 randomization)
#'
#' @author Dan Chaltiel
sample_size_extact_test = function(p1, p2, alpha=0.05, power=0.8,
                                   interval=c(2,500),
                                   alternative=c("two.sided","less","greater"),
                                   ...){
  f=function(x){
    x2=abs(round(x/2))
    p = Exact::power.exact.test(p1, p2, x2, x2, alpha=alpha, 
                                alternative=alternative, ...)$power
    if(p<power) return(.Machine$double.xmax)
    abs(p-power)
  }
  n = optimize(f, interval=interval)
  n = ceiling(n$minimum)
  x = Exact::power.exact.test(p1, p2, n/2, n/2, alpha=alpha, 
                              alternative=alternative, ...)
  if(x$power<power){
    m = c("Optimization failed for {.val p1={p1}}, {.val p2={p2}}, 
                    {.val alpha={alpha}}, {.val power={power}}.",
          i="Increase the interval upper limit (current={interval[2]}).")
    cli::cli_warn(m)
  }
  rtn = tibble::tibble(p1=x$`p1, p2`[1], p2=x$`p1, p2`[2], alpha=x$alpha,
                       power=x$power, n_total=sum(x$`n1, n2`))
  attr(rtn, "method") = x$method
  rtn
}
```

### One-line examples

The function outputs a dataframe with input probabilities and alpha, the
actual attained power, and the total sample size. Each arm will
therefore hold `n_total/2` patients.

``` r
sample_size_extact_test(p1=0.35, p2=0.55, power=0.9, alternative="less")
```

    # A tibble: 1 x 5
         p1    p2 alpha power n_total
      <dbl> <dbl> <dbl> <dbl>   <dbl>
    1  0.35  0.55  0.05 0.902     218

``` r
sample_size_extact_test(p1=0.35, p2=0.55, power=0.8, alternative="less")
```

    # A tibble: 1 x 5
         p1    p2 alpha power n_total
      <dbl> <dbl> <dbl> <dbl>   <dbl>
    1  0.35  0.55  0.05 0.807     158

``` r
# method="fisher"
sample_size_extact_test(p1=0.35, p2=0.55, power=0.8, alternative="less", method="fisher")
```

    # A tibble: 1 x 5
         p1    p2 alpha power n_total
      <dbl> <dbl> <dbl> <dbl>   <dbl>
    1  0.35  0.55  0.05 0.821     182

### Grid search

This example uses `expand_grid()` to calculate the sample size for an
abitrary number of scenarios.

It takes about 1 minute to compute.

``` r
library(tidyverse)
rslt = expand_grid(p=list(c(.35,.50), c(.35,.55), c(.60,.75), c(.60,.80)),
                    power=c(0.8, 0.9)) %>%
  mutate(tbl = map2(p, power, ~{
    sample_size_extact_test(p1=.x[1], p2=.x[2], power=.y, alternative="less")
  }, .progress=TRUE)) %>%
  select(tbl) %>% unnest(tbl) %>%
  arrange(round(power,1), p1)
#timesaveR::to_tribble(rslt) %>% cat
rslt
```

    # A tibble: 8 x 5
         p1    p2 alpha power n_total
      <dbl> <dbl> <dbl> <dbl>   <dbl>
    1  0.35  0.5   0.05 0.804     272
    2  0.35  0.55  0.05 0.807     158
    3  0.6   0.75  0.05 0.802     246
    4  0.6   0.8   0.05 0.804     136
    5  0.35  0.5   0.05 0.900     372
    6  0.35  0.55  0.05 0.902     218
    7  0.6   0.75  0.05 0.901     336
    8  0.6   0.8   0.05 0.902     180

### Warning

As this function uses optimization, it relies on an `interval` argument.
The default is `c(2, 500)`.

- If the computation is too long, your sample size must be high, so you
  should increase the lower limit.

- If the power is not reached, a warning will tell you to increase the
  upper limit.

In the following example, the optimization won’t go above 50 patients,
which is not enough to reach a power of 80%.

``` r
sample_size_extact_test(p1=0.35, p2=0.55, power=0.8, interval=c(2,50))
```

    Warning: Optimization failed for "p1=0.35", "p2=0.55", "alpha=0.05", "power=0.8".
    i Increase the interval upper limit (current=50).

    # A tibble: 1 x 5
         p1    p2 alpha power n_total
      <dbl> <dbl> <dbl> <dbl>   <dbl>
    1  0.35  0.55  0.05 0.260      50

## Base functions

![](https://img.shields.io/badge/East-Different-red.svg)

![](https://img.shields.io/badge/nQuery-Not%20tested-blue.svg)

Base R function `stats::power.prop.test()` gives results very close to
but a bit different from East’s.

However, it is much faster, as the following code is nearly
instantaneous. If you have a lot of scenarios to consider, this method
might be a better choice for starters.

``` r
library(tidyverse)
```

    -- Attaching core tidyverse packages ------------------------ tidyverse 2.0.0 --
    v dplyr     1.1.2     v readr     2.1.4
    v forcats   1.0.0     v stringr   1.5.0
    v ggplot2   3.4.2     v tibble    3.2.1
    v lubridate 1.9.2     v tidyr     1.3.0
    v purrr     1.0.1     
    -- Conflicts ------------------------------------------ tidyverse_conflicts() --
    x dplyr::filter() masks stats::filter()
    x dplyr::lag()    masks stats::lag()
    i Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
expand_grid(p=list(c(.35,.50), c(.35,.55), c(.60,.75), c(.60,.80)),
            power=c(0.8, 0.9)) %>%
  mutate(tbl = map2(p, power, ~{
    t = power.prop.test(p1=.x[1], p2=.x[2], power=.y, alternative="one.sided")
    tibble(p1=t$p1, p2=t$p2, alpha=t$sig.level, power=t$power, n_total=round(2*t$n))
  }, .progress=TRUE)) %>%
  select(tbl) %>% unnest(tbl) %>%
  arrange(round(power,1), p1)
```

    # A tibble: 8 x 5
         p1    p2 alpha power n_total
      <dbl> <dbl> <dbl> <dbl>   <dbl>
    1  0.35  0.5   0.05   0.8     266
    2  0.35  0.55  0.05   0.8     151
    3  0.6   0.75  0.05   0.8     239
    4  0.6   0.8   0.05   0.8     128
    5  0.35  0.5   0.05   0.9     368
    6  0.35  0.55  0.05   0.9     208
    7  0.6   0.75  0.05   0.9     330
    8  0.6   0.8   0.05   0.9     176
