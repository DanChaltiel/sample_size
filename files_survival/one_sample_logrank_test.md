# Sample size for a One-Sample Logrank Test
Dan Chaltiel

This section addresses methods to compute the total sample size (for a
1:1 randomization) for a one-sample logrank test.

## Package `OneArm2stage`

![](https://img.shields.io/badge/East-Untested-blue.svg)

![](https://img.shields.io/badge/nQuery-Untested-blue.svg)

The function `OneArm2stage::Optimal.KJ()` allows to compute the sample
size for a one-sample logrank test.

In the following example, we consider an exponential distribution
(`dist="WB", shape=1`) with a theoric survival of 25% at t=12
(`S0=0.25, x0=12`), a minimal followup time of `tf=12`, an accrual rate
of 6 patient per time unit, the hazard ratio to be proven being
`hr=0.7`.

``` r
x = OneArm2stage::Optimal.KJ(dist="WB", shape=1, S0=0.25, x0=12, 
                             tf=12, rate=6, 
                             hr=0.7, alpha=0.05, beta=0.1)
```

    iter=1&EnH0=10000& Dc1/nbpt=2
    iter=2&EnH0=85.92& Dc1/nbpt=0.18182
    iter=3&EnH0=85.92& Dc1/nbpt=0.13228
    iter=4&EnH0=85.4& Dc1/nbpt=0.10498
    iter=5&EnH0=85.4& Dc1/nbpt=0.07691
    iter=6&EnH0=85& Dc1/nbpt=0.06104
    iter=7&EnH0=84.84& Dc1/nbpt=0.04498
    iter=8&EnH0=84.62& Dc1/nbpt=0.0357
    iter=9&EnH0=84.62& Dc1/nbpt=0.02645
    iter=10&EnH0=84.62& Dc1/nbpt=0.02099
    iter=11&EnH0=84.61& Dc1/nbpt=0.01562
    iter=12&EnH0=84.61& Dc1/nbpt=0.01239
    iter=13&EnH0=84.53& Dc1/nbpt=0.00926
    iter=14&EnH0=84.53& Dc1/nbpt=0.00735
    iter=15&EnH0=84.52& Dc1/nbpt=0.00551
    iter=16&EnH0=84.52& Dc1/nbpt=0.00437
    iter=17&EnH0=84.52& Dc1/nbpt=0.00329
    iter=18&EnH0=84.52& Dc1/nbpt=0.00261
    iter=19&EnH0=84.5& Dc1/nbpt=0.00197
    iter=20&EnH0=84.5& Dc1/nbpt=0.00156
    iter=21&EnH0=84.5& Dc1/nbpt=0.00118

The output contains 3 dataframes: the input parameters, the sample size
for a one-stage design, and the sample size for a two-stage design.

The one-stage design contains the sample size `nsingle`, the accrual
time `tasingle`, and the critical value `csingle`:

``` r
x$Single_stage
```

      nsingle tasingle  csingle
    1      90       15 1.644854

The two-stage design contains the sample size for the interim and final
analyses `n1` and `n`, the critical values `c1` and `c`, the interim
analysis time `t1`, the maximum total study length `MTSL`, the expected
sample size under H0 `ES`, and the probability or early stopping underh
H0 `PS`:

``` r
x$Two_stage
```

      n1      c1  n      c      t1 MTSL      ES     PS
    1 67 -0.4627 93 1.6342 11.0983 27.5 84.5012 0.3218

See `help(OneArm2stage::Optimal.KJ)` for more details.

### Grid search

The following example uses `expand_grid()` to calculate the sample size
for an abitrary number of scenarios.

However,

- `OneArm2stage::Optimal.KJ()` is very verbose (all iterations are
  printed) so we will apply `purrr::quietly()` to mute it

- `OneArm2stage::Optimal.KJ()` is also very long (each call takes about
  a minute) so we will apply `memoise::memoise()` to cache every result
  on the disk.

``` r
library(tidyverse)

quiet_kj = purrr::quietly(OneArm2stage::Optimal.KJ)
compute_kj = memoise::memoise(function(...) quiet_kj(...)$result, 
                              cache=cachem::cache_disk("../cache/OneArm2stage"))

x = expand_grid(accrual=c(5,7), 
                followup=c(6,12), 
                alpha=c(0.05, 0.1), 
                beta=c(0.2,0.3)) %>% 
  rowwise() %>% 
  mutate(tbl = {
    # cli::cat_line(cur_group_id(), "/", n_groups(.))
    d = compute_kj(dist="WB", shape=1, S0=0.30, x0=12, hr=0.7, tf=followup,
                   rate=accrual, alpha=alpha, beta=beta)
    d$Two_stage
  }) %>% 
  unpack(tbl)

x
```

    # A tibble: 16 x 12
       accrual followup alpha  beta    n1     c1     n     c    t1  MTSL    ES    PS
         <dbl>    <dbl> <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
     1       5        6  0.05   0.2    54 -0.231    85  1.62 10.7   23    72.1 0.409
     2       5        6  0.05   0.3    42 -0.233    69  1.61  8.38  19.8  57.9 0.408
     3       5        6  0.1    0.2    44 -0.533    64  1.26  8.73  18.8  58.0 0.297
     4       5        6  0.1    0.3    31 -0.650    50  1.25  6.00  16    44.8 0.258
     5       5       12  0.05   0.2    48 -0.538    71  1.63  9.40  26.2  63.9 0.295
     6       5       12  0.05   0.3    34 -0.620    56  1.62  6.73  23.2  50.0 0.268
     7       5       12  0.1    0.2    35 -0.905    53  1.27  6.90  22.6  49.6 0.183
     8       5       12  0.1    0.3    23 -1.16     39  1.27  4.59  19.8  37.0 0.124
     9       7        6  0.05   0.2    58 -0.333    90  1.62  8.17  18.9  77.9 0.370
    10       7        6  0.05   0.3    44 -0.335    74  1.61  6.20  16.6  62.7 0.369
    11       7        6  0.1    0.2    45 -0.698    68  1.26  6.35  15.7  62.3 0.243
    12       7        6  0.1    0.3    33 -0.724    53  1.25  4.60  13.6  48.1 0.235
    13       7       12  0.05   0.2    48 -0.714    73  1.63  6.79  22.4  67.0 0.238
    14       7       12  0.05   0.3    35 -0.689    58  1.62  4.97  20.3  52.3 0.245
    15       7       12  0.1    0.2    36 -1.09     54  1.27  5.09  19.7  51.5 0.137
    16       7       12  0.1    0.3    24 -1.27     40  1.27  3.31  17.7  38.3 0.102

NB: Using RStudio’s background jobs has been very useful for this
project!

### References

Wu, J, Chen L, Wei J, Weiss H, Chauhan A. (2020). Two-stage phase II
survival trial design. Pharmaceutical Statistics. 2020;19:214-229.
https://doi.org/10.1002/pst.1983
