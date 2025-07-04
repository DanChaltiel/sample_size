# A’Hern design for a 1-sample Test
Dan Chaltiel

This section addresses methods to compute the sample size for a
single-stage one-sample A’Hern design.

## Custom function

![](https://img.shields.io/badge/Validation-95%25-green.svg)

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
                             tolerance=0.00005) {
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
    6   0.6   0.7   0.9   197   130     0.0492      0.903
    7   0.6   0.8   0.8    33    25     0.0444      0.800
    8   0.6   0.8   0.9    45    33     0.0446      0.901

### Tolerance

As you can see, the real type I error rate can be slightly higher than
5%, as in the princeps article. Use the `tolerance` argument if you need
to be strict:

``` r
sample_size_ahern(p0=0.6, p1=0.7, alpha=0.05, power=0.9)
```

        n   r alpha_real power_real
    1 197 130 0.04915422  0.9031079

``` r
sample_size_ahern(p0=0.6, p1=0.7, alpha=0.05, power=0.9, tolerance=0)
```

        n   r alpha_real power_real
    1 197 130 0.04915422  0.9031079

## Validation

By scraping the article, we can validate that the calculation is
correct.

The article used MS Excel BINOMDIST function to compute the
probabilities, which admittedly failed for large N and was then replaced
by a normal approximation. As it never fails in R, this approximation
was not implemented. Most differences are very likely to come from this
fact.

The other difference is the `tolerance` argument. In the original
article, `alpha_real` and `power_real` seem to be rounded, but this is
not detailled. Differences in smaller N can be

``` r
ahern_table = read_lines("ahern_table.txt") %>% 
  as_tibble_col() %>% 
  separate_wider_delim (value, 
                        names=c("p0", "p1", "a0.05_p0.8", "a0.05_p0.9", 
                                "a0.01_p0.8", "a0.01_p0.9"),
                        too_few = "align_end",
                        delim=" ") %>% 
  filter(!is.na(p1)) %>% 
  fill(p0, .direction="down") %>% 
  pivot_longer(-c(p0, p1), names_to=c("alpha", "power"), names_pattern="a(.*)_p(.*)") %>% 
  separate(value, into=c("r", "n"), sep="=") %>% 
  mutate_all(as.numeric)

comparison = ahern_table %>% 
  rowwise() %>%
  mutate(
    ah = list(sample_size_ahern(p0=p0, p1=p1, alpha=alpha, power=power, n_min=n-15, n_max=n+15,
                                tolerance=0.00005) %>% suppressWarnings())
  ) %>% 
  unnest(ah, names_sep="_", keep_empty=TRUE) %>% 
  mutate(diff_n = ah_n-n, diff_r = ah_r-r)
```

Sample sizes are the same in \>95% cases:

``` r
comparison %>% 
  summarise(
    prop_correct__n = mean(!is.na(diff_n) & diff_n==0),
    prop_correct__r = mean(!is.na(diff_r) & diff_r==0),
    mean_diff__n = mean(diff_n[diff_n!=0], na.rm=TRUE),
    mean_diff__r = mean(diff_r[diff_r!=0], na.rm=TRUE)
  ) %>% 
  pivot_longer(everything(), names_sep="__", names_to=c("variable", ".value"))
```

    # A tibble: 2 × 3
      variable         n     r
      <chr>        <dbl> <dbl>
    1 prop_correct 0.955 0.962
    2 mean_diff    2.8   1.24 

Here are the cases with a difference:

``` r
comparison %>% 
  arrange(desc(pmax(diff_n, diff_r))) %>% 
  filter(diff_n!=0 | diff_r!=0 | is.na(diff_n) | is.na(diff_r)) %>% 
  print(n=Inf)
```

    # A tibble: 31 × 12
          p0    p1 alpha power     r     n  ah_n  ah_r ah_alpha_real ah_power_real
       <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <int> <int>         <dbl>         <dbl>
     1  0.3   0.35  0.01   0.9   377  1133  1147   381       0.00998         0.903
     2  0.35  0.4   0.01   0.9   462  1207  1216   465       0.0100          0.900
     3  0.75  0.8   0.01   0.9   713   910   919   720       0.00980         0.901
     4  0.4   0.45  0.01   0.9   548  1266  1274   551       0.00987         0.901
     5  0.55  0.6   0.01   0.8   578   984   991   582       0.00982         0.802
     6  0.3   0.35  0.05   0.9   247   752   758   249       0.0481          0.900
     7  0.8   0.85  0.05   0.8   300   359   365   305       0.0485          0.802
     8  0.8   0.85  0.01   0.9   633   759   764   637       0.00990         0.903
     9  0.25  0.3   0.01   0.9   291  1031  1035   292       0.0100          0.902
    10  0.25  0.35  0.01   0.9    85   270   274    86       0.0100          0.907
    11  0.7   0.8   0.01   0.9   190   247   251   193       0.00913         0.903
    12  0.2   0.5   0.01   0.9    12    30    33    13       0.00819         0.919
    13  0.3   0.4   0.05   0.8    52   141   144    53       0.0473          0.807
    14  0.3   0.4   0.01   0.9   107   293   296   108       0.00971         0.903
    15  0.3   0.55  0.01   0.9    24    51    53    25       0.00628         0.900
    16  0.4   0.55  0.01   0.9    71   142   144    72       0.00951         0.901
    17  0.4   0.6   0.01   0.9    44    82    84    45       0.00811         0.905
    18  0.45  0.55  0.01   0.8   133   253   255   134       0.00928         0.802
    19  0.5   0.55  0.05   0.8   330   618   620   331       0.0498          0.802
    20  0.5   0.55  0.01   0.8   540  1005  1007   541       0.00983         0.801
    21  0.5   0.6   0.01   0.8   144   250   252   145       0.00979         0.806
    22  0.5   0.6   0.01   0.9   183   323   325   184       0.00984         0.903
    23  0.05  0.15  0.05   0.9     8    76    77     8       0.0385          0.907
    24  0.05  0.35  0.05   0.8     3    11    12     3       0.0196          0.849
    25  0.15  0.4   0.05   0.8     7    21    22     7       0.0368          0.842
    26  0.5   0.55  0.01   0.9   693  1300  1301   693       0.00992         0.900
    27  0.55  0.6   0.01   0.9   743  1274  1275   743       0.00998         0.901
    28  0.65  0.7   0.01   0.9   793  1160  1155   789       0.00955         0.900
    29  0.6   0.65  0.01   0.9   779  1230  1223   774       0.00999         0.900
    30  0.7   0.75  0.01   0.9   772  1052  1042   764       0.0100          0.900
    31  0.45  0.5   0.01   0.9   624  1292    NA    NA      NA              NA    
    # ℹ 2 more variables: diff_n <dbl>, diff_r <dbl>

``` r
#differences in small N are due to rounding
sample_size_ahern(p0=0.05, p1=0.15, alpha=0.05, power=0.9,
                  tolerance=0.0005)
```

       n r alpha_real power_real
    1 76 8 0.03599767  0.8999097
