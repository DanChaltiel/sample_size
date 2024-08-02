# Simon’s Two-Stage Design
Dan Chaltiel

## Online calculators

MD Anderson: <https://biostatistics.mdanderson.org/shinyapps/Simon2S/>

UNC Lineberger:
<http://cancer.unc.edu/biostatistics/program/ivanova/SimonsTwoStageDesign.aspx>

## R functions

- [`PhIIdesign::simon2stage()`](https://github.com/IDDI-BE/PhIIdesign):
  not available on CRAN but seems to yield results consistent with
  above-mentioned online calculators.

- [`mtdesign::obtainDesign()`](https://github.com/openpharma/mtdesign):
  available on CRAN but very slow

- [`clinfun::ph2simon()`](https://search.r-project.org/CRAN/refmans/clinfun/html/ph2simon.html):
  available on CRAN but the interface is odd (returns a list)

- [`OneArmPhaseTwoStudy::setupSimon()`](https://cran.r-project.org/web/packages/OneArmPhaseTwoStudy/index.html):
  archived from CRAN on 2021-04-07 (which is usually not a good sign)

## Grid search

![](https://img.shields.io/badge/East-Not%20tested-blue.svg)

![](https://img.shields.io/badge/nQuery-Not%20tested-blue.svg)

The following function is a wrapper around
[`PhIIdesign::simon2stage()`](https://github.com/IDDI-BE/PhIIdesign). I
has not been tested thoroughly yet, but as mentionned the results seems
OK. Final accepted design should probably be checked with a validated
tool.

The cache system uses `memoise`, so that any already computed design
doesn’t need to be computed again.

``` r
library(tidyverse)


#' Grid search for Simon's Two-Stage Design
#'
#' @description
#' Wrapper around [PhIIdesign::simon2stage]
#'
#' @author Dan Chaltiel
grid_search_simon = function(p0, p1, alpha, beta, N_min=25, N_max=50, 
                             verbose=TRUE, cache="cache/simon"){
  simon2stage = PhIIdesign::simon2stage
  simon2stage = purrr::possibly(simon2stage, quiet=TRUE, 
                                otherwise=tibble(MIN="Impossible"))
  simon2stage = memoise::memoise(simon2stage, cache=cachem::cache_disk(cache))
  pc = scales::percent
  f = function(p0, p1, alpha, beta){
    if(verbose) cat(p0, p1, alpha, beta, "\n")
    simon2stage(p0=p0, pa=p1, alpha=alpha, beta=beta, N_min=N_min, N_max=N_max)
  }
  
  rslt = expand_grid(p0=p0, p1=p1, alpha=alpha, beta=beta) %>% 
    rowwise() %>% 
    mutate(simon = list(f(p0=p0, p1=p1, alpha=alpha, beta=beta))) %>% 
    ungroup()
  # browser()
  rslt$simon %>% 
    list_rbind() %>% 
    as_tibble() %>% 
    unite("type", MIN, OPT, ADMISS, sep="") %>%
    transmute(
      design_id=row_number(),
      p0=pc(p0), p1=pc(pa), n=N, n1, r1, r2, 
      alpha=pc(alpha, 1), power=pc(1-beta, 1), 
      type
    )
}
```

### Reporting

The following function uses a single row of the output of the previous
one and give you the protocol section as per the MD Anderson template.

``` r
#' Report a Simon's Two-Stage Design
#'
#' @description
#' Writes the protocol section
#'
#' @source https://biostatistics.mdanderson.org/shinyapps/Simon2S/ 
#' @author Dan Chaltiel
report_simon = function(x){
  stopifnot(is.data.frame(x))
  stopifnot(nrow(x)==1)
  glue::glue("Simon's {x$type} two-stage design (Simon, 1989) will be used for conducting the trial. The null hypothesis is that the true response rate is {x$p0}. and the alternative hypothesis is that the true response rate is {x$p1}. The trial is carried out in two stages. In stage I, a total number of {x$n1} patients is accrued. If there are {x$r1} or fewer responses among these {x$n1} patients, the study will be early stopped. Otherwise, additional {x$n - x$n1} patients will be accrued in stage II, resulting in a total number sample size of {x$n}. If there are {x$r2+1} or more responses among these {x$n} patients, we reject the null hypothesis and claim that the treatment is promising. The design controls the type I error rate at {x$alpha} and yields the power of {x$power}.")
}
```

### Example

``` r
install.packages("tidyverse")
install.packages("memoise")
remotes::install_github("IDDI-BE/PhIIdesign")
```

``` r
#compute all possible designs
gs = grid_search_simon(alpha = c(0.05, 0.1), beta = c(0.2, 0.1),
                       p0 = seq(0.1, 0.2, by=0.05), p1 = seq(0.3, 0.4, by=0.05), 
                       verbose=FALSE)

#choose the right one (for instance here the #75)
result = gs %>% filter(design_id==75)

#write your protocol
report_simon(result)
```

    Simon's Minimax two-stage design (Simon, 1989) will be used for conducting the trial. The null hypothesis is that the true response rate is 20%. and the alternative hypothesis is that the true response rate is 35%. The trial is carried out in two stages. In stage I, a total number of 22 patients is accrued. If there are 4 or fewer responses among these 22 patients, the study will be early stopped. Otherwise, additional 19 patients will be accrued in stage II, resulting in a total number sample size of 41. If there are 12 or more responses among these 41 patients, we reject the null hypothesis and claim that the treatment is promising. The design controls the type I error rate at 10% and yields the power of 80%.
