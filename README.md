# Warming increases the cost of growth in a model vertebrate

This repository contains code and data needed to reproduce the article:

**Barneche DR, Jahn M, Seebacher F**, Warming increases the cost of growth in a model vertebrate. *Functional Ecology* (in press).

## Instructions

All analyses were done in `R`. To compile the paper, including figures and tables we use the [remake](https://github.com/richfitz/remake) package for R. You can install remake using the `devtools` package:

```r
devtools::install_github("richfitz/remake", dependencies = TRUE)
```
(run `install.packages("devtools")` to install devtools if needed.)

The `remake` package also depends on `storr`, install it like this:
```r
devtools::install_github("richfitz/storr", dependencies = TRUE)
```

Next you need to open an R session with working directory set to the root of the project.

We use a number of packages, missing packages can be easily installed by remake:

```r
remake::install_missing_packages()
```

Then, to generate all figures, analyses, and tables, simply run:

```r
remake::make()
```

All output will be automatically placed in a directory called `output` (it is going to be automatically created for you).

Also notice that all the combined Bayesian models in this paper will take a several hours (up to a day) to run on a regular computer.

If you find remake confusing and prefer to run plain R, you can use remake to build a script `build.R` that produces a given output, e.g.

```r
remake::make_script(filename = 'build.R')
```

### This paper was produced using the following software and associated packages:
```
R version 3.5.3 (2019-03-11)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Mojave 10.14.3

Matrix products: default
BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib

locale:
[1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] rstan_2.18.2       StanHeaders_2.18.1 ggplot2_3.1.0      brms_2.8.0         Rcpp_1.0.1         plyr_1.8.4         rmarkdown_1.12     LoLinR_0.0.0.9000 

loaded via a namespace (and not attached):
 [1] Brobdingnag_1.2-6    gtools_3.8.1         threejs_0.3.1        shiny_1.2.0          assertthat_0.2.1     stats4_3.5.3         backports_1.1.3      pillar_1.3.1        
 [9] lattice_0.20-38      glue_1.3.1           digest_0.6.18        promises_1.0.1       colorspace_1.4-1     htmltools_0.3.6      httpuv_1.5.0         Matrix_1.2-15       
[17] dygraphs_1.1.1.6     pkgconfig_2.0.2      purrr_0.3.2          xtable_1.8-3         mvtnorm_1.0-10       scales_1.0.0         processx_3.3.0       later_0.8.0         
[25] tibble_2.1.1         bayesplot_1.6.0      DT_0.5               withr_2.1.2          shinyjs_1.0          lazyeval_0.2.2       cli_1.1.0            magrittr_1.5        
[33] crayon_1.3.4         mime_0.6             evaluate_0.13        ps_1.3.0             nlme_3.1-137         xts_0.11-2           pkgbuild_1.0.3       colourpicker_1.0    
[41] prettyunits_1.0.2    rsconnect_0.8.13     tools_3.5.3          loo_2.1.0            matrixStats_0.54.0   stringr_1.4.0        munsell_0.5.0        callr_3.2.0         
[49] compiler_3.5.3       rlang_0.3.3          grid_3.5.3           ggridges_0.5.1       htmlwidgets_1.3      crosstalk_1.0.0      igraph_1.2.4         miniUI_0.1.1.1      
[57] base64enc_0.1-3      gtable_0.3.0         inline_0.3.15        abind_1.4-5          markdown_0.9         reshape2_1.4.3       R6_2.4.0             gridExtra_2.3       
[65] rstantools_1.5.1     zoo_1.8-5            knitr_1.22           bridgesampling_0.6-0 dplyr_0.8.0.1        shinystan_2.5.0      shinythemes_1.1.2    stringi_1.4.3       
[73] parallel_3.5.3       tidyselect_0.2.5     xfun_0.6             coda_0.19-2          lmtest_0.9-36       
```

### How to download this project for people not familiar with GitHub:  
* on the project main page on GitHub, click on the green button `clone or download` and then click on `Download ZIP`  

## Bug reporting
* Please [report any issues or bugs](https://github.com/dbarneche/zebrafishCostOfGrowth/issues).
