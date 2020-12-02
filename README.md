## Lusk *et al.* (2021)
This is a repository for the simulation code used in Lusk, S. C., C. R. Middaugh, and D. H. Ogle. 2021. Evaluating the performance of methods used to estimate growth parameters from subsampled age data. North American Journal of Fisheries Management XX:XXX-XXX. DOI: XXX.

The pertinent analysis files are described below.

## Repository Structure

* **code/**
    * **Run_Simulation.R** -- Generates a random "entire sample" of data, makes an "aged sample" and constructs an ALK, applies the four methods to expand the aged sample to the entire sample, fits the von Bertalanffy growth models using the five methods, and saves the results outlined in the manuscript to **Simulation_Run_Results.csv**. This script should be sourced as it uses a random number to replicate results and *only needs to be run once* as the resultant CSV file is read into analysis scripts described below for post processing. This script requires many minutes to complete (approaching two hours on my machine ... note this is why the results are saved out and used in subsequent analyses). More details are in comments in the script.
    * **Simulated_Data_Figure.R** -- Produces Figure 2 that illustrate examples of length and age data for each scenario (i.e., by numcap, SDCV, and species).
    * **Simulation_Results.R** -- Analyzes the simulation results in **Simulation_Run_Results.csv** to ultimately produce the results for Tables 2-5 (actual tables were produced outside of R) and Figures 2 and 4.
    * **Visualizing_Fits_Figure.R** -- Produces Figure 5 that illustrates examples of length and age data with model fits.
    * **Common.R** -- Common items and parameters (e.g., Table 1) for the simulations.
    * **EPfuns/**
        * **EPGrowth.R** -- creates two R functions -- `EPdata_prep()` and `EPrun()` used in **Run_Simulation.R** for the *empirical proportions* methods of Perreault *et al.* (2020).
        * **EP_likelihood.cpp** -- C++ code used by **EPGrowth.R** for the maximum likelihood method of Perreault *et al.* (per Perreault). Note that when **Common.R** is run that two new files -- **EP_likelihood.dll** and **EP_likelihood.o** -- will be created. If you get errors when running the empirical proportions portion of the code you may want to delete these two files so that they will be created fresh.

* **data/**
    * Ultimately a file called **Simulation_Run_Results.csv** will be created here by running **Run_Simulation.R**. This file, when created, is used in **Simulation_Results.R**.
    * **Simulation_Run_Results_PUBLISHED.csv** -- This is what the **Simulation_Run_Results.csv** file looked like at the time of publication. It is here for archival purposes, but also as a check against any new version of the file that is created by running **Simulation_Run_Results.csv**. Use the following to compare a newly created **Simulation_Run_Results.csv** with **Simulation_Run_Results_PUBLISHED.csv**. The last three lines should all say `TRUE` (thought `all()` and `identical()` may get caught up on minor differences due to machine infrastructure).

```
dnew <- read.csv("data/Simulation_Run_Results.csv")
dpub <- read.csv("data/Simulation_Run_Results_PUBLISHED.csv")
all(dnew==dpub)
identical(dnew,dpub)
compare::compareEqual(dnew,dpub)
```

* **docs/**
    * Ultimately this will contains figures and tables created by running the analysis scripts described above.
    * **Table2_raw_PUBLISHED.csv** -- This is what the **Table2_raw.csv** file looked like at the time of publication (note that this was used to create Tables 2-5 outside of R). It is here as a check against any new version of the file that is created by running the scripts described above. The last three lines should all say `TRUE` (thought `all()` and `identical()` may get caught up on minor differences due to machine infrastructure).

```
t2new <- read.csv("docs/Table2_Raw.csv")
t2pub <- read.csv("docs/Table2_Raw_PUBLISHED.csv")
all(t2new==t2pub)
identical(t2new,t2pub)
compare::compareEqual(t2new,t2pub)
```

## Suggestions for Implementation

* Use R 4.0.X or greater.
* Install (at least) the following packages from CRAN -- `FSA`, `ggnewscale`, `here`, `patchwork`, `rlang`, `tidyverse`, and `TMB`.
    * `TMB` uses C++ and may require some "finesse" to install. See [this page for help](https://github.com/kaskr/adcomp/wiki/Download).
    * You will also need `compare` if you want to compare data files and summary tables created with this code to what was published.
    * My full session info [is shown below.](#my-session-information)
* The `FSAsim` package is NOT available on CRAN, but can be installed from its GitHub repository [as described here](https://github.com/droglenc/FSAsim#installation).
* We suggest running the code in an RStudio Project environment. If you do not do this then you may need to change some directories (search for `here::here()` code in the scripts).
* By default there are 15 simulation scenarios with 500 iterations each, each of which is computationally heavy. To "test your setup" we suggest changing `numReps` on line 29 of **Run_Simulation.R** from 500 to 5 (or so) at first to make sure that all packages are installed, the file structure is correct, etc. After a successful run with these few iterations then change `numReps` back to 500 to reproduce the analyses of Lusk *et al.* (2021).

## My Session Information

```
> devtools::session_info()
- Session info -----------------------------------------------------------------
 setting  value                       
 version  R version 4.0.2 (2020-06-22)
 os       Windows 10 x64              
 system   x86_64, mingw32             
 ui       RStudio                     
 language (EN)                        
 collate  English_United States.1252  
 ctype    English_United States.1252  
 tz       America/Chicago             
 date     2020-12-01                  

- Packages ---------------------------------------------------------------------
 package     * version     date       lib source                          
 assertthat    0.2.1       2019-03-21 [1] CRAN (R 4.0.0)                  
 callr         3.5.1       2020-10-13 [1] CRAN (R 4.0.3)                  
 cli           2.2.0       2020-11-20 [1] CRAN (R 4.0.3)                  
 colorspace    2.0-0       2020-11-11 [1] CRAN (R 4.0.3)                  
 crayon        1.3.4       2017-09-16 [1] CRAN (R 4.0.0)                  
 desc          1.2.0       2018-05-01 [1] CRAN (R 4.0.0)                  
 devtools      2.3.2       2020-09-18 [1] CRAN (R 4.0.3)                  
 digest        0.6.25      2020-02-23 [1] CRAN (R 4.0.0)                  
 dplyr         1.0.2       2020-08-18 [1] CRAN (R 4.0.3)                  
 ellipsis      0.3.1       2020-05-15 [1] CRAN (R 4.0.0)                  
 fansi         0.4.1       2020-01-08 [1] CRAN (R 4.0.0)                  
 fs            1.5.0       2020-07-31 [1] CRAN (R 4.0.3)                  
 FSA           0.8.31.9000 2020-11-29 [1] local                           
 FSAsim        0.0.6.9999  2020-06-22 [1] Github (droglenc/FSAsim@edcd348)
 generics      0.1.0       2020-10-31 [1] CRAN (R 4.0.3)                  
 ggnewscale    0.4.3       2020-08-27 [1] CRAN (R 4.0.3)                  
 ggplot2       3.3.2       2020-06-19 [1] CRAN (R 4.0.0)                  
 ggpmisc       0.3.7       2020-11-09 [1] CRAN (R 4.0.3)                  
 glue          1.4.2       2020-08-27 [1] CRAN (R 4.0.3)                  
 gtable        0.3.0       2019-03-25 [1] CRAN (R 4.0.0)                  
 here          1.0.0       2020-11-15 [1] CRAN (R 4.0.3)                  
 lifecycle     0.2.0       2020-03-06 [1] CRAN (R 4.0.0)                  
 magrittr      2.0.1       2020-11-17 [1] CRAN (R 4.0.3)                  
 memoise       1.1.0       2017-04-21 [1] CRAN (R 4.0.0)                  
 munsell       0.5.0       2018-06-12 [1] CRAN (R 4.0.0)                  
 pillar        1.4.7       2020-11-20 [1] CRAN (R 4.0.3)                  
 pkgbuild      1.1.0       2020-07-13 [1] CRAN (R 4.0.2)                  
 pkgconfig     2.0.3       2019-09-22 [1] CRAN (R 4.0.0)                  
 pkgload       1.1.0       2020-05-29 [1] CRAN (R 4.0.0)                  
 plyr          1.8.6       2020-03-03 [1] CRAN (R 4.0.0)                  
 prettyunits   1.1.1       2020-01-24 [1] CRAN (R 4.0.0)                  
 processx      3.4.5       2020-11-30 [1] CRAN (R 4.0.3)                  
 ps            1.4.0       2020-10-07 [1] CRAN (R 4.0.3)                  
 purrr         0.3.4       2020-04-17 [1] CRAN (R 4.0.0)                  
 R6            2.5.0       2020-10-28 [1] CRAN (R 4.0.3)                  
 Rcpp          1.0.5       2020-07-06 [1] CRAN (R 4.0.2)                  
 remotes       2.2.0       2020-07-21 [1] CRAN (R 4.0.2)                  
 rlang         0.4.7       2020-07-09 [1] CRAN (R 4.0.0)                  
 rprojroot     2.0.2       2020-11-15 [1] CRAN (R 4.0.3)                  
 rstudioapi    0.13        2020-11-12 [1] CRAN (R 4.0.3)                  
 scales        1.1.1       2020-05-11 [1] CRAN (R 4.0.0)                  
 sessioninfo   1.1.1       2018-11-05 [1] CRAN (R 4.0.0)                  
 stringi       1.5.3       2020-09-09 [1] CRAN (R 4.0.3)                  
 stringr       1.4.0       2019-02-10 [1] CRAN (R 4.0.0)                  
 testthat      3.0.0       2020-10-31 [1] CRAN (R 4.0.3)                  
 tibble        3.0.4       2020-10-12 [1] CRAN (R 4.0.3)                  
 tidyr         1.1.2       2020-08-27 [1] CRAN (R 4.0.3)                  
 tidyselect    1.1.0       2020-05-11 [1] CRAN (R 4.0.0)                  
 usethis       1.6.3       2020-09-17 [1] CRAN (R 4.0.3)                  
 vctrs         0.3.5       2020-11-17 [1] CRAN (R 4.0.3)                  
 withr         2.3.0       2020-09-22 [1] CRAN (R 4.0.3)                  
 yaml          2.2.1       2020-02-01 [1] CRAN (R 4.0.0)                  

[1] C:/Apps/R/R-4.0.2/library
```
