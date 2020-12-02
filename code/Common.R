# This contains common items and parameters for the simulations. Some of these
# are used in the figure making scripts, Thus put into this common file.

#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=
# Prepare Environment ----
#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=
## Clear environment
rm(list=ls(all=TRUE))
## Load packages
library(FSAsim)
library(FSA)
library(tidyverse)
options(dplyr.summarise.inform=FALSE)
library(TMB)
## Set random number seed ... for repeatability
set.seed(63433)

#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=
# Simulation Parameters ----
#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## These may be set by the user
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## Groups to use
groups <- c("Black Crappie","Largemouth Bass","Striped Bass")
numGroups <- length(groups)

## BKC from Jackson and Hurley 2005. Journal of Freshwater Ecology 20:461-467
##    (just made up age distribution based on max age)
## LMB from Fernando et al. 2014. Journal of the Southeastern Association of
##    Fish and Wildlife Agencies 1:75-82 (age structure listed in Table 6)
## STB from Wilson et al. 2013. Biology and Management of Inland Striped Bass
##    and Hybrid Striped Bass (age structure from Lake Ouachita population)
Linfs <- c(286,461,848)
Ks <- c(0.443,0.34,0.31)
t0s <- c(0.106,-0.78,-0.56)
names(Linfs) <- names(Ks) <- names(t0s) <- groups

## Age structure for the groups
ads <- matrix(c(10 , 30,35,25,10, 5,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                144, 90,44,89,28,10, 1, 2,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                1  ,  8, 8,39,18,72,38,23,11, 7, 1, 2, 3, 3, 1, 1,3),
              byrow=FALSE,ncol=numGroups)
colnames(ads) <- groups
rownames(ads) <- 1:nrow(ads)
ads # just to check

## Creates percent table for table in manuscript
ads2 <- ads
ads2[is.na(ads2)] <- 0
round(prop.table(ads2,margin=2)*100,1)

predAges <- matrix(c(2,4,6,
                     3,6,8,
                     4,9,15),
                   byrow=FALSE,ncol=numGroups)

colnames(predAges) <- groups
rownames(predAges) <- 1:nrow(predAges)
predAges # just to check

## Number captured during "sampling"
numCap <- 500
## Variation parameter for simulated capture data
CVs <- c(0.05,0.10,0.15)
## Width of length bin and number to sample per bin when making the ALK
binSize <- 25
nPerBin <- c(5,10,15)

## Model fitting names
mdls <- c("Full Sample","Aged-Only","Weighted Mean from ALK","Assigned Ages",
          "Chih Reweight","Perrault's EP")
mdls.abb <- c("Full","AO","WM","AA","RW","EP")

## Create a typical von B growth function
vbgf <- FSA::vbFuns()

#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=
# Helper Function ----
#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=
vbDataGen <- function(n,minAge,maxAge,ageDist,Linf,K,t0,CV,binsize) {
  dplyr::tibble(id=1:n,
                ageCap=sample(minAge:maxAge,prob=ageDist/sum(ageDist),
                              size=n,replace=TRUE),
                lenCap=vbgf(ageCap,Linf,K,t0)) %>%
    dplyr::mutate(lenCap=lenCap+rnorm(n,mean=0,sd=CV*lenCap),
                  lenCap=round(lenCap,0)) %>%
    dplyr::filter(lenCap>0) %>%
    dplyr::mutate(lcat25=FSA::lencat(lenCap,w=binSize)) %>%
    as.data.frame()
}

## Configure likelihood (TMB) for EP method
### Compile C++ code for likeilhood into a dll file for TMB
TMB::compile("code/EPfuns/EP_likelihood.cpp")
### Loads the dll created above.
dyn.load("code/EPfuns/EP_likelihood")
### Load the EPdata_prep() and EPrun() helper functions
source("code/EPfuns/EPGrowth.R")

#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=
# GGPlot Theme ----
#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=
lblrCV <- function(string) paste("CV:",string)
lblrn <- function(string) paste("n:",string)
lblrGrp <- function(string) stringr::str_replace(string," ","\n")

theme_ALK <- theme_bw() +
  theme(legend.position="none",
        # remove minor grid in each plot
        panel.grid.minor=element_blank(),
        # remove major grid in each plot
        panel.grid.major=element_blank(),
        # change spacing between rows of facets
        panel.spacing.y=unit(2,"mm"),
        # change margins around facet labels
        strip.text.x=element_text(margin=margin(1,0,1,0)),
        strip.text.y=element_text(margin=margin(0,2,0,2)),
        # Remove box around facet labels
        strip.background=element_rect(color=NA),
        # Make tick marks skinnier
        axis.ticks=element_line(size=0.25),
        # Make axis labels black (rather than gray) and control size
        axis.text=element_text(size=10,color="black"),
        # Control size of axis labels
        axis.title = element_text(size=12,face="bold")
  )
