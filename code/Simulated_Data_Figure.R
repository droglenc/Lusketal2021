# This script produces the plot that illustrate examples of length and age data
# for each scenario (i.e., by numcap, SDCV, and group). This is Figure 2.

#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=
# Prepare Environment ----
#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=
## Clear environment, memory, and console ... NOTE THIS!!!
rm(list=ls(all.names=TRUE))
gc()
cat("\014")

## Set working directory ... using here() requires an RStudio project
setwd(here::here())

## Load common items
source("code/Common.R")
library(patchwork)

#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=
# Simulation Results ----
#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=
## Initialize the results data.frame
df <- NULL
## Loop through groups, numPerBin, and SDs ... creating length-age data for
## each combination and saving in results data.frame
for (w in groups) {
  for (j in CVs) {
    ad <- ads[,w]
    ad <- ad[!is.na(ad)]
    minAge <- min(as.numeric(names(ad)))
    maxAge <- max(as.numeric(names(ad)))
    tmp <- vbDataGen(numCap,minAge,maxAge,ad,Linfs[[w]],Ks[[w]],t0s[[w]],
                     j,binsize) %>%
      dplyr::select(id,lenCap,ageCap) %>%
      dplyr::mutate(group=w,CV=j)
    df <- rbind(df,tmp)
  }
}
## This is a trick to order the facets in ggplot below.
df <- dplyr::mutate(df,group=factor(group,levels=groups),CV=factor(CV))

#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=
# Create Figure ----
#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=
simDataPlot <- function(d,xlbl="Age-at-Capture",ylbl="Total Length (mm)") {
  ggplot(data=d,aes(x=ageCap,y=lenCap)) +
    geom_point(shape=1,alpha=0.25,size=1) +
    facet_grid(group~CV,labeller=labeller(CV=lblrCV,group=lblrGrp)) +
    scale_y_continuous(name=ylbl,
                       limits=c(0,NA),expand=expansion(mult=c(0,0.05))) +
    scale_x_continuous(name=xlbl,
                       limits=c(0,NA),expand=expansion(mult=c(0,0.1))) +
    theme_ALK +
    # change spacing between columns of facets
    theme(panel.spacing.x=unit(3,"mm"),
          panel.background=element_rect(fill=NA))
}

BKC <- simDataPlot(filter(df,group=="Black Crappie"),xlbl=NULL,ylbl=NULL) +
  theme(plot.margin=unit(c(0,0,0,0),"mm"))
LMB <- simDataPlot(filter(df,group=="Largemouth Bass"),xlbl=NULL,ylbl=NULL) +
  theme(plot.margin=unit(c(0,0,0,0),"mm"),
        strip.text.x=element_text(color="white"),
        strip.background.x=element_rect(color=NA,fill="white"))
SPB <- simDataPlot(filter(df,group=="Striped Bass"),xlbl=NULL,ylbl=NULL) +
  theme(plot.margin=unit(c(0,0,0,0),"mm"),
        strip.text.x=element_text(color="white"),
        strip.background.x=element_rect(color=NA,fill="white"))

specs <- BKC / LMB / SPB

p <- wrap_elements(grid::textGrob("Total Length (mm)",rot=90,
                                  gp=grid::gpar(fontsize=12))) +
  specs + plot_spacer() + grid::textGrob("Age-at-Capture",
                                         gp=grid::gpar(fontsize=12)) +
  plot_layout(ncol=2,nrow=2,byrow=2,widths=c(0.08,0.92),heights=c(0.96,0.04))
p

ggsave("docs/Figure_2.pdf",p,device="pdf",
       width=3.5,height=3.75,units="in",dpi=1000,scale=1)
ggsave("docs/Figure_2.tiff",p,device="tiff",
       width=3.5,height=3.75,units="in",dpi=300,scale=1)
