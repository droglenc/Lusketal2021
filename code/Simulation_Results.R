# This script is used to analyze the simulation results recorded in
# "data/Simulation_Run_Results.csv" (which was created with "Run_Simulation.R").
# This ultimately produces Tables 2-5 and Figures 3 and 4.

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
## Load packages specific to this script
library(ggplot2)
library(ggnewscale)
library(patchwork)

#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=
# Get Data ----
#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=
## HELPER -- Extracts particular models results from the main data.frame ...
##           will be used to change the format of the data
getModelDF <- function(d,mdl) {
  d %>%
    dplyr::select(group:numLcat,dplyr::ends_with(mdl)) %>%
    dplyr::mutate(model=mdl) %>%
    dplyr::rename_at(dplyr::vars(ends_with(mdl)),
               .funs=function(x) stringr::str_replace_all(x,paste0("_",mdl),""))
}

## HELPER -- Extracts particular vonB parameter results from the main
##           data.frame ... will be used to change the format of the data
getParamMeasureDF <- function(d,param) {
  param <- paste0("_",param)
  d %>%
    dplyr::select(group:CV,model,dplyr::ends_with(param)) %>%
    tidyr::pivot_longer(
      cols=dplyr::ends_with(param),
      names_to="tmp",
      values_to="value") %>%
    tidyr::separate(col=tmp,into=c("measure","param"),sep="_")
}


## Load results data ... from running Run_Simulation.R
## Did not to use raw diffrences so did not include dplyr::starts_with("diff")
## in the last select()
df <- read.csv("data/Simulation_Run_Results.csv") %>%
  dplyr::mutate(n=factor(n),
                CV=factor(CV),
                group=factor(group,levels=groups)) %>%
  dplyr::select(group:numAged,percAged,numLcat,Linf_Full:LenZ_EP,
                dplyr::starts_with("err"))
str(df)

## Reshaping the data.frame
df_full <- getModelDF(df,"Full")
df_AO <- getModelDF(df,"AO")
df_WM <- getModelDF(df,"WM")
df_AA <- getModelDF(df,"AA")
df_RW <- getModelDF(df,"RW")
df_EP <- getModelDF(df,"EP")
df2 <- rbind(df_AO,df_WM,df_AA,df_RW,df_EP) %>%
  dplyr::mutate(model=factor(model,levels=c("AO","WM","AA","RW","EP")))

## Reshaping the data.frame (again)
df_Linf <- getParamMeasureDF(df2,"Linf")
df_K <- getParamMeasureDF(df2,"K")
df_t0 <- getParamMeasureDF(df2,"t0")
df_LenX <- getParamMeasureDF(df2,"LenX")
df_LenY <- getParamMeasureDF(df2,"LenY")
df_LenZ <- getParamMeasureDF(df2,"LenZ")
df3 <- rbind(df_Linf,df_K,df_t0,df_LenX,df_LenY,df_LenZ) %>%
  dplyr::mutate(param=factor(param,levels=c("Linf","K","t0",
                                            "LenX","LenY","LenZ")),
                accurate=case_when(
                  param %in% c("Linf","K","t0") & abs(value)<10 ~ "Accurate",
                  param %in% c("LenX","LenY","LenZ") & abs(value)<10 ~ "Accurate",
                  TRUE ~ "Inaccurate"
                ),
                accurate=factor(accurate,levels=c("Inaccurate","Accurate")))


#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=
# General Summaries ----
#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=
## Compute means of the vonB parameters across replicates to see how closely
## they correspond to the parameters set for the simulations
sumParams <- df2 %>%
  dplyr::group_by(group,n,CV) %>%
  dplyr::summarize(mnLinf=mean(Linf,na.rm=TRUE),
                   mnK=mean(K,na.rm=TRUE),
                   mnt0=mean(t0,na.rm=TRUE)) %>%
  as.data.frame()
sumParams2 <- df2 %>%
  dplyr::group_by(group,CV) %>%
  dplyr::summarize(mnLinf=mean(Linf,na.rm=TRUE),
                   mnK=mean(K,na.rm=TRUE),
                   mnt0=mean(t0,na.rm=TRUE)) %>%
  as.data.frame()
sumParams3 <- df2 %>%
  dplyr::group_by(group) %>%
  dplyr::summarize(mnLinf=mean(Linf,na.rm=TRUE),
                   mnK=mean(K,na.rm=TRUE),
                   mnt0=mean(t0,na.rm=TRUE)) %>%
  as.data.frame()

## !! Distance from chosen Linf increases with increasing CV
cbind(Linfs,Ks,t0s)
sumParams3
sumParams2
sumParams


#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=
# Simulation Summaries Using Mean/SD ----
#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=
## Compute number of fish that were aged (using mean and SD)
sumAged <- df2 %>%
  dplyr::group_by(group,n,CV) %>%
  dplyr::summarize(mnNumAged=mean(numAged,na.rm=TRUE),
                   sdNumAged=sd(numAged,na.rm=TRUE)) %>%
  dplyr::mutate(mnNumAged=formatC(mnNumAged,format="f",digits=0),
                sdNumAged=formatC(sdNumAged,format="f",digits=1),
                numAged=paste0(mnNumAged," (",sdNumAged,")")) %>%
  dplyr::select(group,CV,n,numAged) %>%
  dplyr::arrange(group,CV,n) %>%
  as.data.frame()
sumAged

## Compute number of length categories (using mean and SD)
sumLcat <- df2 %>%
  dplyr::group_by(group,n,CV) %>%
  dplyr::summarize(mnLcat=mean(numLcat,na.rm=TRUE),
                   sdLcat=sd(numLcat,na.rm=TRUE)) %>%
  dplyr::mutate(mnLcat=formatC(mnLcat,format="f",digits=0),
                sdLcat=formatC(sdLcat,format="f",digits=1),
                numLcats=paste0(mnLcat," (",sdLcat,")")) %>%
  dplyr::select(group,CV,n,numLcats) %>%
  dplyr::arrange(group,CV,n) %>%
  as.data.frame()
sumLcat

## Compute number of tries to fit model (more than 1 means some
## convergence error)
sumTries <- df2 %>%
  dplyr::group_by(group,n,CV) %>%
  dplyr::summarize(percConvErrActual=sum(ifelse(triesActual>1,1,0))/max(rep)*100,
                   AA_pConvErr=sum(ifelse(triesAA>1,1,0))/max(rep)*100,
                   WM_pConvErr=sum(ifelse(triesWM>1,1,0))/max(rep)*100,
                   AO_pConvErr=sum(ifelse(triesAO>1,1,0))/max(rep)*100,
                   RW_pConvErr=sum(ifelse(triesRW>1,1,0))/max(rep)*100,
                   EP_pConvErr=sum(ifelse(triesEP>1,1,0))/max(rep)*100) %>%
  dplyr::mutate_at(dplyr::vars(dplyr::ends_with("pConvErr")),
                   formatC,format="f",digits=1) %>%
  dplyr::arrange(group,CV,n) %>%
  as.data.frame()
sumTries


### Percentage error
## Summary statistics for the error variables (using mean and SD)
sumErrs <- df3 %>%
  dplyr::filter(measure=="err") %>%
  dplyr::group_by(group,param,n,CV,model) %>%
  dplyr::summarize(mnErr=mean(value,na.rm=TRUE),
                   sdErr=sd(value,na.rm=TRUE),
                   percAcc=mean(as.numeric(accurate)-1)*100) %>%
  dplyr::mutate(bias=ifelse((mnErr-sdErr)>0|(mnErr+sdErr)<0,"Biased","Unbiased"),
                bias=factor(bias,levels=c("Unbiased","Biased")),
                reliability=ifelse(percAcc>=80,"Reliable","Unreliable"),
                reliability=factor(reliability,levels=c("Reliable","Unreliable")))

sumErrs2 <- sumErrs %>%
  dplyr::mutate_at(vars(dplyr::ends_with("Err")),formatC,format="f",digits=1) %>%
  dplyr::mutate(Err=paste0(mnErr," (",sdErr,")",ifelse(bias=="Biased","^"," ")),
                percAcc=paste0(formatC(percAcc,format="f",digits=1),
                               ifelse(reliability=="Unreliable","*"," "))) %>%
  dplyr::select(group:model,Err,percAcc) %>%
  dplyr::arrange(group,param,model,CV,n) %>%
  tidyr::pivot_wider(names_from=model,values_from=Err:percAcc) %>%
  dplyr::select(group:CV,dplyr::ends_with("AA"),dplyr::ends_with("WM"),
                dplyr::ends_with("AO"),dplyr::ends_with("RW"),
                dplyr::ends_with("EP")) %>%
  tidyr::pivot_wider(names_from=param,values_from=Err_AA:percAcc_EP) %>%
  dplyr::select(group,CV,n,dplyr::ends_with("Linf"),dplyr::ends_with("K"),
                dplyr::ends_with("t0"),dplyr::ends_with("LenX"),
                dplyr::ends_with("LenY"),dplyr::ends_with("LenZ")) %>%
  as.data.frame()
sumErrs2


################################################################################
## Make a summary table and write out to CSV
################################################################################
bys <- c("group","CV","n")
sumTbl <- dplyr::left_join(sumLcat,sumAged,by=bys) %>%
  dplyr::left_join(.,sumTries,by=bys) %>%
  dplyr::left_join(.,sumErrs2,by=bys) %>%
  dplyr::select(-percConvErrActual)
# remove t0, lenZ, numLcats, and convergence error results
sumTbl <- sumTbl %>%
  select(-dplyr::ends_with("t0")) %>%
  select(-numLcats) %>%
  select(-dplyr::ends_with("pConvErr")) %>%
  select(-dplyr::ends_with("LenZ"))
write.csv(sumTbl,"docs/Table2_raw.csv",quote=TRUE,row.names=FALSE)



################################################################################
## Visualizing the bias, reliability, and precision results
################################################################################
## HELPER -- making error plots with violins and mean/SD
errorPlot1 <- function(d1,d2,p,xlbl,ylbl,ylim=NULL,ybrks=waiver(),
                       int=0,err_ref=c(-10,10),trim=0.05) {
  # Isolate selected parameter in d1 and d2
  d1 <- dplyr::filter(d1,param==p)
  d2 <- dplyr::filter(d2,param==p)

  # Used to "automatically" "zoom" the y-axis (to remove some outlying points)
  if (is.null(ylim)) {
    tmp <- d1 %>%
      dplyr::group_by(group,n,CV) %>%
      dplyr::summarize(pn=quantile(value,probs=trim),
                       px=quantile(value,probs=1-trim))
    ylim <- c(min(tmp$pn),max(tmp$px))
  }
  ggplot() +
    geom_hline(yintercept=int,lwd=0.25) +
    geom_hline(yintercept=err_ref,lwd=0.25,lty="dashed",col="gray50") +
    geom_violin(data=d1,aes(x=model,y=value,fill=reliability),
                color="gray40",lwd=0.05,trim=FALSE) +
    scale_fill_manual(values=c("gray45","gray85")) +
    geom_errorbar(data=d2,aes(ymin=mnErr-sdErr,ymax=mnErr+sdErr,x=model),
                  lwd=0.3,width=0) +
    new_scale_fill() +
    geom_point(data=d2,aes(x=model,y=mnErr,fill=bias),
               pch=23,color="black",size=0.75) +
    scale_fill_manual(values=c("black","white")) +
    facet_grid(group~CV+n,labeller=labeller(CV=lblrCV,n=lblrn,group=lblrGrp)) +
    coord_cartesian(ylim=ylim) +
    scale_y_continuous(name=ylbl,labels=function(x) paste0(x,"%"),
                       breaks=ybrks) +
    scale_x_discrete(name=xlbl) +
    theme_ALK +
    # change spacing between columns of facets and angle of axis labels
    theme(panel.spacing.x=unit(c(1,1,3,1,1,3,1,1),"mm"),
          axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
          legend.position="none")
}

## Create a new data.frame that has the "bias" and "reliability" conclusions
## appended to each simulation result ... allows coloring violins based on
## these conclusions
tmp <- select(sumErrs,group,n,CV,model,param,bias,reliability)
df4 <- dplyr::right_join(df3,tmp,by=c("group","n","CV","model","param")) %>%
  select(-measure,-accurate)

## Error plots
peLinf <- errorPlot1(df4,sumErrs,"Linf",xlbl=" ",
                     ylbl=expression(bold(paste("Error in ",L[infinity]))))
peK <- errorPlot1(df4,sumErrs,"K",xlbl="Estimation Method",ylbl="Error in K",
                  ybrks=seq(-50,50,25))
#pet0 <- errorPlot1(df4,sumErrs,"t0",xlbl="Age-Length Key Method",
#                 ylbl=expression(bold(paste("Error in ",t[0]))))
# decided not to include t0 results
ErrPlot1A <- peLinf + peK +
  plot_layout(ncol=1)

ggsave("docs/Figure_3.pdf",ErrPlot1A,device="pdf",
       width=7.25,height=7.5,units="in",dpi=1000,scale=1)
ggsave("docs/Figure_3.tiff",ErrPlot1A,device="tiff",
       width=7.25,height=7.5,units="in",dpi=300,scale=1)


peLenX <- errorPlot1(df4,sumErrs,"LenX",ylim=c(-10,10),
                    xlbl=" ",ylbl="Error in Mean Length at Young Age")
peLenY <- errorPlot1(df4,sumErrs,"LenY",ylim=c(-10,13.5),
                    xlbl="Estimation Method",
                    ylbl="Error in Mean Length at Middle Age")
#peLenZ <- errorPlot1(df4,sumErrs,"LenZ",xlbl="Age-Length Key Method",
#                    ylbl="Error in Mean Length at Old Age")
ErrPlot2A <- peLenX + peLenY +
  plot_layout(ncol=1)

ggsave("docs/Figure_4.pdf",ErrPlot2A,device="pdf",
       width=7.25,height=7.5,units="in",dpi=1000,scale=1)
ggsave("docs/Figure_4.tiff",ErrPlot2A,device="tiff",
       width=7.25,height=7.5,units="in",dpi=300,scale=1)
