# This script produces the plot that illustrates examples of length and age data
# with model fits. This is Figure 5.

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
set.seed(3432)

## HELPER -- big helper function to make the main plot
visFit <- function(w,n,var,ads,ddgw=0.7,ptsz=1,lnsz=1,lntyp="solid",lbl="A",lblsz=1,
                   xmax=NA,ymax=NA,point_guide=FALSE,method_guide=FALSE) {
  ## Get constants
  binsize <- 25
  WT_TOL <- 1000   ## Tolerance to consider weights artificially high in WM
  ad <- ads[,w]
  ad <- ad[!is.na(ad)]
  minAge <- min(as.numeric(names(ad)))
  maxAge <- max(as.numeric(names(ad)))
  
  ## Make data and get model fits
  ## STEP 1 -- Create mock electrofishing data
  df <- vbDataGen(numCap,minAge,maxAge,ad,Linfs[[w]],Ks[[w]],t0s[[w]],
                  var,binSize)
  ## Number per length bin (for Chih method)
  sum.df.by.bin <- df %>%
    group_by(lcat25) %>%
    summarize(Nh=n()) %>%
    mutate(Eh=Nh/sum(Nh))
  
  ## STEP 2 -- Create age sample (i.e., "age" n per 25 mm length group),
  ##           length sample, and age-length key
  df2 <- FSAsim::sample4ALK(ageCap~lenCap,data=df,
                            n.fixed=n,w=binSize)
  agesamp <- dplyr::filter(df2,!is.na(ageCap))
  lensamp <- dplyr::filter(df2,is.na(ageCap))
  alk <- prop.table(xtabs(~lcat25+ageCap,data=agesamp),margin=1)
  ## Number per length bin (for Chih method)
  sum.agesample.by.bin <- agesamp %>%
    group_by(lcat25) %>%
    summarize(nh=n()) %>%
    mutate(eh=nh/sum(nh))
  
  ## STEP 3a -- Estimate growth parameters using all fish in total sample
  suppressWarnings(starts1 <- FSA::vbStarts(lenCap~ageCap,data=df))
  nls1 <- nls(lenCap~vbgf(ageCap,Linf,K,t0),data=df,start=starts1)
  cfs1 <- coef(nls1)
  
  ## STEP 3b -- Estimate growth parameters using only fish in age sample
  suppressWarnings(starts2 <- FSA::vbStarts(lenCap~ageCap,data=agesamp))
  nls2 <- nls(lenCap~vbgf(ageCap,Linf,K,t0),data=agesamp,start=starts2)
  cfs2 <- coef(nls2)
  
  ## STEP 3c -- Estimate growth parameters using mean length-at-age from
  ##            applying the ALK and weighting. Catch  SD values that are
  ##            NA or less than WT_TOL and replace with WT_TOL. This
  ##            will make their weights small but not zero, as is done
  ##            in growth() in fish_methods.
  lfreqL <- xtabs(~lcat25,data=df)
  lfreqA <- xtabs(~lcat25,data=agesamp)
  ageDist <- FSA::alkAgeDist(alk,lfreqA,lfreqL) %>%
    dplyr::mutate(n=round(prop*numCap,1))
  lenAtAge <- FSA::alkMeanVar(alk,lenCap~lcat25+ageCap,
                              data=agesamp,len.n=lfreqL) %>%
    dplyr::full_join(.,ageDist,by="age") %>%
    dplyr::mutate(se=sd/sqrt(n),owt=1/(se^2),
                  wt=ifelse(is.na(owt)|owt>WT_TOL,
                            min(owt,na.rm=TRUE)/10,owt)) %>%
    dplyr::select(age,n,mean,sd,se,owt,wt)
  suppressWarnings(starts3 <- FSA::vbStarts(mean~age,data=lenAtAge))
  nls3 <- nls(mean~vbgf(age,Linf,K,t0),data=lenAtAge,start=starts3,weights=wt)
  cfs3 <- coef(nls3)
  
  ## STEP 3d -- Estimate growth parameters using fish in age sample and
  ##            those in length sample with ages estimated from ALK
  ## Apply observed ALK (ob) using assigned age method (aa)
  mntlobaa <- FSA::alkIndivAge(alk,ageCap~lenCap,data=lensamp)
  ## Combine age sample and assigned age length sample
  comb <- rbind(agesamp,mntlobaa) 
  suppressWarnings(starts4 <- FSA::vbStarts(lenCap~ageCap,data=comb))
  nls4 <- nls(lenCap~vbgf(ageCap,Linf,K,t0),data=comb,start=starts4)
  cfs4 <- coef(nls4)
  
  ## STEP 3e -- Estimate growth parameters using fish in age sample but
  ##            with Chih's reweighting scheme (RW)
  ### Get reweighting value from summaries above
  agesamp <- agesamp %>%
    left_join(select(sum.df.by.bin,-Nh),by="lcat25") %>%
    left_join(select(sum.agesample.by.bin,-nh),by="lcat25") %>%
    mutate(rw=Eh/eh)
  starts5 <- FSA::vbStarts(lenCap~ageCap,data=agesamp)
  nls5 <- nls(lenCap~vbgf(ageCap,Linf,K,t0),data=agesamp,
              start=starts5,weights=rw)
  cfs5 <- coef(nls5)  
  
  ## STEP 3f -- Estimate growth parameters using Perrault's EP method (EP)
  starts6 <- FSA::vbStarts(lenCap~ageCap,data=agesamp)
  starts6$CV <- 0.15
  tmb.info <- EPdata_prep(df,agesamp,lenCap,ageCap,bin_w=binSize,len_mult=1.1)
  cfs6 <- EPrun(tmb.info,starts=starts6,trace=0)
  cfs6 <- cfs6[-4]
  
  ## Combine coefficients from each model fitting
  cfs <- rbind(cfs1,cfs2,cfs3,cfs4,cfs5,cfs6)
  rownames(cfs) <- c("Entire","AO","WM","AA","RW","EP")
   
  ## Plot results
  ## Plot of whether fish was sampled for ageing or not
  ### Add whether fish was aged (or not) and note that this is the Entire sample
  df <- dplyr::mutate(df,sampled=ifelse(!is.na(df2$ageCap),"Yes","No"))
  T <- data.frame(t=seq(-1,max(df$ageCap)+1,0.1))
  
  ggplot() +
    geom_point(data=df,aes(x=ageCap,y=lenCap,group=sampled,color=sampled),
               position=position_dodge(width=ddgw),size=ptsz) +
    scale_color_manual(name="Aged-Sample",
                       breaks=c("Yes","No"),
                       values=c("No"="gray30","Yes"="gray70"),
                       labels=c("Aged","Not Aged"),
                       guide=point_guide) +
    ggnewscale::new_scale_color() +
    stat_function(data=T,aes(x=t,Y=NULL,color="Entire"),
                  fun=vbgf,args=list(Linf=cfs["Entire",]),
                  size=lnsz["Entire"],linetype=lntyp["Entire"]) +
    stat_function(data=T,mapping=aes(x=t,Y=NULL,color="AO"),
                  fun=vbgf,args=list(Linf=cfs["AO",]),
                  size=lnsz["AO"],linetype=lntyp["AO"]) +
    stat_function(data=T,aes(x=t,Y=NULL,color="WM"),
                  fun=vbgf,args=list(Linf=cfs["WM",]),
                  size=lnsz["WM"],linetype=lntyp["WM"]) +
    stat_function(data=T,aes(x=t,Y=NULL,color="AA"),
                  fun=vbgf,args=list(Linf=cfs["AA",]),
                  size=lnsz["AA"],linetype=lntyp["AA"]) +
    stat_function(data=T,aes(x=t,Y=NULL,color="RW"),
                  fun=vbgf,args=list(Linf=cfs["RW",]),
                  size=lnsz["RW"],linetype=lntyp["RW"]) +
    stat_function(data=T,aes(x=t,Y=NULL,color="EP"),
                  fun=vbgf,args=list(Linf=cfs["EP",]),
                  size=lnsz["EP"],linetype=lntyp["EP"]) +
    scale_color_manual(name="Method",
                       breaks=c("AO","WM","AA","RW","EP","Entire"),
                       values=c("AO"="#D55E00","WM"="#E69F00",
                                "AA"="#0072B2","RW"="#009E73",
                                "EP"="#CC79A7","Entire"="black"),
                       guide=method_guide) +
    scale_x_continuous(limits=c(0,xmax),expand=expansion(mult=c(0,0))) +
    scale_y_continuous(limits=c(0,ymax),expand=expansion(mult=c(0,0))) +
    ggpmisc::geom_text_npc(data=data.frame(),aes(npcx=0.02,npcy=0.98,label=lbl),
                           size=lblsz,hjust=0,vjust=1) +
    theme_ALK +
    theme(axis.title=element_blank(),
          legend.title=element_blank(),
          legend.key.size = unit(0.5,"line"),
          legend.margin = margin(0,0,0,0),
          legend.spacing.y = unit(0.07,"line"),
          legend.spacing.x=unit(2,"pt"),
          legend.text=element_text(size=5))
}


ws <- c("Black Crappie","Black Crappie","Striped Bass","Striped Bass")
ns <- c(5,15,5,15)
vars <- c(0.15,0.15,0.15,0.15)

ptsz <- 0.75
lnsz <- c("AO"=0.5,"WM"=0.5,"AA"=0.5,"RW"=0.5,"EP"=0.5,"Entire"=1.5)
lntyp <- c("AO"="84","WM"="84","AA"="8282",
           "RW"="6464","EP"="8484","Entire"="solid")
lblsz <- 3

p1 <- visFit(ws[1],ns[1],vars[1],ads,ddgw=0.5,ptsz=ptsz,lnsz=lnsz,lntyp=lntyp,
             lbl="A",lblsz=lblsz,ymax=400,point_guide="legend") +
  theme(legend.position=c(0.75,0.1))
p2 <- visFit(ws[2],ns[2],vars[2],ads,ddgw=0.5,ptsz=ptsz,lnsz=lnsz,lntyp=lntyp,
             lbl="B",lblsz=lblsz,ymax=400)
p3 <- visFit(ws[3],ns[3],vars[3],ads,ddgw=0.7,ptsz=ptsz,lnsz=lnsz,lntyp=lntyp,
             lbl="C",lblsz=lblsz,ymax=1200,
             method_guide=guide_legend(
               ncol=2,
               override.aes=list(linetype=lntyp,size=lnsz))) +
  theme(legend.position=c(0.7,0.15))
p4 <- visFit(ws[4],ns[4],vars[4],ads,ddgw=0.7,ptsz=ptsz,lnsz=lnsz,lntyp=lntyp,
             lbl="D",lblsz=lblsz,ymax=1200)

p <- p1 + p2 + p3 + p4 + plot_layout(ncol=2,nrow=2)
p <- wrap_elements(grid::textGrob("Total Length (mm)",rot=90,
                                  gp=grid::gpar(fontsize=12)),clip=FALSE) +
  p + plot_spacer() + grid::textGrob("Age-at-Capture",
                                     gp=grid::gpar(fontsize=12)) +
  plot_layout(ncol=2,nrow=2,byrow=2,widths=c(0.02,0.98),heights=c(0.98,0.02))
p

ggsave("docs/Figure_5.pdf",p,device="pdf",
       width=3.75,height=3.5,units="in",dpi=1000,scale=1)
ggsave("docs/Figure_5.tiff",p,device="tiff",
       width=3.75,height=3.5,units="in",dpi=300,scale=1)
