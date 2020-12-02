# This script is used to create age-length key simulation results. Ultimately
# the results of the simulations are recorded in "Simulation_Run_Results.csv"
# which are read into other scripts for analysis. Thus, this script IS ONLY RUN
# ONCE to create the simulation results. Run "Simulated_Data_Figure.R" to see
# what the choices for simulation parameters (which are in "Common.R") look like.
# Run "Simulation_Results.R" to visualize the results.

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


#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=
# Simulation Parameters ----
#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## These may be set by the user ... others set in Common.R
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
numReps <- 500   ## Number of times to run each simulation
KILL_REP <- 10   ## Limit to quit b/c non-convergence
WT_TOL <- 1000   ## Tolerance to consider weights artificially high in WM method

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## Likely don't change anything beneath here
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## Number of simulations (and find combinations of nPerBin and CVs)
numVar <- length(CVs)
numSims <- length(nPerBin)*numVar
simN <- rep(nPerBin,each=numVar)
simVar <- rep(CVs,times=length(nPerBin))
cbind(sim=1:numSims,simN,simVar)  # just to check

## Initiate data.frame to hold the results
results <- NULL


#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=
# Simulation Loop ----
#   Ignore cat()s as these are just progress messages to the console.
#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=
for (w in groups) {
  cat("## Processing",w)
  cat(" using Linf=",Linfs[[w]],", K=",Ks[[w]],", t0=",t0s[[w]]," ##\n",sep="")
  ## Get age structure specific to the group ... must deal with NAs
  ad <- ads[,w]
  ad <- ad[!is.na(ad)]
  minAge <- min(as.numeric(names(ad)))
  maxAge <- max(as.numeric(names(ad)))
  predAge <- predAges[,w]
  for (i in 1:numSims){
    ## Get the simulation parameters specific to the simulation
    n <- simN[i]
    var <- simVar[i]
    cat("== Sim #",i," (of ",numSims," for ",w,") with n = ",n,
        " and CV = ",var," ==\n",sep="")
    cat("-- Rep:")
    for (j in 1:numReps) {
      cat(" ",j,sep="")
      ## (Re)Set a counter for the repeat loop and number of tries by method
      k <- triesActual <- triesAA <- triesWM <- triesAO <- 1
      triesRW <- triesEP <- 1
      repeat {
        cat(letters[k])
        #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
        ## STEP 1 -- Create mock electrofishing data
        #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
        df <- vbDataGen(numCap,minAge,maxAge,ad,Linfs[[w]],Ks[[w]],t0s[[w]],
                        var,binsize)
        ## Number per length bin (for Chih method)
        sum.df.by.bin <- df %>%
          group_by(lcat25) %>%
          summarize(Nh=n()) %>%
          mutate(Eh=Nh/sum(Nh))

        #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
        ## STEP 2 -- Create age sample (i.e., "age" n per 25 mm length group),
        ##           length sample, and age-length key
        #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
        df2 <- FSAsim::sample4ALK(ageCap~lenCap,data=df,n.fixed=n,w=binSize)
        agesamp <- dplyr::filter(df2,!is.na(ageCap))
        lensamp <- dplyr::filter(df2,is.na(ageCap))
        alk <- prop.table(xtabs(~lcat25+ageCap,data=agesamp),margin=1)
        ## Number per length bin (for Chih method)
        sum.agesample.by.bin <- agesamp %>%
          group_by(lcat25) %>%
          summarize(nh=n()) %>%
          mutate(eh=nh/sum(nh))

        #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
        ## STEP 3a -- Estimate growth parameters using all fish in sample
        #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
        e1 <- try({
          suppressWarnings(starts1 <- FSA::vbStarts(lenCap~ageCap,data=df))
          nls1 <- nls(lenCap~vbgf(ageCap,Linf,K,t0),data=df,start=starts1)
          cfs1 <- coef(nls1)
        },silent=TRUE)
        err1 <- class(e1)=="try-error" # TRUE if convergence error
        if (err1) triesActual <- triesActual+1

        #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
        ## STEP 3b -- Estimate parameters using only fish in age sample (AO)
        #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
        if (!err1) {## Only run if STEP 3a did not err (for efficiency).
          e2 <- try({
            suppressWarnings(starts2 <- FSA::vbStarts(lenCap~ageCap,
                                                      data=agesamp))
            nls2 <- nls(lenCap~vbgf(ageCap,Linf,K,t0),data=agesamp,
                        start=starts2)
            cfs2 <- coef(nls2)
          },silent=TRUE)
          err2 <- class(e2)=="try-error" # TRUE if convergence error
          if (err2) triesAO <- triesAO+1
        }

        #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
        ## STEP 3c -- Estimate growth parameters using mean length-at-age from
        ##            applying the ALK and weighting (inverse of SE^2) (WM)
        #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
        if (!err1) {## Only run if STEP 3a did not err.
          ## Estimate n, mean TL, and SE TL for whole sample using ALK (as
          ##   described in Miranda and Bettoli). Weights are inverse of SE^2
          ##   (as described in equation c on page 766 of Kimura 1980). If wt is
          ##   missing, infinity (which means SD was 0 as happens if n==1 or if
          ##   an age occurs in only one length bin), or greater than WT_TOL,
          ##   then wt will be replaced with 1/10th of the smallest "good" wt.
          ##   Note that lfreqL result is n per length category in total sample.
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

          ## Fit vonB to mean length weighted by inverse of SE^2.
          e3 <- try({
            suppressWarnings(starts3 <- FSA::vbStarts(mean~age,data=lenAtAge))
            nls3 <- nls(mean~vbgf(age,Linf,K,t0),data=lenAtAge,
                        start=starts3,weights=wt)
            cfs3 <- coef(nls3)
          },silent=TRUE)
          err3 <- class(e3)=="try-error" # TRUE if convergence error
          if (err3) triesWM <- triesWM+1
        }

        #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
        ## STEP 3d -- Estimate growth parameters using fish combined from the
        ##            age sample and those in length sample that had ages
        ##            estimated from ALK using Isermann-Knight method  (AA)
        #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
        if (!err1) { ## Only run if STEP 3a did not err
          ## Apply observed ALK (ob) using assigned age method (aa)
          ageassigned <- FSA::alkIndivAge(alk,ageCap~lenCap,data=lensamp)
          ## Combine age sample and assigned age length sample
          comb <- rbind(agesamp,ageassigned)
          e4 <- try({
            suppressWarnings(starts4 <- FSA::vbStarts(lenCap~ageCap,data=comb))
            nls4 <- nls(lenCap~vbgf(ageCap,Linf,K,t0),data=comb,start=starts4)
            cfs4 <- coef(nls4)
          },silent=TRUE)
          err4 <- class(e4)=="try-error" # TRUE if convergence error
          if (err4) triesAA <- triesAA+1
        }

        #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
        ## STEP 3e -- Estimate growth parameters using fish in age sample but
        ##            with Chih's re-weighting scheme  (RW)
        #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
        if (!err1) { ## Only run if STEP 3a did not err
          ### Get reweighting value from summaries above
          agesamp <- agesamp %>%
            left_join(select(sum.df.by.bin,-Nh),by="lcat25") %>%
            left_join(select(sum.agesample.by.bin,-nh),by="lcat25") %>%
            mutate(rw=Eh/eh)

          e5 <- try({
            suppressWarnings(starts5 <- FSA::vbStarts(lenCap~ageCap,
                                                      data=agesamp))
            nls5 <- nls(lenCap~vbgf(ageCap,Linf,K,t0),data=agesamp,
                        start=starts5,weights=rw)
            cfs5 <- coef(nls5)
          },silent=TRUE)
          err5 <- class(e5)=="try-error" # TRUE if convergence error
          if (err5) triesRW <- triesRW+1
        }


        #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
        ## STEP 3f -- Estimate growth parameters using Perrault's EP method (EP)
        #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
        if (!err1) { ## Only run if STEP 3a did not err
          e6 <- try({
            suppressWarnings(starts6 <- FSA::vbStarts(lenCap~ageCap,
                                                      data=agesamp))
            starts6$CV <- 0.15
            tmb.info <- EPdata_prep(df,agesamp,lenCap,ageCap,
                                    bin_w=binSize,len_mult=1.1)
            cfs6 <- EPrun(tmb.info,starts=starts6,trace=0)
            cfs6 <- cfs6[-4]
          },silent=TRUE)
          err6 <- class(e6)=="try-error" # TRUE if convergence error
          if (err6) triesEP <- triesEP+1
        }

        #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
        ## If at least one convergence error then break and try the rep again.
        ## However, don't try more than KILL_REP times
        #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
        any_errs <- any(c(err1,err2,err3,err4,err5))
        if (!any_errs) {
          ## no errors!!, goto next rep
          cat("",sep="")
          break
        } else if (k>=KILL_REP) {
          ## tried too many times, goto next rep
          cat("",sep="")
          cfs1 <- cfs2 <- cfs3 <- cfs4 <- cfs5 <- cfs6 <- c(Linf=NA_real_,K=NA_real_,t0=NA_real_)
          break
        } else {
          ## Try rep again
          k <- k+1
          cat(",",sep="")
        }
      } # end repeat loop

      #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
      ## Store the vonB coefficients into the results data.frame
      #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
      tmp1 <- c(data.frame(w,n,var,j,
                           triesActual,triesAA,triesWM,triesAO,triesRW,triesEP,
                           nrow(agesamp),length(lfreqL),stringsAsFactors=FALSE),
                    c(cfs1,lenX1=vbgf(predAge[[1]],cfs1),
                      lenY1=vbgf(predAge[[2]],cfs1),lenZ1=vbgf(predAge[[3]],cfs1),
                      cfs2,lenX2=vbgf(predAge[[1]],cfs2),
                      lenY2=vbgf(predAge[[2]],cfs2),lenZ2=vbgf(predAge[[3]],cfs2),
                      cfs3,lenX3=vbgf(predAge[[1]],cfs3),
                      lenY3=vbgf(predAge[[2]],cfs3),lenZ3=vbgf(predAge[[3]],cfs3),
                      cfs4,lenX4=vbgf(predAge[[1]],cfs4),
                      lenY4=vbgf(predAge[[2]],cfs4),lenZ4=vbgf(predAge[[3]],cfs4),
                      cfs5,lenX5=vbgf(predAge[[1]],cfs5),
                      lenY5=vbgf(predAge[[2]],cfs5),lenZ5=vbgf(predAge[[3]],cfs5),
                      cfs6,lenX6=vbgf(predAge[[1]],cfs6),
                      lenY6=vbgf(predAge[[2]],cfs6),lenZ6=vbgf(predAge[[3]],cfs6)))
      results <- rbind(results,as.data.frame(tmp1))
    } # end of numReps loops (j)
    cat(" --\n")
  } # end of numSims loop (i)
  cat("\n")
} # end of group loop (w)


#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=
# Output Results ----
#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#=-=
## Name variables
##   rep: replicate number within the group, numcap, CV combination
##   tries: number of tries for model to converge for that replicate
##   numaged: number of fish in the age sample
##   percaged: percent of sampled fish that were in the age sample
names(results) <- c("group","n","CV","rep",
                    "triesActual","triesAA","triesWM","triesAO","triesRW","triesEP",
                    "numAged","numLcat",
                    paste(rep(c("Linf","K","t0","LenX","LenY","LenZ"),
                              length(mdls)),
                          rep(mdls.abb,each=6),sep="_"))
## Find raw (diff_) and percent (err_) differences between the estimated
## parameters for a replicate run using one of three different methods and the
## parameter estimates from the complete data set (i.e., before subsampling
## for age). Also the percent of sampled fish that were in the age sample.
results <- results %>%
  dplyr::mutate(percAged=numAged/numCap*100,
                diff_Linf_AO=Linf_AO-Linf_Full,
                diff_K_AO=K_AO-K_Full,
                diff_t0_AO=t0_AO-t0_Full,
                diff_LenX_AO=LenX_AO-LenX_Full,
                diff_LenY_AO=LenY_AO-LenY_Full,
                diff_LenZ_AO=LenZ_AO-LenZ_Full,
                diff_Linf_WM=Linf_WM-Linf_Full,
                diff_K_WM=K_WM-K_Full,
                diff_t0_WM=t0_WM-t0_Full,
                diff_LenX_WM=LenX_WM-LenX_Full,
                diff_LenY_WM=LenY_WM-LenY_Full,
                diff_LenZ_WM=LenZ_WM-LenZ_Full,
                diff_Linf_AA=Linf_AA-Linf_Full,
                diff_K_AA=K_AA-K_Full,
                diff_t0_AA=t0_AA-t0_Full,
                diff_LenX_AA=LenX_AA-LenX_Full,
                diff_LenY_AA=LenY_AA-LenY_Full,
                diff_LenZ_AA=LenZ_AA-LenZ_Full,
                diff_Linf_RW=Linf_RW-Linf_Full,
                diff_K_RW=K_RW-K_Full,
                diff_t0_RW=t0_RW-t0_Full,
                diff_LenX_RW=LenX_RW-LenX_Full,
                diff_LenY_RW=LenY_RW-LenY_Full,
                diff_LenZ_RW=LenZ_RW-LenZ_Full,
                diff_Linf_EP=Linf_EP-Linf_Full,
                diff_K_EP=K_EP-K_Full,
                diff_t0_EP=t0_EP-t0_Full,
                diff_LenX_EP=LenX_EP-LenX_Full,
                diff_LenY_EP=LenY_EP-LenY_Full,
                diff_LenZ_EP=LenZ_EP-LenZ_Full,
                err_Linf_AO=diff_Linf_AO/Linf_Full*100,
                err_K_AO=diff_K_AO/K_Full*100,
                err_t0_AO=diff_t0_AO/t0_Full*100,
                err_LenX_AO=diff_LenX_AO/LenX_Full*100,
                err_LenY_AO=diff_LenY_AO/LenY_Full*100,
                err_LenZ_AO=diff_LenZ_AO/LenZ_Full*100,
                err_Linf_WM=diff_Linf_WM/Linf_Full*100,
                err_K_WM=diff_K_WM/K_Full*100,
                err_t0_WM=diff_t0_WM/t0_Full*100,
                err_LenX_WM=diff_LenX_WM/LenX_Full*100,
                err_LenY_WM=diff_LenY_WM/LenY_Full*100,
                err_LenZ_WM=diff_LenZ_WM/LenZ_Full*100,
                err_Linf_AA=diff_Linf_AA/Linf_Full*100,
                err_K_AA=diff_K_AA/K_Full*100,
                err_t0_AA=diff_t0_AA/t0_Full*100,
                err_LenX_AA=diff_LenX_AA/LenX_Full*100,
                err_LenY_AA=diff_LenY_AA/LenY_Full*100,
                err_LenZ_AA=diff_LenZ_AA/LenZ_Full*100,
                err_Linf_RW=diff_Linf_RW/Linf_Full*100,
                err_K_RW=diff_K_RW/K_Full*100,
                err_t0_RW=diff_t0_RW/t0_Full*100,
                err_LenX_RW=diff_LenX_RW/LenX_Full*100,
                err_LenY_RW=diff_LenY_RW/LenY_Full*100,
                err_LenZ_RW=diff_LenZ_RW/LenZ_Full*100,
                err_Linf_EP=diff_Linf_EP/Linf_Full*100,
                err_K_EP=diff_K_EP/K_Full*100,
                err_t0_EP=diff_t0_EP/t0_Full*100,
                err_LenX_EP=diff_LenX_EP/LenX_Full*100,
                err_LenY_EP=diff_LenY_EP/LenY_Full*100,
                err_LenZ_EP=diff_LenZ_EP/LenZ_Full*100)
## Write results out to a CSV file
write.csv(results,file="data/Simulation_Run_Results.csv",row.names=FALSE)
