#' @name EPGrowth
#' 
#' @title XXX
#'
#' @description XXX
#'
#' @details XXX 
#' @param sample1 First stage sample from population (all with lengths, some with ages).
#' @param sample2 Second stage sample from \code{sample1} (all with lengths and ages).
#' @param vlen Variable that contains observed lengths (either in quotes or not).
#' @param vage Variable that contains observed ages (either in quotes or not).
#' @param bin_w Width of bins for length categories (was \code{bin_size} in Perrault et al.).
#' @param len_mult Multiplier for maximum length when making the full range of lengths.
#' @param info A list returned from \code{EPdata_prep}.
#' @param starts A list of named starting values (must have \code{Linf}, \code{K}, \code{t0}, \code{CV}).
#' @param verbose A logical that indicates whether messages should be returned by \code{MakeADFun}.
#' @param trace A numeric that indicates how often \code{nlminb} should return trace information. Set to 0 to shut-off trace nformation.
#' 
#' @author Derek H. Ogle, \email{derek@@derekogle.com}
#'
#' @return \code{EPdata_prep} returns a list with the following items:
#'   \itemize{
#'     \item \code{age} Vector of observed ages in \code{sample2} (not by individual).
#'     \item \code{len} Vector of observed lengths in \code{sample2} (not by individual). Was \code{vlength} in Perrault et al.
#'     \item \code{vn} Vector of individuals per combination of \code{age} and \code{vlength} in \code{sample2}.
#'     \item \code{lcat} Vector of length categories (for the entire range of possible observed lengths). Was \code{rlength} in Perrault et al.
#'     \item \code{ep} Vector of EP values for each length category in \code{lcat}. Was \code{vp} in Perrault et al.
#'     \item \code{age_ind} Vector of ages seen in \code{sample2}.
#'     \item \code{len_ind} Vector of lengths in entire range of possible observed lengths.
#'     \item \code{sample_size} Number of fish in \code{sample2}
#'     \item \code{bin_w} Width of length bins used.
#'  }
#'
#' @references 
#' 
#' @keywords manip
#'
#' @examples
#' ## None yet
#' 
#' @export
#' @rdname EPGrowth

EPdata_prep <- function(sample1,sample2,vlen,vage,bin_w,len_mult=2) {
  ## Add length categorization variables to both samples
  sample1 <- dplyr::mutate(sample1,lcat=FSA::lencat({{ vlen }},w=bin_w))
  sample2 <- dplyr::mutate(sample2,lcat=FSA::lencat({{ vlen }},w=bin_w))
  
  ## Get number per length bin in both samples
  len.by.bin.1 <- sample1 %>%
    dplyr::group_by(lcat) %>%
    dplyr::summarize(Nh=n()) %>%
    dplyr::ungroup()
  len.by.bin.2 <- sample2 %>%
    dplyr::group_by(lcat) %>%
    dplyr::summarize(nh=n()) %>%
    dplyr::ungroup()
  
  nmlen <- rlang::quo_name(rlang::enquo(vlen))
  
  ## Add Nh and nh to sample2 and then calculate EP values
  sample2 <- dplyr::left_join(sample2,len.by.bin.1,by="lcat") %>%
    dplyr::left_join(len.by.bin.2,by="lcat") %>%
    arrange({{ vlen }})
  tdat <- left_join(len.by.bin.1,len.by.bin.2,by="lcat") %>%
    dplyr::mutate(ep=nh/Nh)
  
  ## Expand EP values to full range of lengths
  len_min <- 0
  len_max <- len_mult*max(tdat$lcat)
  tdat_all <- data.frame(lcat=seq(len_min,len_max,by=bin_w)) %>%
    dplyr::left_join(tdat,by="lcat") %>%
    dplyr::mutate(Nh=ifelse(is.na(Nh),0,Nh),
                  nh=ifelse(is.na(nh),0,nh),
                  ep=ifelse(Nh>0,nh/Nh,1))
  ## ALK by length with length bin attached
  sel.tab <- sample2 %>%
    dplyr::group_by({{ vlen }},{{ vage }}) %>%
    dplyr::summarize(n=dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(lcat=FSA::lencat({{ vlen }},w=bin_w)) %>%
    dplyr::arrange({{ vage }},{{ vlen }})
  
  ## Prepare data list and return it
  list(
    age = unlist(dplyr::select(sel.tab,{{ vage }}),use.names=FALSE),
    len = unlist(dplyr::select(sel.tab,{{ vlen }}),use.names=FALSE),
    vn = sel.tab$n,
    lcat = tdat_all$lcat,
    ep = tdat_all$ep,
    age_ind = 1:max(unlist(dplyr::select(sel.tab,{{ vage }}),use.names=FALSE)),
    len_ind = len_min:(len_max+bin_w),
    sample_size = nrow(sample2),
    bin_w = bin_w
  )
}


#' @export
#' @rdname EPGrowth
#' 
EPrun <- function(info,starts,verbose=FALSE,trace=10) {
  ## Prepare starting values
  params <- list(
    log_Linf = log(starts$Linf),
    log_K = log(starts$K),
    t0 = starts$t0,
    log_CV = log(starts$CV),
    lambda = rep(0,length(info$age_ind)-1)
  )
  ## Pass data and parameters to TMB to make AMDB function
  obj.EP <- TMB::MakeADFun(info,params,DLL="EP_likelihood",
                           inner.control=list(maxit=5000,trace=FALSE),
                           silent=!verbose) 
  ## Run optimizer
  opt.EP <- nlminb(obj.EP$par,obj.EP$fn,obj.EP$gr,
                   control=list(trace=trace,eval.max=1000,iter.max=1000)) 
  ## Get results
  ## Get report from optimiser (2nd part never used)
  #rep.EP <- obj.EP$report()
  #ad_obj_EP <- sdreport(obj.EP,getReportCovariance = FALSE)
  ## Store results
  results <- c(exp(opt.EP$par[1:2]),opt.EP$par[3],exp(opt.EP$par[4]))
  names(results) <- c("Linf","K","t0","CV")
  results
}
