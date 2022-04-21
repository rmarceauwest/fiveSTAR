##=======================##
## minP-based adaptation ##
##=======================##

#simplified from formulas given by Nadarajah and Kotz 2008
#for distribution of min(Z1,Z2) two correlated N(mu,sigma2) variables
minZpdf = function(y,rho){
  fy = 2*dnorm(y)*pnorm((y*(rho-1))/sqrt(1-rho^2))
  return(fy)
}

#' Estimate amalgamated treatment effect using minP approach
#'
#' Estimates overall treatment effect amalgamated from strata-level treatment
#' effect estimates, using optimal weighting (i.e., using the the weighting scheme
#' that gives the best test statistics/minimum p-value, between ni and ni/sqrt(Vi)
#' weights), paying a minor penalty for taking the best of two correlated test
#' statistics
#'
#' @param sf      For survival traits, a \code{\link[survival]{survfit}}
#'  object containing information on number of subjects,
#'  number of events, number at risk, etc. for each strata, pooled across
#'  treatment assignment. For binary and continuous traits, a list of number of
#'  subjects per stratum, pooled across treatment assignment.
#' @param betas   Vector of estimated (log) treatment effect within each strata
#' @param vars    Vector of estimated variances of betas
#' @param cilevel Significance level for adaptive weighted confidence intervals
#' (i.e., (1-cilevel)x100\% CIs for two-tailed tests, and (1-2*cilevel)x100\% CIs
#' for 1-tailed tests) (default = 0.025)
#' @param alternative For tests, whether alternative hypothesis is "less",
#'  "greater", or "two.sided" (default = "less")
#' @param vartype Whether stratum-level variances used to calculate the correlation
#' between sample size and ni/sqrt(Vi) test statistics should be calculated
#' under the null hypothesis ("null") or estimated via the fitted model under the
#' alternative hypothesis ("alt", default). "null" variance may only be used when
#' measure == "HR"
#'
#' @return \itemize{
#'     \item adaptivewtRes: table of amalgamated (log) treatment effect estimate,
#'      variance, ci,and p-value using the optimal weight, as well as
#'      probability the treatment effect is in the desired direction
#'     \item adaptivewts: weights used (e.g., ni or ni/sqrt(Vi))
#'     \item calpha: adjusted critical value for test accounting for selecting
#'     the optimal weighting scheme
#'     \item singlewtres: vector of Z-score test statistics and corresponding
#'     p-values for sample size weights, ni/sqrt(Vi) weights, and inverse variance
#'     weights (1/Vi; output for knowledge though not used in calculations)
#' }
minPadapHR = function(sf,betas,vars,cilevel,alternative="less",vartype="alt"){

  #summary statistics
  ns = sf[[1]]
  Zs = betas/sqrt(vars)

  #=================================================#
  # calculating each beta, Z statistic, and p-value #
  #=================================================#

  #sample size weight values
  sswts = ns/sum(ns)
  beta1 = sum(sswts*betas)
  var1 = sum(sswts^2*vars)
  Z1 = beta1/sqrt(var1)
  if (alternative == "less"){
    p1 = pnorm(Z1,lower.tail=TRUE)
  } else if (alternative == "greater"){
    p1 = pnorm(Z1,lower.tail=FALSE)
  } else if (alternative == "two.sided"){
    p1 = 2*pnorm(abs(Z1),lower.tail=FALSE)
  }

  #"new weight" values (i.e., combining Z statistics over strata)
  newwts = (sswts/sqrt(vars))/(sum(sswts/sqrt(vars)))
  beta2 = sum(newwts*betas)
  var2 = sum(newwts^2*vars)
  Z2 = beta2/sqrt(var2)
  if (alternative == "less"){
    p2 = pnorm(Z2,lower.tail=TRUE)
  } else if (alternative == "greater"){
    p2 = pnorm(Z2,lower.tail=FALSE)
  } else if (alternative == "two.sided"){
    p2 = 2*pnorm(abs(Z2),lower.tail=FALSE)
  }

  #to output (not directly used!): inverse variance weights
  invarwts = (1/vars)/(sum(1/vars))
  beta3 = sum(invarwts*betas)
  var3 = sum(invarwts^2*vars)
  Z3 = beta3/sqrt(var3)
  if (alternative == "less"){
    p3 = pnorm(Z3,lower.tail=TRUE)
  } else if (alternative == "greater"){
    p3 = pnorm(Z3,lower.tail=FALSE)
  } else if (alternative == "two.sided"){
    p3 = 2*pnorm(abs(Z3),lower.tail=FALSE)
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  #==================================#
  # calculating minP and correlation #
  # (between tests 1 and 2 only)     #
  #==================================#

  #calculating correlation between the two test statistics
  #note! null variance formula is only really applicable for measure = "HR"!
  if (vartype == "null"){
    if (family != "cox" | measure != "HR") stop(
      'Please use vartype = "alt" when not doing hazard ratio-based analysis.')
    Vi0 = 4/summary(sf)$table[,"events"]
    rho = sum(ns^2*sqrt(Vi0))/((sqrt(sum(ns^2*Vi0)))*(sqrt(sum(ns^2))))
  } else if (vartype == "alt"){
    rho = sum(ns^2*sqrt(vars))/((sqrt(sum(ns^2*vars)))*(sqrt(sum(ns^2))))
  } else stop("vartype must be 'null' or 'alt'.")

  #results for each weighting scheme if no adaptive weighting is performed
  singlewtres = c(beta1,Z1,p1,beta2,Z2,p2,beta3,Z3,p3,rho)
  names(singlewtres) = c("beta1 (ni weights)","Z1","p1",
                         "beta2 (ni/sqrt(Vi) weights)","Z2","p2",
                         "beta3 (1/Vi weights)","Z3","p3","rho1,2")

  #calculating minP
  minP = min(p1,p2)
  usenewwt = ifelse(minP==p2,1,0)
  adapwts = sswts
  if (minP == p2) adapwts <- newwts

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  #==============================================================#
  # formula-based approach to final critical value for minP test #
  #==============================================================#

  #final estimates
  Zstar = ifelse(usenewwt,Z2,Z1)
  betastar = ifelse(usenewwt,beta2,beta1)
  Vstar = ifelse(usenewwt,var2,var1)

  ### final p-value = minP value = original p-value when only 1 strata
  if (length(sswts)==1){

    pstar = minP
    calpha = ifelse(alternative=="two.sided",qnorm(cilevel/2,0,1),
                    qnorm(cilevel,0,1))
    CIstar = c(betastar + calpha*sqrt(Vstar),
               betastar - calpha*sqrt(Vstar))

  } else { ## multiple strata formed

    if (alternative == "less"){

      #final p-values
      pstar = integrate(function(x) minZpdf(x,rho),lower=-Inf,upper=Zstar)$value

      #calculating critical value for CI
      #note: cilevel is assumed to be desired for correct alternative
      #(e.g., 0.025 for one-sided)
      calpha = uniroot(function(x) integrate(function(y)
        minZpdf(y,rho),lower=-Inf,upper=x)$value - cilevel, lower=-10,
        upper=0)$root

      #"95"% confidence intervals
      CIstar = c(betastar + calpha*sqrt(Vstar),
                 betastar - calpha*sqrt(Vstar))

    } else if (alternative == "greater"){

      #final p-values
      pstar = integrate(function(x) minZpdf(-x,rho),lower=Zstar,upper=Inf)$value

      #calculating critical value for CI
      #note: cilevel is assumed to be desired for correct alternative
      #(e.g., 0.025 for one-sided)
      calpha = uniroot(function(x) integrate(function(y)
        minZpdf(-y,rho),lower=x,upper=Inf)$value - cilevel, lower=0,
        upper=10)$root

      #"95"% confidence intervals
      CIstar = c(betastar - calpha*sqrt(Vstar),
                 betastar + calpha*sqrt(Vstar))

    } else if (alternative == "two.sided"){

      #final p-values
      pstar = 2*integrate(function(x) minZpdf(-x,rho),lower=abs(Zstar),
                          upper=Inf)$value

      #calculating critical value for CI
      #note: cilevel is assumed to be desired for correct alternative
      #(e.g., 0.025 for one-sided)
      calpha = uniroot(function(x) integrate(function(y)
        minZpdf(-y,rho),lower=x,upper=Inf)$value - cilevel/2, lower=0,
        upper=10)$root

      #"95"% confidence intervals
      CIstar = c(betastar - calpha*sqrt(Vstar),
                 betastar + calpha*sqrt(Vstar))

    }

  }

  #probability of test statistic being in the desired direction
  prlt0 = pnorm(0,betastar,sqrt(Vstar))
  prgt0 = pnorm(0,betastar,sqrt(Vstar),lower.tail=FALSE)

  #final amalgamated results
  minPRes = c(betastar,Vstar,CIstar,exp(betastar),exp(CIstar),Zstar,pstar,
              ifelse(alternative=="greater",prgt0,prlt0),calpha)
  names(minPRes) = c("beta","var","ci lower","ci upper","exp(beta)",
                      "exp(ci lower)","exp(ci upper)","minZ","pv",
                     ifelse(alternative=="greater","Pr(beta>0)","Pr(beta<0)"),
                     "calpha")

  #returing all final output
  return(list(adaptivewtRes=minPRes,adaptivewts=adapwts,#calpha=calpha,
              singlewtres=singlewtres))

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# ##=========================================##
# ## Railkar et al 2000 method for obtaining ##
# ## p-value from correlated tests           ##
# ##=========================================##
#
# ###should cilevel be up to user separate from alternative, or the 2-sided significance??
# calcastar = function(rho,nresamp=100000,alternative,cilevel){
#
#   #generating nresamp Z1, Z2 values from bivariate normal w/ correlation rho
#   ZZs = MASS::mvrnorm(n=nresamp,mu=c(0,0),Sigma=matrix(c(1,rho,rho,1),nrow=2))
#
#   #calculating p-values for each pair of Zs
#   if (alternative == "less"){
#     pvs = pnorm(ZZs,lower.tail=TRUE)
#   } else if (alternative == "greater"){
#     pvs = pnorm(ZZs,lower.tail=FALSE)
#   } else if (altnerative == "two.sided"){
#     pvs = 2*pnorm(abs(ZZs),lower.tail=FALSE)
#   } else stop ("Alternative must be one of 'greater', 'less', or 'two.sided'")
#
#   #same as min(x[1],x[2])...
#   #pvfinal = apply(pvs,1,function(x) min(min(x[1],x[2]),max(x[1],x[2])) )
#   phiZ = function(a,b,rho){
#     return(mvtnorm::pmvnorm(lower=c(-Inf,-Inf),upper=c(a,b),mean=c(0,0),
#                             sigma=matrix(c(1,rho,rho,1),nrow=2)))
#   }
#
#   #2-sided
#   calcpv2sided = function(k,cilevel,rho){
#     probftrH0 =
#       2*phiZ(qnorm(1-cilevel/2),qnorm(1-cilevel/(2*k)),rho) -
#       2*phiZ(qnorm(cilevel/2),qnorm(1-cilevel/(2*k)),rho) -
#       2*phiZ(qnorm(1-cilevel/2),qnorm(cilevel/(2*k)),rho) +
#       2*phiZ(qnorm(cilevel/2),qnorm(cilevel/(2*k)),rho) +
#       2*phiZ(qnorm(cilevel/2),qnorm(1-cilevel/2),rho) -
#       phiZ(qnorm(1-cilevel/2),qnorm(1-cilevel/2),rho) -
#       phiZ(qnorm(cilevel/2),qnorm(cilevel/2),rho)
#     rootval = probftrH0 - (1-cilevel)
#     return(rootval)
#   }
#   K2 = uniroot(function(x) calcpv2sided(x,cilevel,rho),interval=c(1,2))$root
#   a2 = cilevel/K2
#
#   #also calculate by formula if alpha = 0.025 (1-sided) or 0.05 (2-sided)
#   a2.fmla = cilevel + (1-rho)^0.55284 - (1+cilevel/2)*(1-rho)^0.53875
#   K2.fmla = cilevel/a2.fmla
#
#   a1.fmla = cilevel + (1-rho)^0.54858 - (1+cilevel/2)*(1-rho)^0.54151
#   K1.fmla = cilevel/a1.fmla
#   #for now..
#   if (cilevel == 0.025){ K1 = K1.fmla; a1 = a1.fmla }
#
#   alphastar = ifelse(alternative=="two.sided",a2,a1)
#   K = ifelse(alternative=="two.sided",K2,K1)
#
#   # #1-sided
#   # calcpv1sided = function(k,cilevel,rho){
#   #   probftrH0 =
#   #     2*phiZ(qnorm(1-cilevel/2),qnorm(1-cilevel/(2*k)),rho) -
#   #     2*phiZ(qnorm(cilevel/2),qnorm(1-cilevel/(2*k)),rho) -
#   #     2*phiZ(qnorm(1-cilevel/2),qnorm(cilevel/(2*k)),rho) +
#   #     2*phiZ(qnorm(cilevel/2),qnorm(cilevel/(2*k)),rho) +
#   #     2*phiZ(qnorm(cilevel/2),qnorm(1-cilevel/2),rho) -
#   #     phiZ(qnorm(1-cilevel/2),qnorm(1-cilevel/2),rho) -
#   #     phiZ(qnorm(cilevel/2),qnorm(cilevel/2),rho)
#   #   rootval = probftrH0 - (1-cilevel)
#   #   return(rootval)
#   # }
#   # K1 = uniroot(function(x) calcp12sided(x,cilevel,rho),interval=c(1,2))$root
#
#   return(list(alphastar,K))
#
# }
#
# adaptiveWtHR2 = function(sf,betas,vars,cilevel,measure,alternative,
#                          corrmethod="pearson",vartype){
#
#   # #currently not allowing other alternatives for binary, continuous traits
#   # #for consistency
#   # if (!(measure %in% c("HR","RMST"))){
#   #   alternative = "two.sided"
#   #   warning("2-sided p-values are reported.")
#   # }
#
#   #sf: a survfit object blinded to trt assignment
#   nns = sf[[1]]
#
#   #sample size weight values
#   sswts = nns/sum(nns)
#   beta1 = sum(sswts*betas)
#   var1 = sum(sswts^2*vars)
#   Z1 = beta1/sqrt(var1)
#   #Z1.1 = sum(nns*betas)/sqrt(sum(nns^2*vars))
#
#   #"new weight" values (i.e., combining Z statistics over strata)
#   newwts = (sswts/sqrt(vars))/(sum(sswts/sqrt(vars)))
#   beta2 = sum(newwts*betas)
#   var2 = sum(newwts^2*vars)
#   Z2 = beta2/sqrt(var2)
#
#   #inverse variance weight values
#   invarwts = (1/vars)/(sum(1/vars))
#   beta3 = sum(invarwts*betas)
#   var3 = sum(invarwts^2*vars)
#   Z3 = beta3/sqrt(var3)
#
#   #calculating asymptotic correlation between Z1 and Z2
#   if (vartype == NULL){
#     Vs = 4/(summary(sf)$table[,"events"])
#   } else Vs = vars
#   rho = sum(nns^2*sqrt(Vs))/(sqrt(sum(nns^2*Vs))*sqrt(sum(nns^2)))
#
#   #finding alphastar by simulation (w/ option for formula if ci = 0.025 (1-sided) or 0.05 (2-sided)?)
#   starvals = calcastar(rho=rho,alternative=alternative,cilevel=cilevel)
#   K = starvals$K
#   alphastar = starvals$alphastar
#
#   #calculate new weights
#   if (isnewwt){
#     newwts = (weights/sqrt(vars))/(sum(weights/sqrt(vars)))
#   } else newwts = weights
#
#   #overall treatment effect estimate, variance, and p-value
#   beta.newwt = sum(newwts*betas)
#   var.newwt = sum(newwts^2*vars)
#   Z.newwt = beta.newwt/sqrt(var.newwt)
#
#   if (alternative == "two.sided"){
#     pv.newwt = 2*pnorm(abs(Z.newwt),lower.tail=FALSE)
#   } else if (alternative == "greater"){
#     pv.newwt = pnorm(Z.newwt,lower.tail=FALSE)
#   } else if (alternative == "less"){
#     pv.newwt = pnorm(Z.newwt,lower.tail=TRUE)
#   } else stop("Alternative must be one of 'greater', 'less', or 'two.sided'.")
#
#   qcrit = qnorm(1-cilevel/2,0,1)
#   ci.newwt = cbind(beta.newwt - qcrit*sqrt(var.newwt),
#                    beta.newwt + qcrit*sqrt(var.newwt))
#
#   #collecting relevant results for return
#   newwtRes = c(beta.newwt,var.newwt,ci.newwt,exp(beta.newwt),exp(ci.newwt),
#                Z.newwt,pv.newwt)
#   names(newwtRes) = c("beta","var","ci lower","ci upper","exp(beta)",
#                       "exp(ci lower)","exp(ci upper)","Zstat","pv")
#
#   return(list(adaptivewtRes=newwtRes,adaptivewts=newwts,adaptivecorr=adaptivecorr))
# }

################################################################################

##============================================##
## Fit Cox Model Within Strata and Make Plots ##
##============================================##

#' Fit CoxPH Model Within Strata
#'
#' Fits a Cox proportional hazards model within each formed strata,
#' providing estimate of log hazard ratio and test of PH within each strata
#'
#' @param time    Follow-up time for right-censored data
#' @param status  Status indicator, 1=event, 0=censored
#' @param arm     Treatment indicator, 1 = test treatment, 0 = control
#' @param treeStrata Vector of strata membership for each subject
#'  (where 1 indicates belonging to the highest risk stratum, and the largest
#'   number indicates belonging to the lowest risk stratum)
#' @param termNodes Vector of names for each tree strata formed (terminal nodes
#' from ctree, representing strata definition in terms of covariates)
#' @param treetype String, whether trees input are "preliminary" (e.g.,
#' from step 3A) or "final" (e.g., from step 3B). Used only in plotting; ignored
#' when plot == FALSE
#' @param alternative For tests, whether alternative hypothesis is "less",
#'  "greater", or "two.sided" (default = "less")
#' @param cilevel Confidence level alpha for overall result and confidence
#' intervals (i.e., produces (1-cilevel)x100\% CIs when alternative="two.sided"
#' and (1-2*cilevel)x100\% CIs otherwise) (default = 0.025)
#' @param inclfrailty A logical variable - whether or not to include a frailty
#' term to each within-strata Cox PH model fit for added robustness against
#' unexplained heterogeneity
#' @param verbose Numeric variable indicating amount of information to print
#' to the terminal (0 = nothing, 1 = notes only, 2+ = notes and intermediate output)
#' @param plot    Logical, whether to create within strata and pooled between
#'  strata Kaplan-Meier plots
#' @param timeunit Optional argument, time unit for survival data
#' (e.g., Months, Years,..); currently only used for plots and ignored if
#' plot is FALSE
#' @param shading Logical variable; whether or not to show confidence bands around
#' Kaplan-Meier curves; ignored when plot != TRUE
#'
#' @return \itemize{
#'     \item fitsummary: summary of cox fits within each strata (estimated log
#'     hazard ratio, variance, (1-cilevel)x100\% CI, and corresponding
#'     exponentiated estimates), along with test statistic and p-value for test
#'     that logHR = 0 and Pr(logHR<0) for each stratum, Grambsch-Therneau (GT)
#'      tests for PH within each strata, and sample size and inverse variance
#'      weights for each stratum
#'     \item stratafit: a \code{\link[survival]{survfit}} object containing
#'     pooled-by-treatment survival information for each strata
#'     \item table: Summary of amalgamated (sample-size weighted) hazard ratio
#'     estimate, variance, (1-cilevel)x100\% confidence interval, test statistic,
#'     and p-value
#'     \item weights: Sample size weights used to construct estimate
#'     \item bystrataKM: Kaplan-Meier survival curve plots within each strata,
#'     returned if plot == TRUE
#'     \item betweenstrataKM: Kaplan-Meier survival curve plots from pooled
#'      treatment assignment data from each strata, returned if plot == TRUE
#' }
#' @import survival
coxbystrata = function(time,status,arm,treeStrata,termNodes=NULL,
                       treetype="final",alternative="less",cilevel=0.025,
                       inclfrailty=FALSE,verbose=0,plot=TRUE,
                       timeunit=NULL,shading=TRUE,shadealpha=0.3){

  #collecting summary data of strata
  dat = data.frame(time,status,arm)
  nstrata = max(treeStrata,na.rm=TRUE)
  if (is.null(termNodes)) termNodes = unique(treeStrata)

  #fitting cox PH model (with or without frailty term) within each formed strata
  if (inclfrailty){
    fsubjid = 1:length(treeStrata)
    coxCtreeStrataFit = lapply(1:nstrata,function(x){
      coxph(Surv(time,status) ~ as.factor(arm)+
              frailty(fsubjid,distribution="gamma"),
            subset=treeStrata==x,data=dat)
    })
  } else { #no frailty term
  coxCtreeStrataFit = lapply(1:nstrata,function(x){
    coxph(Surv(time,status) ~ as.factor(arm),subset=treeStrata==x,data=dat)
  })
  }

  #calculating GT test for PH within each strata
  pv.GTzph = sapply(1:nstrata,function(x)
    tryCatch(cox.zph(coxCtreeStrataFit[[x]], global = FALSE)[[1]][,"p"],
             error=function(e) NA))

  if (verbose >= 3){
    print("Cox model fits within each stratum:")
    print(coxCtreeStrataFit)
    print("Grambsch-Therneau test of PH within each stratum:")
    print(pv.GTzph)
  }

  #pulling out relevant summary statistics from cox fits to calculate coxSS
  #(amalgamated sample-size weighted log hazard ratio estimate)
  coxMat = matrix(unlist(lapply(coxCtreeStrataFit,function(x){
    unlist(x[c("coefficients","var","n")])
  } )),ncol=3,byrow=TRUE)
  colnames(coxMat) = c("bhatj","vhatj","nj")

  cnj = coxMat[,"nj"]
  coxbeta = coxMat[,"bhatj"]
  coxvar = coxMat[,"vhatj"]
  coxwSS = cnj/sum(cnj)

  #calculating coxSS (amalgamated logHR estimate using sample size weights)
  #and corresponding variance, test statistic, and p-value
  coxSS = sum(coxbeta*coxwSS)
  HRcoxSS = exp(coxSS)
  VcoxSS = sum(coxwSS^2*coxvar)
  TcoxSS = coxSS/sqrt(VcoxSS)
  if (alternative == "two.sided"){
    pvcoxSS = 2*pnorm(abs(TcoxSS),lower.tail=FALSE,log.p=FALSE)
  } else if (alternative == "greater"){
    pvcoxSS = pnorm(TcoxSS,lower.tail=FALSE,log.p=FALSE)
  } else if (alternative == "less"){
    pvcoxSS = pnorm(TcoxSS,lower.tail=TRUE,log.p=FALSE)
  } else stop ("Alternative must be one of 'two.sided', 'greater', or 'less'.")

  #calculating probability each betahat is < 0 (or > 0) assuming N(betahat, Vihat) dist'n
  if (alternative == "greater"){
    prlt0 = pnorm(0,coxbeta,sqrt(coxvar), lower.tail = FALSE)
  } else {
    prlt0 = pnorm(0,coxbeta,sqrt(coxvar), lower.tail = TRUE)
  }

  probname = ifelse(alternative=="greater","Pr(beta>0)","Pr(beta<0)")

  #calculating confidence intervals for estimated logHR and exponentiating
  #for estimated amalgamated hazard ratio
  if (alternative == "two.sided"){
    qcrit = qnorm(1-cilevel/2,0,1)
  } else qcrit = qnorm(1-cilevel,0,1)
  coxci = cbind(coxSS - qcrit*sqrt(VcoxSS),coxSS + qcrit*sqrt(VcoxSS))
  expcoxci = exp(coxci)

  #inverse variance weights for reference
  weight.invVar = (1/coxvar)/sum(1/coxvar)

  #amalgamated results using sample size weights
  coxMRSSres = list(table=data.frame(loghrSS=coxSS,vSS=VcoxSS,ciSS=coxci,
                                     hrSS=HRcoxSS,expciSS=expcoxci,
                                     TcoxSS=TcoxSS,pvSS=pvcoxSS),
                    weightSS=coxwSS,weightIV = weight.invVar)

  #confidence interval for logHR within each stratum
  coxstratci = cbind(coxbeta - qcrit*sqrt(coxvar),
                     coxbeta + qcrit*sqrt(coxvar))

  #test statistic and p-value for each stratum
  coxstratT = coxbeta/sqrt(coxvar)
  if (alternative == "two.sided"){
    coxstratpv = 2*pnorm(abs(coxstratT),0,1,lower.tail=FALSE)
  } else if (alternative == "greater"){
    coxstratpv = pnorm(coxstratT,0,1,lower.tail=FALSE)
  } else if (alternative == "less"){
    coxstratpv = pnorm(coxstratT,0,1,lower.tail=TRUE)
  } else stop ("Alternative must be one of 'two.sided', 'greater', or 'less'.")

  #summarized by-stratum results
  coxMRSSmat = cbind(coxbeta,coxvar,coxstratci,exp(coxbeta),exp(coxstratci),
                     coxstratT,coxstratpv,prlt0)
  colnames(coxMRSSmat) = c("bhat","v(bhat)","ci.lower","ci.upper","exp(bhat)",
                           "exp(ci.lower)","exp(ci.upper)","Zstat","pval",
                           probname)

  if (verbose > 2){
    print("Cox fit within each stratum summary:")
    print(coxMRSSmat)
  }

  coxMRSSmat.pt2 = matrix(cbind(unlist(coxwSS),unlist(weight.invVar)),
                          nrow=nstrata,byrow=FALSE)
  colnames(coxMRSSmat.pt2) = c("weight.SS","weight.invVar")

  if (verbose > 2){
    print("Additional weights for strata:")
    print(coxMRSSmat.pt2)
  }

  coxMRSSmat = cbind(coxMRSSmat,pv.GTzph,coxMRSSmat.pt2)

  if (verbose > 1){
    print("Cox fit within each strata")
    print(round(coxMRSSmat,4))
  }

  #KM curves ignoring treatment arm indicator
  Stratum = as.factor(treeStrata)
  s <- survfit(Surv(time,status)~Stratum,data=dat)

  #============================================================================#

  #------------------#
  # Plotting Results #
  #------------------#

  if (plot){

    timelabel = ifelse(is.null(timeunit),"",paste0(" (",timeunit,")"))

    #calculating time breaks for consistency between plots
    brtimeby = NULL
    breaktimes = ggplot2::waiver()
    if (!is.null(timeunit)){
      if (timeunit == "Months"){
        brtimes = c(3,6,12,18,24)
        brtimeby = brtimes[which.min(abs(5.5-max(round(time))/brtimes))]
        breaktimes = seq(0,max(time),brtimeby)
      } else if (timeunit == "Years"){
        brtimes = c(0.25,0.5,1,1.5,2)
        brtimeby = brtimes[which.min(abs(5.5-max(round(time))/brtimes))]
        breaktimes = seq(0,max(time),brtimeby)
      } else if (timeunit == "Days"){
        brtimes =  c(10,20,30,60,90,180,365,540,730)
        brtimeby = brtimes[which.min(abs(5.5-max(round(time))/brtimes))]
        breaktimes = seq(0,max(time),brtimeby)
      }
    }

    #-----------------------------------#
    # plot KM curves within each strata #
    #-----------------------------------#

    #plotting results
    new_df = with(dat,data.frame(arm=c(0,1)))
    #KM survival curves within each strata
    stratFit = lapply(1:nstrata,function(x){
      survfit(Surv(time,status) ~ as.factor(arm),subset=treeStrata==x,data=dat)
    })

    wrapper <- function(x, ...)
    {
      paste(strwrap(x, ...), collapse = "\n")
    }

    labelStart = "S"
    if (treetype == "preliminary") labelStart = "pS"
    if (treetype == "prespecified") labelStart = "S_"

    plotStrataSurvData = function(x){
      withCallingHandlers(survminer::ggsurvplot(
        stratFit[[x]], newdata=new_df, data= dat,subset=treeStrata==x,
        conf.int=FALSE,legend.labs = c("Control","Test"),
        xlab=paste0("Time",timelabel),legend=c(0.75,0.85), ylab="Survival",
        legend.title="Treatment",
        title=wrapper(paste0(labelStart,x," (",round(coxwSS[[x]]*100,1),
                             "% of subjects)"),width=45),
        ggtheme = survminer::theme_survminer(font.main = 14,font.legend=11),
        surv.median.line="hv",break.time.by=brtimeby,risk.table=TRUE,
        risk.table.height=0.2,xlim=c(0,max(time[status==1]))),
        #not printing warning about not reaching median survival time...
        #idea from Duncan Murdoch
        #https://r.789695.n4.nabble.com/Suppress-specific-warnings-td4664591.html
        warning = function(w) {
          if (grepl("Median survival not reached",w$message))
            invokeRestart("muffleWarning")
        })
    }

    splots = lapply(1:nstrata,plotStrataSurvData)
    nr = 1
    nc = nstrata
    p3 = survminer::arrange_ggsurvplots(splots, print = FALSE,ncol=nc, nrow=nr)

    #-------------------------------------------#
    # plots of all fits in the pooled A+B group #
    #-------------------------------------------#

    pseq=seq(.01,.99,by=.01)

    #s=controlStrataFit
    sdata<-data.frame(time=s$time, surv=s$surv, lower=s$lower, upper=s$upper)

    #only one strata
    if (is.null(s$strata)){
      sdata$strata = rep(1,nrow(sdata))
    } else {
      sdata$strata<-rep(names(s$strata), s$strata)
      #fix to ensure it won't mess up the order for > 9 strata
      sdata$strata = as.numeric(gsub(".*=","",sdata$strata))
    }

    s_table <- summary(s)$table
    if (is.null(dim(s_table))) s_table <- t(as.matrix(summary(s)$table,nrow=1))
    sevents <- s_table[,'events']
    #sevents <- summary(s)$table[,'events']
    pctevents <- sevents/sum(sevents)
    sdata$events <- pctevents[sdata$strata]
    sdata$ev5 <- (sdata$events >= 0.05)*1
    sdata$ev10 <- (sdata$events >= 0.10)*1
    sdata$ev10 <- factor(sdata$ev10,levels=c(0,1))

    labelNames = paste0("S",1:nstrata)
    if (treetype == "preliminary") labelNames = paste0("p",labelNames)
    stratLabelName = ifelse(treetype=="preliminary","Preliminary Risk",
                            "Identified Risk")
    if (treetype=="prespecified"){
      labelNames = paste0("S_",1:nstrata)
      stratLabelName = "Design"
    }

    if (!is.null(termNodes)){
      labelNames = paste0(labelNames,": ",termNodes)
      labelNames = sapply(labelNames,function(x) wrapper(x,width=50))
      names(labelNames) = NULL
    }

    if (nstrata <=3){

      spectraledges = RColorBrewer::brewer.pal(5,"Spectral")[
        sort(c(1,5,4)[1:nstrata])]
      p4 = ggplot2::ggplot() +
        ggplot2::geom_step(data=sdata,
                           ggplot2::aes(x=time, y=surv, col=factor(strata)),
                                        #alpha=ev10),
                           lwd=1.2) + ggplot2::theme_bw() +
        ggplot2::scale_x_continuous(breaks=breaktimes) +
        ggplot2::xlab(paste0("Time",timelabel)) + ggplot2::ylab("Survival") +
        ggplot2::scale_color_manual(name=paste0(stratLabelName," Strata"),
                                    labels=labelNames,values=spectraledges) #+
        # ggplot2::scale_alpha_manual(name='Percent Events',values=c(0.2,1.0),
        #                             labels=c('< 10 %','>= 10 %'),drop=FALSE)
      if (shading==TRUE){
        p4 = p4 +
          ggplot2::geom_ribbon(data=sdata,ggplot2::aes(
            x=time,ymin=lower,ymax=upper,fill=factor(strata)),alpha=shadealpha) +
          ggplot2::scale_fill_manual(name=paste0(stratLabelName," Strata"),
                                     labels=labelNames,values=spectraledges)
      }

    }else {

      getPalette = grDevices::colorRampPalette(RColorBrewer::brewer.pal(min(
        4,11),"Spectral"))
      #nstrata,11),"Spectral"))

      p4 = ggplot2::ggplot() +
        ggplot2::geom_step(data=sdata,
                           ggplot2::aes(x=time, y=surv,col=factor(strata)),
                                        #alpha=ev10),
                                        lwd=1.2) +
        ggplot2::theme_bw() +
        ggplot2::scale_x_continuous(breaks=breaktimes) +
        ggplot2::xlab(paste0("Time",timelabel)) + ggplot2::ylab("Survival") +
        ggplot2::scale_color_manual(name=paste0(stratLabelName," Strata"),
                                    labels=labelNames,values=getPalette(nstrata)) #+
        # ggplot2::scale_alpha_manual(name='Percent Events',values=c(0.2,1.0),
        #                             labels=c('< 10 %','>= 10 %'),drop=FALSE)
      if (shading==TRUE){
        p4 = p4 +
          ggplot2::geom_ribbon(data=sdata,ggplot2::aes(
            x=time,ymin=lower,ymax=upper,fill=factor(strata)),alpha=shadealpha) +
          ggplot2::scale_fill_manual(name=paste0(stratLabelName," Strata"),
                                     labels=labelNames,values=getPalette(nstrata))

      }

    }
  } else p3 = p4 = NULL

  #============================================================================#

  return(list(fitsummary=coxMRSSmat,stratafit=s,
              table=coxMRSSres$table,weights=coxMRSSres$weightSS,
              bystrataKM=p3,betweenstrataKM=p4))
}

################################################################################

##=============================================##
## Calculate RMST Within Strata and Make Plots ##
##=============================================##

#' Calculate RMST Within Strata
#'
#' Calculates the restricted mean survival time for each arm and corresponding
#' difference in RMST within each strata
#'
#' @param time    Vector of follow-up times for right-censored data
#' @param status  Status indicator vector, 1=event, 0=censored
#' @param arm     Treatment indicator vector, 1 = test treatment, 0 = control
#' @param tau     End time for RMST calculation; if not specified, the minimum
#' of maximum observed times for each arm is used (calculated separately for
#'  each stratum); see also \code{\link[survRM2]{rmst2}}
#' @param treeStrata Vector of strata membership, as defined by ctree strata
#' @param termNodes Vector of names for each tree strata formed (e.g., number of
#' strata or covariate definition)
#' @param treetype String, whether trees input are "preliminary" (e.g.,
#' from step 3A) or "final" (e.g., from step 3B)
#' @param alternative For tests, whether alternative hypothesis is "less",
#'  "greater", or "two.sided"
#' @param cilevel Confidence level alpha for overall result and confidence
#' intervals (i.e., produces (1-cilevel)x100\% CIs)
#' @param verbose Numeric variable indicating amount of information to print
#' to the terminal (0 = nothing, 1 = notes only, 2 = notes and intermediate output)
#' @param plot    Logical, whether to create within strata and pooled between
#'  strata KM plots; CURRENTLY IGNORED (i.e., no plots are produced)
#' @param timeunit Optional argument, time unit for survival data
#' (e.g., Months, Years,..); currently only used for plots and ignored if
#' plot is FALSE
#'
#' @return \itemize{
#'     \item fitsummary: summary of estimated RMST differences within each
#'      strata (est RMST difference, variance, (1-cilevel)x100\% confidence
#'      interval, test statistic, p-value, tau used, and Pr(RMST>0) for each
#'      stratum)
#'     \item stratafit: a \code{\link[survival]{survfit}} object containing
#'     pooled-by-treatment survival information for each strata
#'     \item weights: Sample size weights
#' }
rmstbystrata = function(time,status,arm,tau=NULL,treeStrata,termNodes=NULL,
                        treetype="final",alternative="two.sided",cilevel=0.05,
                        verbose=0,plot=FALSE,timeunit=NULL){

  #data summaries for computation
  dat = data.frame(time,status,arm)
  nstrata = length(unique(treeStrata))

  #calculating tau as minimum of maximum survival time, separately within each
  #strata if tau is not explicitly input
  if  (is.null(tau)){

    tausub = sapply(1:length(unique(treeStrata)),function(x){

      timesub = time[treeStrata==x]
      statussub= status[treeStrata==x]
      armsub=arm[treeStrata==x]
      datsub = cbind(timesub,statussub,armsub)

      taux = min(max(timesub[armsub==1]),max(timesub[armsub==0]))
      return(taux)

    })
    tau = min(tausub)

    #calculating RMST difference up to time tau within each stratum
    bystrataRMSTfits = lapply(1:length(unique(treeStrata)), function(x){
      rmstres = survRM2::rmst2(time=time[treeStrata==x],
                               status=status[treeStrata==x],
                               arm=arm[treeStrata==x],tau=tau)
      arm1 = rmstres$RMST.arm1$result["RMST",1:2]
      arm0 = rmstres$RMST.arm0$result["RMST",1:2]
      RMSTdiff = arm1[1]-arm0[1]
      RMSTvar = arm1[2]^2 + arm0[2]^2
      RMSTz = RMSTdiff/sqrt(RMSTvar)
      unadjresults = rmstres$unadjusted.result[1,]

      if (alternative == "two.sided"){
        RMSTpv = 2*pnorm(abs(RMSTz),lower.tail=FALSE,log.p=FALSE)#equivalent to unadjresults[4]
      } else if (alternative == "greater"){
        RMSTpv = pnorm(RMSTz,lower.tail=FALSE,log.p=FALSE)
      } else if (alternative == "less"){
        RMSTpv = pnorm(RMSTz,lower.tail=TRUE,log.p=FALSE)
      } else stop("Altnerative must be one of 'greater', 'less', or 'two.sided'.")

      return(data.frame(bhat=RMSTdiff,"v(bhat)"=RMSTvar,
                        ci.lower=unadjresults[2],ci.upper=unadjresults[3],
                        Zstat=RMSTz,pval=RMSTpv,tau=rmstres$tau))
    })
    bystrataRMSTfits = do.call(rbind,bystrataRMSTfits)

  } else if (is.numeric(tau)){ #using tau explicitly if given as input

    bystrataRMSTfits = lapply(1:length(unique(treeStrata)), function(x){
      rmstres = survRM2::rmst2(time=time[treeStrata==x],
                               status=status[treeStrata==x],
                               arm=arm[treeStrata==x],tau=tau,alpha=cilevel)
      arm1 = rmstres$RMST.arm1$result["RMST",1:2]
      arm0 = rmstres$RMST.arm0$result["RMST",1:2]
      RMSTdiff = arm1[1]-arm0[1]
      RMSTvar = arm1[2]^2 + arm0[2]^2
      RMSTz = RMSTdiff/sqrt(RMSTvar)
      unadjresults = rmstres$unadjusted.result[1,]

      if (alternative == "two.sided"){
        RMSTpv = 2*pnorm(abs(RMSTz),lower.tail=FALSE,log.p=FALSE)#equivalent to unadjresults[4]
      } else if (alternative == "greater"){
        RMSTpv = pnorm(RMSTz,lower.tail=FALSE,log.p=FALSE)
      } else if (alternative == "less"){
        RMSTpv = pnorm(RMSTz,lower.tail=TRUE,log.p=FALSE)
      } else stop("Altnerative must be one of 'greater', 'less', or 'two.sided'.")

      return(data.frame(bhat=RMSTdiff,"v(bhat)"=RMSTvar,
                        ci.lower=unadjresults[2],ci.upper=unadjresults[3],
                        Zstat=RMSTz,pval=RMSTpv,tau=rmstres$tau))
    })
    bystrataRMSTfits = do.call(rbind,bystrataRMSTfits,tau)

  } else stop("tau must be a numeric constant")

  #cleaning up results matrix
  labelNames = paste0("S",1:nstrata)
  if (treetype == "preliminary") labelNames = paste0("p",labelNames)
  rownames(bystrataRMSTfits) = labelNames
  wSS = c(table(treeStrata)/length(treeStrata))
  bystrataRMSTfits = data.frame(bystrataRMSTfits,weight.SS=wSS)
  names(bystrataRMSTfits)[2] = "v(bhat)"

  #============================================================================#

  #survival profiles/KM curve data ignoring treatment arm indicator
  Stratum = as.factor(treeStrata)
  s <- survfit(Surv(time,status)~Stratum,data=dat)

  #============================================================================#

  #calc pr(deltaRMST) > 0 for each strata
  prgt0 = sapply(1:nstrata,function(x){
    rmstout = survRM2::rmst2(time[treeStrata==x],status[treeStrata==x],
                             arm[treeStrata==x],tau=tau)
    arm0 = rmstout$RMST.arm0$rmst
    arm1 = rmstout$RMST.arm1$rmst
    deltaRMST = arm1["Est."] - arm0["Est."]
    varRMST = arm1["se"]^2 + arm0["se"]^2
    prgt0 = pnorm(0,deltaRMST,sqrt(varRMST),lower.tail=FALSE)
    return(prgt0)
  })

  bystrataRMSTfits = cbind(bystrataRMSTfits,prgt0)
  colnames(bystrataRMSTfits)[ncol(bystrataRMSTfits)] = "Pr(deltaRMST>0)"

  #============================================================================#

  ##will add plots later

  #output final results
  return(list(fitsummary=bystrataRMSTfits,stratafit=s,weights=wSS))

}

################################################################################

##===========================================================##
## Calculate Model Averaging of AFT Model Fits Within Strata ##
##===========================================================##

#' Calculate Model Averaged AFT Fits Within Strata
#'
#' Fits parametric AFT models within each formed stratum, then performs
#' model averaging to get a stratum-level estimate of time ratio with
#' corresponding estimated variance, confidence interval, and probability the
#' time ratio is > 0
#'
#' @param time    Vector of follow-up times for right-censored data
#' @param status  Vector of status indicators, 1=event, 0=censored
#' @param arm     Vector of treatment indicators, 1 = test treatment, 0 = control
#' @param treeStrata Vector of strata membership, as defined by ctree strata
#' @param termNodes Optional vector of names for each tree strata formed
#' (e.g., number of strata or covariate definition)
#' @param distList Vector of models (survival distributions for parametric AFT
#' fit) to include in model averaging; Each element must be the name of an
#' element from \code{\link[survival]{survreg.distributions}} (see also
#' \code{\link[survival]{survreg}}); default is c("loglogistic","lognormal","weibull")
#' @param ucvar Estimator for the unconditional variance of the model averaging
#' estimate to use. 1 uses Buckland et al. (1997) analytical estimator, and 2 uses
#' the more conservative estimator of Burnham and Anderson (2002) related to
#' Bayesian model averaging. For more details see Turek 2013. Default is "1".
#' @param alternative For tests, whether alternative hypothesis is "less",
#'  "greater", or "two.sided"
#' @param cilevel Confidence level alpha for within-stratum confidence
#' intervals (i.e., produces (1-cilevel)x100\% CIs when alternative="two.sided"
#' and (1-2*cilevel)x100\% CIs otherwise)
#' @param verbose Numeric variable indicating amount of information to print
#' to the terminal (0 = nothing, 1 = notes only, 2 = notes and intermediate output)
#' @param plot    Logical, whether to create within strata and pooled between
#'  strata Kaplan-Meier plots
#' @param treetype String, whether trees input are "preliminary" (e.g.,
#' from step 3A) or "final" (e.g., from step 3B). Used only in plotting; ignored
#' when plot == FALSE
#' @param timeunit Optional argument, time unit for survival data
#' (e.g., Months, Years,..); currently only used for plots and ignored if
#' plot is FALSE
#' @param shading Logical variable; whether or not to show confidence bands around
#' Kaplan-Meier curves; ignored when plot != TRUE
#'
#' @return \itemize{
#'     \item fitsummary: summary of estimated model averaged AFT time
#'     ratio estimate, variance, CI, and Pr(TR > 1) within each strata, as well
#'     as corresponding exponeniated estimates, within-stratum test statistic
#'     and p-value, and weights/AIC information for each model being averaged
#'     \item stratafits: a \code{\link[survival]{survfit}} object containing
#'     pooled-by-treatment survival information for each strata
#'     \item aftfits: list of each AFT model fit for each stratum
#'     \item bystrataKM: Kaplan-Meier survival curve plots within each strata,
#'     returned if plot == TRUE
#'     \item betweenstrataKM: Kaplan-Meier survival curve plots from pooled
#'      treatment assignment data from each strata, returned if plot == TRUE
#' }
#' @import survival
maaftbystrata = function(time,status,arm,treeStrata,termNodes=NULL,
                         distList=c("loglogistic","lognormal","weibull"),
                         ucvar=1,alternative="greater",cilevel=0.025,verbose=0,
                         plot=FALSE,treetype="final",timeunit=NULL,
                         shading=TRUE,shadealpha=0.3){

  #key aspects of input data
  dat = data.frame(time,status,arm)
  nstrata = length(unique(treeStrata))

  #calculating parametric AFT model within each stratum
  srvfits = list()
  for (stratum in 1:nstrata){
    srvfits[[stratum]] = list()

    #explicit in case log logistic model is not in distList
    llogfit = survreg(Surv(time,status)~arm,subset=treeStrata==stratum,
                      dist="loglogistic")

    for (dist in distList){

      if (dist=="loglogistic"){
        srvfits[[stratum]][[dist]] = llogfit

      #using loglogistic model fit as initial values for weibull fit
      #to overcome possible convergence issues
      } else if (dist=="weibull"){
        srvfits[[stratum]][[dist]] = survreg(
          Surv(time,status) ~ arm,subset=treeStrata==stratum,dist=dist,
          init=coef(llogfit),control=survreg.control(maxiter=100))

      } else {
        srvfits[[stratum]][[dist]] = survreg(
          Surv(time,status) ~ arm,subset=treeStrata==stratum,dist=dist)
      }

      # srvfits[[stratum]][[dist]] = tryCatch(
      #   survreg(Surv(time,status) ~ arm,subset=treeStrata==stratum,dist=dist),
      #   error = function(e) e)
      # #if weibull model fails to converge, use loglogistic model as initial values
      # if (inherits(srvfits[[stratum]][[dist]],"error")&(dist=="weibull")){
      #
      #   srvfits[[stratum]][[dist]] =  survreg(
      #     Surv(time,status) ~ arm,subset=treeStrata==stratum,dist=dist,
      #     init=coef(llogfit))
      # }
    }
  }

  #performing model averaging within each stratum
  MAaft = varMAaft1 = varMAaft2 = rep(NA,nstrata)
  MAcis = matrix(NA,ncol=2,nrow=nstrata)
  res = matrix(NA,ncol=4,nrow=nstrata)
  weightList = AIClist = list()
  for (stratum in 1:nstrata){
    #extracting estimate, variance, and AIC from each model fit
    #and performing model averaging within the strata

    modelList = srvfits[[stratum]]
    ests = vars = AICs = rep(NA,length(modelList))

    for (mod in 1:length(modelList)){
      model = modelList[[mod]]
      ests[mod] = model$coefficients[[2]]
      vars[mod] = (summary(model)$table[2,2])^2
      AICs[mod] = AICcmodavg::AICc(model) #same as 2*3-2*model[[stratum]]$loglik[2]
    }

    #need to calculate compared to min AIC so it is computable
    AICs[vars==0] <- Inf # if still not converging within a stratum
    deltaAICs = AICs - min(AICs)
    weights = exp(-0.5*deltaAICs)/sum(exp(-0.5*deltaAICs))
    weightList[[stratum]] = weights
    AIClist[[stratum]] = AICs

    #calculating model averaging estimate
    MAaft[[stratum]] = sum(weights*ests)
    #better to use var2, but for comparison
    varMAaft1[[stratum]] = (sum(weights*sqrt(vars+(ests-MAaft[[stratum]])^2)))^2
    varMAaft2[[stratum]] = sum(weights*(vars+(ests-MAaft[[stratum]])^2))

    #whether to use unconditional variance 1 or 2
    if (ucvar==1){
      varMAaft = varMAaft1[[stratum]]
    } else if (ucvar==2) varMAaft = varMAaft2[[stratum]]

    qq = ifelse(alternative=="two.sided",qnorm(cilevel/2),qnorm(cilevel))
    MAcis[stratum,] = c(MAaft[[stratum]] + qq*sqrt(varMAaft),
                         MAaft[[stratum]] - qq*sqrt(varMAaft))
    res[stratum,] = c(MAaft[[stratum]],varMAaft,MAcis[stratum,])

  }

  #overall results
  colnames(res) = c("bhat","v(bhat)","ci.lower","ci.upper")
  Zs = res[,"bhat"]/sqrt(res[,"v(bhat)"])
  if (alternative=="less"){
    pvals = pnorm(Zs,lower.tail=TRUE)
  } else if (alternative=="greater"){
    pvals = pnorm(Zs,lower.tail=FALSE)
  } else if (alternative=="two.sided"){
    pvals = 2*pnorm(abs(Zs),lower.tail=FALSE)
  }

  #probably true time ratio is > 0 (or < 0)
  if (alternative == "less"){
    prgt0 = pnorm(0,res[,"bhat"],sqrt(res[,"v(bhat)"]),lower.tail=TRUE)
  } else {
    prgt0 = pnorm(0,res[,"bhat"],sqrt(res[,"v(bhat)"]),lower.tail=FALSE)
  }

  probname = ifelse(alternative=="less","Pr(beta<0)","Pr(beta>0)")

  if (nstrata > 1){
    res = cbind(res,exp(res[,c("bhat","ci.lower","ci.upper")]),Zs,pvals,prgt0)
  } else {
    res = matrix(c(res,exp(res[,c("bhat","ci.lower","ci.upper")]),Zs,pvals,
                   prgt0),nrow=1)
  }

  colnames(res) = c("bhat","v(bhat)","ci.lower","ci.upper","exp(bhat)",
                    "exp(ci.lower)","exp(ci.upper)","Zstat","pval",probname)

  weightmat = do.call(rbind,weightList)
  AICmat = do.call(rbind,AIClist)
  colnames(weightmat) = paste0(distList,"_weight")
  colnames(AICmat) = paste0(distList,"_AICc")

  res = cbind(res,AICmat,weightmat)

  #strata fits (km)
  Stratum = as.factor(treeStrata)
  s <- survfit(Surv(time,status)~Stratum,data=dat)

  wSS = (s$n)/sum(s$n)
  res = cbind(res,wSS)
  colnames(res)[ncol(res)] = "weight.SS"

  if (plot){

    timelabel = ifelse(is.null(timeunit),"",paste0(" (",timeunit,")"))

    #calculating time breaks for consistency between plots
    brtimeby = NULL
    breaktimes = ggplot2::waiver()
    if (!is.null(timeunit)){
      if (timeunit == "Months"){
        brtimes = c(3,6,12,18,24)
        brtimeby = brtimes[which.min(abs(5.5-max(round(time))/brtimes))]
        breaktimes = seq(0,max(time),brtimeby)
      } else if (timeunit == "Years"){
        brtimes = c(0.25,0.5,1,1.5,2)
        brtimeby = brtimes[which.min(abs(5.5-max(round(time))/brtimes))]
        breaktimes = seq(0,max(time),brtimeby)
      } else if (timeunit == "Days"){
        brtimes =  c(10,20,30,60,90,180,365,540,730)
        brtimeby = brtimes[which.min(abs(5.5-max(round(time))/brtimes))]
        breaktimes = seq(0,max(time),brtimeby)
      }
    }

    #-----------------------------------#
    # plot KM curves within each strata #
    #-----------------------------------#

    #plotting results
    new_df = with(dat,data.frame(arm=c(0,1)))
    #KM survival curves within each strata
    stratFit = lapply(1:nstrata,function(x){
      survfit(Surv(time,status) ~ as.factor(arm),subset=treeStrata==x,data=dat)
    })

    wrapper <- function(x, ...)
    {
      paste(strwrap(x, ...), collapse = "\n")
    }

    labelStart = "S"
    if (treetype == "preliminary") labelStart = "pS"
    if (treetype == "prespecified") labelStart = "S_"

    plotStrataSurvData = function(x){
      withCallingHandlers(survminer::ggsurvplot(
        stratFit[[x]], newdata=new_df, data= dat,subset=treeStrata==x,
        conf.int=FALSE,legend.labs = c("Control","Test"),
        xlab=paste0("Time",timelabel),legend=c(0.75,0.85), ylab="Survival",
        legend.title="Treatment",
        title=wrapper(paste0(labelStart,x," (",round(wSS[[x]]*100,1),
                             "% of subjects)"),width=45),
        ggtheme = survminer::theme_survminer(font.main = 14,font.legend=11),
        surv.median.line="hv",break.time.by=brtimeby,risk.table=TRUE,
        risk.table.height=0.2,xlim=c(0,max(time[status==1]))),
        #not printing warning about not reaching median survival time...
        #idea from Duncan Murdoch
        #https://r.789695.n4.nabble.com/Suppress-specific-warnings-td4664591.html
        warning = function(w) {
          if (grepl("Median survival not reached",w$message))
            invokeRestart("muffleWarning")
        })
    }

    splots = lapply(1:nstrata,plotStrataSurvData)
    nr = 1
    nc = nstrata
    p3 = survminer::arrange_ggsurvplots(splots, print = FALSE,ncol=nc, nrow=nr)

    #-------------------------------------------#
    # plots of all fits in the pooled A+B group #
    #-------------------------------------------#

    pseq=seq(.01,.99,by=.01)

    #s=controlStrataFit
    sdata<-data.frame(time=s$time, surv=s$surv, lower=s$lower, upper=s$upper)
    #only one strata
    if (is.null(s$strata)){
      sdata$strata = rep(1,nrow(sdata))
    } else {
      sdata$strata<-rep(names(s$strata), s$strata)
      #fix to ensure it won't mess up the order for > 9 strata
      sdata$strata = as.numeric(gsub(".*=","",sdata$strata))
    }

    s_table <- summary(s)$table
    if (is.null(dim(s_table))) s_table <- t(as.matrix(summary(s)$table,nrow=1))
    sevents <- s_table[,'events']
    #sevents <- summary(s)$table[,'events']
    pctevents <- sevents/sum(sevents)
    sdata$events <- pctevents[sdata$strata]
    sdata$ev5 <- (sdata$events >= 0.05)*1
    sdata$ev10 <- (sdata$events >= 0.10)*1
    sdata$ev10 <- factor(sdata$ev10,levels=c(0,1))

    labelNames = paste0("S",1:nstrata)
    if (treetype == "preliminary") labelNames = paste0("p",labelNames)
    stratLabelName = ifelse(treetype=="preliminary","Preliminary Risk",
                            "Identified Risk")
    if (treetype=="prespecified"){
      labelNames = paste0("S_",1:nstrata)
      stratLabelName = "Design"
    }

    if (!is.null(termNodes)){
      labelNames = paste0(labelNames,": ",termNodes)
      labelNames = sapply(labelNames,function(x) wrapper(x,width=50))
      names(labelNames) = NULL
    }

    if (nstrata <=3){

      spectraledges = RColorBrewer::brewer.pal(5,"Spectral")[
        sort(c(1,5,4)[1:nstrata])]
      p4 = ggplot2::ggplot() +
        ggplot2::geom_step(data=sdata,
                           ggplot2::aes(x=time, y=surv, col=factor(strata)),
                                        #alpha=ev10),
                                        lwd=1.2) +
        ggplot2::theme_bw() +
        ggplot2::scale_x_continuous(breaks=breaktimes) +
        ggplot2::xlab(paste0("Time",timelabel)) + ggplot2::ylab("Survival") +
        ggplot2::scale_color_manual(name=paste0(stratLabelName," Strata"),
                                    labels=labelNames,values=spectraledges) #+
        # ggplot2::scale_alpha_manual(name='Percent Events',values=c(0.2,1.0),
        #                             labels=c('< 10 %','>= 10 %'))
      if (shading==TRUE){
        p4 = p4 +
          ggplot2::geom_ribbon(data=sdata,ggplot2::aes(
            x=time,ymin=lower,ymax=upper,fill=factor(strata)),alpha=shadealpha) +
          ggplot2::scale_fill_manual(name=paste0(stratLabelName," Strata"),
                                     labels=labelNames,values=spectraledges)
      }

    }else {

      getPalette = grDevices::colorRampPalette(RColorBrewer::brewer.pal(min(
        nstrata,11),"Spectral"))

      p4 = ggplot2::ggplot() +
        ggplot2::geom_step(data=sdata,
                           ggplot2::aes(x=time, y=surv,col=factor(strata)),
                                        #alpha=ev10),
                                        lwd=1.2) +
        ggplot2::theme_bw() +
        ggplot2::scale_x_continuous(breaks=breaktimes) +
        ggplot2::xlab(paste0("Time",timelabel)) + ggplot2::ylab("Survival") +
        ggplot2::scale_color_manual(name=paste0(stratLabelName," Strata"),
                                    labels=labelNames,values=getPalette(nstrata)) #+
        # ggplot2::scale_alpha_manual(name='Percent Events',values=c(0.2,1.0),
        #                             labels=c('< 10 %','>= 10 %'),drop=FALSE)
      if (shading==TRUE){
        p4 = p4 +
          ggplot2::geom_ribbon(data=sdata,ggplot2::aes(
            x=time,ymin=lower,ymax=upper,fill=factor(strata)),alpha=shadealpha) +
          ggplot2::scale_fill_manual(name=paste0(stratLabelName," Strata"),
                                     labels=labelNames,values=getPalette(nstrata))
      }

    }

  } else p3 = p4 = NULL

  #============================================================================#

  return(list(fitsummary=res,stratafits=s,aftfits=srvfits,
              bystrataKM=p3,betweenstrataKM=p4))

}

################################################################################

##=====================================##
## Make Forest Plot Summary of Results ##
##=====================================##

colfn <- local({
  function(..., clr.line, clr.marker){
    forestplot::fpDrawNormalCI(..., clr.line="darkblue",clr.marker="firebrick2")
  }
})

colfn2 <- local({
  function(..., clr.line, clr.marker){
    forestplot::fpDrawNormalCI(..., clr.line="forestgreen",clr.marker="green3")
  }
})

#ensuring fixed number of digits
printR = function(x,digits){
  return( formatC(round(x,digits), format='f',digits=digits) )
}

printE <- function(x,digits){
  formatC(round(x,digits), digits = digits+1, format = "g",
          drop0trailing = FALSE)
}

roundLabelNums = function(name,digits){
  name = strsplit(name," ")[[1]]
  #matching all numbers of the sort digit[punctuation]digit, except digit-digit
  numpos = grep("\\d+(?!\\-)[[:punct:]]\\d+",name,perl=TRUE)
  numparpos = grep("\\d+[[:punct:]]\\d+[[:punct:]]",name)
  puncttypeL = gsub("\\d+[[:punct:]]\\d+[[:punct:]]+$","",name[numparpos])
  puncttypeL[nchar(puncttypeL)>0] = paste0(sapply(unlist(strsplit(
    puncttypeL[nchar(puncttypeL)>0],"")),function(x)
      gsub("([[:punct:]])",paste0("\\\\","\\1"),x)),collapse="")
  #puncttypeL[nchar(puncttypeL)>0] = paste0("\\",puncttypeL[nchar(puncttypeL)>0])
  puncttypeR = gsub("^c|[[:punct:]]*\\d+[[:punct:]]\\d+","",name[numparpos])
  puncttypeR[nchar(puncttypeR)>0] = paste0("\\",puncttypeR[nchar(puncttypeR)>0])
  if (length(numparpos)!=0){
    name[numparpos] = sapply(1:length(puncttypeL),function(x)
      sub(puncttypeL[x],"",name[numparpos[x]]))
    name[numparpos] = sapply(1:length(puncttypeR),function(x)
      sub(puncttypeR[x],"",name[numparpos[x]]))
  }

  if (length(grep(",",name[numpos]))>0){
    namenums <- strsplit(name[numpos],",")[[1]]
    roundnum <- vector(mode = "list", length = length(namenums))
    for (j in 1:length(namenums)){
      roundnum[j] <- round(as.numeric(namenums[j]),digits)
    }
    name[numpos] <- paste0(unlist(roundnum),collapse=",")
  } else{
    name[numpos] = round(as.numeric(name[numpos]),digits)
  }
  name[numparpos] = paste0(puncttypeL,name[numparpos],puncttypeR)
  name[numparpos] = gsub("\\\\","",name[numparpos])
  name = paste(name,collapse=" ")

  return(name)
}

#' Make forest plot to summarize strata fits
#'
#' Make a forest plot to summarize model fits within each strata
#'
#@inheritParams run5STAR
#' @param MRSSmat matrix summary of model (e.g., cox, glm) fits within each strata
#' with required columns "exp.bhat", "exp.ci.lower.", and "exp.ci.upper." (e.g.,
#' HR/OR and corresponding confidence interval) and "weight.adap"
#' for family = "cox" or "binomial",
#' measure="OR", or columns "ratediff", "ci.lower", "ci.upper" for
#' family = "binomial", measure="RD", or columns "bhat", "ci.lower",
#' "ci.upper" for family = "gaussian"
#' @param MRSSres Row matrix with amalgamated results summary, including estimate
#' and corresponding confidence intervals
#' @param stratafits \code{\link[survival]{survfit}} object with number of
#' subjects, median survival times, etc. within each strata
#' @param family Trait family, current options: "cox", "binomial, "gaussian"
#' @param measure Response of interest; current options
#' are: for survival traits: "HR" (hazard ratio) and "TR"
#' (time ratio from model averaging of AFT models);
#' for binary traits: "RD" (risk difference using Miettinen and Nurminen
#' 1985 method) and "OR" (odds ratio, using GLM model fit);
#' ignored for continuous traits
#' @param labelNames Vector of names for each tree strata formed (e.g., number of
#' strata or covariate definition)
#' @param cilevel Confidence level alpha for overall result and confidence
#' intervals
#' @param alternative Whether to perform two-tailed ("two.tailed) or
#'  one-tailed ("less" or "greater") test
#' @param descfit Matrix summarizing descriptive model fit information (e.g.,
#' model averaged AFT estimate, CI, and Pr(TR > 1))
#' @param descRes Amalgamated results for descriptive (supplemental) model, e.g.,
#' as output from \code{\link{minPadapHR}}
#' @param txtgp Additional forestplot graphical parameters to pass to the
#' function (see also \code{\link[forestplot]{fpTxtGp}})
#' @param wrapwidth Number of characters after which to wrap forest plot label
#' names to next line (to avoid forest plot getting too long)
#' @param ... Optional additional arguments passed into the forestplot function
#' (see also \code{\link[forestplot]{forestplot}} )
#'
#' @return fplot An annotated forest plot showing by-stratum estimates and
#'  confidence intervals, as well as amalgamated result
#' @import forestplot
#' @import grid
strataforestplot = function(MRSSmat,MRSSres,stratafits,family,measure,
                            treetype="final",labelNames=NULL,cilevel,alternative,
                            descfit=NULL,descRes=NULL,fplottype=NULL,
                            txtgp=forestplot::fpTxtGp(
                              cex=0.9,xlab=grid::gpar(cex=0.9),
                              ticks=grid::gpar(cex=0.9),
                              label=grid::gpar(lineheight=0.75)),wrapwidth=70,...){

  avgname <- ifelse(treetype=='prespecified','2-STAR Average','5-STAR Average')

  #extracting summary of strata-level fits, as matrix/data frame
  if (family == "cox"){
    stratafittable = summary(stratafits)$table
    if (is.null(dim(stratafittable))) stratafittable = t(data.frame(stratafittable))
    nn = stratafits$n
  } else nn = stratafits[[1]]

  #defining/cleaning label names for plot
  nstrata = nrow(MRSSmat)
  if (is.null(labelNames)){
    labelNames = paste0("S",1:nstrata)
    if (treetype == "preliminary") labelNames = paste0("p",labelNames)
  }
  labelNames = sapply(labelNames,function(x) roundLabelNums(x,2))
  labelNames = sapply(labelNames,function(x) gsub("\\|","\\|\n",x))
  labelNames = sapply(labelNames,function(x)
    paste(strwrap(x, wrapwidth), collapse = "\n"))
  nlines = max(sapply(labelNames,function(x) stringr::str_count(x, "\n")))

  #annotation text to print on forest plot
  tabletext = paste0(nn," (",printR(MRSSmat$weight.SS*100,1),")")

  tabletext.all = cbind(c(paste0(toupper(substr(treetype,1,1)),
                                 substr(treetype,2,100),
                                 " Strata"),labelNames,avgname),
                        c("No. Subjects (%)",tabletext,
                          paste0(sum(stratafits$n)," (100)")))

  confintlevel = ifelse(alternative=="two.sided",cilevel,2*cilevel)

  #============================================================================#

  #(e.g., where log variable is of interest)
  if ((family == "cox" & measure %in% c("HR","TR"))|
      (family=="binomial" & measure=="OR")){

    nevents = stratafittable[,"events"]
    tabletext.events = paste0(nevents," (",
                              printR((nevents/sum(nevents))*100,1),")")
    tabletext.weights = MRSSmat$weight.adap

    measurename = measure

    tabletext.ci = c(paste0(printR(MRSSmat$exp.bhat.,2)," (",
                            printR(MRSSmat$exp.ci.lower.,2),
                            ", ",printR(MRSSmat$exp.ci.upper.,2),
                            ")"),
                     paste0(printR(MRSSres["exp(beta)"],2),
                            " (",printR(MRSSres["exp(ci lower)"],2),
                            ", ",printR(MRSSres["exp(ci upper)"],2),")") )

    #tabletext.prlt0 = c(printR(MRSSmat[,"Pr.beta.0."]*100,1),"")
    tabletext.prlt0 = c(printR(MRSSmat[,"Pr.beta.0."],3),"")
    tabletext.prlt0[tabletext.prlt0=="1.000"] = ">0.999"
    tabletext.prlt0[tabletext.prlt0=="0.000"] = "<0.001"
                        #printR(MRSSres["Pr(beta<0)"]*100,1))

    # weightName = ifelse(treetype=="prespecified","Inv.Var Weight %",
    #                     "Adap. Wt %")
    weightName <- "Adap. Wt %"

    #---------------------------------------------------#

    gtlt1 = ifelse(alternative=="greater",">1)","<1)")
    ltgt1 = ifelse(alternative=="greater","<1)",">1)")

    if (measure=="HR" & !is.null(descfit)){

      #time ratio information (edit from here!)
      tabletext.tr = c(paste0(printR(descfit[,"exp(bhat)"],2)," (",
                              printR(exp(descfit[,"ci.lower"]),2),", ",
                              printR(exp(descfit[,"ci.upper"]),2),")"),
                       paste0(printR(descRes["exp(beta)"],2)," (",
                              printR(descRes["exp(ci lower)"],2),", ",
                              printR(descRes["exp(ci upper)"],2),")"))

      # when measure is HR, descfit is TR -> want opposite sign from HR direction
      # so if alternative == 'less' or 'two.sided', assume want Pr(beta<0) so
      #           would want Pr(beta>0) for TR
      # if alternative == 'greater', assume want Pr(beta>0) for HR so want
      #           Pr(beta<0) for TR
      tabletext.trgt0 = c(printR(descfit[,ifelse(alternative=="greater",
                                                 "Pr(beta<0)","Pr(beta>0)")],3),"")

      tabletext.trgt0[tabletext.trgt0=="1.000"] = ">0.999"
      tabletext.trgt0[tabletext.trgt0=="0.000"] = "<0.001"

      tabletext.all = list(c(list(paste0(toupper(substr(treetype,1,1)),
                                         substr(treetype,2,100),
                                         " Strata")),labelNames,avgname),
                           c(list("No. Subjects (%)"),tabletext,
                             paste0(sum(stratafits$n)," (100)")),
                           c(list("No. Events (%)"),tabletext.events,
                             paste0(sum(stratafits$n.event)," (100)")),
                           c(list(paste0("Est. HR (",(1-confintlevel)*100,
                                         "% CI)")),tabletext.ci),
                           c(list(paste0("Pr(HR",gtlt1)),tabletext.prlt0),
                           c(list(paste0("Est. TR (",(1-confintlevel)*100,
                                         "% CI)")),tabletext.tr),
                           c(list(paste0("Pr(TR",ltgt1)),tabletext.trgt0))

    } else {
      if (measure == "TR"){

        tabletext.all = list(c(list(paste0(toupper(substr(treetype,1,1)),
                                           substr(treetype,2,100),
                                           " Strata")),labelNames,avgname),
                             c(list("No. Subjects (%)"),tabletext,
                               paste0(sum(stratafits$n)," (100)")),
                             c(list("No. Events (%)"),tabletext.events,
                               paste0(sum(stratafits$n.event)," (100)")),
                             #c(list(weightName),printR(tabletext.weights*100,1),100),
                             c(list(paste0("Est. TR (",(1-confintlevel)*100,
                                           "% CI)")),tabletext.ci),
                             c(list(paste0("Pr(TR",gtlt1)),tabletext.prlt0))

        if (!is.null(descfit)){

          #time ratio information (edit from here!)
          tabletext.hr = c(paste0(printR(descfit[,"exp(bhat)"],2)," (",
                                  printR(exp(descfit[,"ci.lower"]),2),", ",
                                  printR(exp(descfit[,"ci.upper"]),2),")"),
                           ifelse(is.null(descRes),"",paste0(
                             printR(descRes["exp(beta)"],2)," (",
                             printR(descRes["exp(ci lower)"],2),", ",
                             printR(descRes["exp(ci upper)"],2),")")))

          # when measure is TR, descfit is HR -> want opposite sign from TR direction
          # so if alternative == 'greater' or 'two.sided', assume want Pr(beta>0) so
          #           would want Pr(beta<0) for HR
          # if alternative == 'less', assume want Pr(beta<0) for TR so want
          #           Pr(beta>0) for HR
          tabletext.hrgt0 = c(printR(descfit[,ifelse(alternative=="less",
                                                     "Pr(beta>0)","Pr(beta<0)")],3),"")
          tabletext.hrgt0[tabletext.hrgt0=="1.000"] = ">0.999"
          tabletext.hrgt0[tabletext.hrgt0=="0.000"] = "<0.001"

          tabletext.all = list(c(list(paste0(toupper(substr(treetype,1,1)),
                                             substr(treetype,2,100),
                                             " Strata")),labelNames,avgname),
                               c(list("No. Subjects (%)"),tabletext,
                                 paste0(sum(stratafits$n)," (100)")),
                               c(list("No. Events (%)"),tabletext.events,
                                 paste0(sum(stratafits$n.event)," (100)")),
                               c(list(paste0("Est. TR (",(1-confintlevel)*100,
                                             "% CI)")),tabletext.ci),
                               c(list(paste0("Pr(TR",gtlt1)),tabletext.prlt0),
                               c(list(paste0("Est. HR (",(1-confintlevel)*100,
                                             "% CI)")),tabletext.hr),
                               c(list(paste0("Pr(HR",ltgt1)),tabletext.hrgt0))

        }

      } else {
      tabletext.all = list(c(list(paste0(toupper(substr(treetype,1,1)),
                                         substr(treetype,2,100),
                                         " Strata")),labelNames,avgname),
                           c(list("No. Subjects (%)"),tabletext,
                             paste0(sum(stratafits$n)," (100)")),
                           c(list("No. Events (%)"),tabletext.events,
                             paste0(sum(stratafits$n.event)," (100)")),
                           #c(list(weightName),printR(tabletext.weights*100,1),100),
                           c(list(as.expression(
                             bquote(hat(theta) * bold(" (") *
                                      bold(.((1-confintlevel)*100)) *
                                      bold(" % CI)"))) ),tabletext.ci),
                           c(list(as.expression(bquote(bold("Pr(") * bold(theta)*
                                                         bold(.(gtlt1))))),
                             tabletext.prlt0))
      }
    }


    #boxsize: mod from forestplot function
    upper = c(MRSSmat[,"exp.ci.upper."],MRSSres["exp(ci upper)"])
    lower = c(MRSSmat[,"exp.ci.lower."],MRSSres["exp(ci lower)"])
    cwidth = (upper-lower)
    # Set cwidth to min value if the value is invalid
    # this can be the case for reference points
    cwidth[cwidth <= 0 | is.na(cwidth)] <- min(cwidth[cwidth > 0])
    textHeight = 0.8
    info = c(MRSSmat[,"weight.adap"],1)
    info <- info/max(info, na.rm = TRUE)
    # Adjust the dots as it gets ridiculous with small text and huge dots
    if (any(textHeight*(nstrata+.5) * 1.5 < info))
      info <- textHeight*(nstrata+.5) * 1.5 * info/max(info, na.rm=TRUE) +
      textHeight*(nstrata+.5)*1.5/4
    bxsize = 0.2

    clippts = c(-Inf,Inf)
    if (is.infinite(min(lower))) clippts[1] = min(lower[!is.infinite(lower)])-0.1
    if (is.infinite(max(upper))|(max(upper)>1e4)){
      #clippts[2] = max(upper[!is.infinite(upper)])+0.1
      clippts[2] = max(upper[upper < 1e4])+0.1
    }
    clippts[2] = min(clippts[2],1.5*max(MRSSmat[,"exp.bhat."]))

    tickbylist <- c(0.1,0.2,0.5,1)
    clipdiff <- (min(max(upper),clippts[2]) - max(min(lower),clippts[1]))/10
    tickby <- tickbylist[which.min(abs(tickbylist - clipdiff))]
    # xtickpts = seq(plyr::round_any(max(min(lower),clippts[1]),0.2),#,floor),
    #                plyr::round_any(min(max(upper),clippts[2]),0.2),0.2)#,ceiling),0.2)
    xtickpts = seq(plyr::round_any(max(min(lower),clippts[1]),tickby),#,floor),
                   plyr::round_any(min(max(upper),clippts[2]),tickby),tickby)#,ceiling),0.2)


    nc = length(tabletext.all)-1
    xlabname=ifelse(family=="cox","Hazard Ratio (Test / Control)","Odds Ratio")
    if (measure == "TR") xlabname = "Time Ratio (Test / Control)"
    forestplot.default(labeltext=tabletext.all,
                           align=c("l",rep("c",nc)),#"c","c","c","c","c","c"),
                           graph.pos=4,graphwidth=grid::unit(3.5,"inches"),
                           mean=c(NA,MRSSmat[,"exp.bhat."],MRSSres["exp(beta)"]),
                           lower=c(NA,lower),upper=c(NA,upper),
                           is.summary=c(TRUE,rep(FALSE,nstrata),TRUE),
                           xlab=xlabname,zero=1,
                           col=forestplot::fpColors(
                             lines="darkblue",box="darkblue",summary=c("blue")),
                           lwd.zero=grid::gpar(lwd=2),#xticks.digits = 1,
                           fn.ci_sum = colfn,clip=clippts,xticks=xtickpts,
                           colgap=grid::unit(4,"mm"),boxsize=bxsize,
                           lineheight=grid::unit(1.4+0.3*max(0,nlines-3),"cm"),
                           txt_gp=txtgp,...)

    fplot = recordPlot()

    #============================================================================#

  } else if ((family == "gaussian" & measure == "MD")|
             (family=="binomial" & measure=="RD")|
             (family=="cox" & measure=="RMST")){

    tabletext.weights = MRSSmat$weight.adap
    measurename = measure
    if (measure == "RMST") measurename = "RMST Diff."
    if (measure == "RD") measurename = "Risk Diff."
    if (measure == "MD") measurename = "Mean Diff."
    measurenamelong = paste0(gsub("\\.","",measurename),"erence (Test - Control)")

    tabletext.ci = c(paste0(printR(MRSSmat$bhat,2)," (",
                            printR(MRSSmat$ci.lower,2),
                            ", ",printR(MRSSmat$ci.upper,2),")"),
                     paste0(printR(MRSSres[1],2),
                            " (",printR(MRSSres["ci lower"],2),
                            ", ",printR(MRSSres["ci upper"],2),")") )

    weightName = ifelse(treetype=="prespecified","Inv.Var Weight %",
                        "Adap. Weight %")

    tabletext.prltgt0 = c(printR(MRSSmat[,"Pr.beta.0."],3),"")
    tabletext.prltgt0[tabletext.prltgt0=="1.000"] = ">0.999"
    tabletext.prltgt0[tabletext.prltgt0=="0.000"] = "<0.001"

    ltgt0name = paste0("Pr(",measure,ifelse(alternative=="greater",">0)","<0)"))

    tabletext.all = cbind(c(paste0(toupper(substr(treetype,1,1)),
                                   substr(treetype,2,100),
                                   " Strata"),labelNames,avgname),
                          c("No. Subjects (%)",tabletext,
                            paste0(sum(nn)," (100)")),
                          #c(weightName,printR(tabletext.weights*100,1),100),
                          c(paste0(measurename," (",(1-confintlevel)*100,"% CI)"),
                            tabletext.ci),
                          c(ltgt0name,tabletext.prltgt0))

    #adding information on cases if trait is binary
    if (family == "binomial"){

      tabletext.cases = c("No. Cases (%)",paste0(
        stratafits$cj," (",printR(stratafits$cj/sum(stratafits$cj)*100,1),")"),
        paste0(sum(stratafits$cj)," (100)"))

      tabletext.subjA = c("N_A (%)",paste0(
        stratafits$njA," (",printR(stratafits$njA/sum(stratafits$njA)*100,1),")"),
        paste0(sum(stratafits$njA)," (100)"))
      tabletext.subjB = c("N_B (%)",paste0(
        stratafits$njB," (",printR(stratafits$njB/sum(stratafits$njB)*100,1),")"),
        paste0(sum(stratafits$njB)," (100)"))

      tabletext.caseA = c("Ncs_A (%)",paste0(
        stratafits$cjA," (",printR(stratafits$cjA/sum(stratafits$cjA)*100,1),")"),
        paste0(sum(stratafits$cjA)," (100)"))
      tabletext.caseB = c("Ncs_B (%)",paste0(
        stratafits$cjB," (",printR(stratafits$cjB/sum(stratafits$cjB)*100,1),")"),
        paste0(sum(stratafits$cjB)," (100)"))

      tabletext.frac = c("Cases/Subjs (% Subj)",
                         paste0(stratafits$cj," / ",stratafits$nj," (",
                                printR(stratafits$nj/sum(stratafits$nj)*100,1),
                                ")"),paste0(sum(stratafits$cj)," / ",sum(stratafits$nj),
                                                " (100)"))

      if (fplottype == "bytrt"){
        tabletext.all = cbind(tabletext.all[,1],tabletext.subjA,tabletext.subjB,
                              tabletext.caseA,tabletext.caseB,tabletext.all[,3:4])
      } else if (fplottype == "fractional"){
        tabletext.all = cbind(tabletext.all[,1],tabletext.frac,tabletext.all[,3:4])
      } else {
        tabletext.all = cbind(tabletext.all[,1:2],tabletext.cases,tabletext.all[,3:4])
      }

    }

    #boxsize: mod from forestplot function
    upper = c(MRSSmat[,"ci.upper"],MRSSres["ci upper"])
    lower = c(MRSSmat[,"ci.lower"],MRSSres["ci lower"])
    cwidth = (upper-lower)
    # Set cwidth to min value if the value is invalid
    # this can be the case for reference points
    cwidth[cwidth <= 0 | is.na(cwidth)] <- min(cwidth[cwidth > 0])
    textHeight = 0.8
    info = c(MRSSmat[,"weight.adap"],1)
    info <- info/max(info, na.rm = TRUE)
    # Adjust the dots as it gets ridiculous with small text and huge dots
    if (any(textHeight*(nstrata+.5) * 1.5 < info))
      info <- textHeight*(nstrata+.5) * 1.5 * info/max(info, na.rm=TRUE) +
      textHeight*(nstrata+.5)*1.5/4
    bxsize = 0.2

    clippts = c(-Inf,Inf)
    if (is.infinite(min(lower))) clippts[1] = min(lower[!is.infinite(lower)])-0.1
    if (is.infinite(max(upper))|(max(upper)>1e4)){
      #clippts[2] = max(upper[!is.infinite(upper)])+0.1
      clippts[2] = max(upper[upper < 1e4])+0.1
    }

    maxmin <- max(min(lower),clippts[1])
    minmax <- min(max(upper),clippts[2])
    rangefplot <- minmax - maxmin
    stepsizefplot <- round(rangefplot/5,1)
    if (stepsizefplot == 0) stepsizefplot <- round(rangefplot/5,2)
    #stepoptions <- c(0.1, 0.5, 1, 5, 10)
    #stepsizefplot_clean <- stepoptions[which.min(abs(stepoptions-stepsizefplot))]

    xtickpts = seq(plyr::round_any(max(min(lower),clippts[1]),0.1),
                   plyr::round_any(min(max(upper),clippts[2]),0.1),stepsizefplot)#,0.1)

    xtickpts <- sort(c(xtickpts,0))

    ncol = ncol(tabletext.all)
    xlabname=measurenamelong#ifelse(family=="cox","RMST Difference","Risk Difference")
    forestplot::forestplot(
      #labeltext=tabletext.all,align=c("l","c","c","c"),graph.pos=4,
      labeltext=tabletext.all,align=c("l",rep("c",ncol-1)),graph.pos=ncol-1,
      graphwidth=grid::unit(3.5,"inches"),
      mean=c(NA,MRSSmat[,1],MRSSres[1]),
      lower=c(NA,MRSSmat[,"ci.lower"],MRSSres["ci lower"]),
      upper=c(NA,MRSSmat[,"ci.upper"],MRSSres["ci upper"]),
      is.summary=c(TRUE,rep(FALSE,nstrata),TRUE),xlab=xlabname,zero=0,
      col=forestplot::fpColors(
        lines="darkblue",box="darkblue",summary=c("blue")),fn.ci_sum = colfn,
      lineheight=grid::unit(1.4+0.3*max(0,nlines-3),"cm"), #grid::unit(1.4,"cm"),
      lwd.zero=grid::gpar(lwd=2),boxsize=bxsize,clip=clippts,xticks=xtickpts,
      txt_gp=txtgp,...)

    fplot = recordPlot()

  }

  return(fplot)

}


################################################################################

##===========================================================##
## Fit Linear or Logistic Regression GLM Model Within Strata ##
##===========================================================##

#' Fit Linear or Logistic Regression GLM Model Within Strata
#'
#' Fits a linear or logistic regression GLM model within each formed strata,
#' providing estimate of regression coefficient or log odds ratio (EXPERIMENTAL!)
#'
#' @param yy      Continuous response or case control status response vector
#' @param arm     Treatment indicator, 1 = test treatment, 0 = control
#' @param family  Family for glm; options: "binomial" or "gaussian"
#' @param treeStrata Vector of strata membership, as defined by ctree strata
#' @param termNodes Vector of names for each tree strata formed (e.g., number of
#' strata or covariate definition)
#' @param cilevel Confidence level alpha for overall result and confidence
#' intervals
#' @param verbose Numeric variable indicating amount of information to print
#' to the terminal (0 = nothing, 1 = notes only, 2 = notes and intermediate output)
#'
#' @return \itemize{
#'     \item fitsummary: summary of GLM fits within each strata
#'     \item table: Summary of amalgamated regression coefficient or odds ratio
#'      estimate, variance, p-value, and cilevelx100\% confidence interval
#'      assuming sample size weights
#'     \item weights: Sample size weights used to construct estimate
#' }
glmbystrata = function(yy,arm,family,treeStrata,termNodes=NULL,cilevel=0.05,
                       verbose=0,alternative="two.sided"){

  dat = data.frame(yy,arm)
  nstrata = max(treeStrata,na.rm=TRUE)
  if (is.null(termNodes)) termNodes = unique(treeStrata)

  #fitting within each formed strata
  glmCtreeStrataFit = lapply(1:nstrata,function(x){
    glm(yy ~ as.factor(arm),subset=treeStrata==x,data=dat,
        family=family)#binomial(link=logit))
  })

  #pulling out relevant summary statistics from glm fits to calculate glmSS
  glmMat = matrix(unlist(lapply(glmCtreeStrataFit,function(x){
    summary(x)$coef[2,1:2]
  } )),ncol=2,byrow=TRUE)
  glmMat[,2] = glmMat[,2]^2

  cnj = sapply(1:nstrata,function(x) sum(treeStrata==x))
  glmMat = cbind(glmMat,cnj)

  colnames(glmMat) = c("bhatj","vhatj","nj")

  #calculating glmSS
  glmbeta = glmMat[,"bhatj"]
  glmvar = glmMat[,"vhatj"]

  glmwSS = cnj/sum(cnj)
  glmSS = sum(glmbeta*glmwSS)
  ORglmSS = exp(glmSS)
  VglmSS = sum(glmwSS^2*glmvar)
  TglmSS = glmSS/sqrt(VglmSS)
  pvglmSS = 2*pnorm(abs(TglmSS),lower.tail=FALSE,log.p=FALSE)

  #two-sided test
  if (alternative == "two.sided"){
    qcrit = qnorm(1-cilevel/2,0,1)
  } else qcrit = qnorm(1-cilevel,0,1)

  #overall and within-strata cis
  glmci = cbind(glmSS - qcrit*sqrt(VglmSS),glmSS + qcrit*sqrt(VglmSS))
  expglmci = exp(glmci)
  glmstratci = cbind(glmbeta - qcrit*sqrt(glmvar),
                     glmbeta + qcrit*sqrt(glmvar))

  glmstratpv = 2*pnorm(abs(glmbeta/sqrt(glmvar)),0,1,lower.tail=FALSE)

  #compiling within-strata effect estimate summary matrix
  if (family == "binomial"){

    glmMRSSmat = cbind(glmbeta,glmvar,glmstratci,exp(glmbeta),exp(glmstratci),
                       glmstratpv)
    colnames(glmMRSSmat) = c("bhat","v(bhat)","ci.lower","ci.upper","exp(bhat)",
                             "exp(ci.lower)","exp(ci.upper)","pval")

    glmMRSSres = list(table=data.frame(logorSS=glmSS,vSS=VglmSS,ciSS=glmci,
                                       orSS=ORglmSS,expciSS=expglmci,
                                       pvSS=pvglmSS),weightSS=glmwSS)

  } else if (family == "gaussian"){

    #(no need for exp() terms for quantitative traits)
    glmMRSSmat = cbind(glmbeta,glmvar,glmstratci,glmstratpv)
    colnames(glmMRSSmat) = c("bhat","v(bhat)","ci.lower","ci.upper","pval")

    glmMRSSres = list(table=data.frame(betaSS=glmSS,vSS=VglmSS,ciSS=glmci,
                                       pvSS=pvglmSS),weightSS=glmwSS)
  } else stop("Family must be 'binomial' or 'gaussian'.")

  glmMRSSmat.pt2 = matrix(unlist(glmwSS),nrow=nstrata,byrow=FALSE)
  colnames(glmMRSSmat.pt2) = c("weight.SS")
  glmMRSSmat = cbind(glmMRSSmat,glmMRSSmat.pt2)

  if (verbose > 1){
    print("GLM fit within each strata")
    print(round(glmMRSSmat,4))
  }

  #============================================================================#

  return(list(fitsummary=glmMRSSmat,table=glmMRSSres$table,
              weights=glmMRSSres$weightSS))
}

################################################################################

##========================================##
## Estimate Mean Difference Within Strata ##
##========================================##

#' Estimate mean difference Within Strata
#'
#' Estimates difference in means and significance within each formed strata
#'
#' @param yy      Continuous response vector
#' @param arm     Treatment indicator, 1 = test treatment, 0 = control
#' @param treeStrata Vector of strata membership, as defined by ctree strata
#' @param treetype String, whether trees input are "preliminary" or "final"
#' @param termNodes Vector of names for each tree strata formed (e.g., number of
#' strata or covariate definition)
#' @param cilevel Confidence level alpha for overall result and confidence
#' intervals
#' @param verbose Numeric variable indicating amount of information to print
#' to the terminal (0 = nothing, 1 = notes only, 2 = notes and intermediate output)
#' @param plot    Logical, whether to create risk difference forest plots
#'
#' @return \itemize{
#'     \item fitsummary: summary of mean difference within each strata
#'     \item table: Summary of amalgamated mean difference estimate, variance,
#'     p-value, and cilevelx100\% confidence interval assuming sample size weights
#'     \item weights: Sample size weights used to construct estimate
#' }
mdbystrata = function(yy,arm,treeStrata,treetype='final',termNodes=NULL,cilevel=0.05,
                       verbose=0,alternative="two.sided",plot=FALSE){

  dat = data.frame(yy,arm)
  nstrata = max(treeStrata,na.rm=TRUE)
  if (is.null(termNodes)) termNodes = unique(treeStrata)

  # if (measure == "MD"){

    #performing t-test within each stratum
    #assume unequal variance per group
    #setting levels to ensure always subtract arm 1 - arm 0 (vs 0 - 1 by default)
    mdCtreeStrataFit = lapply(1:nstrata,function(x){
      t.test(yy ~ factor(arm, levels = c(1,0)), alternative = alternative,
             var.equal = FALSE, conf.level = 1-cilevel, subset = treeStrata==x,
             data = dat)
    })

    if (verbose > 2) print(mdCtreeStrataFit)

    mdtestMat = matrix(unlist(lapply(mdCtreeStrataFit, function(x){
      meandiff = x$estimate["mean in group 1"] - x$estimate["mean in group 0"]
      stderr = x$stderr
      c(meandiff, stderr)
    })), ncol = 2, byrow = TRUE)
    mdtestMat[,2] = mdtestMat[,2]^2

  # } else if (measure == "MDnp"){
  #
  #   #nonparametric mean difference from Wilcoxon rnak sum test
  #   mdCtreeStrataFit = lapply(1:nstrata,function(x){
  #     wilcox.test(yy ~ factor(arm, levels = c(1,0)), alternative = alternative,
  #                 mu = 0, paired = FALSE, conf.int = TRUE,
  #                 conf.level = 1-cilevel, subset = treeStrata==x, data = dat)
  #   })
  # }

  cnj = sapply(1:nstrata,function(x) sum(treeStrata==x))
  mdtestMat = cbind(mdtestMat,cnj)

  colnames(mdtestMat) = c("bhatj","vhatj","nj")

  #calculating glmSS
  mdbeta = mdtestMat[,"bhatj"]
  mdvar = mdtestMat[,"vhatj"]

  mdwSS = cnj/sum(cnj)
  mdSS = sum(mdbeta*mdwSS)
  #mdSS = exp(mdSS)
  VmdSS = sum(mdwSS^2*mdvar)
  TmdSS = mdSS/sqrt(VmdSS)

  if (alternative == "two.sided"){
    pvmdS = 2*pnorm(abs(TmdSS),lower.tail=FALSE,log.p=FALSE)
  } else if (alternative == "greater"){
    pvmdSS = pnorm(TmdSS,lower.tail=FALSE,log.p=FALSE)
  } else if (alternative == "less"){
    pvmdSS = pnorm(TmdSS,lower.tail=TRUE,log.p=FALSE)
  } else stop ("Alternative must be one of 'two.sided', 'greater', or 'less'.")

  #calculating probability each betahat is < 0 assuming N(betahat, Vihat) dist'n
  prlt0 = pnorm(0,mdbeta,sqrt(mdvar))
  prgt0 = pnorm(0,mdbeta,sqrt(mdvar),lower.tail=FALSE)

  #two-sided test
  if (alternative == "two.sided"){
    qcrit = qnorm(1-cilevel/2,0,1)
  } else qcrit = qnorm(1-cilevel,0,1)

  #overall ci
  mdci = cbind(mdSS - qcrit*sqrt(VmdSS),mdSS + qcrit*sqrt(VmdSS))

  #inverse variance weights for reference
  weight.invVar = (1/mdvar)/sum(1/mdvar)

  #amalgamated results using sample size weights
  mdMRSSres = list(table=data.frame(mdSS=mdSS,vSS=VmdSS,ciSS=mdci,
                                    TmdSS=TmdSS,pvSS=pvmdSS),
                    weightSS=mdwSS,weightIV = weight.invVar)

  #confidence intervals within each stratum
  mdstratci = cbind(mdbeta - qcrit*sqrt(mdvar),
                    mdbeta + qcrit*sqrt(mdvar))

  #test statistic and p-value for each stratum
  mdstratT = mdbeta/sqrt(mdvar)
  if (alternative == "two.sided"){
    mdstratpv = 2*pnorm(abs(mdstratT),0,1,lower.tail=FALSE)
  } else if (alternative == "greater"){
    mdstratpv = pnorm(mdstratT,0,1,lower.tail=FALSE)
  } else if (alternative == "less"){
    mdstratpv = pnorm(mdstratT,0,1,lower.tail=TRUE)
  } else stop ("Alternative must be one of 'two.sided', 'greater', or 'less'.")

  #summarized by-stratum results
  mdMRSSmat = cbind(mdbeta,mdvar,mdstratci,mdstratT,mdstratpv,
                    ifelse(rep(alternative=="greater",nstrata),prgt0,prlt0))
  colnames(mdMRSSmat) = c("bhat","v(bhat)","ci.lower","ci.upper","Zstat","pval",
                           ifelse(alternative=="greater","Pr(beta>0)",
                                  "Pr(beta<0)"))

  if (verbose > 2) print(mdMRSSmat)

  mdMRSSmat.pt2 = matrix(cbind(unlist(mdwSS),unlist(weight.invVar)),
                          nrow=nstrata,byrow=FALSE)
  colnames(mdMRSSmat.pt2) = c("weight.SS","weight.invVar")

  mdMRSSmat = cbind(mdMRSSmat,mdMRSSmat.pt2)

  if (verbose > 1){
    print("Esimated mean difference fit within each strata")
    print(round(mdMRSSmat,4))
  }

  if (plot == TRUE){

    getPalette = grDevices::colorRampPalette(RColorBrewer::brewer.pal(min(
      nstrata,11),"Spectral"))

    wrapper <- function(x, ...)
    {
      paste(strwrap(x, ...), collapse = "\n")
    }

    labelStart = "S"
    if (treetype == "preliminary") labelStart = "pS"
    if (treetype == "prespecified") labelStart = "S_"

    labelNames = paste0("S",1:nstrata)
    if (treetype == "preliminary") labelNames = paste0("p",labelNames)
    stratLabelName = ifelse(treetype=="preliminary","Preliminary Risk",
                            "Identified Risk")
    if (treetype=="prespecified"){
      labelNames = paste0("S_",1:nstrata)
      stratLabelName = "Design" #"Prespecified"
    }

    if (!is.null(termNodes)){
      labelNames = paste0(labelNames,": ",termNodes)
      labelNames = sapply(labelNames,function(x) wrapper(x,width=50))
      names(labelNames) = NULL
    }

    mdplotdat <- data.frame(yy, treeStrata, arm)

    mdplotdat$treeStrata <- as.factor(mdplotdat$treeStrata)
    mdplotdat$arm <- as.factor(mdplotdat$arm)
    mdplotdat$Treatment <- as.character(mdplotdat$arm)
    mdplotdat$Treatment[mdplotdat$Treatment == 0] <- 'Control'
    mdplotdat$Treatment[mdplotdat$Treatment == 1] <- 'Test'

    facetstrattitle <- list()
    for (x in 1:length(termNodes)){
      facetstrattitle[[x]] <- wrapper(paste0(labelStart,x," (",
                                             round(mdwSS[[x]]*100,1),
                                             "% of subjects)"),width=45)
    }

    strata_labeller <- function(variable,value){
      return(facetstrattitle[value])
    }

    p3 <- ggplot2::ggplot(mdplotdat) +
      ggplot2::geom_boxplot(aes(y = yy, x = Treatment, group = Treatment,
                                fill = Treatment)) +
      ggplot2::facet_grid(.~treeStrata, labeller=strata_labeller) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank())

    p4 <- ggplot2::ggplot(mdplotdat) +
      ggplot2::geom_boxplot(aes(y = yy, x = treeStrata, group = treeStrata,
                       fill = treeStrata)) +
      ggplot2::theme_bw() +
      ggplot2::scale_fill_manual(name=paste0(stratLabelName," Strata"),#" Risk Strata"),
                                  labels = labelNames,values = getPalette(nstrata)) +
      ggplot2::theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank())

  } else p3 <- p4 <- NULL

  #KM curves ignoring treatment arm indicator
  # Stratum = as.factor(treeStrata)
  # s <- list(sapply(unique(Stratum), function(x) sum(Stratum==x)))
  #need order to match order of the strata being output
  s <- list(cnj)
  #list(unname(table(Stratum)))

  #============================================================================#

  return(list(fitsummary=mdMRSSmat, table=mdMRSSres$table, stratafit=s,
              weights=mdMRSSres$weightSS, p3=p3, p4=p4))
}


################################################################################

##========================================##
## Estimate Risk Difference Within Strata ##
##========================================##

#' Estimate Risk Difference Within Strata
#'
#' Estimates risk difference and significance within each formed strata
#'
#' @param yy      Case control status response vector
#' @param arm     Treatment indicator, 1 = test treatment, 0 = control
#' @param treeStrata Vector of strata membership, as defined by ctree strata
#' @param treetype String, whether trees input are "preliminary" or "final"
#' @param termNodes Vector of names for each tree strata formed (e.g., number of
#' strata or covariate definition)
#' @param alternative For tests, whether alternative hypothesis is "less",
#'  "greater", or "two.sided" (default = "two.sided")
#' @param cilevel Confidence level alpha for overall result and confidence
#' intervals
#' @param verbose Numeric variable indicating amount of information to print
#' to the terminal (0 = nothing, 1 = notes only, 2 = notes and intermediate output)
#' @param plot    Logical, whether to create risk difference forest plots
#'
#' @return \itemize{
#'     \item fitsummary: summary of risk difference within each strata
#'     \item table: Summary of amalgamated risk difference estimate, variance,
#'     p-value, and cilevelx100\% confidence interval assuming sample size weights
#'     \item weights: Sample size weights used to construct estimate
#' }
ratediffbystrata = function(yy,arm,treeStrata,treetype="final",
                            termNodes=NULL,alternative="two.sided",cilevel=0.05,
                            verbose=0,plot=TRUE){

  dat = data.frame(yy,arm)
  nstrata = max(treeStrata,na.rm=TRUE)
  if (is.null(termNodes)) termNodes = unique(treeStrata)

  #alt to glm: risk difference test within each strata (Wald-based)
  rdCtreeStrataFit = lapply(1:nstrata, function(x) {

    #data for strata overall and strata for each treatment arm
    datx = dat[treeStrata == x,]
    datx0 = datx[datx$arm == 0,]
    datx1 = datx[datx$arm == 1,]

    #number of cases in each arm
    c0 = sum(datx0$yy == 1)
    c1 = sum(datx1$yy == 1)

    #number of subjects in each arm
    n0 = nrow(datx0)
    n1 = nrow(datx1)

    #proportion of cases in each arm
    p0 = c0/n0
    p1 = c1/n1

    #risk difference and corresponding wald-based variance (bias-adjusted)
    RD = p1 - p0
    varRD = p0*(1-p0)/(n0-1) + p1*(1-p1)/(n1-1)
    return(c(RD, varRD, n0+n1, n0, n1, c0+c1, c0, c1))

  })

  rdMat = matrix(unlist(rdCtreeStrataFit),ncol=8,byrow=TRUE)
  colnames(rdMat) = c("rdj","vhatj","nj", "njB", "njA", "cj", "cjB", "cjA")

  # #alt to glm: risk difference test within each strata (MN/score-based)
  # rdCtreeStrataFit = lapply(1:nstrata,function(x) {
  #   datx = dat[treeStrata==x,]
  #   c1 = sum(datx$yy[datx$arm==1])
  #   S1 = sum(datx$arm==1)
  #   c0 = sum(datx$yy[datx$arm==0])
  #   S0 = sum(datx$arm==0)
  #
  #   r1 = c1/S1; r0 = c0/S0
  #   c=c1+c0; S=S1+S0; r = c/S
  #
  #   #proper chisq stat for RD = 0 (R1 = R0):
  #   RDall=ratesci::scoreci(c1,S1,c0,S0,distrib="bin",measure="RD",level=1-cilevel,
  #                          skew=FALSE,weighting="MN")
  #   Ts0 = ( (r1-r0)^2 )/( r*(1-r)*(S/(S-1))*(1/S1+1/S0) )
  #
  #   #########calc variance
  #   RD = 0
  #   L0 = c0*RD*(1-RD)
  #   L1 = (S0*RD - S - 2*c0)*RD+c
  #   L2 = (S1+2*S0)*RD-S-c
  #   L3 = S
  #
  #   q = L2^3/(3*L3)^3 - L1*L2/(6*L3^2) + L0/(2*L3)
  #   sgn = sign(q)
  #   p = sgn*(L2^2/(3*L3)^2 - L1/(3*L3))^(1/2)
  #   if (sign(p) != sign(q)) p = -p
  #   a = (1/3)*(pi + acos(q/p^3))
  #   R0tilde = 2*p*cos(a) - L2/(3*L3)
  #   R1tilde = R0tilde + RD
  #
  #   VarRD = (R1tilde*(1-R1tilde)/S1 + R0tilde*(1-R0tilde)/S0)*(S/(S+1))
  #   ################################
  #   RDall$var = VarRD
  #   RDall
  #
  # })

  # #pulling out relevant summary statistics from glm fits to calculate glmSS
  # rdMat = matrix(unlist(lapply(rdCtreeStrataFit,function(x){
  #   cbind(x$estimates[,"MLE"],x$var)
  # } )),ncol=2,byrow=TRUE)
  # cnj = sapply(1:nstrata,function(x) sum(treeStrata==x))
  # rdMat = cbind(rdMat,cnj)
  # colnames(rdMat) = c("rdj","vhatj","nj")

  #calculating rdS
  rdest = rdMat[,"rdj"]
  rdvar = rdMat[,"vhatj"]
  cnj = rdMat[,"nj"]

  rdwSS = cnj/sum(cnj)
  rdSS = sum(rdest*rdwSS)
  rdSSp = rdSS*100
  VrdSS = sum(rdwSS^2*rdvar)
  TrdSS = rdSS/sqrt(VrdSS)

  if (alternative == "two.sided"){
    pvrdSS = 2*pnorm(abs(TrdSS),lower.tail=FALSE,log.p=FALSE)
  } else if (alternative == "greater"){
    pvrdSS = pnorm(TrdSS,lower.tail=FALSE,log.p=FALSE)
  } else if (alternative == "less"){
    pvrdSS = pnorm(TrdSS,lower.tail=TRUE,log.p=FALSE)
  } else stop ("Alternative must be one of 'two.sided', 'greater', or 'less'.")

  # #two-sided test
  # qcrit = qnorm(1-cilevel/2,0,1)
  #
  # #overall and within-strata cis
  # rdci = cbind(rdSS - qcrit*sqrt(VrdSS),rdSS + qcrit*sqrt(VrdSS))
  # rdcip = rdci*100
  # rdstratci2 = cbind(rdest - qcrit*sqrt(rdvar),
  #                    rdest + qcrit*sqrt(rdvar))
  # rdstratci = matrix(unlist(lapply(rdCtreeStrataFit,function(x){
  #   x$estimates[,c("Lower","Upper")]
  # })),ncol=2,byrow=TRUE)
  #
  # rdstratpv = sapply(rdCtreeStrataFit,function(x){
  #   x$pval[,"pval2sided"]
  # })
  # rdstratpv2 = 2*pnorm(abs(rdest/sqrt(rdvar)),0,1,lower.tail=FALSE)
  #
  # #compiling within-strata effect estimate summary matrix
  # rdMRSSmat = cbind(rdest,rdvar,rdstratci,rdstratpv)
  # colnames(rdMRSSmat) = c("ratediff","v(ratediff)","ci.lower","ci.upper","pval")
  #
  # rdMRSSres = list(table=data.frame(ratediffSS=rdSS,vSS=VrdSS,ciSS=rdci,
  #                                   pvSS=pvrdSS),weightSS=rdwSS)
  #
  # rdMRSSmat.pt2 = matrix(unlist(rdwSS),nrow=nstrata,byrow=FALSE)
  # colnames(rdMRSSmat.pt2) = c("weight.SS")
  # rdMRSSmat = cbind(rdMRSSmat,rdMRSSmat.pt2)
  # labelNames = paste0("S",1:nstrata)
  # if (treetype == "preliminary"){
  #   labelNames = paste0("p",labelNames)
  # }
  # rownames(rdMRSSmat) = labelNames
  #
  # if (verbose > 1){
  #   print("Estimated risk difference within each strata")
  #   print(round(rdMRSSmat,4))
  # }
  #
  # #============================================================================#
  #
  # return(list(fitsummary=rdMRSSmat,table=rdMRSSres$table,
  #             weights=rdMRSSres$weightSS))

  #calculating probability each betahat is < 0 assuming N(betahat, Vihat) dist'n
  prlt0 = pnorm(0,rdest,sqrt(rdvar))
  prgt0 = pnorm(0,rdest,sqrt(rdvar),lower.tail=FALSE)

  #two-sided test
  if (alternative == "two.sided"){
    qcrit = qnorm(1-cilevel/2,0,1)
  } else qcrit = qnorm(1-cilevel,0,1)

  #overall ci
  rdci = cbind(rdSS - qcrit*sqrt(VrdSS),rdSS + qcrit*sqrt(VrdSS))

  #inverse variance weights for reference
  weight.invVar = (1/rdvar)/sum(1/rdvar)

  #amalgamated results using sample size weights
  rdMRSSres = list(table=data.frame(rdSS=rdSS,vSS=VrdSS,ciSS=rdci,
                                    TrdSS=TrdSS,pvSS=pvrdSS),
                   weightSS=rdwSS,weightIV = weight.invVar)

  #confidence intervals within each stratum
  rdstratci = cbind(rdest - qcrit*sqrt(rdvar),
                    rdest + qcrit*sqrt(rdvar))

  #test statistic and p-value for each stratum
  rdstratT = rdest/sqrt(rdvar)
  if (alternative == "two.sided"){
    rdstratpv = 2*pnorm(abs(rdstratT),0,1,lower.tail=FALSE)
  } else if (alternative == "greater"){
    rdstratpv = pnorm(rdstratT,0,1,lower.tail=FALSE)
  } else if (alternative == "less"){
    rdstratpv = pnorm(rdstratT,0,1,lower.tail=TRUE)
  } else stop ("Alternative must be one of 'two.sided', 'greater', or 'less'.")

  #summarized by-stratum results
  rdMRSSmat = cbind(rdest,rdvar,rdstratci,rdstratT,rdstratpv,
                    ifelse(alternative==rep("greater",nstrata),prgt0,prlt0))
  colnames(rdMRSSmat) = c("bhat","v(bhat)","ci.lower","ci.upper","Zstat","pval",
                          ifelse(alternative=="greater","Pr(beta>0)",
                                 "Pr(beta<0)"))

  rdMRSSmat.pt2 = matrix(cbind(unlist(rdwSS),unlist(weight.invVar)),
                         nrow=nstrata,byrow=FALSE)
  colnames(rdMRSSmat.pt2) = c("weight.SS","weight.invVar")

  rdMRSSmat = cbind(rdMRSSmat,rdMRSSmat.pt2)

  if (verbose > 1){
    print("Esimated mean difference fit within each strata")
    print(round(rdMRSSmat,4))
  }

  #No. subjects and No. cases by strata overall, and by treatment assignment
  casesubjMat <- rdMat[,3:8]
  if (is.null(ncol(casesubjMat))){
    casesubjMat <- matrix(casesubjMat, nrow=1,
                          dimnames = list(1,colnames(rdMat)[3:8]))
  }
  s <- split(casesubjMat, rep(1:ncol(casesubjMat), each = nrow(casesubjMat)))
  names(s) <- colnames(casesubjMat)

  #============================================================================#

  return(list(fitsummary=rdMRSSmat, table=rdMRSSres$table, stratafit=s,
              weights=rdMRSSres$weightSS))

}

################################################################################

##===============================##
## Check Events Per Strata x Arm ##
##===============================##

#' Summarize events per strata x arm
#'
#' Summarizes number of events, subjects, and percent censoring within each
#' strata and treatment arm from formed risk strata
#'
#' @param treeStrata    Vector of subject-level strata membership
#' @param arm           Treatment indicator, 1 = test treatment, 0 = control
#' @param status        For family="cox", the censoring status (1 = event,
#' 0 = censored); for family="binomial" the case/control status
#' @param family        Trait family, current options: "cox" or "binomial"
#'
#' @return Table summarizing number of events, subjects, and percent censoring
#' (for family = "cox") or percent controls (for family = "binomial")
#' within each strata and treatment arm
eventsbystrata = function(treeStrata,arm,status,family="cox"){

  #number of strata
  nstrata = length(unique(treeStrata))

  #number of events per strata for control and treatment arms
  neventsPerTreeStrataControl = sapply(1:nstrata,function(x)
    sum(status[(treeStrata==x)&(arm==0)]))
  neventsPerTreeStrataTrt = sapply(1:nstrata,function(x)
    sum(status[(treeStrata==x)&(arm==1)]))

  #number of subjects per strata for control and treatment arms
  nsubjPerTreeStrataControl = sapply(1:nstrata,function(x)
    sum(treeStrata==x & arm == 0))
  nsubjPerTreeStrataTrt = sapply(1:nstrata,function(x)
    sum(treeStrata==x & arm == 1))

  #percent of subjects censored (or controls for binary traits) for control
  #and treatment arms
  pctCensorControl = 1-neventsPerTreeStrataControl/nsubjPerTreeStrataControl
  pctCensorTrt = 1-neventsPerTreeStrataTrt/nsubjPerTreeStrataTrt

  #summary of all above values
  perStrataEventSummary = matrix(c(neventsPerTreeStrataControl,
                                   neventsPerTreeStrataTrt,
                                   nsubjPerTreeStrataControl,
                                   nsubjPerTreeStrataTrt,
                                   pctCensorControl,pctCensorTrt),
                                 ncol=2*nstrata,byrow=TRUE)
  colnames(perStrataEventSummary) = paste(rep(paste0("S",1:nstrata),2),
                                          rep(c("Control","Test"),
                                              each=nstrata),sep="-")
  if (family == "cox"){
    rownames(perStrataEventSummary) = c("nevents","nsubj","pctcensor")
  } else if (family == "binomial"){
    rownames(perStrataEventSummary) = c("ncases","nsubj","pctcontrols")
  }

  return(perStrataEventSummary)

}


################################################################################
