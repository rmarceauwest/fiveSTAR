##============================##
## Perform adaptive weighting ##
##============================##

#' Estimate amalgamated treatment effect
#'
#' Estimates overall treatment effect amalgamated from adaptively weighted
#' strata-level treatment effects
#'
#' @param weights Original, non-adaptive weights (e.g., sample size weights)
#' @param betas   Vector of estimated treatment effect within each strata
#' @param vars    Vector of estimated variances of betas
#' @param cilevel Significance level for adaptive weighted confidence intervals
#' (i.e., (1-cilevel)x100\% CIs)
#'
#' @return \itemize{
#'     \item adaptivewtRes: table of amalgamated log-hazard estimate, variance,
#'     cilevelx100\% ci,and p-value using the adaptive weight
#'     \item adaptivewts: adaptive weights used
#'     \item adaptivecorr: rank correlation between betas and vars used to determine
#'     where original (i.e., combining estimates) or new weights (i.e., combining
#'     test statistics) are used
#' }
adaptivewtHR = function(weights,betas,vars,cilevel){

  #whether to use the new weight or original input weight
  adaptivecorr = cor(betas,vars,method="kendall")
  isnewwt = (adaptivecorr > 0)

  #original "weights" of 1 if no strata found
  if (length(betas)==1) isnewwt = FALSE

  #calculate new weights
  if (isnewwt){
    newwts = (weights/sqrt(vars))/(sum(weights/sqrt(vars)))
  } else newwts = weights

  #overall treatment effect estimate, variance, and p-value
  beta.newwt = sum(newwts*betas)
  var.newwt = sum(newwts^2*vars)
  pv.newwt = 2*pnorm(abs(beta.newwt/sqrt(var.newwt)),lower.tail=FALSE)

  qcrit = qnorm(1-cilevel/2,0,1)
  ci.newwt = cbind(beta.newwt - qcrit*sqrt(var.newwt),
                   beta.newwt + qcrit*sqrt(var.newwt))

  #collecting relevant results for return
  newwtRes = c(beta.newwt,var.newwt,ci.newwt,exp(beta.newwt),exp(ci.newwt),
               pv.newwt)
  names(newwtRes) = c("beta","var","ci lower","ci upper","exp(beta)",
                      "exp(ci lower)","exp(ci upper)","pv")

  return(list(adaptivewtRes=newwtRes,adaptivewts=newwts,adaptivecorr=adaptivecorr))
}

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
#' @param treeStrata Vector of strata membership, as defined by ctree strata
#' @param termNodes Vector of names for each tree strata formed (e.g., number of
#' strata or covariate definition)
#' @param treetype String, whether trees input are "preliminary" (e.g.,
#' from step 3A) or "final" (e.g., from step 3B)
#' @param cilevel Confidence level alpha for overall result and confidence
#' intervals (i.e., produces (1-cilevel)x100\% CIs)
#' @param verbose Numeric variable indicating amount of information to print
#' to the terminal (0 = nothing, 1 = notes only, 2 = notes and intermediate output)
#' @param plot    Logical, whether to create within strata and pooled between
#'  strata KM plots
#' @param timeunit Optional argument, time unit for survival data
#' (e.g., Months, Years,..); currently only used for plots and ignored if
#' plot is FALSE
#'
#' @return \itemize{
#'     \item fitsummary: summary of cox fits within each strata, along with
#'     Grambsch-Therneau (GT) tests for PH within each strata
#'     \item stratafit: list of full coxph fits within each strata
#'     \item table: Summary of amalgamated hazard ratio estimate, variance,
#'     p-value, and cilevelx100\% confidence interval assuming sample size weights
#'     \item weights: Sample size weights used to construct estimate
#'     \item bystrataKM: Kaplan-Meier survival curve plots within each strata,
#'     returned if plot = TRUE
#'     \item betweenstrataKM: Kaplan-Meier survival curve plots from pooled
#'      treatment assignment data from each strata, returned if plot = TRUE
#' }
#' @import survival
#'
coxbystrata = function(time,status,arm,treeStrata,termNodes=NULL,
                       treetype="final",cilevel=0.05,verbose=0,plot=TRUE,
                       timeunit=NULL){

  dat = data.frame(time,status,arm)
  nstrata = max(treeStrata,na.rm=TRUE)
  if (is.null(termNodes)) termNodes = unique(treeStrata)

  #fitting within each formed strata
  coxCtreeStrataFit = lapply(1:nstrata,function(x){
    coxph(Surv(time,status) ~ as.factor(arm),subset=treeStrata==x,data=dat)
  })

  #calculating GT test for PH within each strata
  pv.GTzph = sapply(1:nstrata,function(x)
    cox.zph(coxCtreeStrataFit[[x]])[[1]][,"p"])

  #pulling out relevant summary statistics from cox fits to calculate coxSS
  coxMat = matrix(unlist(lapply(coxCtreeStrataFit,function(x){
    unlist(x[c("coefficients","var","n")])
  } )),ncol=3,byrow=TRUE)
  colnames(coxMat) = c("bhatj","vhatj","nj")

  cnj = coxMat[,"nj"]
  coxbeta = coxMat[,"bhatj"]
  coxvar = coxMat[,"vhatj"]
  coxwSS = cnj/sum(cnj)

  #calculating coxSS (amalgamated HR estimate using sample size weights)
  #and corresponding variance, test statistic, and p-value
  coxSS = sum(coxbeta*coxwSS)
  HRcoxSS = exp(coxSS)
  VcoxSS = sum(coxwSS^2*coxvar)
  TcoxSS = coxSS/sqrt(VcoxSS)
  pvcoxSS = 2*pnorm(abs(TcoxSS),lower.tail=FALSE,log.p=FALSE)

  #calculating confidence intervals for estimated logHR and exponentiating
  #for estimated amalgamated hazard ratio
  qcrit = qnorm(1-cilevel/2,0,1)
  coxci = cbind(coxSS - qcrit*sqrt(VcoxSS),coxSS + qcrit*sqrt(VcoxSS))
  expcoxci = exp(coxci)

  #inverse variance weights for reference
  weight.invVar = (1/coxvar)/sum(1/coxvar)

  #amalgamated results using sample size weights
  coxMRSSres = list(table=data.frame(loghrSS=coxSS,vSS=VcoxSS,ciSS=coxci,
                                     hrSS=HRcoxSS,expciSS=expcoxci,
                                     pvSS=pvcoxSS),weightSS=coxwSS,
                    weightIV = weight.invVar)

  coxstratci = cbind(coxbeta - qcrit*sqrt(coxvar),
                     coxbeta + qcrit*sqrt(coxvar))
  coxstratpv = 2*pnorm(abs(coxbeta/sqrt(coxvar)),0,1,lower.tail=FALSE)

  #by-stratum results
  coxMRSSmat = cbind(coxbeta,coxvar,coxstratci,exp(coxbeta),exp(coxstratci),
                     coxstratpv)
  colnames(coxMRSSmat) = c("bhat","v(bhat)","ci.lower","ci.upper","exp(bhat)",
                           "exp(ci.lower)","exp(ci.upper)","pval")
  coxMRSSmat.pt2 = matrix(cbind(unlist(coxwSS),unlist(weight.invVar)),
                          nrow=nstrata,byrow=FALSE)
  colnames(coxMRSSmat.pt2) = c("weight.SS","weight.invVar")
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
        brtimes = c(3,6,12)
        brtimeby = brtimes[which.min(abs(5.5-max(round(time))/brtimes))]
        breaktimes = seq(0,max(time),brtimeby)
      } else if (timeunit == "Years"){
        brtimes = c(0.25,0.5,1)
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
    sdata$strata<-rep(names(s$strata), s$strata)
    #fix to ensure it won't mess up the order for > 9 strata
    sdata$strata = as.numeric(gsub(".*=","",sdata$strata))
    #only one strata
    if (is.null(s$strata)) sdata$strata = rep(1,nrow(sdata))

    labelNames = paste0("S",1:nstrata)
    if (treetype == "preliminary") labelNames = paste0("p",labelNames)
    stratLabelName = ifelse(treetype=="preliminary","Preliminary","Final")
    if (treetype=="prespecified"){
      labelNames = paste0("S_",1:nstrata)
      stratLabelName = "Prespecified"
    }

    if (!is.null(termNodes)){
      labelNames = paste0(labelNames,": ",termNodes)
      labelNames = sapply(labelNames,function(x) wrapper(x,width=40))
      names(labelNames) = NULL
    }

    if (nstrata <=3){

      spectraledges = RColorBrewer::brewer.pal(5,"Spectral")[
        sort(c(1,5,4)[1:nstrata])]
      p4 = ggplot2::ggplot() +
        ggplot2::geom_step(data=sdata, ggplot2::aes(x=time, y=surv,
                                                    col=factor(sdata$strata)),
                           lwd=1.2) + ggplot2::theme_bw() +
        ggplot2::scale_x_continuous(breaks=breaktimes) +
        ggplot2::xlab(paste0("Time",timelabel)) + ggplot2::ylab("Survival") +
        ggplot2::scale_color_manual(name=paste0(stratLabelName," Risk Strata"),
                                    labels=labelNames,values=spectraledges) +
        ggplot2::geom_ribbon(data=sdata,ggplot2::aes(
          x=time,ymin=lower,ymax=upper,fill=factor(sdata$strata)),alpha=0.2) +
        ggplot2::scale_fill_manual(name=paste0(stratLabelName," Risk Strata"),
                                   labels=labelNames,values=spectraledges)

    }else {

      getPalette = grDevices::colorRampPalette(RColorBrewer::brewer.pal(min(
        nstrata,11),"Spectral"))

      p4 = ggplot2::ggplot() +
        ggplot2::geom_step(data=sdata, ggplot2::aes(x=time, y=surv,
                                                    col=factor(sdata$strata)),
                           lwd=1.2) + ggplot2::theme_bw() +
        ggplot2::scale_x_continuous(breaks=breaktimes) +
        ggplot2::xlab(paste0("Time",timelabel)) + ggplot2::ylab("Survival") +
        ggplot2::scale_color_manual(name=paste0(stratLabelName," Risk Strata"),
                                    labels=labelNames,values=getPalette(nstrata)) +
        ggplot2::geom_ribbon(data=sdata,ggplot2::aes(
          x=time,ymin=lower,ymax=upper,fill=factor(sdata$strata)),alpha=0.2) +
        ggplot2::scale_fill_manual(name=paste0(stratLabelName," Risk Strata"),
                                    labels=labelNames,values=getPalette(nstrata))
    }

  } else p3 = p4 = NULL

  #============================================================================#

  return(list(fitsummary=coxMRSSmat,stratafit=s,
              table=coxMRSSres$table,weights=coxMRSSres$weightSS,
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

roundLabelNums = function(name,digits){
  name = strsplit(name," ")[[1]]
  numpos = grep("\\d+[[:punct:]]\\d+",name)
  numparpos = grep("\\d+[[:punct:]]\\d+[[:punct:]]",name)
  puncttypeL = gsub("\\d+[[:punct:]]\\d+[[:punct:]]+$","",name[numparpos])
  puncttypeL[nchar(puncttypeL)>0] = paste0("\\",puncttypeL[nchar(puncttypeL)>0])
  puncttypeR = gsub("^[[:punct:]]*\\d+[[:punct:]]\\d+","",name[numparpos])
  puncttypeR[nchar(puncttypeR)>0] = paste0("\\",puncttypeR[nchar(puncttypeR)>0])
  if (length(numparpos)!=0){
    name[numparpos] = sapply(1:length(puncttypeL),function(x)
      sub(puncttypeL[x],"",name[numparpos[x]]))
    name[numparpos] = sapply(1:length(puncttypeR),function(x)
      sub(puncttypeR[x],"",name[numparpos[x]]))
  }
  name[numpos] = printR(as.numeric(name[numpos]),digits)
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
#' measure="OR, or columns "ratediff", "ci.lower", "ci.upper" for
#' family = "binomial", measure="RD", or columns "bhat", "ci.lower",
#' "ci.upper" for family = "gaussian"
#' @param MRSSres row matrix with amalgamated results summary, including estimate
#' and corresponding confidence intervals
#' @param stratafits survfit object with number of subjects, median survival
#' times, etc. within each strata
#' @param family Trait family, current options: "cox", "binomial, "gaussian"
#' @param measure Response of interest for binary traits; current options
#' are: "RD" (rate difference using Miettinen and Nurminen 1985 method), or
#' "OR" (odds ratio, using GLM model fit); ignored for family!="binomial"
#' @param treetype String, whether trees input are "preliminary" or "final"
#' @param labelNames Vector of names for each tree strata formed (e.g., number of
#' strata or covariate definition)
#' @param cilevel Confidence level alpha for overall result and confidence
#' intervals
#' @param ... Optional additional arguments passed into the forestplot function
#' (@seealso \code{\link[forestplot]{forestplot}})
#'
#' @return fplot An annotated forest plot showing by-stratum estimates and
#'  confidence intervals, as well as amalgamated result
strataforestplot = function(MRSSmat,MRSSres,stratafits,family,measure,
                            treetype="final",labelNames=NULL,cilevel,...){

  nstrata = nrow(MRSSmat)
  if (is.null(labelNames)){
    labelNames = paste0("S",1:nstrata)
    if (treetype == "preliminary") labelNames = paste0("p",labelNames)
  }
  labelNames = sapply(labelNames,function(x) roundLabelNums(x,1))
  labelNames = sapply(labelNames,function(x) gsub("\\|","\\|\n",x))
  labelNames = sapply(labelNames,function(x)
    paste(strwrap(x, 70), collapse = "\n"))
  nlines = max(sapply(labelNames,function(x) stringr::str_count(x, "\n")))

  tabletext = paste0(stratafits$n," (",printR(MRSSmat$weight.SS*100,1),")")

  tabletext.all = cbind(c(paste0(toupper(substr(treetype,1,1)),
                                 substr(treetype,2,100),
                                 " Strata"),labelNames,"5-STAR Average"),
                        c("No. Subjects (%)",tabletext,
                          paste0(sum(stratafits$n)," (100)")))

  #============================================================================#

  if (family == "cox"|(family=="binomial" & measure=="OR")){

    nevents = summary(stratafits)$table[,"events"]
    tabletext.events = paste0(nevents," (",
                              printR((nevents/sum(nevents))*100,1),")")
    tabletext.weights = MRSSmat$weight.adap

    measurename = measure
    if (family == "cox") measurename = "HR"

    tabletext.ci = c(paste0(printR(MRSSmat$exp.bhat.,2)," (",
                            printR(MRSSmat$exp.ci.lower.,2),
                            ", ",printR(MRSSmat$exp.ci.upper.,2),
                            ")"),
                     paste0(printR(MRSSres["exp(beta)"],2),
                            " (",printR(MRSSres["exp(ci lower)"],2),
                            ", ",printR(MRSSres["exp(ci upper)"],2),")") )

    weightName = ifelse(treetype=="prespecified","Inv.Var Weight %",
                        "Adap. Weight %")

    tabletext.all = cbind(c(paste0(toupper(substr(treetype,1,1)),
                                   substr(treetype,2,100),
                                   " Strata"),labelNames,"5-STAR Average"),
                          c("No. Subjects (%)",tabletext,
                            paste0(sum(stratafits$n)," (100)")),
                          c("No. Events (%)",tabletext.events,
                            paste0(sum(stratafits$n.event)," (100)")),
                          c(paste0(measurename," (",(1-cilevel)*100,"% CI)"),
                            tabletext.ci),
                          c(weightName,printR(tabletext.weights*100,1),100))

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
    if (is.infinite(max(upper))) clippts[2] = max(upper[!is.infinite(upper)])+0.1

    xlabname=ifelse(family=="cox","Hazard Ratio","Odds Ratio")
    forestplot::forestplot(labeltext=tabletext.all,align=c("l","c","c","c","c"),
                           graph.pos=4,graphwidth=grid::unit(4,"inches"),
                           mean=c(NA,MRSSmat[,"exp.bhat."],MRSSres["exp(beta)"]),
                           lower=c(NA,lower),upper=c(NA,upper),
                           is.summary=c(TRUE,rep(FALSE,nstrata),TRUE),xlab=xlabname,zero=1,
                           col=forestplot::fpColors(lines="darkblue",box="darkblue",
                                                    summary=c("blue")),
                           lwd.zero=grid::gpar(lwd=2),xticks.digits = 2,
                           fn.ci_sum = colfn,clip=clippts, colgap=grid::unit(4,"mm"),
                           lineheight=grid::unit(1.4+0.3*max(0,nlines-3),"cm"),boxsize=bxsize,
                           #line.margin=0.5*nlines,
                           txt_gp=forestplot::fpTxtGp(cex=0.9,xlab=grid::gpar(cex=0.9),
                                                      ticks=grid::gpar(cex=0.9),
                                                      label=grid::gpar(lineheight=0.75)),...)
    # txt_gp=forestplot::fpTxtGp(xlab=grid::gpar(cex=0.8),
    #                            ticks=grid::gpar(cex=0.8),
    #                            label=grid::gpar(cex=0.8),
    #                            summary=grid::gpar(cex=0.8),
    #                            title=grid::gpar(cex=0.8)),...)

    fplot = recordPlot()

    #============================================================================#

  } else if (family=="binomial" & measure=="RD"){

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

    forestplot::forestplot(
      labeltext=c(labelNames,"5-STAR Average"),align="c",
      graphwidth=grid::unit(4,"inches"),
      mean=c(MRSSmat[,"ratediff"],MRSSres["ratediff"]),
      lower=c(MRSSmat[,"ci.lower"],MRSSres["ci lower"]),
      upper=c(MRSSmat[,"ci.upper"],MRSSres["ci upper"]),
      is.summary=c(rep(FALSE,nstrata),TRUE),xlab="Rate Difference",zero=0,
      col=forestplot::fpColors(lines="forestgreen",
                               box="forestgreen",
                               summary="green3"),fn.ci_sum = colfn2,
      lineheight=grid::unit(1.4,"cm"),lwd.zero=2,boxsize=bxsize,
      txt_gp=forestplot::fpTxtGp(cex=0.9,xlab=grid::gpar(cex=0.9),
                                 ticks=grid::gpar(cex=0.9)),...)
    # txt_gp=forestplot::fpTxtGp(xlab=grid::gpar(cex=0.8),
    #                            ticks=grid::gpar(cex=0.8),
    #                            label=grid::gpar(cex=0.8),
    #                            summary=grid::gpar(cex=0.8),
    #                            title=grid::gpar(cex=0.8)),...)

    fplot = recordPlot()

  } else { #family = "gaussian"

    forestplot::forestplot(
      labeltext=c(labelNames,"5-STAR Average"),align="c",
      graphwidth=grid::unit(4,"inches"),
      mean=c(MRSSmat[,"bhat"],MRSSres["beta"]),
      lower=c(MRSSmat[,"ci.lower"],MRSSres["ci lower"]),
      upper=c(MRSSmat[,"ci.upper"],MRSSres["ci upper"]),
      is.summary=c(rep(FALSE,nstrata),TRUE),xlab="Mean Difference",zero=0,
      col=forestplot::fpColors(lines="firebrick3",
                               box="firebrick3",
                               summary=c("firebrick1")),
      lineheight=grid::unit(1.4,"cm"),
      txt_gp=forestplot::fpTxtGp(cex=0.9,xlab=grid::gpar(cex=0.9),
                                 ticks=grid::gpar(cex=0.9)),...)
    # txt_gp=forestplot::fpTxtGp(xlab=grid::gpar(cex=0.9),
    #                            ticks=grid::gpar(cex=0.9),
    #                            label=grid::gpar(cex=0.9),
    #                            summary=grid::gpar(cex=0.9),
    #                            title=grid::gpar(cex=0.9)),...)

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
#' providing estimate of regression coefficient or log odds ratio
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
#'      estimate, variance, p-value, and cilevelx100% confidence interval
#'      assuming sample size weights
#'     \item weights: Sample size weights used to construct estimate
#' }
glmbystrata = function(yy,arm,family,treeStrata,termNodes=NULL,cilevel=0.05,
                       verbose=0){

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
  qcrit = qnorm(1-cilevel/2,0,1)

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
## Estimate Rate Difference Within Strata ##
##========================================##

#' Estimate rate difference Within Strata
#'
#' Estimates rate difference and significance within each formed strata using
#' methodology of Miettinen and Nurminen (1985)
#'
#' @param yy      Case control status response vector
#' @param arm     Treatment indicator, 1 = test treatment, 0 = control
#' @param family  Family for glm; options: "binomial" or "gaussian"
#' @param treeStrata Vector of strata membership, as defined by ctree strata
#' @param treetype String, whether trees input are "preliminary" or "final"
#' @param termNodes Vector of names for each tree strata formed (e.g., number of
#' strata or covariate definition)
#' @param cilevel Confidence level alpha for overall result and confidence
#' intervals
#' @param verbose Numeric variable indicating amount of information to print
#' to the terminal (0 = nothing, 1 = notes only, 2 = notes and intermediate output)
#' @param plot    Logical, whether to create rate difference forest plots
#'
#' @return \itemize{
#'     \item fitsummary: summary of rate difference within each strata
#'     \item table: Summary of amalgamated rate difference estimate, variance,
#'     p-value, and cilevelx100% confidence interval assuming sample size weights
#'     \item weights: Sample size weights used to construct estimate
#' }
ratediffbystrata = function(yy,arm,family,treeStrata,treetype="final",
                            termNodes=NULL,cilevel=0.05,
                            verbose=0,plot=TRUE){

  dat = data.frame(yy,arm)
  nstrata = max(treeStrata,na.rm=TRUE)
  if (is.null(termNodes)) termNodes = unique(treeStrata)

  #alt to glm: rate difference test within each strata
  rdCtreeStrataFit = lapply(1:nstrata,function(x) {
    datx = dat[treeStrata==x,]
    c1 = sum(datx$yy[datx$arm==1])
    S1 = sum(datx$arm==1)
    c0 = sum(datx$yy[datx$arm==0])
    S0 = sum(datx$arm==0)

    r1 = c1/S1; r0 = c0/S0
    c=c1+c0; S=S1+S0; r = c/S

    #proper chisq stat for RD = 0 (R1 = R0):
    RDall=ratesci::scoreci(c1,S1,c0,S0,distrib="bin",measure="RD",level=1-cilevel,
                           skew=FALSE,weighting="MN")
    Ts0 = ( (r1-r0)^2 )/( r*(1-r)*(S/(S-1))*(1/S1+1/S0) )

    #########calc variance
    RD = 0
    L0 = c0*RD*(1-RD)
    L1 = (S0*RD - S - 2*c0)*RD+c
    L2 = (S1+2*S0)*RD-S-c
    L3 = S

    q = L2^3/(3*L3)^3 - L1*L2/(6*L3^2) + L0/(2*L3)
    sgn = sign(q)
    p = sgn*(L2^2/(3*L3)^2 - L1/(3*L3))^(1/2)
    if (sign(p) != sign(q)) p = -p
    a = (1/3)*(pi + acos(q/p^3))
    R0tilde = 2*p*cos(a) - L2/(3*L3)
    R1tilde = R0tilde + RD

    VarRD = (R1tilde*(1-R1tilde)/S1 + R0tilde*(1-R0tilde)/S0)*(S/(S+1))
    ################################
    RDall$var = VarRD
    RDall

  })

  #pulling out relevant summary statistics from glm fits to calculate glmSS
  rdMat = matrix(unlist(lapply(rdCtreeStrataFit,function(x){
    cbind(x$estimates[,"MLE"],x$var)
  } )),ncol=2,byrow=TRUE)
  cnj = sapply(1:nstrata,function(x) sum(treeStrata==x))
  rdMat = cbind(rdMat,cnj)
  colnames(rdMat) = c("rdj","vhatj","nj")

  #calculating rdS
  rdest = rdMat[,"rdj"]
  rdvar = rdMat[,"vhatj"]

  rdwSS = cnj/sum(cnj)
  rdSS = sum(rdest*rdwSS)
  rdSSp = rdSS*100
  VrdSS = sum(rdwSS^2*rdvar)
  TrdSS = rdSS/sqrt(VrdSS)
  pvrdSS = 2*pnorm(abs(TrdSS),lower.tail=FALSE,log.p=FALSE)

  #two-sided test
  qcrit = qnorm(1-cilevel/2,0,1)

  #overall and within-strata cis
  rdci = cbind(rdSS - qcrit*sqrt(VrdSS),rdSS + qcrit*sqrt(VrdSS))
  rdcip = rdci*100
  rdstratci2 = cbind(rdest - qcrit*sqrt(rdvar),
                     rdest + qcrit*sqrt(rdvar))
  rdstratci = matrix(unlist(lapply(rdCtreeStrataFit,function(x){
    x$estimates[,c("Lower","Upper")]
  })),ncol=2,byrow=TRUE)

  rdstratpv = sapply(rdCtreeStrataFit,function(x){
    x$pval[,"pval2sided"]
  })
  rdstratpv2 = 2*pnorm(abs(rdest/sqrt(rdvar)),0,1,lower.tail=FALSE)

  #compiling within-strata effect estimate summary matrix
  rdMRSSmat = cbind(rdest,rdvar,rdstratci,rdstratpv)
  colnames(rdMRSSmat) = c("ratediff","v(ratediff)","ci.lower","ci.upper","pval")

  rdMRSSres = list(table=data.frame(ratediffSS=rdSS,vSS=VrdSS,ciSS=rdci,
                                    pvSS=pvrdSS),weightSS=rdwSS)

  rdMRSSmat.pt2 = matrix(unlist(rdwSS),nrow=nstrata,byrow=FALSE)
  colnames(rdMRSSmat.pt2) = c("weight.SS")
  rdMRSSmat = cbind(rdMRSSmat,rdMRSSmat.pt2)
  labelNames = paste0("S",1:nstrata)
  if (treetype == "preliminary"){
    labelNames = paste0("p",labelNames)
  }
  rownames(rdMRSSmat) = labelNames

  if (verbose > 1){
    print("Estimated rate difference within each strata")
    print(round(rdMRSSmat,4))
  }

  #============================================================================#

  return(list(fitsummary=rdMRSSmat,table=rdMRSSres$table,
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
                                          rep(c("Control","Test"),each=nstrata),sep="-")
  if (family == "cox"){
    rownames(perStrataEventSummary) = c("nevents","nsubj","pctcensor")
  } else if (family == "binomial"){
    rownames(perStrataEventSummary) = c("ncases","nsubj","pctcontrols")
  }

  return(perStrataEventSummary)

}


################################################################################
