#packrat:::recursivePackageDependencies("fiveSTAR",lib.loc = .libPaths()[1])

#from Hadley Wickham, http://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-td859876.html
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

################################################################################

##================##
## Filtering Step ##
##================##


#' Filter covariates for 5STAR
#'
#' Performs Step 2 of the 5-STAR Algorithm:
#' Forms dummy variable matrix for factors, and fits an elastic net (ENET) or
#' random forest (RF) model to determine which covariates to keep for building
#' trees in 5-STAR
#'
#' @param yy      Response - either Surv() object for time to event data or 1
#' column matrix of case/control status for binary data
#' @param X       Data frame of all possible stratification covariates
#' @param family  Trait family, current options: "cox", "binomial", or "gaussian"
#' @param plot    Whether to make plots for filter results (e.g., variable importance
#' plot for RF-based filtering or solution path for ENET-based filtering)
#' @param verbose Numeric variable indicating amount of information to print
#' to the terminal (0 = nothing, 1 = notes only, 2+ = notes and intermediate output)
#' @param filter.hyper List of control parameters for filtering step
#' (see also \code{\link{filter_control}}); key agruments include: \itemize{
#' \item method: filtering method; current options are: "ENET" (Elastic Net),
#' "RF" (Random Forest with default or input parameters, using delete-d
#' jackknife confidence intervals for VIMP for variable selection), "RFbest"
#' (Random Survival Forest, performing tuning for mtry and nodesize, and using
#' double bootstrap confidence intervals for VIMP for variable selection; more
#' accurate but slower than RF option). See \code{\link[glmnet]{glmnet}} for more
#' details on the elastic net and \code{\link[randomForestSRC]{rfsrc}} for more
#' details on the random forest filtering
#'  \item lambdatype: Optional elastic net parameter; whether to use the tuning
#'  parameter lambda that minimizes cross validation error
#'  (lambdatype="min", default) or the largest lambda that gives error within 1
#'  standard error of the minimum error (lambdatype="1se"). Ignored when method
#'  is "RF" or "RFbest"
#' \item mixparm: Optional elastic net mixing parameter alpha or grid of alpha
#' values to search over. If nothing is entered, the algorithm will perform a
#' grid search for the best mixing parameter between 0.05 and 0.95.
#' Ignored when method is "RF" or "RFbest"
#'  \item vimpalpha: Optional significance level for RF VIMP confidence intervals
#'  \item nfolds: number of folds to use for cross validation tuning of lambda
#'  parameter for elastic net filtering. Default = 10. Ignored when method="RF"
#'  or "RFbest"
#'  \item filterseed: optional seed for filtering step
#'  \item ... : Optional arguments for glmnet or rfsrc
#' }
#' @param vars2keep: List of variable names (matching column names in X) of
#' variables to be passed through filtering step without penalization, etc.
#' This is ignored in the main 5-STAR algorithm but may be used when filtering
#' step is used as a stand alone function. Currently only used when method = ENET
#'
#' @return cov2keep: List of all covariates kept after filtering step,
#' selected by elastic net or random forest
#' @return For method=ENET, additionally outputs \itemize{
#'    \item cvout: vector containing fraction of null deviance explained,
#'    mean cross validation error, and optimal tuning parameter values
#'     (see \code{\link[glmnet]{cv.glmnet}} and \code{\link[glmnet]{glmnet}}
#'      for more details)
#'    \item beta: the coefficients of glmnet fit given tuned parameters
#'    \item ENETplot: plot containing the deviance over different combinations of
#'    the mixing and tuning parameters alpha and lambda, with optimal
#'    alpha shown in blue and minimum deviance point for each alpha in black
#'    (left panel), and solution path for the best mixing parameter (right panel).
#'    Output when plot = TRUE
#' }
#' @return For method=RF or RFbest, additionally outputs:
#' \itemize{
#'    \item varselectmat:  matrix of VIMP confidence intervals, p-values,
#'    and selection decision for each variable
#'    \item VIMPplot: default variable importance CI plot, output if plot = TRUE
#'    \item varselect2plot: variable selection information from subsample.rfsrc for
#' making customizable VIMP CI plots
#'    \item forest: rfsrc object
#' }
#'
#' @export
filter5STAR = function(yy,X,family="cox",plot=FALSE,verbose=0,
                       filter.hyper=filter_control(
                         method="ENET",lambdatype="min",mixparm=NULL,
                         vimpalpha=0.05,nfolds=10,filterseed=2019,...),
                       vars2keep=NULL){

  method = filter.hyper$method

  X = droplevels(X)

  #setting seed
  set.seed(filter.hyper$filterseed)

  ##=================##
  ## METHOD = "ENET" ##
  ##=================##
  if (method == "ENET"){

    #----------------------------------#
    # pulling out important parameters #
    #----------------------------------#
    nfolds = filter.hyper$nfolds
    lambdatype = filter.hyper$lambdatype
    mixparm = filter.hyper$mixparm

    optargs = filter.hyper[!(names(filter.hyper) %in%
                               c("method","lambdatype","mixparm","vimpalpha",
                                 "filterseed","nfolds"))]

    if (!(lambdatype %in% c("1se","min"))) stop("lambdatype must be one of
                                                '1se',`min'.")

    #-----------------------------------------------------------------#
    # data cleaning/putting data in correct format for enet filtering #
    #-----------------------------------------------------------------#

    #ensuring yy is a column
    if (is.null(nrow(yy))) yy = as.matrix(yy)

    #for now, complete case analysis for filtering
    #(i.e., removing individuals with any missing covariate information)
    XC.nomiss  = X[rowSums(is.na(X))==0,]
    yyC.nomiss = yy[rowSums(is.na(X))==0,]

    # removing also any 0's from the dataset for time to event data
    if (family == 'cox'){
      XC.nomiss <- XC.nomiss[yyC.nomiss[,1] > 0,]
      yyC.nomiss <- yyC.nomiss[yyC.nomiss[,1] > 0,]
    }

    #continuous (and possibly binary; non-factor) variables
    XC.numeric = XC.nomiss[,which(lapply(XC.nomiss,is.numeric)==TRUE)]
    if (is.null(nrow(XC.numeric))){
      XC.numeric = data.frame(XC.numeric)
      names(XC.numeric) = names(which(lapply(XC.nomiss,is.numeric)==TRUE))
    }

    #creating dummy variables for all categorical variables
    dummyNames = names(which(lapply(XC.nomiss,is.factor)==TRUE))
    #checking if any levels need to be dropped after data cleaning
    #(if any levels have no obs)
    XCd.dummy = XC.nomiss[,dummyNames]
    XCd.dummy = do.call(data.frame,lapply(dummyNames,function(x)
      droplevels(XC.nomiss[,x])))

    #removing covariates that have only 1 level after removing missing data
    polymorphicDummy = sapply(XCd.dummy,function(x) length(levels(x)))>1
    monoDummyNames = dummyNames[!polymorphicDummy]
    if (length(monoDummyNames >0) & verbose > 0) print(paste0(
      "Note: removing covariates ",paste0(monoDummyNames,collapse=", "),
      " from analysis (only 1 level after removing missing subjects)"))
    dummyNames = dummyNames[!(dummyNames %in% monoDummyNames)]
    if (length(dummyNames)>0) dummyNames = paste0(dummyNames,".d")
    XCd.dummy = XCd.dummy[,sapply(XCd.dummy,function(x) length(levels(x)))>1]

    if (is.null(nrow(XCd.dummy))) XCd.dummy = data.frame(XCd.dummy)
    if (ncol(XCd.dummy)>0) colnames(XCd.dummy) = dummyNames

    if (ncol(X)!= ncol(XC.numeric)+ncol(XCd.dummy)+length(monoDummyNames))
      warning("No. columns in numeric and factor matrices not equal to number of covariates")

    if (length(dummyNames)>0){
      dummyFormula = paste(dummyNames,collapse=" + ")
      XC.dummy = model.matrix(as.formula(paste0("~ ",dummyFormula)), XCd.dummy)
      XC.enet = data.matrix(cbind(XC.dummy,XC.numeric)[,-1])
    } else XC.enet = data.matrix(XC.numeric)

    #penalty factor - allowing it to be 0 for those vars to keep
    penaltyfac = rep(1,ncol(XC.enet))
    if (!is.null(vars2keep)&length(vars2keep > 0)){
      vars2keepList = paste0(vars2keep,collapse="|")
      penaltyfac[grep(vars2keepList,colnames(XC.enet))] = 0
    }

    #calculating cross validation folds to keep constant across
    #grid search of alpha
    if (family == "cox"){
      foldid=quiet(c060::balancedFolds(class.column.factor=yyC.nomiss[,2],
                                       cross.outer=nfolds))
    } else if (family=="binomial") {
      foldid=quiet(c060::balancedFolds(class.column.factor=yyC.nomiss,
                                       cross.outer=nfolds))
    } else if (family == "gaussian"){
      foldid = sample(1:nfolds,size=length(yyC.nomiss),replace=TRUE)
    }

    #setting mixing parameter to default if not input into function
    if (is.null(mixparm)){
      a <- seq(0.05, 0.95, 0.05)
    } else a = mixparm

    #-------------------------------------------------------#
    # fitting elastic net, doing grid search for best alpha #
    #-------------------------------------------------------#

    plotdata = list()
    search = NULL
    kk = 0
    for(i in a) {
      kk = kk + 1
      cvfit = do.call(glmnet::cv.glmnet,c(list(
        XC.enet, yyC.nomiss,family = family,alpha = i, foldid=foldid,
        penalty.factor=penaltyfac),optargs))

      plotdata[[kk]] = data.frame(lambda=cvfit$lambda,deviance=cvfit$cvm,
                                  lower=cvfit$cvlo,upper=cvfit$cvup,alpha=i)

      ##here ... is: type.measure = "deviance", standardize=TRUE
      # cvfit = glmnet::cv.glmnet(XC.enet, yyC.nomiss, family = family,
      #                           alpha = i, foldid=foldid, ...)
      optlambda = ifelse(lambdatype=="1se",cvfit$lambda.1se,cvfit$lambda.min)
      search = rbind(search,data.frame(cvm = cvfit$cvm[
      cvfit$lambda == optlambda],lambda = optlambda, alpha = i))
    }
    cv <- search[search$cvm == min(search$cvm), ]
    md <- do.call(glmnet::glmnet,c(list(
      XC.enet, yyC.nomiss, family = family,lambda = cv$lambda,
      alpha = cv$alpha, penalty.factor=penaltyfac),optargs))

    cvoutput = cbind(md$dev.ratio,cv)
    colnames(cvoutput) = c("dev.ratio","cvm","lambda","alpha")

    #-----------------------------#
    # ENET plots (when plot=TRUE) #
    #-----------------------------#

    plotdata = do.call(rbind,plotdata)
    best = ifelse(plotdata[,"alpha"]==cv$alpha,TRUE,FALSE)
    if (plot){

      #deviance over tuning paramters plot
      devalphaPlot = ggplot2::ggplot(plotdata) +
        ggplot2::geom_point(ggplot2::aes(x=log(lambda),y=deviance,group=alpha,
                                col=best),size=1) +
        ggplot2::geom_line(ggplot2::aes(x=log(lambda),y=deviance,group=alpha,
                                        color=best)) +
        ggplot2::geom_point(ggplot2::aes(x=log(lambda),y=cvm,group=alpha),
                            data=search,pch=19,color="black",size=1.5) +
        ggplot2::theme_classic() + ggplot2::theme(legend.position="none") +
        ggplot2::scale_color_manual(values=c("gray","blue")) +
        ggplot2::ylab("Deviance") +
        ggplot2::theme(axis.text=ggplot2::element_text(size=ggplot2::rel(1.4)),
                       axis.title=ggplot2::element_text(size=ggplot2::rel(1.4))) +
        ggplot2::xlab(expression(paste("log(",lambda,")")))

      #solution path plot
      cvfit.plot = do.call(glmnet::glmnet,c(list(
        XC.enet, yyC.nomiss,family = family,alpha = cv$alpha,
        penalty.factor=penaltyfac),optargs))

      solnpathdata = data.frame(beta=as.matrix(cvfit.plot$beta),
                                cov=rownames(cvfit.plot$beta))
      colnames(solnpathdata) = c(cvfit.plot$lambda,"cov")
      solnpathdata = tidyr::gather(solnpathdata,key=lambda,value=beta,-cov)
      solnpathdata$covgroup = sapply(solnpathdata$cov,function(x)
        trimws(strsplit(as.character(x),"\\.")[[1]][1]))
      solnpathdata$lambda = as.numeric(solnpathdata$lambda)

      #find coefficients corresponding to optimal lambda by cv tuning
      #note: currently some issue with mismatch between plot and cov2keep
      #when best beta value is small but non-zero (e.g., 6e-16)
      best.data = solnpathdata[(abs(solnpathdata$lambda - cv$lambda) < 1e-6),]
      tocolor = unique(best.data$covgroup[abs(best.data$beta) > 0])
      covgroup2 = solnpathdata$covgroup
      covgroup2[!(covgroup2 %in% tocolor)] = "notselected"
      covgroupTF = covgroup2 == "notselected"

      colvals = rep("darkgray",length(unique(solnpathdata$covgroup)))
      #if all are selected, need the full number of covariate groups; otherwise
      #need one less (i.e., -1 for the gray "not selected" category)
      nuniquecolvals = ifelse(mean(covgroupTF)==0,length(unique(covgroup2)),
                                   length(unique(covgroup2))-1)
      colvals[gtools::mixedsort(unique(solnpathdata$covgroup)) %in%
                 intersect(solnpathdata$covgroup,covgroup2)] = rainbow(nuniquecolvals)

      solnpathPlot = ggplot2::ggplot(solnpathdata) +
        ggplot2::geom_line(ggplot2::aes(
          x=log(lambda),y=beta,group=cov,
          col=factor(covgroup,levels=gtools::mixedsort(unique(covgroup)))),lwd=0.8) +
        ggplot2::theme_classic() + ggplot2::geom_hline(yintercept=0) +
        ggplot2::geom_vline(xintercept=log(cv$lambda),col="black",lty=2) +
        ggplot2::scale_color_manual(name="Covariate",labels=paste0(
          "X",seq(1:length(unique(solnpathdata$covgroup)))),
          values = colvals) +
        ggplot2::theme(axis.text=ggplot2::element_text(size=ggplot2::rel(1.4)),
                       axis.title=ggplot2::element_text(size=ggplot2::rel(1.4)),
                       legend.text = ggplot2::element_text(size=ggplot2::rel(1.1)),
                       legend.title = ggplot2::element_text(size=ggplot2::rel(1.2))) +
        ggplot2::ylab(expression(beta)) +
        ggplot2::xlab(expression(paste("log(",lambda,")")))

      enetplot = gridExtra::arrangeGrob(grobs=list(devalphaPlot,solnpathPlot),
                                        ncol=2,widths=c(0.65,1))
    } else enetplot = NULL

    enetcoef = coef(md)

    #extracting covariates to keep from filtering step
    cov2keep = names(enetcoef[,1])[abs(enetcoef[,1])>0]
    if (length(dummyNames)>0){
      cov2keep.cat = dummyNames[sapply(dummyNames,function(x)
        length(grep(x,cov2keep))>0)]
      cov2keep.cat = sub(".d","",cov2keep.cat)
    } else cov2keep.cat = NULL
    cov2keep.num = cov2keep[cov2keep %in% names(XC.numeric)]
    cov2keep.all = c(cov2keep.cat,cov2keep.num)

  ##===============##
  ## METHOD = "RF" ##
  ##===============##

  } else if (method == "RF"){

    vimpalpha = filter.hyper$vimpalpha
    optargs = filter.hyper[!(names(filter.hyper) %in%
                               c("method","lambdatype","mixparm","vimpalpha",
                                 "filterseed","nfolds"))]

    if (family=="cox"){
      dat = data.frame(yy[,1],yy[,2],X)
      colnames(dat) = c("time","status",colnames(X))
      srv.o = do.call(randomForestSRC::rfsrc,c(list(
        formula=Surv(time,status)~.,data=dat),optargs))
    } else if (family=="binomial"|family=="gaussian") {
      if (family=="binomial") yy = as.factor(yy)
      dat = data.frame(yy,X)
      srv.o = do.call(randomForestSRC::rfsrc,c(list(
        formula=yy~.,data=dat),optargs))
    }

    srv.smp.o = randomForestSRC::subsample(srv.o)
    if (plot){
      plot(srv.smp.o)
      pvimp = recordPlot()
    } else pvimp = NULL
    srv.varselect.o=randomForestSRC::extract.subsample(
      srv.smp.o,alpha=vimpalpha)$var.jk.sel.Z
    cov2keep.all = rownames(srv.varselect.o)[which(srv.varselect.o$signif)]

  ##===================##
  ## METHOD = "RFbest" ##
  ##===================##
  } else if (method == "RFbest"){

    vimpalpha = filter.hyper$vimpalpha
    #note: no optional arguments here?

    if (family=="cox"){
      dat = data.frame(yy[,1],yy[,2],X)
      colnames(dat) = c("time","status",colnames(X))
      srv.tune.o = randomForestSRC::tune.rfsrc(Surv(time,status)~.,dat)
      srv.o = randomForestSRC::rfsrc(Surv(time,status)~.,dat,
                                     nodesize=srv.tune.o$optimal[1],
                                     mtry=srv.tune.o$optimal[2])
    } else if (family=="binomial"|family=="gaussian") {
      if (family=="binomial") yy = as.factor(yy)
      dat = data.frame(yy,X)
      srv.tune.o = randomForestSRC::tune.rfsrc(yy~.,dat)
      srv.o = randomForestSRC::rfsrc(yy~.,dat,
                                     nodesize=srv.tune.o$optimal[1],
                                     mtry=srv.tune.o$optimal[2])
    }

    srv.smp.o = randomForestSRC::subsample(srv.o,B=100,bootstrap=TRUE)

    if (plot){
      plot(srv.smp.o)
      pvimp = recordPlot()
    } else pvimp = NULL
    srv.varselect.o=randomForestSRC::extract.subsample(
      srv.smp.o,alpha=vimpalpha)$var.sel.Z
    cov2keep.all = rownames(srv.varselect.o)[which(srv.varselect.o$signif)]

  } else{
    warning("Invalid filter type. Please select one of 'ENET','RF', or
                 'RFbest'. Returning all covariates.")
    cov2keep.all = colnames(X)
  }

  ##==============================================##
  ## collecting final function results for output ##
  ##==============================================##
  covOrder=match(colnames(X),cov2keep.all)[
    !is.na(match(colnames(X),cov2keep.all))]
  cov2keep.all = cov2keep.all[covOrder]

  if (method=="ENET"){
    return(list(cov2keep=cov2keep.all,cvout=cvoutput,beta=enetcoef,
                ENETplot=enetplot))
  } else if (method == "RF"|method=="RFbest"){
    return(list(cov2keep=cov2keep.all,varselectmat=srv.varselect.o,
                  VIMPplot=pvimp,varselect2plot=srv.smp.o,forest=srv.o))
  }
}

################################################################################

##===========================================================================##
## simplifyand: Used to simplify the output from 5-STAR tree step to easier  ##
## to read/interpret. Currently combines numeric results for same variable   ##
##===========================================================================##
simplifyand = function(node){

  #sepand = do.call(rbind,strsplit(node,"&"))
  #regular expression from Wiktor Stribizew (https://stackoverflow.com/questions/39733645/split-string-by-space-except-whats-inside-parentheses)
  nodeelements = trimws(strsplit(node,"&")[[1]])
  septab = do.call(rbind,sapply(nodeelements,function(x)
    strsplit(x, "(\\((?:[^()]++|(?1))*\\))(*SKIP)(*F)| ", perl=TRUE)))
    #strsplit(x," |(?>\\(.*?\\).*?\\K(, |$))",perl=TRUE)))

  possiblecombos = names(table(septab[,1]) > 1)
  cnd = lapply(possiblecombos,function(x){

    condition = NULL
    septabcombo = septab[septab[,1]==x,]
    if (is.null(dim(septabcombo)) ) septabcombo = matrix(septabcombo,nrow=1)

    eq = (septabcombo[,2] == "=")
    if (sum(septabcombo[,2]=="%in%")>0){

      valin = septabcombo[septabcombo[,2]=="%in%",3]
      valin = sapply(valin,function(x) gsub("c\\(","",x))
      valin = sapply(valin,function(x) gsub("\\)","",x))
      valin = sapply(valin,function(x) strsplit(x,", "))
      condition = c(condition,paste0(x," %in% c(",paste0(
        Reduce(intersect,valin),collapse=","),")"))

    } else if (sum(eq) > 0){

      valeq = septabcombo[eq,3]
      if ( length(setdiff(unique(eval(parse(text=x))),unique(valeq))) > 0 ){
        condition = c(condition,paste(x,"=",valeq))
      }

    } else {

      lt = (septabcombo[,2] == "<=" | septabcombo[,2] == "<")
      if (sum(lt)>0){
        vallt = min(septabcombo[lt,3])
        signchar = ifelse(any(septabcombo[,2]=="<"),"<","<=")
        condition = c(condition,paste(x,signchar,vallt))
      }

      gt = ( septabcombo[,2] == ">="|  septabcombo[,2] == ">")
      if (sum(gt)>0){
        valgt = max(septabcombo[gt,3])
        signchar2 = ifelse(any(septabcombo[,2]==">"),">",">=")
        condition = c(condition,paste(x,signchar2,valgt))
      }

    }

    condition = paste(condition,collapse=" & ")

  })

  fullcnd = paste(cnd,collapse=" & ")
  return(fullcnd)

}

##=============================##
## IGNORING SIMPLIFYOR FOR NOW ##
##=============================##
# #for now just checking simplest case: if node is of the form
# #(stuff & Xcondition) | (stuff & Xconditioncomplement)
# simplifyor = function(node){
#
#   septab = strsplit(node,"\\|")
#
#   #separating ands within each side of the or's
#   #(currently assumes only one "&" per subnode?)
#   sepand = lapply(septab,function(x) strsplit(x,"\\&"))
#   if (length(sepand[[1]])==2){
#
#     Rside = sub("^\\(","",trimws(unlist(lapply(sepand[[1]],function(x) x[1]))))
#     Rside.unique = unique(Rside)
#     Lside = sub("\\)$","",trimws(unlist(lapply(sepand[[1]],function(x) x[2]))))
#     Lside.unique = unique(Lside)
#
#     if (length(Rside.unique)==1) {
#
#       septab2 = do.call(rbind,sapply(Lside.unique,function(x) strsplit(x," ")))
#       if (length(unique(septab2[,1])==1) &
#           (septab2[1,2]=="<="&septab2[2,2]==">")|(septab2[1,2]==">="&
#                                                   septab2[2,2]=="<") &
#           length(unique(septab2[,3])==1)) fullcnd = Rside.unique
#
#     } else if (length(Lside.unique)==1){
#
#       septab2 = do.call(rbind,sapply(Rside.unique,function(x) strsplit(x," ")))
#       if (length(unique(septab2[,1])==1) &
#           (septab2[1,2]=="<="&septab2[2,2]==">")|(septab2[1,2]==">="&
#                                                   septab2[2,2]=="<") &
#           length(unique(septab2[,3])==1)) fullcnd = Lside.unique
#     } else fullcnd = node
#
#   } else {
#     fullcnd = node
#   }
#   return(fullcnd)
# }

################################################################################

##======================##
## Fit Trees for 5-STAR ##
##======================##

#' Fit Trees for 5-STAR
#'
#' Fits preliminary (3A) and pooled final (3B) trees to determine homogenous risk
#' strata for use with 5-STAR algorithm
#'
#' @param yy      Trait/response (either a binary or continuous
#'  covariate or a Surv() object summarizing follow-up time for right-censored
#'  data and status indicator where 1=dead, 0=censored)
#' @param X       Data frame of all possible stratification covariates
#' @param family  Trait family, current options: "cox", "binomial", or "gaussian"
#' @param verbose Numeric variable indicating amount of information to print
#' to the terminal (0 = nothing, 1 = notes only, 2+ = notes and intermediate output)
#' @param tree.hyper List of control variables for tree fitting (see
#' \code{\link{tree_control}} for details), many of which will be passed into
#'  \code{\link[partykit]{ctree_control}}.
#'
#' @return \itemize{
#'     \item strataids: vector of subject-level strata membership (1 = highest risk,
#'     2 = 2nd highest risk, etc.)
#'     \item stratadefn: definition of final formed strata, in terms of covariates
#'     \item finaltree: ctree object; final, pruned tree built by ctree,
#'     in terms of preliminary strata membership
#'     \item prelimstrataids: vector of subject-level strata membership for
#'      preliminary, unpruned tree
#'     \item prelimstratadefn: definition of final formed strata, in terms of
#'     covariates for preliminary, unpruned tree
#'     \item prelimtree: ctree object; preliminary, unpruned tree built by
#'     ctree, in terms of covariates
#' }
#'
#' @import partykit
#' @export
#'
fittrees = function(yy,X,family="cox",verbose=0,tree.hyper = tree_control()){

  ##====================================================##
  ## Step 3: Conditional Inference Trees to Find Strata ##
  ##====================================================##

  #ensuring X is in the proper format if only a single covariate is passed in
  if (is.null(dim(X))) X = data.frame(X)

  #new data frame with just filtered covariates
  if (family == "binomial" & !is.factor(yy)) yy = factor(yy)
  datdf = data.frame(yy,X)

  #extracting number of subjects
  nsubj = ifelse(is.null(nrow(yy)),length(yy),nrow(yy))

  #extracting key control parameters from tree_control: minbuckets and alpha,
  #as these are allowed to be different between preliminary (3A) and final (3B)
  #trees
  minbuckets   = tree.hyper$minbucket
  alpha        = tree.hyper$alpha

  #allowing different minbucket for original and pruning steps
  if (length (minbuckets)==1){
    minbucket1 = minbucket2 = minbuckets
  } else {
    minbucket1 = minbuckets[1]
    minbucket2 = minbuckets[2]
  }

  #allowing differing alphas for original and pruning steps
  if (length (alpha)==1){
    alpha1 = alpha2 = alpha
  } else {
    alpha1 = alpha[1]
    alpha2 = alpha[2]
  }

  #saving the rest of the control parameters to be passed into ctree_control()
  #function directly
  ctreeparms = tree.hyper[!(names(tree.hyper) %in%
                              c("minbucket","alpha"))]

  ##=========##
  ## Step 3A ##
  ##=========##

  #fitting ctree if there are any variables left after filtering
  if (length(colnames(X))>0){

    #running preliminary ctree algorithm
    if (family=="cox"){
      stree = do.call(partykit::ctree,c(list(Surv(yy[,1],yy[,2]) ~ ., data=datdf,
                                             minbucket=minbucket1,alpha=alpha1),
                                        ctreeparms))
    } else {
      stree = do.call(partykit::ctree,c(list(yy ~ ., data=datdf,
                                             minbucket=minbucket1,alpha=alpha1),
                                        ctreeparms))
    }

    #extract strata defintions (i.e., rules at terminal nodes)
    termNodes = .list.rules.party(stree)
    if (verbose > 1) print(paste0("unordered termNode ",1:length(termNodes),": ",
                              termNodes))
    termNodes = sapply(termNodes,simplifyand)
    if (verbose > 1) print(paste0("unordered simplified termNode ",
                              1:length(termNodes),": ",termNodes))

  } else stree=streeprune=NULL

  #continuing if any strata were found in step 3A
  if (exists("termNodes")){

    if (length(termNodes)>1){

      #determining which preliminary strata each subject belongs in
      #(even subjects with missing values for key covariates will be assigned a
      #node using this approach)
      pS = predict(stree,type="node")
      datdf2 = datdf
      if (length(pS) < nrow(datdf)){
        datdf2 = datdf[(rownames(datdf) %in% names(pS)),]
      }

      #ordering preliminary strata to input into step 3B
      if (family == "cox"){

        #calculating survival km within each strata
        sf=survfit(yy~pS,data=datdf2)

        #calculating RMST up to minimax time for each strata
        #(blinded to treatment assignment)
        ##GOES TO MAX TIME, NOT MAX OBSERVED TIME
        ##(but this seems to be the same as is used by survRM2)
        tau = min(sapply(unique(pS),function(x){
          max(datdf2$yy[,1][pS==x]) }))
        rmeans = summary(sf,rmean=tau)$table[,"*rmean"]

        #ordering risk order by RMST, saving for input into step 3B
        rmeanriskorder = order(rmeans)
        orderedstratalevels = sort(unique(pS))[rmeanriskorder]
        pS = ordered(pS,levels=orderedstratalevels)
        pS = ordered(as.numeric(pS))

      } else if (family == "binomial"){

        #below has issues when 2 nodes have same predicted probabilities
        #rmeans = unique(predict(stree,type="prob",simplify=TRUE))[,1]
        rmeans = predict(stree,type="prob",simplify=TRUE)[which(
          !duplicated(predict(stree, type = "node")))]
        #assuming higher probability = higher risk**
        rmeanriskorder = order(rmeans, decreasing = TRUE)#order(rmeans)
        orderedstratalevels = unique(pS)[rmeanriskorder]
        pS = ordered(pS,levels=orderedstratalevels)
        pS = ordered(as.numeric(pS))

      } else if (family == "gaussian"){

        rmeans = unique(predict(stree,type="response",simplify=TRUE))
        rmeanriskorder = order(rmeans)
        orderedstratalevels = unique(pS)[rmeanriskorder]
        pS = ordered(pS,levels=orderedstratalevels)
        pS = ordered(as.numeric(pS))

      }

      #reordering terminal nodes so output is given in terms of risk order (high -> low)
      termNodes = termNodes[as.character(orderedstratalevels)]

      ##=========##
      ## Step 3B ##
      ##=========##

      #post-pruning strata: fit ctree again w/ ordered strata as input variables
      if (family=="cox"){
        streeprune = do.call(partykit::ctree,c(list(
          Surv(yy[,1],yy[,2]) ~ pS, data=datdf2,minbucket=minbucket2,alpha=alpha2),
          ctreeparms))
      } else {
        streeprune = do.call(partykit::ctree,c(list(
          yy ~ pS, data=datdf2,minbucket=minbucket2,alpha=alpha2),ctreeparms))
      }

      #terminal node definition in terms of preliminary strata
      termNodes2 = .list.rules.party(streeprune)

      #continuing if any strata were found in step 3B
      if (length(termNodes2)>1){

        #----------------------------------------------------------------#
        # writing 3B terminal nodes in terms of initial input covariates #
        #----------------------------------------------------------------#
        #writing terminal nodes in a way to extract subjects fitting those criteria
        termNodesSep2 = strsplit(termNodes2,"&")
        termNodesSep2num = lapply(termNodesSep2,function(x)
          do.call(rbind,strsplit(x,"%in%"))[,2])
        termNodesSep2num = lapply(termNodesSep2num,function(b)
          if (length(b)>1){
            Reduce(intersect,lapply(b,function(x) eval(parse(text=x))))
          } else eval(parse(text=b))
        )
        termNodesSep2 = paste0("pS %in% ",termNodesSep2num)

        newNodeVars = lapply(lapply(termNodesSep2,function(x) gsub("\\D","-",x)),
                             function(y) strsplit(y,"-"))
        newNodeVars = lapply(newNodeVars,function(x)
          as.numeric(unlist(x)[unlist(x)!=""]))

        #adding robustness: pulling out by termNode names rather than position
        prunedTermNodes = unlist(lapply(newNodeVars,function(x)
          paste(paste0("(",termNodes[x],")"),collapse = " | ")))
          #paste(paste0("(",termNodes[names(termNodes) %in% x],")"),collapse = " | ")))
        prunedTermNodesFinal = gsub("(^c*)\\(","\\(datdf$",prunedTermNodes)
        prunedTermNodesFinal = gsub("& ","& datdf$",prunedTermNodesFinal)
        prunedTermNodesFinal = gsub("\\| \\(","\\| \\(datdf$",
                                    prunedTermNodesFinal)

        #extracting final node for each subject
        prunedpS = predict(streeprune,type="node")
        datdf3 = datdf2
        if (length(prunedpS) < nrow(datdf2)){
          datdf3 = datdf2[(rownames(datdf2) %in% names(prunedpS)),]
        }

      } else { # case where only one strata after pruning

        prunedTermNodes = prunedTermNodesFinal = "datdf"

      }

    } else{ #case where no initial subgroups are detected
      streeprune = stree#NULL
      termNodes = termNodesFinal = "datdf"
      prunedTermNodes = prunedTermNodesFinal = "datdf"
      pS = as.factor(rep(1,nsubj))
      prunedpS = rep(1,nsubj)

    }

  } else{ #case where no subgroups are detected or no tree was fit due to
          #nothing passing filtering step

    termNodes = termNodesFinal = "datdf"
    prunedTermNodes = prunedTermNodesFinal = "datdf"
    pS = as.factor(rep(1,nsubj))
    prunedpS = rep(1,nsubj)

  }

  #============================================================================#

  if (!is.data.frame(datdf)) datdf = as.data.frame(datdf)

  #extracting final strata (as list of strata each subj belongs to)
  treeStrata = as.numeric(factor(pS,ordered=FALSE,
                                 levels=sort(as.numeric(levels(pS)))))
  treeStrataPruned = as.numeric(factor(prunedpS,ordered=FALSE))

  #number of found strata
  nstrata = length(unique(treeStrata))
  nstrataPruned = length(unique(treeStrataPruned))

  #============================================================================#

  #-----------------------------------------------------------------------#
  # ensuring tree strata/terminal nodes are ordered high risk -> low risk #
  #-----------------------------------------------------------------------#

  if (family == "cox"){
    #order preliminary strata by risk determined by restricted mean (AUC of KM curve)
    sf=survfit(yy~treeStrata)
    if (length(unique(treeStrata))>1){
      tau.pre = min(sapply(unique(treeStrata),function(x){
        max(yy[,1][treeStrata==x]) }))
      rmeans = summary(sf,rmean=tau.pre)$table[,"*rmean"]
    } else rmeans = summary(sf,rmean="common")$table["*rmean"]
    riskOrder = order(rmeans)

    #order pruned strata by risk determined by restricted mean (AUC of KM curve)
    prunedsf=survfit(yy~treeStrataPruned)
    if (length(unique(treeStrataPruned))>1){
      tau.prune = min(sapply(unique(treeStrataPruned),function(x){
        max(yy[,1][treeStrataPruned==x]) }))
      prunedrmeans = summary(prunedsf,rmean=tau.prune)$table[,"*rmean"]
    } else prunedrmeans = summary(sf,rmean="common")$table["*rmean"]
    prunedRiskOrder = order(prunedrmeans)

  } else if (family == "binomial"|family=="gaussian"){

    #risk order by probability of case
    # below: same as: predict(glm(as.numeric(yy)~strata(treeStrata)-1,
    #              family=binomial(link=logit)),type="response") !!
    ts = treeStrata
    if (length(unique(ts))>1) ts = factor(ts)
    rmeans = glm(as.numeric(yy)~ts-1)$coef
    riskOrder = order(rmeans,decreasing=TRUE)

    tsp = treeStrataPruned
    if (length(unique(tsp))>1) tsp = factor(tsp)
    prunedrmeans = glm(as.numeric(yy)~tsp-1)$coef
    prunedRiskOrder = order(prunedrmeans,decreasing=TRUE)

  }

  #----------------------------------------------------------------#
  # ordering preliminary and final terminal nodes high -> low risk #
  #----------------------------------------------------------------#
  termNodes = termNodes[riskOrder]
  prunedTermNodes = prunedTermNodes[prunedRiskOrder]

  riskGroups = lapply(1:nstrata,function(x) which(treeStrata==riskOrder[x]))
  for (x in 1:length(riskGroups)){
    treeStrata[riskGroups[[x]]]=x
  }

  prunedRiskGroups = lapply(1:nstrataPruned,function(x)
    which(treeStrataPruned==prunedRiskOrder[x]))
  for (x in 1:length(prunedRiskGroups)){
    treeStrataPruned[prunedRiskGroups[[x]]]=x
  }

  #-------------------------#
  # outputing final results #
  #-------------------------#
  return(list(strataids=treeStrataPruned,stratadefn=prunedTermNodes,
              finaltree=streeprune,
              prelimstrataids=treeStrata,prelimstratadefn=termNodes,
              prelimtree=stree))

}


################################################################################


##==========================##
## 5-STAR Control Functions ##
##==========================##

#' Control for filtering step of 5-STAR algorithm
#'
#' Parameters for control of the elastic net or random forest filtering step
#' (step 2 of the 5-STAR algorithm)
#' @param method specifying method used for filtering; either "ENET" for
#' elastic net, "RF" for random forest using defaults in randomForestSRC package
#' (with possible modifications through "..."), or "RFbest" for tuned random forest
#' (slower but more optimized)
#' @param lambdatype Optional elastic net parameter; whether to use lambda that
#' minimizes cross validation error (lambdatype="min", default) or the largest
#' lambda that gives error within 1 standard error of minimum error
#' (lambdatype="1se"). Ignored when method is "RF" or "RFbest"
#' @param mixparm Optional elastic net mixing parameter alpha or grid of alpha
#' values to search over. If nothing is entered, will search for best value
#' between 0.05 and 0.95. Ignored when method is "RF" or "RFbest"
#' @param vimpalpha For "RF" or "RFbest", the significance level for variable
#' importance (VIMP) confidence intervals (CIs) to determine which covariates
#' are passed through the filtering stage (default = 0.05, ignored for
#' method == "ENET")
#' @param nfolds Number of folds used for cross validation when tuning the elastic
#' net tuning parameter lambda (Default = 10; ignored when method is "RF" or
#' "RFbest")
#' @param filterseed Seed for the filtering step to control variability in
#' cross validation step of the elastic net filtering
#' @param ... Optional additional arguments passed into
#' \code{\link[glmnet]{glmnet}} or \code{\link[randomForestSRC]{rfsrc}} functions
#' @return A list of control parameters for filtering step
#' @export
filter_control = function(method="ENET",lambdatype="min",mixparm=NULL,
                          vimpalpha=0.05,nfolds=10,filterseed=2019,...){
  list(method=method, lambdatype=lambdatype, mixparm=mixparm,
       vimpalpha=vimpalpha, nfolds=nfolds, filterseed=filterseed,...)
}


#' Control for strata-finding tree step of 5-STAR algorithm
#'
#' Parameters for control of conditional inference tree steps of 5-STAR algorithm
#' (3A and 3B) - mainly for passing into \code{\link[partykit]{ctree}} and
#' \code{\link[partykit]{ctree_control}} functions
#' @param minbucket Vector of minimum set of weights per terminal node for
#' initial and pruning steps (e.g., minimum number of patients/terminal node for
#' steps 3A and 3B) - if a single number is given, the same value is
#' used for both preliminary and final trees (Default = max(50,0.05*nsubj) for both steps)
#' @param alpha vector of significance level for variable selection for tree
#'  splits in preliminary and final trees (3A and 3B) - if a single number is
#'  given, the same value is used for both preliminary and final trees
#'  (Default is (0.1,0.1))
#' @param testtype from \code{\link[partykit]{ctree_control}}:
#' "a character specifying how to compute the distribution of the test
#' statistic" (default = "Bonferroni")
#' @param majority ctree control parameter, specifying whether to randomly
# assign subjects with missing information in splitting variable (FALSE) or to
# go with the majority (TRUE) (default = FALSE)
#' @param maxsurrogate ctree control parameter defining number of surrogate
# splits to evaluate for missing covariates (see ctree_control)
#' @param maxdepth Maximum tree depth (default for 5-STAR is 2)
#' @param ... additional parameters to be passed into
#' \code{\link[partykit]{ctree}} function
#' @return A list of control parameters for strata formation step
#' @export
tree_control = function(minbucket=NULL,alpha=c(0.1,0.1),testtype="Bonferroni",
                        majority=FALSE,maxsurrogate=3,maxdepth=2,...){
  list(minbucket=minbucket,alpha=alpha,testtype=testtype,
       majority=majority,maxsurrogate=maxsurrogate,maxdepth=maxdepth,...)
}

################################################################################

#' Prepares data for 5-STAR algorithm
#'
#' Cleans covariate matrix, removing covariates with too much missingess or those
#' that can't be split on (i.e., due to too few obs/minor levels),
#' and coverts all character covariates into factors
#' @param yy Trait/response (a Surv() object summarizing follow-up time for
#' right-censored data and status indicator where 1=dead, 0=censored).
#' Ignored for family != "cox"
#' @param X Data frame of all possible stratification covariates
#' @param family Trait family, current options: "cox", "binomial", or "gaussian"
#' @param missthreshold Vector of lower and upper bound of acceptable missingness
#' levels for each covariate, such that covariates with less than the first
#' missthreshold value will be passed to step 2 (filtering step), those with
#' missingess greater than the second missthreshold value will be removed from
#' analysis, and those with missingness between these two values will be included
#' only if they are significantly correlated with the outcome scores
#' (e.g., logrank scores). If a scalar is entered, covariates with less than
#' that amount of missingness will be included and those with greater will be
#' removed. For family other than "cox", only the first value is used.
#' @param minbucket Minimum number of subjects/terminal node
#' @param verbose Numeric variable indicating amount of information to print
#' to the terminal (0 = nothing, 1 = notes only, 2 = notes and intermediate output)
#' @return \itemize{
#' \item X: cleaned covariate matrix
#' }
#'
#' @export
prep5STAR = function(yy,X,family="cox",missthreshold=c(0.1,0.2),verbose=0,
                     minbucket){

  #------------------------------------------#
  # converting character variables to factor #
  #------------------------------------------#

  ischar = sapply(X,is.character)
  if (sum(ischar)>0){
    if (verbose > 0) print(paste0("Note: character variable(s) ",
                 paste(names(X[ischar]),collapse=", "),
                 " being converted to factors"))
    if (sum(ischar)==1){
      X[,ischar] = sapply(X[,ischar],factor)
    } else X[,ischar] = lapply(X[,ischar],factor)
  }

  #---------------------------------------------------#
  # check that valid missingness threshold is entered #
  #---------------------------------------------------#

  #checking correct length
  if (length(missthreshold) > 1 & family!="cox"){
    if (verbose > 0) print(
      "Note: Only the first value of missthreshold will be used.")
    missthreshold = missthreshold[1]
  }
  if (length(missthreshold) > 2){
    if (verbose > 0) print(
      "Note: Only the first two values of missthreshold are used.")
    missthreshold = missthreshold[1:2]
  }
  if (length(missthreshold) == 1) missthreshold = rep(missthreshold,2)

  if (any(missthreshold < 0)) stop("Missingness threshold must be non-negative")
  if (all(missthreshold > 1)){
    if (all(missthreshold <= 100)){
      print(paste0("Setting missinnesss to proportion ",missthreshold/100," (",
            missthreshold,"/100)"))
      missthreshold = misshthreshold/100
    } else stop("Missingness threshold must be a proportion in [0,1]")
  }

  #-----------------------------------------------#
  # removing covariates with too much missingness #
  #-----------------------------------------------#

  #removing covariates with > upper boundary missingness
  covs2rm = which(colMeans(is.na(X)) > missthreshold[2])
  if (length(covs2rm)>0) X = X[,-covs2rm]

  #for survival data, using rank correlation to determine which covariates
  #with missingness between lower and upper limits to keep
  if (family == "cox"){

    #calculating logrank scores
    LRscores = coin::logrank_trafo(yy)
    #matrix of covariates to check (i.e., missingness in
    #  (missthreshold[1],missthreshold[2]))
    Xmiss = X[,colMeans(is.na(X)) > missthreshold[1]]

    if (length(Xmiss) > 0){
      #(more than one covariates with missingness between boundaries)
      if (!is.null(ncol(Xmiss))){

        #converting ordinal, binary variables to numeric in test matrix
        #(same for rank-based correlation)
        isOrdered = sapply(Xmiss,is.ordered)
        Xmiss[isOrdered]=lapply(Xmiss[isOrdered],as.numeric)
        isBinary = sapply(Xmiss,function(x) is.factor(x)&(length(levels(x))==2))
        Xmiss[isBinary] = sapply(Xmiss[isBinary],function(x)
          (x==levels(x)[1])*1 )

        #performing test of rank correlation
        Xmiss.num = Xmiss[,sapply(Xmiss,is.numeric)]

        if (!is.null(ncol(Xmiss.num))){
          cortest = sapply(Xmiss[,sapply(Xmiss,is.numeric)],function(x)
            cor.test(LRscores,x,method="kendall",alternative="two.sided")$p.value)
        } else {
          cortest = cor.test(LRscores,Xmiss.num,method="kendall",
                             alternative="two.sided")$p.value
          names(cortest) = names(which(sapply(Xmiss,is.numeric)))
        }

        #performing kruskal wallis h test (rank one-way ANOVA) for nominal variables
        Xmiss.factor = Xmiss[,sapply(Xmiss,function(y) !is.numeric(y))]
        if (!is.null(ncol(Xmiss.factor))){
          if (ncol(Xmiss.factor)>0){
            cortest2 = sapply(Xmiss.factor,function(x){
              LRscores.nomiss = LRscores[!is.na(x)]
              x.nomiss = x[!is.na(x)]
              kruskal.test(x=LRscores.nomiss,g=x.nomiss)$p.value
            })
          #still shows as non-null if an empty matrix (e.g., 0x6...)
          } else cortest2 = NULL
        } else {
          cortest2 = kruskal.test(x=LRscores[!is.na(Xmiss.factor)],
                                  g=Xmiss.factor[!is.na(Xmiss.factor)])$p.value
          names(cortest2) = setdiff(names(Xmiss),names(cortest))
        }

        cortest.all = c(cortest,cortest2)

        #removing those not "significantly" correlated
        torm = names(cortest.all)[cortest.all >= 0.05]

      } else { #case where only one covariate with missingness between boundaries

        #converting ordinal, binary variables to numeric
        if (is.factor(Xmiss)&(length(levels(Xmiss))==2)){
          Xmiss = (Xmiss==levels(Xmiss)[1])*1
        } else if (is.ordered(Xmiss)){
          Xmiss = as.numeric(Xmiss)
        }

        if (is.numeric(Xmiss)){
          #performing test of rank correlation
          cortest = cor.test(LRscores,Xmiss,method="kendall",
                             alternative="two.sided")$p.value
        } else { #nomial covariate
          #performing kruskal wallis h test (rank one-way ANOVA)
          cortest = kruskal.test(LRscores,Xmiss)$p.value
        }

        torm = ifelse(cortest >= 0.05,names(which(colMeans(
          is.na(X)) > missthreshold[1])),NA)

      }

      if (length(torm[!is.na(torm)]) > 0) X = X[,!(colnames(X) %in% torm)]
    } else torm = NA #nothing with > missthreshold[1] missingness -- nothing to do here

  } else torm = NA

  #names of covariates to remove
  covs2rm.all = c(names(covs2rm),torm[!is.na(torm)])

  if ((length(covs2rm.all) > 0) & (verbose > 0)){
      print(paste0("Removing covariate(s) ",paste(covs2rm.all,collapse=", "),
                   " from analysis due to too much missingness"))
  }

  # also removing covariates that couldn't possibly be split on
  # (e.g., those with < minbucket in all the less frequent groups combined)

  #removing variables w/ no hope of splitting (e.g., w/ < minbucket in all
  #non-major groups)
  covs.cansplit = sapply(X,function(x){
    tab=table(x)
    sum(tab[names(tab) != names(which.max(tab))])>=minbucket
  })
  covs2rm.cantsplit = names(X)[!covs.cansplit]
  X = X[,covs.cansplit]
  if ((length(covs2rm.cantsplit) > 0) & (verbose > 0)){
    print(paste0("Removing covariate(s) ",paste(covs2rm.cantsplit,collapse=", "),
                 " from analysis as no splits are possible with these covariates",
                 " (number of subjects in non-major groups <",minbucket,")"))
  }

  #removing variables w/ > 31 categories that are UNordered factors
  #(CTree can't handle this type of data)
  covs.manyunordlevels = sapply(X,function(x){
    (is.factor(x))&(!is.ordered(x))&(length(levels(x))>31)
  })
  covs2rm.manyunordlevels = names(X)[covs.manyunordlevels]
  X = X[,!(names(X) %in% covs2rm.manyunordlevels)]
  if ((length(covs2rm.manyunordlevels) > 0) & (verbose > 0)){
    print(paste0("Removing covariate(s) ",paste(covs2rm.manyunordlevels,
                                              collapse=", "),
                 " from analysis as CTree can't find splits in unordered",
                 " variables with > 31 levels."))
  }

  # returning cleaned covariate matrix
  return(X)
}

################################################################################

#---------------------#
# cleaning node names #
#---------------------#
#(to add more cleaning later)
cleanNodeNames = function(nodename){

  #R code to human-readable text
  newname = gsub('%in% c\\(\"',"in {",nodename)
  newname = gsub('\"\\)',"}",newname)
  newname = gsub('\\"',"",newname)

  #if only one, change "in {}" to "="
  equalTerms = unlist(regmatches(newname, gregexpr("in \\{.+?\\}", newname)))
  if (length(equalTerms) > 0){

    equalTermsReplace = sapply(equalTerms, function(x){
      if (length(strsplit(x,",")[[1]]) == 1){
        x = stringr::str_replace(x,"in \\{(.+?)\\}","= \\1")
      }
      return(x)
    })

    for (i in 1:length(equalTerms)){
      newname = sub(equalTerms[i],equalTermsReplace[i],newname,fixed=TRUE)
    }

  }

  #finally, replacing "&" with "," for cleaner naming
  newname = gsub(" & ",", ",newname)

  return(newname)
}


################################################################################

##==================##
## 5-STAR Main Code ##
##==================##

#' Run 5-STAR Algorithm
#'
#' Runs 5-STAR algorithm for forming homogeneous risk strata and calculating
#' an overall treatment effect estimate amalgamated from the formed strata
#'
#' @param yy      Trait/response (currently must be either continous, binary, or
#' a Surv() object summarizing follow-up time for right-censored data and status
#'  indicator where 1=dead, 0=censored)
#' @param arm     Binary treatment indicator, 1 = test treatment, 0 = control
#' @param X       Data frame of all covariates (as specified in Step 1 of 5-STAR)
#' @param family  Trait family, current options: "cox", "binomial", and
#' "gaussian" (note: binomial and gaussian are still experimental!)
#' @param measure Response of interest; current options
#' are: for survival traits: "HR" (hazard ratio, default when family = "cox")
#' and "TR" (time ratio from model averaging of AFT models);
#' for binary traits: "RD" (risk difference, default when family = "binomial");
#' for continuous traits: "MD" (mean difference, default when family = "gaussian")
#' @param cilevel Confidence level alpha for overall result and confidence
#' intervals (default = 0.025, for one-tailed tests)
#' @param alternative For tests, whether alternative hypothesis is "less",
#'  "greater", or "two.sided". Currently ignored for family!="cox" (all others
#'  assume a 2-sided test for now)
#' @param vartype Whether stratum-level variances used to calculate the correlation
#' between sample size ni and ni/sqrt(Vi) test statistics should be calculated
#' under the null ("null") or estimated via the Cox PH model ("alt"). "alt"
#' is default, and must be used when measure != "HR"
#' @param timeunit  Optional parameter for plots defining scale of time to event
#'  data (e.g., "Months", "Years"); ignored for family!="cox" and plot == FALSE
#' @param tau Optional numeric variable specifying end time point for rmst
#' calculations; if NULL, will calculate the minimum of maximum times over arm
#' within each formed stratum; ignored when measure != "RMST"
#' @param inclfrailty A logical variable - whether or not to include a frailty
#' term to each within-strata Cox PH model fit for added robustness against
#' unexplained heterogeneity. Ignored when measure != "HR".
#' @param verbose   Numeric variable indicating amount of information to print
#' to the terminal (0 = nothing, 1 = notes only, 2+ = notes and intermediate output)
#' @param plot  Logical, whether to return summary plots (filtering plots,
#' within strata and pooled KM plots, and forest plots summarizing strata-level
#' and overall results)
#' @param filter.hyper A list of control parameters for filtering step (Step 2),
#' (see \code{\link{filter_control}})
#' @seealso \code{\link{filter_control}}
#' @param tree.hyper A list of control parameters for strata building step (Step 3)
#' (see \code{\link{tree_control}})
#' @seealso \code{\link{tree_control}}
#' @param distList Vector of models (survival distributions for parametric AFT
#' fit) to include in model averaging; Each element must be the name of an
#' element from \code{\link[survival]{survreg.distributions}} (see also
#' \code{\link[survival]{survreg}}); default is c("weibull","lognormal","loglogistic");
#' ignored when family != "cox"
#' @param ucvar Estimator for the unconditional variance of the model averaging
#' estimate to use. 1 uses Buckland et al. (1997) analytical estimator, and 2 uses
#' the more conservative estimator of Burnham and Anderson (2002) related to
#' Bayesian model averaging. For more details see Turek 2013. Default is "1".
#' @param shading For plotting, whether to add shaded confidence intervals
#' around the pooled KM curves to better differentiate how similar the survival
#' curves are. Default is FALSE. Ignored when plot = FALSE and/or family != "cox"
#' @param missthreshold
#' @param fplottype
#'
#' @details Filtering step is performed on complete case basis (e.g., removing
#' all individuals with any missing covariate data)
#' @return \itemize{
#'     \item inputcovs: vector containing names of all covariates input into
#'     the algorithm (after removing those with too much missingness, as defined
#'     by missthreshold, or can't be split on) (see \code{\link{prep5STAR}} for details)
#'     \item filteredcovs: vector containing names of all covariates left after
#'     filtering step (step 2)
#'     \item filterdetails: list containing additional output from filter step;
#'     for method=ENET, this includes: \itemize{
#'    \item cvout: containing fraction of null deviance explained, mean cross validation
#' error, and optimal lambda and alpha value (see cv.glmnet and glmnet for
#' more details)
#'    \item beta: the coefficients of glmnet fit given tuned choice of alpha and lambda
#'    \item ENETplot: plots of deviance and solution path over alpha/lambda parameters
#' } For method=RF or RFbest, this includes:
#' \itemize{
#'    \item varselectmat: matrix of VIMP confidence intervals, p-values, and selection
#' decision for each variable
#'    \item VIMPplot: default variable importance CI plot, output if plot = TRUE
#'    \item varselect2plot: variable selection information from subsample.rfsrc for
#' making customizable VIMP CI plots
#'    \item forest: rfsrc object
#' }
#'     \item prelimtree: ctree object; preliminary tree (3A) built by ctree,
#'     in terms of covariates
#'     \item prelimstratadefn: vector of rules defining the preliminary strata
#'     (step 3A), as defined by covariates chosen in initial ctree run, ordered from
#'     highest to lowest risk
#'     \item prelimstrataids: vector of preliminary (3A) strata membership for
#'     each subject (where 1 indicates highest risk, and the largest number
#'     indicates lowest risk, i.e., matching order of prelimstratadefn)
#'     \item prelimstratasurvfits: a \code{\link[survival]{survfit}} object,
#'     summarizing survival curves in each preliminary (3A) strata, pooled by
#'     treatment assignment; returned when family=="cox"
#'     \item prelimbystratasummary: matrix summarizing  model fits within
#'     each preliminary (3A) strata, including estimated (log) hazard ratio
#'     (or TR, OR, RD, etc. when measure!="HR"), corresponding (1-cilevel)x100\%
#'     CI, Gramsch and Therneau test of non-PH p-value, and sample size,
#'     inverse variance, and adaptive weights
#'     \item prelimsinglewtres: Preliminary results using only a single
#'     weighting scheme (e.g., just ni sample size weights or just ni/Vi weights)
#'     \item prelimres5star: Summary of amalgamated estimand, variance,
#'     p-value, and (1-cilevel)x100\% confidence interval for adaptive weights
#'     from preliminary strata (3A)
#'     \item finaltree: ctree object; final, pruned tree (3B) built by ctree,
#'     in terms of preliminary strata
#'     \item stratadefn: vector of rules defining the final strata
#'     (step 3B), as defined by covariates translated from pooling ctree run,
#'     ordered from highest to lowest risk
#'     \item strataids:  vector of final (3B) strata membership for
#'     each subject (where 1 indicates highest risk, and the largest number
#'     indicates lowest risk, i.e., matching order of stratadefn)
#'     \item stratasurvfits: a \code{\link[survival]{survfit}} object,
#'     summarizing survival curves in each final (3B) strata, pooled by
#'     treatment assignment; returned when family=="cox"
#'     \item bystratasummary: matrix summarizing  model fits within
#'     each final (3B) strata, including estimated (log) hazard ratio (or OR,
#'     RD, etc. when measure!="HR"), corresponding (1-cilevel)x100\% CI,
#'     Gramsch and Therneau test of non-PH p-value, and sample size,
#'     inverse variance, and adaptive weights
#'     \item singlewtres: Final results using only a single weighting scheme
#'     (e.g., just ni sample size weights or just ni/Vi weights)
#'     \item res5star: Summary of amalgamated estimand, variance,
#'     p-value, and cilevelx100% confidence interval for adaptive weights
#'     \item prelimbetweenstrataKM: Kaplan-Meier survival curve plots from pooled
#'      treatment assignment data from each preliminary (3A) strata,
#'      returned if plot == TRUE and family == "cox"
#'     \item prelimbystrataKM: Kaplan-Meier survival curve plots within each
#'     preliminary (3A) strata, returned if plot == TRUE and family == "cox"
#'     \item prelimforestplot: a forest plot summarizing estimated coefficient and
#'     confidence interval for each preliminary (3A) strata, as well as amalgamated
#'     result, returned when plot == TRUE and family == "cox"
#'     \item betweenstrataKM: Kaplan-Meier survival curve plots from pooled
#'      treatment assignment data from each final (3B) strata,
#'      returned if plot == TRUE and family == "cox"
#'     \item bystrataKM: Kaplan-Meier survival curve plots within each final
#'     (3B) strata, returned if plot == TRUE and family == "cox"
#'     \item forestplot: a forest plot summarizing estimated coefficient and
#'     confidence interval for each final (3B) strata, as well as amalgamated
#'     result, returned if plot == TRUE and family == "cox"
#' }
#'
#' @import survival
#' @export
run5STAR = function(yy,arm,X,family="cox",#measure="HR",
                    measure = ifelse(family=="cox","HR",
                                     ifelse(family=="binomial","RD","MD")),
                    alternative=ifelse(measure %in% c("HR","RD"),"less","greater"),
                    cilevel=ifelse(alternative=="two.sided",0.05,0.025),
                    vartype="alt",missthreshold=c(0.1,0.2),timeunit=NULL,
                    tau=NULL,inclfrailty=FALSE,verbose=0,plot=TRUE,
                    filter.hyper=filter_control(),
                    tree.hyper = tree_control(),
                    distList=c("weibull","lognormal","loglogistic"),
                    ucvar=1,shading=FALSE,fplottype="overall",
                    vars2keep=NULL){

  # if minbucket is not set through tree_control(), setting it to the default
  # value here to allow the default to be data-dependent
  if (is.null(tree.hyper$minbucket)){
    tree.hyper$minbucket <- max(ceiling(0.05*nrow(X)),50)
  }

  #setting measure if it is missing
  #(currently ignored as a default is set)
  # if (is.null(measure)){
  #
  #   if (family == "cox"){
  #     measure = "HR"
  #   } else if (family == "gaussian"){
  #     measure = "MD"
  #   } else if (family == "binomial"){
  #     measure = "RD"
  #   } else stop("family must be one of 'cox', 'gaussian', or 'binomial'.")
  #
  # }
  if (!(family %in% c("cox","binomial","gaussian"))) stop(
    "family must be one of 'cox','binomial', or 'gaussian'.")

  message(paste0("Running 5-STAR algorithm with ",family," family and ",measure,
                 " measure. Alternative hypothesis direction: ",
                 ifelse(alternative == 'two.sided',"two.sided",paste0("test ",
                 alternative," than control."))))

  #checking family/measure compatibility
  if (family!="cox"){
    warning("Currently family!='cox' is still experimental!")
  } else if (!(measure %in% c("HR","TR"))) {
    warning("Currently methodology is highly experimental for method!='TR' or 'HR'!")
  }

  #need to specify measure if family is binomial or cox
  if (family=="binomial" & !(measure %in% c("RD","OR"))){
    stop("For binary traits, measure must be either `RD` or `OR`.")
  }
  if (family=="cox" & !(measure %in% c("HR","TR","RMST"))){
    stop(paste0("For survival traits, measure must be either 'HR' ",
                "(hazard ratio), 'TR' (time ratio from AFT model) or 'RMST' ",
                "(restricted mean survival time)."))
  }

  #making sure "alt" variance is selected when measure != "HR"
  if (measure!="HR" & vartype=="null"){
    stop(paste0('Correlation between test statistics must be calculated using',
                'the alternative variance (vartype="alt") unless measure == "HR".'))
  }

  #============================================================================#

  #------------------------------------------------------#
  # extracting/comining dat of interest to usable format #
  #------------------------------------------------------#

  #removing missing values in response
  obs2rm = is.na(yy)
  yy = yy[!obs2rm]
  X = X[!obs2rm,]
  arm = arm[!obs2rm]

  if (verbose > 1 & sum(obs2rm > 0)){
    print(paste0("Removing ",sum(obs2rm)," missing observations."))
  }

  #converting yy to a 1 column matrix if it is a vector
  if (is.null(nrow(yy))) yy = matrix(yy,ncol=1)

  #extracting all data, number of subjects, potential covariate matrix
  dat = cbind(yy,arm,X)
  nsubj = nrow(yy)
  X = as.data.frame(X)

  #============================================================================#

  #----------------------------------------#
  # checking for missingness in covariates #
  #----------------------------------------#

  X = prep5STAR(yy=yy,X=X,family=family,missthreshold=missthreshold,
                minbucket=tree.hyper$minbucket[1],verbose=verbose)

  #input covariate names
  inputcovs = names(X)

  #============================================================================#

  #---------------------------#
  # Step 2: Filter Covariates #
  #---------------------------#

  filterRes = filter5STAR(yy=yy,X=X,family=family,plot=plot,
                          verbose=verbose,filter.hyper=filter.hyper,vars2keep=NULL)

  cov2keep.all = filterRes$cov2keep
  X = X[,colnames(X) %in% cov2keep.all]

  if (length(cov2keep.all)==1){
    X = data.frame(X)
    colnames(X) = cov2keep.all
  }

  if (verbose > 1){
    print(paste0("Covariates kept after filtering step: ",
                 paste(colnames(X),collapse=",")))
  }

  #other details, plots, etc. to output
  filterdetails = filterRes[names(filterRes) != "cov2keep"]

  #============================================================================#

  #----------------------------------------------------#
  # Step 3: Conditional Inference Trees to Find Strata #
  #----------------------------------------------------#

  trees5STAR = fittrees(yy=yy,X=X,family=family,verbose=verbose,
                        tree.hyper=tree.hyper)

  treeStrataPruned = trees5STAR$strataids
  prunedTermNodes = trees5STAR$stratadefn
  streeprune = trees5STAR$finaltree
  treeStrata = trees5STAR$prelimstrataids
  termNodes = trees5STAR$prelimstratadefn
  stree = trees5STAR$prelimtree

  #cleaning node names
  prunedTermNodes = sapply(prunedTermNodes,cleanNodeNames)
  termNodes = sapply(termNodes,cleanNodeNames)
  names(prunedTermNodes) = names(termNodes) = NULL

  if (verbose > 1){

    if (length (termNodes) > 1){
      print(paste0("Preliminary Strata: ",termNodes))
    } else print("No strata were formed at end of step 3A. Preliminary stratum consists of the whole data set.")

    if (length (prunedTermNodes) > 1){
      print(paste0("Final Pruned Strata: ",prunedTermNodes))
    } else print("No strata were formed at end of step 3B. Final pruned stratum consists of the whole data set.")

  }

  #============================================================================#

  #----------------------------------------------------------------------#
  # calculating number of events per strataxarm                          #
  # can't calculate an estimate if not at least one event per strataxarm #
  #----------------------------------------------------------------------#

  statusVec = yy
  if (is.factor(statusVec)) statusVec = as.numeric(as.character(statusVec))

  if (family == "cox") statusVec = yy[,2]

  if (family!="gaussian"){

    #for survival traits, need at least one event per strata x arm
    perStrataEventSummary = eventsbystrata(treeStrata=treeStrata,arm=arm,
                                           status=statusVec,family=family)
    perStrataEventSummaryPruned = eventsbystrata(treeStrata=treeStrataPruned,
                                                 arm=arm,status=statusVec,
                                                 family=family)

    ###skipping if not at least 1 event per strata x arm combination
    minEventsPerStratArm = min(perStrataEventSummary[1,])
    minEventsPerStratArmPruned = min(perStrataEventSummaryPruned[1,])

    #(should have at least as many events/group as in initial tree fit)
    if (minEventsPerStratArmPruned < 1)
      stop(paste0("Need at least one event per strata x arm",
                  " (events per preliminary strata: ",
                  paste0(sapply(1:length(unique(treeStrata)),function(x)
                    sum(statusVec[treeStrata==x])),collapse="-"),
                  " and per final strata: ",
                  paste0(sapply(1:length(unique(treeStrataPruned)),function(x)
                    sum(statusVec[treeStrataPruned==x])),collapse="-"),")"))


    if (family == "binomial"){
      if (min(perStrataEventSummary[3,])==0)
        stop("Need at least one control per strata x arm")
    }

  } else {
    minEventsPerStratArm = NULL
    minEventsPerStratArmPruned = NULL
  }

  #============================================================================#

  #---------------------------------------------#
  # Step 4: Fit Model Within Each Formed Strata #
  #---------------------------------------------#

  if (family == "cox" & measure == "HR"){

    time = yy[,1]; status = yy[,2]
    ByStrataFits = coxbystrata(time,status,arm,treeStrata,termNodes,
                               treetype="preliminary",alternative=alternative,
                               cilevel,inclfrailty=inclfrailty,verbose,plot,
                               timeunit,shading=shading)

    MRSSmat = ByStrataFits$fitsummary
    p3 = ByStrataFits$bystrataKM
    p4 = ByStrataFits$betweenstrataKM

    prunedByStrataFits = coxbystrata(time=time,status=status,arm=arm,
                                     treeStrata=treeStrataPruned,
                                     termNodes=prunedTermNodes,
                                     treetype="final",alternative=alternative,
                                     cilevel=cilevel,inclfrailty=inclfrailty,
                                     verbose,plot,timeunit,shading=shading)

    prunedMRSSmat = prunedByStrataFits$fitsummary
    p3prune = prunedByStrataFits$bystrataKM
    p4prune = prunedByStrataFits$betweenstrataKM

  } else if (family == "cox" & measure == "RMST"){

    time = yy[,1]; status = yy[,2]
    ByStrataFits = rmstbystrata(time,status,arm,tau=tau,treeStrata,termNodes,
                               treetype="preliminary",alternative=alternative,
                               cilevel=cilevel,verbose=verbose,plot=plot,
                               timeunit=timeunit)

    MRSSmat = ByStrataFits$fitsummary
    p3 = NULL#ByStrataFits$bystrataKM
    p4 = NULL#ByStrataFits$betweenstrataKM

    prunedByStrataFits = rmstbystrata(time=time,status=status,arm=arm,tau=tau,
                                     treeStrata=treeStrataPruned,
                                     termNodes=prunedTermNodes,
                                     treetype="final",alternative=alternative,
                                     cilevel=cilevel,verbose=verbose,
                                     plot=plot,timeunit=timeunit)

    prunedMRSSmat = prunedByStrataFits$fitsummary
    p3prune = NULL#prunedByStrataFits$bystrataKM
    p4prune = NULL#prunedByStrataFits$betweenstrataKM

  } else if (family=="cox" & measure=="TR"){

    time = yy[,1]; status = yy[,2]
    ByStrataFits = maaftbystrata(time,status,arm,treeStrata=treeStrata,
                                 termNodes=termNodes,distList=distList,
                                 ucvar=ucvar,alternative=alternative,
                                 cilevel=cilevel,verbose=verbose,
                                 plot=plot,treetype="preliminary",
                                 timeunit=timeunit,shading=shading)

    MRSSmat = ByStrataFits$fitsummary
    p3 = ByStrataFits$bystrataKM
    p4 = ByStrataFits$betweenstrataKM

    prunedByStrataFits = maaftbystrata(time=time,status=status,arm=arm,
                                       treeStrata=treeStrataPruned,
                                       termNodes=prunedTermNodes,
                                       distList=distList,ucvar=ucvar,
                                       alternative=alternative,
                                       cilevel=cilevel,verbose,
                                       plot=plot,treetype="final",
                                       timeunit=timeunit,shading=shading)

    prunedMRSSmat = prunedByStrataFits$fitsummary
    p3prune = prunedByStrataFits$bystrataKM
    p4prune = prunedByStrataFits$betweenstrataKM

  } else if (family == "gaussian" & measure == "MD"){

    ByStrataFits = mdbystrata(yy=yy,arm=arm,treeStrata=treeStrata,
                              treetype='preliminary',
                              termNodes=termNodes,cilevel=cilevel,
                              verbose=verbose,alternative=alternative,
                              plot=plot)

    MRSSmat = ByStrataFits$fitsummary

    prunedByStrataFits = mdbystrata(yy=yy,arm=arm,treeStrata=treeStrataPruned,
                                    treetype='final',
                                    termNodes=prunedTermNodes,cilevel=cilevel,
                                    verbose=verbose,alternative=alternative,
                                    plot=plot)

    prunedMRSSmat = prunedByStrataFits$fitsummary

    p3 = ByStrataFits$p3
    p4 = ByStrataFits$p4
    p3prune = prunedByStrataFits$p3
    p4prune = prunedByStrataFits$p4

    ## no longer allowing "OR"/glm-based option for measure
  # } else if (family=="gaussian"|(family=="binomial" & measure=="OR")){
  #
  #   ByStrataFits = glmbystrata(yy,arm,family,treeStrata,termNodes,cilevel,
  #                              verbose,alternative=alternative)
  #   MRSSmat = ByStrataFits$fitsummary
  #
  #   prunedByStrataFits = glmbystrata(yy=yy,arm=arm,family=family,
  #                                    treeStrata=treeStrataPruned,
  #                                    termNodes=prunedTermNodes,
  #                                    cilevel=cilevel,verbose=verbose,
  #                                    alternative=alternative)
  #
  #   prunedMRSSmat = prunedByStrataFits$fitsummary
  #
  #   p3 = p4 = p3prune = p4prune = fplot = fplotprune = NULL

  } else if (family == "binomial" & measure == "RD"){

    ByStrataFits = ratediffbystrata(yy=yy,arm=arm,
                                    treeStrata=treeStrata,
                                    treetype="preliminary",termNodes=termNodes,
                                    alternative=alternative,
                                    cilevel=cilevel,verbose,plot)
    MRSSmat = ByStrataFits$fitsummary

    prunedByStrataFits = ratediffbystrata(yy=yy,arm=arm,
                                          treeStrata=treeStrataPruned,
                                          # treetype="final",
                                          termNodes=prunedTermNodes,
                                          alternative=alternative,
                                          cilevel=cilevel,verbose,plot)

    prunedMRSSmat = prunedByStrataFits$fitsummary

    #currently no plots output for family != "cox"
    p3 = p4 = p3prune = p4prune = NULL

  }

  #============================================================================#

  #-----------------------------------------------------#
  # Step 5: Amalgamate Results, using minP approach     #
  # accounting for correlation between the two possible #
  # test statistics (ni, ni/sqrtVi weights)             #
  #-----------------------------------------------------#

  adaptiveWtResSS = minPadapHR(sf=ByStrataFits$stratafit, betas=MRSSmat[,1],
                               vars=MRSSmat[,2], cilevel=cilevel,
                               alternative=alternative, vartype=vartype)

  adaptiveWtResSSprune = minPadapHR(sf=prunedByStrataFits$stratafit,
                                    betas=prunedMRSSmat[,1],
                                    vars=prunedMRSSmat[,2],
                                    cilevel=cilevel,alternative=alternative,
                                    vartype=vartype)

  #combining new weights in coxmrss output
  MRSSmat = data.frame(MRSSmat,weight.adap=adaptiveWtResSS$adaptivewts)
  prunedMRSSmat = data.frame(prunedMRSSmat,
                             weight.adap=adaptiveWtResSSprune$adaptivewts)

  #storing results from adaptive weight test
  res5starprelim = adaptiveWtResSS$adaptivewtRes
  res5star = adaptiveWtResSSprune$adaptivewtRes

  #single-weight results
  singlewtresprelim = adaptiveWtResSS$singlewtres
  singlewtres = adaptiveWtResSSprune$singlewtres

  #don't need exp results if trait is quantitative or measure is RD
  if (family == "gaussian"|(family=="binomial" & measure=="RD")|(
    family == "cox" & measure == "RMST")){
    res5starprelim = res5starprelim[!(names(res5starprelim) %in%
                                        c("exp(beta)", "exp(ci lower)",
                                          "exp(ci upper)"))]
    res5star = res5star[!(names(res5star) %in%
                            c("exp(beta)", "exp(ci lower)","exp(ci upper)"))]
  }

  #"beta" is actually the risk difference if measure is "RD"
  if (measure=="RD"){
    names(res5star)[1] = names(res5starprelim)[1] = c("riskdiff")
  }
  if (measure=="RMST"){
    names(res5star)[1] = names(res5starprelim)[1] = c("rmstdiff")
  }
  if (measure=="MD"){
    names(res5star)[1] = names(res5starprelim)[1] = c("meandiff")
  }

  #============================================================================#

  #------------------------------------------#
  # making forest plots to summarize results #
  #------------------------------------------#

  if (plot){

    #calc rmst within each strata also to include in forest plots
    if (family=="cox" & measure == "HR"){

      #fit model averaged AFT fit within each stratum to add in forest plot

      #-------------------------------------------------#
      # descriptive fits - model averaged aft by strata #
      #-------------------------------------------------#
      descaft = maaftbystrata(time=time,status=status,arm=arm,
                              treeStrata=treeStrata,
                              distList=distList,ucvar=ucvar,
                              alternative=ifelse(
                                 alternative=="greater","less","greater"),
                              cilevel=ifelse(
                                  alternative=="two.sided",cilevel/2,cilevel),
                              verbose=0)

      descfit = descaft$fitsummary

      descRes = minPadapHR(sf=descaft$stratafit,betas=descfit[,"bhat"],
                           vars=descfit[,"v(bhat)"],
                           cilevel=ifelse(alternative=="two.sided",cilevel/2,
                                          cilevel),
                           alternative=ifelse(alternative=="greater",
                                               "less","greater"),
                           vartype="alt")$adaptivewtRes

      pruneddescaft = maaftbystrata(time=time,status=status,arm=arm,
                                    treeStrata=treeStrataPruned,
                              distList=distList,ucvar=ucvar,
                              alternative=ifelse(
                                 alternative=="greater","less","greater"),
                              cilevel=ifelse(alternative=="two.sided",cilevel/2,
                                             cilevel),
                              verbose=0)

      pruneddescfit = pruneddescaft$fitsummary

      pruneddescRes = minPadapHR(sf=pruneddescaft$stratafit,
                                 betas=pruneddescfit[,"bhat"],
                                 vars=pruneddescfit[,"v(bhat)"],
                                 cilevel=ifelse(alternative=="two.sided",cilevel/2,
                                                cilevel),
                                 alternative=ifelse(alternative=="greater",
                                                     "less","greater"),
                                 vartype="alt")$adaptivewtRes


    } else if (family=="cox" & measure=="TR"){

      #-------------------------------------------#
      # descriptive fits - cox ph model by strata #
      #-------------------------------------------#
      desccox = coxbystrata(time,status,arm,treeStrata,termNodes,
                            treetype="preliminary",
                            alternative=ifelse(
                              alternative=="less","greater","less"),
                            cilevel=ifelse(
                               alternative=="two.sided",cilevel/2,cilevel),
                            inclfrailty=inclfrailty,verbose=0,plot=FALSE,
                            timeunit,shading)
      descfit = desccox$fitsummary

      descRes = minPadapHR(sf=desccox$stratafit,betas=descfit[,"bhat"],
                           vars=descfit[,"v(bhat)"],
                           cilevel=ifelse(alternative=="two.sided",cilevel/2,
                                          cilevel),
                           alternative=ifelse(alternative=="less",
                                              "greater","less"),
                           vartype="alt")$adaptivewtRes

      pruneddesccox = coxbystrata(time=time,status=status,arm=arm,
                                  treeStrata=treeStrataPruned,
                                  termNodes=prunedTermNodes,
                                  treetype="final",alternative=ifelse(
                                    alternative=="less","greater","less"),
                                  cilevel=ifelse(
                                     alternative=="two.sided",cilevel/2,cilevel),
                                  inclfrailty=inclfrailty,
                                  verbose=0,plot=FALSE,timeunit,shading)

      pruneddescfit = pruneddesccox$fitsummary

      pruneddescRes = minPadapHR(sf=pruneddesccox$stratafit,
                                 betas=pruneddescfit[,"bhat"],
                                 vars=pruneddescfit[,"v(bhat)"],
                                 cilevel=ifelse(alternative=="two.sided",cilevel/2,
                                                cilevel),
                                 alternative=ifelse(alternative=="less",
                                                    "greater","less"),
                                 vartype="alt")$adaptivewtRes

    } else descfit = pruneddescfit = descRes = prunedDescRes = NULL

    #----------------------------------------------#
    # Forest plots of preliminar and final results #
    #----------------------------------------------#

    #forest plot of preliminary strata
    fplot = strataforestplot(MRSSmat = MRSSmat,MRSSres=res5starprelim,
                             stratafits = ByStrataFits$stratafit,family=family,
                             measure=measure,treetype="preliminary",
                             labelNames=termNodes,cilevel=cilevel,
                             alternative=alternative,descfit=descfit,
                             descRes=descRes,fplottype=fplottype)

    #forest plot of final strata
    fplotprune = strataforestplot(MRSSmat = prunedMRSSmat,MRSSres=res5star,
                                  stratafits = prunedByStrataFits$stratafit,
                                family=family,measure=measure,treetype="final",
                                labelNames=prunedTermNodes,cilevel=cilevel,
                                alternative=alternative,descfit=pruneddescfit,
                                descRes=pruneddescRes,fplottype=fplottype)

  } else fplot = fplotprune = NULL

  #============================================================================#

  #----------------------------------------------------#
  # Returning final list of all 5-STAR function output #
  #----------------------------------------------------#
  return(list(inputcovs=inputcovs,filteredcovs=cov2keep.all,
              filterdetails=filterdetails,
              prelimtree=stree,prelimstratadefn=termNodes,
              prelimstrataids=treeStrata,
              prelimstratasurvfits=ByStrataFits$stratafit,
              prelimbystratasummary=MRSSmat,
              prelimsinglewtres=singlewtresprelim,
              prelimres5star=res5starprelim,
              finaltree=streeprune,stratadefn=prunedTermNodes,
              strataids=treeStrataPruned,
              stratasurvfits=prunedByStrataFits$stratafit,
              bystratasummary=prunedMRSSmat,
              singlewtres=singlewtres,res5star=res5star,
              prelimbystrataPlot=p3,prelimbetweenstrataPlot=p4,prelimforestplot=fplot,
              bystrataPlot=p3prune,betweenstrataPlot=p4prune,forestplot=fplotprune))

} #end of main function
