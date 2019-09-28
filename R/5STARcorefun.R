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
#' Forms dummy variable matrix for factors, and fits an elastic net or random
#' forest model to determine which covariates to keep for building trees in 5-STAR
#'
#' @param yy      Response - either Surv() object for time to event data or 1
#' column matrix of case/control status for binary data
#' @param X       Data frame of all possible stratification covariates
#' @param family  Trait family, current options: "cox", "binomial", or "gaussian"
#' @param cilevel Optional significance level for RF VIMP confidence intervals
#' @param plot    Whether to make include VIMP CI plots (ignored if
#' method = "ENET")
#' @param verbose Numeric variable indicating amount of information to print
#' to the terminal (0 = nothing, 1 = notes only, 2 = notes and intermediate output)
#' @param filter_control List of contorl parameters for filtering step, @seealso
#' \code{\link{filter_control}}; key agruments include: \itemize{
#' \item method: filtering method; current options are: "ENET" (Elastic Net),
#' "RF" (Random Forest with default or input parameters, using delete-d
#' jackknife confidence intervals for VIMP for variable selection), "RFbest"
#' (Random Survival Forest, performing tuning for mtry and nodesize, and using
#'  double bootstrap confidence intervals for VIMP for variable selection;
#'  more accurate but slower than RF option)
#'  \item mixparm: Optional elastic net mixing parameter alpha or grid of alpha
#' values to search over. If nothing is entered, will search for best value
#' between 0.05 and 0.95. Ignored when method is "RF" or "RFbest"
#'  \item lambdatype: Optional elastic net parameter; whether to use lambda that
#' minimizes cross validation error (lambdatype="min") or the largest lambda
#' that gives error within 1 standard error of minimum error (lambdatype="1se",
#' default). Ignored when method is "RF" or "RFbest"
#'  \item ... : Optional arguments for glmnet or rfsrc
#' }
#'
#' @return cov2keep: List of all covariates to keep after filtering step, selected by
#' elastic net or random forest
#' @return For method=ENET, additionally outputs \itemize{
#'    \item cvout: containing fraction of null deviance explained, mean cross validation
#' error, and optimal lambda and alpha value (see cv.glmnet and glmnet for
#' more details)
#'    \item beta: the coefficients of glmnet fit given tuned choice of alpha and lambda
#' }
#' @return For method=RF or RFbest, additionally outputs:
#' \itemize{
#'    \item varselectmat:  matrix of VIMP confidence intervals, p-values, and selection
#' decision for each variable
#'    \item VIMPplot: default variable importance CI plot, output if plot = TRUE
#'    \item varselect2plot: variable selection information from subsample.rfsrc for
#' making customizable VIMP CI plots
#'    \item forest: rfsrc object
#' }
#' @export
filter5STAR = function(yy,X,family="cox",cilevel=0.05,plot=FALSE,verbose=0,
                       filter.hyper=filter_control(
                         method="ENET",lambdatype="min",mixparm=NULL,...)){

  method = filter.hyper$method
  lambdatype = filter.hyper$lambdatype
  mixparm = filter.hyper$mixparm

  optargs = filter.hyper[!(names(filter.hyper) %in%
                             c("method","lambdatype","mixparm"))]

  X = droplevels(X)

  if (method == "ENET"){

    if (!(lambdatype %in% c("1se","min"))) stop("lambdatype must be one of
                                                '1se',`min'.")

    #ensuring yy is a column
    if (is.null(nrow(yy))) yy = as.matrix(yy)

    #for now, complete case analysis for filtering
    #(i.e., removing individuals with any missing covariate information)
    XC.nomiss  = X[rowSums(is.na(X))==0,]
    yyC.nomiss = yy[rowSums(is.na(X))==0,]

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

    polymorphicDummy = sapply(XCd.dummy,function(x) length(levels(x)))>1
    monoDummyNames = dummyNames[!polymorphicDummy]
    if (length(monoDummyNames >0) & verbose > 0) print(paste0(
      "Note: removing covariates ",paste0(monoDummyNames,collapse=", "),
      " from analysis (only 1 level after removing missing subjects)"))
    dummyNames = dummyNames[!(dummyNames %in% monoDummyNames)]
    dummyNames = paste0(dummyNames,".d")
    XCd.dummy = XCd.dummy[,sapply(XCd.dummy,function(x) length(levels(x)))>1]

    if (is.null(nrow(XCd.dummy))) XCd.dummy = data.frame(XCd.dummy)
    colnames(XCd.dummy) = dummyNames

    if (ncol(X)!= ncol(XC.numeric)+ncol(XCd.dummy)+length(monoDummyNames))
      warning("No. columns in numeric and factor matrices not equal to number of covariates")

    if (length(dummyNames)>0){
      dummyFormula = paste(dummyNames,collapse=" + ")
      XC.dummy = model.matrix(as.formula(paste0("~ ",dummyFormula)), XCd.dummy)
      XC.enet = data.matrix(cbind(XC.dummy,XC.numeric)[,-1])
    } else XC.enet = data.matrix(XC.numeric)

    if (family == "cox"){
      foldid=quiet(c060::balancedFolds(class.column.factor=yyC.nomiss[,2],
                                       cross.outer=10))
    } else if (family=="binomial") {
      foldid=quiet(c060::balancedFolds(class.column.factor=yyC.nomiss,
                                       cross.outer=10))
    } else if (family == "gaussian"){
      foldid = sample(1:10,size=length(yyC.nomiss),replace=TRUE)
    }

    #fitting survival elastic net searching for best alpha
    if (is.null(mixparm)){
      a <- seq(0.05, 0.95, 0.05)
    } else a = mixparm

    search = NULL
    for(i in a) {
      cvfit = do.call(glmnet::cv.glmnet,c(list(
        XC.enet, yyC.nomiss,family = family,alpha = i, foldid=foldid),optargs))
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
      alpha = cv$alpha),optargs))

    cvoutput = cbind(md$dev.ratio,cv)
    colnames(cvoutput) = c("dev.ratio","cvm","lambda","alpha")

#     library(magrittr)
      enetcoef = coef(md)
#     #standardized coefficient
#     enetcoef.std = enetcoef*apply(XC.enet,2,sd)
#     df = data.frame(XC.enet)
#     sds <- df %>% sapply(sd)
#     enetcoef.std2 <- enetcoef * c(1, c(sds[names(enetcoef)][-1]))
#     std_coeff_lr2[1] <- coeff_lr2[1] + sum(coeff_lr2[-1] * mus[names(coeff_lr2)][-1])
    # XC.enetscale = scale(XC.enet)

    # #plots of lambda solution path and/or coefficients
    # if (plot){
    #   enetplotmat = data.frame(covariate=rownames(enetcoef.std),beta=enetcoef.std[,"s0"])
    #   enetplotmat$covariate = gsub(".d",".",enetplotmat$covariate)
    #   enetplotmat = enetplotmat[order(abs(enetplotmat[,2]),decreasing=FALSE),]
    #   pvimp = ggplot2::ggplot(enetplotmat[enetplotmat[,2]!=0,]) +
    #     ggplot2::geom_col(ggplot2::aes(y=abs(beta),x=covariate),
    #                       col="black",fill="forestgreen") +
    #     ggplot2::ylab(expression(abs(beta))) + ggplot2::coord_flip() +
    #     ggplot2::theme_bw() +ggplot2::ggtitle("ENET Important Variables") +
    #     ggplot2::xlab("") +ggplot2::scale_x_discrete(limits=enetplotmat[
    #       enetplotmat[,2]!=0,]$covariate)
    # } else pvimp = NULL

    #extracting covariates to keep from filtering step
    cov2keep = names(enetcoef[,1])[abs(enetcoef[,1])>0]
    if (length(dummyNames)>0){
      cov2keep.cat = dummyNames[sapply(dummyNames,function(x)
        length(grep(x,cov2keep))>0)]
      cov2keep.cat = sub(".d","",cov2keep.cat)
    } else cov2keep.cat = NULL
    cov2keep.num = cov2keep[cov2keep %in% names(XC.numeric)]
    cov2keep.all = c(cov2keep.cat,cov2keep.num)

  } else if (method == "RF"){

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
    srv.varselect.o=randomForestSRC::extract.subsample(srv.smp.o,alpha=cilevel)$
      var.jk.sel.Z
    cov2keep.all = rownames(srv.varselect.o)[which(srv.varselect.o$signif)]

  } else if (method == "RFbest"){

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
      srv.smp.o,alpha=cilevel)$var.sel.Z
    cov2keep.all = rownames(srv.varselect.o)[which(srv.varselect.o$signif)]

  } else{
    warning("Invalid filter type. Please select one of 'ENET','RF', or
                 'RFbest'. Returning all covariates.")
    cov2keep.all = colnames(X)
  }

  covOrder=match(colnames(X),cov2keep.all)[!is.na(match(colnames(X),cov2keep.all))]
  cov2keep.all = cov2keep.all[covOrder]

  if (method=="ENET"){
    return(list(cov2keep=cov2keep.all,cvout=cvoutput,beta=enetcoef))#,
                #VIMPplot=pvimp))
  } else if (method == "RF"|method=="RFbest"){
    return(list(cov2keep=cov2keep.all,varselectmat=srv.varselect.o,
                  VIMPplot=pvimp,varselect2plot=srv.smp.o,forest=srv.o))
  }
}

################################################################################

simplifyand = function(node){

  #sepand = do.call(rbind,strsplit(node,"&"))
  septab = do.call(rbind,sapply(trimws(strsplit(node,"&")[[1]]),function(x)
    strsplit(x," |(?>\\(.*?\\).*?\\K(, |$))",perl=TRUE)))

  #septab = do.call(rbind,strsplit(node," "))

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
#' Fits preliminary (3A) and pruned (3B) trees to determine homogenous risk
#' strata for use with 5-STAR algorithm
#'
#' @param yy      Trait/response (currently must be either a binary or continuous
#'  covariate or a Surv() object summarizing follow-up time for right-censored
#'  data and status indicator where 1=dead, 0=censored)
# @param arm     Treatment indicator, 1 = test treatment, 0 = control
#' @param X       Data frame of all possible stratification covariates
#' @param family  Trait family, current options: "cox", "binomial", or "gaussian"
#' @param cilevel Confidence level alpha for overall result and confidence
#' intervals
#' @param verbose Numeric variable indicating amount of information to print
#' to the terminal (0 = nothing, 1 = notes only, 2 = notes and intermediate output)
#' @param tree.hyper List of control variables for tree fitting, @seealso
#' \code{\link{tree_control}}, many of which will be passed into
#'  \code{\link[partykit]{ctree_control}}
#'
#' @return \itemize{
#'     \item strataids: vector of subject-level strata membership
#'     \item stratadef: definition of final formed strata, in terms of covariates
#'     \item finaltree: final, pruned tree built by ctree, in terms of
#'      preliminary strata membership
#'     \item prelimstrataids: vector of subject-level strata membership for
#'      preliminary, unpruned tree
#'     \item prelimstratadef: definition of final formed strata, in terms of
#'     covariates for preliminary, unpruned tree
#'     \item prelimtree: preliminary, unpruned tree built by ctree, in terms
#'      of covariates
#' }
#'
#' @import partykit
#' @export
#'
fittrees = function(yy,X,family="cox",cilevel=0.05,verbose=0,
                     tree.hyper = tree_control()){
                       # minbucket=40,alpha=c(0.1,0.2),testtype="Bonferroni",
                       # majority=FALSE,maxsurrogate=3,...)){

  #----------------------------------------------------#
  # Step 3: Conditional Inference Trees to Find Strata #
  #----------------------------------------------------#

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
                              c("minbucket","alpha"))] #"splitweights",
  #control = do.call(ctree_control,ctreeparms)

  #---------#
  # step 3A #
  #---------#

  #fitting ctree if there are any variables left after filtering
  if (length(colnames(X))>0){

    #if modifying list, need to change logmincriterion as alpha is not retained
    #from ctree_control() function
    ###control$minbucket = minbucket1
    ###control$logmincriterion = log(1-alpha1)

    #running ctree algorithm
    ###stree = partykit::ctree(yy ~ ., data=datdf,control=control)
    stree = do.call(partykit::ctree,c(list(yy ~ ., data=datdf,
                                           minbucket=minbucket1,alpha=alpha1),
                                      ctreeparms))

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

      # #writing terminal nodes in a way to extract subjects fitting those criteria
      # termNodesSep = strsplit(termNodes,"&")
      # dataName = "datdf"
      # termNodesSep = lapply(termNodesSep,function(x) paste0(dataName,"$",trimws(x)))
      # termNodesFinal = lapply(termNodesSep,function(x) paste(x,collapse=" & "))
      #
      # #list of indicator of all subjects in each stratum/node
      # inNode = lapply(1:length(termNodes),function(x){
      #   eval(parse(text=termNodesFinal[x]))
      # } )
      # inNodeMat = do.call(cbind,inNode)

      #determining which preliminary strata each subject belongs in
      #(this does not lead to issues with NAs like above approach, as
      #even subjects with missing values will be assigned a node)
      pS = predict(stree,type="node")
      datdf2 = datdf
      if (length(pS) < nrow(datdf)){
        datdf2 = datdf[(rownames(datdf) %in% names(pS)),]
      }

      #ordering preliminary strata to input into step 3B
      if (family == "cox"){

        #calculating survival km within each strata
        #inputting preliminary strata into tuning step as ordered factor
        sf=survfit(yy~pS,data=datdf2)
        rmeans = summary(sf,rmean="common")$table[,"*rmean"]
        rmeanriskorder = order(rmeans)
        orderedstratalevels = sort(unique(pS))[rmeanriskorder]
        pS = ordered(pS,levels=orderedstratalevels)
        pS = ordered(as.numeric(pS))

      } else if (family == "binomial"){

        rmeans = unique(predict(stree,type="prob",simplify=TRUE))[,1]
        rmeanriskorder = order(rmeans)
        orderedstratalevels = unique(pS)[rmeanriskorder]
        pS = ordered(pS,levels=orderedstratalevels)

      } else if (family == "gaussian"){

        rmeans = unique(predict(stree,type="response",simplify=TRUE))
        rmeanriskorder = order(rmeans)
        orderedstratalevels = unique(pS)[rmeanriskorder]
        pS = ordered(pS,levels=orderedstratalevels)

      }

      #reordering terminal nodes so output is given in terms of risk order (high -> low)
      termNodes = termNodes[as.character(orderedstratalevels)]

      #---------#
      # step 3B #
      #---------#

      #post-pruning strata: fit ctree again w/ ordered strata as input variables
      ###control$minbucket = minbucket2
      ###control$logmincriterion = log(1-alpha2)
      ###streeprune = partykit::ctree(yy ~ pS,data=datdf2,control=control)
      # streeprune = do.call(partykit::ctree(yy ~ pS,data=datdf2),
      #                      c(list(minbucket=minbucket2,alpha=alpha2),ctreeparms))
      streeprune = do.call(partykit::ctree,c(list(yy ~ pS, data=datdf2,
                                             minbucket=minbucket2,alpha=alpha2),
                                        ctreeparms))
      termNodes2 = .list.rules.party(streeprune)

      #continuing if any strata were found in step 3B
      if (length(termNodes2)>1){

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

        prunedTermNodes = unlist(lapply(newNodeVars,function(x)
          paste(paste0("(",termNodes[x],")"),collapse = " | ")))
        prunedTermNodesFinal = gsub("(^c*)\\(","\\(datdf$",prunedTermNodes)
        prunedTermNodesFinal = gsub("& ","& datdf$",prunedTermNodesFinal)
        prunedTermNodesFinal = gsub("\\| \\(","\\| \\(datdf$",
                                    prunedTermNodesFinal)

        # #list of indicator of all subjects in each stratum/node
        # inNodePruned = lapply(1:length(prunedTermNodes),function(x){
        #   eval(parse(text=prunedTermNodesFinal[x]))
        # } )
        # inNodeMatPruned = do.call(cbind,inNodePruned)

        #extracting final node for each subject
        prunedpS = predict(streeprune,type="node")
        datdf3 = datdf2
        if (length(prunedpS) < nrow(datdf2)){
          datdf3 = datdf2[(rownames(datdf2) %in% names(prunedpS)),]
        }

      } else { # case where only one strata after pruning

        prunedTermNodes = prunedTermNodesFinal = "datdf"
        #inNodeMatPruned = matrix(1,nrow=nsubj,ncol=1)

      }

    } else{ #case where no initial subgroups are detected
      streeprune = NULL
      termNodes = termNodesFinal = "datdf"
      #inNodeMat = inNodeMatPruned = matrix(1,nrow=nsubj,ncol=1)
      prunedTermNodes = prunedTermNodesFinal = "datdf"
      pS = as.factor(rep(1,nsubj))
      prunedpS = rep(1,nsubj)

    }

  } else{ #case where no subgroups are detected or no tree was fit due to
          #nothing passing filtering step

    termNodes = termNodesFinal = "datdf"
    #inNodeMat = matrix(1,nrow=nsubj,ncol=1)
    prunedTermNodes = prunedTermNodesFinal = "datdf"
    #inNodeMatPruned = matrix(1,nrow=nsubj,ncol=1)
    pS = as.factor(rep(1,nsubj))
    prunedpS = rep(1,nsubj)

  }

  #============================================================================#

  #--------------------------------------------------------------------------#
  # removing observations, if any, that don't fit into any strata defintions #
  #--------------------------------------------------------------------------#

  if (!is.data.frame(datdf)) datdf = as.data.frame(datdf)
  #datPruned = datdf; armPruned=arm; timePruned = time; statusPruned = status;

  #extracting final strata (as list of strata each subj belongs to)
  treeStrata = as.numeric(factor(pS,ordered=FALSE,
                                 levels=sort(as.numeric(levels(pS)))))
  treeStrataPruned = as.numeric(factor(prunedpS,ordered=FALSE))

  #number of found strata
  nstrata = length(unique(treeStrata))
  nstrataPruned = length(unique(treeStrataPruned))

  #============================================================================#

  #-----------------------------------------------------------------------#
  # ensuring tree strata/terminal nodes are ordered high risk -- low risk #
  #-----------------------------------------------------------------------#

  if (family == "cox"){
    #order preliminary strata by risk determined by restricted mean (AUC of KM curve)
    sf=survfit(yy~treeStrata)
    if (length(unique(treeStrata))>1){
      rmeans = summary(sf,rmean="common")$table[,"*rmean"]
    } else rmeans = summary(sf,rmean="common")$table["*rmean"]
    riskOrder = order(rmeans)

    #order pruned strata by risk determined by restricted mean (AUC of KM curve)
    prunedsf=survfit(yy~treeStrataPruned)
    if (length(unique(treeStrataPruned))>1){
      prunedrmeans = summary(prunedsf,rmean="common")$table[,"*rmean"]
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
#' (with possible modifications through ...), or "RFbest" for tuned random forest
#' @param lambdatype Optional elastic net parameter; whether to use lambda that
#' minimizes cross validation error (lambdatype="min", default) or the largest
#' lambda that gives error within 1 standard error of minimum error
#' (lambdatype="1se"). Ignored when method is "RF" or "RFbest"
#' @param mixparm Optional elastic net mixing parameter alpha or grid of alpha
#' values to search over. If nothing is entered, will search for best value
#' between 0.05 and 0.95. Ignored when method is "RF" or "RFbest"
#' @param ... Optional additional arguments passed into
#' \code{\link[glmnet]{glmnet}} or \code{\link[randomForestSRC]{rfsrc}} functions
#' @return A list of control parameters for filtering step
#' @export
filter_control = function(method="ENET",lambdatype="min",mixparm=NULL,...){
  list(method=method, lambdatype=lambdatype, mixparm=mixparm,...)
}


#' Control for strata-finding tree step of 5-STAR algorithm
#'
#' Parameters for control of conditional inference tree steps of 5-STAR algorithm
#' (3A and 3B) - mainly for passing into ctree and ctree_control functions
#' @param minbucket Vector of minimum set of weights per terminal node for
#' initial and pruning steps - if a single number is given, the same value is
#' used for both preliminary and final trees
#' @param alpha vector of significance level for variable selection for tree
#'  splits in preliminary and final trees (3A and 3B) - if a single number is
#'  given, the same value is used for both preliminary and final trees
#' @param testtype from \code{\link[partykit]{ctree_control}}:
#' "a character specifying how to compute the distribution of the test
#' statistic" (default = "Bonferroni")
#' @param majority from \code{\link[partykit]{ctree_control}}:
#' whether to place all missing observations with the most common node, or
#' assign them randomly
#' @param maxsurrogate from \code{\link[partykit]{ctree_control}}: "number of
#' surrogate splits to evaluate"
#' @param testtype  ctree control parameter, specifying how to compute the
# distribution of the test statistic (see ctree_control)
#' @param majority  ctree control parameter, specifying whether to randomly
# assign subjects with missing information in splitting variable (FALSE) or to
# go with the majority (TRUE) (see ctree_control)
#' @param maxsurrogate ctree control parameter defining number of surrogate
# splits to evaluate for missing covariates (see ctree_control)
#' @param ... additional parameters to be passed into
#' \code{\link[partykit]{ctree}} function
#' @return A list of control parameters for strata formation step
#' @export
tree_control = function(minbucket=40,alpha=c(0.1,0.2),testtype="Bonferroni",
                        majority=FALSE,maxsurrogate=3,maxdepth=3,...){
  list(minbucket=minbucket,alpha=alpha,testtype=testtype,
       majority=majority,maxsurrogate=maxsurrogate,maxdepth=maxdepth,...)
}

################################################################################

#' Prepares data for 5-STAR algorithm
#'
#' Cleans covariate matrix, removing covariates with too much missingess,
#' and coverts all character covariates into factors
#' @param yy Trait/response (currently must be either continous, binary, or
#' a Surv() object summarizing follow-up time for right-censored data and status
#'  indicator where 1=dead, 0=censored). Ignored for family != "cox"
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
#' @param verbose Numeric variable indicating amount of information to print
#' to the terminal (0 = nothing, 1 = notes only, 2 = notes and intermediate output)
#' @return \itemize{
#' \item X: cleaned covariate matrix
#' \item input covs: names of all covariates which will be passed into filtering
#' step of 5-STAR algorithm
#' }
#' @export
prep5STAR = function(yy,X,family="cox",missthreshold=c(0.1,0.2),verbose=0){

  #------------------------------------------#
  # converting character variables to factor #
  #------------------------------------------#

  ischar = sapply(X,is.character)
  if (sum(ischar)>0){
    if (verbose > 0) print(paste0("Note: character variables ",
                 paste(names(XX[ischar]),collapse=", "),
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

        #removing those not significantly correlated
        torm = names(cortest.all)[cortest.all >= 0.05]

       #(only one covariate with missingness between boundaries)
      } else {

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

        # tokeep = ifelse(cor.test(LRscores,Xmiss,method="kendall",
        #                          alternative="two.sided")$p.value < 0.05,
        #                 names(which(colMeans(is.na(XX)) > 0.1 &
        #                               colMeans(is.na(XX)) < 0.2)),NA)
      }

      if (length(torm[!is.na(torm)]) > 0) X = X[,!(colnames(X) %in% torm)]
      #XX = XX[,(colMeans(is.na(XX)) <= missThreshold)|(colnames(XX) %in% tokeep)]
    } else torm = NA #nothing with > missthreshold[1] missingness -- nothing to do here
    ##else XX = XX[,colMeans(is.na(XX)) <= missThreshold]

  } else torm = NA

  #printing out which covariates were removed
  covs2rm.all = c(names(covs2rm),torm[!is.na(torm)])

  if ((length(covs2rm.all) > 0) & (verbose > 0)){
      print(paste0("Removing covariates ",paste(covs2rm.all,collapse=", "),
                   " from analysis due to too much missingness"))
  }

  # #covariates to remove due to too much missingness
  # covs2rm = which(colMeans(is.na(X)) > missthreshold)
  # if (length(covs2rm)>0){
  #   print(paste0("Removing covariates ",paste(names(covs2rm),collapse=", "),
  #                " from analysis due to too much missingness (>",
  #                missthreshold,")"))
  #   X = X[,-covs2rm]
  # }

  return(X)
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
#' @param X       Data frame of all possible stratification covariates
#' @param family  Trait family, current options: "cox", "binomial", and
#' "gaussian" (note: binomial and gaussian are still experimental!)
#' @param measure Response of interest for binary traits; current options
#' are: "RD" (rate difference using Miettinen and Nurminen 1985 method), or
#' "OR" (odds ratio, using GLM model fit); ignored for family!="binomial"
#' @param cilevel Confidence level alpha for overall result and confidence
#' intervals
#' @param missthrehold Optional parameter specifying proportion of missingness
#' allowed for each covariate. See \code{\link{prep5STAR}} for more details.
#' @param verbose   Numeric variable indicating amount of information to print
#' to the terminal (0 = nothing, 1 = notes only, 2 = notes and intermediate output)
#' @param plot       Logical, whether to return summary plots (within strata and
#' pooled KM plots, forest plots summarizing strata-level and overall results)
#' @param timeunit  Optional parameter for plots defining scale of time to event
#'  data (e.g., "Months", "Years"); ignored for family!="cox"
#' @param filter.hyper A list of control parameters for filtering step,
#' (see \code{\link{filter_control}})
#' @seealso \code{\link{filter_control}}
#' @param tree.hyper A list of control parameters for strata building step
#' (see \code{\link{tree_control}})
#' @seealso \code{\link{tree_control}}
#'
#' @details Filtering step is performed on complete case basis (e.g., removing
#' all individuals with any missing covariate data)
#' @return \itemize{
#'     \item inputcovs: vector containing names of all covariates input into
#'     the algorithm (after removing those with too much missingness, as defined
#'     by missthreshold)
#'     \item filteredcovs: vector containing all covariates left after
#'     filtering step
#'     \item filterdetails: list containing additional output from filter step;
#'     for method=ENET, this includes: \itemize{
#'    \item cvout: containing fraction of null deviance explained, mean cross validation
#' error, and optimal lambda and alpha value (see cv.glmnet and glmnet for
#' more details)
#'    \item beta: the coefficients of glmnet fit given tuned choice of alpha and lambda
#' } For method=RF or RFbest, this includes:
#' \itemize{
#'    \item varselectmat: matrix of VIMP confidence intervals, p-values, and selection
#' decision for each variable
#'    \item VIMPplot: default variable importance CI plot, output if plot = TRUE
#'    \item varselect2plot: variable selection information from subsample.rfsrc for
#' making customizable VIMP CI plots
#'    \item forest: rfsrc object
#' }
#'     \item prelimtree: preliminary tree (3A) built by ctree, in terms
#'      of covariates
#'     \item prelimstratadefn: vector of rules defining the preliminary strata
#'     (step 3A), as defined by covariates chosen in initial ctree run, ordered from
#'     highest to lowest risk
#'     \item prelimstrataids: vector of preliminary (3A) strata membership for
#'     each subject (where 1 indicates highest risk, and the largest number
#'     indicates lowest risk, i.e., matching order of prelimstratadefn)
#'     \item prelimstratasurvfits: a \code{\link[survival]{survfit}} object,
#'     summarizing survival curves in each preliminary (3A) strata, pooled by
#'     treatment assignment; returned when family="cox"
#'     \item prelimbystratasummary: matrix summarizing  cox/glm fits within
#'     each preliminary (3A) strata, including estimated (log) hazard ratio
#'     (or OR, RD, etc. when family != "cox"), corresponding (1-cilevel)x100\%
#'     CI, Gramsch and Therneau test of non-PH p-value, and sample size,
#'     inverse variance, and adaptive weights
#'     \item prelimadaptivecorr: rank correlation between estimated coefficients
#'     beta and corresponding variances from each preliminary (3A) strata
#'     \item prelimres5star: Summary of amalgamated hazard ratio (family="cox") or
#'     odds ratio (family="binomial") estimate, variance,
#'     p-value, and (1-cilevel)x100\% confidence interval for adaptive weights
#'     from preliminary strata (3A)
#'     \item finaltree: final, pruned tree (3B) built by ctree, in terms of
#'      preliminary strata
#'     \item stratadefn: vector of rules defining the final strata
#'     (step 3B), as defined by covariates translated from pruning ctree run,
#'     ordered from highest to lowest risk
#'     \item strataids:  vector of final (3B) strata membership for
#'     each subject (where 1 indicates highest risk, and the largest number
#'     indicates lowest risk, i.e., matching order of stratadefn)
#'     \item stratasurvfits: a \code{\link[survival]{survfit}} object,
#'     summarizing survival curves in each final (3B) strata, pooled by
#'     treatment assignment; returned when family="cox"
#'     \item bystratasummary: matrix summarizing  cox/glm fits within
#'     each final (3B) strata, including estimated (log) hazard ratio (or OR,
#'     RD, etc. when family != "cox"), corresponding (1-cilevel)x100\% CI,
#'     Gramsch and Therneau test of non-PH p-value, and sample size,
#'     inverse variance, and adaptive weights
#'     \item adaptivecorr: rank correlation between estimated coefficients
#'     beta and corresponding variances from each final (3B) strata; used
#'     to determine if adaptive weights will combine betas or test statistics
#'     \item res5star: Summary of amalgamated hazard ratio (family="cox") or
#'     odds ratio (family="binomial") estimate, variance,
#'     p-value, and cilevelx100% confidence interval for adaptive weights
#'     \item prelimbetweenstrataKM: Kaplan-Meier survival curve plots from pooled
#'      treatment assignment data from each preliminary (3A) strata,
#'      returned if plot = TRUE and family = "cox"
#'     \item prelimbystrataKM: Kaplan-Meier survival curve plots within each
#'     preliminary (3A) strata, returned if plot = TRUE and family = "cox"
#'     \item prelimforestplot: a forest plot summarizing estimated coefficient and
#'     confidence interval for each final (3B) strata, as well as amalgamated
#'     result
#'     \item betweenstrataKM: Kaplan-Meier survival curve plots from pooled
#'      treatment assignment data from each final (3B) strata,
#'      returned if plot = TRUE and family = "cox"
#'     \item bystrataKM: Kaplan-Meier survival curve plots within each final
#'     (3B) strata, returned if plot = TRUE and family = "cox"
#'     \item forestplot: a forest plot summarizing estimated coefficient and
#'     confidence interval for each final (3B) strata, as well as amalgamated
#'     result
#' }
#'
#' @import partykit
#' @import survival
#' @export
run5STAR = function(yy,arm,X,family="cox",measure="OR",cilevel=0.05,
                    missthreshold=c(0.1,0.2),verbose=0,plot=TRUE,
                    timeunit=NULL,filter.hyper = filter_control(),
                    tree.hyper = tree_control()){


  #need to specify measure if family is binomial
  if (family=="binomial" & !(measure %in% c("RD","OR"))){
    stop("For binary traits, measure must be either `RD` or `OR`.")
  }

  #------------------------------------------------------#
  # extracting/comining dat of interest to usable format #
  #------------------------------------------------------#

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
                verbose=verbose)

  #input covariate names
  inputcovs = names(X)

  #============================================================================#

  #---------------------------#
  # Step 2: Filter Covariates #
  #---------------------------#

  filterRes = filter5STAR(yy=yy,X=X,family=family,cilevel=cilevel,plot=plot,
                          verbose=verbose,filter.hyper=filter.hyper)

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
  #VIMPplot = filterRes$VIMPplot
  filterdetails = filterRes[names(filterRes) != "cov2keep"]

  #============================================================================#

  #----------------------------------------------------#
  # Step 3: Conditional Inference Trees to Find Strata #
  #----------------------------------------------------#

  trees5STAR = fittrees(yy=yy,X=X,family=family,cilevel=cilevel,verbose=verbose,
                        tree.hyper=tree.hyper)

  treeStrataPruned = trees5STAR$strataids
  prunedTermNodes = trees5STAR$stratadefn
  streeprune = trees5STAR$finaltree
  treeStrata = trees5STAR$prelimstrataids
  termNodes = trees5STAR$prelimstratadefn
  stree = trees5STAR$prelimtree

  if (verbose > 1){
    print(paste0("Preliminary Strata: ",termNodes))
    print(paste0("Final Pruned Strata: ",prunedTermNodes))
  }

  #============================================================================#

  #-------------------------------------------------------------------------#
  # calculating number of events per strataxarm                             #
  # skipping the rest of simu loop if not at least one event per strataxarm #
  #-------------------------------------------------------------------------#

  statusVec = yy
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
                    sum(status[treeStrata==x])),collapse="-"),
                  " and per final strata: ",
                  paste0(sapply(1:length(unique(treeStrataPruned)),function(x)
                    sum(status[treeStrataPruned==x])),collapse="-"),")"))


    if (family == "binomial"){
      if (min(perStrataEventSummary[3,])==0)
        stop("Need at least one control per strata x arm")
    }

  } else {
    minEventsPerStratArm = NULL
    minEventsPerStratArmPruned = NULL
  }

  #============================================================================#

  #-------------------------------------------------#
  # Step 4: Fit Cox Model Within Each Formed Strata #
  #-------------------------------------------------#

  if (family == "cox"){

    time = yy[,1]; status = yy[,2]
    ByStrataFits = coxbystrata(time,status,arm,treeStrata,termNodes,
                               treetype="preliminary",cilevel,verbose,plot,
                               timeunit)

    MRSSmat = ByStrataFits$fitsummary
    p3 = ByStrataFits$bystrataKM
    p4 = ByStrataFits$betweenstrataKM

    prunedByStrataFits = coxbystrata(time=time,status=status,arm=arm,
                                     treeStrata=treeStrataPruned,
                                     termNodes=prunedTermNodes,
                                     treetype="final",cilevel=cilevel,verbose,
                                     plot,timeunit)

    prunedMRSSmat = prunedByStrataFits$fitsummary
    p3prune = prunedByStrataFits$bystrataKM
    p4prune = prunedByStrataFits$betweenstrataKM

  } else if (family=="gaussian"|(family=="binomial" & measure=="OR")){

    ByStrataFits = glmbystrata(yy,arm,family,treeStrata,termNodes,cilevel,
                               verbose)
    MRSSmat = ByStrataFits$fitsummary

    prunedByStrataFits = glmbystrata(yy=yy,arm=arm,family=family,
                                     treeStrata=treeStrataPruned,
                                     termNodes=prunedTermNodes,
                                     cilevel=cilevel,verbose)

    prunedMRSSmat = prunedByStrataFits$fitsummary

    p3 = p4 = p3prune = p4prune = fplot = fplotprune = NULL

  } else if (family == "binomial" & measure == "RD"){

    ByStrataFits = ratediffbystrata(yy=yy,arm=arm,family=family,
                                    treeStrata=treeStrata,
                                    treetype="preliminary",termNodes=termNodes,
                                    cilevel=cilevel,verbose,plot)
    MRSSmat = ByStrataFits$fitsummary

    prunedByStrataFits = ratediffbystrata(yy=yy,arm=arm,family=family,
                                          treeStrata=treeStrataPruned,
                                          # treetype="final",
                                          termNodes=prunedTermNodes,
                                          cilevel=cilevel,verbose,plot)

    prunedMRSSmat = prunedByStrataFits$fitsummary

    #currently no plots output for family != "cox"
    p3 = p4 = p3prune = p4prune = NULL

  }

  #============================================================================#

  #----------------------------------------------------------------------#
  # Step 5: Amalgamate Results, Using New Weight if cor(betahati,Vi) > 0 #
  #----------------------------------------------------------------------#

  adaptiveWtResSS = adaptivewtHR(weights=ByStrataFits$weights,
                                 betas=MRSSmat[,1],#[,"bhat"],
                                 vars=MRSSmat[,2],#[,"v(bhat)"],
                                 cilevel=cilevel)

  adaptiveWtResSSprune = adaptivewtHR(weights=prunedByStrataFits$weights,
                                      betas=prunedMRSSmat[,1],#[,"bhat"],
                                      vars=prunedMRSSmat[,2],#[,"v(bhat)"],
                                      cilevel=cilevel)
  adaptivecorr = adaptiveWtResSS$adaptivecorr
  adaptivecorrpruned = adaptiveWtResSSprune$adaptivecorr

  #combining new weights in coxmrss output
  MRSSmat = data.frame(MRSSmat,weight.adap=adaptiveWtResSS$adaptivewts)
  prunedMRSSmat = data.frame(prunedMRSSmat,weight.adap=adaptiveWtResSSprune$adaptivewts)

  #storing results from adaptive weight test
  res5starprelim = adaptiveWtResSS$adaptivewtRes
  res5star = adaptiveWtResSSprune$adaptivewtRes

  #don't need exp results if trait is quantitative or measure is RD
  if (family == "gaussian"|(family=="binomial" & measure=="RD")){
    res5starprelim = res5starprelim[!(names(res5starprelim) %in%
                                        c("exp(beta)", "exp(ci lower)",
                                          "exp(ci upper)"))]
    res5star = res5star[!(names(res5star) %in%
                            c("exp(beta)", "exp(ci lower)","exp(ci upper)"))]
  }

  #"beta" is actually the rate difference if measure is "RD"
  if (measure=="RD"){
    names(res5star)[1] = names(res5starprelim)[1] = c("ratediff")
  }

  #============================================================================#

  #------------------------------------------#
  # making forest plots to summarize results #
  #------------------------------------------#

  if (plot){

    #forest plot of preliminary strata
    fplot = strataforestplot(MRSSmat = MRSSmat,MRSSres=res5starprelim,
                             stratafits = ByStrataFits$stratafit,family=family,
                             measure=measure,treetype="preliminary",
                             labelNames=termNodes,cilevel=cilevel)

    #forest plot of final strata
    fplotprune = strataforestplot(MRSSmat = prunedMRSSmat,MRSSres=res5star,
                                  stratafits = prunedByStrataFits$stratafit,
                                family=family,measure=measure,treetype="final",
                                labelNames=prunedTermNodes,cilevel=cilevel)

  } else fplot = fplotprune = NULL

  #============================================================================#

  return(list(inputcovs=inputcovs,filteredcovs=cov2keep.all,
              filterdetails=filterdetails,#VIMPplot=VIMPplot,
              prelimtree=stree,prelimstratadefn=termNodes,
              prelimstrataids=treeStrata,
              prelimstratasurvfits=ByStrataFits$stratafit,
              #prelimminevents=minEventsPerStratArm,
              prelimbystratasummary=MRSSmat,prelimadaptivecorr=adaptivecorr,
              prelimres5star=res5starprelim,
              finaltree=streeprune,stratadefn=prunedTermNodes,
              strataids=treeStrataPruned,
              stratasurvfits=prunedByStrataFits$stratafit,
              #minevents=minEventsPerStratArmPruned,
              bystratasummary=prunedMRSSmat,adaptivecorr=adaptivecorrpruned,
              res5star=res5star,
              prelimbystrataKM=p3,prelimbetweenstrataKM=p4,prelimforestplot=fplot,
              bystrataKM=p3prune,betweenstrataKM=p4prune,forestplot=fplotprune))

} #end of main function
