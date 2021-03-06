% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/5STARcorefun.R
\name{filter5STAR}
\alias{filter5STAR}
\title{Filter covariates for 5STAR}
\usage{
filter5STAR(yy, X, family = "cox", plot = FALSE, verbose = 0,
  filter.hyper = filter_control(method = "ENET", lambdatype = "min",
  mixparm = NULL, vimpalpha = 0.05, nfolds = 10, filterseed = 2019, ...),
  vars2keep = NULL)
}
\arguments{
\item{yy}{Response - either Surv() object for time to event data or 1
column matrix of case/control status for binary data}

\item{X}{Data frame of all possible stratification covariates}

\item{family}{Trait family, current options: "cox", "binomial", or "gaussian"}

\item{plot}{Whether to make plots for filter results (e.g., variable importance
plot for RF-based filtering or solution path for ENET-based filtering)}

\item{verbose}{Numeric variable indicating amount of information to print
to the terminal (0 = nothing, 1 = notes only, 2+ = notes and intermediate output)}

\item{filter.hyper}{List of control parameters for filtering step
(see also \code{\link{filter_control}}); key agruments include: \itemize{
\item method: filtering method; current options are: "ENET" (Elastic Net),
"RF" (Random Forest with default or input parameters, using delete-d
jackknife confidence intervals for VIMP for variable selection), "RFbest"
(Random Survival Forest, performing tuning for mtry and nodesize, and using
double bootstrap confidence intervals for VIMP for variable selection; more
accurate but slower than RF option). See \code{\link[glmnet]{glmnet}} for more
details on the elastic net and \code{\link[randomForestSRC]{rfsrc}} for more
details on the random forest filtering
 \item lambdatype: Optional elastic net parameter; whether to use the tuning
 parameter lambda that minimizes cross validation error
 (lambdatype="min", default) or the largest lambda that gives error within 1
 standard error of the minimum error (lambdatype="1se"). Ignored when method
 is "RF" or "RFbest"
\item mixparm: Optional elastic net mixing parameter alpha or grid of alpha
values to search over. If nothing is entered, the algorithm will perform a
grid search for the best mixing parameter between 0.05 and 0.95.
Ignored when method is "RF" or "RFbest"
 \item vimpalpha: Optional significance level for RF VIMP confidence intervals
 \item nfolds: number of folds to use for cross validation tuning of lambda
 parameter for elastic net filtering. Default = 10. Ignored when method="RF"
 or "RFbest"
 \item filterseed: optional seed for filtering step
 \item ... : Optional arguments for glmnet or rfsrc
}}

\item{vars2keep:}{List of variable names (matching column names in X) of
variables to be passed through filtering step without penalization, etc.
This is ignored in the main 5-STAR algorithm but may be used when filtering
step is used as a stand alone function. Currently only used when method = ENET}
}
\value{
cov2keep: List of all covariates kept after filtering step,
selected by elastic net or random forest

For method=ENET, additionally outputs \itemize{
   \item cvout: vector containing fraction of null deviance explained,
   mean cross validation error, and optimal tuning parameter values
    (see \code{\link[glmnet]{cv.glmnet}} and \code{\link[glmnet]{glmnet}}
     for more details)
   \item beta: the coefficients of glmnet fit given tuned parameters
   \item ENETplot: plot containing the deviance over different combinations of
   the mixing and tuning parameters alpha and lambda, with optimal
   alpha shown in blue and minimum deviance point for each alpha in black
   (left panel), and solution path for the best mixing parameter (right panel).
   Output when plot = TRUE
}

For method=RF or RFbest, additionally outputs:
\itemize{
   \item varselectmat:  matrix of VIMP confidence intervals, p-values,
   and selection decision for each variable
   \item VIMPplot: default variable importance CI plot, output if plot = TRUE
   \item varselect2plot: variable selection information from subsample.rfsrc for
making customizable VIMP CI plots
   \item forest: rfsrc object
}
}
\description{
Performs Step 2 of the 5-STAR Algorithm:
Forms dummy variable matrix for factors, and fits an elastic net (ENET) or
random forest (RF) model to determine which covariates to keep for building
trees in 5-STAR
}
