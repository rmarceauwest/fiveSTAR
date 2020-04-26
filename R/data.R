#' 5-STAR Example Hypothetical Dataset
#'
#' Survival and covariate data for 5-STAR simulated hypothetical example
#'
#' @format A data frame with 600 subjects (rows) and 53 variables:
#' \describe{
#'   \item{time}{observed survival time, in months}
#'   \item{status}{binary censoring indicator where 1 indicates the true
#'   survival time is observed and 0 indicates the subject's time is censored}
#'   \item{arm}{binary vector of treatment assignment where 1 indicates
#'   receiving test treatment and 0 indicates receiving control treatment}
#'   \item{X1-X25}{25 candidate prognostic binary factor covariates}
#'   \item{X26-X50}{25 candidate prognostic continuous covariates}
#' }
"survdata"
