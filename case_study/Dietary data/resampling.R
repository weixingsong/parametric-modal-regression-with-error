#' bootstrap critical value of statistic
#'
#' @param samplesize  number of observations in the sample
#' @param statistic   test statistic to be used
#' @param tuning      tuning parameter used for the test statistic (\code{NULL} stands for no tuning parameter needed)
#' @param k_estimator value of the estimated \code{shape} parameter
#' @param boot.param  number of bootstrap iterations
#' @param alpha       significance level of the test
#'
#' @return returns the critical value for the goodness-of-fit test using the \code{statistic}.
#'
#' @examples
#' crit.values(samplesize=20,statistic=HME1,tuning=1,k_estimator=2,boot.param=100,alpha=0.05)
#'
#'@export
crit.values <- function(samplesize,statistic,tuning=NULL,k_estimator,boot.param=500,alpha=0.05)
{
  n=samplesize
  y <- rep(0,boot.param)
  if(!is.null(tuning))
  {
  for (j in 1:boot.param)
    {
    BS <- stats::rgamma(n,k_estimator,1)
    BS_k_estimator <- gamma_est(BS)[1]
    BS_lambda_estimator <- mean(BS)/BS_k_estimator
    x_BS=BS/BS_lambda_estimator
    y[j]=statistic(x_BS,BS_k_estimator,tuning)
    }
  } else {
  for (j in 1:boot.param)
      {
      BS <- stats::rgamma(n,k_estimator,1)
      BS_k_estimator <- gamma_est(BS)[1]
      BS_lambda_estimator <- mean(BS)/BS_k_estimator
      x_BS=BS/BS_lambda_estimator
      y[j]=statistic(x_BS,BS_k_estimator)
      }
  }
  
  Tn <- sort(y)
  k=floor((1-alpha)*boot.param)
  p_n <-Tn[k]+(1-alpha)*(Tn[k+1]-Tn[k])
  return(list(p_n=p_n,Tn=Tn))
}

