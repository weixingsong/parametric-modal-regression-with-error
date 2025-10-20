#' The Betsch-Ebner goodness-of-fit test for the gamma family
#'
#' @description
#' This function computes the goodness-of-fit test for the gamma family due to Betsch and Ebner (2019).
#'
#' @param data  a vector of positive numbers.
#' @param a     positive tuning parameter.
#' @param boot  number of bootstrap iterations used to obtain critical value.
#' @param alpha level of significance of the test.
#'
#' @return a list containing the value of the test statistic, the approximated critical value and a test decision on the significance level \code{alpha}: \cr
#' \describe{
#'         \item{\code{$T.value}}{the value of the test statistic.}
#'         \item{\code{$cv}}{the approximated critical value.}
#'         \item{\code{$par.est}}{number of points used in approximation.}
#'         \item{\code{$Decision}}{the comparison of the critical value and the value of the test statistic.}
#'         \item{\code{$sig.level}}{level of significance chosen.}
#'         \item{\code{$boot.run}}{number of bootstrap iterations.}
#'}
#'
#'
#' @details
#' The test is of weighted \eqn{L^2} type and uses a characterization of the distribution function of the gamma distribution. Critical values are obtained by a parametric bootstrap procedure, see \code{\link{crit.values}}.
#'
#' @references
#' Betsch, S., Ebner, B. (2019) "A new characterization of the Gamma distribution and associated goodness of fit tests", Metrika, 82(7):779-806. \href{https://doi.org/10.1007/s00184-019-00708-7}{DOI}
#'
#' @examples
#' test.BE(stats::rgamma(20,3,6),boot=100)
#'
#' @export
test.BE<-function(data,a=1,boot=500,alpha=0.05)
{
  if (!is.vector(data)){warning("data must be of type vector")}
  para=gamma_est(data)
  Test.value=BE(data/para[2],para[1],a)
  cv=crit.values(samplesize=length(data),statistic=BE,tuning=a,k_estimator=para[1],boot.param=boot,alpha=alpha)
  result <- list("Test" = "BE","parameter" = a, "T.value"=Test.value,"cv"=cv,"par.est"=para,"Decision"=Test.value>cv,"sig.level"=alpha,"boot.run"=boot)
  attr(result, "class") <- "gofgamma"
  return(result)
}


#' The Kolmogorov-Smirnov goodness-of-fit test for the gamma family
#'
#' @description
#' This function computes the goodness-of-fit test for the gamma family in the spirit of Kolmogorov and Smirnov. Note that this tests the composite hypothesis of fit to the family of gamma distributions, i.e. a bootstrap procedure is implemented to perform the test.
#'
#' @param data  a vector of positive numbers.
#' @param boot  number of bootstrap iterations used to obtain critical value.
#' @param alpha level of significance of the test.
#'
#' @return a list containing the value of the test statistic, the approximated critical value and a test decision on the significance level \code{alpha}: \cr
#' \describe{
#'         \item{\code{$T.value}}{the value of the test statistic.}
#'         \item{\code{$cv}}{the approximated critical value.}
#'         \item{\code{$par.est}}{number of points used in approximation.}
#'         \item{\code{$Decision}}{the comparison of the critical value and the value of the test statistic.}
#'         \item{\code{$sig.level}}{level of significance chosen.}
#'         \item{\code{$boot.run}}{number of bootstrap iterations.}
#'}
#'
#'
#' @details
#' The Kolmogorov Smirnov test is computed as described in Henze et. al. (2012). Critical values are obtained by a parametric bootstrap procedure, see \code{\link{crit.values}}.
#'
#' @references
#' Henze, N., Meintanis, S.G., Ebner, B. (2012) "Goodness-of-fit tests for the Gamma distribution based on the empirical Laplace transform". Communications in Statistics - Theory and Methods, 41(9): 1543-1556. \href{https://doi.org/10.1080/03610926.2010.542851}{DOI}
#'
#' @examples
#' test.KS(stats::rgamma(20,3,6),boot=100)
#'
#' @export
test.KS <- function(data,boot=500,alpha=0.05){
  if (!is.vector(data)){warning("data must be of type vector")}
  para=gamma_est(data)
  cv=crit.values(samplesize=length(data),statistic=KS,tuning=NULL,k_estimator=para[1],boot.param=boot,alpha=alpha)$p_n
  Tn=crit.values(samplesize=length(data),statistic=KS,tuning=NULL,k_estimator=para[1],boot.param=boot,alpha=alpha)$Tn
  Test.value=KS(data/para[2],para[1])
  pvalue=mean(Tn>=Test.value)
  result <- list("Test" = "KS","parameter" = NULL, "T.value"=Test.value,"cv"=cv,"par.est"=para,"Decision"=Test.value>cv,"pvalue"=pvalue,"sig.level"=alpha,"boot.run"=boot)
  attr(result, "class") <- "gofgamma"
  return(result)
}



#' The Cramer-von Mises goodness-of-fit test for the gamma family
#'
#' @description
#' This function computes the goodness-of-fit test for the gamma family in the spirit of Cramer and von Mises. Note that this tests the composite hypothesis of fit to the family of gamma distributions, i.e. a bootstrap procedure is implemented to perform the test.
#'
#' @param data  a vector of positive numbers.
#' @param boot  number of bootstrap iterations used to obtain critical value.
#' @param alpha level of significance of the test.
#'
#' @return a list containing the value of the test statistic, the approximated critical value and a test decision on the significance level \code{alpha}: \cr
#' \describe{
#'         \item{\code{$T.value}}{the value of the test statistic.}
#'         \item{\code{$cv}}{the approximated critical value.}
#'         \item{\code{$par.est}}{number of points used in approximation.}
#'         \item{\code{$Decision}}{the comparison of the critical value and the value of the test statistic.}
#'         \item{\code{$sig.level}}{level of significance chosen.}
#'         \item{\code{$boot.run}}{number of bootstrap iterations.}
#'}
#'
#'
#' @details
#' The Cramér-von Mises test is computed as described in Henze et. al. (2012). Critical values are obtained by a parametric bootstrap procedure, see \code{\link{crit.values}}.
#'
#' @references
#' Henze, N., Meintanis, S.G., Ebner, B. (2012) "Goodness-of-fit tests for the Gamma distribution based on the empirical Laplace transform". Communications in Statistics - Theory and Methods, 41(9): 1543-1556. \href{https://doi.org/10.1080/03610926.2010.542851}{DOI}
#'
#' @examples
#' test.CM(stats::rgamma(20,3,6),boot=100)
#'
#' @export
test.CM <- function(data,boot=500,alpha=0.05){
  para=gamma_est(data)
  cv=crit.values(samplesize=length(data),statistic=CM,tuning=NULL,k_estimator=para[1],boot.param=boot,alpha=alpha)$p_n
  Tn=crit.values(samplesize=length(data),statistic=CM,tuning=NULL,k_estimator=para[1],boot.param=boot,alpha=alpha)$Tn
  Test.value=CM(data/para[2],para[1])
  pvalue=mean(Tn>=Test.value)
  result <- list("Test" = "CM","parameter" = NULL, "T.value"=Test.value,"cv"=cv,"par.est"=para,"Decision"=Test.value>cv,"pvalue"=pvalue,"sig.level"=alpha,"boot.run"=boot)
  attr(result, "class") <- "gofgamma"
  return(result)
}



#' The Anderson-Darling goodness-of-fit test for the gamma family
#'
#' @description
#' This function computes the goodness-of-fit test for the gamma family in the spirit of Anderson and Darling. Note that this tests the composite hypothesis of fit to the family of gamma distributions, i.e. a bootstrap procedure is implemented to perform the test.
#'
#' @param data  a vector of positive numbers.
#' @param boot  number of bootstrap iterations used to obtain critical value.
#' @param alpha level of significance of the test.
#'
#' @return a list containing the value of the test statistic, the approximated critical value and a test decision on the significance level \code{alpha}: \cr
#' \describe{
#'         \item{\code{$T.value}}{the value of the test statistic.}
#'         \item{\code{$cv}}{the approximated critical value.}
#'         \item{\code{$par.est}}{number of points used in approximation.}
#'         \item{\code{$Decision}}{the comparison of the critical value and the value of the test statistic.}
#'         \item{\code{$sig.level}}{level of significance chosen.}
#'         \item{\code{$boot.run}}{number of bootstrap iterations.}
#'}
#'
#'
#' @details
#' The Anderson-Darling  test is computed as described in Henze et. al. (2012). Critical values are obtained by a parametric bootstrap procedure, see \code{\link{crit.values}}.
#'
#' @references
#' Henze, N., Meintanis, S.G., Ebner, B. (2012) "Goodness-of-fit tests for the Gamma distribution based on the empirical Laplace transform". Communications in Statistics - Theory and Methods, 41(9): 1543-1556. \href{https://doi.org/10.1080/03610926.2010.542851}{DOI}
#'
#' @examples
#' test.AD(stats::rgamma(20,3,6),boot=100)
#'
#' @export
test.AD <- function(data,boot=500,alpha=0.05){
  if (!is.vector(data)){warning("data must be of type vector")}
  para=gamma_est(data)
  cv=crit.values(samplesize=length(data),statistic=AD,tuning=NULL,k_estimator=para[1],boot.param=boot,alpha=alpha)
  Test.value=AD(data/para[2],para[1])
  result <- list("Test" = "AD","parameter" = NULL, "T.value"=Test.value,"cv"=cv,"par.est"=para,"Decision"=Test.value>cv,"sig.level"=alpha,"boot.run"=boot)
  attr(result, "class") <- "gofgamma"
  return(result)
}



#' The Watson goodness-of-fit test for the gamma family
#'
#' @description
#' This function computes the goodness-of-fit test for the gamma family in the spirit of Watson. Note that this tests the composite hypothesis of fit to the family of gamma distributions, i.e. a bootstrap procedure is implemented to perform the test.
#'
#' @param data  a vector of positive numbers.
#' @param boot  number of bootstrap iterations used to obtain critical value.
#' @param alpha level of significance of the test.
#'
#' @return a list containing the value of the test statistic, the approximated critical value and a test decision on the significance level \code{alpha}: \cr
#' \describe{
#'         \item{\code{$T.value}}{the value of the test statistic.}
#'         \item{\code{$cv}}{the approximated critical value.}
#'         \item{\code{$par.est}}{number of points used in approximation.}
#'         \item{\code{$Decision}}{the comparison of the critical value and the value of the test statistic.}
#'         \item{\code{$sig.level}}{level of significance chosen.}
#'         \item{\code{$boot.run}}{number of bootstrap iterations.}
#'}
#'
#'
#' @details
#' The Watson test is computed as described in Henze et. al. (2012). Critical values are obtained by a parametric bootstrap procedure, see \code{\link{crit.values}}.
#'
#' @references
#' Henze, N., Meintanis, S.G., Ebner, B. (2012) "Goodness-of-fit tests for the Gamma distribution based on the empirical Laplace transform". Communications in Statistics - Theory and Methods, 41(9): 1543-1556. \href{https://doi.org/10.1080/03610926.2010.542851}{DOI}
#'
#' @examples
#' test.WA(stats::rgamma(20,3,6),boot=100)
#'
#' @export
test.WA <- function(data,boot=500,alpha=0.05){
  para=gamma_est(data)
  cv=crit.values(samplesize=length(data),statistic=WA,tuning=NULL,k_estimator=para[1],boot.param=boot,alpha=alpha)
  Test.value=WA(data/para[2],para[1])
  result <- list("Test" = "WA","parameter" = NULL, "T.value"=Test.value,"cv"=cv,"par.est"=para,"Decision"=Test.value>cv,"sig.level"=alpha,"boot.run"=boot)
  attr(result, "class") <- "gofgamma"
  return(result)
}



#' The first Henze-Meintanis-Ebner goodness-of-fit test for the gamma family
#'
#' @description
#' This function computes the first goodness-of-fit test for the gamma family due to Henze, Meintanis and Ebner (2012).
#'
#' @param data  a vector of positive numbers.
#' @param a     positive tuning parameter.
#' @param boot  number of bootstrap iterations used to obtain critical value.
#' @param alpha level of significance of the test.
#'
#' @return a list containing the value of the test statistic, the approximated critical value and a test decision on the significance level \code{alpha}: \cr
#' \describe{
#'         \item{\code{$T.value}}{the value of the test statistic.}
#'         \item{\code{$cv}}{the approximated critical value.}
#'         \item{\code{$par.est}}{number of points used in approximation.}
#'         \item{\code{$Decision}}{the comparison of the critical value and the value of the test statistic.}
#'         \item{\code{$sig.level}}{level of significance chosen.}
#'         \item{\code{$boot.run}}{number of bootstrap iterations.}
#'}
#'
#'
#' @details
#' The test is of weighted \eqn{L^2} type and uses a characterization of the distribution function of the gamma distribution. Critical values are obtained by a parametric bootstrap procedure, see \code{\link{crit.values}}.
#'
#' @references
#' Henze, N., Meintanis, S.G., Ebner, B. (2012) "Goodness-of-fit tests for the Gamma distribution based on the empirical Laplace transform". Communications in Statistics - Theory and Methods, 41(9): 1543-1556. \href{https://doi.org/10.1080/03610926.2010.542851}{DOI}
#'
#' @examples
#' test.HME1(stats::rgamma(20,3,6),boot=100)
#'
#' @export
test.HME1 <-function(data,a=1,boot=500,alpha=0.05){
  if (!is.vector(data)){warning("data must be of type vector")}
  para=gamma_est(data)
  cv=crit.values(samplesize=length(data),statistic=HME1,tuning=a,k_estimator=para[1],boot.param=boot,alpha=alpha)
  Test.value=HME1(data/para[2],para[1],a)
  result <- list("Test" = "HME1","parameter" = a, "T.value"=Test.value,"cv"=cv,"par.est"=para,"Decision"=Test.value>cv,"sig.level"=alpha,"boot.run"=boot)
  attr(result, "class") <- "gofgamma"
  return(result)
}



#' The second Henze-Meintanis-Ebner goodness-of-fit test for the gamma family
#'
#' @description
#' This function computes the second goodness-of-fit test for the gamma family due to Henze, Meintanis and Ebner (2012).
#'
#' @param data  a vector of positive numbers.
#' @param a     positive tuning parameter.
#' @param boot  number of bootstrap iterations used to obtain critical value.
#' @param alpha level of significance of the test.
#'
#' @return a list containing the value of the test statistic, the approximated critical value and a test decision on the significance level \code{alpha}: \cr
#' \describe{
#'         \item{\code{$T.value}}{the value of the test statistic.}
#'         \item{\code{$cv}}{the approximated critical value.}
#'         \item{\code{$par.est}}{number of points used in approximation.}
#'         \item{\code{$Decision}}{the comparison of the critical value and the value of the test statistic.}
#'         \item{\code{$sig.level}}{level of significance chosen.}
#'         \item{\code{$boot.run}}{number of bootstrap iterations.}
#'}
#'
#'
#' @details
#' The test is of weighted \eqn{L^2} type and uses a characterization of the distribution function of the gamma distribution. Critical values are obtained by a parametric bootstrap procedure, see \code{\link{crit.values}}.
#'
#' @references
#' Henze, N., Meintanis, S.G., Ebner, B. (2012) "Goodness-of-fit tests for the Gamma distribution based on the empirical Laplace transform". Communications in Statistics - Theory and Methods, 41(9): 1543-1556. \href{https://doi.org/10.1080/03610926.2010.542851}{DOI}
#'
#' @examples
#' test.HME2(stats::rgamma(20,3,6),boot=100)
#'
#' @export
test.HME2 <-function(data,a=4,boot=500,alpha=0.05){
  if (!is.vector(data)){warning("data must be of type vector")}
  para=gamma_est(data)
  cv=crit.values(samplesize=length(data),statistic=HME2,tuning=a,k_estimator=para[1],boot.param=boot,alpha=alpha)
  Test.value=HME2(data/para[2],para[1],a)
  result <- list("Test" = "HME2","parameter" = a, "T.value"=Test.value,"cv"=cv,"par.est"=para,"Decision"=Test.value>cv,"sig.level"=alpha,"boot.run"=boot)
  attr(result, "class") <- "gofgamma"
  return(result)
}
