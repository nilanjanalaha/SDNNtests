
################### Gives confidence interval for Hellinger distance  ##################
#' Confidence interval for Hellinger distance estimators
#'
#' Calculates the confidence interval of the squared Hellinger distance, as
#' discussed in Laha et al. (2021). This function employs  shape
#' constrained methods to estimate the underlying densities. See below for more details.
#' which include the log-concave density estimator of Dumbgen and Rufibatch (2009),
#' and the smoothed log-concave density estimator
#' given by ,
#'
#'
#'
#' @param x Vector of m independent and identically distributed random variables;
#'           corresponds to the first sample.
#' @param  y Vector of n independent and identically distributed random variables;
#'           corresponds to the second sample.
#' @param alpha Level of significance; returns a \eqn{1-\alpha\%} confidence interval.
#'
#'@details The confidence intervals follow the construnction of Laha et al. (2021). One
#'  confidence interval is based on the unimodal density estimator of Birges (1997), and
#'  provably requires the underlying densities to be unimodal to work.
#' The other confidence intervals are based on the log-concave density estimators, i.e.
#' the log-concave MLE of Dumbgen and Rufibatch (2009), and the smoothed log-concave
#' MLE of Chen and Samworth (2013) .
#' Both are computed using the function \link[logcondens]{logConDens} of logcondens package.
#'
#'
#'@return  \itemize{
#'\item lc.ci - Confidence intervals based on the log-concave MLE.
#'\item sm.lc.ci - Confidence intervals based on the smoothed log-concave MLE.
#'\item um.ci - Confidence intervals based on unimodal density estimator.
#'}
#'
#'@references Laha, N., Moodie, Z., Huang, Y., and Luedtke (2020), A.
#' \emph{ Improved inference for vaccine-induced immune responses
#'        via shape-constrained methods}. Submitted.
#' @references Dumbgen, L. and Rufibatch, K. (2009). \emph{Maximum likelihood estimation of a logconcave density and
#'           its distribution function: Basic properties and uniform
#'           consistency}, Bernoulli, 15, 40–68.
#' @references Chen, Y. and Samworth, R. J. (2013). \emph{Smoothed
#'           log-concave maximum likelihood estimation with applications},
#'           Statistica Sinica, 23, 1303-1398.
#' @references Birge, L. (1997). \emph{Estimation of unimodal densities
#'          without smoothness assumptions}, Ann. Statist., 25, 970–981.
#'@author \href{https://connects.catalyst.harvard.edu/Profiles/display/Person/184207}{Nilanjana Laha}
#' (maintainer), \email{nlaha@@hsph.harvard.edu},
#'
#' Alex Luedtke, \email{aluedtke@@uw.edu}.
#' @seealso \code{\link{hd.uni}}, \code{\link{hd.lc.sm}}, \code{\link{hd.uni}}
#' @examples
#' x <- sort(rnorm(100)); y <- sort(rgamma(50, shape=1));
#' hell.ci(x,y)$lc.ci
#' @export
hell.ci <- function(x, y, alpha=0.05)
{
  #length of vectors
  m <- length(x)
  n <- length(y)

  #The confidence intervals
  lc.ci <- ci(hd.lc(x, y), m, n, alpha)
  sm.lc.ci <- ci(hd.lc.sm(x,y), m, n, alpha)
  um.ci <- ci(hd.uni(x,y), m, n, alpha)

  data.frame(lc.ci, sm.lc.ci, um.ci)
}

#The inner function
ci <- function(H, m, n, alpha)
{
  N <- m+n
  #Calculate the sd
  lambda <- m/N
  var <- (2*H-H^2)/(4*lambda*(1-lambda))

  #Wald type Confidence interval
  sd2 <- sqrt(var/N)*qnorm(1-alpha/2)
   c(H-sd2,H+sd2)
}


