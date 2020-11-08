library(logcondens)
library(pracma)




############### calculates the Hellinger distance between KDE estimators from ################
###################################### two samples  #######################################

# low and up are the lower and upper bounds of f1, can be infinity

Hell <- function(f1, f2, low, up)
{
  h <- function(x) (sqrt(max(0,f1(x)))-sqrt(max(0,f2(x))))^2/2
  pracma::integral(Vectorize(h), low, up,method="Kronrod")
}



##################### Gives Hellinger distance for Birge's estimator ###########################
#' Hellinger distance between unimodal densities
#'
#' Provides an estimate and a confidence interval for
#' the hellinger distance between two unimodal densities.
#' This function uses the density estimator of Birges (1997), given by
#' \code{\link{calc_mode}}. See Laha et al. (2020) for more details.
#'
#' @details  This function calls \code{link{calc_mode}} where the parameter
#'           t is taken to be one.
#'
#' @param x Vector of m independent and identically distributed random variables;
#'           corresponds to the first sample.
#' @param  y Vector of n independent and identically distributed random variables;
#'           corresponds to the second sample.
#' @param alpha A positive number between 0 and 0.50. The function returns an
#'              (1-alpha)\eqn{\%} confidence interval. The default value is 0.05.
#'
#'@return  A vector of three elements. The first element is the point estimate
#'         of the Hellinger distance. The second and the third elements give the
#'          left and right endpoints of the confidence interval.
#'
#'@references Laha, N., Moodie, Z., Huang, Y., and Luedtke, A. (2020).
#' \emph{ Improved inference for vaccine-induced immune responses
#'        via shape-constrained methods}. Submitted.
#'@references Birge, L. (1997). \emph{Estimation of unimodal densities
#'          without smoothness assumptions}, Ann. Statist., 25, 970–981.
#'@author \href{https://connects.catalyst.harvard.edu/Profiles/display/Person/184207}{Nilanjana Laha}
#' (maintainer), \email{nlaha@@hsph.harvard.edu},
#'
#' Alex Luedtke, \email{aluedtke@@uw.edu}.
#' @seealso \code{\link{calc_mode}}, \code{\link{hd.lc}}, \code{\link{hd.lc.sm}}
#' @examples
#' x <- sort(rnorm(100)); y <- sort(rgamma(50, shape=1));
#' hd.uni(x, y, 0.05)
#' @export
hd.uni <- function(x,y, alpha=0.05)
{
  x <- sort(x)
  y <- sort(y)
  N <- length(x)+length(y)
  eta <- 1/N

  s1 <- calc_mode(x,eta)
  s2 <- calc_mode(y,eta)

  pool <- sort(c(s1$x.knots, s2$x.knots))
  v1 <-  grenander.inter.density(pool,s1$x.knots, s1$F.knots, s1$f.knots)
  v2 <-  grenander.inter.density(pool,s2$x.knots, s2$F.knots, s2$f.knots)


  f1x <-  grenander.inter.density(x,s1$x.knots, s1$F.knots, s1$f.knots)
  f1y <-  grenander.inter.density(y,s1$x.knots, s1$F.knots, s1$f.knots)

  f2x <-  grenander.inter.density(x,s2$x.knots, s2$F.knots, s2$f.knots)
  f2y <-  grenander.inter.density(y,s2$x.knots, s2$F.knots, s2$f.knots)

  m <- length(pool)
  s <- 0

  for ( i in 1:(m-1))
  {
    a <- pool[i]
    b <- pool[i+1]
    s=s+ sqrt(v1[i+1]*v2[i+1])*(b-a)

  }
  H <- 1-s # the hellinger distance
  ci <- ci.birge(f1x,f1y,f2x,f2y,v1,v2,pool,m,N,H, alpha)
  c(H,ci)
}

################## gives the the confidence interval for Birge's estimator #####################

ci.birge <- function(f1x,f1y,f2x,f2y,v1,v2,pool,m,N,H, alpha)
{
  scf <- 0
  scg <- 0
  for ( i in 1:(m-1))
  {
    a <- pool[i]
    b <- pool[i+1]
    scf <- scf + ifelse(v1[i+1]>0,(1-sqrt(v2[i+1]/v1[i+1])-H)^2*v1[i+1],0)*(b-a)/4
    scg <- scg + ifelse(v2[i+1]>0,(1-sqrt(v1[i+1]/v2[i+1])-H)^2*v2[i+1],0)*(b-a)/4
  }

  l <- length(f1x)/N
  var <- scf/l+scg/(1-l)
  sd <- -sqrt(var/N)*qnorm(alpha/2)

  sm.ci <- c(H-sd,H+sd) #ci constructed using f-hat to smooth out
  sm.ci
}


################### Gives Hellinger distance for log-concave density estimator ##################
#' Hellinger distance between log-concave densities
#'
#' Provides an estimate of
#' the hellinger distance between two log-concave densities.
#' This function uses the log-concave density estimator of Dumbgen and Rufibatch (2009),
#' given by \link[logcondens]{logConDens} of logcondens package.
#'
#'
#'
#' @param x Vector of m independent and identically distributed random variables;
#'           corresponds to the first sample.
#' @param  y Vector of n independent and identically distributed random variables;
#'           corresponds to the second sample.
#'
#'
#'@return  A point estimator of the Hellinger distance.
#'
#'@references Laha, N., Moodie, Z., Huang, Y., and Luedtke (2020), A.
#' \emph{ Improved inference for vaccine-induced immune responses
#'        via shape-constrained methods}. Submitted.
#' @references Dumbgen, L. and Rufibatch, K. (2009). \emph{Maximum likelihood estimation of a logconcave density and
#'           its distribution function: Basic properties and uniform
#'           consistency}, Bernoulli, 15, 40–68.
#'@author \href{https://connects.catalyst.harvard.edu/Profiles/display/Person/184207}{Nilanjana Laha}
#' (maintainer), \email{nlaha@@hsph.harvard.edu},
#'
#' Alex Luedtke, \email{aluedtke@@uw.edu}.
#' @seealso \code{\link{hd.uni}}, \code{\link{hd.lc.sm}}
#' @examples
#' x <- sort(rnorm(100)); y <- sort(rgamma(50, shape=1));
#' hd.lc(x,y)
#' @export
hd.lc <- function(x,y)
{

  mlex <- logcondens::logConDens(x, smoothed=FALSE)
  mley <- logcondens::logConDens(y, smoothed =FALSE)

  f1 <- function(x) logcondens::evaluateLogConDens(x,mlex,which=2)[,3]
  f2 <- function(x) logcondens::evaluateLogConDens(x,mley,which=2)[,3]

  Hell(f1, f2 ,min(x,y), max(x,y))
}
##################### smoothed log-concave estimator ###################################
#' Hellinger distance between log-concave densities
#'
#' Provides an estimate of
#' the Hellinger distance between two log-concave densities.
#' This function uses the smoothed log-concave density estimator of Dumbgen and Rufibatch (2009),
#' given by \link[logcondens]{logConDens} of logcondens package.
#'
#'
#'
#' @param x Vector of m independent and identically distributed random variables;
#'           corresponds to the first sample.
#' @param  y Vector of n independent and identically distributed random variables;
#'           corresponds to the second sample.
#'
#'
#'@return  A point estimator of the Hellinger distance.
#'
#'@references Laha, N., Moodie, Z., Huang, Y., and Luedtke (2020), A.
#' \emph{ Improved inference for vaccine-induced immune responses
#'        via shape-constrained methods}. Submitted.
#' @references Dumbgen, L. and Rufibatch, K. (2009). \emph{Maximum likelihood estimation of a logconcave density and
#'           its distribution function: Basic properties and uniform
#'           consistency}, Bernoulli, 15, 40–68.
#'@author \href{https://connects.catalyst.harvard.edu/Profiles/display/Person/184207}{Nilanjana Laha}
#' (maintainer), \email{nlaha@@hsph.harvard.edu},
#'
#' Alex Luedtke, \email{aluedtke@@uw.edu}.
#' @seealso \code{\link{hd.lc}}, \code{\link{hd.uni}}
#' @examples
#' x <- sort(rnorm(100)); y <- sort(rgamma(50, shape=1));
#' hd.lc.sm(x,y)
#' @export
hd.lc.sm <- function(x,y)
{

  mlex <- logcondens::logConDens(x, smoothed = TRUE)
  mley <- logcondens::logConDens(y)

  f1 <- function(x) logcondens::evaluateLogConDens(x,mlex,which=4)[,5]
  f2 <- function(x) logcondens::evaluateLogConDens(x,mley,which=4)[,5]
  t <- max(sd(x),sd(y))
  Hell(f1, f2 ,min(x,y)-5*t, max(x,y)+5*t)
}
