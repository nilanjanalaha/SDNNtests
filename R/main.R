################# The asymptotic tests #######################################################


#################### gives the asymptotic test statistics for log-concave tests based on distance #######

logcon_a<- function(x,y, pr)
{

  m <- length(x)
  n <- length(y)
  N <- m+n
  a <- ceiling(pr*N) #Truncation
  b <- ceiling((1-pr)*N) #Truncation
  pool1 <- sort(c(x,y))

  mlef1 <- logcondens::logConDens(x,smoothed = FALSE)
  mlef2 <- logcondens::logConDens(y,smoothed = FALSE)

  pool <- sort(c(mlef1$knots,mlef2$knots))

  ph1 <-  logcondens::evaluateLogConDens(pool,mlef1,1)[,2]
  ph2 <-  logcondens::evaluateLogConDens(pool,mlef2,1)[,2]
  php1 <- diff(ph1)/diff(pool)
  php2 <- diff(ph2)/diff(pool)

  # checking the knots between where F2-F1 attains maximum in some intermediate points, not at the kbnots
  ind <- which(php1*php2<0)
  # Combining the points where to have checked
  poolt1 <- pool[pool>=pool[a] & pool<=pool[b]]
  poolt <- c(poolt1,seq(pool1[a],pool1[b], by=0.05))
  poolt <- sort(poolt)

  #T1
  F2 <- logcondens::evaluateLogConDens(poolt,mlef2,3)[,4]
  F1 <- logcondens::evaluateLogConDens(poolt,mlef1,3)[,4]
  diff.lc <- F2 - F1
  t1 <- (F2-F1)/sqrt(F2*(1-F2)/n+F1*(1-F1)/m)
  t1 <- min(t1)


  z <- seq(pr,1-pr,0.05)
  z1 <- c(1:N)/N
  z1 <- z1[z1>pr & z1 < 1-pr]
  z <- sort(unique(c(z,z1)))

  #Calculating T4
  pool <- quantile(sort(c(x,y)),z)
  F1 <- logcondens::evaluateLogConDens(pool,mlef1,3)[,4]
  F2 <- logcondens::evaluateLogConDens(pool,mlef2,3)[,4]
  t4 <- sqrt(m*n/(m+n))*(F2-F1)/sqrt(z*(1-z))
  t4 <- min(t4)
  # t2 and t4 stand for the truncated ones
  c(t1,t4)

  c(t1,t4)

}

############# gives statistic of the unimodal tests based on the difference #######
unimod_a <- function(x,y,p, pr)
{


  m <- length(x)
  n <- length(y)
  N <- m+n
  eta <- 1/N^p
  a <- ceiling(pr*N) #truncation
  b <- ceiling((1-pr)*N) #truncation

  s1 <- calc_mode(sort(x),eta)
  s2 <- calc_mode(sort(y),eta)

  poolin <- sort(c(x,y))
  pool <- c(poolin[a:b], seq(poolin[a], poolin[b], 0.05))
  pool <- sort(pool)
  F1 <-  grenander.inter(pool,s1$x.knots, s1$F.knots, s1$f.knots)
  F2 <-  grenander.inter(pool,s2$x.knots, s2$F.knots, s2$f.knots)
  tmp1 <- (F2-F1)/sqrt(F2*(1-F2)/n+F1*(1-F1)/m)
  t1 <- min(tmp1)

  z <- seq(pr,1-pr,0.05)
  z1 <- c(1:N)/N
  z1 <- z1[a:b]
  z <- sort(unique(c(z,z1)))

  pool <- quantile(poolin,z)
  F1 <-  grenander.inter(pool,s1$x.knots, s1$F.knots, s1$f.knots)
  F2 <-  grenander.inter(pool,s2$x.knots, s2$F.knots, s2$f.knots)
  t4 <- sqrt(m*n/(m+n))*(F2-F1)/sqrt(z*(1-z))
  t4 <- min(t4)
  # t2 and t4 stand for the truncated ones
  c(t1,t4)
}

############# gives statistic of the nonparametric tests based on the difference #######
nonparam_a <- function(x,y, pr)
{

  m <- length(x)
  n <- length(y)
  N <- m+n
  a <- ceiling(pr*N) #truncation
  b <- ceiling((1-pr)*N) #truncation
  pool <- sort(c(x,y))


  F1 <-  ecdf(x)(pool)
  F2 <-  ecdf(y)(pool)
  tmp1 <-(F2-F1)/sqrt(F2*(1-F2)/n+F1*(1-F1)/m)
  t1 <- min(tmp1[a:b])


  z <- seq(pr, 1-pr, 0.01)
  z1 <- c(1:N)/N
  z1 <- z1[a:b]
  z <- sort(unique(c(z,z1)))

  pool <- quantile(pool,z)
  F1 <-  ecdf(x)(pool)
  F2 <- ecdf(y)(pool)
  t4 <- sqrt(m*n/(m+n))*(F2-F1)/sqrt(z*(1-z))
  t4 <- min(t4)
  # t2 and t4 stand for the truncated ones
  c(t1,t4)
}

############## combining the two shape constrained tests ###################################

boot_test_a <- function(y1,y2, p, pr,Method)
{
  if(Method=="NP") temp <- nonparam_a(y1, y2, pr)
  if(Method=="UM") temp <- unimod_a(y1, y2, p, pr)
  if(Method=="LC") temp <- logcon_a(y1, y2, pr)
  temp
}

############## SDNN ################################
#' One-sided test of  stochastic dominance against the null of non-dominance
#'
#' Calculates the p-values of one-sided
#' tests of restricted stochastic dominance against the null of non-dominance.
#' The concerned tests are the minimum t-statistic test
#' and the two sample empirical process (TSEP)
#' test of  \emph{Laha et al., 2020}. Each test can be either nonparametric,
#' or semiparametric, ie.  using unimodality or log-concavity assumption on the underlying
#' densities.
#'
#'
#' @details Suppose \eqn{X_1,\ldots, X_m} and \eqn{Y_1,\ldots, y_n} are independent random
#' variables with distribution \eqn{F} and \eqn{G}, respectively. Denote by
#' \eqn{D_{p,m,n}} the set \deqn{[H_n^{-1}(p), H_n^{-1}(1-p)]} where
#' \eqn{p\in[0,0.5)} and \eqn{H_n} is the empirical distribution
#' function of the combined sample \deqn{\{X_1,\ldots, x_m,Y_1,\ldots,Y_n\}.} The function SDNN tests
#'  \eqn{H_0: F(x)\geq G(x)} for some \eqn{x\in D_{p,m,n}} vs
#'  \eqn{H_1: F(x)<G(x)} for all \eqn{x\in D_{p,m,n}}.
#'  For more details, see Laha et al., 2020.
#'
#' @details  \code{Method}:
#'             "NP" corresponds to the nonparametric tests. "UM" corresponds
#'              to tests which use the function \code{\link{calc_mode}}
#'              to estimate the densities of \eqn{X_i}'s and \eqn{Y_i}'s. This function
#'              estimates the unimodal density estimator of Birge (1997).
#'               "LC" corresponds to tests which use  the log-concave MLE of Dumbgen and Rufibatch (2009)
#'               to estimate the latter densities. For more detail, see Laha et al., 2020.
#' @details  \code{t}:   The parameter \eqn{t}
#'                 corresponds to the parameter \eqn{\tau}
#'                 in Birge (1997). Higher values of \eqn{t} leads to more accurate
#'                  estimation of the unimodal densities.  This value ontrols the accuracy in
#'            unimodal density estimation upto the term \eqn{(m+n)^{-t}}. A value greater
#'              than or equal to one is recommended.  See Birge (1997)
#'            for more information.
#' @details \code{p}:   p corresponds to the set \eqn{D_{p,m,n}} in \eqn{H_0} and \eqn{H_1}.
#'            To overcome the difficulties arising from the tail region,
#'           \eqn{100p} percent of data is trimmed from both sides of the
#'           combined sample.
#'
#' @param x Vector of m independent and identically distributed random variables;
#'           corresponds to the first sample.
#' @param  y Vector of n independent and identically distributed random variables;
#'           corresponds to the second sample.
#' @param Method Must be one among "NP" (nonparametric), "UM" (unimodal), and "LC" (log-concave). See 'Details'.
#' @param t  A positive real number. Only required when Method={"UM"}, default value is 1.
#'           See Details.
#' @param p The proportion of combined data to be trimmed from each end prior to testing,
#'           should take value in \eqn{[0,0.50)}, the default is set to 0.05.
#'
#'@return  A list of two numbers.
#'\itemize{
#'\item T1 - The p-value of the test based on minimum T-statistic of
#'            Laha \emph{et al.}, (2020)
#'\item T2 - The p-value of the TSEP test of Laha \emph{et al.}, (2020).
#' }
#'@references Laha, N., Moodie, Z., Huang, Y., and Luedtke, A.
#' \emph{ Improved inference for vaccine-induced immune responses
#'        via shape-constrained methods}. Submitted.
#'@references Dumbgen, L. and Rufibatch, K. (2009). \emph{Maximum likelihood estimation of a logconcave density and
#'           its distribution function: Basic properties and uniform
#'           consistency}, Bernoulli, 15, 40–68.
#'@references Birge, L. (1997). \emph{Estimation of unimodal densities
#'          without smoothness assumptions}, Ann. Statist., 25, 970–981.
#'@author \href{https://connects.catalyst.harvard.edu/Profiles/display/Person/184207}{Nilanjana Laha}
#' (maintainer), \email{nlaha@@hsph.harvard.edu},
#'
#' Zoe Moodie, \email{zoe@@fredhutch.org},
#'
#' Y Huang, \email{yhuang@@fredhutch.org},
#'
#' Alex Luedtke, \email{aluedtke@@uw.edu}.
#' @seealso \code{\link{calc_mode}}
#' @examples
#' x <- rnorm(100); y <- rgamma(50, shape=1);
#'  SDNN(x, y, Method="UM", t=1, p=0.01)
#' @export
SDNN <- function(x, y, Method,  t, p)
{
  #setting default values of pr and p
  if(missing(p)) p <- 0.05
  if(missing(t)) t <- 1

  x <- sort(x)
  y <- sort(y)
  test_stat <- boot_test_a(x, y, t, p,Method)
  #Calculating the p-value
  pval <- 1-pnorm(test_stat)
  list(T1=pval[1],T2=pval[2])
}



