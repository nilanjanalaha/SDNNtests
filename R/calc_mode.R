
# Calculating unimodal regulerization about r where x is the dataset

# x is the datapoints
# Linear interpolation function
# r is a point which is not a sample point, r is the mode
unimodal_Grenander <- function(r,x)
{
  n <- length(x)
  xp <- x[x>=r]
  xm <- x[x<r]
  al <- length(xp)/length(x)
  Fm <-  my_grenander(xm,"increasing",r,1-al)
  Fp <- my_grenander(xp,"decreasing",r,al)
  return(list(mode=r, x.knots=c(Fm$x.knots,Fp$x.knots[-1]), F.knots = c(Fm$F.knots,Fp$F.knots[-1]+1-al), f.knots = c(Fm$f.knots,Fp$f.knots)))
}

# Calculating ai and bi at y
# If calculated correctly, an should increase and bn should decrease
# ind says if r is adatapoint, ind =0 if its not
calcab <- function(r,x,flag)
{
  # x is sorted
  n <- length(x)

    xp <- x[x>=r]
    xm <- x[x<r]
    al <- length(xp)/length(x) # sample proportion as in Prakash Rao's paper greater than r
    imax <- length(which(x<=r))


    #finding ai
    a <- 1/n #when r=x[1]
    if(r>x[1])
    {
    Fm <-  my_grenander(xm,"increasing",r,1-al)
   # It contains 1-a as the end-point,  Why are you doing that???
    xmn <- c(xm,r)#
    Fhatm <- grenander.inter(xmn,Fm$x.knots, Fm$F.knots, Fm$f.knots) #interpolating

    if(flag==0)
    {
     v <- c(c(1:imax),imax)/n-Fhatm # if r is not an observation, does not need it
    } else    v <- c(1:imax)/n-Fhatm
    a <- max(max(v),0)
    }

    #finding bi
    b <- 1/n  #when r=x[n]
    if(r<x[n])
    {
    Fp <- my_grenander(xp,"decreasing",r,al)
     Fp.knot <- Fp$F.knots+1-al
       Fhatp <-  grenander.inter(xp,Fp$x.knots, Fp.knot, Fp$f.knots)
    if(flag==0)
     {  v <- Fhatp-c(imax:(n-1))/n } else # if r is not an observation, does not need it
   { v <- Fhatp[-1]-c(imax:(n-1))/n}
    b <- max(v,0)
    }
    c(a,b)
}

# Gives the unimodal density estimator for unknown mode
############## calc_mode ################################
#' Unimodal density estimator when the mode is unknown
#'
#' Estimates the density of a given sample under the assumption that
#' the underlying density is unimodal. The mode is estimated from
#' the data.
#' The method is based on the unimodal regularization of Birge (1997).
#'
#'Birge(1997)'s  estimator gives a a pieciewise constant unimodal
#' density. The discontinuity points of the respective density are called the knots, which belong
#' to the set of datapoints. The density estimator is constant between two
#' consecutive knots.
#'
#'
#' @details  \code{t}:   The parameter \eqn{t}
#'                 corresponds to the parameter \eqn{\tau}
#'                 in Birge (1997). Higher values of \eqn{t} leads to more accurate
#'                  estimation of the unimodal densities.  This value ontrols the accuracy in
#'            unimodal density estimation upto the term \eqn{n^{-t}} where \eqn{n} is the sample
#'            size. We recommend a value greater than or equal to one.  See Birge (1997)
#'            for more details.
#'
#'
#' @param x  Vector of independent and identically distributed random variables;
#'           must be sorted.
#' @param t  A positive real number. Default value is one.
#'
#'@return
#'\itemize{
#'\item mode - The estimator of the mode
#'\item x.knots - The vector of the knots of the estimated density.
#'\item F.knots - A vector consisting the values of the estimated
#'                 distribution function evaluated at the knots.
#'\item f.knots - A vector whose i-th element gives the value of the
#'                estimated density on the segment
#'                joining the i-th and (i+1)-th knot.
#'                Recall that the estimated density is piecewise constant between two knots.
#' }

#'@references Birge, L. (1997). \emph{Estimation of unimodal densities
#'          without smoothness assumptions}, Ann. Statist., 25, 970–981.
#'@author \href{https://connects.catalyst.harvard.edu/Profiles/display/Person/184207}{Nilanjana Laha}
#' (maintainer), \email{nlaha@@hsph.harvard.edu},
#' Alex Luedtke, \email{aluedtke@@uw.edu}.
#' @examples
#' x <- sort(rnorm(100)); calc_mode(x, 1)
#' @export
calc_mode<- function(x,t)
{
  y <- x
    res <- sapply(y,calcab,x=x,flag=1) #calculating ai and bi for each x
  an <- res[1,]
  bn <- res[2,]
  k=0
  while(an[k+1]<bn[k+1])
  {
   k=k+1
  }

  # dichotomy argument
  #setting the parameters
  as <- an[k]
  an <- an[k+1]
  bs <- bn[k]
  bn <- bn[k+1]
  st <- x[k]
  end <- x[k+1]
  d <- max(an,bs)-max(as,bn)
  count <- 0
  while(d>t)
  {
    count <- count+1
    y <- (st+end)/2
    res <- calcab(y,x,0)
    newa <- res[1]
    newb <- res[2]

    #Update parameters
    if(newa>newb)
    {
      an <- newa
      bn <- newb
      end <- y
    } else { as <- newa
      bs <- newb
      st <- y
    }
    d <- max(an,bs)-max(as,bn)
  }
  # Selecting the mode
  #we know the number of iteration
  m <- (st+end)/2
  unimodal_Grenander(m,x)
}

########################### infrimum of Difference of two samples ##################################################################



