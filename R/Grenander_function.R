# x is the dataset
# a is where the support starts
Grenander_general<- function(x,a)
{
  n <- length(x) #n=sample size
  F.knots <- array(0) #Array for y.knots of lcm of ecdf

  if(x[1]>a)# This means that end-point is not included in the data
  {
  y <- c(a,unique(x)) # x knots from empirical cdf
  y1 <- c(0,cumsum(tabulate(match(x,unique(x))))/n) # y knots from empirical cdf
  F.knots[1] = 0
  } else {
    y <- unique(x)
    y1 <- cumsum(tabulate(match(x,unique(x))))/n
    F.knots[1] = y1[1]
  }

  i = 1
  x.knots <- array(0) #Array for x.knots of lcm of ecdf

  f.knots <- array(0) #Array for density estimate
  x.knots[1] = y[1]

  count=2
  while(i != length(y))
  {
    slope <- array(0)
    for(j in (i+1):length(y)) #Computing all possible slopes
    {
      slope[j-i] <- ((y1[j]-y1[i])/((y[j]-y[i])))
    }

    for(j in (i+1):length(y))
    {
      if(slope[j-i] == max(slope)) #checking which slope is maximum
      {
        x.knots[count] = y[j]
        F.knots[count] = y1[j]
        f.knots[count-1] = max(slope)
        count = count+1
        i = j
        break
      }
    }
  }
  return(list(sample.values=x,x.knots=x.knots,F.knots=F.knots,f.knots=f.knots))
  }


#-----------------------------------------------------------------------------------------------------------------------

# param is F^(-1)(0) for the decreasing density and F^(-1)(1) for the other part
# a is a fraction that is the data fraction, equals to 1 in general if whole data is used
my_grenander <- function(x, type, param,a)
{

  if (type=="decreasing")
  {
    res <- Grenander_general(x,param)
    res$f.knots <- res$f.knots*a
    res$F.knots <- res$F.knots*a
    return(res)
  }

   res <- Grenander_general(sort(-x),-param)
    res$sample.values <- -res$sample.values
    n <- length( res$x.knots)
    turn <- n-c(0,1:(n-1))
    res$x.knots <- -res$x.knots[turn]
    res$f.knots <- res$f.knots[turn][-1]*a
    res$F.knots <- (1-res$F.knots[turn])*a
    res
}

# Interpolating the grenander function to get all the values of x
# x should be sorted

grenander.inter <- function(x,knot,F.knot,f.knot)
{
  xx <-x
  #print(knot)
  n <- length(x)
  m <- length(knot)
  bg <- length(x[x<knot[1]])
   if(bg==n)
    return(rep(F.knot[1],n))
  # cases where all the x's are less than knots
  en <- length(x[x>=knot[m]])
  if(en==n)
    return(rep(F.knot[m],n))         # cases where all the x's are greater than knots

  x <- x[(bg+1):(n-en)]
 beg <- rep(F.knot[1],bg) # Taking care of x-i's outside support
 end <- rep(F.knot[m],en)
 if(length(x)==0)
   return(c(beg,end))

  count=0
  Fhat <- 1:length(x)
  fhat <- Fhat


  for(i in 1:length(x))
  {    if(x[i]==knot[1])
  {
    count <- 1
    Fhat[i] <- F.knot[count]
    fhat[i] <- f.knot[count]
  } else {
    #print(c(en,x[i],count,knot[count+1],knot))
    #print(c(i,count+1,length(knot)))
    while(x[i]> knot[count+1])
      count <- count+1

      Fhat[i] <- F.knot[count]+(x[i]-knot[count])*f.knot[count]
      fhat[i] <- f.knot[count]
  }
  }

  c(beg,Fhat,end)
  }

####################################################################
# This function gives fhat(x-)
grenander.inter.density <- function(x,knot,F.knot,f.knot)
{
  xx <-x

  n <- length(x)
  m <- length(knot)
  bg <- length(x[x<knot[1]])
  if(bg==n)
    return(rep(0,n))         # cases where all the x's are less than knots
  en <- length(x[x>knot[m]])
  if(en==n)
    return(rep(0,n))         # cases where all the x's are greater than knots

  x <- x[(bg+1):(n-en)]
  beg <- rep(0, bg) # Taking care of x-i's outside support
  end <- rep(0, en)
  if(length(x)==0)
    return(c(beg,end))

  count=0
  Fhat <- 1:length(x)
  fhat <- Fhat


  for(i in 1:length(x))
      {
   # print(c(i,x,knot[1]))
    if(x[i]==knot[1])
  {
    count <- 1
    Fhat[i] <- F.knot[count]
    fhat[i] <- 0  }
 else {
    #print(c(en,x[i],count,knot[count+1],knot))
    #print(c(length(knot),count+1))
    while(x[i]> knot[count+1])
      count <- count+1

    Fhat[i] <- F.knot[count]+(x[i]-knot[count])*f.knot[count]
    fhat[i] <- f.knot[count]

  }
  }

  c(beg,fhat,end)
}
