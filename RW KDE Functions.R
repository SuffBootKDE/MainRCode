#### R( . ) ####
# Calculate squared integration of function #
r.int <- function(f, ..., lower = -Inf, upper = Inf){
  require(pracma)
  f2 <- function(x,...) f(x,...)^2
  
  integrate(f2,..., lower, upper,subdivisions = 1000)$value
}


#### Calculate the R.int for known SkewNorm distn ####
snorm.r.int <- function(x, omega =1, alpha = 0){
  require(sn)
  
  n <- length(x)
  
  #CDF
  DF <- function(x) psn(x, omega = omega, alpha = alpha)
  
  #pdf
  f <- function(x) 2/omega*dnorm(x/omega)*pnorm(alpha*x/omega)
  
  #first derivative
  fp <- function(x) 2/omega^2*((-x/omega)*dnorm(x/omega)*pnorm(alpha*x/omega)
                               + alpha*dnorm(x/omega)*dnorm(alpha*x/omega))
  
  #second derivative
  fpp <- function(x) 2/omega^3*dnorm(x/omega)*((x^2/omega^2 - 1)*pnorm(alpha*x/omega)  
                                               - alpha*(2 - alpha^2)*x/omega*dnorm(alpha*x/omega))
  
  #inner part of term of R
  gx <- function(x){
    fpp(x) + (n-1)*(fpp(x)*DF(x) + 3*fp(x)*f(x)) 
    }
  
  
  return(r.int(gx))
  
  
}

#### Pseudo Rule of Thumb/Plugin ####
# Takes data, calculates MLE estimates of alpha and omega in Azzalini skewnorm.
# Calculates R(G(x)) for the skew norm with estimated parameters
# Uses that to calculate h_AMISE for the assumed skewnormal distribution
rw.plugin2 <- function(x){
  
  n <- length(x)
  RK <- 1/(2*sqrt(pi))
  sdK <- 1
  #Rf <- integrate(r.inner(x), min(x)-3*sd(x), max(x)+3*sd(x), subdivisions = 1000)$value
  
  # Estimate mean, sd, and skew of sample.
  fit <- sn.mple(y = x)
  
  # Convert values estimated to direct SN paramters and store.
  parm <- cp2dp(fit$cp, family = "SN")
  
  xi <- parm[1]; omega <- parm[2]; alpha <- parm[3]
  
  Rf <- snorm.r.int(x, omega = omega, alpha = alpha)
  
  ( ( RK*sum(1/(sqrt(1:n))) )/( 2*n^2*sdK^4*Rf )  )^(1/5)
}

#### Classic Rule of THumb Bandwidth ####
# s * ( 4 / (3 * n) )^(1/5)
null.bw0 <- function(x){
  n <- length(x)
  s <- sd(x)
  bw <- s * ( 4 / (3 * n) )^(1/5)
  return(bw)
}

#### h_opt calculation functions ####
### Function calculates the integral term of numerator in h_opt for ranks
kofn <- function(n) 1/n*sum(1/sqrt(1:n))


#### Rank Weighted Rule of Thumb BW 

#### Rank Weighted Rule of Thumb BW 
rw.bw0. <- function(x, kernel = "gaussian"){
  n <- length(x)
  s <- sd(x)
  
  if(kernel == "epanechnikov"){ rk <- 3/5 ; sk <- 1/5}
  else if(kernel == "gaussian"){ rk <- 1/(2*sqrt(pi)); sk <- 1}
  else {stop("Invalid kernel specification.")}
  
  # Special gaussian integrals. used for rule of thumb
  #I1 <- r.int(function(x) {(x^2 - 1)*dnorm(x)})$value
  I1 <- 0.2115711
  #I2 <- r.int(function(x){ (x^2 - 1)*dnorm(x)*pnorm(x)-3*x*dnorm(x)^2} )$value
  I2 <- 0.1454667
  #I3 <- integrate(function(x) {(x^2-1)*dnorm(x)*((x^2-1)*dnorm(x)*pnorm(x)
  #                                               -3*x*dnorm(x)^2)}, 
  #                lower = -Inf, upper = Inf)$value
  I3 <- 0.1057855
  
  bw = s*( ( rk*sum(1/sqrt(1:n))) /
             (2 * n^2 * sk * (I1 + (n-1)^2 * I2 + 2*(n-1)*I3) ) )^(1/5)
  return(bw)
}

rw.bw0 <- rw.bw0.


# Replicate of rank.lscv2
#### LSCV Rankweighted Bandwidth (Prototype) ####
# Functions takes data and calculates RankLSCV bandwidth.
rw.lscv<- function(x, lower  = 1E-6, upper = NULL, tol = 1E-10, kernel = dnorm){
  if(is.null(upper)){ 
    upper = 1.144 * sqrt(var(x)) * n^(-1/5)
  }
  
  # Sample size
  n <- length(x)
  
  # Get all pairwise indices i < j
  pair.ind <- combn(1:n,2)
  i <- pair.ind[1,] ; j <- pair.ind[2,]
  
  # Get all differences, i < j
  diffs <- x[i] - x[j]
  
  # constants from expression
  kk1 <- 1/n^2 * 1/sqrt(2*pi) * sum(1/sqrt(2*(1:n)))
  nn2 <- 2/n^2
  
  # create LSCV h function
  lscv <- function(h){
    
    kk1/h + nn2*sum(1/(h*sqrt(i+j)) * dnorm(diffs/(h*sqrt(i+j))) ) - (4/(n*(n-1)*h))*sum(1/sqrt(j)*dnorm(diffs/(h*sqrt(j))))
    
  } 
  h <- optimize(lscv, c(lower, upper), tol = tol)$minimum
  return(h)
}

#### Transformation Estimator ####
#### LSCV Rankweighted Transformation Bandwidth ####
# Trans is the transformation function being used
rw.trans.lscv3 <- function(x, trans = log, transdx = function(x) 1/x,
                     lower  = 1E-6, upper = NULL, tol = 1E-10, kernel = dnorm){
  
  #transformed values
  y <- transdx(x)
  
  if(is.null(upper)){ 
    upper = 1.144 * sqrt(var(x)) * n^(-1/5)
  }
  
  # Sample size
  n <- length(x)
  
  # Get all pairwise indices i < j
  pair.ind <- combn(1:n,2)
  i <- pair.ind[1,] ; j <- pair.ind[2,]
  # Get all differences
  diffs <- x[i] - x[j]
  
  
  # All possible pairs for use in first term
  i. <- rep(1:n, n); j. <- rep(1:n, each = n)
  k. <- sqrt(i.*j.)
  w. <- (x[i.] - x[j.])/k.
  
  
  # create LSCV h function
  lscv <- function(h){
    
    
    1/n^2*1/h*sum(1/k.*exp(-(w./h)^2/4 )/(sqrt(4*pi) ) ) - (4/(n*(n-1)*h))*sum(1/sqrt(j)*dnorm(diffs/(h*sqrt(j))))
  } 
  h <- optimize(lscv, c(lower, upper), tol = tol)$minimum
  return(h)
}


#### PLUG IN ESTIMATOR STUFF ####

### Function for estimating derivatives ###
# fp.hat returns a function that computes derivate values
# requires
fp.hat <- function(x, order = 0){
  require(kedd)
  
  r <- order
  
  # Get bandwidth for derivative order.
  bw <- h.ucv(x, deriv.order = order)$h
  
  
  
  fphat <- function(x0, h = bw){
    
    m <- length(x0)
    n <- length(x)
    
    #Get kernel estimates for each point in x
    y <- rowSums(matrix((-1)^r*kernel.fun((rep(x, times = m) - rep(x0, each = n))/h, 
                                          deriv.order = order, kernel = "gaussian")$kx/(n*h^(order + 1)),
                        nrow = m, ncol = n, byrow = T))
    return(y)
    
  }
  return(fphat)
  
}

###R Function for R(f) = int(f^2, -inf, inf) ###
# Function for helping in computing R term. Integrate to get actual R term
# This function just computes f(x)^2
# Much faster.
r.inner <- function(x){
  
  n <- length(x)
  
  # Create density order estimating functions.
  tempcdf <- ecdf(x) # CDF
  f <- fp.hat(x, order = 0) # Regular Density
  fp <- fp.hat(x, order = 1) # First deriv
  fpp <- fp.hat(x, order = 2) # Second Deriv
  
  
  R <- function(y){
    ( fpp(y) + 
        (n-1)*(fpp(y)*tempcdf(y) + 3*f(y)*fp(y)) )^2
  }
  
  return(R)
}


### Plugin Type Bandwidth for rankweighted Estimator ###
# Computes plug in bandwidth for rank weight kde by estimating h_amise with
# kdes for functional of F in the denominator of h_amise: R(f'' + (n-1)(f''*F +3*f*f'))
# the functional is calculated using LSCV kdes in the R function
# Needs modification for anything besides gaussian Kernel w/ sd of 1
rw.plugin1 <- function(x){
  
  n <- length(x)
  RK <- 1/(2*sqrt(pi))
  sdK <- 1
  #Rf <- integrate(r.inner(x), min(x)-3*sd(x), max(x)+3*sd(x), subdivisions = 1000)$value
  Rf <- integral(r.inner(x), -Inf, Inf)
  
  ( ( RK*sum(1/(sqrt(1:n))) )/( 2*n^2*sdK^4*Rf )  )^(1/5)
}

#### END PLUG IN ESTIMATOR STUFF ####

#### Rank Weighted Density Calculator ####
# x is data, bw = bandwidth, sk = standard deviation of Kernel, nref = number of points to estimate
# at, xref = desired points to estimate at. 
# K is the desired kernel function
rw.density <- function(x, bw = "plugin",  sk = 1, nref = 512, xref = NA, y.only = F, 
                       kernel = "gaussian"){
  require(gplm)
  
  x <- sort(x)
  n <- length(x)
  r <- 1:n # Ranks since x is now sorted
  s <- sd(x)
  
  #set Kernel function
  K <- function(x) kernel.function2(x, kernel = kernel)
  
  # These are for getting values above and below the range of the data.
  if(all(is.na(xref))) xref <- seq(min(x) - s, max(x) + s, length.out = nref)
  
  # 
  nref <- length(xref)
  
  # Default bw
  if(bw == "rot"){
    
    bw <- rw.bw0(x, kernel = kernel )
    
  }else if(bw=="plugin"){
  
  bw <- rw.plugin1(x)
  
  }else if(bw == "lscv"){
    
    bw <- rw.lscv(x, kernel = kernel)
    
  }else if(bw == "plugin2"){
    
    bw <- rw.plugin2(x)
    
  }else if(is.numeric){}
  
  # Had an else stop() here, maybe it shouldn't be here?
  
  dens.temp <- function(t) (1/n) * (1/bw) * sum( (1/sqrt(r)) * K( (t - x)/(bw*sqrt(r))))
  yref <- as.vector(apply(as.matrix(xref), 1, dens.temp))
  
  if(!y.only) return(data.frame(xref, yref))
  if(y.only) return(yref)
}



#### Get rankweighted bandwidth ONLY ####

rw.bw <- function(x, bw = "rot", kernel = "gaussian"){
  require(gplm)
  
  x <- sort(x)
  n <- length(x)
  r <- 1:n # Ranks since x is now sorted
  s <- sd(x)
  
  #set Kernel function
  K <- function(x) kernel.function2(x, kernel = kernel)
  
  # Default bw
  if(bw == "rot"){
    
    bw <- rw.bw0(x, kernel = kernel )
    
  }else if(bw=="plugin"){
    
    bw <- rw.plugin1(x)
    
  }else if(bw == "lscv"){
    
    bw <- rw.lscv(x, kernel = kernel)
    
  }else if(bw == "plugin2"){
    
    bw <- rw.plugin2(x)
    
  }else if(is.numeric(bw)){}
  else stop("Bad bw")
  
  return(bw)

}

#### Equivalent function for calculating Std Bws. Equivalent to rw.bws directly above ####
std.bw <- function(x, bw = "nrd0"){
  require(gplm)
  
  x <- sort(x)
  n <- length(x)
  r <- 1:n # Ranks since x is now sorted
  s <- sd(x)
  
  #set Kernel function
  #K <- function(x) kernel.function2(x, kernel = kernel)
  
  # Default bw
  if(bw == "nrd0"){
    
    bw <- bw.nrd0(x)
    
  }else if(bw=="SJ"){
    
    bw <- bw.SJ(x)
    
  }else if(bw == "ucv"){
    
    bw <- bw.ucv(x)
    
  }else if(is.numeric(bw)){}
  
   else stop("Bad bw")
  
  return(bw)
  
}
#### Rank Weighted Density Function ####
# Returns function that can calculate fhat at a given vector x
# x is data, bw = bandwidth, sk = standard deviation of Kernel, nref = number of points to estimate
# at, xref = desired points to estimate at. 
# K is the desired kernel function
rw.fhat <- function(x, h, kernel = "gaussian"){
  require(gplm)
  
  #x <- sort(x)
  n <- length(x)
  #r <- 1:n # Ranks since x is now sorted
  r <- rank(x)
  
  #set Kernel function
  K <- function(x) dnorm(x)
  
  # These are for getting values above and below the range of the data.
  fhat <-function(x0){
    k = length(x0)
    x. = matrix(rep(x, k), k, n, byrow = T); r. <- matrix(rep(r, k), k, n, byrow = T)
    sqrt.r. <- sqrt(r.)
    x0. = matrix(rep(x0, each = n), k, n, byrow = T)
    z = (x0.-x.)/(sqrt.r.*h)
    W = (1/n)*K(z)*1/(h*sqrt.r.)
    return(rowSums(W))
  }
  return(fhat)
}

#### LSCV fhat that takes h as arugment ####
# x is data, bw = bandwidth, sk = standard deviation of Kernel, nref = number of points to estimate
# at, xref = desired points to estimate at. 
# K is the desired kernel function
rw.fhath <- function(x, kernel = "gaussian"){
  require(gplm)
  
  x <- sort(x)
  n <- length(x)
  r <- 1:n # Ranks since x is now sorted
  
  #set Kernel function
  K <- function(x) kernel.function2(x, kernel = kernel)
  
  # These are for getting values above and below the range of the data.
  fhat <-function(h){
    x0 = x
    k = length(x0)
    x. = matrix(rep(x, k), k, n, byrow = T); r. <- matrix(rep(r, k), k, n, byrow = T)
    sqrt.r. <- sqrt(r.)
    x0. = matrix(rep(x0, each = n), k, n, byrow = T)
    W = K((x0.-x)/(sqrt.r.*h))*1/(h*sqrt.r.)
    return((1/n)*rowSums(W))
  }
  return(fhat)
}

#### Function for getting just y values out of density() ####
# suply number of points, from, and to.
density.yonly <- function(x, bw = "nrd0", kernel = "gaussian", n, from, to){
  density(x, kernel = kernel, bw = bw, n = n, from = from, to = to)$y
}



#### 2 Normal Mixture generation and Density Calculation ####
# Normal Mixture 
# Mixture probability of p
# p is for first norm. Function creates two.
rnorm.2mix <- function(n,mu, sd, p){
  if(length(mu) + length(sd) !=4) return(cat("mu and sd need to be vectors of size two"))
  else{
    q <- rbinom(n,1,p)
    x <- q*rnorm(n,mu[1],sd[1]) + (1-q)*rnorm(n,mu[2],sd[2])
    return(x)
  }
}

dnorm.2mix <- function(x, mu, sd, p){
  if(length(mu) + length(sd) !=4) return(cat("mu and sd need to be vectors of size two"))
  else{
    y <- p*dnorm(x,mu[1],sd[1]) + (1-p)*dnorm(x,mu[2],sd[2])
    return(y)
  }
}





#### norm2.mix.sim ####
# This function creates random data from a normal mixture.
# Plots a histogram, true density (black dash), RankWeight Density (red), and classic density (blue)
# Will also simulate from single normal if mix = F
# ref is for how small the step size is when making x.ref vector which
# is for plotting true density. 
norm2.mix.sim <- function(n, mu, sd, p, mix = T, ref = .01){
  # DATA
  if(mix){
    x <- sort(rnorm.2mix(n,mu, sd , p))
    
    #x.ref is for plotting true density. Important at small sample sizes.
    x.ref <- seq(min(x), max(x), ref)
    y <- dnorm.2mix(x.ref, mu, sd, p)
  }
  
  else if(!mix){
    x <- sort(rnorm(n,mu, sd))
    
    x.ref <- seq(min(x), max(x), ref)
    y <- dnorm(x.ref, mu, sd)
  }
  # Classic Density
  d.null.0 <- density(x, bw = null.bw0(x))
  
  # RW density
  d.rw.0 <- rw.density(x)
  
  #Plot it
  hist(x, freq = F, ylim =c(0,max(y, d.null.0$y, d.rw.0$y) ),
       main = paste("n = ",n, ", p = ",p, ", mu = ", paste(mu, collapse = ", "), ", sd = ", 
                    paste(sd, collapse = ", ")))
  
  lines(x.ref, y, lty = 2)
  lines(d.rw.0, col = "red")
  lines(d.null.0$x, d.null.0$y, col = "blue")
}

#### Generate k normal  mixtures ####
# p is mixture probabilities. It could also be relative size of parts
# p is normalized so that it adds to 1. p = c(1, 2, 2) -> c(.2, .4, .4)
# must provide vectors of mu and sd

rnorm.mix <- function(n, mu, sd, p){
  # number of mixtures
  k <- length(p)
  
  p <- p/sum(p)
  
  # data vector and density vector.
  x <- rep(NA, n) 
  y <- rep(0, n)
  
  # Multinomial variable that creates mixture indicators
  mix <- rmultinom(n, 1, p) == 1
  
  for(i in 1:k){
    x[mix[i, ]] <- rnorm(sum(mix[i,]), mu[i], sd[i])
  }
  
  return(x)
}



dnorm.mix <- function(x, mu, sd, p){
  p <- p/sum(p)
  k <- length(p)
  y <- rep(0, length(x))
  
  for(i in 1:k){
    y <- y + p[i]*dnorm(x, mu[i], sd[i])
  }
  return(y)
  
}






# Transformation STUFF ----------------------------------------------------

#### Transformed Rank Weighted Density Calculator ####
# x is data, bw = bandwidth, sk = standard deviation of Kernel, nref = number of points to estimate
# at, xref = desired points to estimate at. 
# K is the desired kernel function
trw.density <- function(x, bw = "plugin", tt = log, dtt = function(x) 1/x,  
                        sk = 1, nref = 512, xref = NA, y.only = F, kernel = "gaussian"){
  require(gplm)
  
  
  
  x <- sort(x)
  n <- length(x)
  r <- 1:n # Ranks since x is now sorted
  s <- sd(x)
  
  #set Kernel function
  K <- function(x) kernel.function2(x, kernel = kernel)
  
  # These are for getting values above and below the range of the data.
  if(all(is.na(xref))) xref <- seq(max(c(min(x) - s, 0.0001)), max(x) + s, length.out = nref)
  
  # 
  nref <- length(xref)
  
  if(bw == "lscv"){
    
    bw <- trw.lscv(x)
    
  }else if(is.numeric(bw)){}
  
  # Had an else stop() here, maybe it shouldn't be here?
  
  # transform data
  ttx = tt(x)
  
  
  dens.temp <- function(t) {
    
    
    (abs(dtt(t))/n) * (1/bw) * sum( (1/sqrt(r)) * K( (tt(t) - tt(x))/(bw*sqrt(r))))
  }
  yref <- as.vector(apply(as.matrix(xref), 1, dens.temp))
  
  if(!y.only) return(data.frame(xref, yref))
  if(y.only) return(yref)
}


#### Transformed KDE Density Calculator ####
# x is data, bw = bandwidth, sk = standard deviation of Kernel, nref = number of points to estimate
# at, xref = desired points to estimate at. 
# K is the desired kernel function
t.density <- function(x, bw = "lscv", tt = log, dtt = function(x) 1/x,  
                      sk = 1, nref = 512, xref = NA, y.only = F, kernel = "gaussian"){
  require(gplm)
  
  
  
  x <- sort(x)
  n <- length(x)
  s <- sd(x)
  
  #set Kernel function
  K <- function(x) kernel.function2(x, kernel = kernel)
  
  # These are for getting values above and below the range of the data.
  if(all(is.na(xref))) xref <- seq(max(c(min(x) - s, 0.0001)), max(x) + s, length.out = nref)
  
  # 
  nref <- length(xref)
  
  if(bw == "lscv"){
    
    bw <- t.lscv(x)
    
  }else if(is.numeric(bw)){}
  
  # Had an else stop() here, maybe it shouldn't be here?
  
  # transform data
  ttx = tt(x)
  
  
  dens.temp <- function(t) {
    
    
    (abs(dtt(t))/n) * (1/bw) * sum( K( (tt(t) - tt(x))/(bw)))
  }
  yref <- as.vector(apply(as.matrix(xref), 1, dens.temp))
  
  if(!y.only) return(data.frame(xref, yref))
  if(y.only) return(yref)
  
}

#### Transformation LSCV Rankweighted  Bandwidth ####
# Trans is the transformation function being used
trw.lscv <- function(x, tt = log, tdx = function(x) 1/x,
                     lower  = 1E-6, upper = NULL, tol = 1E-10, kernel = dnorm){
  
  #save regular values
  x.bk <- x
  
  #transformed values
  x<- tt(x)
  dx <- tdx(x)
  
  n <- length(x)
  
  if(is.null(upper)){ 
    upper = 1.144 * sqrt(var(x)) * n^(-1/5)
  }
  
  # Sample size
  n <- length(x)
  
  # Get all pairwise indices i < j
  pair.ind <- combn(1:n,2)
  i <- pair.ind[1,] ; j <- pair.ind[2,]
  root.j <- sqrt(j)
  
  # Get all differences
  diffs <- x[i] - x[j]
  
  #other vales
  over.xi.root.j <- 1/x.bk[i]/root.j
  
  # All possible pairs for use in first term
  i. <- rep(1:n, n); j. <- rep(1:n, each = n)
  c. <- j.*x[i.] + i.*x[j.]; d. <- j.*x[i.]^2 + i.*x[j.]^2
  ij. <- i.*j.; root.ipj <- sqrt(i. + j.); i.root.ipj <- 1/root.ipj
  root.ij <- sqrt(ij.)
  
  
  
  # create LSCV h function
  lscv <- function(h){
    
    (1/n^2)*(1/h)*sum(  i.root.ipj * dnorm( (d. - ( c. - ij.*h^2 )*i.root.ipj )/( root.ij*h ) ) ) - 
      4/(h*n*(n-1))*sum( over.xi.root.j * dnorm( (x[i] - x[j])/(root.j*h) )) 
    
  } 
  h <- optimize(lscv, c(lower, upper), tol = tol)$minimum
  return(h)
  
  
  
}

#### Transformation LSCV Bandwidth ####
# Trans is the transformation function being used
t.lscv <- function(x, tt = log, tdx = function(x) 1/x,
                   lower  = 1E-6, upper = NULL, tol = 1E-10, kernel = dnorm){
  
  #save regular values
  x.bk <- x
  
  #transformed values
  x<- tt(x)
  dx <- tdx(x)
  
  n <- length(x)
  
  if(is.null(upper)){ 
    upper = 1.144 * sqrt(var(x)) * n^(-1/5)
  }
  
  # Sample size
  n <- length(x)
  
  # Get all pairwise indices i < j
  pair.ind <- combn(1:n,2)
  i <- pair.ind[1,] ; j <- pair.ind[2,]
  
  
  # Get all differences
  diffs <- x[i] - x[j]
  
  
  # All possible pairs for use in first term
  i. <- rep(1:n, n); j. <- rep(1:n, each = n)
  
  # some constant stuff 
  d. <- x[i.] + x[j.]; c. <- x[i.]^2 + x[j.]^2
  
  
  
  
  # create LSCV h function
  lscv <- function(h){
    
    (1/(2*sqrt(pi)*h*n^2))*sum( exp( (1/(4*h^2))*(d. - h^2)^2 - c./(2*h^2) ) ) -
      4/(n*(n-1)*h)*sum( 1/x.bk[i]*dnorm( (diffs)/h ) )
    
  } 
  h <- optimize(lscv, c(lower, upper), tol = tol)$minimum
  return(h)
  
  
  
}



#### plotting/simulation function for rule of thumb ####
# rsample <- RNG function for desired distribution
# dsamp <- density function for desired distribution
# ... parameters for distribution 
# lower and upper are bounds for plot of true density.
plot.density.sim0 <- function(n, rsamp, dsamp, lower = NA, upper = NA,...){
  test <- rsamp(n, ...)
  
  # Set lower and upper if undefined
  if(is.na(lower)) lower = min(test); if(is.na(upper)) upper = max(test)
  
  # Reference values for plotting true density.
  x <- seq(lower, upper, .001)
  y <- dsamp(x, ...)
  
  # Classic Density Estimate
  null.test.density <- density(test, bw = null.bw0(test))
  
  # PLOTS
  plot(x,y, type = "l")
  lines(rw.density(test), col = "blue", lty = 2)
  lines(rw.density(test, bw = bw.ucv(test)), col = "green", lty = 2)
  lines(null.test.density$x, null.test.density$y, col = "red", lty = 2)
}



#### Create Var, Bias, MSE at each x values ####
# Yhat is data matrix containing several density estimate values at given xref points
# yref is the true density values at xref points. Needed for bias.
# Each row is for a single value of xref. So if there are 100 columns, this implies
# there are 100 different density estimates at a given point. 
# Function returns the variance at each xref, the average bias, and the estimated MSE
density.stats <- function(fhat, xref, yref){
  varf <- apply(fhat, 1, var)
  meanfhat <- apply(fhat, 1, mean)
  biasf <- meanfhat - yref
  msef <- varf + biasf ^2
  
  return(data.frame(xref, yref, meanfhat, varf, biasf, msef))
}




#### General Simulations for Bandwidths Function ####
# M = Simulation Size, n = sample size(s), rsamp and dsamp are sample and ensity function for
# desired distirbution of simulations, lower, upper are the bounds for where the reference x points
# are to be made. nref is how many points to estimate the density at.
# ... are the distribution parameters for rsamp and dsamp
# rw.bws and std.bws are vectors of strings that list bandwidth selection methods to be used.
# if there are no specified bandwidth methods, rule of thumb procedures are automatically applied
#CHECK
bw.density.sim <- function(M, n, rsamp, dsamp, kernel = "gaussian", rw.bws, std.bws, lower, upper, nref = 512,
                           debug = T, ...){
  distsamples <- list()
  for(i in 1:length(n)){
    distsamples[[i]] <- matrix(rsamp(n[i]*M, ...), nrow = n[i], ncol = M)
  }
  
  # Number of reference points for density estimate
  xref <- seq(lower, upper, length.out = nref)
  yref <- dsamp(xref,...)
  
  # Density estimates at each value of xref for each sample
  # rank weighted rule of thumb
  disthats <- list()
  disthats$rw <- list()
  
    # standard case rule of thumb
  disthats$std <- list()
  
  # Store h! # HARD TO DO! NEED TO REFORMAT EVERYTHING

  
  # Calculate density estimates for rank weighted and standard estimators
  # Do so for each bandwidth type specified [j]
  # And at each sample size [i]
  for(j in 1:length(rw.bws)){
    
    disthats$rw[[j]] <- list()
    
    for(i in 1:length(n)){
      
      disthats$rw[[j]][[i]] <- apply(distsamples[[i]], MARGIN = 2,
                                rw.density, bw = rw.bws[j], xref = xref, nref = nref,
                                kernel = kernel, y.only = T)      
    }
  }
  
  if(debug) print(" Density Values Calculated:RW")
    
  for(j in 1:length(std.bws)){
    
    disthats$std[[j]] <- list()
    
    for(i in 1:length(n)){
      
            disthats$std[[j]][[i]] <- apply(distsamples[[i]], MARGIN = 2,
                                 density.yonly, kernel = kernel, bw = std.bws[j],
                                 n = nref, from = lower, to = upper)
    }
  }
  
  if(debug) print(" Density Values Calculated:STD")
  
  # Name respective estimates by bw type selected
  names(disthats$rw) <- rw.bws
  names(disthats$std) <- std.bws
  
  # Get stats for each sample size
  # Create similar list nesting structure as with disthats
  rw <- list()
  std<- list()
  
  #vector samples size names vector
  n.names <- paste("n", n, sep = "")
  
  # Stats calculations for rankweighted case
  for(j in 1:length(rw.bws)){
    rw[[j]] <- list()
    
    for(i in 1:length(n)){
      rw[[j]][[i]] <- density.stats(disthats$rw[[j]][[i]], xref = xref,
                               yref = yref)
    }
    
    names(rw[[j]]) <- n.names
    
  }
  
  if(debug) print(" density stats for RW.")
  
  names(rw) <- rw.bws
  
  #vector samples size names vector
  n.names <- paste("n", n, sep = "")
  
  # Stats calculations for standard case
  for(j in 1:length(std.bws)){
    std[[j]] <- list()
    
    for(i in 1:length(n)){
      std[[j]][[i]] <- density.stats(disthats$std[[j]][[i]], xref = xref,
                                    yref = yref)
    }
    
    names(std[[j]]) <- n.names
    
  }
  
  # assign bw type names
  names(std) <- std.bws #### NEED TO DOUBLE CHECK UP TO THIS POINT ####
  #### THEN MODIFY DENSITY CALCULATION FUNCTIONS TO ACCEPT BANDWIDTH SELECTIONS ####

  
  return(list(rw = rw, std = std))
}

#### Tranformed KDE Simulation ####
# M = Simulation Size, n = sample size(s), rsamp and dsamp are sample and ensity function for
# desired distirbution of simulations, lower, upper are the bounds for where the reference x points
# are to be made. nref is how many points to estimate the density at.
# ... are the distribution parameters for rsamp and dsamp
# rw.bws and std.bws are vectors of strings that list bandwidth selection methods to be used.
# if there are no specified bandwidth methods, rule of thumb procedures are automatically applied
#CHECK
t.bw.density.sim <- function(M, n, rsamp, dsamp, kernel = "gaussian", rw.bws = "lscv", std.bws = "lscv", 
                             lower, upper, nref = 512, debug = T, ...){
  distsamples <- list()
  for(i in 1:length(n)){
    distsamples[[i]] <- matrix(rsamp(n[i]*M, ...), nrow = n[i], ncol = M)
  }
  
  # Number of reference points for density estimate
  xref <- seq(lower, upper, length.out = nref)
  yref <- dsamp(xref,...)
  
  # Density estimates at each value of xref for each sample
  # rank weighted rule of thumb
  disthats <- list()
  disthats$rw <- list()
  
  # standard case rule of thumb
  disthats$std <- list()
  
  # Store h! # HARD TO DO! NEED TO REFORMAT EVERYTHING
  
  
  # Calculate density estimates for rank weighted and standard estimators
  # Do so for each bandwidth type specified [j]
  # And at each sample size [i]
  for(j in 1:length(rw.bws)){
    
    disthats$rw[[j]] <- list()
    
    for(i in 1:length(n)){
      
      disthats$rw[[j]][[i]] <- apply(distsamples[[i]], MARGIN = 2,
                                     trw.density, bw = rw.bws[j], xref = xref, nref = nref,
                                     kernel = kernel, y.only = T)      
    }
  }
  
  if(debug) print(" Density Values Calculated:RW")
  
  for(j in 1:length(std.bws)){
    
    disthats$std[[j]] <- list()
    
    for(i in 1:length(n)){
      
      disthats$std[[j]][[i]] <- apply(distsamples[[i]], MARGIN = 2,
                                      t.density, bw = std.bws[j], xref = xref, nref = nref,
                                      kernel = kernel, y.only = T)   
    }
  }
  
  if(debug) print(" Density Values Calculated:STD")
  
  # Name respective estimates by bw type selected
  names(disthats$rw) <- rw.bws
  names(disthats$std) <- std.bws
  
  # Get stats for each sample size
  # Create similar list nesting structure as with disthats
  rw <- list()
  std<- list()
  
  #vector samples size names vector
  n.names <- paste("n", n, sep = "")
  
  # Stats calculations for rankweighted case
  for(j in 1:length(rw.bws)){
    rw[[j]] <- list()
    
    for(i in 1:length(n)){
      rw[[j]][[i]] <- density.stats(disthats$rw[[j]][[i]], xref = xref,
                                    yref = yref)
    }
    
    names(rw[[j]]) <- n.names
    
  }
  
  if(debug) print(" density stats for RW.")
  
  names(rw) <- rw.bws
  
  #vector samples size names vector
  n.names <- paste("n", n, sep = "")
  
  # Stats calculations for standard case
  for(j in 1:length(std.bws)){
    std[[j]] <- list()
    
    for(i in 1:length(n)){
      std[[j]][[i]] <- density.stats(disthats$std[[j]][[i]], xref = xref,
                                     yref = yref)
    }
    
    names(std[[j]]) <- n.names
    
  }
  
  # assign bw type names
  names(std) <- std.bws #### NEED TO DOUBLE CHECK UP TO THIS POINT ####
  #### THEN MODIFY DENSITY CALCULATION FUNCTIONS TO ACCEPT BANDWIDTH SELECTIONS ####
  
  
  return(list(rw = rw, std = std))
}


#### BW Only Simulations####
# M = Simulation Size, n = sample size(s), rsamp and dsamp are sample and ensity function for
# desired distirbution of simulations, 
# ... are the distribution parameters for rsamp and dsamp
# rw.bws and std.bws are vectors of strings that list bandwidth selection methods to be used.
# if there are no specified bandwidth methods, rule of thumb procedures are automatically applied
#CHECK
bw.density.sim_bwonly <- function(M, n, rsamp, dsamp, kernel = "gaussian", rw.bws, std.bws, ...){
  distsamples <- list()
  for(i in 1:length(n)){
    distsamples[[i]] <- matrix(rsamp(n[i]*M, ...), nrow = n[i], ncol = M)
  }
  

  
  # Density estimates at each value of xref for each sample
  # rank weighted rule of thumb
  h <- list()
  h$rw <- list()
  
  # standard case rule of thumb
  h$std <- list()
  
  # Store h! # HARD TO DO! NEED TO REFORMAT EVERYTHING
  
  #vector samples size names vector
  n.names <- paste("n", n, sep = "")
  
  time <- proc.time()[3]
  # Calculate density estimates for rank weighted and standard estimators
  # Do so for each bandwidth type specified [j]
  # And at each sample size [i]
  for(j in 1:length(rw.bws)){
    
    h$rw[[j]] <- data.frame(rep(NA, M))
    
    for(i in 1:length(n)){
      
      h$rw[[j]][,i] <- apply(distsamples[[i]], MARGIN = 2,
                                     rw.bw, bw = rw.bws[j])
    }
    names(h$rw[[j]]) <- n.names
  }
  
  print(paste("Density Values Calculated:RW, Time = ", round(proc.time()[3] - time, 3)) )
  
  time <- proc.time()[3]
  for(j in 1:length(std.bws)){
    
    h$std[[j]] <- list()
    
    for(i in 1:length(n)){
      
      h$std[[j]][[i]] <- apply(distsamples[[i]], MARGIN = 2,
                                      std.bw, bw = std.bws[j])
    }
  }
  
  print(paste("Density Values Calculated:STD", round(proc.time()[3] - time, 3) ))
  

  
  return(h)
}

#### Multiplot Function ####
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#### Create plots for comparing bias, var, and MSE. ####
# Function takes in the list for a set of simulations done for given distribution.
# must specify that n layer. That is, specify the index of the 'n' values used
# if n = c(10, 50, 300) and we want graphs for size 50 sample sizes, use nind = 2.
# 'distribution' is a distribution name 
rot.density.plots <- function(X,nind, distribution = "", abs.bias = T){
  require(ggplot2); require(cowplot)
  # p0 is the g
  xref <- X$rw.rot[[nind]]$xref
  yref <- X$rw.rot[[nind]]$yref
  p0 = qplot(xref, geom = 'blank') +  # Create common set of x values
    geom_line(aes(y = yref, colour = 'True Distribution')) +  # True distribution
    geom_line(aes(y = X$rw.rot[[nind]]$meanfhat, colour = 'Rank Weighted')) + # rank weighted                    
    geom_line(aes(y = X$std.rot[[nind]]$meanfhat, colour = 'Standard')) + # Standard Estimator
    ggtitle(paste("Avg Density Estimates: ", distribution, ", ", names(X$rw.rot)[nind] )) +
    scale_colour_manual(name = 'Density', values = c('blue', 'red', 'black'))
  
  # Bias Plot
  if(!abs.bias) {
    p1 = qplot(xref, geom = 'blank') +  # Create common set of x values
      geom_line(aes(y = X$rw.rot[[nind]]$biasf, colour = 'Rank Weighted')) + # rank weighted                    
      geom_line(aes(y = X$std.rot[[nind]]$biasf, colour = 'Standard')) +
      ggtitle("Bias") +
      scale_colour_manual(name = 'Bias', values = c( 'blue', 'red'))
  }
  else {
    p1 = qplot(xref, geom = 'blank') +  # Create common set of x values
      geom_line(aes(y = abs(X$rw.rot[[nind]]$biasf), colour = 'Rank Weighted')) + # rank weighted                    
      geom_line(aes(y = abs(X$std.rot[[nind]]$biasf), colour = 'Standard')) +
      ggtitle("Absolute Bias") +
      scale_colour_manual(name = 'Abs. Bias', values = c( 'blue', 'red'))
  }
  # Variance Plot
  p2 = qplot(xref, geom = 'blank') +  # Create common set of x values
    geom_line(aes(y = X$rw.rot[[nind]]$varf, colour = 'Rank Weighted')) + # rank weighted                    
    geom_line(aes(y = X$std.rot[[nind]]$varf, colour = 'Standard')) +
    ggtitle("Variance") +
    scale_colour_manual(name = 'Var.', values = c( 'blue', 'red'))
  
  # Variance Plot
  p3 = qplot(xref, geom = 'blank') +  # Create common set of x values
    geom_line(aes(y = X$rw.rot[[nind]]$msef, colour = 'Rank Weighted')) + # rank weighted                    
    geom_line(aes(y = X$std.rot[[nind]]$msef, colour = 'Standard')) +
    ggtitle("Mean Square Error") +
    scale_colour_manual(name = 'MSE', values = c( 'blue', 'red'))
  
  multiplot(p0,p1,p2,p3)
  
}

#### Create plots for comparing bias, var, and MSE. ####
# Function takes in the list for a set of simulations done for given distribution.
# must specify that n layer. That is, specify the index of the 'n' values used
# if n = c(10, 50, 300) and we want graphs for size 50 sample sizes, use nind = 2.
# 'distribution' is a distribution name 
density.plots <- function(X,nind, distribution = "", abs.bias = T){
  require(ggplot2); require(cowplot)
  # p0 is the g
  xref <- X$rw[[nind]]$xref
  yref <- X$rw[[nind]]$yref
  p0 = qplot(xref, geom = 'blank') +  # Create common set of x values
    geom_line(aes(y = yref, colour = 'True Distribution')) +  # True distribution
    geom_line(aes(y = X$rw[[nind]]$meanfhat, colour = 'Rank Weighted')) + # rank weighted                    
    geom_line(aes(y = X$std[[nind]]$meanfhat, colour = 'Standard')) + # Standard Estimator
    ggtitle(paste("Avg Density Estimates: ", distribution, ", ", names(X$rw)[nind] )) +
    scale_colour_manual(name = 'Density', values = c('blue', 'red', 'black'))
  
  # Bias Plot
  if(!abs.bias) {
    p1 = qplot(xref, geom = 'blank') +  # Create common set of x values
      geom_line(aes(y = X$rw[[nind]]$biasf, colour = 'Rank Weighted')) + # rank weighted                    
      geom_line(aes(y = X$std[[nind]]$biasf, colour = 'Standard')) +
      ggtitle("Bias") +
      scale_colour_manual(name = 'Bias', values = c( 'blue', 'red'))
  }
  else {
    p1 = qplot(xref, geom = 'blank') +  # Create common set of x values
      geom_line(aes(y = abs(X$rw[[nind]]$biasf), colour = 'Rank Weighted')) + # rank weighted                    
      geom_line(aes(y = abs(X$std[[nind]]$biasf), colour = 'Standard')) +
      ggtitle("Absolute Bias") +
      scale_colour_manual(name = 'Abs. Bias', values = c( 'blue', 'red'))
  }
  # Variance Plot
  p2 = qplot(xref, geom = 'blank') +  # Create common set of x values
    geom_line(aes(y = X$rw[[nind]]$varf, colour = 'Rank Weighted')) + # rank weighted                    
    geom_line(aes(y = X$std[[nind]]$varf, colour = 'Standard')) +
    ggtitle("Variance") +
    scale_colour_manual(name = 'Var.', values = c( 'blue', 'red'))
  
  # Variance Plot
  p3 = qplot(xref, geom = 'blank') +  # Create common set of x values
    geom_line(aes(y = X$rw[[nind]]$msef, colour = 'Rank Weighted')) + # rank weighted                    
    geom_line(aes(y = X$std[[nind]]$msef, colour = 'Standard')) +
    ggtitle("Mean Square Error") +
    scale_colour_manual(name = 'MSE', values = c( 'blue', 'red'))
  
  multiplot(p0,p1,p2,p3)
  
}

#### MISE Estimator ####
# Takes set of points for estimated MSE at several different x values and uses integration rules to
# approximate the area underneath the curve. ideally... this would have several different types of
# rules. Default is trapezoidal. 
# The X is the data input. THe data should come from the list object created in the rot.density.sim
# function. The function outputs a dataframe with the approximated MISE for for standard case and
# rank weighted case at each sample size.
# Current numerical integration rules being used
# "trapezoid" for trapezoidal rule
rot.mise.approx <- function(X, rule ="simpson"){
  require(caTools) #needed for trapz function (trapezoidal rule)
  
  xref <- X$rw.rot[[1]]$xref # Get x values since they are all the same by how sim function works
  n <- names(X[[1]]) # Used from rownames of output dataframe.
  
  mise <- data.frame(rw.rot = rep(NA, length(n)), std.rot = rep(NA, length(n)),
                     diff = rep(NA, length(n)))
  
  for(i in 1:length(n)){
    mise$rw.rot[i] = trapz(xref, X$rw.rot[[i]]$msef)
    mise$std.rot[i] = trapz(xref, X$std.rot[[i]]$msef)
    mise$diff[i] = mise$rw.rot[i] - mise$std.rot[i]
  }
  row.names(mise) <- n
  
  return(mise)
}

#### MISE Estimator ####
# Takes set of points for estimated MSE at several different x values and uses integration rules to
# approximate the area underneath the curve. ideally... this would have several different types of
# rules. Default is trapezoidal. 
# The X is the data input. THe data should come from the list object created in the rot.density.sim
# function. The function outputs a dataframe with the approximated MISE for for standard case and
# rank weighted case at each sample size.
# Current numerical integration rules being used
# "trapezoid" for trapezoidal rule
mise.approx <- function(X, rule ="simpson"){
  require(caTools) #needed for trapz function (trapezoidal rule)
  
  xref <- X$rw[[1]]$xref # Get x values since they are all the same by how sim function works
  n <- names(X[[1]]) # Used from rownames of output dataframe.
  
  mise <- data.frame(rw = rep(NA, length(n)), std = rep(NA, length(n)),
                     diff = rep(NA, length(n)))
  
  for(i in 1:length(n)){
    mise$rw[i] = trapz(xref, X$rw[[i]]$msef)
    mise$std[i] = trapz(xref, X$std[[i]]$msef)
    mise$diff[i] = mise$rw[i] - mise$std[i]
  }
  row.names(mise) <- n
  
  return(mise)
}

#### MISE Estimator for generalized simulationer ####
mise.approx2 <- function(X, rule ="simpson"){
  require(caTools) #needed for trapz function (trapezoidal rule)
  
  xref <- X$rw[[1]][[1]]$xref # Get x values since they are all the same by how sim function works
  n <- names(X[[1]][[1]]) # Used from rownames of output dataframe.
  
  mise <- list(rw = data.frame(), std = data.frame())
  
  #Get names of rank wieghted bandwidth methods
  rw.names <- names(X$rw)
  
  # Get names of standard bandwidth methods
  std.names <- names(X$std)
  
  for(j in 1:length(names(X$rw))){
    
    for(i in 1:length(n)){
      
      mise$rw[i,j] = trapz(xref, X$rw[[j]][[i]]$msef)

    }
  }
  
  names(mise$rw) <- rw.names
  row.names(mise$rw) <- n
  
  for(j in 1:length(names(X$std))){
    
    for(i in 1:length(n)){
      
      mise$std[i,j] = trapz(xref, X$std[[j]][[i]]$msef)
      
    }
  }
  
  names(mise$std) <- std.names
  row.names(mise$std) <- n

  return(mise)
}



#### Kernel FUnction ####
kernel.function2 <- function (u, kernel = "biweight", product = TRUE) 
{
  if (kernel == "triangular") {
    kernel <- "triangle"
  }
  if (kernel == "rectangle" || kernel == "rectangular") {
    kernel <- "uniform"
  }
  if (kernel == "quartic") {
    kernel <- "biweight"
  }
  if (kernel == "normal") {
    kernel <- "gaussian"
  }
  kernel.names <- c("triangle", "uniform", "epanechnikov", 
                    "biweight", "triweight", "gaussian")
  
  d<-1
  c1 <- c(1, 0.5, 0.75, 0.9375, 1.09375, NA)
  pp <- c(1, 2, 2, 2, 2, 0)
  qq <- c(1, 0, 1, 2, 3, NA)
  names(c1) <- names(pp) <- names(qq) <- kernel.names
  p <- pp[kernel]
  q <- qq[kernel]
  volume.d <- pi^(d/2)/gamma(d/2 + 1)
  r1 <- c(d + 1, 1, (d + 2), (d + 2) * (d + 4), (d + 2) * (d + 
                                                             4) * (d + 6), NA)
  r2 <- c(1, 1, 2, 8, 48, NA)
  names(r1) <- names(r2) <- kernel.names
  if (p > 0) {
    if (product) {
      x <- 1 - sqrt(u * u)^p
      c <- c1[kernel]
      k <- (c^d) * x^q *(x >= 0)
    }
    else {
      x <- 1 - sqrt(rowSums(u * u))^p
      c <- r1[kernel]/(r2[kernel] * volume.d)
      k <- c * x^q * (x >= 0)
    }
  }
  else {
    k <- dnorm(u)
  }
  return(k)
}