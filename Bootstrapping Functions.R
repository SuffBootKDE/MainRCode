# gplm, pracma, tictoc
library(pracma);library(gplm)


#### Boot Sample Function ####
# Generates a list where each list item is a vector of indicies for a bootstrap sample
# Allows for Sufficient Bootstrapping by setting "suff = T"
# Inputs are:
# n = the original sample size
# m = the size of the bootstrap sample
# B = number of bootstrap samples desired.
# suff = the indicator for whether sufficient bootstrapping is to be used

boot.indices <- function(n, m = NA, B, suff = F){
	
	m = as.integer(m)
	# Future plans... Make sure B is positive and that suff is BOOLEAN.
	# Also, check that everthing is integer type
	
	# Check if m is specified and a valid size
	# If either check fails, set m = n
	if(is.na(m) || m > n) m = n
	
	# Vector of potential indicices
	i = 1:n
	
	# For regular bootstrap, sampling with replacement.
	if(!suff){
		
		#This is a matrix of the indices for all bootstrap samples
		# Each ROW represents the indices for a boot strap sample
		
		m = matrix(sample(i, size = m*B, replace = T),
					  nrow = B, ncol = m, byrow = T)
		
		# Create the list from the matrix.
		# List format is used so that any functions that use this setup can deal
		##with vectors of indices that are not the same length which occurs in the 
		##sufficient bootstrap case.
		ind = lapply(seq_len(nrow(m)), function(i) m[i,])
		
		return(ind)
		
		
	}
	
	
	# For sufficient bootstrap, sample without replacement, then only take unique indices.
	if(suff){
		
		#This is a matrix of the indices for all bootstrap samples
		# Each ROW represents the indices for a boot strap sample
		
		s = matrix(sample(i, size = m*B, replace = T),
					  nrow = B, ncol = m, byrow = T)
		
		# Create the list from the matrix.
		# List format is used so that any functions that use this setup can deal
		##with vectors of indices that are not the same length which occurs in the 
		##sufficient bootstrap case. 
		# THe unique function only selects unique indices in bootstrap sample -> Sufficient Boot
		ind = lapply(seq_len(nrow(s)), function(i) unique(s[i,]))
		
		return(ind)		
	
	}
	
}


#### Bootstrap Sample ####
# Create a list of bootstrap samples: regular or sufficient
# a data vector x is given and other inputs are same as boot.indices
# x = data (no missing values please.) Must be numeric
# m = the size of the bootstrap sample. If unspecified (or too large), you get m = n
# B = number of bootstrap samples desired.
# suff = the indicator for whether sufficient bootstrapping is to be used
boot.samp <- function(x, m = NA, B = 100, suff = F){
	
	n = length(x)
	
	# Check if m is specified and a valid size
	# If either check fails, set m = n
	if(is.na(m) || m > n) m = n
	
	ind = boot.indices(n, m, B, suff)
	samp = lapply(ind, function(i) x[i])
	
	return(samp)
	
}	

#### Smooth Bootstrap Sample ####
# Create a list of smoothed bootstrap samples
# a data vector x is given and other inputs are same as boot.indices
# x = data (no missing values please.) Must be numeric
# m = the size of the bootstrap sample. If unspecified (or too large), you get m = n
# B = number of bootstrap samples desired.
# suff = the indicator for whether sufficient bootstrapping is to be used
# h = a baseline bandwidth used for kde. default to Sheather Jones
# Note the "smoothing" is done assuming a Gaussian kernel with var = 1
boot.smooth <- function(x, B, m = NA, h = NA, suff = F){
	
	# Default to Sheather Jones bandwidth
	if(is.na(h)) h = std.bw(x, "SJ")
	
	n = length(x)
	xbar = mean(x)
	varx = var(x)
	
	# Check if m is specified and a valid size
	# If either check fails, set m = n
	if(is.na(m) || m > n) m = n
	
	ind = boot.indices(n, m, B, suff)
	samp = lapply(ind, function(i) xbar + ( x[i] - xbar +h*rnorm(1))/sqrt(1 + h^2/varx) )
	
	return(samp)
	
}	



#### SECOND BOOTSTRAP BW OPTIMIZATION TRY####
# Calculate bootstrap based bandwidth.
# x is sample data
# h0 is pilot bandwidth
# B is bootstrap sample size.
# upper.limit is the upperbound for the optimize call.

# This currently only works for regular smoothed bootstrap. Need to fix sufficient bootstrap code
boot.FJ.bw <- function(x, h0 = NA, B = 1000, suff = F, upper.limit = 2){
	
	n <-length(x)
	# Initial Smoothing Parameter
	if(is.na(h0)) h0 <- std.bw(x, "SJ")
	
	# Smoothed bootstrap sample. Smoothing based on intial h.
	xx <- boot.smooth(x, B = B, h = h0, suff = suff)
	
	# we need a function to optimize
	
	
### NEED CODE ###
	
	
}


#################################
#### SUFFICIENT BOOTSTRAP BW ####
#################################
# x is the data
# m is reduced bootstrap sample size, m = NA if resampling size is to be n
# B is number of bootstrap samples
# suff is boolean for doing sufficient boot (T) or not (F)
# lower and upper are the bounds for MISE integral
# opt.limit is the upperbound in optimization for bandwidth

suff.bw <-function(x, m = NA, B = 100, suff = T, lower = -5, upper = 5, opt.limit = 10){
	require(pracma)
   
   
	m = as.integer(m)
	n <-length(x)
	n0 = n
	
	h0 <- std.bw(x, "SJ")
		# Smoothed bootstrap sample)
	xx <- boot.samp(x, m = m, B = B, suff = suff)
	### NEED CODE ###
	
}









#### Function for getting just y values out of density() ####
# supply number of points, from, and to.
density.yonly <- function(x, bw = "nrd0", kernel = "gaussian", n, from, to){
	density(x, kernel = kernel, bw = bw, n = n, from = from, to = to)$y
}


### fhat function (CHANGED) #### 
# Creates a density estimator that will give the estimated density at given points
# x is the sample data
# the function returned has two arguments. x0 and h
# x0 is the set of estimation points desired
# h is bandwidth
fhat <- function(x){
	
	n = length(x)
	fh <- function(x0, h){
		
		# Matrix of differences for estimation points and sample data created by outer();
		# rows indicate which estimation point x_0 we are estimating the density at;
		# columns are each sample point.;
		# Apply kernel to each difference
		# use apply over matrix to get vector of estimated density values
		d = apply((n*h)^(-1)*dnorm(outer(x0, x, '-')/h), 1, sum)
		return(d)
		
	}
	return(fh)
	
}

#### BMISE Sq. Diff. ####
# This functions takes the 
bmise.g <- function(x, xx, h){
	

	
}




##################
# Kernel Bandwidth Function
##################
# Function for calculating Std Bws. Equivalent to rw.bws directly above #
# Coordinates all baseline bandwidth functions into one
# nrd0 (normal rule of thumb)
# SJ (Sheather-Jones plug-in bw)
# ucv (Unbiased or Least Squares Cross Validation)
std.bw <- function(x, bw = "nrd0"){
	require(gplm)
	
	#I don't think I need these here. Let's call it legacy code that I'm too scared to remove.
#	x <- sort(x)
#	n <- length(x)
#	r <- 1:n # Ranks since x is now sorted
#	s <- sd(x)
	
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



##################
# Optimal Bandwidth Calculator
##################
# Take reference density and calculate optimal bandwidth for that density
# 
# K is kernel function used, normal by default
# n is sample size used
# dd is second derivative of estimated density, default is normal
# m and v are mean and variance of kernel,
# Lower and upper bounds of integration are Inf by default
h.opt.fn <- function(K = dnorm, n, dd = function(x){(x^2-1)*dnorm(x)}, m = NA, v = NA, 
							l = Inf, u = Inf){
	
	if(is.na(m)){m = integrate(function(x){ x*K(x)}, l, u)$value}
	if(is.na(v)){v = integrate(function(x){ (x-m)^2*K(x)}, l, u)$value} 
	hamise = ( r.int(K)/(v^2*r.int(dd)*n) )^(1/5)
	return(hamise)
}





#### R( . ) ####
# Calculate squared integration of function #
r.int <- function(f, ..., lower = -Inf, upper = Inf){
	require(pracma)
	f2 <- function(x,...) f(x,...)^2
	
	integrate(f2,..., lower, upper,subdivisions = 1000)$value
}


#### Equivalent function for calculating  standard bandwidths####
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
