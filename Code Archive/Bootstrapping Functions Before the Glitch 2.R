
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
boot.samp <- function(x, m = NA, B, suff = F){
	
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
boot.FJ.bw <-function(x, h0 = NA, B = 1000, suff = F, upper.limit = 2){
	
	n <-length(x)
	# Initial Smoothing Parameter
	if(is.na(h0)) h0 <- std.bw(x, "SJ")
	
	# Smoothed bootstrap sample)
	xx <- boot.smooth(x, B = B, h = h0, suff = suff)
	fx <- fhat(x, h0)
	fxx <- lapply(xx, fhat, h = h0) # creating this list is being verrrrry problematic
	
	
	#This function returns the sum of squared bootstrap density deviations from original sample
	# density estimate at a given set of points, i.e., SUM (f_j - f)^2
	# CURRENTLY THIS IS HORRIBLY STRUCTURED BUT IT WORKS
	bigfun<-function(x,h){
		con <- fx(x,h0)
		dxx <- colSums(matrix(unlist(
			lapply(xx, function(c) ((1/n)*colSums((1/h)*dnorm(outer(c, x, '-')/h)) - con)^2 )
		),B, length(x), T) )
		return(dxx)
	}
	
	
	bmise.fj <- Vectorize(function(h){ 
		
		integrate(bigfun, lower = -Inf, upper = Inf, h)$value #SLOW PART ALSO HERE
		
	}, "h")
	
	optimize(bmise.fj, c(0,upper.limit))$minimum # SLOW PART IS HERE
	
	
}


#################################
#### SUFFICIENT BOOTSTRAP BW ####
#################################
# Function for taking a bandwidth of one sample size and rescaling it for that of another sample 
# size. This would be useful in the base of boostrap bandwidth estimation where smaller samples 
# sizes are used.
# See Hall's 1990 paper "Using Bootstrap in estimating the MISE and selecting h in nonpar problems"
# h = h_1(n_1/n)^(1/(2r + 1))
# r is the kernel order. Default order is 2


suff.bw <-function(x, m = NA, B = 1000, suff = F, upper.limit = 2){
	
	m = as.integer(m)
	n <-length(x)
	
	h0 <- std.bw(x, "SJ")
		# Smoothed bootstrap sample)
	xx <- boot.samp(x, m = m, B = B, suff = suff)
	fx <- fhat(x, h0)
	#fxx <- lapply(xx, fhat, h = h0) # creating this list is being verrrrry problematic
	
	
	#This function returns the sum of squared bootstrap density deviations from original sample
	# density estimate at a given set of points, i.e., SUM (f_j - f)^2
	# CURRENTLY THIS IS HORRIBLY STRUCTURED BUT IT WORKS
	bigfun<-function(x,h){
		con <- fx(x,h0)
		dxx <- colSums(matrix(unlist(
			lapply(xx, function(c) {((1/n)*colSums((1/h)*dnorm(outer(c, x, '-')/h)) - con)^2 } )
		),B, length(x), T) )
		return(dxx)
	}
	
	
	bmise.fj <- Vectorize(function(h){ 
		
		integrate(bigfun, lower = -Inf, upper = Inf, h)$value #SLOW PART ALSO HERE
		
	}, "h")
	
	optimize(bmise.fj, c(0,upper.limit))$minimum # SLOW PART IS HERE
	
	
}










#### Function for getting just y values out of density() ####
# supply number of points, from, and to.
density.yonly <- function(x, bw = "nrd0", kernel = "gaussian", n, from, to){
	density(x, kernel = kernel, bw = bw, n = n, from = from, to = to)$y
}


### fhat function ####
# Creates a density estimator that will give the estimated density at given points
# x is the sample data
# h is the desired bandwidth to be used for kde, default is sheather jones bw
# fhat will return a function with one mandatory input and one voluntary input. The mandatory input
# is the set of points that the density is to be estimated at. 
fhat <- function(x, h){
	
	n = length(x)
	ff <- function(xx, hh = h){
		
		# Matrix of differences for estimation points and sample data created by outer();
		# rows indicate which estimation point x_0 we are estimating the density at;
		# columns are each sample point.;
		# Apply kernel to each difference
		# use apply over matrix to get vector of estimated density values
		d = apply((n*hh)^(-1)*dnorm(outer(xx, x, '-')/hh), 1, sum)
		return(d)
		
	}
	return(ff)
	
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
# Bootstrap Bandwidth Function
##################
# Function for taking a bandwidth of one sample size and rescaling it for that of another sample 
# size. This would be useful in the base of boostrap bandwidth estimation where smaller samples 
# sizes are used.
# See Hall's 1990 paper "Using Bootstrap in estimating the MISE and selecting h in nonpar problems"
# h = h_1(n_1/n)^(1/(2r + 1))
# r is the kernel order. Default order is 2

#boot.bw <- function(h1, n1, n, r = 2) return(h1*(n1/n)^(1/(2*r + 1)))



