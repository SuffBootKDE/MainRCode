 
# General Simulations for Bandwidths Function -------------------------------------------------


# M = Simulation Size, n = sample size(s), rsamp and dsamp are sample and density function for
# desired distirbution of simulations, dds is the second derivative of distribution for h_opt
#lower, upper are the bounds for where the reference x points
# are to be made. nref is how many points to estimate the density at.
# ... are the distribution parameters for rsamp and dsamp
# bws is a string indicating bandwidth types to use in the simulation
# possible choices involve bootstrapping, plugin, nrd0, etc...
### NEED TO FIGURE OUT HOW TO GET CIs from bootstrapping.
bw.density.sim <- function(M, n, rsamp, dsamp, dds, kernel = "gaussian", bws, 
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
	
	# Store h! # HARD TO DO! NEED TO REFORMAT EVERYTHING
	
	
	# # Calculate density estimates for rank weighted and standard estimators
	# # Do so for each bandwidth type specified [j]
	# # And at each sample size [i]
	# for(j in 1:length(rw.bws)){
	# 
	# 	disthats$rw[[j]] <- list()
	# 
	# 	for(i in 1:length(n)){
	# 
	# 		disthats$rw[[j]][[i]] <- apply(distsamples[[i]], MARGIN = 2,
	# 												 rw.density, bw = rw.bws[j], xref = xref, nref = nref,
	# 												 kernel = kernel, y.only = T)
	# 	}
	# }
	
	# if(debug) print(" Density Values Calculated:RW")
	
	for(j in 1:length(bws)){
		
		disthats[[j]] <- list()
		
		for(i in 1:length(n)){
			
			disthats[[j]][[i]] <- apply(distsamples[[i]], MARGIN = 2,
													  density.yonly, kernel = kernel, bw = bws[j],
													  n = nref, from = lower, to = upper)
		}
	}
	
	if(debug) print(" Density Values Calculated:STD")
	
	# Name respective estimates by bw type selected
	names(disthats) <- bws
	
	# Get stats for each sample size
	# Create similar list nesting structure as with disthats
	std<- list()
	
	# #vector samples size names vector
	# n.names <- paste("n", n, sep = "")
	# 
	# # Stats calculations for rankweighted case
	# for(j in 1:length(rw.bws)){
	# 	rw[[j]] <- list()
	# 	
	# 	for(i in 1:length(n)){
	# 		rw[[j]][[i]] <- density.stats(disthats$rw[[j]][[i]], xref = xref,
	# 												yref = yref)
	# 	}
	# 	
	# 	names(rw[[j]]) <- n.names
	# 	
	# }
	# 
	# if(debug) print(" density stats for RW.")
	# 
	# names(rw) <- rw.bws
	
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
	names(std) <- bws #### NEED TO DOUBLE CHECK UP TO THIS POINT ####
	#### THEN MODIFY DENSITY CALCULATION FUNCTIONS TO ACCEPT BANDWIDTH SELECTIONS ####
	
	
	return(std)
}



# Normal Mixtures -----------------------------------------------------------------------------

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


